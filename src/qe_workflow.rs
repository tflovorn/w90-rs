use qe::pw::input::Input as PwInput;
use qe::pw::input::{Calculation, Ibrav, KPoints, Occupations, Smearing, SpinType};
use qe::pw::input::LatticeUnits as PwLatticeUnits;
use qe::pw::input::PositionCoordinateType as PwCoord;
use qe::pw::input::AtomCoordinate as PwAtomCoordinate;
use input::Input as W90Input;
use input::{Disentanglement, MLWFIterationMode, Projection};
use input::LatticeUnits as W90LatticeUnits;
use input::PositionCoordinateType as W90Coord;
use input::AtomCoordinate as W90AtomCoordinate;
use input::Positions as W90Positions;
use input::Cell as W90Cell;

pub fn nscf_input(
    scf: &PwInput,
    diago_thr_init: f64,
    num_bands: u64,
    smearing_type: &Smearing,
    smearing_size: f64,
    nscf_nk: [u64; 3],
) -> Result<PwInput, Error> {
    if let Calculation::Scf { .. } = scf.calculation {
        let mut nscf = scf.clone();
        nscf.calculation = Calculation::Nscf {
            diago_thr_init,
            nbnd: Some(num_bands),
            nosym: Some(true),
        };
        nscf.system.occupations = Occupations::Smearing(smearing_type.clone(), smearing_size);
        nscf.k_points = KPoints::CrystalUniform(nscf_nk);

        Ok(nscf)
    } else {
        Err(Error::WrongCalculation)
    }
}

pub fn bands_input(nscf: &PwInput, bands_kpoints: &KPoints) -> Result<PwInput, Error> {
    let (diago_thr_init, nbnd, nosym) = match nscf.calculation {
        Calculation::Nscf {
            diago_thr_init,
            nbnd,
            nosym,
        } => Ok((diago_thr_init, nbnd, nosym)),
        _ => Err(Error::WrongCalculation),
    }?;
    match nosym {
        Some(nosym) => {
            if !nosym {
                return Err(Error::NoSym);
            }
        }
        None => {
            return Err(Error::NoSym);
        }
    };

    let mut bands = nscf.clone();
    bands.calculation = Calculation::Bands {
        diago_thr_init,
        nbnd,
        nosym,
    };

    bands.k_points = match *bands_kpoints {
        KPoints::CrystalBands { .. } => Ok(bands_kpoints.clone()),
        _ => Err(Error::WrongKPointsBands),
    }?;

    Ok(bands)
}

pub fn w90_input(
    nscf: &PwInput,
    num_wann: u64,
    mlwf_iteration_mode: &MLWFIterationMode,
    disentanglement: &Disentanglement,
    projection_units: Option<W90LatticeUnits>,
    projections: Vec<Projection>,
) -> Result<W90Input, Error> {
    let num_bands = match nscf.calculation {
        Calculation::Nscf { nbnd, .. } => {
            if let Some(nbnd) = nbnd {
                Ok(nbnd)
            } else {
                Err(Error::NumBands)
            }
        }
        _ => Err(Error::WrongCalculation),
    }?;

    let spinors = match nscf.system.spin_type {
        Some(ref spin_type) => match spin_type {
            &SpinType::NonPolarized | &SpinType::CollinearPolarized => false,
            &SpinType::Noncollinear { .. } => true,
        },
        None => false,
    };

    let (lattice_units, cell) = match nscf.system.ibrav {
        Ibrav::Free(ref cell) => match cell.units {
            PwLatticeUnits::Alat => {
                let alat = nscf.system.alat;
                let lat_vecs = scale_cell(cell.cell, alat);
                (W90LatticeUnits::Bohr, lat_vecs)
            }
            PwLatticeUnits::Bohr => (W90LatticeUnits::Bohr, cell.cell),
            PwLatticeUnits::Angstrom => (W90LatticeUnits::Angstrom, cell.cell),
        }, // TODO - support other Ibrav cases.
           // Here we don't have a simple way to extract the lattice vectors from PwInput.
           // May need to generate it by hand, or just leave unsupported for w90_input.
           // Another possibilty: extract from scf output. Then must run this step after
           // scf finishes, though.
    };
    let unit_cell_cart = W90Cell {
        units: lattice_units,
        cell,
    };

    let (coordinate_type, coordinates) = match nscf.atomic_positions.coordinate_type {
        PwCoord::AlatCartesian => {
            let alat = nscf.system.alat;
            let coordinates = scale_coords(&nscf.atomic_positions.coordinates, alat);
            Ok((W90Coord::BohrCartesian, coordinates))
        }
        PwCoord::BohrCartesian => Ok((
            W90Coord::BohrCartesian,
            map_coords(&nscf.atomic_positions.coordinates),
        )),
        PwCoord::AngstromCartesian => Ok((
            W90Coord::AngstromCartesian,
            map_coords(&nscf.atomic_positions.coordinates),
        )),
        PwCoord::Crystal => Ok((
            W90Coord::Crystal,
            map_coords(&nscf.atomic_positions.coordinates),
        )),
        PwCoord::CrystalSG => Err(Error::CrystalSG),
    }?;
    let positions = W90Positions {
        coordinate_type,
        coordinates,
    };

    let k_points = match nscf.k_points {
        KPoints::CrystalUniform(k_points) => Ok(k_points),
        _ => Err(Error::WrongKPointsNscf),
    }?;

    Ok(W90Input {
        num_bands,
        num_wann,
        write_hr: Some(true),
        mlwf_iteration_mode: mlwf_iteration_mode.clone(),
        disentanglement: Some(disentanglement.clone()),
        spinors,
        projection_units,
        projections,
        unit_cell_cart,
        positions,
        k_points,
    })
}

#[derive(Fail, Debug)]
pub enum Error {
    #[fail(display = "Unexpected type of calculation input.")] WrongCalculation,
    #[fail(display = "Must have `nosym = true` in nscf calculation.")] NoSym,
    #[fail(display = "Must specify `nbnd` in nscf calculation.")] NumBands,
    #[fail(display = "Must have `KPoints::CrystalUniform` in nscf calculation.")] WrongKPointsNscf,
    #[fail(display = "Must input `KPoints::CrystalBands`.")] WrongKPointsBands,
    #[fail(display = "`CrystalSG` positions unsupported.")] CrystalSG,
}

fn scale_cell(cell: [[f64; 3]; 3], alat: f64) -> [[f64; 3]; 3] {
    // TODO - clean way to do this? Can't do collect::<[f64; 3]>() etc.
    // Would be simple if using an Array2 instead of nested array, but there
    // can't specify static size.
    // Const generics would fix this problem.
    let a0 = [alat * cell[0][0], alat * cell[0][1], alat * cell[0][2]];
    let a1 = [alat * cell[1][0], alat * cell[1][1], alat * cell[1][2]];
    let a2 = [alat * cell[2][0], alat * cell[2][1], alat * cell[2][2]];

    [a0, a1, a2]
}

fn scale_coords(coordinates: &Vec<PwAtomCoordinate>, alat: f64) -> Vec<W90AtomCoordinate> {
    coordinates
        .iter()
        .map(|c| W90AtomCoordinate {
            species: c.species.clone(),
            r: [alat * c.r[0], alat * c.r[1], alat * c.r[2]],
        })
        .collect()
}

fn map_coords(coordinates: &Vec<PwAtomCoordinate>) -> Vec<W90AtomCoordinate> {
    coordinates
        .iter()
        .map(|c| W90AtomCoordinate {
            species: c.species.clone(),
            r: c.r,
        })
        .collect()
}
