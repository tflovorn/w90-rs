#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Input {
    pub num_bands: u64,
    pub num_wann: u64,
    pub write_hr: Option<bool>,
    pub mlwf_iteration_mode: MLWFIterationMode,

    pub disentanglement: Option<Disentanglement>,

    pub spinors: bool,
    pub projection_units: Option<LatticeUnits>,
    pub projections: Vec<Projection>,

    pub unit_cell_cart: Cell,
    pub positions: Positions,

    pub kpoints: [u64; 3],
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum MLWFIterationMode {
    /// Use projected wavefunctions only: do not mix to acheive maximal localization.
    /// `num_iter` is set to 0 in Wannier90 input file.
    ProjectionOnly,
    /// Mix projected wavefunctions to achieve maximal localization, using a maximum
    /// of `num_iter` iterations.
    MLWF { num_iter: u64 },
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Disentanglement {
    pub dis_win_min: f64,
    pub dis_win_max: f64,
    pub dis_froz_min: f64,
    pub dis_froz_max: f64,
    pub dis_num_iter: u64,
    pub dis_mix_ratio: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum Projection {
    Random,
    Site {
        site: ProjectionSite,
        ang_mtm: Vec<AngularMomentum>,
        zaxis: Option<[f64; 3]>,
        xaxis: Option<[f64; 3]>,
        radial: Option<u64>,
        zona: Option<f64>,
    },
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum ProjectionSite {
    /// Projection centered on all atoms of the given species.
    Species(String),
    /// `c = x, y, z` type of projection centered at Cartesian position (x, y, z).
    CenterCartesian([f64; 3]),
    /// `f = r1, r2, r3` type of projection centered at lattice coordinate position (r1, r2, r3).
    CenterCrystal([f64; 3]),
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum AngularMomentum {
    S,
    P,
    D,
    F,
    // TODO: hybrid orbitals and individual l=l,mr=mr orbitals.
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Cell {
    pub units: LatticeUnits,
    pub cell: [[f64; 3]; 3],
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum LatticeUnits {
    Bohr,
    Angstrom,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Positions {
    pub coordinate_type: PositionCoordinateType,
    pub coordinates: Vec<AtomCoordinate>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum PositionCoordinateType {
    BohrCartesian,
    AngstromCartesian,
    Crystal,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AtomCoordinate {
    pub species: String,
    pub r: [f64; 3],
}

pub fn validate(input: &Input) -> Result<(), ErrorList> {
    let mut errs = Vec::new();

    // TODO: Check that all atom-centered projections have Species that exist in the coordinates.

    // Check that `Random` does not appear more than once in the list of projections.
    let random_count = input
        .projections
        .iter()
        .filter(|&p| *p == Projection::Random)
        .count();
    if random_count > 1 {
        errs.push(Error::RandomCount);
    }

    // TODO: Check that the number of projections given is compatible with `num_wann`, or that
    // `Random` is in the list of projections.
    //if random_count == 0 {
    //
    //}

    if errs.len() == 0 {
        Ok(())
    } else {
        Err(ErrorList { errs })
    }
}

#[derive(Fail, Debug)]
pub enum Error {
    #[fail(display = "`Random` may appear at most once in the list of projections.")] RandomCount,
    //#[fail(display = "Number of projections is incompatible with `num_wann`.")]
    //ProjectionNumber,
}

pub type ErrorList = ::qe::error::ErrorList<Error>;
