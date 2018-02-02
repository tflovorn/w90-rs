use std::path::Path;
use std::io;
use std::io::Write;
use std::fs::File;
use qe::pw::input::generate_uniform_kpoints;
use input;
use input::{AngularMomentum, Disentanglement, Input, LatticeUnits, MLWFIterationMode,
            PositionCoordinateType, Projection, ProjectionSite};

pub fn make_input_file(input: &Input) -> Result<String, Error> {
    input::validate(&input)?;

    let header = make_header(&input);

    let mut input_sections = vec![header];

    if let Some(ref disentanglement) = input.disentanglement {
        input_sections.push(make_disentanglement(disentanglement));
    }

    let projections = make_projections(&input);
    let cell = make_unit_cell(&input);
    let positions = make_positions(&input);
    let k_points = make_kpoints(&input);

    input_sections.extend(vec![projections, cell, positions, k_points]);

    let input_text = input_sections.join("\n");

    return Ok(input_text);
}

fn make_header(input: &Input) -> String {
    let mut lines = Vec::new();

    lines.push(format!("num_bands = {}", input.num_bands));
    lines.push(format!("num_wann = {}", input.num_wann));
    lines.push(format!("num_iter = {}", input.mlwf_iteration_mode.value()));

    push_bool_field(&mut lines, "num_iter", input.write_hr);

    lines.join("\n")
}

fn push_bool_field(lines: &mut Vec<String>, name: &str, b: Option<bool>) {
    if let Some(b) = b {
        let val = if b {
            String::from(".true.")
        } else {
            String::from(".false.")
        };

        lines.push(format!("{}={}", name, val));
    };
}

fn make_disentanglement(dis: &Disentanglement) -> String {
    let mut lines = Vec::new();

    lines.push(format!("dis_win_min = {}", dis.dis_win_min));
    lines.push(format!("dis_win_max = {}", dis.dis_win_max));
    lines.push(format!("dis_froz_min = {}", dis.dis_froz_min));
    lines.push(format!("dis_froz_max = {}", dis.dis_froz_max));
    lines.push(format!("dis_num_iter = {}", dis.dis_num_iter));
    lines.push(format!("dis_mix_ratio = {}", dis.dis_mix_ratio));

    lines.join("\n")
}

fn make_projections(input: &Input) -> String {
    let mut lines = Vec::new();

    push_bool_field(&mut lines, "spinors", Some(input.spinors));

    lines.push(String::from("begin projections"));

    if let Some(ref units) = input.projection_units {
        lines.push(units.value());
    }
    for proj in &input.projections {
        lines.push(proj.value());
    }

    lines.push(String::from("end projections"));

    lines.join("\n")
}

fn make_unit_cell(input: &Input) -> String {
    let cell = &input.unit_cell_cart;
    let mut lines = vec![String::from("begin unit_cell_cart"), cell.units.value()];

    for i in 0..3 {
        let a = cell.cell[i];
        lines.push(format!("  {}  {}  {}", a[0], a[1], a[2]));
    }

    lines.push(String::from("end unit_cell_cart"));

    lines.join("\n")
}

fn make_positions(input: &Input) -> String {
    let mut lines = Vec::new();

    let pos = &input.positions;
    match pos.coordinate_type {
        PositionCoordinateType::BohrCartesian => {
            lines.push(String::from("begin atoms_cart"));
            lines.push(String::from("bohr"));
        }
        PositionCoordinateType::AngstromCartesian => {
            lines.push(String::from("begin atoms_cart"));
            lines.push(String::from("ang"));
        }
        PositionCoordinateType::Crystal => {
            lines.push(String::from("begin atoms_frac"));
        }
    }

    for coord in &pos.coordinates {
        let r = &coord.r;
        lines.push(format!(" {} {} {} {}", coord.species, r[0], r[1], r[2]));
    }

    match pos.coordinate_type {
        PositionCoordinateType::BohrCartesian | PositionCoordinateType::AngstromCartesian => {
            lines.push(String::from("end atoms_cart"));
        }
        PositionCoordinateType::Crystal => {
            lines.push(String::from("end atoms_frac"));
        }
    }

    lines.join("\n")
}

fn make_kpoints(input: &Input) -> String {
    let nk = input.k_points;

    let mut lines = vec![
        format!("mp_grid = {} {} {}", nk[0], nk[1], nk[2]),
        String::from("begin kpoints"),
    ];

    for k in generate_uniform_kpoints(nk) {
        lines.push(format!("{} {} {}", k[0], k[1], k[2]));
    }

    lines.push(String::from("end kpoints"));
    lines.join("\n")
}

pub fn write_input_file<P: AsRef<Path>>(input: &Input, file_path: P) -> Result<(), Error> {
    let input_text = make_input_file(input)?;

    let mut file = File::create(file_path)?;
    file.write_all(input_text.as_bytes())?;

    Ok(())
}

#[derive(Fail, Debug)]
pub enum Error {
    #[fail(display = "{}", _0)] Input(input::ErrorList),
    #[fail(display = "{}", _0)] Io(#[cause] io::Error),
}

impl From<input::ErrorList> for Error {
    fn from(errs: input::ErrorList) -> Error {
        Error::Input(errs)
    }
}

impl From<io::Error> for Error {
    fn from(e: io::Error) -> Error {
        Error::Io(e)
    }
}

/// A `Field` has a method `value()` which returns its textual representation on the
/// right-hand side of a `field_name = value` expression in the Wannier90 input file.
pub trait Field {
    fn value(&self) -> String;
}

impl Field for MLWFIterationMode {
    fn value(&self) -> String {
        format!(
            "{}",
            match *self {
                MLWFIterationMode::ProjectionOnly => 0,
                MLWFIterationMode::MLWF { num_iter } => num_iter,
            }
        )
    }
}

impl Field for LatticeUnits {
    fn value(&self) -> String {
        String::from(match *self {
            LatticeUnits::Bohr => "bohr",
            LatticeUnits::Angstrom => "ang",
        })
    }
}

impl Field for Projection {
    fn value(&self) -> String {
        match *self {
            Projection::Random => String::from("random"),
            Projection::Site {
                ref site,
                ref ang_mtm,
                zaxis,
                xaxis,
                radial,
                zona,
            } => {
                let mut proj = format!("{}:", site.value());
                for (i, ang_mtm) in ang_mtm.iter().enumerate() {
                    if i != 0 {
                        proj.push_str(";");
                    }
                    proj.push_str(&ang_mtm.value());
                }
                if let Some(zaxis) = zaxis {
                    proj.push_str(&format!(":z={},{},{}", zaxis[0], zaxis[1], zaxis[2]));
                };
                if let Some(xaxis) = xaxis {
                    proj.push_str(&format!(":x={},{},{}", xaxis[0], xaxis[1], xaxis[2]));
                };
                if let Some(radial) = radial {
                    proj.push_str(&format!(":r={}", radial));
                };
                if let Some(zona) = zona {
                    proj.push_str(&format!(":zona={}", zona));
                };
                proj
            }
        }
    }
}

impl Field for ProjectionSite {
    fn value(&self) -> String {
        match *self {
            ProjectionSite::Species(ref species) => species.clone(),
            ProjectionSite::CenterCartesian(r) => format!("c = {}, {}, {}", r[0], r[1], r[2]),
            ProjectionSite::CenterCrystal(r) => format!("f = {}, {}, {}", r[0], r[1], r[2]),
        }
    }
}

impl Field for AngularMomentum {
    fn value(&self) -> String {
        String::from(match *self {
            AngularMomentum::S => "l=0",
            AngularMomentum::P => "l=1",
            AngularMomentum::D => "l=2",
            AngularMomentum::F => "l=3",
        })
    }
}
