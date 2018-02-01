extern crate w90;

use w90::input;
use w90::input::{AngularMomentum, AtomCoordinate, Cell, Disentanglement, LatticeUnits,
                 MLWFIterationMode, PositionCoordinateType, Positions, Projection, ProjectionSite};
use w90::serialize;

#[test]
fn generate_input() {
    let disentanglement = Some(Disentanglement {
        dis_win_min: -6.5582,
        dis_win_max: 8.4418,
        dis_froz_min: -4.5582,
        dis_froz_max: 6.4418,
        dis_num_iter: 1000,
        dis_mix_ratio: 0.5,
    });

    let projections = vec![
        Projection::Site {
            site: ProjectionSite::Species(String::from("Se")),
            ang_mtm: vec![AngularMomentum::P],
            zaxis: None,
            xaxis: None,
            radial: None,
            zona: None,
        },
        Projection::Site {
            site: ProjectionSite::Species(String::from("W")),
            ang_mtm: vec![AngularMomentum::D],
            zaxis: None,
            xaxis: None,
            radial: None,
            zona: None,
        },
    ];

    let unit_cell_cart = Cell {
        units: LatticeUnits::Bohr,
        cell: [
            [3.13603975949, -5.43178019799, 0.0],
            [3.13603975949, 5.43178019799, 0.0],
            [0.0, 0.0, 68.6629186029],
        ],
    };

    let positions = Positions {
        coordinate_type: PositionCoordinateType::Crystal,
        coordinates: vec![
            AtomCoordinate {
                species: String::from("Se"),
                r: [0.0, 0.0, 0.275217856494],
            },
            AtomCoordinate {
                species: String::from("W"),
                r: [0.333333333333, 0.666666666667, 0.321438654707],
            },
            AtomCoordinate {
                species: String::from("Se"),
                r: [0.0, 0.0, 0.36765945292],
            },
        ],
    };

    let test_input = input::Input {
        num_bands: 44,
        num_wann: 22,
        write_hr: Some(true),
        mlwf_iteration_mode: MLWFIterationMode::ProjectionOnly,
        disentanglement,
        spinors: true,
        projection_units: None,
        projections,
        unit_cell_cart,
        positions,
        kpoints: [9, 9, 1],
    };

    let input_text = serialize::make_input_file(&test_input).unwrap();

    println!("{}", input_text);
}
