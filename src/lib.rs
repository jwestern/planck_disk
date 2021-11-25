#![allow(unused)]
use pyo3::{exceptions::PyValueError, prelude::*, PyErrArguments};
use std::ops::Range;
use std::f64::consts::PI;

// Physical constants in CGS units
pub const SPEED_OF_LIGHT: f64 = 2.99792458e10;
pub const K_BOLTZMANN: f64 = 1.38065812e-16;
pub const M_PROTON: f64 = 1.6726e-24;
pub const G_NEWTON: f64 = 6.6725985e-8;
pub const SIGMA_STEFAN_BOLTZMANN: f64 = 5.6705119e-5;
pub const PARSEC: f64 = 3085678000000000000.0;
pub const SOLAR_MASS: f64 = 1.989e33;
pub const PLANCK_CONSTANT: f64 = 6.6261e-27;
pub const NPLANCK: i32 = 15; // how many terms to use when approximating the antiderivative of the Planck spectrum

// Auxiliary function used to approximate the antiderivative of the Planck spectrum
fn planck_aux(x: f64, n: i64) -> f64 {
    ( x.powi(3) / n + 3.0 * x.powi(2) / n.powi(2) + 6.0 * x / n.powi(3) + 6.0 / n.powi(4) ) * f64::exp( -n * x )
}

// Approximation of the antiderivative of the Planck spectrum
fn planck_anti_approx(x: f64, m: i64} -> f64 {
    let mut result: f64 = 0.0;
    for n in 0..m {
        result = result + planck_aux(x, n);
    }
    result
}

enum EquationOfState {
    Gas { adiabatic_index: f64 },
    Radiation,
    GasPlusRadiation,
}

impl EquationOfState {
    /// TODO
    fn temperature(&self, density: f64, pressure: f64) -> f64 {
        match self {
            Self::Gas { .. } => pressure / density * mproton / kboltzmann,
            Self::Radiation => 0.0,
            Self::GasPlusRadiation => 0.0,
        }
    }

    /// TODO
    fn pressure(&self, density: f64, temperature: f64) -> f64 {
        match self {
            Self::Gas { .. } => density * temperature * kboltzmann / mproton,
            Self::Radiation => 0.0,
            Self::GasPlusRadiation => 0.0,
        }
    }
}

/// Returns the disk photosphere temperature of a disk, given the vertically
/// integrated mass density, the vertically integrated pressure, and the
/// opacity, kappa.
fn photosphere_temperature(surface_density: f64, surface_pressure: f64, kappa: f64) -> f64 {
    0.0 //need temperature and density, then use (4/3)*sigma*T^4/(kappa*Sigma)    
}

/// Returns the luminosity per unit surface area (radiant flux), of a
/// blackbody with the given temperature, over a frequency range.
fn surface_luminosity(photosphere_temperature: f64, frequency_range: &Range<f64>) -> f64 {
    // This is the expensive function which integrates the Planck spectrum.
    // When computing multiple bands, it's more efficient to do them all at once,
    // because neighboring bands share an antiderivative at their mutual boundary.
    let xl = energy_ratio(frequency_range[0], photosphere_temperature);
    let xr = energy_ratio(frequency_range[1], photosphere_temperature);
    let planck_antideriv_l = planck_anti_approx(xl, NPLANCK);
    let planck_antideriv_r = planck_anti_approx(xr, NPLANCK);
    let prefactor = 2.0 / SPEED_OF_LIGHT.powi(2) * ( K_BOLTZMANN * photosphere_temperature ).powi(4) / PLANCK_CONSTANT.powi(3)
    PI * prefactor * ( planck_antideriv_r - planck_antideriv_l )
}

#[pyfunction(
    surface_pressure = "vec![]",
    surface_density = "vec![]",
    frequency_range = "(1e14, 1e15)",
    eos = "\"gamma_law\"",
    kappa = "0.4",
)]
#[pyo3(text_signature = "(surface_pressure, surface_density, frequency_range=(1e14, 1e15), eos='gamma_law', kappa=0.4)")]
/// Compute the blackbody radiant flux for a patch on a thin disk as a
/// function of the vertically integrated gas pressure and the surface
/// density. Arguments are assumed to be in CGS units.
fn radiant_flux(
    surface_pressure: Vec<f64>,
    surface_density: Vec<f64>,
    frequency_range: (f64, f64),
    eos: &str,
    kappa: f64,
) -> PyResult<Vec<f64>> {

    if surface_pressure.len() != surface_density.len() {
        return Err(PyValueError::new_err(
            "input arrays must have the same length ",
        ));
    }

    let eos = match eos {
        "gamma_law" => EquationOfState::Gas {
            adiabatic_index: 5.0 / 3.0,
        },
        _ => {
            return Err(PyValueError::new_err("eos must be 'gamma_law'"))
        }        
    };

    let scale_height = 1.0; // TODO
    let frequency_range = frequency_range.0..frequency_range.1;

    let num_zones = surface_pressure.len();
    let mut flux = vec![0.0; num_zones];

    for i in 0..num_zones {
        let rho = surface_density[i] / scale_height;
        let pre = surface_pressure[i] / scale_height;
        let midplane_temp = eos.temperature(rho, pre);
        let photosphere_temp =
            photosphere_temperature(surface_density[i], surface_pressure[i], kappa);
        flux[i] = surface_luminosity(photosphere_temp, &frequency_range);
    }
    Ok(flux)
}

#[pymodule]
fn planck_disk(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(radiant_flux, m)?)?;
    Ok(())
}
