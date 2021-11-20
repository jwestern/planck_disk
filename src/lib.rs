#![allow(unused)]
use std::ops::Range;
use pyo3::{PyErrArguments, exceptions::PyValueError, prelude::*};

enum EquationOfState {
    Gas { adiabatic_index: f64 },
    Radiation,
    GasPlusRadiation,
}

impl EquationOfState {
    fn temperature(&self, density: f64, pressure: f64) -> f64 {
        match self {
            Self::Gas { .. } => 0.0,
            Self::Radiation => 0.0,
            Self::GasPlusRadiation => 0.0,
        }
    }

    fn pressure(&self, density: f64, temperature: f64) -> f64 {
        match self {
            Self::Gas { .. } => 0.0,
            Self::Radiation => 0.0,
            Self::GasPlusRadiation => 0.0,
        }
    }
}

/// Returns the disk photosphere temperature of a disk, given the vertically
/// integrated mass density, the vertically integrated pressure, and the
/// opacity, kappa.
fn photosphere_temperature(surface_density: f64, surface_pressure: f64, kappa: f64) -> f64 {
    0.0
}

/// Returns the luminosity per unit surface area (radiant flux), of a
/// blackbody with the given temperature, over a frequency range.
fn surface_luminosity(photosphere_temperature: f64, frequency_range: &Range<f64>) -> f64 {
    0.0
}

#[pyfunction]
/// Compute the blackbody radiant flux for a patch on a thin disk as a
/// function of the vertically integrated gas pressure and the surface
/// density.
fn radiant_flux(surface_pressure: Vec<f64>, surface_density: Vec<f64>) -> PyResult<Vec<f64>> {
    if surface_pressure.len() != surface_density.len() {
        return Err(PyValueError::new_err("input arrays must have the same length"))
    }

    let eos = EquationOfState::Gas {
        adiabatic_index: 5.0 / 3.0,
    };

    let scale_height = 1.0;
    let kappa = 1.0;
    let frequency_range = 1e12..1e13;

    let num_zones = surface_pressure.len();
    let mut flux = vec![0.0; num_zones];

    for i in 0..num_zones {
        let rho = surface_density[i] / scale_height;
        let pre = surface_pressure[i] / scale_height;
        let midplane_temp = eos.temperature(rho, pre);
        let photosphere_temp = photosphere_temperature(surface_density[i], surface_pressure[i], kappa);
        flux[i] = surface_luminosity(photosphere_temp, &frequency_range);
    }
    Ok(flux)
}

#[pymodule]
fn planck_disk(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(radiant_flux, m)?)?;
    Ok(())
}
