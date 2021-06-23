use statrs::distribution::{ContinuousCDF, Normal};

use crate::newton::*;

pub fn sn_cdf_powi(x: f64, n: i32) -> f64 {
    let sn = Normal::new(0.0, 1.0).unwrap();
    let sn_cdf_x = sn.cdf(x);
    sn_cdf_x.powi(n)
}

pub fn sn_icdf_powi(x: f64, n: i32) -> f64 {
    if n == 0 {
        panic!();
    } else if n == 1 {
        let sn = Normal::new(0.0, 1.0).unwrap();
        sn.inverse_cdf(x)
    } else {
        let f = |y: f64| sn_cdf_powi(y, n) - x;
        newton(f, 0.0, 1e-7, 1e-7)
    }
}