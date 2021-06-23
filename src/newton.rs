pub fn newton<F: Fn(f64)->f64>(f: F, x0: f64, dx: f64, tol: f64) -> f64 {
    let mut x = x0;
    loop {
        let fx = f(x);
        if fx.abs() < tol {
            return x;
        }
        let fx_prime = (f(x + dx) - fx) / dx;
        x -= fx / fx_prime;
    }
}