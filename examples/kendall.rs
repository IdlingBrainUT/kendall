use kendall::kendall::*;
use chrono::*;
use rand::distributions::{Uniform, Distribution};

fn main() {
    let mut rng = rand::thread_rng();
    let ud = Uniform::new(0.0, 1.0);
    let a = (0..1000000).map(|_| if ud.sample(&mut rng) < 0.5 { 0.0 } else { ud.sample(&mut rng) }).collect::<Vec<f64>>();
    let b = (0..1000000).map(|_| if ud.sample(&mut rng) < 0.5 { 0.0 } else { ud.sample(&mut rng) }).collect::<Vec<f64>>();

    /*
    let start = Local::now();
    println!("{}", tau(&a, &b).unwrap());
    println!("{}", Local::now() - start);
    */
    
    let start = Local::now();
    println!("=== SDTAU ===");
    println!("tau = {:.5}", sdtau(&a, &b).unwrap());
    println!("time = {:.5}\n", Local::now() - start);

    let start = Local::now();
    println!("=== NZTAU ===");
    println!("tau = {:.5}", nztau(&a, &b).unwrap());
    println!("time = {:.5}\n", Local::now() - start);

    let start = Local::now();
    println!("=== NDNZTAU ===");
    println!("tau = {:.5}", ndnztau(&a, &b).unwrap());
    println!("time = {:.5}\n", Local::now() - start);

    println!("tau = {:.5}", ndnztau_scale(&a, &b).unwrap());
}