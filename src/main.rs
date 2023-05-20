use clap::Parser;
use rand::{distributions::Standard, Rng};
use std::{error::Error, fs, path::PathBuf, str::FromStr};
use serde::Deserialize;

fn main() -> Result<(), Box<dyn Error>> {
    let args = Arguments::parse();
    let config_file = fs::read_to_string(args.config_file)?;
    let config: Config = toml::from_str(&config_file)?;
    
    let j_beta = config.coupling / config.beta;
    let h_beta = config.field / config.beta;
    
    // Spin state of the system, the value 1 = spin up, and the value -1 = spin down.
    let mut z = vec![1; config.sites];
    for _ in 0..config.init_iters {
        gibbs_sample(&mut z, j_beta, h_beta, config.sites);
    }
    
    let mut count: usize = 0;
    let mut mean = 0.0;
    let mut variance = 0.0;
    let mut corr = vec![0.0; config.sites / 2];
    for i in 0..config.sample_iters {
        gibbs_sample(&mut z, j_beta, h_beta, config.sites);

        if i % config.sample_freq == 0 {
            // Compute the magnetization of the system and update the mean and variance.
            count += 1;
            let m = z.iter().sum::<i32>() as f64;
            let delta = m - mean;
            mean += delta / count as f64;
            variance += delta * (m - mean);

            // Compute the inter-site correlation and update its average value.
            let temp_corr = correlation(&z, config.sites);
            for i in 0..config.sites / 2 {
                let delta = temp_corr[i] - corr[i];
                corr[i] += delta / count as f64;
            }
        }
    }
    variance /= (count - 1) as f64;

    println!("observations = {count}, M = {mean}, std. dev. = {}", variance.sqrt());
    
    let mut corr_string = String::new();
    for (i, value) in corr.iter().enumerate() {
        corr_string.push_str(&format!("{i:>6} {value:>15.7E}\n"));
    }

    let correlation_file = config.correlation_file.unwrap_or(PathBuf::from_str("corr.dat").unwrap());
    fs::write(correlation_file, corr_string)?;

    Ok(())
}

/// Update the state of the system using Gibbs sampling.
fn gibbs_sample(z: &mut Vec<i32>, j_beta: f64, h_beta: f64, sites: usize) {
    let p = prob(z[sites - 1], z[1], j_beta, h_beta);
    let u: f64 = rand::thread_rng().sample(Standard);
    z[0] = if u < p { 1 } else { -1 };

    for i in 1..sites - 1 {
        let p = prob(z[i - 1], z[i + 1], j_beta, h_beta);
        let u: f64 = rand::thread_rng().sample(Standard);
        z[i] = if u < p { 1 } else { -1 };
    }

    let p = prob(z[sites - 2], z[0], j_beta, h_beta);
    let u: f64 = rand::thread_rng().sample(Standard);
    z[sites - 1] = if u < p { 1 } else { -1 };
}

/// The probability that a spin is in the spin up state, conditional on the state of its neighbors.
fn prob(z_m: i32, z_p: i32, j_beta: f64, h_beta: f64) -> f64 {
    let arg = z_m + z_p;
    let arg = j_beta * arg as f64 + h_beta;
    arg.exp() / (2.0 * arg.cosh())
}

/// Compute the expectation value of the inter-site correlation z[i] * z[j] as a function of inter-site
/// distance |i - j| for distances up to half the size of the periodic box.
fn correlation(z: &Vec<i32>, sites: usize) -> Vec<f64> {
    let mut corr: Vec<f64> = Vec::with_capacity(sites / 2);
    let m2 = z.iter().map(|z| (z * z) as f64).sum::<f64>();
    corr.push(m2 / sites as f64);
    for i in 1..sites / 2 {
        let m2 = (0..sites).map(|j| (z[j] * z[(j + i) % sites]) as f64).sum::<f64>();
        corr.push(m2 / sites as f64);
    }
    corr
}

#[derive(Deserialize)]
struct Config {
    sites: usize,
    init_iters: u32,
    sample_iters: u32,
    sample_freq: u32,
    coupling: f64,
    field: f64,
    beta: f64,
    correlation_file: Option<PathBuf>
}

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Arguments {
    config_file: PathBuf
}