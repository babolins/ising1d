use clap::Parser;
use rand::{distributions::Standard, Rng};
use std::{error::Error, fs, path::PathBuf, str::FromStr};
use serde::Deserialize;

fn main() -> Result<(), Box<dyn Error>> {
    let args = Arguments::parse();
    let config_file = fs::read_to_string(args.config_file)?;
    let config: Config = toml::from_str(&config_file)?;

    let j = config.coupling;
    let h = config.field;
    
    // Spin state of the system, the value 1 = spin up, and the value -1 = spin down.
    let mut sites = vec![1; config.sites];
    assert!(sites.len() >= 3);
    
    for _ in 0..config.init_iters {
        gibbs_sample(&mut sites[..], j, h);
    }
    
    let mut count: usize = 0;
    let mut mean = 0.0;
    let mut variance = 0.0;
    let mut corr = vec![0.0; config.sites / 2];
    for i in 0..config.sample_iters {
        gibbs_sample(&mut sites[..], j, h);

        if i % config.sample_freq == 0 {
            // Compute the magnetization of the system and update the mean and variance.
            count += 1;
            let m = sites.iter().sum::<i32>() as f64;
            let delta = m - mean;
            mean += delta / count as f64;
            variance += delta * (m - mean);

            // Compute the inter-site correlation and update its average value.
            let temp_corr = correlation(&sites);
            for i in 0..corr.len() {
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
fn gibbs_sample(sites: &mut [i32], j: f64, h: f64) {
    let p = prob(sites[sites.len() - 1], sites[1], j, h);
    let u: f64 = rand::thread_rng().sample(Standard);
    sites[0] = if u < p { 1 } else { -1 };

    for i in 1..sites.len() - 1 {
        let p = prob(sites[i - 1], sites[i + 1], j, h);
        let u: f64 = rand::thread_rng().sample(Standard);
        sites[i] = if u < p { 1 } else { -1 };
    }

    let p = prob(sites[sites.len() - 2], sites[0], j, h);
    let u: f64 = rand::thread_rng().sample(Standard);
    sites[sites.len() - 1] = if u < p { 1 } else { -1 };
}

/// The probability that a spin is in the spin up state, conditional on the state of its neighbors.
fn prob(z_m: i32, z_p: i32, j: f64, h: f64) -> f64 {
    let arg = z_m + z_p;
    let arg = j * arg as f64 + h;
    arg.exp() / (2.0 * arg.cosh())
}

/// Compute the expectation value of the inter-site correlation z[i] * z[j] as a function of inter-site
/// distance |i - j| for distances up to half the size of the periodic box.
fn correlation(sites: &[i32]) -> Vec<f64> {
    let mut corr: Vec<f64> = Vec::with_capacity(sites.len() / 2);
    let m2 = sites.iter().map(|z| (z * z) as f64).sum::<f64>();
    corr.push(m2 / sites.len() as f64);
    for i in 1..sites.len() / 2 {
        let m2 = (0..sites.len()).map(|j| (sites[j] * sites[(j + i) % sites.len()]) as f64).sum::<f64>();
        corr.push(m2 / sites.len() as f64);
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
    correlation_file: Option<PathBuf>
}

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Arguments {
    config_file: PathBuf
}