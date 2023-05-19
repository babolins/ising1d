use rand::Rng;
use rand::distributions::Standard;

fn main() {
    let sites: usize = 1024;
    let init_iters: usize = 10000;
    let sample_iters: usize = 100000;
    let sample_freq: usize = 100;
    let j = 1.0;
    let h = 0.0;
    let beta = 1.0;
    let j_beta = j / beta;
    let h_beta = h / beta;
    
    let mut z = vec![1; sites];
    for _ in 0..init_iters {
        gibbs_sample(&mut z, j_beta, h_beta, sites);
    }
    
    let mut count: usize = 0;
    let mut mean = 0.0;
    let mut variance = 0.0;
    for i in 0..sample_iters {
        gibbs_sample(&mut z, j_beta, h_beta, sites);

        if i % sample_freq == 0 {
            count += 1;
            let m = z.iter().sum::<i32>() as f64;
            let delta = m - mean;
            mean += delta / count as f64;
            variance += delta * (m - mean);
        }
    }
    variance /= (count - 1) as f64;

    println!("observations = {count}, M = {mean}, std. dev. = {}", variance.sqrt());
}

fn gibbs_sample(z: &mut Vec<i32>, j_beta: f64, h_beta: f64, sites: usize) {
    let arg = z[sites - 1] + z[1];
    let arg = -j_beta * arg as f64 + h_beta;
    let p = prob(arg);
    let u: f64 = rand::thread_rng().sample(Standard);
    z[0] = if u < p { 1 } else { -1 };


    for i in 1..sites - 1 {
        let arg = z[i - 1] + z[i + 1];
        let arg = -j_beta * arg as f64 + h_beta;
        let p = prob(arg);
        let u: f64 = rand::thread_rng().sample(Standard);
        z[i] = if u < p { 1 } else { -1 };
    }

    let arg = z[sites - 2] + z[0];
    let arg = -j_beta * arg as f64 + h_beta;
    let p = prob(arg);
    let u: f64 = rand::thread_rng().sample(Standard);
    z[sites - 1] = if u < p { 1 } else { -1 };
}

fn prob(arg: f64) -> f64 {
    arg.exp() / (2.0 * arg.cosh())
}