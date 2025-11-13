pub mod data;
pub mod fasta;
pub mod translator;
pub mod tests;

use flate2::read::GzDecoder;
use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

pub fn load_genome_gz(path: &Path) -> (String, Vec<f32>) {
    let file = File::open(path).unwrap();
    let decoder = GzDecoder::new(file);
    let reader = BufReader::new(decoder);
    let mut sequence = String::new();
    let mut profile: Vec<f32> = Vec::new();

    for (i, line) in reader.lines().enumerate() {
        if i == 0 {
            continue;
        }

        if let Ok(content) = line {
            let split: Vec<&str> = content.split_whitespace().collect();

            sequence += &split[1].to_ascii_uppercase().replace("T", "U");
            let value = if split.len() > 27 {
                split[27].parse().unwrap_or(f32::NAN)
            } else {
                f32::NAN
            };
            profile.push(value);
        }
    }

    (sequence, profile)
}
