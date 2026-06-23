use std::path::Path;
use std::str::FromStr;

use anyhow::{anyhow, Result};
use bigtools::BigWigRead;
use clap::ValueEnum;
use log::{debug, warn};

use crate::translator::{Translator, YeastChromosome};

pub struct BigWigProfile {
    pub name: String,
    pub values: Vec<f32>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum MissingMode {
    /// Emit `0` for positions with no bigWig data (correct for Ribo-Seq count tracks).
    Zero,
    /// Emit an empty cell for missing values (preserves gap vs. zero distinction).
    Empty,
    /// Emit the literal string `nan`.
    Nan,
    /// Emit `?` to mark missing/unknown values.
    Question,
}

pub fn format_value(value: f32) -> String {
    if value == value.trunc() {
        format!("{}", value as i64)
    } else {
        format!("{}", value)
    }
}

pub fn format_bigwig_cell(value: f32, mode: MissingMode) -> String {
    if value.is_nan() {
        match mode {
            MissingMode::Zero => "0".to_string(),
            MissingMode::Empty => String::new(),
            MissingMode::Nan => "nan".to_string(),
            MissingMode::Question => "?".to_string(),
        }
    } else {
        format_value(value)
    }
}

fn parse_chrom(name: &str) -> Option<YeastChromosome> {
    let trimmed = name.strip_prefix("chr").unwrap_or(name);
    if trimmed.eq_ignore_ascii_case("M") || trimmed.eq_ignore_ascii_case("Mito") {
        return YeastChromosome::from_str("Mito").ok();
    }
    YeastChromosome::from_str(trimmed).ok()
}

pub fn load(path: &Path, name: String, translator: &Translator, genome_len: usize) -> Result<BigWigProfile> {
    let mut bigwig = BigWigRead::open_file(path)
        .map_err(|e| anyhow!("failed to open bigWig {}: {}", path.display(), e))?;

    let chroms: Vec<(String, u32)> = bigwig
        .chroms()
        .iter()
        .map(|c| (c.name.clone(), c.length))
        .collect();

    let mut values = vec![f32::NAN; genome_len];
    let mut filled = 0usize;

    for (chrom_name, length) in &chroms {
        let yeast_chrom = match parse_chrom(chrom_name) {
            Some(c) => c,
            None => {
                debug!("bigWig {}: skipping unknown contig {}", path.display(), chrom_name);
                continue;
            }
        };

        let base = match translator.chromosome_offset(&yeast_chrom) {
            Some(b) => b,
            None => {
                warn!(
                    "bigWig {}: chromosome {} not present in genome reference; skipping",
                    path.display(),
                    chrom_name
                );
                continue;
            }
        };

        let interval = match bigwig.get_interval(chrom_name, 0, *length) {
            Ok(it) => it,
            Err(e) => {
                warn!(
                    "bigWig {}: failed to read interval for {} ({}): {}",
                    path.display(),
                    chrom_name,
                    length,
                    e
                );
                continue;
            }
        };

        for value in interval {
            let value = match value {
                Ok(v) => v,
                Err(e) => {
                    warn!(
                        "bigWig {}: error reading value on {}: {}",
                        path.display(),
                        chrom_name,
                        e
                    );
                    continue;
                }
            };
            let start = value.start as usize;
            let end = value.end as usize;
            if start >= end {
                continue;
            }
            let lo = base + start;
            let hi = base + end;
            if hi > genome_len {
                warn!(
                    "bigWig {}: value on {} exceeds genome length ({} > {}); clamping",
                    path.display(),
                    chrom_name,
                    hi,
                    genome_len
                );
            }
            let clamp_hi = hi.min(genome_len);
            if lo >= clamp_hi {
                continue;
            }
            for v in &mut values[lo..clamp_hi] {
                *v = value.value;
            }
            filled += clamp_hi - lo;
        }
    }

    debug!(
        "bigWig {}: filled {} of {} positions",
        path.display(),
        filled,
        genome_len
    );

    Ok(BigWigProfile { name, values })
}

pub fn aggregate_dir(
    dir: &Path,
    name: String,
    translator: &Translator,
    genome_len: usize,
) -> Result<BigWigProfile> {
    use std::fs;

    let mut paths: Vec<std::path::PathBuf> = match fs::read_dir(dir) {
        Ok(rd) => rd
            .filter_map(|e| e.ok())
            .map(|e| e.path())
            .filter(|p| p.extension().and_then(|s| s.to_str()) == Some("bw"))
            .collect(),
        Err(e) => return Err(anyhow!("failed to read bigWig dir {}: {}", dir.display(), e)),
    };
    paths.sort();

    if paths.is_empty() {
        return Err(anyhow!(
            "no .bw files found in directory {}",
            dir.display()
        ));
    }

    let mut values: Vec<f32> = vec![f32::NAN; genome_len];
    let mut files_summed = 0usize;

    for path in &paths {
        let profile = load(path, String::new(), translator, genome_len)?;
        for (acc, new) in values.iter_mut().zip(profile.values.iter()) {
            *acc = if acc.is_nan() {
                *new
            } else if new.is_nan() {
                *acc
            } else {
                *acc + *new
            };
        }
        files_summed += 1;
    }

    let filled = values.iter().filter(|v| !v.is_nan()).count();
    log::info!(
        "bigWig dir {}: summed {} files, {} of {} positions non-NaN",
        dir.display(),
        files_summed,
        filled,
        genome_len
    );

    Ok(BigWigProfile { name, values })
}

#[cfg(test)]
mod tests {
    use super::{parse_chrom, MissingMode};
    use crate::translator::YeastChromosome;
    use bigtools::{beddata::BedParserStreamingIterator, BigWigRead, BigWigWrite, Value};
    use std::collections::HashMap;

    fn write_test_bigwig(path: &std::path::Path, intervals: Vec<(u32, u32, f32)>) {
        let mut chrom_sizes: HashMap<String, u32> = HashMap::new();
        chrom_sizes.insert("chrI".to_string(), 230218);
        let write = BigWigWrite::create_file(path, chrom_sizes).unwrap();
        let values = intervals
            .into_iter()
            .map(|(start, end, value)| {
                (
                    "chrI".to_string(),
                    Value { start, end, value },
                )
            })
            .collect::<Vec<_>>();
        let data = BedParserStreamingIterator::wrap_infallible_iter(values.into_iter(), true);
        let runtime = tokio::runtime::Builder::new_multi_thread()
            .enable_all()
            .build()
            .unwrap();
        write.write(data, runtime).unwrap();
    }

    #[test]
    fn parse_chrom_names() {
        assert_eq!(parse_chrom("chrI"), Some(YeastChromosome::I));
        assert_eq!(parse_chrom("chrXVI"), Some(YeastChromosome::XVI));
        assert_eq!(parse_chrom("chrM"), Some(YeastChromosome::Mito));
        assert_eq!(parse_chrom("chrMito"), Some(YeastChromosome::Mito));
        assert_eq!(parse_chrom("I"), Some(YeastChromosome::I));
        assert_eq!(parse_chrom("Mito"), Some(YeastChromosome::Mito));
        assert_eq!(parse_chrom("2-micron"), None);
        assert_eq!(parse_chrom("chrFoo"), None);
    }

    #[test]
    fn format_bigwig_cell_modes() {
        let nan = f32::NAN;
        assert_eq!(super::format_bigwig_cell(nan, MissingMode::Zero), "0");
        assert_eq!(super::format_bigwig_cell(nan, MissingMode::Empty), "");
        assert_eq!(super::format_bigwig_cell(nan, MissingMode::Nan), "nan");
        assert_eq!(super::format_bigwig_cell(nan, MissingMode::Question), "?");

        assert_eq!(super::format_bigwig_cell(42.0, MissingMode::Zero), "42");
        assert_eq!(super::format_bigwig_cell(42.0, MissingMode::Empty), "42");
        assert_eq!(super::format_bigwig_cell(7113.0, MissingMode::Nan), "7113");
        assert_eq!(super::format_bigwig_cell(1.5, MissingMode::Question), "1.5");
        assert_eq!(super::format_bigwig_cell(0.0, MissingMode::Zero), "0");
    }

    #[test]
    fn load_maps_values_to_genome_indices() {
        let (genome, _) = crate::load_genome_gz(std::path::Path::new("tests/genome.txt.gz"));
        let translator = crate::translator::Translator::new(&genome);
        let offset = translator
            .chromosome_offset(&YeastChromosome::I)
            .expect("chr I present in reference");

        let tmp = tempfile::NamedTempFile::new().unwrap();
        write_test_bigwig(tmp.path(), vec![(100, 110, 42.0)]);

        let profile = super::load(
            tmp.path(),
            "test".to_string(),
            &translator,
            genome.len(),
        )
        .unwrap();

        assert_eq!(profile.name, "test");
        assert_eq!(profile.values.len(), genome.len());
        for i in 100..110 {
            assert_eq!(profile.values[offset + i], 42.0, "position {}", i);
        }
        assert!(profile.values[offset + 99].is_nan());
        assert!(profile.values[offset + 110].is_nan());
        let _ = BigWigRead::open_file(tmp.path());
    }

    #[test]
    fn aggregate_dir_sums_and_preserves_all_nan() {
        let (genome, _) = crate::load_genome_gz(std::path::Path::new("tests/genome.txt.gz"));
        let translator = crate::translator::Translator::new(&genome);
        let offset = translator
            .chromosome_offset(&YeastChromosome::I)
            .expect("chr I present in reference");

        let dir = tempfile::TempDir::new().unwrap();
        // File A: [100, 110) value 2.0
        write_test_bigwig(
            &dir.path().join("a.bw"),
            vec![(100, 110, 2.0)],
        );
        // File B: [105, 115) value 3.0
        write_test_bigwig(
            &dir.path().join("b.bw"),
            vec![(105, 115, 3.0)],
        );

        let profile = super::aggregate_dir(
            dir.path(),
            "agg".to_string(),
            &translator,
            genome.len(),
        )
        .unwrap();

        assert_eq!(profile.name, "agg");
        assert_eq!(profile.values.len(), genome.len());

        // only A
        for i in 100..105 {
            assert_eq!(profile.values[offset + i], 2.0, "position {} (A only)", i);
        }
        // A + B
        for i in 105..110 {
            assert_eq!(profile.values[offset + i], 5.0, "position {} (A+B)", i);
        }
        // only B
        for i in 110..115 {
            assert_eq!(profile.values[offset + i], 3.0, "position {} (B only)", i);
        }
        // neither
        assert!(profile.values[offset + 99].is_nan(), "all-gap position should be NaN");
        assert!(profile.values[offset + 200].is_nan(), "all-gap position should be NaN");
    }

    #[test]
    fn aggregate_dir_errors_on_empty_directory() {
        let (genome, _) = crate::load_genome_gz(std::path::Path::new("tests/genome.txt.gz"));
        let translator = crate::translator::Translator::new(&genome);
        let dir = tempfile::TempDir::new().unwrap();

        let result = super::aggregate_dir(
            dir.path(),
            "agg".to_string(),
            &translator,
            genome.len(),
        );
        let err = match result {
            Ok(_) => panic!("expected error for empty directory"),
            Err(e) => e,
        };

        assert!(
            err.to_string().contains("no .bw files"),
            "expected 'no .bw files' in error, got: {}",
            err
        );
    }
}
