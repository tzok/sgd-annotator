use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

use flate2::read::GzDecoder;
use regex::Regex;

use crate::translator::GenomicRange;

#[derive(Debug)]
pub struct Fasta {
    header: String,
    sequence: String,
}

#[derive(Debug, Eq, PartialEq)]
pub enum FastaType {
    Chromosome,
    Gene,
    UTR,
}

impl Fasta {
    fn new(header: &str, sequence: &str) -> Self {
        Self {
            header: header.to_string(),
            sequence: sequence
                .chars()
                .map(|c| c.to_ascii_uppercase())
                .map(|c| match c {
                    'T' => 'U',
                    _ => c,
                })
                .collect(),
        }
    }

    fn reversed(&self) -> Self {
        Self {
            header: self.header.to_string(),
            sequence: self
                .sequence
                .chars()
                .rev()
                .map(|c| match c {
                    'A' => 'U',
                    'U' => 'A',
                    'C' => 'G',
                    'G' => 'C',
                    _ => c,
                })
                .collect(),
        }
    }

    pub fn fasta_type(&self) -> FastaType {
        if self.header.starts_with(">sacCer3") {
            return FastaType::UTR;
        } else if self.header.starts_with(">tpg") || self.header.starts_with(">ref") {
            return FastaType::Chromosome;
        }
        FastaType::Gene
    }

    pub fn genomic_range(&self) -> GenomicRange {
        match self.fasta_type() {
            FastaType::Chromosome => self.genomic_range_for_chromosome(),
            FastaType::Gene => self.genomic_range_for_gene(),
            FastaType::UTR => self.genomic_range_for_utr(),
        }
    }

    fn genomic_range_for_chromosome(&self) -> GenomicRange {
        let regex = Regex::new(r"\[chromosome=([IVX]+|mitochondrion)\]").unwrap();
        let capture = regex.captures(&self.header).unwrap();
        let chromosome = capture[1].to_string();
        GenomicRange {
            chromosome: if chromosome == "mitochondrion" {
                "Mito".to_string()
            } else {
                chromosome
            },
            range: 1..self.sequence.len() + 1,
        }
    }

    fn genomic_range_for_gene(&self) -> GenomicRange {
        let regex = Regex::new(r"Chr ([IVX]+|Mito) from (\d+)-(\d+)").unwrap();
        let captures = regex.captures(&self.header).unwrap();
        let chromosome = captures[1].to_string();
        let from = captures[2].parse().unwrap();
        let to = captures[3].parse().unwrap();

        if from < to {
            GenomicRange {
                chromosome: chromosome,
                range: from..to,
            }
        } else {
            GenomicRange {
                chromosome: chromosome,
                range: to..from,
            }
        }
    }

    fn genomic_range_for_utr(&self) -> GenomicRange {
        let regex = Regex::new(r"range=chr([IVX]+):(\d+)-(\d+)").unwrap();
        let captures = regex.captures(&self.header).unwrap();
        let chromosome = captures[1].to_string();
        let from = captures[2].parse().unwrap();
        let to = captures[3].parse().unwrap();

        if from < to {
            GenomicRange {
                chromosome: chromosome,
                range: from..to,
            }
        } else {
            GenomicRange {
                chromosome: chromosome,
                range: to..from,
            }
        }
    }

    pub fn coding_ranges(&self) -> Option<Vec<GenomicRange>> {
        match self.fasta_type() {
            FastaType::Chromosome => None,
            FastaType::Gene => self.coding_ranges_for_gene(),
            FastaType::UTR => None,
        }
    }

    fn coding_ranges_for_gene(&self) -> Option<Vec<GenomicRange>> {
        let main_regex = Regex::new(r"Chr ([IVX]+|Mito) from ((\d+-\d+,)+)").unwrap();
        let main_match = main_regex.captures(&self.header)?;
        let chromosome = &main_match[1];
        let regex = Regex::new(r"(\d+)-(\d+)").unwrap();
        let mut result = Vec::new();

        for captures in regex.captures_iter(&main_match[2]) {
            let from: usize = captures[1].parse().unwrap();
            let to: usize = captures[2].parse().unwrap();

            if from < to {
                result.push(GenomicRange {
                    chromosome: chromosome.to_string(),
                    range: from..to,
                });
            } else {
                result.push(GenomicRange {
                    chromosome: chromosome.to_string(),
                    range: to..from,
                });
            }
        }

        Some(result)
    }

    pub fn noncoding_ranges(&self) -> Option<Vec<GenomicRange>> {
        match self.fasta_type() {
            FastaType::Chromosome => None,
            FastaType::Gene => self.noncoding_ranges_for_gene(),
            FastaType::UTR => None,
        }
    }

    fn noncoding_ranges_for_gene(&self) -> Option<Vec<GenomicRange>> {
        if let Some(coding) = self.coding_ranges() {
            if coding.len() > 1 {
                let mut result = Vec::new();
                for i in 1..coding.len() {
                    result.push(GenomicRange {
                        chromosome: coding[i].chromosome.clone(),
                        range: coding[i - 1].range.end..coding[i].range.start,
                    })
                }
                return Some(result);
            }
        }
        None
    }

    pub fn systematic_name(&self) -> String {
        match self.fasta_type() {
            FastaType::Chromosome => self.systematic_name_for_chromosome(),
            FastaType::Gene => self.systematic_name_for_gene(),
            FastaType::UTR => self.systematic_name_for_utr(),
        }
    }
    fn systematic_name_for_chromosome(&self) -> String {
        format!("chr{}", self.genomic_range_for_chromosome().chromosome)
    }

    fn systematic_name_for_gene(&self) -> String {
        self.header[1..]
            .split_whitespace()
            .nth(0)
            .map(|s| s.to_string())
            .unwrap()
    }

    fn systematic_name_for_utr(&self) -> String {
        self.header
            .split('_')
            .nth(4)
            .map(|s| s.to_string())
            .unwrap()
    }

    pub fn standard_name(&self) -> String {
        match self.fasta_type() {
            FastaType::Chromosome => self.standard_name_for_chromosome(),
            FastaType::Gene => self.standard_name_for_gene(),
            FastaType::UTR => self.standard_name_for_utr(),
        }
    }

    fn standard_name_for_chromosome(&self) -> String {
        self.systematic_name_for_chromosome()
    }

    fn standard_name_for_gene(&self) -> String {
        self.header[1..]
            .split_whitespace()
            .nth(1)
            .map(|s| s.to_string())
            .unwrap()
    }

    fn standard_name_for_utr(&self) -> String {
        self.systematic_name_for_utr()
    }
}

pub fn load_fasta_gz(path: &Path) -> HashMap<String, Fasta> {
    let file = File::open(path).unwrap();
    let decoder = GzDecoder::new(file);
    let reader = BufReader::new(decoder);
    let mut header = String::new();
    let mut sequence = String::new();
    let mut result = HashMap::new();

    for line in reader.lines() {
        if let Ok(content) = line {
            if content.starts_with(">") {
                if !header.is_empty() {
                    let fasta = Fasta::new(&header, &sequence);
                    result.insert(fasta.systematic_name(), fasta);
                    header.clear();
                    sequence.clear();
                }
                header += content.trim();
            } else {
                sequence += content.trim();
            }
        }
    }

    if !header.is_empty() {
        let fasta = Fasta::new(&header, &sequence);
        result.insert(fasta.systematic_name(), fasta);
    }

    result
}

pub fn load_utr_fasta_gz(path: &Path) -> HashMap<String, Fasta> {
    load_fasta_gz(path)
        .into_iter()
        .map(|(_, fasta)| (fasta.systematic_name_for_utr(), fasta))
        .collect()
}
