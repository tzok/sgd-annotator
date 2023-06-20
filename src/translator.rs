use std::{collections::HashMap, fmt::Display, ops::Range, path::PathBuf, str::FromStr};

use crate::fasta::{load_fasta_gz, Fasta};
use glob::glob;
use log::debug;
use rayon::prelude::{IntoParallelRefIterator, ParallelExtend, ParallelIterator};

#[derive(Clone, Debug, Hash, PartialEq, Eq)]
pub enum YeastChromosome {
    I,
    II,
    III,
    IV,
    V,
    VI,
    VII,
    VIII,
    IX,
    X,
    XI,
    XII,
    XIII,
    XIV,
    XV,
    XVI,
    Mito,
}

impl FromStr for YeastChromosome {
    type Err = ();
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "I" => Ok(Self::I),
            "II" => Ok(Self::II),
            "III" => Ok(Self::III),
            "IV" => Ok(Self::IV),
            "V" => Ok(Self::V),
            "VI" => Ok(Self::VI),
            "VII" => Ok(Self::VII),
            "VIII" => Ok(Self::VIII),
            "IX" => Ok(Self::IX),
            "X" => Ok(Self::X),
            "XI" => Ok(Self::XI),
            "XII" => Ok(Self::XII),
            "XIII" => Ok(Self::XIII),
            "XIV" => Ok(Self::XIV),
            "XV" => Ok(Self::XV),
            "XVI" => Ok(Self::XVI),
            "Mito" => Ok(Self::Mito),
            _ => Err(()),
        }
    }
}

impl Display for YeastChromosome {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::I => write!(f, "I"),
            Self::II => write!(f, "II"),
            Self::III => write!(f, "III"),
            Self::IV => write!(f, "IV"),
            Self::V => write!(f, "V"),
            Self::VI => write!(f, "VI"),
            Self::VII => write!(f, "VII"),
            Self::VIII => write!(f, "VIII"),
            Self::IX => write!(f, "IX"),
            Self::X => write!(f, "X"),
            Self::XI => write!(f, "XI"),
            Self::XII => write!(f, "XII"),
            Self::XIII => write!(f, "XIII"),
            Self::XIV => write!(f, "XIV"),
            Self::XV => write!(f, "XV"),
            Self::XVI => write!(f, "XVI"),
            Self::Mito => write!(f, "Mito"),
        }
    }
}

#[derive(Debug, PartialEq, Eq)]
pub struct GenomicRange {
    pub chromosome: YeastChromosome,
    pub start: usize,
    pub end: usize,
}

pub struct Translator {
    mapping: HashMap<YeastChromosome, usize>,
}

impl Translator {
    pub fn new(genome: &str) -> Self {
        let paths: Vec<PathBuf> = glob("../data/chr??.fsa.gz")
            .unwrap()
            .filter_map(|path| {
                if path.is_ok() {
                    Some(path.unwrap())
                } else {
                    None
                }
            })
            .collect();

        let mut vec: Vec<(YeastChromosome, usize)> = Vec::new();
        vec.par_extend(paths.par_iter().filter_map(|path| {
            debug!("Loading chromosome {}", path.display());
            let chromosome = load_fasta_gz(&path);
            let fasta = chromosome.iter().map(|(_, fasta)| fasta).next()?;
            let index = find_sequence(genome, &fasta)?;
            Some((fasta.genomic_range().chromosome.clone(), index))
        }));

        let mapping: HashMap<YeastChromosome, usize> = vec.into_iter().collect();
        Self { mapping }
    }

    pub fn translate_fasta(&self, fasta: &Fasta) -> Option<(usize, usize)> {
        let range = fasta.genomic_range();
        let chromosome = &range.chromosome;
        Some((
            self.translate_nt(chromosome, range.start)?,
            self.translate_nt(chromosome, range.end)?,
        ))
    }

    pub fn translate_nt(&self, chromosome: &YeastChromosome, index: usize) -> Option<usize> {
        let base = self.mapping.get(chromosome)?;
        Some(base + index - 1)
    }
}

fn find_sequence(genome: &str, fasta: &Fasta) -> Option<usize> {
    if let Some(index) = genome.find(fasta.sequence().as_str()) {
        return Some(index);
    }
    debug!(
        "Failed to find sequence {} in genome",
        fasta.systematic_name()
    );
    None
}
