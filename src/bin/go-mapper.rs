use anyhow::{Context, Result};
use serde::Serialize;
use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
    str::FromStr,
};

#[derive(Debug)]
enum Aspect {
    C,
    F,
    P,
}

impl FromStr for Aspect {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "C" => Ok(Aspect::C),
            "F" => Ok(Aspect::F),
            "P" => Ok(Aspect::P),
            _ => Err(format!("Invalid aspect: {}", s)),
        }
    }
}

#[derive(Debug)]
struct GeneOntologySlimLine {
    orf: String,
    gene: String,
    aspect: Aspect,
    term: String,
}

impl GeneOntologySlimLine {
    fn new(orf: &str, gene: &str, aspect: &str, term: &str) -> Self {
        Self {
            orf: orf.to_string(),
            gene: gene.to_string(),
            aspect: Aspect::from_str(aspect).unwrap(),
            term: term.to_string(),
        }
    }
}

#[derive(Debug, Serialize)]
struct Row {
    #[serde(rename(serialize = "Systematic name"))]
    orf: String,
    #[serde(rename(serialize = "Standard name"))]
    gene: String,
    #[serde(rename(serialize = "GO-P (Process)"))]
    go_p: String,
    #[serde(rename(serialize = "GO-F (Function)"))]
    go_f: String,
    #[serde(rename(serialize = "GO-C (Component)"))]
    go_c: String,
    #[serde(rename(serialize = "mRNA abundance (Csardi)"))]
    mrna_abundance_csardi: f64,
    #[serde(rename(serialize = "Ribosome footprint density (Csardi)"))]
    ribosome_footprint_density_csardi: f64,
    #[serde(rename(serialize = "Translational efficiency (Csardi)"))]
    translational_efficiency_csardi: f64,
    #[serde(rename(serialize = "mRNA abundance (Weinberg)"))]
    mrna_abundance_weinberg: f64,
    #[serde(rename(serialize = "Ribosome footprint density (Weinberg)"))]
    ribosome_footprint_density_weinberg: f64,
    #[serde(rename(serialize = "Translational efficiency (Weinberg)"))]
    translational_efficiency_weinberg: f64,
}

impl Row {
    fn new(
        orf: &str,
        gene: &str,
        go_p: &str,
        go_f: &str,
        go_c: &str,
        csardi: Option<&TranslationalEfficiency>,
        weinberg: Option<&TranslationalEfficiency>,
    ) -> Self {
        let (mrna_csardi, rd_csardi, te_csardi) = if let Some(csardi) = csardi {
            (
                csardi.mrna_abundance,
                csardi.ribosome_footprint_density,
                csardi.translational_efficiency,
            )
        } else {
            (f64::NAN, f64::NAN, f64::NAN)
        };
        let (mrna_weinberg, rd_weinberg, te_weinberg) = if let Some(weinberg) = weinberg {
            (
                weinberg.mrna_abundance,
                weinberg.ribosome_footprint_density,
                weinberg.translational_efficiency,
            )
        } else {
            (f64::NAN, f64::NAN, f64::NAN)
        };

        Self {
            orf: orf.to_string(),
            gene: gene.to_string(),
            go_p: go_p.to_string(),
            go_f: go_f.to_string(),
            go_c: go_c.to_string(),
            mrna_abundance_csardi: mrna_csardi,
            ribosome_footprint_density_csardi: rd_csardi,
            translational_efficiency_csardi: te_csardi,
            mrna_abundance_weinberg: mrna_weinberg,
            ribosome_footprint_density_weinberg: rd_weinberg,
            translational_efficiency_weinberg: te_weinberg,
        }
    }
}

#[derive(Debug)]
struct TranslationalEfficiency {
    mrna_abundance: f64,
    ribosome_footprint_density: f64,
    translational_efficiency: f64,
}

impl TranslationalEfficiency {
    fn new(
        mrna_abundance: f64,
        ribosome_footprint_density: f64,
        translational_efficiency: f64,
    ) -> Self {
        Self {
            mrna_abundance,
            ribosome_footprint_density,
            translational_efficiency,
        }
    }
}

fn read_translational_efficiency_csv(
    path: &Path,
) -> Result<HashMap<String, TranslationalEfficiency>> {
    let file: File = File::open(path).unwrap();
    let reader: BufReader<File> = BufReader::new(file);
    let mut result = HashMap::new();

    for line in reader.lines().skip(1) {
        if let Ok(line) = line {
            let fields: Vec<&str> = line.split(',').collect();
            result.insert(
                fields[0].to_string(),
                TranslationalEfficiency::new(
                    fields[1]
                        .parse::<f64>()
                        .with_context(|| format!("Failed to parse {}", fields[1]))?,
                    fields[2]
                        .parse::<f64>()
                        .with_context(|| format!("Failed to parse {}", fields[2]))?,
                    fields[3]
                        .parse::<f64>()
                        .with_context(|| format!("Failed to parse {}", fields[3]))?,
                ),
            );
        }
    }

    return Ok(result);
}

fn read_gene_ontology_slim(path: &Path) -> Vec<GeneOntologySlimLine> {
    let file: File = File::open(path).unwrap();
    let reader: BufReader<File> = BufReader::new(file);
    let mut result = Vec::new();

    for line in reader.lines() {
        if let Ok(line) = line {
            let fields: Vec<&str> = line.split('\t').collect();
            result.push(GeneOntologySlimLine::new(
                fields[0], fields[1], fields[3], fields[4],
            ));
        }
    }

    result
}

fn combine_terms(
    slims: &Vec<GeneOntologySlimLine>,
    csardi: &HashMap<String, TranslationalEfficiency>,
    weinberg: &HashMap<String, TranslationalEfficiency>,
) -> Vec<Row> {
    let mut c_map: HashMap<String, Vec<String>> = HashMap::new();
    let mut f_map: HashMap<String, Vec<String>> = HashMap::new();
    let mut p_map: HashMap<String, Vec<String>> = HashMap::new();

    for slim in slims.iter() {
        match slim.aspect {
            Aspect::C => {
                c_map
                    .entry(slim.orf.clone())
                    .or_insert(Vec::new())
                    .push(slim.term.clone());
            }
            Aspect::F => {
                f_map
                    .entry(slim.orf.clone())
                    .or_insert(Vec::new())
                    .push(slim.term.clone());
            }
            Aspect::P => {
                p_map
                    .entry(slim.orf.clone())
                    .or_insert(Vec::new())
                    .push(slim.term.clone());
            }
        }
    }

    let mut results = Vec::new();
    let mut processed = HashSet::new();

    for slim in slims.iter() {
        if !processed.contains(&slim.orf) {
            processed.insert(slim.orf.clone());

            let go_p = if let Some(p) = p_map.get(&slim.orf) {
                p.join("; ")
            } else {
                "".to_string()
            };
            let go_f = if let Some(f) = f_map.get(&slim.orf) {
                f.join("; ")
            } else {
                "".to_string()
            };
            let go_c = if let Some(c) = c_map.get(&slim.orf) {
                c.join("; ")
            } else {
                "".to_string()
            };

            results.push(Row::new(
                &slim.orf,
                &slim.gene,
                &go_p,
                &go_f,
                &go_c,
                csardi.get(&slim.orf),
                weinberg.get(&slim.orf),
            ));
        }
    }

    results
}

fn main() {
    let slim_lines = read_gene_ontology_slim(Path::new("../data/go_slim_mapping.tab"));
    let csardi =
        read_translational_efficiency_csv(Path::new("../data/translational-efficiency-csardi.csv"))
            .unwrap();
    let weinberg = read_translational_efficiency_csv(Path::new(
        "../data/translational-efficiency-weinberg.csv",
    ))
    .unwrap();

    let rows = combine_terms(&slim_lines, &csardi, &weinberg);
    let mut writer = csv::Writer::from_path(Path::new("go_mapper.csv")).unwrap();

    for row in rows {
        writer.serialize(row).unwrap();
    }
}
