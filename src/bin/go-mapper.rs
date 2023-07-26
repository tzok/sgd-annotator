use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
    str::FromStr,
};

use serde::Serialize;

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
}

impl Row {
    fn new(orf: &str, gene: &str, go_p: &str, go_f: &str, go_c: &str) -> Self {
        Self {
            orf: orf.to_string(),
            gene: gene.to_string(),
            go_p: go_p.to_string(),
            go_f: go_f.to_string(),
            go_c: go_c.to_string(),
        }
    }
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

fn combine_terms(slims: &Vec<GeneOntologySlimLine>) -> Vec<Row> {
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
            results.push(Row::new(&slim.orf, &slim.gene, &go_p, &go_f, &go_c));
        }
    }

    results
}

fn main() {
    let slim_lines = read_gene_ontology_slim(Path::new("../data/go_slim_mapping.tab"));
    let rows = combine_terms(&slim_lines);
    let mut writer = csv::Writer::from_path(Path::new("go_mapper.csv")).unwrap();

    for row in rows {
        writer.serialize(row).unwrap();
    }
}
