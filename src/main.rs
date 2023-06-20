use std::{
    cmp::{max, min},
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    path::Path,
};

use clap::Parser;
use fasta::{load_fasta_gz, load_utr_fasta_gz, Fasta};
use flate2::read::GzDecoder;

use translator::Translator;

mod fasta;
mod tests;
mod translator;

#[derive(Parser)]
#[command(version)]
struct Args {
    #[arg(short, long)]
    input: String,

    #[arg(short, long)]
    output: String,
}

fn load_genome_gz(path: &Path) -> String {
    let file = File::open(path).unwrap();
    let decoder = GzDecoder::new(file);
    let reader = BufReader::new(decoder);
    let mut result = String::new();

    for (i, line) in reader.lines().enumerate() {
        if i == 0 {
            continue;
        }

        if let Ok(content) = line {
            result += &content
                .split_whitespace()
                .nth(1)
                .map(|s| s.to_ascii_uppercase().replace("T", "U"))
                .unwrap();
        }
    }

    result
}

fn translate_all(
    all_genomic: &Vec<(&String, &Fasta)>,
    utr5p: &HashMap<String, Fasta>,
    utr3p: &HashMap<String, Fasta>,
    translator: &Translator,
) -> HashMap<String, (usize, usize)> {
    let mut ranges = HashMap::new();

    for (name, fasta) in all_genomic {
        if let Some((mut start, mut end)) =
            translator.translate_genomic_range(&fasta.genomic_range())
        {
            if let Some(utr) = utr5p.get(*name) {
                if let Some((utr_start, utr_end)) =
                    translator.translate_genomic_range(&utr.genomic_range())
                {
                    start = min(start, utr_start);
                    end = max(end, utr_end);
                }
            }

            if let Some(utr) = utr3p.get(*name) {
                if let Some((utr_start, utr_end)) =
                    translator.translate_genomic_range(&utr.genomic_range())
                {
                    start = min(start, utr_start);
                    end = max(end, utr_end);
                }
            }

            ranges.insert((*name).to_string(), (start, end));
        }
    }

    ranges
}

fn create_graph(
    all_genomic: &Vec<(&String, &Fasta)>,
    ranges: &HashMap<String, (usize, usize)>,
) -> HashMap<String, Vec<String>> {
    let names: Vec<&String> = all_genomic.iter().map(|(name, _)| *name).collect();
    let mut graph = HashMap::new();

    for i in 0..names.len() {
        let name_i = names.get(i).unwrap().to_string();

        if let Some((start_i, end_i)) = ranges.get(&name_i) {
            for j in (i + 1)..names.len() {
                let name_j = names.get(j).unwrap().to_string();

                if let Some((start_j, end_j)) = ranges.get(&name_j) {
                    if (start_i < start_j && start_j < end_i)
                        || (start_i < end_j && end_j < end_i)
                        || (start_j < start_i && start_i < end_j)
                        || (start_j < end_i && end_i < end_j)
                    {
                        if !graph.contains_key(&name_i) {
                            graph.insert(name_i.clone(), Vec::new());
                        }
                        graph.get_mut(&name_i).unwrap().push(name_j.clone());

                        if !graph.contains_key(&name_j) {
                            graph.insert(name_j.clone(), Vec::new());
                        }
                        graph.get_mut(&name_j).unwrap().push(name_i.clone());
                    }
                }
            }
        }
    }

    graph
}

fn determine_order(
    all_genomic: &Vec<(&String, &Fasta)>,
    graph: &HashMap<String, Vec<String>>,
) -> HashMap<String, usize> {
    let mut colors: HashMap<String, usize> = HashMap::new();

    for (current, _) in all_genomic.iter() {
        let mut available = [true, true, true, true, true, true, true, true, true, true];

        if let Some(adjacent) = graph.get(*current) {
            for next in adjacent {
                if let Some(color) = colors.get(next) {
                    available[*color] = false;
                }
            }
        }

        for i in 0..available.len() {
            if available[i] {
                colors.insert(current.to_string(), i);
                break;
            }
        }
    }

    colors
}

fn fill_annotations(
    genome: &str,
    all_genomic: &Vec<(&String, &Fasta)>,
    orf_genomic: &HashMap<String, Fasta>,
    rna_genomic: &HashMap<String, Fasta>,
    utr5p: &HashMap<String, Fasta>,
    utr3p: &HashMap<String, Fasta>,
    orf_coding: &HashMap<String, Fasta>,
    rna_coding: &HashMap<String, Fasta>,
    ranges: &HashMap<String, (usize, usize)>,
    orders: &HashMap<String, usize>,
    translator: &Translator,
) -> Vec<Vec<String>> {
    let max = orders.values().max().unwrap();
    let mut annotations: Vec<Vec<String>> = Vec::with_capacity(genome.len());
    for _ in 0..genome.len() {
        let mut v = Vec::with_capacity((*max + 1) * 4);
        for _ in 0..(*max + 1) * 4 {
            v.push(String::new());
        }
        annotations.push(v);
    }

    for (name, fasta) in all_genomic.iter() {
        let category = if orf_genomic.contains_key(*name) {
            "ORF"
        } else if rna_genomic.contains_key(*name) {
            "RNA"
        } else {
            "Other"
        };

        let order = orders.get(*name).unwrap();
        let (start, end) = ranges.get(*name).unwrap();
        for i in *start..=*end {
            annotations[i][*order * 4] = category.to_string();
            annotations[i][*order * 4 + 1] = "?".to_string();
            annotations[i][*order * 4 + 2] = fasta.systematic_name().to_string();
            annotations[i][*order * 4 + 3] = fasta.standard_name().to_string();
        }

        if category == "ORF" {
            if let Some(fasta) = utr5p.get(*name) {
                let range = fasta.genomic_range();
                let (start, end) = translator.translate_genomic_range(&range).unwrap();
                for i in start..=end {
                    annotations[i][*order * 4 + 1] = "UTR 5'".to_string();
                }
            }

            if let Some(fasta) = utr3p.get(*name) {
                let range = fasta.genomic_range();
                let (start, end) = translator.translate_genomic_range(&range).unwrap();
                for i in start..=end {
                    annotations[i][*order * 4 + 1] = "UTR 3'".to_string();
                }
            }

            if let Some(fasta) = orf_coding.get(*name) {
                for range in fasta.coding_ranges().unwrap().iter() {
                    let (start, end) = translator.translate_genomic_range(range).unwrap();
                    for i in start..=end {
                        annotations[i][*order * 4 + 1] = "Exon".to_string();
                    }
                }
                if fasta.noncoding_ranges().is_some() {
                    for range in fasta.noncoding_ranges().unwrap().iter() {
                        let (start, end) = translator.translate_genomic_range(range).unwrap();
                        for i in start..=end {
                            annotations[i][*order * 4 + 1] = "Intron".to_string();
                        }
                    }
                }
            }
        } else if category == "RNA" {
            if let Some(fasta) = rna_coding.get(*name) {
                for range in fasta.coding_ranges().unwrap().iter() {
                    let (start, end) = translator.translate_genomic_range(range).unwrap();
                    for i in start..=end {
                        annotations[i][*order * 4 + 1] = "Exon".to_string();
                    }
                }
                if fasta.noncoding_ranges().is_some() {
                    for range in fasta.noncoding_ranges().unwrap().iter() {
                        let (start, end) = translator.translate_genomic_range(range).unwrap();
                        for i in start..=end {
                            annotations[i][*order * 4 + 1] = "Intron".to_string();
                        }
                    }
                }
            }
        }
    }

    annotations
}

fn store_result(input: &Path, output: &Path, annotations: Vec<Vec<String>>) {
    let fin = File::open(input).unwrap();
    let fout = File::create(output).unwrap();
    let decoder = GzDecoder::new(fin);
    let reader = BufReader::new(decoder);
    let mut writer = BufWriter::new(fout);

    let count = annotations.get(0).unwrap().len();

    for (i, line) in reader.lines().enumerate() {
        if let Ok(line) = line {
            let _ = writer.write(line.as_bytes());

            if i == 0 {
                for j in 0..count / 4 {
                    let _ = writer.write(
                        format!(
                            "\tType {}\tSubtype {}\tSystematic name {}\tStandard name {}",
                            j + 1,
                            j + 1,
                            j + 1,
                            j + 1
                        )
                        .as_bytes(),
                    );
                }
                let _ = writer.write("\n".as_bytes());
            } else {
                for annotation in annotations.get(i - 1).unwrap().iter() {
                    let _ = writer.write(b"\t");
                    let _ = writer.write(annotation.as_bytes());
                }
                let _ = writer.write(b"\n");
            }
        }
    }
}

fn main() {
    let args = Args::parse();

    let genome = load_genome_gz(Path::new(&args.input));
    let translator = Translator::new(&genome);

    let orf_genomic = load_fasta_gz(Path::new("../data/orf_genomic.fasta.gz"));
    let rna_genomic = load_fasta_gz(Path::new("../data/rna_genomic.fasta.gz"));
    let other_genomic = load_fasta_gz(Path::new("../data/other_features_genomic.fasta.gz"));

    let all_genomic: Vec<(&String, &Fasta)> = orf_genomic
        .iter()
        .chain(rna_genomic.iter())
        .chain(other_genomic.iter())
        .filter(|(_, fasta)| {
            translator
                .translate_genomic_range(&fasta.genomic_range())
                .is_some()
        })
        .collect();

    let orf_coding = load_fasta_gz(Path::new("../data/orf_coding.fasta.gz"));
    let rna_coding = load_fasta_gz(Path::new("../data/rna_coding.fasta.gz"));

    let utr5p = load_utr_fasta_gz(Path::new("../data/SGD_all_ORFs_5prime_UTRs.fsa.gz"));
    let utr3p = load_utr_fasta_gz(Path::new("../data/SGD_all_ORFs_3prime_UTRs.fsa.gz"));

    let ranges = translate_all(&all_genomic, &utr5p, &utr3p, &translator);
    let graph = create_graph(&all_genomic, &ranges);
    let orders = determine_order(&all_genomic, &graph);

    let annotations = fill_annotations(
        &genome,
        &all_genomic,
        &orf_genomic,
        &rna_genomic,
        &utr5p,
        &utr3p,
        &orf_coding,
        &rna_coding,
        &ranges,
        &orders,
        &translator,
    );
    store_result(Path::new(&args.input), Path::new(&args.output), annotations);
}
