use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader},
    ops::Range,
    path::Path,
};

use clap::Parser;
use flate2::read::GzDecoder;
use rayon::prelude::*;
use regex::Regex;

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

fn load_genome(path: &Path) -> String {
    let file = File::open(path).unwrap();
    let reader = BufReader::new(file);
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

// fn find_sequence(
//     genome: &str,
//     genomic_1000: &HashMap<String, Fasta>,
//     genomic: &HashMap<String, Fasta>,
//     name: &str,
//     reverse: bool,
// ) -> Option<FileRange> {
//     let fasta_1000 = genomic_1000.get(name)?;
//     let fasta = genomic.get(name)?;

//     if reverse {
//         let index_1000 = genome.find(fasta_1000.reversed().sequence.as_str())?;
//         let index = genome[index_1000..].find(fasta.reversed().sequence.as_str())?;
//         Some(FileRange {
//             range: index_1000 + index..index_1000 + index + fasta.sequence.len(),
//         })
//     } else {
//         let index_1000 = genome.find(fasta_1000.sequence.as_str())?;
//         let index = genome[index_1000..].find(fasta.sequence.as_str())?;
//         Some(FileRange {
//             range: index_1000 + index..index_1000 + index + fasta.sequence.len(),
//         })
//     }
// }

// fn collect_found(
//     genome: &str,
//     genomic_1000: &HashMap<String, Fasta>,
//     genomic: &HashMap<String, Fasta>,
// ) -> HashMap<String, FileRange> {
//     genomic_1000
//         .par_iter()
//         .filter_map(|(name, _)| {
//             if let Some(range) = find_sequence(genome, genomic_1000, genomic, name, false) {
//                 return Some((name.to_string(), range));
//             }
//             if let Some(range) = find_sequence(genome, genomic_1000, genomic, name, true) {
//                 return Some((name.to_string(), range));
//             }
//             None
//         })
//         .collect()
// }

// fn extend_file_range_for_utr(
//     found: &mut HashMap<String, FileRange>,
//     orf: &HashMap<String, Fasta>,
//     utrs: &HashMap<String, Fasta>,
// ) {
//     for (name, file_range) in found.iter_mut() {
//         if let Some(gene) = orf.get(name) {
//             if let Some(gene_range) = gene.genomic_range() {
//                 if let Some(utr) = utrs.get(name) {
//                     if let Some(utr_range) = utr.genomic_range_for_utr() {
//                         assert!(gene_range.chromosome == utr_range.chromosome);

//                         if gene_range.range.start == utr_range.range.end {
//                             file_range.range.start -= utr_range.range.end - utr_range.range.start;
//                         } else if gene_range.range.end == utr_range.range.start {
//                             file_range.range.end += utr_range.range.end - utr_range.range.start;
//                         } else {
//                             panic!(
//                                 "{}: {} {} {} {}",
//                                 name,
//                                 gene_range.range.start,
//                                 gene_range.range.end,
//                                 utr_range.range.start,
//                                 utr_range.range.end
//                             );
//                         }
//                     }
//                 }
//             }
//         }
//     }
// }

// fn name_ranges(
//     found: &HashMap<String, FileRange>,
//     genomic: &HashMap<String, Fasta>,
//     coding: Option<&HashMap<String, Fasta>>,
//     utr5: Option<&HashMap<String, Fasta>>,
//     utr3: Option<&HashMap<String, Fasta>>,
// ) -> HashMap<String, Vec<NamedFileRange>> {
//     for (name, file_range) in found.iter() {
//         if let Some(gene) = genomic.get(name) {
//             if let Some(coding) = coding {
//                 if let Some(gene_coding) = coding.get(name) {
//                     if let Some(coding_ranges) = gene_coding.coding_ranges() {
//                         // TODO!!!
//                     }
//                 }
//             }
//         }
//     }

//     HashMap::new()
// }

// fn update_utrs(
//     genome: &str,
//     orf: &HashMap<String, FileRange>,
//     utr: &HashMap<String, Fasta>,
//     tag: &str,
// ) -> HashMap<String, (String, usize, usize)> {
//     let mut result: HashMap<String, (String, usize, usize)> = HashMap::new();

//     for (name, range) in orf.iter_mut() {
//         let from = range.from;
//         let to = range.to;

//         if let Some(fasta) = utr.get(name) {
//             let sequence = fasta.sequence.to_owned();
//             if genome[from - sequence.len() + 1..from + 1] == sequence {
//                 range.from -= sequence.len() + 1;
//                 for subrange in range.subranges.iter_mut() {
//                     subrange.from += sequence.len() + 1;
//                     subrange.to += sequence.len() + 1;
//                 }
//                 range.subranges.push(Subrange {
//                     from: 0,
//                     to: sequence.len(),
//                     tag: tag.to_string(),
//                 });
//                 continue;
//             }
//             if genome[to - 1..to + fasta.sequence.len() - 1] == sequence {
//                 range.subranges.push(Subrange {
//                     from: -sequence.len() + 1,
//                     to: 0,
//                     tag: tag.to_string(),
//                 });
//                 result.insert(name.to_string(), (tag.to_owned(), to, to + sequence.len()));
//                 continue;
//             }

//             let sequence = fasta.reversed().sequence;
//             if genome[from - sequence.len() + 1..from + 1] == sequence {
//                 result.insert(
//                     name.to_string(),
//                     (tag.to_owned(), from - sequence.len(), from),
//                 );
//                 continue;
//             }
//             if genome[to - 1..to + fasta.sequence.len() - 1] == sequence {
//                 result.insert(name.to_string(), (tag.to_owned(), to, to + sequence.len()));
//                 continue;
//             }

//             println!("Failed to find UTR position for {}", name);
//         }
//     }

//     result
// }

// fn determine_order(
//     orf: &HashMap<String, (usize, usize)>,
//     rna: &HashMap<String, (usize, usize)>,
//     other: &HashMap<String, (usize, usize)>,
// ) -> HashMap<String, usize> {
//     let all: Vec<(&String, &(usize, usize))> =
//         orf.iter().chain(rna.iter()).chain(other.iter()).collect();
//     let mut graph: HashMap<String, Vec<String>> = HashMap::new();

//     for i in 0..all.len() {
//         let (name_i, (from_i, to_i)) = all[i];

//         for j in i + 1..all.len() {
//             let (name_j, (from_j, to_j)) = all[j];

//             if (from_i <= from_j && from_j <= to_i)
//                 || (from_i <= to_j && to_j <= to_i)
//                 || (from_j <= from_i && from_i <= to_j)
//                 || (from_j <= to_i && to_i <= to_j)
//             {
//                 if !graph.contains_key(name_i) {
//                     graph.insert(name_i.to_string(), Vec::new());
//                 }
//                 graph.get_mut(name_i).unwrap().push(name_j.to_string());

//                 if !graph.contains_key(name_j) {
//                     graph.insert(name_j.to_string(), Vec::new());
//                 }
//                 graph.get_mut(name_j).unwrap().push(name_i.to_string());
//             }
//         }
//     }

//     let mut order: Vec<(String, usize)> = graph
//         .iter()
//         .map(|(name, vec)| (name.to_string(), vec.len()))
//         .collect();
//     order.sort_by_key(|(_, len)| len.to_owned());

//     let mut colors: HashMap<String, usize> = HashMap::new();

//     for (current, _) in order {
//         let mut available = [true, true, true, true, true, true, true, true];

//         if let Some(adjacent) = graph.get(&current) {
//             for next in adjacent {
//                 if let Some(color) = colors.get(next) {
//                     available[*color] = false;
//                 }
//             }
//         }

//         for i in 0..available.len() {
//             if available[i] {
//                 colors.insert(current.to_string(), i);
//                 break;
//             }
//         }
//     }

//     for (name, _) in all.iter() {
//         if !colors.contains_key(*name) {
//             colors.insert(name.to_string(), 0);
//         }
//     }

//     colors
// }

// fn fill_annotations(
//     annotations: &mut Vec<Vec<String>>,
//     category: &str,
//     found: &HashMap<String, (usize, usize)>,
//     order: &HashMap<String, usize>,
//     genomic: &HashMap<String, Fasta>,
//     coding: Option<&HashMap<String, Fasta>>,
// ) {
//     for (name, (from, to)) in found.iter() {
//         let index = order.get(name).unwrap();

//         for i in *from..*to {
//             let annotation = annotations.get_mut(i).unwrap();
//             annotation[4 * index] = category.to_string();
//             annotation[4 * index + 1] = "Intron".to_string();
//             annotation[4 * index + 2] = genomic.get(name).unwrap().systematic_name();
//             annotation[4 * index + 3] = genomic.get(name).unwrap().standard_name();
//         }

//         if let Some(coding) = coding {
//             if let Some(vec) = coding.get(name).unwrap().coding_ranges() {
//                 for (i, j) in vec.iter() {
//                     for k in *i..=*j {
//                         let annotation = annotations.get_mut(*from + k).unwrap();
//                         annotation[4 * index + 1] = "Exon".to_string();
//                     }
//                 }
//             }
//         }
//     }
// }

// fn store_result(input: &Path, output: &Path, annotations: Vec<Vec<String>>) {
//     let fin = File::open(input).unwrap();
//     let fout = File::create(output).unwrap();
//     let reader = BufReader::new(fin);
//     let mut writer = BufWriter::new(fout);

//     let count = annotations.get(0).unwrap().len();

//     for (i, line) in reader.lines().enumerate() {
//         if let Ok(line) = line {
//             writer.write(line.as_bytes());

//             if i == 0 {
//                 for j in 0..count / 4 {
//                     writer.write(
//                         format!(
//                             "\tType {}\tSubtype {}\tSystematic name {}\tStandard name {}\n",
//                             j + 1,
//                             j + 1,
//                             j + 1,
//                             j + 1
//                         )
//                         .as_bytes(),
//                     );
//                 }
//             } else {
//                 for annotation in annotations.get(i - 1).unwrap().iter() {
//                     writer.write(b"\t");
//                     writer.write(annotation.as_bytes());
//                 }
//                 writer.write(b"\n");
//             }
//         }
//     }
// }

// fn update_exon_intron(found: &HashMap<String, Range>, coding: &HashMap<String, Fasta>) {
//     for (name, range) in found.iter_mut() {
//         if let Some(coding) = coding.get(name) {
//             range.subranges = Vec::new();

//             let mut marks = Vec::new();
//             let mut exons = HashSet::new();

//             if let Some(vec) = coding.coding_ranges() {
//                 for (i, j) in vec.iter() {
//                     marks.push(*i);
//                     marks.push(*j);
//                     exons.insert((*i, *j));
//                 }
//             }

//             marks.sort();

//             for i in 1..marks.len() {
//                 let tuple = (marks[i - 1], marks[i]);

//                 if exons.contains(&tuple) {
//                     range.subranges.push(Subrange {
//                         from: tuple.0,
//                         to: tuple.1,
//                         tag: "Exon".to_string(),
//                     })
//                 } else {
//                     range.subranges.push(Subrange {
//                         from: tuple.0,
//                         to: tuple.1,
//                         tag: "Intron".to_string(),
//                     })
//                 }
//             }
//         }
//     }
// }

// fn load_utr(path: &Path) -> HashMap<String, Fasta> {
//     load_fasta_gz(path)
//         .into_iter()
//         .map(|(_, fasta)| (fasta.systematic_name_for_utr(), fasta))
//         .collect()
// }

fn main() {
    let args = Args::parse();

    let genome = load_genome(Path::new(&args.input));

    // let orf_genomic_1000 = load_fasta_gz(Path::new("../data/orf_genomic_1000.fasta.gz"));
    // let orf_genomic = load_fasta_gz(Path::new("../data/orf_genomic.fasta.gz"));
    // let orf_coding = load_fasta_gz(Path::new("../data/orf_coding.fasta.gz"));
    // let mut orf_found = collect_found(&genome, &orf_genomic_1000, &orf_genomic);

    // let utr5p = load_utr(Path::new("../data/SGD_all_ORFs_5prime_UTRs.fsa.gz"));
    // let utr3p = load_utr(Path::new("../data/SGD_all_ORFs_3prime_UTRs.fsa.gz"));
    // extend_file_range_for_utr(&mut orf_found, &orf_genomic, &utr5p);
    // extend_file_range_for_utr(&mut orf_found, &orf_genomic, &utr3p);

    // let file_range = orf_found.get("YAL064W-B").unwrap();
    // let orf_named = name_ranges(&orf_found, Some(&orf_coding), Some(&utr5p), Some(&utr3p));

    // return;

    //  update_exon_intron(&orf_found, &orf_coding);
    // update_utrs(&genome, &orf_found, &utr5p, "UTR 5'");
    // update_utrs(&genome, &orf_found, &utr3p, "UTR 3'");

    // let rna_genomic_1000 = load_fasta_gz(Path::new("../data/rna_genomic_1000.fasta.gz"));
    // let rna_genomic = load_fasta_gz(Path::new("../data/rna_genomic.fasta.gz"));
    // let rna_coding = load_fasta_gz(Path::new("../data/rna_coding.fasta.gz"));
    // let rna_found = collect_found(&genome, &rna_genomic_1000, &rna_genomic);
    // update_exon_intron(&orf_found, &orf_coding);

    // let other_genomic_1000 =
    //     load_fasta_gz(Path::new("../data/other_features_genomic_1000.fasta.gz"));
    // let other_genomic = load_fasta_gz(Path::new("../data/other_features_genomic.fasta.gz"));

    // let other_found = collect_found(&genome, &other_genomic_1000, &other_genomic);

    // let utr5p_found = update_utrs(&genome, &orf_found, &utr5p, "UTR 5'");
    // let utr3p_found = update_utrs(&genome, &orf_found, &utr3p, "UTR 3'");

    // let order = determine_order(&orf_found, &rna_found, &other_found);
    // let max = order.iter().map(|(_, order)| order).max().unwrap();

    // let mut annotations = Vec::with_capacity(genome.len());
    // for _ in 0..genome.len() {
    //     let arr = vec!["?".to_string(); 4 * (max + 1)];
    //     annotations.push(arr);
    // }

    // fill_annotations(
    //     &mut annotations,
    //     "ORF",
    //     &orf_found,
    //     &order,
    //     &orf_genomic,
    //     Some(&orf_coding),
    // );
    // fill_annotations(
    //     &mut annotations,
    //     "RNA",
    //     &rna_found,
    //     &order,
    //     &rna_genomic,
    //     Some(&rna_coding),
    // );
    // fill_annotations(
    //     &mut annotations,
    //     "Other",
    //     &other_found,
    //     &order,
    //     &other_genomic,
    //     None,
    // );

    // store_result(Path::new(&args.input), Path::new(&args.output), annotations);
}
