use anyhow::{Context, Result};
use clap::Parser;
use std::{
    collections::HashMap,
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};

use sgd_annotator::data::ensure_all_data;
use sgd_annotator::fasta::{load_fasta_gz, load_utr_fasta_gz, Fasta};
use sgd_annotator::translator::GenomicRange;

#[derive(Parser)]
#[command(version)]
struct Args {
    #[arg(short, long)]
    output: String,
}

fn format_locus(range: &GenomicRange) -> String {
    format!(
        "Chr {} from {}-{}",
        range.chromosome, range.start, range.end
    )
}

fn format_locus_ranges(ranges: &[GenomicRange]) -> String {
    let chromosome = ranges
        .get(0)
        .expect("Expected at least one range")
        .chromosome
        .clone();
    let joined = ranges
        .iter()
        .map(|range| format!("{}-{}", range.start, range.end))
        .collect::<Vec<_>>()
        .join(",");
    format!("Chr {} from {}", chromosome, joined)
}

fn write_record(
    writer: &mut BufWriter<File>,
    name: &str,
    category: &str,
    locus: &str,
    sequence: &str,
) -> Result<()> {
    writeln!(writer, ">{name}|{category}|{locus}")?;
    writeln!(writer, "{sequence}")?;
    Ok(())
}

fn sorted_entries(map: &HashMap<String, Fasta>) -> Vec<(&String, &Fasta)> {
    let mut entries: Vec<(&String, &Fasta)> = map.iter().collect();
    entries.sort_by(|(a, _), (b, _)| a.cmp(b));
    entries
}

fn intron_sequence(
    genomic_fasta: &Fasta,
    genomic_range: &GenomicRange,
    intron_range: &GenomicRange,
    oriented_genomic_sequence: &str,
) -> Result<String> {
    if intron_range.start < genomic_range.start || intron_range.end > genomic_range.end {
        anyhow::bail!(
            "Intron range {}-{} outside genomic range {}-{} for {}",
            intron_range.start,
            intron_range.end,
            genomic_range.start,
            genomic_range.end,
            genomic_fasta.systematic_name()
        );
    }

    let is_reverse = genomic_fasta.header.contains("reverse complement");
    let start_idx = if is_reverse {
        genomic_range.end - intron_range.end
    } else {
        intron_range.start - genomic_range.start
    };
    let end_idx = if is_reverse {
        genomic_range.end - intron_range.start
    } else {
        intron_range.end - genomic_range.start
    };

    let bytes = oriented_genomic_sequence.as_bytes();
    if end_idx >= bytes.len() {
        anyhow::bail!(
            "Computed intron slice {}-{} out of bounds for {}",
            start_idx,
            end_idx,
            genomic_fasta.systematic_name()
        );
    }
    let slice = &bytes[start_idx..=end_idx];
    Ok(String::from_utf8(slice.to_vec())?)
}

fn emit_genomic(
    writer: &mut BufWriter<File>,
    entries: &HashMap<String, Fasta>,
) -> Result<()> {
    for (_, fasta) in sorted_entries(entries) {
        let locus = format_locus(&fasta.genomic_range());
        write_record(
            writer,
            &fasta.systematic_name(),
            "Genomic",
            &locus,
            &fasta.sequence(),
        )?;
    }
    Ok(())
}

fn emit_utr(
    writer: &mut BufWriter<File>,
    entries: &HashMap<String, Fasta>,
    category: &str,
) -> Result<()> {
    for (_, fasta) in sorted_entries(entries) {
        let locus = format_locus(&fasta.genomic_range());
        write_record(
            writer,
            &fasta.systematic_name(),
            category,
            &locus,
            &fasta.sequence(),
        )?;
    }
    Ok(())
}

fn emit_coding_and_introns(
    writer: &mut BufWriter<File>,
    entries: &HashMap<String, Fasta>,
    genomic_lookup: &HashMap<String, &Fasta>,
) -> Result<()> {
    for (name, fasta) in sorted_entries(entries) {
        let coding_ranges = fasta
            .coding_ranges()
            .context("Expected coding ranges for CDS FASTA")?;
        let locus = format_locus_ranges(&coding_ranges);
        write_record(writer, name, "CDS", &locus, &fasta.sequence())?;

        if let Some(noncoding_ranges) = fasta.noncoding_ranges() {
            let genomic_fasta = genomic_lookup
                .get(name)
                .with_context(|| format!("Missing genomic FASTA for {name}"))?;
            let genomic_range = genomic_fasta.genomic_range();
            let oriented_genomic_sequence = genomic_fasta.sequence();
            let is_reverse = genomic_fasta.header.contains("reverse complement");
            let mut intron_ranges = noncoding_ranges;
            if is_reverse {
                intron_ranges.reverse();
            }

            for (index, intron_range) in intron_ranges.iter().enumerate() {
                let locus = format_locus(intron_range);
                let sequence = intron_sequence(
                    genomic_fasta,
                    &genomic_range,
                    intron_range,
                    &oriented_genomic_sequence,
                )?;
                let intron_name = format!("{}.intron{}", name, index + 1);
                write_record(writer, &intron_name, "Intron", &locus, &sequence)?;
            }
        }
    }
    Ok(())
}

fn main() -> Result<()> {
    let args = Args::parse();

    ensure_all_data()?;

    let orf_genomic = load_fasta_gz(Path::new("data/orf_genomic.fasta.gz"));
    let rna_genomic = load_fasta_gz(Path::new("data/rna_genomic.fasta.gz"));
    let other_genomic = load_fasta_gz(Path::new("data/other_features_genomic.fasta.gz"));
    let orf_coding = load_fasta_gz(Path::new("data/orf_coding.fasta.gz"));
    let rna_coding = load_fasta_gz(Path::new("data/rna_coding.fasta.gz"));
    let utr5p = load_utr_fasta_gz(Path::new("data/5prime_utr.fsa.gz"));
    let utr3p = load_utr_fasta_gz(Path::new("data/3prime_utr.fsa.gz"));

    let mut genomic_lookup: HashMap<String, &Fasta> = HashMap::new();
    for (name, fasta) in orf_genomic.iter() {
        genomic_lookup.insert(name.clone(), fasta);
    }
    for (name, fasta) in rna_genomic.iter() {
        genomic_lookup.insert(name.clone(), fasta);
    }
    for (name, fasta) in other_genomic.iter() {
        genomic_lookup.insert(name.clone(), fasta);
    }

    let output = File::create(&args.output)
        .with_context(|| format!("Failed to create output file {}", args.output))?;
    let mut writer = BufWriter::new(output);

    emit_genomic(&mut writer, &orf_genomic)?;
    emit_genomic(&mut writer, &rna_genomic)?;
    emit_genomic(&mut writer, &other_genomic)?;
    emit_coding_and_introns(&mut writer, &orf_coding, &genomic_lookup)?;
    emit_coding_and_introns(&mut writer, &rna_coding, &genomic_lookup)?;
    emit_utr(&mut writer, &utr5p, "5UTR")?;
    emit_utr(&mut writer, &utr3p, "3UTR")?;

    Ok(())
}
