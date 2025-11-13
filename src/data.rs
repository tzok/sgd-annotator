use std::{fs, path::PathBuf};
use anyhow::Result;

pub fn data_dir() -> PathBuf {
    let mut path = std::env::current_dir().unwrap();
    path.push("data");
    if !path.exists() {
        std::fs::create_dir_all(&path).unwrap();
    }
    path
}

pub fn ensure_file(url: &str, filename: &str) -> Result<()> {
    let data_dir = data_dir();
    let file_path = data_dir.join(filename);

    if file_path.exists() {
        return Ok(());
    }

    let response = reqwest::blocking::get(url)?;
    let bytes = response.bytes()?;
    fs::write(&file_path, bytes)?;
    Ok(())
}

pub fn ensure_chromosome(url: &str, name: &str) -> Result<()> {
    let data_dir = data_dir();
    let fsa_path = data_dir.join(name).with_extension("fsa");
    let gz_path = data_dir.join(name).with_extension("fsa.gz");

    if gz_path.exists() {
        return Ok(());
    }

    if !fsa_path.exists() {
        let response = reqwest::blocking::get(url)?;
        let bytes = response.bytes()?;
        fs::write(&fsa_path, bytes)?;
    }

    // Compress to gz
    let content = fs::read_to_string(&fsa_path)?;
    let gz_file = fs::File::create(&gz_path)?;
    let mut encoder = flate2::write::GzEncoder::new(gz_file, flate2::Compression::default());
    encoder.write_all(content.as_bytes())?;
    encoder.finish()?;

    fs::remove_file(fsa_path)?;
    Ok(())
}

pub fn ensure_utr_downloaded(url: &str, output_filename: &str) -> Result<()> {
    let data_dir = data_dir();
    let zip_path = data_dir.join(output_filename).with_extension("zip");
    let gz_path = data_dir.join(output_filename).with_extension("fsa.gz");

    if gz_path.exists() {
        return Ok(());
    }

    // Download the zip if it doesn't exist.
    if !zip_path.exists() {
        let response = reqwest::blocking::get(url)?;
        let bytes = response.bytes()?;
        fs::write(&zip_path, bytes)?;
    }

    // Extract the zip
    let zip_file = fs::File::open(zip_path)?;
    let mut archive = zip::ZipArchive::new(zip_file)?;
    let mut first_file = archive.by_index(0)?;
    let content = first_file.bytes().collect::<Result<Vec<_>, _>>()?;

    // Write as gzipped
    let gz_file = fs::File::create(&gz_path)?;
    let mut encoder = flate2::write::GzEncoder::new(gz_file, flate2::Compression::default());
    encoder.write_all(&content)?;
    encoder.finish()?;

    // We can remove the zip file to save space?
    fs::remove_file(zip_path)?;

    Ok(())
}

pub fn ensure_all_data() -> Result<()> {
    // Chromosomes
    let chromosomes = [
        ("chr01", "http://sgd-archive.yeastgenome.org/sequence/S288C_reference/chromosomes/fasta/chr01.fsa"),
        ("chr02", "http://sgd-archive.yeastgenome.org/sequence/S288C_reference/chromosomes/fasta/chr02.fsa"),
        ("chr03", "http://sgd-archive.yeastgenome.org/sequence/S288C_reference/chromosomes/fasta/chr03.fsa"),
        ("chr04", "http://sgd-archive.yeastgenome.org/sequence/S288C_reference/chromosomes/fasta/chr04.fsa"),
        ("chr05", "http://sgd-archive.yeastgenome.org/sequence/S288C_reference/chromosomes/fasta/chr05.fsa"),
        ("chr06", "http://sgd-archive.yeastgenome.org/sequence/S288C_reference/chromosomes/fasta/chr06.fsa"),
        ("chr07", "http://sgd-archive.yeastgenome.org/sequence/S288C_reference/chromosomes/fasta/chr07.fsa"),
        ("chr08", "http://sgd-archive.yeastgenome.org/sequence/S288C_reference/chromosomes/fasta/chr08.fsa"),
        ("chr09", "http://sgd-archive.yeastgenome.org/sequence/S288C_reference/chromosomes/fasta/chr09.fsa"),
        ("chr10", "http://sgd-archive.yeastgenome.org/sequence/S288C_reference/chromosomes/fasta/chr10.fsa"),
        ("chr11", "http://sgd-archive.yeastgenome.org/sequence/S288C_reference/chromosomes/fasta/chr11.fsa"),
        ("chr12", "http://sgd-archive.yeastgenome.org/sequence/S288C_reference/chromosomes/fasta/chr12.fsa"),
        ("chr13", "http://sgd-archive.yeastgenome.org/sequence/S288C_reference/chromosomes/fasta/chr13.fsa"),
        ("chr14", "http://sgd-archive.yeastgenome.org/sequence/S288C_reference/chromosomes/fasta/chr14.fsa"),
        ("chr15", "http://sgd-archive.yeastgenome.org/sequence/S288C_reference/chromosomes/fasta/chr15.fsa"),
        ("chr16", "http://sgd-archive.yeastgenome.org/sequence/S288C_reference/chromosomes/fasta/chr16.fsa"),
        ("chrmt", "http://sgd-archive.yeastgenome.org/sequence/S288C_reference/chromosomes/fasta/chrmt.fsa"),
    ];
    for (name, url) in &chromosomes {
        ensure_chromosome(url, name)?;
    }

    // Other files
    let files = [
        ("orf_genomic.fasta.gz", "http://sgd-archive.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_genomic.fasta.gz"),
        ("orf_coding.fasta.gz", "http://sgd-archive.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_coding.fasta.gz"),
        ("rna_genomic.fasta.gz", "http://sgd-archive.yeastgenome.org/sequence/S288C_reference/rna/rna_genomic.fasta.gz"),
        ("rna_coding.fasta.gz", "http://sgd-archive.yeastgenome.org/sequence/S288C_reference/rna/rna_coding.fasta.gz"),
        ("other_features_genomic.fasta.gz", "http://sgd-archive.yeastgenome.org/sequence/S288C_reference/other_features/other_features_genomic.fasta.gz"),
    ];
    for (filename, url) in &files {
        ensure_file(url, filename)?;
    }

    // UTRs
    let utrs = [
        ("5prime_utr", "http://sgd-archive.yeastgenome.org/sequence/S288C_reference/SGD_all_ORFs_5prime_UTRs.fsa.zip"),
        ("3prime_utr", "http://sgd-archive.yeastgenome.org/sequence/S288C_reference/SGD_all_ORFs_3prime_UTRs.fsa.zip"),
    ];
    for (name, url) in &utrs {
        ensure_utr_downloaded(url, name)?;
    }

    // Also go_slim_mapping.tab
    ensure_file("http://sgd-archive.yeastgenome.org/curation/literature/go_slim_mapping.tab", "go_slim_mapping.tab")?;

    // Translational efficiency files (from the code in go-mapper.rs)
    ensure_file("../data/translational-efficiency-csardi.csv", "translational-efficiency-csardi.csv")?;
    ensure_file("../data/translational-efficiency-weinberg.csv", "translational-efficiency-weinberg.csv")?;

    Ok(())
}
