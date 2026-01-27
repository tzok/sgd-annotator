# sgd-annotator

Tools for working with Saccharomyces cerevisiae (SGD) reference data. The project downloads reference FASTA and annotation data, then provides CLI utilities for annotation export, GO slim mapping, and normalized FASTA generation.

## Data sources

On first run, binaries download data from the SGD archive into `data/` (chromosomes, ORF/RNA FASTA, UTRs, and GO slim mapping). This is handled by `sgd_annotator::data::ensure_all_data`.

## Binaries

### `sgd-annotator`

Reads a genome input file (gzipped, tab-delimited per nucleotide), overlays annotations for ORFs/RNAs/UTRs/introns, and writes a TSV-like output that augments the input with annotation columns.

Example:

```bash
cargo run --bin sgd-annotator -- --input path/to/genome.txt.gz --output annotated.tsv
```

### `go-mapper`

Builds a CSV of GO slim terms and translational efficiency metrics for ORFs.

Example:

```bash
cargo run --bin go-mapper
```

Output: `go_mapper.csv`

### `normalize-fasta`

Creates a normalized FASTA that includes entries from all non-genome data files (genomic, CDS, intron, 5UTR, 3UTR). Headers use the format:

```
>systematic name|category|locus
```

Introns are inferred by subtracting coding ranges from genomic ranges and are suffixed as `.intron1`, `.intron2`, etc.

Example:

```bash
cargo run --bin normalize-fasta -- --output normalized.fasta
```

## Notes

- The code normalizes sequences to RNA (T -> U) and orients sequences according to strand when applicable.
- Mitochondrial entries use `Chr Mito` in loci.
