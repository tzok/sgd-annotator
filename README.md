# sgd-annotator

Tools for working with Saccharomyces cerevisiae (SGD) reference data. The project downloads reference FASTA and annotation data, then provides CLI utilities for annotation export, GO slim mapping, and normalized FASTA generation.

## Data sources

On first run, binaries download data from the SGD archive into `data/` (chromosomes, ORF/RNA FASTA, UTRs, and GO slim mapping). This is handled by `sgd_annotator::data::ensure_all_data`.

## Binaries

### `sgd-annotator`

Reads a genome input file (gzipped, tab-delimited per nucleotide), overlays annotations for ORFs/RNAs/UTRs/introns, and writes a TSV-like output that augments the input with annotation columns.

Optional `--bigwig` flags append one per-nucleotide column per bigWig file to the output, sampled at each genome position. Chromosome names are reconciled between the bigWig (e.g. `chrI`, `chrM`) and the SGD reference; 0-based half-open bigWig coordinates are mapped onto the 1-based genome index.

Positions with no data in the bigWig are emitted according to `--bigwig-missing` (default `zero`, emitting `0`), since Ribo-Seq coverage tracks conventionally represent zero-read positions as gaps. Other modes: `empty` (blank cell, preserves the gap-vs-zero distinction), `nan` (literal `nan`), `question` (`?`).

Example:

```bash
cargo run --bin sgd-annotator -- --input path/to/genome.txt.gz --output annotated.tsv \
    --bigwig Global_Scanning.26_07_2016.bw --bigwig-name ribo_global \
    --bigwig Sen15_All.RiboCov.bw --bigwig-missing zero
```

`--bigwig-name` (optional, paired by position with `--bigwig`) overrides the output column header, which otherwise defaults to the bigWig file stem.

`--bigwig-dir <DIR>` (repeatable) aggregates all `*.bw` files in a directory into a single summed column. Aggregation uses NaN-aware summation: a position is `NaN` in the aggregate only when every file has a gap there; otherwise the available values are summed. This preserves the "no data anywhere" signal so `--bigwig-missing` still applies meaningfully to aggregate columns. `--bigwig-dir-name <NAME>` (optional, paired by position with `--bigwig-dir`) overrides the aggregate column header, which otherwise defaults to the directory's trailing path component (e.g. `--bigwig-dir /tmp/download` → column `download`). Column order in the output: all `--bigwig` columns first (in flag order), then all `--bigwig-dir` aggregate columns (in flag order).

Example:

```bash
cargo run --bin sgd-annotator -- --input path/to/genome.txt.gz --output annotated.tsv \
    --bigwig Global_Scanning.26_07_2016.bw --bigwig-name ribo_global \
    --bigwig-dir /tmp/download --bigwig-dir-name ribo_total --bigwig-missing zero
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
