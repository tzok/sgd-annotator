#[cfg(test)]
mod tests {
    use crate::{fasta::load_fasta_gz, translator::GenomicRange};
    use std::path::Path;

    #[test]
    fn fasta_genomic() {
        let fasta = load_fasta_gz(Path::new("../data/tests/genomic.fasta.gz"));

        assert!(fasta.contains_key("YAL068C"));
        let gene = fasta.get("YAL068C").unwrap();
        assert_eq!(gene.systematic_name(), "YAL068C");
        assert_eq!(gene.standard_name(), "PAU8");
        assert_eq!(
            gene.genomic_range(),
            GenomicRange {
                chromosome: "I".to_string(),
                range: 1807..2169,
            }
        );

        assert!(fasta.contains_key("YAL067W-A"));
        let gene = fasta.get("YAL067W-A").unwrap();
        assert_eq!(gene.systematic_name(), "YAL067W-A");
        assert_eq!(gene.standard_name(), "YAL067W-A");
        assert_eq!(
            gene.genomic_range(),
            GenomicRange {
                chromosome: "I".to_string(),
                range: 2480..2707,
            }
        );
    }
}
