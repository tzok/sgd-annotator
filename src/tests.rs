#[cfg(test)]
mod tests {
    use crate::{
        fasta::load_fasta_gz,
        fasta::{load_utr_fasta_gz, FastaType},
        translator::GenomicRange,
    };
    use std::path::Path;

    #[test]
    fn fasta_genomic() {
        let fasta = load_fasta_gz(Path::new("../data/tests/genomic.fasta.gz"));

        assert!(fasta.contains_key("YAL068C"));
        let gene = fasta.get("YAL068C").unwrap();
        assert_eq!(gene.fasta_type(), FastaType::Gene);
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
        assert_eq!(gene.fasta_type(), FastaType::Gene);
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

    #[test]
    fn fasta_coding() {
        let fasta = load_fasta_gz(Path::new("../data/tests/coding.fasta.gz"));

        assert!(fasta.contains_key("YAL068C"));
        let gene = fasta.get("YAL068C").unwrap();
        assert_eq!(gene.fasta_type(), FastaType::Gene);
        assert_eq!(gene.systematic_name(), "YAL068C");
        assert_eq!(gene.standard_name(), "PAU8");
        let coding = gene.coding_ranges();
        assert!(coding.is_some());
        let coding = coding.unwrap();
        assert_eq!(coding.len(), 1);
        assert_eq!(
            coding[0],
            GenomicRange {
                chromosome: "I".to_string(),
                range: 1807..2169,
            }
        );
        let noncoding = gene.noncoding_ranges();
        assert!(noncoding.is_none());

        assert!(fasta.contains_key("YAL003W"));
        let gene = fasta.get("YAL003W").unwrap();
        assert_eq!(gene.fasta_type(), FastaType::Gene);
        assert_eq!(gene.systematic_name(), "YAL003W");
        assert_eq!(gene.standard_name(), "EFB1");
        let coding = gene.coding_ranges();
        assert!(coding.is_some());
        let coding = coding.unwrap();
        assert_eq!(coding.len(), 2);
        assert_eq!(
            coding[0],
            GenomicRange {
                chromosome: "I".to_string(),
                range: 142174..142253,
            }
        );
        assert_eq!(
            coding[1],
            GenomicRange {
                chromosome: "I".to_string(),
                range: 142620..143160,
            }
        );
        let noncoding = gene.noncoding_ranges();
        assert!(noncoding.is_some());
        let noncoding = noncoding.unwrap();
        assert_eq!(noncoding.len(), 1);
        assert_eq!(
            noncoding[0],
            GenomicRange {
                chromosome: "I".to_string(),
                range: 142253..142620,
            }
        );
    }

    #[test]
    fn fasta_utr() {
        let fasta = load_utr_fasta_gz(Path::new("../data/tests/utr.fasta.gz"));

        assert!(fasta.contains_key("YAL067C"));
        let gene = fasta.get("YAL067C").unwrap();
        assert_eq!(gene.fasta_type(), FastaType::UTR);
        assert_eq!(gene.systematic_name(), "YAL067C");
        assert_eq!(
            gene.genomic_range(),
            GenomicRange {
                chromosome: "I".to_string(),
                range: 9016..9049,
            }
        );

        assert!(fasta.contains_key("YAL066W"));
        let gene = fasta.get("YAL066W").unwrap();
        assert_eq!(gene.fasta_type(), FastaType::UTR);
        assert_eq!(gene.systematic_name(), "YAL066W");
        assert_eq!(
            gene.genomic_range(),
            GenomicRange {
                chromosome: "I".to_string(),
                range: 9807..10091,
            }
        );
    }

    #[test]
    fn fasta_chromosome() {
        let fasta = load_fasta_gz(Path::new("../data/chr01.fsa.gz"));

        assert!(fasta.contains_key("chrI"));
        let gene = fasta.get("chrI").unwrap();
        assert_eq!(gene.fasta_type(), FastaType::Chromosome);
        assert_eq!(gene.systematic_name(), "chrI");
        assert_eq!(
            gene.genomic_range(),
            GenomicRange {
                chromosome: "I".to_string(),
                range: 1..230219,
            }
        );
    }
}
