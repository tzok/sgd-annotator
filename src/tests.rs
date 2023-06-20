#[cfg(test)]
mod tests {
    use crate::{
        fasta::load_fasta_gz,
        fasta::{load_utr_fasta_gz, FastaType},
        load_genome_gz,
        translator::{GenomicRange, Translator, YeastChromosome},
    };
    use std::path::Path;

    fn init() {
        let _ = env_logger::builder().is_test(true).try_init();
    }

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
                chromosome: YeastChromosome::I,
                start: 1807,
                end: 2169,
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
                chromosome: YeastChromosome::I,
                start: 2480,
                end: 2707,
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
                chromosome: YeastChromosome::I,
                start: 1807,
                end: 2169,
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
                chromosome: YeastChromosome::I,
                start: 142174,
                end: 142253,
            }
        );
        assert_eq!(
            coding[1],
            GenomicRange {
                chromosome: YeastChromosome::I,
                start: 142620,
                end: 143160,
            }
        );
        let noncoding = gene.noncoding_ranges();
        assert!(noncoding.is_some());
        let noncoding = noncoding.unwrap();
        assert_eq!(noncoding.len(), 1);
        assert_eq!(
            noncoding[0],
            GenomicRange {
                chromosome: YeastChromosome::I,
                start: 142253,
                end: 142620,
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
                chromosome: YeastChromosome::I,
                start: 9016,
                end: 9049,
            }
        );

        assert!(fasta.contains_key("YAL066W"));
        let gene = fasta.get("YAL066W").unwrap();
        assert_eq!(gene.fasta_type(), FastaType::UTR);
        assert_eq!(gene.systematic_name(), "YAL066W");
        assert_eq!(
            gene.genomic_range(),
            GenomicRange {
                chromosome: YeastChromosome::I,
                start: 9807,
                end: 10091,
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
                chromosome: YeastChromosome::I,
                start: 1,
                end: 230218,
            }
        );
    }

    #[test]
    fn translator() {
        init();

        let genome = load_genome_gz(Path::new("../data/tests/genome.txt.gz"));
        let sample = load_fasta_gz(Path::new("../data/tests/sample.fasta.gz"));
        let translator = Translator::new(&genome);

        for name in vec![
            "YAL068C",
            "YBL113C",
            "YCL075W",
            "YDL248W",
            "YEL077C",
            "YFL068W",
            "YGL263W",
            "YHL050C",
            "YIL177C",
            "YJL225C",
            "YKL224C",
            "YLL067C",
            "YML133C",
            "YNL339C",
            "YOL166W-A",
            "YPL283C",
            "YAL067C",
            "YEL073C",
            "YJL218W",
            "YBL107C",
            "YDR227W",
            "YIL169C",
            "YFL062W",
            "YKL169C",
            "YOL164W",
            "YCL064C",
            "YGL263W",
            "YLR125W",
            "YNL331C",
            "YPL279C",
            "YHL048W",
            "YML131W",
        ] {
            assert!(sample.contains_key(name));
            let fasta = sample.get(name).unwrap();
            let range = fasta.genomic_range();
            let start = translator
                .translate_nt(&range.chromosome, range.start)
                .unwrap();
            let end = translator
                .translate_nt(&range.chromosome, range.end)
                .unwrap();
            assert_eq!(&genome[start..=end], fasta.sequence());
        }
    }
}
