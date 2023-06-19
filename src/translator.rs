use std::ops::Range;

#[derive(Debug, PartialEq, Eq)]
pub struct GenomicRange {
    pub chromosome: String,
    pub range: Range<usize>,
}

pub struct FileRange {
    pub range: Range<usize>,
}
