//! A variant file read into memory once, so that a cohort can take each
//! sample's alleles from it without reading the file again for every sample.
//!
//! Reading the file was the smaller part of the cost. Every record of a cohort
//! VCF carries a column per sample, and splitting the whole line to reach one
//! sample's column made each sample pay for all the others: the work grew with
//! sites x samples x samples. Here the columns of each line are found once, and
//! a sample's column is a lookup.
//!
//! The two readers keep their own rules. The plain-text reader splits a line as
//! it is written; the one for `.vcf.gz` and BCF trims it first. Each line keeps
//! the split its own reader makes, so a parser sees from memory exactly the
//! columns it would have read from the file.

use super::vcf::is_column_header;
use crate::error::AppResult;
use std::io::{BufRead, BufReader};
use std::sync::Arc;

/// How a reader splits a line into columns.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Splitting {
    /// The plain-text reader splits the line as it is written.
    AsWritten,
    /// The reader for compressed and binary files trims the line first, and
    /// skips it when what is left is blank or a comment.
    Trimmed,
}

/// The columns of a record that a parser reads: the first nine, how many the
/// line has in all, and the selected sample's own.
pub(crate) struct RecordColumns<'a> {
    /// The fixed columns, `CHROM` to `FORMAT`; empty past the end of the line.
    pub(crate) fixed: [&'a str; 9],
    /// How many columns splitting the line on tabs gives.
    pub(crate) count: usize,
    /// The column asked for, if the line has one there.
    pub(crate) sample: Option<&'a str>,
}

impl<'a> RecordColumns<'a> {
    /// Split `text` on tabs, keeping the fixed columns and the one at
    /// `sample_column`.
    pub(crate) fn split(text: &'a str, sample_column: Option<usize>) -> Self {
        let mut fixed = [""; 9];
        let mut count = 0;
        let mut sample = None;
        for (index, column) in text.split('\t').enumerate() {
            if let Some(slot) = fixed.get_mut(index) {
                *slot = column;
            }
            if sample_column == Some(index) {
                sample = Some(column);
            }
            count = index + 1;
        }
        RecordColumns {
            fixed,
            count,
            sample,
        }
    }
}

/// A line as the parsers read it: its text, its columns, and for a `#CHROM`
/// line the sample names it gives.
pub(crate) trait VcfLine {
    fn text(&self) -> &str;

    /// The record's columns, split by the rule of the reader the line came
    /// from, with the column at `sample_column` picked out.
    fn record(&self, sample_column: Option<usize>) -> RecordColumns<'_>;

    /// The names a `#CHROM` line gives the samples, from its tenth column on,
    /// split as the line is written. `None` when the line has nine columns or
    /// fewer.
    fn header_samples(&self) -> Option<Arc<Vec<String>>>;
}

/// A line just read from the file.
pub(crate) struct ReadLine<'a> {
    pub(crate) text: &'a str,
    pub(crate) splitting: Splitting,
}

impl VcfLine for ReadLine<'_> {
    fn text(&self) -> &str {
        self.text
    }

    fn record(&self, sample_column: Option<usize>) -> RecordColumns<'_> {
        match self.splitting {
            Splitting::AsWritten => RecordColumns::split(self.text, sample_column),
            Splitting::Trimmed => RecordColumns::split(self.text.trim(), sample_column),
        }
    }

    fn header_samples(&self) -> Option<Arc<Vec<String>>> {
        samples_named_by(self.text)
    }
}

/// The sample names a `#CHROM` line gives, split as the line is written.
pub(crate) fn samples_named_by(header: &str) -> Option<Arc<Vec<String>>> {
    let columns: Vec<&str> = header.split('\t').collect();
    (columns.len() > 9).then(|| Arc::new(columns[9..].iter().map(|s| s.to_string()).collect()))
}

/// A line of a file held in memory, with its columns found when it was read.
pub(crate) struct HeldLine {
    text: String,
    /// Where each column of the line's split part starts, then one past the
    /// end of that part, so column `i` is `text[starts[i]..starts[i + 1] - 1]`.
    /// Empty for a line its reader never splits into a record.
    starts: Vec<u32>,
    /// What [`VcfLine::header_samples`] answers for a `#CHROM` line.
    samples: Option<Arc<Vec<String>>>,
}

impl HeldLine {
    /// `None` when the line is too long for its offsets to be held.
    fn new(text: String, splitting: Splitting) -> Option<Self> {
        let samples = if is_column_header(&text) {
            samples_named_by(&text)
        } else {
            None
        };
        // The part of the line its reader splits into a record, if any.
        let part = match splitting {
            Splitting::AsWritten => {
                (!text.starts_with("##") && !is_column_header(&text) && !text.trim().is_empty())
                    .then_some(text.as_str())
            }
            Splitting::Trimmed => {
                let trimmed = text.trim();
                (!trimmed.is_empty() && !trimmed.starts_with('#')).then_some(trimmed)
            }
        };
        let mut starts = Vec::new();
        if let Some(part) = part {
            let offset = part.as_ptr() as usize - text.as_ptr() as usize;
            u32::try_from(offset + part.len() + 1).ok()?;
            starts.push(offset as u32);
            starts.extend(
                part.bytes()
                    .enumerate()
                    .filter(|&(_, byte)| byte == b'\t')
                    .map(|(at, _)| (offset + at + 1) as u32),
            );
            starts.push((offset + part.len() + 1) as u32);
        }
        Some(HeldLine {
            text,
            starts,
            samples,
        })
    }

    fn column(&self, index: usize) -> &str {
        &self.text[self.starts[index] as usize..self.starts[index + 1] as usize - 1]
    }
}

impl VcfLine for HeldLine {
    fn text(&self) -> &str {
        &self.text
    }

    fn record(&self, sample_column: Option<usize>) -> RecordColumns<'_> {
        if self.starts.is_empty() {
            // Not a line its reader splits into a record, so not one a parser
            // asks this of; answer the way splitting it would.
            return RecordColumns::split(&self.text, sample_column);
        }
        let count = self.starts.len() - 1;
        let mut fixed = [""; 9];
        for (index, slot) in fixed.iter_mut().enumerate().take(count) {
            *slot = self.column(index);
        }
        RecordColumns {
            fixed,
            count,
            sample: sample_column
                .filter(|&index| index < count)
                .map(|index| self.column(index)),
        }
    }

    fn header_samples(&self) -> Option<Arc<Vec<String>>> {
        self.samples.clone()
    }
}

/// The most a cohort holds of a variant file, text and column offsets
/// together. A file larger than this is read again by every sample, as it was
/// before files were held, rather than take the memory of a workstation.
const MOST_HELD: u64 = 4 << 30;

/// A variant file read whole, the way its parser reads it.
pub(crate) struct VcfText {
    path: String,
    splitting: Splitting,
    lines: Vec<HeldLine>,
    /// What stopped the reading before the end, as the reader reported it. A
    /// parser reaching this point fails with it, as it would have reading the
    /// file.
    unread: Option<(std::io::ErrorKind, String)>,
}

impl VcfText {
    /// Read `path` through the reader its parser uses. `None` when the file is
    /// too large to hold, or has a line too long for its offsets, which leaves
    /// the parser to read the file itself.
    pub(crate) fn read(path: &str) -> AppResult<Option<Self>> {
        Self::read_holding_at_most(path, MOST_HELD)
    }

    fn read_holding_at_most(path: &str, most_held: u64) -> AppResult<Option<Self>> {
        let splitting = if super::vcf_fast::use_fast_parser(path) {
            Splitting::AsWritten
        } else {
            Splitting::Trimmed
        };
        let reader: Box<dyn BufRead> = match splitting {
            Splitting::AsWritten => Box::new(BufReader::with_capacity(
                64 * 1024,
                std::fs::File::open(path)
                    .map_err(|e| format!("Cannot open VCF file '{}': {}", path, e))?,
            )),
            Splitting::Trimmed => super::vcf::open_vcf_text(path)?,
        };
        let mut lines = Vec::new();
        let mut unread = None;
        let mut held_bytes = 0u64;
        for line in reader.lines() {
            match line {
                Ok(text) => match HeldLine::new(text, splitting) {
                    Some(held) => {
                        held_bytes += (held.text.len() + 4 * held.starts.len()) as u64;
                        if held_bytes > most_held {
                            return Ok(None);
                        }
                        lines.push(held);
                    }
                    None => return Ok(None),
                },
                Err(error) => {
                    unread = Some((error.kind(), error.to_string()));
                    break;
                }
            }
        }
        Ok(Some(VcfText {
            path: path.to_string(),
            splitting,
            lines,
            unread,
        }))
    }

    pub(crate) fn path(&self) -> &str {
        &self.path
    }

    pub(crate) fn splitting(&self) -> Splitting {
        self.splitting
    }

    pub(crate) fn lines(&self) -> &[HeldLine] {
        &self.lines
    }

    /// The error the reading stopped at, made again so it prints the same.
    pub(crate) fn unread(&self) -> Option<std::io::Error> {
        self.unread
            .as_ref()
            .map(|(kind, text)| std::io::Error::new(*kind, text.clone()))
    }

    /// The lines as the file gave them, ending with the error that stopped
    /// the reading, if one did.
    pub(crate) fn replay(&self) -> impl Iterator<Item = std::io::Result<&str>> + '_ {
        self.lines
            .iter()
            .map(|line| Ok(line.text.as_str()))
            .chain(self.unread().map(Err))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn held(text: &str, splitting: Splitting) -> HeldLine {
        HeldLine::new(text.to_string(), splitting).unwrap()
    }

    /// From memory a parser must see the columns it would have split from
    /// the line, whatever the line looks like: short, padded, with empty
    /// columns or a trailing tab.
    #[test]
    fn held_columns_are_the_split_columns() {
        let lines = [
            "c\t1\t.\tA\tG\t.\tPASS\tDP=3\tGT\t1\t0\t1",
            "c\t1\t.\tA\tG\t.\tPASS\tDP=3",
            "c\t1\t.\tA",
            "c\t1\t.\tA\tG\t.\tPASS\tDP=3\tGT\t1\t\t",
            "  c\t1\t.\tA\tG\t.\tPASS\t.\tGT\t0/1 \t1 ",
            "\tc\t1\t.\tA\tG\t.\tPASS\t.\tGT\t1",
            "c 1 . A G",
        ];
        for text in lines {
            for splitting in [Splitting::AsWritten, Splitting::Trimmed] {
                let read = ReadLine { text, splitting };
                let held = held(text, splitting);
                for sample_column in [None, Some(9), Some(10), Some(11), Some(12), Some(40)] {
                    let (a, b) = (read.record(sample_column), held.record(sample_column));
                    assert_eq!(a.fixed, b.fixed, "{text:?} {splitting:?}");
                    assert_eq!(a.count, b.count, "{text:?} {splitting:?}");
                    assert_eq!(
                        a.sample, b.sample,
                        "{text:?} {splitting:?} {sample_column:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn a_header_line_names_its_samples() {
        let header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2";
        let names = held(header, Splitting::AsWritten).header_samples().unwrap();
        assert_eq!(*names, vec!["S1".to_string(), "S2".to_string()]);
        assert_eq!(
            ReadLine {
                text: header,
                splitting: Splitting::Trimmed
            }
            .header_samples(),
            held(header, Splitting::Trimmed).header_samples()
        );
        let bare = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
        assert!(held(bare, Splitting::AsWritten).header_samples().is_none());
    }

    /// A file that stops being readable part way keeps the lines before the
    /// failure, and hands the failure back where the reader met it.
    #[test]
    fn a_read_failure_is_kept_for_the_parser() {
        let dir = std::env::temp_dir().join(format!("get_mnv_vcf_text_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("broken.vcf");
        let mut bytes =
            b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n".to_vec();
        bytes.extend_from_slice(b"c\t1\t.\tA\tG\t.\tPASS\t.\n");
        bytes.extend_from_slice(b"c\t2\t.\tA\tG\t.\tPASS\t\xff\n");
        bytes.extend_from_slice(b"c\t3\t.\tA\tG\t.\tPASS\t.\n");
        std::fs::write(&path, bytes).unwrap();
        let text = VcfText::read(path.to_str().unwrap()).unwrap().unwrap();
        assert_eq!(text.lines().len(), 3);
        let error = text.unread().expect("the bad line stops the reading");
        let direct = BufReader::new(std::fs::File::open(&path).unwrap())
            .lines()
            .find_map(Result::err)
            .unwrap();
        assert_eq!(error.to_string(), direct.to_string());
        assert_eq!(error.kind(), direct.kind());
        assert!(text.replay().last().unwrap().is_err());
        assert!(
            VcfText::read_holding_at_most(path.to_str().unwrap(), 64)
                .unwrap()
                .is_none(),
            "a file larger than the limit is left to be read again"
        );
        std::fs::remove_dir_all(&dir).ok();
    }
}
