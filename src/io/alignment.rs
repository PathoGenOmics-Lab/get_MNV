//! Indexed alignment input: BAM or CRAM, read through one interface.
//!
//! Every read-support path works on `bam::Record`, so a CRAM record is handed
//! over in that shape: decoded with the reference sequence, written with the
//! BAM record encoder and read back. That costs a copy per read, small next to
//! decoding CRAM itself, and keeps a single reading of CIGARs, bases and
//! qualities instead of a second one that could drift from the first.

use noodles::sam::alignment::io::Write as _;
use noodles::{bam, bgzf, cram, fasta, sam};
use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::path::Path;
use std::sync::{Mutex, OnceLock};

type BamReader = bam::io::IndexedReader<bgzf::io::Reader<File>>;

/// Whether a path names a CRAM file rather than a BAM.
pub fn is_cram_path(path: &str) -> bool {
    path.to_ascii_lowercase().ends_with(".cram")
}

/// A coordinate-sorted, indexed BAM (`.bai`/`.csi`) or CRAM (`.crai`).
pub enum AlignmentReader {
    Bam(Box<BamReader>),
    Cram(Box<cram::io::IndexedReader<File>>),
}

impl AlignmentReader {
    /// Opens `path` and reads its header. A CRAM is decoded against
    /// `reference_fasta`, the FASTA the run annotates with: its reads are
    /// stored as differences from that sequence and mean nothing without it.
    pub fn open(path: &str, reference_fasta: &str) -> io::Result<(Self, sam::Header)> {
        if is_cram_path(path) {
            let mut reader = cram::io::indexed_reader::Builder::default()
                .set_reference_sequence_repository(reference_repository(reference_fasta)?)
                .build_from_path(path)?;
            let header = reader.read_header()?;
            Ok((Self::Cram(Box::new(reader)), header))
        } else {
            let mut reader = bam::io::indexed_reader::Builder::default().build_from_path(path)?;
            let header = reader.read_header()?;
            Ok((Self::Bam(Box::new(reader)), header))
        }
    }

    /// The records overlapping `region`, read one at a time into a
    /// `bam::Record` with [`AlignmentQuery::read_record`].
    pub fn query<'r>(
        &'r mut self,
        header: &'r sam::Header,
        region: &noodles::core::Region,
    ) -> io::Result<AlignmentQuery<'r>> {
        match self {
            Self::Bam(reader) => Ok(AlignmentQuery::Bam(reader.query(header, region)?)),
            Self::Cram(reader) => Ok(AlignmentQuery::Cram(Box::new(CramQuery {
                query: reader.query(header, region)?,
                header,
                record: sam::alignment::RecordBuf::default(),
                encoded: Vec::new(),
            }))),
        }
    }
}

pub enum AlignmentQuery<'r> {
    Bam(bam::io::reader::Query<'r, bgzf::io::Reader<File>>),
    Cram(Box<CramQuery<'r>>),
}

pub struct CramQuery<'r> {
    query: cram::io::reader::Query<'r, 'r, 'r, File>,
    header: &'r sam::Header,
    record: sam::alignment::RecordBuf,
    encoded: Vec<u8>,
}

impl AlignmentQuery<'_> {
    /// Reads the next overlapping record into `record`, returning 0 at the end,
    /// as `bam::io::reader::Query::read_record` does.
    pub fn read_record(&mut self, record: &mut bam::Record) -> io::Result<usize> {
        match self {
            Self::Bam(query) => query.read_record(record),
            Self::Cram(cram) => {
                if cram.query.read_record_buf(&mut cram.record)? == 0 {
                    return Ok(0);
                }
                cram.encoded.clear();
                bam::io::Writer::from(&mut cram.encoded)
                    .write_alignment_record(cram.header, &cram.record)?;
                bam::io::Reader::from(cram.encoded.as_slice()).read_record(record)
            }
        }
    }
}

/// The reference a CRAM is decoded against, built once per FASTA and shared by
/// every reader that opens one: worker threads each open their own reader, and
/// loading the sequence once per thread would multiply its memory.
fn reference_repository(fasta_path: &str) -> io::Result<fasta::Repository> {
    static REPOSITORIES: OnceLock<Mutex<HashMap<String, fasta::Repository>>> = OnceLock::new();
    let cache = REPOSITORIES.get_or_init(|| Mutex::new(HashMap::new()));
    let mut cache = cache
        .lock()
        .map_err(|_| io::Error::other("reference cache lock poisoned"))?;
    if let Some(repository) = cache.get(fasta_path) {
        return Ok(repository.clone());
    }

    // With a `.fai` index the sequence is read on demand; without one the
    // whole FASTA is loaded, which is what the run does with it anyway.
    let index_path = format!("{fasta_path}.fai");
    let repository = if Path::new(&index_path).exists() {
        let reader = fasta::io::indexed_reader::Builder::default().build_from_path(fasta_path)?;
        fasta::Repository::new(fasta::repository::adapters::IndexedReader::new(reader))
    } else {
        let mut reader = fasta::io::reader::Builder.build_from_path(fasta_path)?;
        let records = reader.records().collect::<io::Result<Vec<_>>>()?;
        fasta::Repository::new(records)
    };
    cache.insert(fasta_path.to_string(), repository.clone());
    Ok(repository)
}
