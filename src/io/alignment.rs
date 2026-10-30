//! Indexed alignment input: BAM or CRAM, read through one interface.
//!
//! Every read-support path works on `bam::Record`, so a CRAM record is handed
//! over in that shape: decoded with the reference sequence, written with the
//! BAM record encoder and read back. That costs a copy per read, small next to
//! decoding CRAM itself, and keeps a single reading of CIGARs, bases and
//! qualities instead of a second one that could drift from the first.

use noodles::sam::alignment::io::Write as _;
use noodles::{bam, bgzf, cram, fasta, sam};
use std::collections::{HashMap, HashSet};
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
    Cram(Box<CramReader>),
}

/// A CRAM and its `.crai`, queried through the entries that matter.
///
/// noodles' own indexed query reads every container the index lists for the
/// contig and decodes it whole, once for each of its slices: a query for one
/// gene decoded the whole contig (60,000 reads took 132 s as CRAM and 3 s as
/// BAM), and a container of three slices returned every read three times.
pub struct CramReader {
    reader: cram::io::Reader<File>,
    index: cram::crai::Index,
    /// The index entries of the current query, kept here so the query can
    /// borrow them.
    selected: cram::crai::Index,
}

impl AlignmentReader {
    /// Opens `path` and reads its header. A CRAM is decoded against
    /// `reference_fasta`, the FASTA the run annotates with: its reads are
    /// stored as differences from that sequence and mean nothing without it.
    pub fn open(path: &str, reference_fasta: &str) -> io::Result<(Self, sam::Header)> {
        if is_cram_path(path) {
            let index = cram::crai::fs::read(format!("{path}.crai"))?;
            let mut reader = cram::io::reader::Builder::default()
                .set_reference_sequence_repository(reference_repository(reference_fasta)?)
                .build_from_path(path)?;
            let header = reader.read_header()?;
            let cram = CramReader {
                reader,
                index,
                selected: Vec::new(),
            };
            Ok((Self::Cram(Box::new(cram)), header))
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
            Self::Cram(cram) => {
                let CramReader {
                    reader,
                    index,
                    selected,
                } = &mut **cram;
                select_containers(index, header, region, selected);
                let selected: &'r cram::crai::Index = selected;
                Ok(AlignmentQuery::Cram(Box::new(CramQuery {
                    query: reader.query(header, selected, region)?,
                    contig: header.reference_sequences().get_index_of(region.name()),
                    header,
                    record: sam::alignment::RecordBuf::default(),
                    encoded: Vec::new(),
                })))
            }
        }
    }
}

/// The index entries a query for `region` has to read: those on its contig
/// whose span overlaps it, and one per container, since each slice of a
/// container has an entry of its own and every entry decodes the whole
/// container.
fn select_containers(
    index: &cram::crai::Index,
    header: &sam::Header,
    region: &noodles::core::Region,
    selected: &mut cram::crai::Index,
) {
    selected.clear();
    // A contig the header does not name selects nothing; the query itself
    // then fails with its own message.
    let Some(contig) = header.reference_sequences().get_index_of(region.name()) else {
        return;
    };
    let interval = region.interval();
    let mut containers = HashSet::new();
    for entry in index {
        if entry.reference_sequence_id() != Some(contig) {
            continue;
        }
        let overlaps = match entry.alignment_start() {
            Some(start) => {
                let end = start
                    .checked_add(entry.alignment_span().saturating_sub(1))
                    .unwrap_or(start);
                interval.intersects((start..=end).into())
            }
            // An entry without a start says nothing about where its reads
            // are, so it is read.
            None => true,
        };
        if overlaps && containers.insert(entry.offset()) {
            selected.push(entry.clone());
        }
    }
}

pub enum AlignmentQuery<'r> {
    Bam(bam::io::reader::Query<'r, bgzf::io::Reader<File>>),
    Cram(Box<CramQuery<'r>>),
}

pub struct CramQuery<'r> {
    query: cram::io::reader::Query<'r, 'r, 'r, File>,
    /// The contig queried, by its index in the header.
    contig: Option<usize>,
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
                // noodles decodes whole containers and keeps the reads that
                // overlap the query's positions, on any contig. A container
                // can hold the reads of several contigs, which samtools writes
                // for short ones such as plasmids, so a query on one contig
                // also counted another's reads at the same positions.
                loop {
                    if cram.query.read_record_buf(&mut cram.record)? == 0 {
                        return Ok(0);
                    }
                    if cram.record.reference_sequence_id() == cram.contig {
                        break;
                    }
                }
                // An unnamed read decodes with the name `*`, which the BAM
                // encoder refuses, so a CRAM of unnamed reads stopped the run
                // that the same reads as BAM completed. Without a name the
                // read is counted as unnamed, as from a BAM.
                if cram
                    .record
                    .name()
                    .is_some_and(|name| <[u8]>::eq(name.as_ref(), b"*"))
                {
                    *cram.record.name_mut() = None;
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

#[cfg(test)]
mod tests {
    use super::select_containers;
    use noodles::core::Position;
    use noodles::cram::crai;
    use noodles::sam::header::record::value::{map::ReferenceSequence, Map};
    use std::num::NonZeroUsize;

    fn entry(contig: usize, start: usize, span: usize, offset: u64) -> crai::Record {
        crai::Record::new(Some(contig), Position::new(start), span, offset, 0, 0)
    }

    fn offsets(region: &str, index: &crai::Index) -> Vec<u64> {
        let length = NonZeroUsize::new(1000).unwrap();
        let header = noodles::sam::Header::builder()
            .add_reference_sequence("chr1", Map::<ReferenceSequence>::new(length))
            .add_reference_sequence("chr2", Map::<ReferenceSequence>::new(length))
            .build();
        let mut selected = Vec::new();
        select_containers(index, &header, &region.parse().unwrap(), &mut selected);
        selected.iter().map(crai::Record::offset).collect()
    }

    /// Only the containers that overlap the query are read, and each of them
    /// once however many of its slices do.
    #[test]
    fn a_query_reads_each_overlapping_container_once() {
        let index = vec![
            entry(0, 1, 100, 10),
            entry(0, 101, 100, 10),
            entry(0, 300, 50, 500),
            entry(1, 1, 400, 900),
        ];
        assert_eq!(offsets("chr1:50-150", &index), vec![10]);
        assert_eq!(offsets("chr1:150-320", &index), vec![10, 500]);
        assert_eq!(offsets("chr1:200-200", &index), vec![10]);
        assert_eq!(offsets("chr1:201-299", &index), Vec::<u64>::new());
        assert_eq!(offsets("chr2:10-20", &index), vec![900]);
        assert_eq!(offsets("chr3:10-20", &index), Vec::<u64>::new());
    }
}
