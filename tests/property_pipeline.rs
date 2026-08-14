//! Pipeline-level invariants over generated inputs.
//!
//! The example-based tests in `integration.rs` and the scenario suite say what
//! the output should be for inputs somebody thought of. These say what has to
//! hold for *any* input, and they are here because the defects this branch fixed
//! were almost never a single wrong answer: they were two answers to the same
//! question. The same change written two ways, the same run's TSV and VCF, the
//! same gene given as a TSV or as a GFF, the same records in a different order.
//!
//! Each property runs the whole pipeline, so the case counts are deliberately
//! modest; they are a net, not a benchmark. A counterexample is written to
//! `tests/proptest-regressions/` and checked in, so it is replayed for ever after.

use clap::Parser;
use get_mnv::cli::Args;
use get_mnv::pipeline;
use proptest::prelude::*;
use std::collections::{BTreeSet, HashMap};
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

const GENE_START: usize = 101;
const GENE_END: usize = 400;
const CONTIG_LEN: usize = 600;

fn case_dir() -> PathBuf {
    static COUNTER: AtomicUsize = AtomicUsize::new(0);
    let id = COUNTER.fetch_add(1, Ordering::Relaxed);
    let dir = std::env::temp_dir().join(format!("get_mnv_prop_{}_{}", std::process::id(), id));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("create case dir");
    dir
}

/// Intron of the spliced gene model, in genomic coordinates.
const INTRON_START: usize = 251;
const INTRON_END: usize = 350;

fn write_reference(dir: &Path, sequence: &str) -> PathBuf {
    // A real intron opens with GT and closes with AG; without them there is no
    // splice site to call and the spliced model is untested.
    let mut bases: Vec<char> = sequence.chars().collect();
    bases[INTRON_START - 1] = 'G';
    bases[INTRON_START] = 'T';
    bases[INTRON_END - 2] = 'A';
    bases[INTRON_END - 1] = 'G';
    let sequence: String = bases.into_iter().collect();
    let path = dir.join("ref.fasta");
    fs::write(&path, format!(">chr1\n{sequence}\n")).expect("write reference");
    path
}

/// The sequence as written to disk, which is what the records must agree with.
fn stored_sequence(sequence: &str) -> Vec<char> {
    let mut bases: Vec<char> = sequence.chars().collect();
    bases[INTRON_START - 1] = 'G';
    bases[INTRON_START] = 'T';
    bases[INTRON_END - 2] = 'A';
    bases[INTRON_END - 1] = 'G';
    bases
}

/// The same gene said three ways: as a TSV annotation, as a single GFF CDS row,
/// and as that row carrying a Parent so it is read through the spliced-transcript
/// model. One CDS row is one segment, so the spliced model degenerates to the
/// same feature and all three must agree. The insertion boundary rule did not,
/// and the only difference between the inputs was the Parent attribute.
fn write_annotations(dir: &Path) -> (PathBuf, PathBuf) {
    let tsv = dir.join("genes.txt");
    fs::write(&tsv, format!("geneA\t{GENE_START}\t{GENE_END}\t+\n")).expect("write genes tsv");
    let gff = dir.join("genes.gff3");
    fs::write(
        &gff,
        format!(
            "##gff-version 3\nchr1\tsyn\tCDS\t{GENE_START}\t{GENE_END}\t.\t+\t0\tID=c1;Name=geneA\n"
        ),
    )
    .expect("write gff");
    let transcript = dir.join("transcript.gff3");
    fs::write(
        &transcript,
        format!(
            "##gff-version 3\n\
chr1\tsyn\tgene\t{GENE_START}\t{GENE_END}\t.\t+\t.\tID=g1;Name=geneA\n\
chr1\tsyn\tmRNA\t{GENE_START}\t{GENE_END}\t.\t+\t.\tID=m1;Parent=g1;Name=geneA\n\
chr1\tsyn\tCDS\t{GENE_START}\t{GENE_END}\t.\t+\t0\tID=c1;Parent=m1;Name=geneA\n"
        ),
    )
    .expect("write transcript gff");
    let _ = transcript;
    (tsv, gff)
}

/// The transcript-model spelling of the same single-exon gene.
fn transcript_annotation(dir: &Path) -> PathBuf {
    dir.join("transcript.gff3")
}

/// The same single-exon gene read on the minus strand, where the transcript runs
/// from higher coordinates down. Upstream and downstream swap there, and so does
/// the side of an anchor its inserted bases land on, which the boundary rules had
/// to be told twice.
fn write_minus_annotation(dir: &Path) -> PathBuf {
    let tsv = dir.join("minus.txt");
    fs::write(&tsv, format!("geneM\t{GENE_START}\t{GENE_END}\t-\n")).expect("write minus tsv");
    let path = dir.join("minus.gff3");
    fs::write(
        &path,
        format!(
            "##gff-version 3\n\
chr1\tsyn\tgene\t{GENE_START}\t{GENE_END}\t.\t-\t.\tID=g3;Name=geneM\n\
chr1\tsyn\tmRNA\t{GENE_START}\t{GENE_END}\t.\t-\t.\tID=m3;Parent=g3;Name=geneM\n\
chr1\tsyn\tCDS\t{GENE_START}\t{GENE_END}\t.\t-\t0\tID=c3;Parent=m3;Name=geneM\n"
        ),
    )
    .expect("write minus gff");
    path
}

/// A two-exon transcript over the same span, so the invariants also run against
/// a gene that has an intron, a splice donor and a splice acceptor. Every defect
/// keyed on a record's POS rather than on the bases it changes needed one to show.
fn write_spliced_annotation(dir: &Path) -> PathBuf {
    let path = dir.join("spliced.gff3");
    let exon1_end = INTRON_START - 1;
    let exon2_start = INTRON_END + 1;
    fs::write(
        &path,
        format!(
            "##gff-version 3\n\
chr1\tsyn\tgene\t{GENE_START}\t{GENE_END}\t.\t+\t.\tID=g2;Name=geneB\n\
chr1\tsyn\tmRNA\t{GENE_START}\t{GENE_END}\t.\t+\t.\tID=m2;Parent=g2;Name=geneB\n\
chr1\tsyn\tCDS\t{GENE_START}\t{exon1_end}\t.\t+\t0\tID=e1;Parent=m2;Name=geneB\n\
chr1\tsyn\tCDS\t{exon2_start}\t{GENE_END}\t.\t+\t0\tID=e2;Parent=m2;Name=geneB\n"
        ),
    )
    .expect("write spliced gff");
    path
}

fn write_vcf(path: &Path, records: &[String]) {
    let body = records.join("\n");
    fs::write(
        path,
        format!(
            "##fileformat=VCFv4.2\n\
##contig=<ID=chr1,length={CONTIG_LEN}>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{body}\n"
        ),
    )
    .expect("write vcf");
}

fn run_pipeline(dir: &Path, vcf: &Path, annotation: &Path, prefix: &str, both: bool) -> PathBuf {
    let reference = dir.join("ref.fasta");
    let mut argv = vec![
        "get_mnv".to_string(),
        "--vcf".to_string(),
        vcf.to_string_lossy().into_owned(),
        "--fasta".to_string(),
        reference.to_string_lossy().into_owned(),
    ];
    if annotation.extension().and_then(|e| e.to_str()) == Some("gff3") {
        argv.push("--gff".to_string());
    } else {
        argv.push("--genes".to_string());
    }
    argv.push(annotation.to_string_lossy().into_owned());
    if both {
        argv.push("--both".to_string());
    }

    let mut args = Args::parse_from(argv);
    args.output_dir = Some(dir.to_string_lossy().into_owned());
    args.output_prefix = Some(prefix.to_string());
    pipeline::run(&args).expect("pipeline should annotate generated input");
    dir.join(format!("{prefix}.MNV.tsv"))
}

fn read_rows(path: &Path) -> Vec<HashMap<String, String>> {
    let content = fs::read_to_string(path).expect("read TSV");
    let mut lines = content.lines();
    let header: Vec<String> = lines
        .next()
        .expect("TSV header")
        .split('\t')
        .map(str::to_string)
        .collect();
    lines
        .map(|line| {
            header
                .iter()
                .cloned()
                .zip(line.split('\t').map(str::to_string))
                .collect()
        })
        .collect()
}

/// What a row says, with nothing in it that depends on how the record was
/// written down: the bases that changed and the consequence drawn from them.
fn row_summary(row: &HashMap<String, String>) -> Vec<String> {
    [
        "Gene",
        "Positions",
        "Reference Bases",
        "Base Changes",
        "AA Changes",
        "Variant Type",
        "Change Type",
        "SO Term",
        "Impact",
        "Event Components",
        "HGVS g.",
    ]
    .iter()
    .map(|column| format!("{column}={}", row.get(*column).cloned().unwrap_or_default()))
    .collect()
}

fn summaries(path: &Path) -> BTreeSet<Vec<String>> {
    read_rows(path).iter().map(row_summary).collect()
}

/// What a run said about one gene, ignoring every other row it produced.
fn summaries_for(path: &Path, gene: &str) -> BTreeSet<Vec<String>> {
    read_rows(path)
        .iter()
        .filter(|row| row.get("Gene").map(String::as_str) == Some(gene))
        .map(row_summary)
        .collect()
}

/// The three gene models over one span, written together in one file.
///
/// Distinct names, so a row can be attributed to the model that produced it.
fn write_combined_annotation(dir: &Path) -> PathBuf {
    let path = dir.join("combined.gff3");
    let mut content = String::from("##gff-version 3\n");
    for source in ["genes.gff3", "spliced.gff3", "minus.gff3"] {
        let text = fs::read_to_string(dir.join(source)).expect("read annotation");
        for line in text.lines().filter(|line| !line.starts_with("##")) {
            content.push_str(line);
            content.push('\n');
        }
    }
    fs::write(&path, content).expect("write combined gff");
    path
}

fn other_base(base: char) -> char {
    match base {
        'A' => 'C',
        'C' => 'G',
        'G' => 'T',
        _ => 'A',
    }
}

/// Sites far enough apart that padding one cannot reach another, so each
/// generated record describes exactly one independent change. The range runs
/// from well before the gene to well after it, so a case can land inside the
/// coding span, outside every feature, or straddling the edge.
fn spaced_sites() -> impl Strategy<Value = Vec<usize>> {
    // The grid is offset by the gene's own start, so a slot lands exactly on the
    // first and last coding base as well as on either side of them. Every
    // boundary defect this branch fixed lived on one of those bases.
    proptest::collection::vec(0usize..44, 1..6).prop_map(|slots| {
        let mut sites: Vec<usize> = slots
            .into_iter()
            .map(|slot| GENE_START - 40 + slot * 10)
            .collect();
        sites.push(GENE_START);
        sites.push(GENE_END);
        // Codon boundaries hang off the gene's start on the plus strand and off
        // its end on the minus, so a grid anchored at one end never lands on the
        // other's. Both are needed: a record anchored exactly on a codon's edge is
        // where the boundary rules are decided.
        sites.push(GENE_END - 27);
        sites.push(GENE_END - 28);
        sites.sort_unstable();
        sites.dedup();
        sites
    })
}

/// Several changes inside one codon, or in two neighbouring ones.
///
/// The spaced grid above keeps records independent, which is what the padding
/// property needs and exactly what hides an interaction: almost every defect
/// this branch fixed needed two changes close enough to meet. A codon is three
/// bases, so "close enough" means within three.
fn clustered_sites() -> impl Strategy<Value = Vec<usize>> {
    (0usize..96, proptest::collection::vec(0usize..6, 2..4)).prop_map(|(codon, offsets)| {
        let codon_start = GENE_START + codon * 3;
        let mut sites: Vec<usize> = offsets
            .into_iter()
            .map(|offset| codon_start + offset)
            .collect();
        sites.sort_unstable();
        sites.dedup();
        sites
    })
}

/// What a generated record does at its site.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Change {
    Substitution,
    /// Two adjacent bases changed in one record, which is how a caller spells an
    /// MNP. Outside a gene this is the multi-position row whose second base the
    /// VCF writer used to drop.
    Mnp,
    Insertion,
    Deletion,
}

fn change_kinds(len: usize) -> impl Strategy<Value = Vec<Change>> {
    proptest::collection::vec(
        prop::sample::select(vec![
            Change::Substitution,
            Change::Substitution,
            Change::Mnp,
            Change::Insertion,
            Change::Deletion,
        ]),
        len,
    )
}

/// One VCF record for the change at `site`, spelled the plain way.
fn record_for(bases: &[char], site: usize, change: Change) -> String {
    let reference = bases[site - 1];
    match change {
        Change::Substitution => format!(
            "chr1\t{site}\t.\t{reference}\t{}\t100\tPASS\tDP=30",
            other_base(reference)
        ),
        Change::Mnp => {
            let reference_pair: String = bases[site - 1..site + 1].iter().collect();
            let alternate_pair: String = reference_pair.chars().map(other_base).collect();
            format!("chr1\t{site}\t.\t{reference_pair}\t{alternate_pair}\t100\tPASS\tDP=30")
        }
        Change::Insertion => format!(
            "chr1\t{site}\t.\t{reference}\t{reference}{}{}\t100\tPASS\tDP=30",
            other_base(reference),
            other_base(other_base(reference))
        ),
        Change::Deletion => {
            let deleted: String = bases[site - 1..site + 2].iter().collect();
            format!("chr1\t{site}\t.\t{deleted}\t{reference}\t100\tPASS\tDP=30")
        }
    }
}

fn records_for(bases: &[char], sites: &[usize], kinds: &[Change]) -> Vec<String> {
    sites
        .iter()
        .enumerate()
        .map(|(index, &site)| {
            record_for(
                bases,
                site,
                kinds.get(index).copied().unwrap_or(Change::Substitution),
            )
        })
        .collect()
}

fn reference_sequence() -> impl Strategy<Value = String> {
    proptest::collection::vec(prop::sample::select(vec!['A', 'C', 'G', 'T']), CONTIG_LEN)
        .prop_map(|bases| bases.into_iter().collect())
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 200, .. ProptestConfig::default() })]

    /// Padding a record with reference bases does not change what it says.
    ///
    /// `146 A>C` and `144 GTA>GTC` describe the same single substitution; the
    /// second merely starts two bases earlier. Every consequence the tool draws
    /// must come from the change, so the two runs must agree row for row. This is
    /// the invariant behind the splice-donor call that came out `intergenic` for
    /// the padded spelling of a change it called HIGH when written bare.
    #[test]
    fn padding_a_record_does_not_change_what_it_says(
        sequence in reference_sequence(),
        annotation_index in 0usize..5,
        sites in spaced_sites(),
        pad in 1usize..4,
    ) {
        let dir = case_dir();
        write_reference(&dir, &sequence);
        let (genes, _) = write_annotations(&dir);
        write_spliced_annotation(&dir);
        write_minus_annotation(&dir);
        let annotation = [
            genes,
            dir.join("genes.gff3"),
            dir.join("transcript.gff3"),
            dir.join("spliced.gff3"),
            dir.join("minus.gff3"),
        ][annotation_index]
            .clone();
        let bases = stored_sequence(&sequence);

        let mut bare = Vec::new();
        let mut padded = Vec::new();
        for &site in &sites {
            let reference = bases[site - 1];
            let alternate = other_base(reference);
            bare.push(format!(
                "chr1\t{site}\t.\t{reference}\t{alternate}\t100\tPASS\tDP=30"
            ));

            let start = site - pad;
            let ref_run: String = bases[start - 1..site].iter().collect();
            let mut alt_run: Vec<char> = ref_run.chars().collect();
            let last = alt_run.len() - 1;
            alt_run[last] = alternate;
            let alt_run: String = alt_run.into_iter().collect();
            padded.push(format!(
                "chr1\t{start}\t.\t{ref_run}\t{alt_run}\t100\tPASS\tDP=30"
            ));
        }

        let bare_vcf = dir.join("bare.vcf");
        let padded_vcf = dir.join("padded.vcf");
        write_vcf(&bare_vcf, &bare);
        write_vcf(&padded_vcf, &padded);

        let bare_tsv = run_pipeline(&dir, &bare_vcf, &annotation, "bare", false);
        let padded_tsv = run_pipeline(&dir, &padded_vcf, &annotation, "padded", false);

        prop_assert_eq!(
            summaries(&bare_tsv),
            summaries(&padded_tsv),
            "padding changed the annotation; sites {:?}, pad {}",
            sites,
            pad
        );
        let _ = fs::remove_dir_all(&dir);
    }

    /// The order of the records in the file does not change the output.
    ///
    /// Nothing about a variant depends on which line of the input it arrived on,
    /// so the tool sorts. An answer that moves when the input is shuffled is a
    /// dependence on arrival order that nobody asked for, and the frameshift
    /// linkage map had one: its key named an indel by coordinate alone, so which
    /// record was read last decided the other's answer.
    #[test]
    fn record_order_does_not_change_the_output(
        sequence in reference_sequence(),
        annotation_index in 0usize..5,
        (sites, kinds) in prop_oneof![spaced_sites(), clustered_sites()].prop_flat_map(|sites| {
            let len = sites.len();
            (Just(sites), change_kinds(len))
        }),
    ) {
        let dir = case_dir();
        write_reference(&dir, &sequence);
        let (_genes, _) = write_annotations(&dir);
        write_spliced_annotation(&dir);
        write_minus_annotation(&dir);
        let annotation = [
            dir.join("genes.txt"),
            dir.join("genes.gff3"),
            dir.join("transcript.gff3"),
            dir.join("spliced.gff3"),
            dir.join("minus.gff3"),
        ][annotation_index]
            .clone();
        let bases = stored_sequence(&sequence);

        let records = records_for(&bases, &sites, &kinds);
        let mut reversed = records.clone();
        reversed.reverse();

        let forward_vcf = dir.join("forward.vcf");
        let reversed_vcf = dir.join("reversed.vcf");
        write_vcf(&forward_vcf, &records);
        write_vcf(&reversed_vcf, &reversed);

        let forward_tsv = run_pipeline(&dir, &forward_vcf, &annotation, "forward", false);
        let reversed_tsv = run_pipeline(&dir, &reversed_vcf, &annotation, "reversed", false);

        prop_assert_eq!(
            fs::read_to_string(&forward_tsv).unwrap(),
            fs::read_to_string(&reversed_tsv).unwrap(),
            "reversing the records changed the output; sites {:?}",
            sites
        );
        let _ = fs::remove_dir_all(&dir);
    }

    /// One gene described two ways is one gene.
    ///
    /// `geneA 101 400 +` in a TSV annotation and a single GFF CDS row over the
    /// same span with phase 0 are the same feature, so they must produce the same
    /// rows. The insertion boundary rule differed between the two annotation
    /// models until this branch, which is exactly what this catches.
    #[test]
    fn a_gene_given_two_ways_annotates_the_same(
        sequence in reference_sequence(),
        (sites, kinds) in prop_oneof![spaced_sites(), clustered_sites()].prop_flat_map(|sites| {
            let len = sites.len();
            (Just(sites), change_kinds(len))
        }),
    ) {
        let dir = case_dir();
        write_reference(&dir, &sequence);
        let (genes, gff) = write_annotations(&dir);
        write_spliced_annotation(&dir);
        write_minus_annotation(&dir);
        let bases = stored_sequence(&sequence);

        let records = records_for(&bases, &sites, &kinds);
        let vcf = dir.join("variants.vcf");
        write_vcf(&vcf, &records);

        let from_tsv = run_pipeline(&dir, &vcf, &genes, "from_tsv", false);
        let from_gff = run_pipeline(&dir, &vcf, &gff, "from_gff", false);
        let from_transcript = run_pipeline(
            &dir,
            &vcf,
            &transcript_annotation(&dir),
            "from_transcript",
            false,
        );

        prop_assert_eq!(
            summaries(&from_tsv),
            summaries(&from_gff),
            "the TSV and GFF spellings disagree; sites {:?}",
            sites
        );
        prop_assert_eq!(
            summaries(&from_gff),
            summaries(&from_transcript),
            "adding a Parent changed the annotation; sites {:?}",
            sites
        );

        // And the same pair on the minus strand, where the transcript reads from
        // higher coordinates down. Upstream and downstream swap there, and so
        // does the side of an anchor its inserted bases land on, so the two
        // models had to be told the boundary rules twice and disagreed until they
        // were. A property that compares two spellings can only see a defect the
        // two do not share, which is exactly what this pair is for.
        let minus_tsv = run_pipeline(&dir, &vcf, &dir.join("minus.txt"), "minus_tsv", false);
        let minus_gff = run_pipeline(&dir, &vcf, &dir.join("minus.gff3"), "minus_gff", false);
        prop_assert_eq!(
            summaries(&minus_tsv),
            summaries(&minus_gff),
            "the two models disagree on the minus strand; sites {:?}",
            sites
        );
        let _ = fs::remove_dir_all(&dir);
    }

    /// The TSV and the VCF of one run describe the same changes.
    ///
    /// They are two renderings of one answer, so neither may carry a base the
    /// other leaves out. A multi-base substitution outside every gene reached the
    /// TSV and not the VCF, and an intergenic indel was deleted from the VCF by a
    /// threshold the TSV judged it to pass; both were two outputs of a single run
    /// disagreeing about what had been called.
    #[test]
    fn the_tsv_and_the_vcf_agree_on_what_changed(
        sequence in reference_sequence(),
        annotation_index in 0usize..5,
        (sites, kinds) in prop_oneof![spaced_sites(), clustered_sites()].prop_flat_map(|sites| {
            let len = sites.len();
            (Just(sites), change_kinds(len))
        }),
    ) {
        let dir = case_dir();
        write_reference(&dir, &sequence);
        let (genes, _) = write_annotations(&dir);
        write_spliced_annotation(&dir);
        write_minus_annotation(&dir);
        let annotation = [
            genes,
            dir.join("genes.gff3"),
            dir.join("transcript.gff3"),
            dir.join("spliced.gff3"),
            dir.join("minus.gff3"),
        ][annotation_index]
            .clone();
        let bases = stored_sequence(&sequence);

        let records = records_for(&bases, &sites, &kinds);
        let vcf = dir.join("variants.vcf");
        write_vcf(&vcf, &records);

        let tsv_path = run_pipeline(&dir, &vcf, &annotation, "both", true);
        let vcf_path = dir.join("both.MNV.vcf");

        let tsv_positions: BTreeSet<usize> = read_rows(&tsv_path)
            .iter()
            .flat_map(|row| {
                row.get("Positions")
                    .cloned()
                    .unwrap_or_default()
                    .split(", ")
                    .filter_map(|value| value.trim().parse::<usize>().ok())
                    .collect::<Vec<_>>()
            })
            .collect();
        let vcf_body = fs::read_to_string(&vcf_path).expect("read VCF");
        let records: Vec<(usize, usize)> = vcf_body
            .lines()
            .filter(|line| !line.starts_with('#') && !line.trim().is_empty())
            .filter_map(|line| {
                let mut fields = line.split('\t');
                let start = fields.nth(1)?.parse::<usize>().ok()?;
                let reference = fields.nth(1)?;
                Some((start, reference.len()))
            })
            .collect();
        // A record may render several bases at once, `301 AA>CC` for an MNP, so
        // the VCF covers a base when some record's REF spans it. What must hold
        // is that neither output carries a base the other leaves out.
        let covered: BTreeSet<usize> = records
            .iter()
            .flat_map(|&(start, len)| start..start + len.max(1))
            .collect();
        let record_starts: BTreeSet<usize> = records.iter().map(|&(start, _)| start).collect();

        let missing_from_vcf: Vec<usize> = tsv_positions
            .iter()
            .copied()
            .filter(|position| !covered.contains(position))
            .collect();
        prop_assert!(
            missing_from_vcf.is_empty(),
            "the TSV reports {:?} and the VCF does not carry them; sites {:?}",
            missing_from_vcf,
            sites
        );

        let invented_by_vcf: Vec<usize> = record_starts
            .iter()
            .copied()
            .filter(|position| !tsv_positions.contains(position))
            .collect();
        prop_assert!(
            invented_by_vcf.is_empty(),
            "the VCF has records at {:?} that no TSV row reports; sites {:?}",
            invented_by_vcf,
            sites
        );
        let _ = fs::remove_dir_all(&dir);
    }

    /// Annotating several genes at once says about each what naming it alone says.
    ///
    /// A gene's account of a base is its own. Three defects on this branch broke
    /// that: a gene lost its intron row once an overlapping gene annotated the
    /// base coding, two units sharing a Name suppressed one another, and a base a
    /// slippage join reads twice reached only one of its codons. Each made a run
    /// naming several genes report fewer rows for one of them than a run naming it
    /// alone, which is the shape this catches without having to think of the case.
    #[test]
    fn annotating_genes_together_says_what_each_says_alone(
        sequence in reference_sequence(),
        (sites, kinds) in prop_oneof![spaced_sites(), clustered_sites()].prop_flat_map(|sites| {
            let len = sites.len();
            (Just(sites), change_kinds(len))
        }),
    ) {
        let dir = case_dir();
        write_reference(&dir, &sequence);
        write_annotations(&dir);
        write_spliced_annotation(&dir);
        write_minus_annotation(&dir);
        let combined = write_combined_annotation(&dir);
        let bases = stored_sequence(&sequence);

        let records = records_for(&bases, &sites, &kinds);
        let vcf = dir.join("variants.vcf");
        write_vcf(&vcf, &records);

        let together = run_pipeline(&dir, &vcf, &combined, "together", false);

        for (gene, annotation) in [
            ("geneA", "genes.gff3"),
            ("geneB", "spliced.gff3"),
            ("geneM", "minus.gff3"),
        ] {
            let alone = run_pipeline(
                &dir,
                &vcf,
                &dir.join(annotation),
                &format!("alone_{gene}"),
                false,
            );
            prop_assert_eq!(
                summaries_for(&alone, gene),
                summaries_for(&together, gene),
                "{} says something different when it is not alone; sites {:?}",
                gene,
                sites
            );
        }
        let _ = fs::remove_dir_all(&dir);
    }
}
