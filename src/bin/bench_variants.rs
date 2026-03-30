use get_mnv::error::{AppError, AppResult};
use get_mnv::io::{self, Reference, VcfPosition};
use get_mnv::variants::{get_mnv_variants_for_gene, Gene, Strand};
use rayon::prelude::*;
use rayon::{ThreadPool, ThreadPoolBuilder};
use std::fs::OpenOptions;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::{Instant, SystemTime, UNIX_EPOCH};

#[derive(Debug, Clone)]
struct BenchConfig {
    warmup: usize,
    iterations: usize,
    dataset: Option<PathBuf>,
    contig: Option<String>,
    csv_out: Option<PathBuf>,
    max_avg_ms: Option<f64>,
    threads: Option<usize>,
    synthetic_scale: usize,
}

#[derive(Debug, Clone)]
struct BenchResult {
    mode: String,
    warmup: usize,
    iterations: usize,
    threads: usize,
    variants_per_iteration: usize,
    elapsed_ms: f64,
    avg_ms: f64,
}

fn parse_arg_usize(args: &[String], flag: &str, default: usize) -> usize {
    args.windows(2)
        .find(|pair| pair[0] == flag)
        .and_then(|pair| pair[1].parse::<usize>().ok())
        .unwrap_or(default)
        .max(1)
}

fn parse_arg_string(args: &[String], flag: &str) -> Option<String> {
    args.windows(2)
        .find(|pair| pair[0] == flag)
        .map(|pair| pair[1].clone())
}

fn parse_arg_optional_usize(args: &[String], flag: &str) -> Option<usize> {
    args.windows(2)
        .find(|pair| pair[0] == flag)
        .and_then(|pair| pair[1].parse::<usize>().ok())
}

fn parse_arg_f64(args: &[String], flag: &str) -> Option<f64> {
    args.windows(2)
        .find(|pair| pair[0] == flag)
        .and_then(|pair| pair[1].parse::<f64>().ok())
}

fn parse_config() -> BenchConfig {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    BenchConfig {
        warmup: parse_arg_usize(&args, "--warmup", 5),
        iterations: parse_arg_usize(&args, "--iters", 30),
        dataset: parse_arg_string(&args, "--dataset").map(PathBuf::from),
        contig: parse_arg_string(&args, "--contig"),
        csv_out: parse_arg_string(&args, "--csv").map(PathBuf::from),
        max_avg_ms: parse_arg_f64(&args, "--max-avg-ms"),
        threads: parse_arg_optional_usize(&args, "--threads"),
        synthetic_scale: parse_arg_usize(&args, "--synthetic-scale", 1),
    }
}

fn build_thread_pool(threads: Option<usize>) -> AppResult<Option<ThreadPool>> {
    match threads {
        Some(0) => Err(AppError::config("--threads must be >= 1")),
        Some(1) | None => Ok(None),
        Some(count) => ThreadPoolBuilder::new()
            .num_threads(count)
            .build()
            .map(Some)
            .map_err(|e| {
                AppError::config(format!("Failed to configure benchmark thread pool: {e}"))
            }),
    }
}

fn run_iterations<F>(iterations: usize, thread_pool: Option<&ThreadPool>, task: F) -> Vec<usize>
where
    F: Fn() -> usize + Sync,
{
    if let Some(pool) = thread_pool {
        pool.install(|| (0..iterations).into_par_iter().map(|_| task()).collect())
    } else {
        (0..iterations).map(|_| task()).collect()
    }
}

fn append_csv(path: &Path, result: &BenchResult) -> AppResult<()> {
    let existed = path.exists();
    let mut file = OpenOptions::new().create(true).append(true).open(path)?;
    if !existed {
        writeln!(
            file,
            "timestamp_unix,mode,warmup,iters,threads,variants_per_iteration,elapsed_ms,avg_ms"
        )?;
    }
    let timestamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_err(|e| AppError::msg(format!("System time error: {e}")))?
        .as_secs();
    writeln!(
        file,
        "{},{},{},{},{},{},{:.3},{:.3}",
        timestamp,
        result.mode,
        result.warmup,
        result.iterations,
        result.threads,
        result.variants_per_iteration,
        result.elapsed_ms,
        result.avg_ms
    )?;
    Ok(())
}

fn run_synthetic(
    config: &BenchConfig,
    thread_pool: Option<&ThreadPool>,
    threads: usize,
) -> BenchResult {
    let chrom = "chr1";
    let scale = config.synthetic_scale.max(1);
    let codon_repeats = 20_000usize.saturating_mul(scale);
    let gene_end = 60_000usize.saturating_mul(scale);
    let ref_sequence = "ATG".repeat(codon_repeats);
    let reference = Reference {
        sequence: &ref_sequence,
    };
    let gene = Gene {
        name: "gene1".to_string(),
        start: 1,
        end: gene_end,
        strand: Strand::Plus,
    };
    let snp_positions: Vec<VcfPosition> = (10..=(gene_end - 10))
        .step_by(10)
        .map(|pos| VcfPosition {
            position: pos,
            ref_allele: "A".to_string(),
            alt_allele: "G".to_string(),
            original_dp: Some(30),
            original_freq: Some(0.5),
            original_info: None,
        })
        .collect();

    let _ = run_iterations(config.warmup, thread_pool, || {
        get_mnv_variants_for_gene(&gene, &snp_positions, &reference, chrom).len()
    });

    let start = Instant::now();
    let produced_per_iteration = run_iterations(config.iterations, thread_pool, || {
        get_mnv_variants_for_gene(&gene, &snp_positions, &reference, chrom).len()
    });
    let elapsed = start.elapsed().as_secs_f64() * 1000.0;
    let produced = produced_per_iteration.last().copied().unwrap_or(0);

    BenchResult {
        mode: format!("synthetic:scale={scale}"),
        warmup: config.warmup,
        iterations: config.iterations,
        threads,
        variants_per_iteration: produced,
        elapsed_ms: elapsed,
        avg_ms: elapsed / config.iterations as f64,
    }
}

fn infer_annotation_path(dataset_dir: &Path) -> AppResult<PathBuf> {
    let gff = dataset_dir.join("MTB_ancestor.gff3");
    if gff.exists() {
        return Ok(gff);
    }
    let tsv = dataset_dir.join("anot_genes.txt");
    if tsv.exists() {
        return Ok(tsv);
    }
    Err(AppError::validation(format!(
        "Could not find annotation file in dataset '{}'. Expected MTB_ancestor.gff3 or anot_genes.txt",
        dataset_dir.display()
    )))
}

fn run_dataset(
    config: &BenchConfig,
    dataset_dir: &Path,
    thread_pool: Option<&ThreadPool>,
    threads: usize,
) -> AppResult<BenchResult> {
    let fasta_path = dataset_dir.join("MTB_ancestor.fas");
    let vcf_path = dataset_dir.join("G35894.var.snp.vcf");
    if !fasta_path.exists() {
        return Err(AppError::validation(format!(
            "Dataset FASTA not found: {}",
            fasta_path.display()
        )));
    }
    if !vcf_path.exists() {
        return Err(AppError::validation(format!(
            "Dataset VCF not found: {}",
            vcf_path.display()
        )));
    }
    let annotation_path = infer_annotation_path(dataset_dir)?;

    let references = io::load_references(fasta_path.to_string_lossy().as_ref())?;
    let snp_by_contig = io::load_vcf_positions_by_contig(
        vcf_path.to_string_lossy().as_ref(),
        None,
        false,
        false,
        false,
    )?;

    let target_contigs = if let Some(contig) = config.contig.as_ref() {
        vec![contig.clone()]
    } else {
        let mut contigs = snp_by_contig.keys().cloned().collect::<Vec<_>>();
        contigs.sort();
        contigs
    };
    if target_contigs.is_empty() {
        return Err(AppError::validation("No contigs available in dataset VCF"));
    }

    let mut contig_jobs = Vec::new();
    for contig in &target_contigs {
        let reference = io::reference_for_chrom(&references, contig)?;
        let snp_list = snp_by_contig
            .get(contig)
            .ok_or_else(|| {
                AppError::validation(format!("Dataset missing contig '{contig}' in VCF"))
            })?
            .clone();
        let default_features = vec!["gene".to_string(), "pseudogene".to_string()];
        let genes = io::load_genes(
            annotation_path.to_string_lossy().as_ref(),
            &snp_list,
            Some(contig),
            &default_features,
        )?;
        contig_jobs.push((contig.clone(), reference, snp_list, genes));
    }

    let run_one_iteration = || {
        contig_jobs
            .iter()
            .map(|(contig, reference, snp_list, genes)| {
                genes
                    .iter()
                    .map(|gene| get_mnv_variants_for_gene(gene, snp_list, reference, contig).len())
                    .sum::<usize>()
            })
            .sum::<usize>()
    };

    let _ = run_iterations(config.warmup, thread_pool, run_one_iteration);

    let start = Instant::now();
    let produced_per_iteration = run_iterations(config.iterations, thread_pool, run_one_iteration);
    let elapsed = start.elapsed().as_secs_f64() * 1000.0;
    let produced = produced_per_iteration.last().copied().unwrap_or(0);

    Ok(BenchResult {
        mode: format!("dataset:{}", dataset_dir.display()),
        warmup: config.warmup,
        iterations: config.iterations,
        threads,
        variants_per_iteration: produced,
        elapsed_ms: elapsed,
        avg_ms: elapsed / config.iterations as f64,
    })
}

fn print_result(result: &BenchResult) {
    println!("Benchmark mode: {}", result.mode);
    println!("Warmup iterations: {}", result.warmup);
    println!("Measured iterations: {}", result.iterations);
    println!("Threads: {}", result.threads);
    println!(
        "Variants produced per iteration: {}",
        result.variants_per_iteration
    );
    println!("Total measured time: {:.3} ms", result.elapsed_ms);
    println!("Average time per iteration: {:.3} ms", result.avg_ms);
}

fn run() -> AppResult<()> {
    let config = parse_config();
    let thread_pool = build_thread_pool(config.threads)?;
    let threads = thread_pool
        .as_ref()
        .map_or(1, rayon::ThreadPool::current_num_threads);

    let result = if let Some(dataset_dir) = config.dataset.as_ref() {
        run_dataset(&config, dataset_dir, thread_pool.as_ref(), threads)?
    } else {
        run_synthetic(&config, thread_pool.as_ref(), threads)
    };

    print_result(&result);
    if let Some(csv_path) = config.csv_out.as_ref() {
        append_csv(csv_path, &result)?;
        println!("CSV appended to {}", csv_path.display());
    }
    if let Some(max_avg_ms) = config.max_avg_ms {
        if result.avg_ms > max_avg_ms {
            return Err(AppError::validation(format!(
                "Benchmark regression: avg_ms {:.3} exceeded threshold {:.3}",
                result.avg_ms, max_avg_ms
            )));
        }
        println!(
            "Threshold check passed: avg_ms {:.3} <= {:.3}",
            result.avg_ms, max_avg_ms
        );
    }

    Ok(())
}

fn main() {
    if let Err(err) = run() {
        eprintln!("ERROR: {err}");
        std::process::exit(1);
    }
}
