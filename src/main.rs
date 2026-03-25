//! CLI entry point for get_MNV.

use env_logger::Env;
use get_mnv::cli;
use get_mnv::{error, pipeline};

fn main() {
    env_logger::Builder::from_env(Env::default().default_filter_or("info")).init();
    let args = cli::parse_args();
    if let Err(error) = pipeline::run(&args) {
        if let Some(error_json_path) = args.error_json.as_deref() {
            if let Err(write_err) = error::write_error_json(error_json_path, &error) {
                eprintln!(
                    "[E100] Failed to write error JSON '{}': {}",
                    error_json_path, write_err
                );
            }
        }
        eprintln!("{}", error);
        std::process::exit(error.code.exit_code());
    }
}
