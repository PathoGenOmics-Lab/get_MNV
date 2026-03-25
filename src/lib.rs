pub mod cli;
pub mod error;
pub mod io;
pub mod output;
pub mod pipeline;
pub mod read_count;
pub mod utils;
pub mod variants;

/// Core library version, sourced from Cargo.toml at compile time.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
