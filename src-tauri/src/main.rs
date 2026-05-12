// Prevents additional console window on Windows in release, DO NOT REMOVE!!
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

mod commands;

use commands::{
    check_output_conflicts, detect_variant_input_format, ensure_fasta_index, get_bam_view,
    get_core_version, get_gff_features, read_tsv_file, resolve_output_paths, run_analysis,
    write_text_file,
};

fn main() {
    env_logger::init();

    tauri::Builder::default()
        .plugin(tauri_plugin_dialog::init())
        .plugin(tauri_plugin_shell::init())
        .invoke_handler(tauri::generate_handler![
            run_analysis,
            get_core_version,
            detect_variant_input_format,
            get_gff_features,
            ensure_fasta_index,
            check_output_conflicts,
            resolve_output_paths,
            read_tsv_file,
            get_bam_view,
            write_text_file,
        ])
        .run(tauri::generate_context!())
        .expect("error while running get_MNV application");
}
