# scripts/run_neat_config_simulate_reads.py

import subprocess
import sys
import logging
import json
import yaml
from pathlib import Path

# --- Snakemake objects ---
try:
    snk_input = snakemake.input
    snk_output = snakemake.output
    snk_params = snakemake.params
    snk_log = snakemake.log
    snk_wildcards = snakemake.wildcards
    snk_threads = snakemake.threads
    snk_config = snakemake.config
except AttributeError as e:
    raise AttributeError(f"Missing expected Snakemake object attribute: {e}")

# --- Define File Paths ---
ref_fasta = snk_input.ref_fasta
params_json = snk_input.params_json
output_r1_gz = snk_output.r1
output_r2_gz = snk_output.r2
output_dir = Path(output_r1_gz).parent
output_prefix = output_dir / "sim"
neat_config_yaml_path = snk_output.neat_config_yaml
log_file = snk_log.path

# --- Get Config Parameters ---
coverage = snk_config.get("SIMULATION_COVERAGE", 10)

# --- Setup Logging ---
Path(log_file).parent.mkdir(parents=True, exist_ok=True)
logging.basicConfig(filename=log_file, filemode='w', level=logging.INFO,
                    format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

# Log parameters
logging.info("--- Snakemake Parameters ---")
logging.info(f"Reference FASTA: {ref_fasta}")
logging.info(f"Input Params JSON: {params_json}")
logging.info(f"Output Prefix for NEAT (derived): {output_prefix}")
logging.info(f"Expected R1: {output_r1_gz}")
logging.info(f"Expected R2: {output_r2_gz}")
logging.info(f"Log File: {log_file}")
logging.info(f"Threads: {snk_threads}")
logging.info(f"Wildcards: {snk_wildcards}")
logging.info(f"Coverage (from config): {coverage}")
logging.info("----------------------------")

# --- Read Parameters from Input JSON ---
logging.info(f"Reading parameters from {params_json}")
try:
    with open(params_json, 'r', encoding='utf-8') as f:
        sim_params = json.load(f)
    mean_ins = sim_params.get("mean_insert_size")
    stdev_ins = sim_params.get("std_dev_insert_size")
    avg_read_len = sim_params.get("average_read_length")
    if mean_ins is None or stdev_ins is None or avg_read_len is None:
        raise ValueError("Essential parameters missing from JSON.")
    read_len_mean = int(round(avg_read_len)) # Use integer read length
    logging.info(f"Found Data Params: ReadLen={read_len_mean}, InsertSize={mean_ins},{stdev_ins}")
except Exception as e:
    logging.error(f"Error reading or processing parameters from {params_json}: {e}", exc_info=True)
    sys.exit(1)

# --- Generate NEAT Config YAML ---
logging.info(f"Generating NEAT config file: {neat_config_yaml_path}")
neat_config_dict = {
    'reference': str(ref_fasta),
    'read_len': read_len_mean,
    'coverage': int(coverage),
    'paired_ended': True,
    'fragment_mean': float(mean_ins),
    'fragment_st_dev': float(stdev_ins),
    'avg_seq_error': 0,
    'mutation_rate': 0,
    'rescale_qualities': True,
    'threads': int(snk_threads)
    # Add 'rng_seed' here if reproducible simulation needed
}

try:
    Path(neat_config_yaml_path).parent.mkdir(parents=True, exist_ok=True)
    with open(neat_config_yaml_path, 'w', encoding='utf-8') as yf:
        yaml.dump(neat_config_dict, yf, default_flow_style=False)
except Exception as e:
    logging.error(f"Failed to write NEAT config YAML to {neat_config_yaml_path}: {e}", exc_info=True)
    sys.exit(1)

# --- Construct NEAT Command ---
cmd_neat = [
    "/usr/bin/time", "-v",
    "neat",
    "--log-name", f'{log_file}.neat',  # Log file for NEAT (for some reason NEAT outputs ending .log instead of .neat)
    "read-simulator",
    "-c", str(neat_config_yaml_path),
    "-o", str(output_prefix)
]

# --- Run NEAT ---
logging.info(f"Running command: {' '.join(cmd_neat)}")
try:
    with open(log_file, "a") as log_f:
        process = subprocess.run(cmd_neat, check=True, text=True, encoding='utf-8',
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        log_f.write(f"\n--- neat read-simulator stdout ---\n")
        log_f.write(process.stdout)
        log_f.write(f"\n--- neat read-simulator stderr ---\n")
        log_f.write(process.stderr)
    logging.info(f"neat read-simulator completed successfully.")

    # --- Verification ---
    # Check if the FINAL gzipped files Snakemake expects were created
    if not Path(output_r1_gz).is_file():
        logging.error(f"Final output R1 file not found: {output_r1_gz}")
        sys.exit(1)
    if not Path(output_r2_gz).is_file():
        logging.error(f"Final output R2 file not found: {output_r2_gz}")
        sys.exit(1)
    logging.info("Expected output files found.")

except subprocess.CalledProcessError as e:
    logging.error(f"neat read-simulator failed with exit code {e.returncode}")
    if e.stdout: logging.error(f"neat stdout:\n{e.stdout}")
    if e.stderr: logging.error(f"neat stderr:\n{e.stderr}")
    sys.exit(1)
except Exception as e:
    logging.error(f"An unexpected error occurred running NEAT: {e}", exc_info=True)
    sys.exit(1)

logging.info("Read simulation script finished successfully.")