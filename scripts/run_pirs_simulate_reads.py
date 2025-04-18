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
final_output_r1 = Path(snk_output.r1)
final_output_r2 = Path(snk_output.r2)
output_dir = Path(final_output_r1).parent
output_prefix = output_dir / "pirs"
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
logging.info(f"Output Prefix for pIRS: {output_prefix}")
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
    
    # Convert to integers
    mean_ins = int(round(mean_ins))
    stdev_ins = int(round(stdev_ins))
    avg_read_len = int(round(avg_read_len))
    logging.info(f"Found Data Params: ReadLen={avg_read_len}, InsertSize={mean_ins},{stdev_ins}")

except Exception as e:
    logging.error(f"Error reading or processing parameters from {params_json}: {e}", exc_info=True)
    sys.exit(1)

# pIRS will produce output fastq-files that are formatted like so:
# {output_prefix}_{read_len}_{insert_size}_1.fq.gz
expected_pirs_output_r1 = Path(f'{output_prefix}_{avg_read_len}_{mean_ins}_1.fq.gz')
expected_pirs_output_r2 = Path(f'{output_prefix}_{avg_read_len}_{mean_ins}_2.fq.gz')

# --- Construct pIRS Command ---
# No substitution errors, no indel errors, no GC bias
# Will produce output fastq-files that are formatted like so:
# {output_prefix}_{read_len}_{insert_size}_1.fq.gz
cmd_pirs = [
    "/usr/bin/time", "-v",
    "pirs",
    "simulate",
    "--read-len", str(avg_read_len),
    "--coverage", str(coverage),
    "--insert-len-mean", str(mean_ins),
    "--insert-len-sd", str(stdev_ins),
    "--no-substitution-errors",
    "--no-indel-errors",
    "--no-gc-bias",
    "--output-prefix", str(output_prefix),
    "--compress",
    "--no-logs",
    "--threads", str(snk_threads),
    str(ref_fasta)
]

# --- Run pIRS ---
logging.info(f"Running command: {' '.join(cmd_pirs)}")
try:
    with open(log_file, "a") as log_f:
        process = subprocess.run(cmd_pirs, check=True, text=True, encoding='utf-8',
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        log_f.write(f"\n--- pIRS read-simulator stdout ---\n")
        log_f.write(process.stdout)
        log_f.write(f"\n--- pIRS read-simulator stderr ---\n")
        log_f.write(process.stderr)
    logging.info(f"pIRS read-simulator completed successfully.")

    # --- Verification and renaming of output files ---
    try:
        # Rename the files to the final output names expected by Snakemake
        expected_pirs_output_r1.rename(final_output_r1)
        expected_pirs_output_r2.rename(final_output_r2)
    except FileNotFoundError as e:
        logging.error(f"Couldn't find expected raw pIRS output file: {e}")
        sys.exit(1)
    except FileExistsError as e:
        # Don't overwrite existing files
        logging.error(f"Error: Output file already exists, exiting to not overwrite: {e}")
        sys.exit(1)
    
    logging.info(f"Successfully renamed '{expected_pirs_output_r1}' to '{final_output_r1}'")
    logging.info(f"Successfully renamed '{expected_pirs_output_r2}' to '{final_output_r2}'")

except subprocess.CalledProcessError as e:
    logging.error(f"pIRS read-simulator failed with exit code {e.returncode}")
    if e.stdout: logging.error(f"pIRS stdout:\n{e.stdout}")
    if e.stderr: logging.error(f"pIRS stderr:\n{e.stderr}")
    sys.exit(1)
except Exception as e:
    logging.error(f"An unexpected error occurred running pIRS: {e}", exc_info=True)
    sys.exit(1)

logging.info("Read simulation script finished successfully.")