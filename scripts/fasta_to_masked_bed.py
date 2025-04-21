import subprocess
import sys
import logging
import re
from pathlib import Path
import gzip

# --- Snakemake objects ---
try:
    snk_input = snakemake.input
    snk_output = snakemake.output
    snk_params = snakemake.params
    snk_log = snakemake.log
    snk_wildcards = snakemake.wildcards
except AttributeError as e:
    raise AttributeError(f"Missing expected Snakemake object attribute: {e}")

# --- File Paths ---
input_fasta = snk_input.ref_fasta
output_bed = snk_output.masked_bed
temp_masked_fasta = snk_output.masked_fasta # masker output goes here
log_file = snk_log.path
masker_executable = snk_params.masker_exec

# --- Setup Logging ---
Path(log_file).parent.mkdir(parents=True, exist_ok=True)
logging.basicConfig(filename=log_file, filemode='w', level=logging.INFO,
                    format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

logging.info("--- Snakemake Parameters ---")
logging.info(f"Input FASTA: {input_fasta}")
logging.info(f"Output BED: {output_bed}")
logging.info(f"Temp Masked FASTA: {temp_masked_fasta}")
logging.info(f"Log File: {log_file}")
logging.info(f"Masker Executable: {masker_executable}")
logging.info(f"Wildcards: {snk_wildcards}")
logging.info("----------------------------")

# --- Step 1: Run masker to generate lowercase masked fasta ---
logging.info(f"Running masker on {input_fasta}...")
cmd_mask = [
    masker_executable,
    "-in", str(input_fasta),
    "-outfmt", "fasta",
    "-out", str(temp_masked_fasta)
]
try:
    # Ensure output directory exists for temp file
    Path(temp_masked_fasta).parent.mkdir(parents=True, exist_ok=True)
    with open(log_file, "a") as log_f:
        process = subprocess.run(cmd_mask, check=True, text=True, encoding='utf-8',
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        log_f.write("\n--- masker stdout ---\n")
        log_f.write(process.stdout)
        log_f.write("\n--- masker stderr ---\n")
        log_f.write(process.stderr)
    logging.info("masker completed successfully.")

except subprocess.CalledProcessError as e:
    logging.error(f"masker failed with exit code {e.returncode}")
    if e.stdout: logging.error(f"masker stdout:\n{e.stdout}")
    if e.stderr: logging.error(f"masker stderr:\n{e.stderr}")
    sys.exit(1)
except FileNotFoundError:
    logging.error(f"Error: '{masker_executable}' command not found. Check KRAKEN2_BIN_DIR or PATH.")
    sys.exit(1)
except Exception as e:
    logging.error(f"An unexpected error occurred running masker: {e}", exc_info=True)
    sys.exit(1)


# --- Step 2: Parse lowercase masked fasta to BED ---
logging.info(f"Parsing {temp_masked_fasta} to generate BED file {output_bed}")
lc_pattern = re.compile(r"[a-z]+") # Find runs of lowercase letters (masker masks by lowercasing)

try:
    # Ensure output directory exists for BED file
    Path(output_bed).parent.mkdir(parents=True, exist_ok=True)

    with open(temp_masked_fasta, "r", encoding='utf-8') as fin, \
         open(output_bed, "w", encoding='utf-8') as fout:

        current_contig = None
        current_pos_on_contig = 0 # Use 0-based coordinates internally

        for line in fin:
            line = line.strip()
            if not line: continue # Skip empty lines

            if line.startswith(">"):
                # Extract contig name (part after > up to first space); ">Contig_1 some information" -> "Contig_1"
                current_contig = line.split(None, 1)[0][1:]
                current_pos_on_contig = 0
                logging.debug(f"Processing contig: {current_contig}")
            
            elif current_contig:
                seq_len = len(line)
                # Find all lowercase stretches within this line
                for match in lc_pattern.finditer(line):
                    # Get match start/end relative to line start
                    match_start_on_line = match.start()
                    match_end_on_line = match.end()

                    # Convert to 0-based genomic coordinates
                    bed_start = current_pos_on_contig + match_start_on_line
                    bed_end = current_pos_on_contig + match_end_on_line

                    # Write BED record: chrom, start (0-based), end (exclusive)
                    fout.write(f"{current_contig}\t{bed_start}\t{bed_end}\n")

                # Update genomic position for next line
                current_pos_on_contig += seq_len
            else:
                # Sequence line before first header? Ignore? Log warning?
                logging.warning(f"Skipping sequence line before first header: {line[:50]}...")

    logging.info(f"Finished parsing. Masked regions saved to {output_bed}")

except FileNotFoundError:
    logging.error(f"Masked fasta file not found for parsing: {temp_masked_fasta}")
    sys.exit(1)
except Exception as e:
    logging.error(f"An error occurred during parsing or writing BED: {e}", exc_info=True)
    sys.exit(1)

logging.info("Masked region identification script finished successfully.")
