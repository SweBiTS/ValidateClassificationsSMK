import sys
import logging
import re
from pathlib import Path
from collections import defaultdict

# --- Snakemake objects ---
try:
    snk_input = snakemake.input
    snk_output = snakemake.output
    snk_params = snakemake.params
    snk_log = snakemake.log
except AttributeError as e:
    raise AttributeError(f"Missing expected Snakemake object attribute: {e}")

# --- File Paths & Params ---
combined_kraken_file = snk_input.combined_kraken_out
# snk_output is a list of all expected output file paths
expected_output_files = snk_output # Snakemake.OutputFiles object behaves like list
log_file = snk_log.path
header_regex_pattern = snk_params.header_regex

# --- Setup Logging ---
Path(log_file).parent.mkdir(parents=True, exist_ok=True)
logging.basicConfig(filename=log_file, filemode='w', level=logging.INFO,
                    format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

logging.info("--- Snakemake Parameters ---")
logging.info(f"Combined Kraken Input: {combined_kraken_file}")
logging.info(f"Number of Expected Output Files: {len(expected_output_files)}")
logging.info(f"Log File: {log_file}")
logging.info(f"Header Regex: {header_regex_pattern}")
logging.info("----------------------------")

# --- Pre-compile regex and prepare output file mapping ---
logging.info("Compiling regex and preparing output file map...")
output_handles = {} # Dictionary to hold open file handles { (taxid, genome): handle }
try:
    header_regex = re.compile(header_regex_pattern)
    # Create a mapping from (taxid, genome) tuple to output file path
    output_file_map = {}
    output_dirs = set()

    # Regex to extract wildcards from expected output paths
    # Assumes pattern ".../{tax_id}/{genome_basename}/sim.kraken.out"
    out_pattern_regex = re.compile(r".*/(\d+)/([a-zA-Z0-9_.\-]+)/sim\.kraken\.out$")

    for outfile_path_str in expected_output_files:
        outfile_path = Path(outfile_path_str)
        match = out_pattern_regex.match(outfile_path_str)
        if match:
            tax_id = match.group(1)
            genome_basename = match.group(2)
            output_file_map[(tax_id, genome_basename)] = outfile_path_str
            output_dirs.add(outfile_path.parent)
        else:
            # Treat as fatal error if expected path cannot be parsed
            logging.error(f"FATAL: Could not parse expected wildcards from required output path: {outfile_path_str}")
            logging.error("This likely indicates a mismatch between rule output patterns and the parsing regex in the script, or an issue with the get_all_kraken_outs function.")
            sys.exit(1)

    if not output_file_map:
         # If NO output files could be parsed, then something is wrong
         raise ValueError("Could not create mapping from wildcards to output files. Check patterns/helper function.")

    # Create all output directories beforehand
    for out_dir in output_dirs:
         out_dir.mkdir(parents=True, exist_ok=True)

    # Open all output files (consider resource limits if # files is huge)
    for key, filepath in output_file_map.items():
        output_handles[key] = open(filepath, "w", encoding='utf-8')

    logging.info(f"Prepared and opened {len(output_handles)} output file destinations.")

except Exception as e:
    logging.error(f"Error during setup: {e}", exc_info=True)
    # Close any files that were opened
    for handle in output_handles.values():
         handle.close()
    sys.exit(1)


# --- Process Combined Kraken File and Split ---
logging.info(f"Processing combined kraken file: {combined_kraken_file}")
processed_lines = 0
try:
    with open(combined_kraken_file, "r", encoding='utf-8') as fin:
        for line in fin:
            processed_lines += 1
            parts = line.strip().split('\t')

            # Kraken output format: C/U, ReadID, TaxID, Length, k-mer pairs
            if len(parts) >= 5:
                read_id = parts[1]
                
                # Try to parse embedded TaxID and GenomeBaseName
                match = header_regex.search(read_id)
                if match:
                    tax_id = match.group(1)
                    genome_basename = match.group(2)
                    key = (tax_id, genome_basename)

                    # Check if we have a corresponding output file handle
                    if key in output_handles:
                        output_handles[key].write(line) # Write original line
                    else:
                        # ERROR: Read ID parsed but doesn't match any expected output
                        logging.error(f"Read ID '{read_id}' contained TaxID/Genome '{key}' but no matching output file was expected based on rule definition. Stopping.")
                        sys.exit(1)
                else:
                    # ERROR: Cannot parse required info from read ID format
                    logging.error(f"Could not parse TaxID/Genome (using regex '{header_regex_pattern}') from read ID '{read_id}'. Stopping.")
                    sys.exit(1)
            else:
                # ERROR: Kraken output line format incorrect
                logging.error(f"Malformed Kraken output line (expected >=5 tab-separated fields, got {len(parts)}): {line.strip()}")
                sys.exit(1)

            # Progress indicator
            if processed_lines % 5000000 == 0:
                logging.info(f"Processed {processed_lines} lines...")

except FileNotFoundError:
    logging.error(f"Combined Kraken input file not found: {combined_kraken_file}")
    sys.exit(1)
except Exception as e:
    logging.error(f"Error processing Kraken file {combined_kraken_file} around line {processed_lines}: {e}", exc_info=True)
    sys.exit(1)
finally:
    # Ensure all output files are closed
    logging.info("Closing output file handles.")
    for handle in output_handles.values():
         handle.close()

logging.info(f"Finished splitting successfully. Processed Lines: {processed_lines}")
logging.info("Kraken output splitting script finished.")