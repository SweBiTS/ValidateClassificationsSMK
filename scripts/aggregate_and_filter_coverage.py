import sys
import logging
import pandas as pd
from pathlib import Path

# --- Snakemake objects ---
try:
    snk_input = snakemake.input
    snk_output = snakemake.output
    snk_params = snakemake.params
    snk_log = snakemake.log
except AttributeError as e:
    raise AttributeError(f"Missing expected Snakemake object attribute: {e}")

# --- File Paths & Params ---
# snk_input is a list-like object containing coverage file paths
coverage_files = snk_input.coverage_files
output_passed_list_tsv = snk_output.passed_samples_tsv
log_file = snk_log.path
min_coverage_threshold = float(snk_params.cov_threshold)

# --- Setup Logging ---
Path(log_file).parent.mkdir(parents=True, exist_ok=True)
logging.basicConfig(filename=log_file, filemode='w', level=logging.INFO,
                    format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

logging.info("--- Snakemake Parameters ---")
logging.info(f"Number of Input Coverage Files: {len(coverage_files)}")
logging.info(f"Output Passing Samples TSV: {output_passed_list_tsv}")
logging.info(f"Log File: {log_file}")
logging.info(f"Min Coverage Threshold: {min_coverage_threshold}")
logging.info("----------------------------")

# --- Process input files, filter, collect passing samples ---
passing_samples = [] # List to store dicts {'TaxID': tid, 'GenomeBasename': gbase}
processed_files = 0
error_files = 0

# --- Process input files, filter, collect passing samples ---
passing_samples = [] # List to store dicts {'TaxID': tid, 'GenomeBasename': gbase}
processed_files = 0
logging.info("Processing input coverage files...")
for cov_file in coverage_files:
    processed_files += 1
    logging.debug(f"Reading {cov_file}")
    try:
        with open(cov_file, 'r') as f:
            line = f.readline().strip()
            if not line:
                logging.error(f"FATAL: Coverage file contains only whitespace or is empty: {cov_file}")
                sys.exit(1)

            parts = line.split('\t')
            if len(parts) != 3:
                logging.error(f"FATAL: Malformed line in {cov_file} (expected 3 tab-separated fields, found {len(parts)}): '{line}'")
                sys.exit(1)

            tax_id = parts[0]
            genome_basename = parts[1]
            coverage_value_str = parts[2]
            try:
                coverage_value = float(coverage_value_str)
            except ValueError:
                logging.error(f"FATAL: Could not convert coverage value in {cov_file} ('{coverage_value_str}') to float.")
                sys.exit(1)

            # Check threshold
            if coverage_value >= min_coverage_threshold:
                logging.info(f"PASS: {tax_id}/{genome_basename} (Coverage={coverage_value:.4f})")
                passing_samples.append({'TaxID': tax_id, 'GenomeBasename': genome_basename})
            else:
                logging.info(f"FAIL: {tax_id}/{genome_basename} (Coverage={coverage_value:.4f})")

            # Check for unexpected extra lines (should only be one line)
            extra_line = f.readline()
            if extra_line:
                logging.error(f"FATAL: Expected only one line in {cov_file} but found more. Content: '{line}' followed by '{extra_line.strip()}'...")
                sys.exit(1)

    except FileNotFoundError:
            logging.error(f"FATAL: Input coverage file not found: {cov_file}")
            sys.exit(1)
    except Exception as e:
            logging.error(f"FATAL: Error processing file {cov_file}: {e}", exc_info=True)
            sys.exit(1)

logging.info(f"Finished processing {processed_files} input coverage files successfully.")
logging.info(f"Found {len(passing_samples)} samples passing threshold.")

# --- Write Output TSV ---
logging.info(f"Writing passing samples to {output_passed_list_tsv}")
try:
    if passing_samples:
        summary_df = pd.DataFrame(passing_samples)
        summary_df = summary_df[["TaxID", "GenomeBasename"]]
        Path(output_passed_list_tsv).parent.mkdir(parents=True, exist_ok=True)
        summary_df.to_csv(output_passed_list_tsv, sep='\t', index=False, header=True)
    else:
        logging.warning("No samples passed the coverage threshold. Creating empty output file with header.")
        Path(output_passed_list_tsv).parent.mkdir(parents=True, exist_ok=True)
        with open(output_passed_list_tsv, 'w') as f_out:
             f_out.write("TaxID\tGenomeBasename\n")

    logging.info(f"Successfully wrote list of {len(passing_samples)} passing samples to {output_passed_list_tsv}")

except Exception as e:
    logging.error(f"Failed to write output TSV {output_passed_list_tsv}: {e}", exc_info=True)
    sys.exit(1)

# --- Final Exit ---
logging.info("Aggregation and filtering script finished successfully.")
sys.exit(0)