import subprocess
import sys
import logging
import json
import re
from pathlib import Path

# --- Snakemake inputs/outputs ---
# Access attributes directly - script will fail if they aren't defined in the rule
try:
    input_bam = snakemake.input.bam
    output_stats_raw = snakemake.output.stats_raw
    output_json = snakemake.output.params_json
    log_file = snakemake.log.path
    threads = snakemake.threads
    wildcards = snakemake.wildcards
except AttributeError as e:
    # If Snakemake object doesn't have an expected attribute
    raise AttributeError(f"Missing expected Snakemake input/output/log attribute: {e}")

# --- Setup Logging ---
# Ensure log directory exists before setting up logging
try:
    Path(log_file).parent.mkdir(parents=True, exist_ok=True)
except Exception as e:
    print(f"Error: Could not create log directory {Path(log_file).parent}. Exiting. Details: {e}", file=sys.stderr)
    sys.exit(1)

logging.basicConfig(filename=log_file, filemode='w', level=logging.INFO,
                    format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

# Log Snakemake parameters received (useful for debugging)
logging.info("--- Snakemake Parameters ---")
logging.info(f"Input BAM: {input_bam}")
logging.info(f"Output Raw Stats: {output_stats_raw}")
logging.info(f"Output JSON Params: {output_json}")
logging.info(f"Log File: {log_file}")
logging.info(f"Threads: {threads}")
logging.info(f"Wildcards: {wildcards}")
logging.info("----------------------------")
logging.info(f"Starting samtools stats for {input_bam}")

# --- Construct samtools command ---
cmd = [
    "samtools", "stats",
    f"--threads", str(threads),
    str(input_bam)
]

# --- Run samtools stats ---
logging.info(f"Running command: {' '.join(cmd)}")
try:
    # Run command and redirect stdout to the raw output file
    # Ensure output directory for stats_raw exists
    Path(output_stats_raw).parent.mkdir(parents=True, exist_ok=True)
    with open(output_stats_raw, "w", encoding='utf-8') as stats_f, open(log_file, "a") as log_f:
        process = subprocess.run(cmd, check=True, text=True, encoding='utf-8',
                                 stdout=stats_f, stderr=subprocess.PIPE)
        # Write stderr to log file
        log_f.write("\n--- samtools stats stderr ---\n")
        log_f.write(process.stderr)
    logging.info("samtools stats completed successfully.")

except subprocess.CalledProcessError as e:
    logging.error(f"samtools stats failed with exit code {e.returncode}")
    if e.stderr: logging.error(f"samtools stats stderr:\n{e.stderr}")
    sys.exit(1)
except FileNotFoundError:
    logging.error("Error: 'samtools' command not found. Is it in the environment PATH?")
    sys.exit(1)
except Exception as e:
    logging.error(f"An unexpected error occurred running samtools stats: {e}", exc_info=True)
    sys.exit(1)


# --- Parse samtools stats output ---
logging.info(f"Parsing samtools stats output file: {output_stats_raw}")
extracted_params = {
    "mean_insert_size": None,
    "std_dev_insert_size": None,
    "average_read_length": None,
    "source_bam": str(input_bam),
    "source_stats_file": str(output_stats_raw)
}
found_count = 0
expected_count = 3 # How many params we are looking for

try:
    with open(output_stats_raw, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith("SN"):
                # Use regex for more robust parsing of potential extra whitespace/comments
                match = re.match(r"SN\s+(.+?):\s+([\d\.\-eE]+)", line)
                if not match: continue

                key = match.group(1).strip()
                value_str = match.group(2).strip()

                try:
                    value = float(value_str) # Try converting captured numeric string to float
                except ValueError:
                    logging.warning(f"Could not convert value '{value_str}' to float for key '{key}'. Skipping.")
                    continue # Skip if value captured wasn't actually numeric

                if key == "insert size average":
                    extracted_params["mean_insert_size"] = value
                    found_count += 1
                elif key == "insert size standard deviation":
                    extracted_params["std_dev_insert_size"] = value
                    found_count += 1
                elif key == "average length": # Get average read length too
                    extracted_params["average_read_length"] = value
                    found_count +=1

            # Stop if we found all needed params
            if found_count >= expected_count:
                 break

    # Validation after parsing
    parsing_successful = True
    if extracted_params["mean_insert_size"] is None:
        logging.warning("Could not find 'insert size average' in stats output.")
        parsing_successful = False
    if extracted_params["std_dev_insert_size"] is None:
        if extracted_params["mean_insert_size"] is not None:
            logging.warning("Could not find 'insert size standard deviation' in stats output. Assigning 0.0 as mean was found.")
            extracted_params["std_dev_insert_size"] = 0.0 # Assign a default only if mean exists
        else:
             logging.warning("Could not find 'insert size standard deviation' in stats output.")
             parsing_successful = False # Mark failure if mean is also missing

    if extracted_params["average_read_length"] is None:
        logging.warning("Could not find 'average length' in stats output.")
        parsing_successful = False

    if not parsing_successful:
         logging.error("Failed to parse one or more required parameters from samtools stats output.")
         # Exit if essential parameters are missing
         if extracted_params["mean_insert_size"] is None or extracted_params["std_dev_insert_size"] is None:
             logging.critical("Exiting: Missing essential simulation parameters (Mean/StdDev Insert Size).")
             sys.exit(1)
         # Otherwise, just log error and continue to write partial JSON if needed


    logging.info(f"Extracted params: Mean Insert={extracted_params['mean_insert_size']}, "
                 f"StdDev Insert={extracted_params['std_dev_insert_size']}, "
                 f"Avg Length={extracted_params['average_read_length']}")

    # --- Create Output JSON ---
    logging.info(f"Writing simulation parameters to {output_json}")
    # Ensure output directory exists
    Path(output_json).parent.mkdir(parents=True, exist_ok=True)
    with open(output_json, 'w', encoding='utf-8') as jf:
        json.dump(extracted_params, jf, indent=4)

    logging.info("Script finished successfully.")

except FileNotFoundError:
     logging.error(f"Samtools stats file not found for parsing: {output_stats_raw}")
     sys.exit(1)
except Exception as e:
    logging.error(f"An error occurred during parsing or writing JSON: {e}", exc_info=True)
    sys.exit(1)