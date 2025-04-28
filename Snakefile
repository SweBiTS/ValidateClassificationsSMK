from pathlib import Path
import pandas as pd
import yaml
import logging
import sys

# TODO: Implement dynamic resource allocation that depend on the size of the input files.
# TODO: Implement resubmission of jobs that fail, that increases the memory and/or time resources based on the "attempts" variable
# TODO: Things for final report:
    # - For each tax_id and genome_basename:
        # - Number of reads to start with
        # - Mapping information (check the stats file from BBMap):
            # - Percentage reads mapped
            # - Percentage unambiguous
            # - Percentage ambiguous
            # - Percentage unmapped
        # - Percentage duplicates (removed)
        # - Number of generated reads and coverage
        # - Percentage correctly classified (to tax_id in mapping spec)


# ============================== #
# --- Housekeeping and Setup --- #
# ============================== #

# --- Configuration file --- #
configfile: "config.yaml"

# --- Define Core Paths and Patterns--- #
LOG_DIR = Path("logs")
OUTPUT_DIR = Path("output")
INPUT_DIR = Path("input")
GENOMES_DIR = Path("supporting_files/genomes")
ENVS_DIR = Path("envs")
MAPPING_SPEC = Path(config["MAPPING_SPECIFICATION"])
STARTING_FASTQ_PATTERN = str(INPUT_DIR / config["FASTQ_FILE_PATTERN"])

# --- Validate Kraken2 Path (conditinal on workflow mode in the main config) --- #
if config.get("RUN_VALIDATION", True) and config.get("KRAKEN2_BIN_DIR"):
    kraken2_bin_dir = Path(config["KRAKEN2_BIN_DIR"]).resolve()
    if not kraken2_bin_dir.is_dir():
        sys.exit(f"ERROR: KRAKEN2_BIN_DIR is not a valid directory: {kraken2_bin_dir}")
    kraken2_exe_path = kraken2_bin_dir / "kraken2"
    if not kraken2_exe_path.is_file() or not os.access(kraken2_exe_path, os.X_OK):
         sys.exit(f"ERROR: 'kraken2' not found or not executable in {kraken2_bin_dir}")
    print(f"WORKFLOW_INFO: Will use Kraken 2 binary in directory {kraken2_bin_dir}")

# --- Environment Detection --- #
# Detect WSL vs native Linux. BBMap JNI acceleration (`usejni=t`) failed on WSL
# during development due to linking issues. Set flag to `usejni=f` (WSL) or 
# `usejni=t` (Linux) for compatibility/performance.
IS_WSL = 'WSL_DISTRO_NAME' in os.environ or 'WSL_INTEROP' in os.environ

if IS_WSL:
    print("WORKFLOW_INFO: WSL environment detected. Explicitly disabling BBMap JNI acceleration (usejni=f) due to known compatibility issues.")
    BBMAP_JNI_FLAG = "usejni=f"
else:
    print("WORKFLOW_INFO: Non-WSL environment detected. Enabling BBMap JNI acceleration (usejni=t).")
    BBMAP_JNI_FLAG = "usejni=t"

# ==================== #
# --- HELPER FUNCS --- #
# ==================== #

# --- Load Mapping Specification --- #
def load_mapping_spec(config_mapping_spec_path):
    """
    Load the mapping specification from the provided path.
    Expected columns: 'tax_id', 'name', and 'fasta_file'.
    Validates fasta filenames (must end in .fa/.fasta/.fna).
    Removes duplicate rows.
    Adds a 'genome_basename' column derived from 'fasta_file'.
    Returns a pandas DataFrame.
    """
    mapping_spec_file = Path(config_mapping_spec_path)
    required_columns = ['tax_id', 'name', 'fasta_file']
    allowed_suffixes = {'.fasta', '.fa', '.fna'}

    if not mapping_spec_file.is_file():
        sys.exit(f"WORKFLOW_ERROR: Mapping specification file not found at '{mapping_spec_file}'. Exiting.")

    try:
        mapping_df = pd.read_csv(
            mapping_spec_file,
            sep='\t',
            comment='#',
            dtype={'tax_id': str, 'name': str, 'fasta_file': str})

        # Validate that the mapping spec file is formatted correctly and is non-empty
        missing_cols = [col for col in required_columns if col not in mapping_df.columns]
        if missing_cols:
            sys.exit(f"WORKFLOW_ERROR: Mapping specification TSV '{mapping_spec_file}' is missing required columns: {missing_cols}. Exiting.")

        mapping_df.dropna(subset=required_columns, inplace=True)
        if mapping_df.empty:
            print(f"WORKFLOW_WARNING: Mapping specification file '{mapping_spec_file}' is empty or contains only comments/header.")
            return pd.DataFrame(columns=required_columns + ['genome_basename'])

        # Validate fasta_file suffixes
        invalid_files = []
        for fname in mapping_df['fasta_file']:
            if Path(fname).suffix.lower() not in allowed_suffixes:
                invalid_files.append(fname)
        if invalid_files:
            error_msg = (f"WORKFLOW_ERROR: Found rows in '{mapping_spec_file}' with invalid 'fasta_file' suffixes "
                         f"(must end in {', '.join(allowed_suffixes)}):\n" + "\n".join(invalid_files))
            sys.exit(error_msg)

        # Drop duplicated rows
        initial_rows = len(mapping_df)
        mapping_df.drop_duplicates(inplace=True)
        dropped_count = initial_rows - len(mapping_df)
        if dropped_count > 0:
            print(f"WORKFLOW_INFO: Removed {dropped_count} duplicate rows from mapping specification.")

        # Add the genome_basename column (remove suffix)
        mapping_df['genome_basename'] = mapping_df['fasta_file'].apply(lambda f: Path(f).stem)
        
        print(f"WORKFLOW_INFO: Found {len(mapping_df)} unique and valid entries from '{mapping_spec_file}'.")

    except Exception as e:
        sys.exit(f"WORKFLOW_ERROR: An error occurred reading or processing mapping specification TSV '{mapping_spec_file}': {e}. Exiting.")

    return mapping_df

# --- Generate Target Output Files Function --- #
def get_target_outputs(mapping_spec, pattern_template):
    """
    Generates a list of output files based on the mapping specification and the
    provided file name pattern. Finds the wildcards from the mapping specification
    and applies them to the provided pattern template.
    """
    target_files = []
    for tax_id, genome_list in mapping_spec.items():
        for genome_fasta in genome_list:
            genome_basename = Path(genome_fasta).stem
            target_path = Path(pattern_template.format(
                tax_id=tax_id,
                genome_basename=genome_basename))
            target_files.append(str(target_path))
    
    return target_files

# --- Target Function for Rule All --- #
def get_workflow_target():
    """
    Determines the final target file for the workflow based on the RUN_VALIDATION flag in the config.
    
    Target file is always at least the mapping branch target file.

    If validation:
        Add a flag file that triggers the coverage checkpoint to run, which determines the combinations of
        tax_id/genome_basename values to validate.
    """
    final_target = MAPPING_BRANCH_TARGET
    if config.get("RUN_VALIDATION", True):
        logger.info("RUN_VALIDATION=TRUE. Running the validation branch of the workflow.")
        final_target = [MAPPING_BRANCH_TARGET, VALIDATION_BRANCH_TARGET]
    return final_target


# ================================== #
# --- WORKFLOWS AND TARGET FILES --- #
# ================================== #

wildcard_constraints:
    tax_id="[0-9]+"

# --- Final Targets --- #
MAPPING_BRANCH_TARGET = str(OUTPUT_DIR / "1_mapping.done")
VALIDATION_BRANCH_TARGET = str(OUTPUT_DIR / "2_validation.done")

# --- Load the mapping specification --- #
MAPPING_SPEC_DATA = load_mapping_spec(MAPPING_SPEC)

# --- Include Rule Modules --- #
include: "rules/map.smk"
include: "rules/simulation_validation.smk"

localrules: all

# --- Define the final targets of the workflow --- #
rule all:
    input:
        # Depends on the RUN_VALIDATION flag in the config
        get_workflow_target()