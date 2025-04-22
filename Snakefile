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
    k2mask_exe_path = kraken2_bin_dir / "k2mask"
    if not kraken2_exe_path.is_file() or not os.access(kraken2_exe_path, os.X_OK):
         sys.exit(f"ERROR: 'kraken2' not found or not executable in {kraken2_bin_dir}")
    if not k2mask_exe_path.is_file() or not os.access(k2mask_exe_path, os.X_OK):
         sys.exit(f"ERROR: 'k2mask' not found or not executable in {kraken2_bin_dir}")
    print(f"WORKFLOW_INFO: Will use Kraken 2 binaries in directory {kraken2_bin_dir}")

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
    Load the mapping specification from the provided path. The mapping specification
    is a YAML file that contains taxID entries and their corresponding genome FASTA files.
    """
    try:
        with open(config_mapping_spec_path, 'r') as f:
            raw_mapping_spec = yaml.safe_load(f)
            if raw_mapping_spec is None:
                mapping_spec_data = {}
                print(f"WORKFLOW_INFO: Mapping specification file '{config_mapping_spec_path}' is empty or invalid.")
            else:
                mapping_spec_data = {str(k): v for k, v in raw_mapping_spec.items()}
                print(f"WORKFLOW_INFO: Loaded {len(mapping_spec_data)} taxID entries from mapping specification file.")
    except FileNotFoundError:
        sys.exit(f"WORKFLOW_INFO: Mapping specification file not found at '{config_mapping_spec_path}'. Exiting.")
    except yaml.YAMLError as e:
        sys.exit(f"WORKFLOW_INFO: Error parsing YAML file '{config_mapping_spec_path}': {e}. Exiting.")
    except Exception as e:
        sys.exit(f"WORKFLOW_INFO: An unexpected error occurred reading '{config_mapping_spec_path}': {e}. Exiting.")
    
    return mapping_spec_data

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

# --- Final Targets --- #
MAPPING_BRANCH_TARGET = str(OUTPUT_DIR / "mapping.done")
VALIDATION_BRANCH_TARGET = str(OUTPUT_DIR / "validation.done")

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