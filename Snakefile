from pathlib import Path
import logging
import sys

# --- Configuration ---
configfile: "config.yaml"

# --- Setup Logging ---
# Get log level string from config, default to INFO, convert to upper case
log_level_str = config.get("LOG_LEVEL", "INFO").upper()

# Convert string level name to the corresponding logging constant
log_level_int = logging.getLevelName(log_level_str)

# Validate the level
if not isinstance(log_level_int, int):
    print(f"Warning: Invalid LOG_LEVEL '{log_level_str}' provided in config.yaml. Defaulting to INFO.", file=sys.stderr)
    log_level_int = logging.INFO

# Configure root logger
logging.basicConfig(
    level=log_level_int,
    format='%(asctime)s %(levelname)-8s %(name)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    stream=sys.stderr)

# Log the effective level being used
logging.info(f"Logging level set to: {logging.getLevelName(logging.getLogger().getEffectiveLevel())}")


# --- Define Core Paths ---
# These Path objects are globally accessible in included files
INPUT_DIR_P = Path(config["INPUT_DIR"])
OUTPUT_DIR_P = Path(config["OUTPUT_DIR"])
LOG_DIR_P = Path(config["LOG_DIR"])
INDEX_DIR_BBMAP_P = OUTPUT_DIR_P / config["INDEX_DIR_BBMAP"]  # Special case
GENOMES_DIR_P = Path(config["GENOMES_DIR"])
ENVS_DIR_P = Path(config["ENVS_DIR"])

logging.info(f"Input directory: {INPUT_DIR_P}")
logging.info(f"Input directory: {INPUT_DIR_P}")
logging.info(f"Output directory: {OUTPUT_DIR_P}")
logging.info(f"Log directory: {LOG_DIR_P}")
logging.info(f"Index directory (BBMap): {INDEX_DIR_BBMAP_P}")
logging.info(f"Genomes directory: {GENOMES_DIR_P}")
logging.info(f"Conda environments directory: {ENVS_DIR_P}")

# --- Validate Kraken2 Path (CONDITIONAL on workflow mode in the main config) ---
# Only validate if validation is requested AND a specific binary dir is given
if config.get("RUN_VALIDATION", True) and config.get("KRAKEN2_BIN_DIR"):
    kraken2_bin_dir = Path(config["KRAKEN2_BIN_DIR"]).resolve()
    logger.info(f"RUN_VALIDATION is true and KRAKEN2_BIN_DIR specified. Validating provided Kraken 2 path: {kraken2_bin_dir}")
    if not kraken2_bin_dir.is_dir():
        sys.exit(f"ERROR: KRAKEN2_BIN_DIR is not a valid directory: {kraken2_bin_dir}")
    kraken2_exe_path = kraken2_bin_dir / "kraken2"
    k2mask_exe_path = kraken2_bin_dir / "k2mask"
    if not kraken2_exe_path.is_file() or not os.access(kraken2_exe_path, os.X_OK):
         sys.exit(f"ERROR: 'kraken2' not found or not executable in {kraken2_bin_dir}")
    if not k2mask_exe_path.is_file() or not os.access(k2mask_exe_path, os.X_OK):
         sys.exit(f"ERROR: 'k2mask' not found or not executable in {kraken2_bin_dir}")
    logger.info(f"Validated that Kraken 2 executables are found in the specified directory.")
elif config.get("RUN_VALIDATION", True):
    logger.info("RUN_VALIDATION is true, workflow will run mapping-based classification validation.")
    logger.info("KRAKEN2_BIN_DIR not set. Rules will rely on conda environment for Kraken 2 tools.")
else:
    logger.info("RUN_VALIDATION is false, workflow will only run the mapping pipeline.")

def get_all_rewritten_fastqs_r1(wildcards):
    """
    Gets list of all header-rewritten R1 FASTQ files. This is for merging all simulated 
    reads before classification so save compute resources (need only load the Kraken database once).
    """
    return get_all_target_outputs(mapping_spec_data, CLEANED_SIM_READS_R1_PATTERN)

def get_all_rewritten_fastqs_r2(wildcards):
    """Dito for R2."""
    return get_all_target_outputs(mapping_spec_data, CLEANED_SIM_READS_R2_PATTERN)

def get_all_kraken_outs(wildcards):
    """
    Gets list of all expected per-sample kraken output files (after the 
    classification of aggregated R1 and R2 files).
    """
    return get_all_target_outputs(mapping_spec_data, KRAKEN_OUT_SIM_PATTERN)


# --- Include Rule Modules ---
include: "rules/map.smk"
include: "rules/simulation_validation.smk"

# --- Define Conditional Targets for Rule All ---
def get_final_targets():
    # Always include the main mapping output indices (map.smk)
    targets = get_all_target_indices()

    # Conditionally add the single final validation report file
    if config.get("RUN_VALIDATION", True):
        logger.info("Adding final validation summary reports to rule all targets.")
        final_validation_outputs = get_all_target_outputs(
            mapping_spec_data, PER_SAMPLE_SUMMARY_PATTERN
        )
        targets.extend(final_validation_outputs)
    else:
         logger.info("Skipping validation pipeline targets for rule all.")

    return targets

# --- Rule All: Defines the final targets of the workflow ---
localrules: all

rule all:
    input:
        get_final_targets()