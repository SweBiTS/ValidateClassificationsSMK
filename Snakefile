from pathlib import Path
import yaml
import logging
import sys

# TODO: For de-cluttering of the Snakefile, move all patterns and helper functions into their own files, and import them here.
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

# --- Setup Logging --- #
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

# --- Define Core Paths --- #
LOG_DIR = Path("logs")
OUTPUT_DIR = Path("output")
INPUT_DIR = Path("input")
BBMAP_INDEX_DIR = OUTPUT_DIR / "bbmap_indices"
GENOMES_DIR = Path("supporting_files/genomes")
ENVS_DIR = Path("envs")

# --- Pattern definitions --- #


# --- Validate Kraken2 Path (CONDITIONAL on workflow mode in the main config) --- #
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

# --- Environment Detection --- #
# Detect WSL vs native Linux. BBMap JNI acceleration (`usejni=t`) failed on WSL
# during development due to linking issues. Set flag to `usejni=f` (WSL) or 
# `usejni=t` (Linux) for compatibility/performance.
IS_WSL = 'WSL_DISTRO_NAME' in os.environ or 'WSL_INTEROP' in os.environ

if IS_WSL:
    logger.warning("WSL environment detected. Explicitly disabling BBMap JNI acceleration (usejni=f) due to known compatibility issues.")
    BBMAP_JNI_FLAG = "usejni=f"
else:
    logger.info("Non-WSL environment detected. Enabling BBMap JNI acceleration (usejni=t).")
    BBMAP_JNI_FLAG = "usejni=t"

# --- Load Mapping Specification --- #
config_mapping_spec_path = Path(config["MAPPING_SPECIFICATION"])
try:
    logger.info(f"Loading mapping specification from: {config_mapping_spec_path}")
    with open(config_mapping_spec_path, 'r') as f:
        raw_mapping_spec = yaml.safe_load(f)
        if raw_mapping_spec is None:
            mapping_spec_data = {}
            logger.warning(f"Mapping specification file '{config_mapping_spec_path}' is empty or invalid.")
        else:
            mapping_spec_data = {str(k): v for k, v in raw_mapping_spec.items()}
            logger.info(f"Loaded {len(mapping_spec_data)} taxID entries from spec file.")
except FileNotFoundError:
    logger.error(f"Mapping specification file not found at '{config_mapping_spec_path}'. Exiting.")
    sys.exit(1)
except yaml.YAMLError as e:
    logger.error(f"Error parsing YAML file '{config_mapping_spec_path}': {e}. Exiting.")
    sys.exit(1)
except Exception as e:
    logger.error(f"An unexpected error occurred reading '{config_mapping_spec_path}': {e}. Exiting.")
    sys.exit(1)


# ==================== #
# --- HELPER FUNCS --- #
# ==================== #

# --- Get Genome Basename --- #
def get_genome_basename(fasta_filename):
    if fasta_filename and isinstance(fasta_filename, str):
        return Path(fasta_filename).stem
    logger.warning(f"Invalid or non-string genome filename '{fasta_filename}' encountered in spec file.")
    return None

# --- Generate Target Output Files Function --- #
def get_all_target_outputs(mapping_spec, pattern_template):
    """
    Generates a list of output files based on the mapping specification and the
    provided file name pattern. Finds the wildcards from the mapping specification
    and applies them to the provided pattern template.
    """
    target_files = []
    logger.info(f"Generating target output file list using pattern: {pattern_template}")
    count = 0
    warning_count = 0
    for tax_id, genome_list in mapping_spec.items():
        if not isinstance(genome_list, list):
             logger.warning(f"Expected list of genomes for tax_id {tax_id}, found {type(genome_list)}. Skipping this entry.")
             warning_count += 1
             continue
        for genome_fasta in genome_list:
            genome_basename = get_genome_basename(genome_fasta)
            if genome_basename:
                target_path = Path(pattern_template.format(
                    tax_id=tax_id,
                    genome_basename=genome_basename))
                target_files.append(str(target_path))
                count += 1
            else:
                warning_count += 1
                pass
    logger.info(f"Generated {count} target output files.")
    if warning_count > 0:
        logger.warning(f"Encountered {warning_count} issues while processing mapping specification.")
    return target_files

# --- Define Functions for Rule All --- #
# Final output files from the mapping branch of the workflow (map.smk)
def get_all_target_mapping_outputs():
    return get_all_target_outputs(mapping_spec_data, DEDUP_BAM_INDEX_OUT_PATTERN)

# Final output files from the validation branch of the workflow (simulation_validation.smk)
def get_all_target_validation_outputs():
    return get_all_target_outputs(mapping_spec_data, PLOT_SCATTER_HTML_PATTERN)

# ================================== #
# --- WORKFLOWS AND TARGET FILES --- #
# ================================== #

# --- Include Rule Modules --- #
include: "rules/map.smk"
include: "rules/simulation_validation.smk"

# --- Define Conditional Targets for Rule All --- #
def get_final_targets():
    # Always include the main mapping output indices (map.smk)
    targets = get_all_target_mapping_outputs()

    # Conditionally add the single final validation report file
    if config.get("RUN_VALIDATION", True):
        logger.info("Adding final validation summary reports to rule all targets.")
        final_validation_outputs = get_all_target_validation_outputs()
        targets.extend(final_validation_outputs)
    else:
        logger.info("Skipping validation pipeline targets for rule all.")

    return targets

# --- Define the final targets of the workflow --- #
localrules: all

rule all:
    input:
        get_final_targets()