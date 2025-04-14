from pathlib import Path
import yaml
import logging
import sys

# TODO: For de-cluttering of the Snakefile, move all patterns and helper functions into their own files, and import them here.
# TODO: Implement dynamic resource allocation that depend on the size of the input files.
# TODO: Implement resubmission of jobs that fail, that increases the memory and/or time resources based on the "attempts" variable

# ==================================================
#      Housekeeping and Setup
# ==================================================

# --- Configuration file ---
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

# --- PATTERN DEFINITIONS ---
# --- map.smk output patterns ---
RAW_BAM_OUT_PATTERN = str(OUTPUT_DIR_P / "mapping/{tax_id}/{genome_basename}/mapping.bam")
SORTED_BAM_OUT_PATTERN = str(OUTPUT_DIR_P / "mapping/{tax_id}/{genome_basename}/mapping.sorted.bam")
DEDUP_BAM_OUT_PATTERN = str(OUTPUT_DIR_P / "mapping/{tax_id}/{genome_basename}/mapping_dedup.bam")
DEDUP_BAM_INDEX_OUT_PATTERN = str(OUTPUT_DIR_P / "mapping/{tax_id}/{genome_basename}/mapping_dedup.bam.bai")
STATS_OUT_PATTERN = str(OUTPUT_DIR_P / "mapping/{tax_id}/{genome_basename}/mapping_stats.txt")

# --- map.smk log patterns ---
LOG_BBMAP_INDEX_PATTERN = str(LOG_DIR_P / "bbmap_index/{genome_basename}.log")
LOG_MAP_PATTERN = str(LOG_DIR_P / "bbmap_map/{tax_id}_{genome_basename}.log")
LOG_SORT_PATTERN = str(LOG_DIR_P / "samtools_sort/{tax_id}_{genome_basename}.log")
LOG_MARKDUP_PATTERN = str(LOG_DIR_P / "sambamba_markdup/{tax_id}_{genome_basename}.log")

# --- simulation_validation.smk output patterns ---
SAMTOOLS_STATS_OUT_PATTERN = str(OUTPUT_DIR_P / "validation/stats/{tax_id}/{genome_basename}/samtools_stats.txt")
SIM_PARAMS_OUT_PATTERN = str(OUTPUT_DIR_P / "validation/params/{tax_id}/{genome_basename}/sim_params.json")
SIM_READS_R1_PATTERN = str(OUTPUT_DIR_P / "validation/simulated_reads/{tax_id}/{genome_basename}/sim_r1.fastq.gz")
SIM_READS_R2_PATTERN = str(OUTPUT_DIR_P / "validation/simulated_reads/{tax_id}/{genome_basename}/sim_r2.fastq.gz")
NEAT_CONFIG_YAML_PATTERN = str(OUTPUT_DIR_P / "validation/simulated_reads/{tax_id}/{genome_basename}/neat_config.yml")
CLEANED_SIM_READS_R1_PATTERN = str(OUTPUT_DIR_P / "validation/simulated_reads/{tax_id}/{genome_basename}/sim_clean_R1.fq.gz")
CLEANED_SIM_READS_R2_PATTERN = str(OUTPUT_DIR_P / "validation/simulated_reads/{tax_id}/{genome_basename}/sim_clean_R2.fq.gz")
AGG_SIM_R1_PATTERN = str(OUTPUT_DIR_P / "validation/simulated_reads/aggregated/all_sim_clean_R1.fq.gz")
AGG_SIM_R2_PATTERN = str(OUTPUT_DIR_P / "validation/simulated_reads/aggregated/all_sim_clean_R2.fq.gz")
ALL_KRAKEN_OUT_PATTERN = str(OUTPUT_DIR_P / "validation/kraken2_sim/aggregated/all_sim.kraken.out")
ALL_KRAKEN_REPORT_PATTERN = str(OUTPUT_DIR_P / "validation/kraken2_sim/aggregated/all_sim.kraken.report")
KRAKEN_OUT_SIM_PATTERN = str(OUTPUT_DIR_P / "validation/kraken2_sim/{tax_id}/{genome_basename}/sim.kraken.out")
KRAKEN_REPORT_SIM_PATTERN = str(OUTPUT_DIR_P / "validation/kraken2_sim/{tax_id}/{genome_basename}/sim.kraken.report")
CORRECT_READ_IDS_PATTERN = str(OUTPUT_DIR_P / "validation/filtered_reads/{tax_id}/{genome_basename}/sim_correct_read_ids.txt")
FILTERED_SIM_READS_R1_PATTERN = str(OUTPUT_DIR_P / "validation/filtered_reads/{tax_id}/{genome_basename}/sim_correct_R1.fq.gz")
FILTERED_SIM_READS_R2_PATTERN = str(OUTPUT_DIR_P / "validation/filtered_reads/{tax_id}/{genome_basename}/sim_correct_R2.fq.gz")
MAPPED_SIM_BAM_PATTERN = str(OUTPUT_DIR_P / "validation/mapped_sim/{tax_id}/{genome_basename}/mapped_correct_sim.bam")
MAPPED_SIM_STATS_PATTERN = str(OUTPUT_DIR_P / "validation/mapped_sim/{tax_id}/{genome_basename}/mapped_correct_sim_stats.txt")
SORTED_MAPPED_SIM_BAM_PATTERN = str(OUTPUT_DIR_P / "validation/mapped_sim/{tax_id}/{genome_basename}/mapped_correct_sim.sorted.bam")
MASKED_REGIONS_BED_PATTERN = str(OUTPUT_DIR_P / "validation/masked_regions/{genome_basename}.masked.bed")
COVERED_REGIONS_BED_PATTERN = str(OUTPUT_DIR_P / "validation/coverage_sim/{tax_id}/{genome_basename}/covered_regions.bed")
MASKED_FASTA_PATTERN = str(OUTPUT_DIR_P / "validation/masked_regions/{genome_basename}.masked_lowercase.fasta")
CLASSIFIABLE_REGIONS_BED_PATTERN = str(OUTPUT_DIR_P / "validation/classifiable_regions/{tax_id}/{genome_basename}/classifiable_regions.bed")
OVERLAPPING_READS_BAM_PATTERN = str(OUTPUT_DIR_P / "validation/final_analysis/{tax_id}/{genome_basename}/overlapping_reads.bam")
PER_SAMPLE_SUMMARY_PATTERN = str(OUTPUT_DIR_P / "validation/summary_reports/{tax_id}/{genome_basename}/summary.tsv")
PLOT_SCATTER_HTML_PATTERN  = str(OUTPUT_DIR_P / "validation/plots/{tax_id}/{genome_basename}/scatter_reads_vs_length.html")

# --- simulation_validation.smk log patterns ---
LOG_ESTIMATE_PARAMS_PATTERN = str(LOG_DIR_P / "estimate_sim_params/{tax_id}/{genome_basename}.log")
LOG_SIMULATE_READS_PATTERN = str(LOG_DIR_P / "simulate_reads/{tax_id}/{genome_basename}.log")
LOG_REWRITE_HEADERS_PATTERN = str(LOG_DIR_P / "rewrite_headers/{tax_id}/{genome_basename}.log")
LOG_AGG_FASTQ_PATTERN = str(LOG_DIR_P / "aggregate_sim_fastqs/aggregate_sim_fastqs.log")
LOG_KRAKEN2_BATCH_PATTERN = str(LOG_DIR_P / "kraken2_sim/classify_all_simulated.log")
LOG_SPLIT_KRAKEN_PATTERN = str(LOG_DIR_P / "split_kraken/split_kraken_output.log")
LOG_EXTRACT_IDS_PATTERN = str(LOG_DIR_P / "extract_ids/{tax_id}/{genome_basename}.log")
LOG_FILTER_FASTQ_PATTERN = str(LOG_DIR_P / "filter_fastq/{tax_id}/{genome_basename}.log")
LOG_MAP_SIM_PATTERN = str(LOG_DIR_P / "map_simulated/{tax_id}/{genome_basename}.log")
LOG_SORT_SIM_PATTERN = str(LOG_DIR_P / "sort_simulated/{tax_id}/{genome_basename}.log")
LOG_CALC_SIM_COV_PATTERN = str(LOG_DIR_P / "calculate_sim_coverage/{tax_id}/{genome_basename}.log")
LOG_MASKED_REGIONS_PATTERN = str(LOG_DIR_P / "find_masked_regions/{genome_basename}.log")
LOG_DEFINE_CLASSIFIABLE_PATTERN = str(LOG_DIR_P / "define_classifiable_regions/{tax_id}/{genome_basename}.log")
LOG_INTERSECT_READS_PATTERN = str(LOG_DIR_P / "intersect_reads/{tax_id}/{genome_basename}.log")
LOG_SUMMARIZE_VALIDATION_PATTERN = str(LOG_DIR_P / "summarize_validation/{tax_id}/{genome_basename}.log")
LOG_PLOT_SCATTER_PATTERN = str(LOG_DIR_P / "plot_scatter/{tax_id}/{genome_basename}.log")

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

# --- Environment Detection ---
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

# --- Load Mapping Specification ---
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

# --- Helper Function: Get Genome Basename ---
def get_genome_basename(fasta_filename):
    if fasta_filename and isinstance(fasta_filename, str):
        return Path(fasta_filename).stem
    logger.warning(f"Invalid or non-string genome filename '{fasta_filename}' encountered in spec file.")
    return None

# --- Helper Function: Find Actual Genome File Path ---
def get_actual_fasta_path(wildcards):
    """Finds the full path to the genome file, checking common extensions."""
    basename = wildcards.genome_basename
    genomes_dir = GENOMES_DIR_P
    
    # Define the extensions to check in order of preference
    extensions_to_check = [".fasta", ".fa", ".fna"]

    for ext in extensions_to_check:
        potential_path = genomes_dir / f"{basename}{ext}"
        if potential_path.exists():
            logger.debug(f"Found genome file for basename '{basename}' at: {potential_path}")
            return str(potential_path) # Return the path string immediately if found

    # If the loop finishes without returning, no file was found
    raise FileNotFoundError(
        f"Genome file for basename '{basename}' not found with any of the extensions "
        f"{extensions_to_check} in directory '{genomes_dir}'")

# --- Generate Target Output Files Function ---
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

# --- Define Functions for Rule All ---
def get_all_target_indices():
    """The final output files from the mapping branch of the workflow (map.smk)."""
    return get_all_target_outputs(mapping_spec_data, DEDUP_BAM_INDEX_OUT_PATTERN)

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
    """Gets list of all expected per-sample kraken output files."""
    return get_all_target_outputs(mapping_spec_data, KRAKEN_OUT_SIM_PATTERN)

# --- Generate the list of expected Kraken output files (after splitting them) ---
# This is needed because I couldn't use a function in the output section of the rule
# But since we know what the output files will be called, we can generate a list of them
# here and use it in the rule.
ALL_EXPECTED_KRAKEN_OUTS = get_all_kraken_outs(None)
if not ALL_EXPECTED_KRAKEN_OUTS:
    logger.warning("No expected Kraken output files generated based on mapping spec. Check patterns/spec file.")

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
            mapping_spec_data, PLOT_SCATTER_HTML_PATTERN)
        targets.extend(final_validation_outputs)
    else:
        logger.info("Skipping validation pipeline targets for rule all.")

    return targets

# --- Rule All: Defines the final targets of the workflow ---
localrules: all

rule all:
    input:
        get_final_targets()