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
logging.info(f"Output directory: {OUTPUT_DIR_P}")
logging.info(f"Log directory: {LOG_DIR_P}")
logging.info(f"Index directory (BBMap): {INDEX_DIR_BBMAP_P}")
logging.info(f"Genomes directory: {GENOMES_DIR_P}")
logging.info(f"Conda environments directory: {ENVS_DIR_P}")

# --- Include Rule Modules ---
include: "rules/map.smk"

# --- Rule All: Defines the final targets of the workflow ---
localrules: all

rule all:
    input:
        # Call the function(s) defined in map.smk
        get_all_target_indices()
        # get_all_target_stats()