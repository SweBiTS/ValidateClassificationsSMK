# Snakefile

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
        get_all_target_bams()
        # get_all_target_stats()


# # Snakefile

# from pathlib import Path
# import logging
# import sys

# # --- Configuration ---
# configfile: "config.yaml"

# # --- Setup Logging ---
# # Get log level string from config, default to INFO, convert to upper case
# log_level_str = config.get("log_level", "INFO").upper()

# # Convert string level name (e.g., "INFO") to the corresponding logging constant (e.g., logging.INFO)
# log_level_int = logging.getLevelName(log_level_str)

# # Check if the conversion was successful (returns an int for valid names)
# # If not, default to INFO and issue a warning.
# if not isinstance(log_level_int, int):
#     print(f"Warning: Invalid log_level '{log_level_str}' provided in config.yaml. Defaulting to INFO.", file=sys.stderr)
#     log_level_int = logging.INFO

# # Configure root logger
# logging.basicConfig(
#     level=log_level_int,
#     format='%(asctime)s %(levelname)-8s %(name)s: %(message)s',
#     datefmt='%Y-%m-%d %H:%M:%S',
#     stream=sys.stderr)

# # Log the effective level being used
# logging.info(f"Logging level set to: {logging.getLevelName(logging.getLogger().getEffectiveLevel())}")


# # --- Define Core Paths ---
# INPUTDIR = Path(config["inputdir"])
# OUTDIR = Path(config["outdir"])
# LOGDIR = Path(config["logdir"])
# # Place indexdir within OUTDIR for better organization
# INDEXDIR = OUTDIR / config.get("indexdir", "bbmap_indices") # Use config value as subdir name
# GENOMESDIR = Path(config["genomesdir"])

# logging.info(f"Input directory: {INPUTDIR}")
# logging.info(f"Output directory: {OUTDIR}")
# logging.info(f"Log directory: {LOGDIR}")
# logging.info(f"Index directory: {INDEXDIR}")
# logging.info(f"Genomes directory: {GENOMESDIR}")

# # --- Include Rule Modules ---
# include: "rules/map.smk"

# # --- Rule All: Defines the final targets of the workflow ---
# localrules: all

# rule all:
#     input:
#         get_all_target_bams()
#         # get_all_target_stats()

# # from pathlib import Path
# import logging

# # --- Configuration ---
# configfile: "config.yaml"

# # --- Setup Logging ---
# logging.basicConfig(
#     level=config.get("log_level", "INFO"), # Control level via config?
#     format='%(asctime)s %(levelname)-8s %(name)s: %(message)s',
#     datefmt='%Y-%m-%d %H:%M:%S',
#     stream=sys.stderr # Or configure file handler etc.
# )
# logging.getLogger("snakemake").setLevel(logging.WARNING) # Example: Quiet snakemake internals

# # --- Define Core Paths ---
# # These are globally accessible in included files
# INPUTDIR = Path(config["inputdir"])
# OUTDIR = Path(config["outdir"])
# LOGDIR = Path(config["logdir"])
# # Use get for indexdir in case it's optional or named differently
# # Ensure indexdir is relative to OUTDIR if desired by config structure
# INDEXDIR = OUTDIR / config.get("indexdir", "bbmap_indices") # Example: Place under OUTDIR, default name
# GENOMESDIR = Path(config["genomesdir"])

# logging.info(f"Input directory: {INPUTDIR}")
# logging.info(f"Output directory: {OUTDIR}")
# logging.info(f"Log directory: {LOGDIR}")
# logging.info(f"Index directory: {INDEXDIR}")
# logging.info(f"Genomes directory: {GENOMESDIR}")

# # --- Include Rule Modules ---
# # Functions and rules within map.smk can now use INPUTDIR, OUTDIR, etc.
# include: "rules/map.smk"

# # --- Rule All: Defines the final targets of the workflow ---
# # Calls the function from map.smk to get the list of target files
# localrules: all

# rule all:
#     input:
#         # Call the function(s) defined in map.smk
#         get_all_target_bams()
#         # get_all_target_stats() # Uncomment if you have stats targets defined

# # from pathlib import Path


# # configfile: "config.yaml"

# # # The paths to the input and output directories
# # INPUTDIR = Path(config["inputdir"])
# # OUTDIR = Path(config["outdir"])
# # LOGDIR = Path(config["logdir"])
# # INDEXDIR= Path(config["indexdir"])
# # GENOMESDIR = Path(config["genomesdir"])

# # # Make sure the output directories exist
# # LOGDIR.mkdir(parents=True, exist_ok=True)
# # OUTDIR.mkdir(parents=True, exist_ok=True)
# # INDEXDIR.mkdir(parents=True, exist_ok=True)

# # all_outputs = []

# # # SAMPLES = set(glob_wildcards(INPUTDIR/config["sample_pattern"]).sample)
# # # if len(SAMPLES) < 1:
# # #     raise WorkflowError("Found no samples! Check input file pattern and path in config.yaml")


# # ###################################
# # # Extract reads from a clade
# # ###################################

# # include: "rules/map.smk"


# # localrules: all

# # rule all:
# #     input:
# #         all_outputs
