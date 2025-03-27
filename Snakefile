from pathlib import Path


configfile: "config.yaml"

# The paths to the input and output directories
INPUTDIR = Path(config["inputdir"])
OUTDIR = Path(config["outdir"])
LOGDIR = Path(config["logdir"])
INDEXDIR= Path(config["indexdir"])
GENOMESDIR = Path(config["genomesdir"])

# Make sure the output directories exist
LOGDIR.mkdir(parents=True, exist_ok=True)
OUTDIR.mkdir(parents=True, exist_ok=True)
INDEXDIR.mkdir(parents=True, exist_ok=True)

all_outputs = []

# SAMPLES = set(glob_wildcards(INPUTDIR/config["sample_pattern"]).sample)
# if len(SAMPLES) < 1:
#     raise WorkflowError("Found no samples! Check input file pattern and path in config.yaml")


###################################
# Extract reads from a clade
###################################

include: "rules/map.smk"


localrules: all

rule all:
    input:
        all_outputs
