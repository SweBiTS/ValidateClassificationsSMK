import sys
import os
from pathlib import Path
import logging
import json
import re
import subprocess

# ==================================================
#      Housekeeping and setup
# ==================================================

# --- Logging ---
logger = logging.getLogger(__name__) # Use module-specific logger

# --- Define Output Patterns ---
# DEDUP_BAM_OUT_PATTERN is already defined in map.smk and is available to us here since map.smk is included in the Snakefile before this file
SAMTOOLS_STATS_OUT_PATTERN = "{outdir}/validation/stats/{tax_id}/{genome_basename}/samtools_stats.txt"
SIM_PARAMS_OUT_PATTERN = "{outdir}/validation/params/{tax_id}/{genome_basename}/sim_params.json"
LOG_ESTIMATE_PARAMS_PATTERN = "{logdir}/estimate_sim_params/{tax_id}/{genome_basename}.log"

# --- Rule: Estimate Simulation Parameters using samtools stats ---
rule estimate_simulation_params:
    input:
        # The deduplicated BAM file from the mapping step
        bam = DEDUP_BAM_OUT_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    output:
        # Raw output from samtools stats
        stats_raw = SAMTOOLS_STATS_OUT_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        ),
        # JSON file with parsed parameters for the read simulation
        params_json = SIM_PARAMS_OUT_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    log:
        path = LOG_ESTIMATE_PARAMS_PATTERN.format(
            logdir=LOG_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    conda:
        ".." / ENVS_DIR_P / "map.yaml"
    threads: config.get("SAMTOOLS_STATS_THREADS", 1)
    resources:
        mem_mb=config.get("SAMTOOLS_STATS_MEM_MB", 2000),
        runtime=config.get("SAMTOOLS_STATS_RUNTIME", "20m")
    script:
        "../scripts/run_samtools_stats_parse.py"