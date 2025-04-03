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

# --- Define Output and Log Patterns ---
# DEDUP_BAM_OUT_PATTERN is already defined in map.smk and is available to us here since map.smk is included in the Snakefile before this file
SAMTOOLS_STATS_OUT_PATTERN = "{outdir}/validation/stats/{tax_id}/{genome_basename}/samtools_stats.txt"
SIM_PARAMS_OUT_PATTERN = "{outdir}/validation/params/{tax_id}/{genome_basename}/sim_params.json"
LOG_ESTIMATE_PARAMS_PATTERN = "{logdir}/estimate_sim_params/{tax_id}/{genome_basename}.log"
SIM_READS_R1_PATTERN = "{outdir}/validation/simulated_reads/{tax_id}/{genome_basename}/sim_r1.fastq.gz"
SIM_READS_R2_PATTERN = "{outdir}/validation/simulated_reads/{tax_id}/{genome_basename}/sim_r2.fastq.gz"
NEAT_CONFIG_YAML_PATTERN = "{outdir}/validation/simulated_reads/{tax_id}/{genome_basename}/neat_config.yml"
LOG_SIMULATE_READS_PATTERN = "{logdir}/simulate_reads/{tax_id}/{genome_basename}.log"

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
        runtime=config.get("SAMTOOLS_STATS_RUNTIME", "20m"),
        cpus_per_task=config.get("SAMTOOLS_STATS_THREADS", 1)
    script:
        "../scripts/run_samtools_stats_parse.py"

rule simulate_reads_neat:
    input:
        # The FASTA file to simulate reads from
        ref_fasta = get_actual_fasta_path,
        # Parameters from previous step
        params_json = SIM_PARAMS_OUT_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    output:
        # Output simulated reads (fastq.gz files)
        r1 = SIM_READS_R1_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        ),
        r2 = SIM_READS_R2_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        ),
        neat_config_yaml = NEAT_CONFIG_YAML_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    log:
        path = LOG_SIMULATE_READS_PATTERN.format(
            logdir=LOG_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    conda:
        ".." / ENVS_DIR_P / "simulate_reads.yaml"
    threads: config.get("NEAT_SIM_THREADS", 1)
    resources:
        runtime=config.get("NEAT_SIM_RUNTIME", "2h"),
        mem_mb=config.get("NEAT_SIM_MEM_MB", 4000),
        cpus_per_task=config.get("NEAT_SIM_THREADS", 1)
    script:
        "../scripts/run_neat_simulate_reads.py"