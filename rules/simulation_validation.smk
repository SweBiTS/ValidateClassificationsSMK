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
SIM_READS_R1_PATTERN = "{outdir}/validation/simulated_reads/{tax_id}/{genome_basename}/sim_raw_r1.fastq.gz"
SIM_READS_R2_PATTERN = "{outdir}/validation/simulated_reads/{tax_id}/{genome_basename}/sim_raw_r2.fastq.gz"
NEAT_CONFIG_YAML_PATTERN = "{outdir}/validation/simulated_reads/{tax_id}/{genome_basename}/neat_config.yml"
LOG_SIMULATE_READS_PATTERN = "{logdir}/simulate_reads/{tax_id}/{genome_basename}.log"
CLEANED_SIM_READS_R1_PATTERN = "{outdir}/validation/filtered_reads/{tax_id}/{genome_basename}/sim_R1.fq.gz"
CLEANED_SIM_READS_R2_PATTERN = "{outdir}/validation/filtered_reads/{tax_id}/{genome_basename}/sim_R2.fq.gz"
LOG_CLEAN_FASTQ_PATTERN = "{logdir}/clean_fastq/{tax_id}/{genome_basename}.log"
KRAKEN_OUT_SIM_PATTERN = "{outdir}/validation/kraken2_sim/{tax_id}/{genome_basename}/sim.kraken.out"
KRAKEN_REPORT_SIM_PATTERN = "{outdir}/validation/kraken2_sim/{tax_id}/{genome_basename}/sim.kraken.report"
LOG_KRAKEN2_SIM_PATTERN = "{logdir}/kraken2_sim/{tax_id}/{genome_basename}.log"
CORRECT_READ_IDS_PATTERN = "{outdir}/validation/filtered_reads/{tax_id}/{genome_basename}/sim_correct_read_ids.txt"
LOG_EXTRACT_IDS_PATTERN = "{logdir}/extract_ids/{tax_id}/{genome_basename}.log"
FILTERED_SIM_READS_R1_PATTERN = "{outdir}/validation/filtered_reads/{tax_id}/{genome_basename}/sim_correct_R1.fq.gz"
FILTERED_SIM_READS_R2_PATTERN = "{outdir}/validation/filtered_reads/{tax_id}/{genome_basename}/sim_correct_R2.fq.gz"
LOG_FILTER_FASTQ_PATTERN = "{logdir}/filter_fastq/{tax_id}/{genome_basename}.log"


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
        r1 = temp(SIM_READS_R1_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )),
        r2 = temp(SIM_READS_R2_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )),
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

# --- Rule: Clean FASTQ Headers (Remove /1, /2) using sed  ---
rule clean_fastq_headers:
    input:
        # Input are the raw simulated reads from NEAT (that have appended "/1" or "/2" to the headers)
        r1 = SIM_READS_R1_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        ),
        r2 = SIM_READS_R2_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    output:
        # Output are cleaned FASTQ files, marked as temp
        r1 = CLEANED_SIM_READS_R1_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        ),
        r2 = CLEANED_SIM_READS_R2_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    log:
        path = LOG_CLEAN_FASTQ_PATTERN.format(
            logdir=LOG_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    conda:
        ".." / ENVS_DIR_P / "extract_reads.yaml"
    threads: config.get("CLEAN_FASTQ_THREADS", 1)
    resources:
        runtime=config.get("CLEAN_FASTQ_RUNTIME", "60m"),
        mem_mb=config.get("CLEAN_FASTQ_MEM_MB", 1000),
        cpus_per_task=config.get("CLEAN_FASTQ_THREADS", 1)
    shell:
        # Use pigz for decompression/compression, sed for substitution
        # Redirect stderr for the whole block to the log file
        "mkdir -p $(dirname {output.r1}) && "
        "mkdir -p $(dirname {log.path}) && "
        "{{ "
            "(pigz -dc {input.r1} | sed -E 's/^(@[^[:space:]]+)\/[12]([[:space:]].*|$)/\1\2/' | pigz -c > {output.r1}) && "
            "(pigz -dc {input.r2} | sed -E 's/^(@[^[:space:]]+)\/[12]([[:space:]].*|$)/\1\2/' | pigz -c > {output.r2}) ; "
        "}} 2> {log.path}"

# --- Rule: Classify Simulated Reads using Kraken2 ---
rule classify_simulated:
    input:
        r1 = CLEANED_SIM_READS_R1_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        ),
        r2 = CLEANED_SIM_READS_R2_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    output:
        kraken_out = KRAKEN_OUT_SIM_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        ),
        report = KRAKEN_REPORT_SIM_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    params:
        db_path = config["KRAKEN2_DB_PATH"],
        # Determine executable path: use specified path or default 'kraken2' command
        executable = config.get("KRAKEN2_EXEC", "kraken2")
    log:
        path = LOG_KRAKEN2_SIM_PATTERN.format(
            logdir=LOG_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    conda:
        ".." / ENVS_DIR_P / "classify.yaml" if not config.get("KRAKEN2_EXEC") else None
    threads: config.get("KRAKEN2_THREADS", 8)
    resources:
        runtime=config.get("KRAKEN2_RUNTIME", "2h"),
        mem_mb=config.get("KRAKEN2_MEM_MB", 32000),
        cpus_per_task=config.get("KRAKEN2_THREADS", 8)
    shell:
        "mkdir -p $(dirname {output.kraken_out}) && "
        "mkdir -p $(dirname {log.path}) && "
        "/usr/bin/time -v "
        "{params.executable} "
            "--db {params.db_path} "
            "--threads {threads} "
            "--paired "
            "--output {output.kraken_out} "
            "--report {output.report} "
            "--gzip-compressed "
            "{input.r1} {input.r2} "
            "2> {log.path}" # Redirect stderr to log

# --- Rule: Extract Correctly Classified Read IDs using selmeout.py ---
rule extract_correct_read_ids:
    input:
        # Kraken output from simulated reads classification
        kraken_out = KRAKEN_OUT_SIM_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        ),
        # Taxonomy files (paths from config)
        nodes = config["TAXONOMY_NODES_PATH"],
        names = config["TAXONOMY_NAMES_PATH"]
    output:
        # List of read IDs (classified_taxid, root_taxid, read_id)
        read_ids = CORRECT_READ_IDS_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    params:
        selmeout_path = "scripts/selmeout.py",
        tax_id = "{tax_id}",
        mode = "clade"
    log:
        path = LOG_EXTRACT_IDS_PATTERN.format(
            logdir=LOG_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    conda:
        ".." / ENVS_DIR_P / "extract_reads.yaml"
    threads: config.get("EXTRACT_IDS_THREADS", 1)
    resources:
        runtime=config.get("EXTRACT_IDS_RUNTIME", "60m"),
        mem_mb=config.get("EXTRACT_IDS_MEM_MB", 2000),
        cpus_per_task=config.get("EXTRACT_IDS_THREADS", 1)
    shell:
        # Run the selmeout.py script
        "mkdir -p $(dirname {output.read_ids}) && "
        "mkdir -p $(dirname {log.path}) && "
        "/usr/bin/time -v "
        "python {params.selmeout_path} "
            "--input {input.kraken_out} "
            "--output {output.read_ids} "
            "--tax_id {params.tax_id} "
            "--mode {params.mode} "
            "--nodes {input.nodes} "
            "--names {input.names} "
            "> {log.path} 2>&1" # Capture script stdout/stderr

# --- Rule: Filter FASTQ files using seqtk ---
rule filter_fastq_with_seqtk:
    input:
        # Original simulated FASTQ files
        r1 = CLEANED_SIM_READS_R1_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        ),
        r2 = CLEANED_SIM_READS_R2_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        ),
        # List of read IDs (classified_taxid, root_taxid, read_id) from previous rule
        read_ids = CORRECT_READ_IDS_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    output:
        # Final filtered FASTQ files (gzipped, marked temp)
        r1 = FILTERED_SIM_READS_R1_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        ),
        r2 = FILTERED_SIM_READS_R2_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    log:
        path = LOG_FILTER_FASTQ_PATTERN.format(
            logdir=LOG_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    conda:
        ".." / ENVS_DIR_P / "extract_reads.yaml"
    threads: config.get("SEQTK_THREADS", 1)
    resources:
        runtime=config.get("SEQTK_RUNTIME", "60m"),
        mem_mb=config.get("SEQTK_MEM_MB", 2000),
        cpus_per_task=config.get("SEQTK_THREADS", 1)
    shell:
        "mkdir -p $(dirname {output.r1}) && "
        "mkdir -p $(dirname {log.path}) && "
        # Start brace group to capture stderr from all commands within
        "{{ "
            "( /usr/bin/time -v seqtk subseq {input.r1} <(cut -f 3 {input.read_ids}) | pigz -c > {output.r1} ) && "
            "( /usr/bin/time -v seqtk subseq {input.r2} <(cut -f 3 {input.read_ids}) | pigz -c > {output.r2} ) ; "
        # Redirect stderr for the entire brace group to the log file
        "}} 2> {log.path}"