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
SIM_READS_R1_PATTERN = "{outdir}/validation/simulated_reads/{tax_id}/{genome_basename}/sim_r1.fastq.gz"
SIM_READS_R2_PATTERN = "{outdir}/validation/simulated_reads/{tax_id}/{genome_basename}/sim_r2.fastq.gz"
NEAT_CONFIG_YAML_PATTERN = "{outdir}/validation/simulated_reads/{tax_id}/{genome_basename}/neat_config.yml"
CLEANED_SIM_READS_R1_PATTERN = "{outdir}/validation/simulated_reads/{tax_id}/{genome_basename}/sim_clean_r1.fq.gz"
CLEANED_SIM_READS_R2_PATTERN = "{outdir}/validation/simulated_reads/{tax_id}/{genome_basename}/sim_clean_r2.fq.gz"
KRAKEN_OUT_SIM_PATTERN = "{outdir}/validation/kraken2_sim/{tax_id}/{genome_basename}/sim.kraken.out"
KRAKEN_REPORT_SIM_PATTERN = "{outdir}/validation/kraken2_sim/{tax_id}/{genome_basename}/sim.kraken.report"
CORRECT_READ_IDS_PATTERN = "{outdir}/validation/filtered_reads/{tax_id}/{genome_basename}/sim_correct_read_ids.txt"
FILTERED_SIM_READS_R1_PATTERN = "{outdir}/validation/filtered_reads/{tax_id}/{genome_basename}/sim_correct_R1.fq.gz"
FILTERED_SIM_READS_R2_PATTERN = "{outdir}/validation/filtered_reads/{tax_id}/{genome_basename}/sim_correct_R2.fq.gz"
MAPPED_SIM_BAM_PATTERN = "{outdir}/validation/mapped_sim/{tax_id}/{genome_basename}/mapped_correct_sim.bam"
MAPPED_SIM_STATS_PATTERN = "{outdir}/validation/mapped_sim/{tax_id}/{genome_basename}/mapped_correct_sim_stats.txt"
SORTED_MAPPED_SIM_BAM_PATTERN = "{outdir}/validation/mapped_sim/{tax_id}/{genome_basename}/mapped_correct_sim.sorted.bam"
MASKED_REGIONS_BED_PATTERN = "{outdir}/validation/masked_regions/{genome_basename}.masked.bed"
COVERED_REGIONS_BED_PATTERN = "{outdir}/validation/coverage_sim/{tax_id}/{genome_basename}/covered_regions.bed"
MASKED_FASTA_PATTERN = "{outdir}/validation/masked_regions/{genome_basename}.masked_lowercase.fasta"
CLASSIFIABLE_REGIONS_BED_PATTERN = "{outdir}/validation/classifiable_regions/{tax_id}/{genome_basename}/classifiable_regions.bed"
OVERLAPPING_READS_BAM_PATTERN = "{outdir}/validation/final_analysis/{tax_id}/{genome_basename}/overlapping_reads.bam"
PER_SAMPLE_SUMMARY_PATTERN = "{outdir}/validation/summary_reports/{tax_id}/{genome_basename}/summary.tsv"
PLOT_SCATTER_HTML_PATTERN  = "{outdir}/validation/plots/{tax_id}/{genome_basename}/scatter_reads_vs_length.html"

# --- Define Log Patterns ---
LOG_ESTIMATE_PARAMS_PATTERN = "{logdir}/estimate_sim_params/{tax_id}/{genome_basename}.log"
LOG_SIMULATE_READS_PATTERN = "{logdir}/simulate_reads/{tax_id}/{genome_basename}.log"
LOG_CLEAN_FASTQ_PATTERN = "{logdir}/clean_fastq/{tax_id}/{genome_basename}.log"
LOG_KRAKEN2_SIM_PATTERN = "{logdir}/kraken2_sim/{tax_id}/{genome_basename}.log"
LOG_EXTRACT_IDS_PATTERN = "{logdir}/extract_ids/{tax_id}/{genome_basename}.log"
LOG_FILTER_FASTQ_PATTERN = "{logdir}/filter_fastq/{tax_id}/{genome_basename}.log"
LOG_MAP_SIM_PATTERN = "{logdir}/map_simulated/{tax_id}/{genome_basename}.log"
LOG_SORT_SIM_PATTERN = "{logdir}/sort_simulated/{tax_id}/{genome_basename}.log"
LOG_CALC_SIM_COV_PATTERN = "{logdir}/calculate_sim_coverage/{tax_id}/{genome_basename}.log"
LOG_MASKED_REGIONS_PATTERN = "{logdir}/find_masked_regions/{genome_basename}.log"
LOG_DEFINE_CLASSIFIABLE_PATTERN = "{logdir}/define_classifiable_regions/{tax_id}/{genome_basename}.log"
LOG_INTERSECT_READS_PATTERN = "{logdir}/intersect_reads/{tax_id}/{genome_basename}.log"
LOG_SUMMARIZE_VALIDATION_PATTERN = "{logdir}/summarize_validation/{tax_id}/{genome_basename}.log"
LOG_PLOT_SCATTER_PATTERN = "{logdir}/plot_scatter/{tax_id}/{genome_basename}.log"

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

# --- Rule: Simulate error-free reads from original (unmasked) genome using NEAT ---
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

# --- Rule: Clean FASTQ Headers (Remove /1, /2) using awk  ---
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
        # Output are cleaned FASTQ files
        r1 = temp(CLEANED_SIM_READS_R1_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )),
        r2 = temp(CLEANED_SIM_READS_R2_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        ))
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
        r""" # <-- Add 'r' here to make it a raw string
        mkdir -p $(dirname {output.r1}) &&
        mkdir -p $(dirname {log.path}) &&
        {{
            (pigz -dc {input.r1} | awk 'NR % 4 == 1 {{sub(/\/[12]$/, "", $1)}} {{print}}' | pigz -c > {output.r1}) &&
            (pigz -dc {input.r2} | awk 'NR % 4 == 1 {{sub(/\/[12]$/, "", $1)}} {{print}}' | pigz -c > {output.r2}) ;
        }} 2> {log.path}
        """

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
        kraken_out = temp(KRAKEN_OUT_SIM_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )),
        report = KRAKEN_REPORT_SIM_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    params:
        db_path = config["KRAKEN2_DB_PATH"],
        # Determine executable path: use specified path or default 'kraken2' command
        # Path already validated in main Snakefile, if no path is given we assume conda execution
        executable = str(Path(config["KRAKEN2_BIN_DIR"]) / "kraken2") if config.get("KRAKEN2_BIN_DIR") else "kraken2"
    log:
        path = LOG_KRAKEN2_SIM_PATTERN.format(
            logdir=LOG_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    conda:
        ".." / ENVS_DIR_P / "classify.yaml" if not config.get("KRAKEN2_BIN_DIR") else None
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
        read_ids = temp(CORRECT_READ_IDS_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        ))
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

# --- Rule: Filter correctly classified simulated reads from FASTQ files using seqtk ---
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
        r1 = temp(FILTERED_SIM_READS_R1_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )),
        r2 = temp(FILTERED_SIM_READS_R2_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        ))
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

# --- Rule: Map Correctly Classified Simulated Reads using BBMap ---
rule map_correct_simulated:
    input:
        # Use the filtered FASTQ files from seqtk step
        r1 = FILTERED_SIM_READS_R1_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        ),
        r2 = FILTERED_SIM_READS_R2_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        ),
        # Reference the index created by bbmap_index rule in map.smk
        idx_dir = rules.bbmap_index.output.idx_dir # Uses {genome_basename} wildcard
    output:
        # Output BAM and mapping stats
        bam = temp(MAPPED_SIM_BAM_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )),
        stats = MAPPED_SIM_STATS_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    params:
        build_id = 1, # Assuming build_id is always 1 for this workflow
        mem = config.get("BBMAP_MAP_MEM", "20g"),
        minid = config["BBMAP_MINID"],
        pairedonly = config["BBMAP_PAIREDONLY"]
    log:
        path = LOG_MAP_SIM_PATTERN.format(
            logdir=LOG_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    conda:
        ".." / ENVS_DIR_P / "map.yaml"
    threads: config.get("BBMAP_MAP_THREADS", 4)
    resources:
        runtime=config.get("BBMAP_MAP_RUNTIME", "1h"),
        mem_mb=config.get("BBMAP_MAP_MEM_MB", 8000),
        cpus_per_task=config.get("BBMAP_MAP_THREADS", 4)
    shell:
        "mkdir -p $(dirname {output.bam}) && "
        "mkdir -p $(dirname {log.path}) && "
        "/usr/bin/time -v "
        "bbmap.sh "
            "in1={input.r1} "
            "in2={input.r2} "
            "path={input.idx_dir} "
            "build={params.build_id} "
            "out={output.bam} "
            "statsfile={output.stats} "
            "mappedonly=t "
            "minid={params.minid} "
            "ambiguous=random "
            "pairedonly={params.pairedonly} "
            "-Xmx{params.mem} "
            "threads={threads} "
            "pigz=t unpigz=t "
            "overwrite=t "
            f"{BBMAP_JNI_FLAG} "
            "> {log.path} 2>&1" # Redirect stderr to log

# --- Rule: Sort Mapped Simulated Reads BAM file by coordinate ---
rule sort_mapped_simulated:
    input:
        # Input is the unsorted BAM from map_correct_simulated
        bam = MAPPED_SIM_BAM_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    output:
        bam_sorted = SORTED_MAPPED_SIM_BAM_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    params:
        mem_per_thread = config.get("SAMTOOLS_SORT_MEM_PER_THREAD", "2G"),
    log:
        path = LOG_SORT_SIM_PATTERN.format(
            logdir=LOG_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    conda:
        ".." / ENVS_DIR_P / "map.yaml"
    threads: config.get("SAMTOOLS_SORT_THREADS", 2)
    resources:
        mem_mb=lambda wildcards, threads: int(threads) * int(config.get("SAMTOOLS_SORT_MEM_PER_THREAD_MB", 2048)),
        runtime=config.get("SAMTOOLS_SORT_RUNTIME", "2h"),
        cpus_per_task=config.get("SAMTOOLS_SORT_THREADS", 2)
    shell:
        "mkdir -p $(dirname {output.bam_sorted}) && "
        "mkdir -p $(dirname {log.path}) && "
        "/usr/bin/time -v samtools sort "
            "-@ {threads} "
            "-m {params.mem_per_thread} "
            "-o {output.bam_sorted} "
            "{input.bam} "
            "2> {log.path}"

# --- Rule: Calculate Coverage from Mapped Simulated Reads ---
rule calculate_simulated_coverage:
    input:
        # The sorted BAM file containing mapped simulated reads
        bam = SORTED_MAPPED_SIM_BAM_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    output:
        # BED file of regions covered >= threshold (marked temp)
        covered_bed = COVERED_REGIONS_BED_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    params:
        threshold = 1  # Coverage threshold (site is kept if coverage >= threshold)
    log:
        path = LOG_CALC_SIM_COV_PATTERN.format(
            logdir=LOG_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    conda:
        ".." / ENVS_DIR_P / "analysis.yaml"
    threads: config.get("CALC_SIM_COV_THREADS", 1)
    resources:
        runtime=config.get("CALC_SIM_COV_RUNTIME", "30m"),
        mem_mb=config.get("CALC_SIM_COV_MEM_MB", 2000),
        cpus_per_task=config.get("CALC_SIM_COV_THREADS", 1)
    shell:
        # Use bedtools genomecov -bga, coverage filter with awk, merge resulting bed-regions
        "mkdir -p $(dirname {output.covered_bed}) && "
        "mkdir -p $(dirname {log.path}) && "
        "bedtools genomecov -bga -ibam {input.bam} "
            "| awk -v threshold={params.threshold} '$4 >= threshold' "
            "| bedtools merge -i stdin "
            "> {output.covered_bed} "
            "2> {log.path}"

# --- Rule: Identify Kraken 2 Masked Regions using k2mask ---
rule find_masked_regions:
    input:
        # Original reference genome
        ref_fasta = get_actual_fasta_path
    output:
        # Final BED file of masked regions
        masked_bed = MASKED_REGIONS_BED_PATTERN.format(
            outdir=OUTPUT_DIR_P, genome_basename="{genome_basename}"
        ),
        # Intermediate lowercase masked fasta
        masked_fasta = MASKED_FASTA_PATTERN.format(
            outdir=OUTPUT_DIR_P, genome_basename="{genome_basename}"
        )
    params:
        # Determine k2mask executable path (uses KRAKEN2_BIN_DIR if set, else 'k2mask')
        k2mask_exec = str(Path(config["KRAKEN2_BIN_DIR"]) / "k2mask") if config.get("KRAKEN2_BIN_DIR") else "k2mask"
    log:
        path = LOG_MASKED_REGIONS_PATTERN.format(
            logdir=LOG_DIR_P, genome_basename="{genome_basename}"
        )
    conda:
        # fasta_to_masked_bed.py needs k2mask (Kraken 2)
        ".." / ENVS_DIR_P / "classify.yaml" if not config.get("KRAKEN2_BIN_DIR") else None
    threads: config.get("K2MASK_THREADS", 1)
    resources:
        runtime=config.get("K2MASK_RUNTIME", "60m"),
        mem_mb=config.get("K2MASK_MEM_MB", 2000),
        cpus_per_task=config.get("K2MASK_THREADS", 1)
    script:
        "../scripts/fasta_to_masked_bed.py"

# --- Rule: Define Classifiable Regions (subtract the masked regions (k2mask)) ---
rule define_classifiable_regions:
    input:
        # Regions covered by correctly classified simulated reads
        coverage_bed = COVERED_REGIONS_BED_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        ),
        # BED file of low-complexity regions (output of find_masked_regions)
        masked_bed = MASKED_REGIONS_BED_PATTERN.format(
            outdir=OUTPUT_DIR_P, genome_basename="{genome_basename}"
        )
    output:
        # Final BED file of classifiable regions, i.e., regions covered by 
        # simulated and correctly kraken classified reads AND that are not masked
        classifiable_bed = CLASSIFIABLE_REGIONS_BED_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    log:
        path = LOG_DEFINE_CLASSIFIABLE_PATTERN.format(
            logdir=LOG_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    conda:
        ".." / ENVS_DIR_P / "analysis.yaml"
    threads: config.get("BEDTOOLS_SUBTRACT_THREADS", 1)
    resources:
        runtime=config.get("BEDTOOLS_SUBTRACT_RUNTIME", "15m"),
        mem_mb=config.get("BEDTOOLS_SUBTRACT_MEM_MB", 2000),
        cpus_per_task=config.get("BEDTOOLS_SUBTRACT_THREADS", 1)
    shell:
        # Subtract masked regions (-b) from covered regions (-a)
        "mkdir -p $(dirname {output.classifiable_bed}) && "
        "mkdir -p $(dirname {log.path}) && "
        "bedtools subtract -a {input.coverage_bed} -b {input.masked_bed} "
        "> {output.classifiable_bed} "
        "2> {log.path}"

# --- Rule: Intersect Real Reads with Classifiable Regions ---
rule intersect_real_reads_with_classifiable_regions:
    input:
        # Use the deduplicated BAM of REAL reads from map.smk
        bam = DEDUP_BAM_OUT_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        ),
        # Use the final classifiable regions BED
        regions = CLASSIFIABLE_REGIONS_BED_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    output:
        # Output BAM containing unique reads overlapping classifiable regions
        overlapping_bam = OVERLAPPING_READS_BAM_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    log:
        path = LOG_INTERSECT_READS_PATTERN.format(
            logdir=LOG_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    conda:
        ".." / ENVS_DIR_P / "analysis.yaml"
    threads: config.get("INTERSECT_READS_THREADS", 1)
    resources:
        runtime=config.get("INTERSECT_READS_RUNTIME", "30m"),
        mem_mb=config.get("INTERSECT_READS_MEM_MB", 2000),
        cpus_per_task=config.get("INTERSECT_READS_THREADS", 1)
    shell:
        # Use bedtools intersect -u; output is BAM. Use -sorted for speed.
        "mkdir -p $(dirname {output.overlapping_bam}) && "
        "mkdir -p $(dirname {log.path}) && "
        "bedtools intersect "
            "-a {input.bam} "               # Input BAM (real reads)
            "-b {input.regions} "           # Classifiable regions BED
            "-u "                           # Report each alignment record in -a once if *any* overlap exists
            "-sorted "                      # Assume inputs are sorted for speed
            "-wa "                          # Write the original BAM records
            "> {output.overlapping_bam} "
            "2> {log.path}"

# --- Rule: Summarize Validation Coverage Per Contig ---
rule summarize_validation_coverage:
    input:
        # BAM containing only alignment records overlapping classifiable regions
        overlapping_bam = OVERLAPPING_READS_BAM_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        ),
        # BED file defining the classifiable regions (to calculate length)
        classifiable_bed = CLASSIFIABLE_REGIONS_BED_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    output:
        summary_tsv = PER_SAMPLE_SUMMARY_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    log:
        path = LOG_SUMMARIZE_VALIDATION_PATTERN.format(
            logdir=LOG_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    conda:
        ".." / ENVS_DIR_P / "analysis.yaml"
    threads: config.get("SUMMARIZE_VALIDATION_THREADS", 1)
    resources:
        runtime=config.get("SUMMARIZE_VALIDATION_RUNTIME", "20m"), # Integer minutes
        mem_mb=config.get("SUMMARIZE_VALIDATION_MEM_MB", 2000),
        cpus_per_task=config.get("SUMMARIZE_VALIDATION_THREADS", 1)
    script:
        "../scripts/summarize_validation.py"

# --- Rule: Create Validation Scatter Plot ---
rule plot_validation_scatter:
    input:
        # Per-sample summary TSV from summarize_validation_coverage
        summary_tsv = PER_SAMPLE_SUMMARY_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    output:
        plot_html = PLOT_SCATTER_HTML_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    log:
        path = LOG_PLOT_SCATTER_PATTERN.format(
            logdir=LOG_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    conda:
        ".." / ENVS_DIR_P / "analysis.yaml"
    threads: config.get("PLOT_THREADS", 1)
    resources:
        runtime=config.get("PLOT_RUNTIME", "15m"),
        mem_mb=config.get("PLOT_MEM_MB", 1000),
        cpus_per_task=config.get("PLOT_THREADS", 1)
    script:
        "../scripts/plot_scatter_plotly.py"