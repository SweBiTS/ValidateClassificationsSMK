import yaml
import sys
import os
from pathlib import Path
import logging

# TODO: Implement dynamic resource allocation that depend on the size of the input files.
# TODO: Implement resubmission of jobs that fail, that increases the memory and/or time resources based on the "attempts" variable

# ==================================================
#      Housekeeping
# ==================================================

# --- Logging ---
# Logger specific to this module - inherits config from root logger
logger = logging.getLogger(__name__)


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


# ==================================================
#      Mapping setup
# ==================================================

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


# --- Helper Function: Get Genome Basename (Stem) ---
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
        f"{extensions_to_check} in directory '{genomes_dir}'"
    )

# --- Generate Target Output Files Function (for rule all) ---
def get_all_target_outputs(mapping_spec, outdir_pattern_template):
    """Generates list of expected final output files based on the mapping spec."""
    target_files = []
    logger.info(f"Generating target output file list using pattern: {outdir_pattern_template}")
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
                target_path = Path(outdir_pattern_template.format(
                    outdir=OUTPUT_DIR_P,
                    tax_id=tax_id,
                    genome_basename=genome_basename
                ))
                target_files.append(str(target_path))
                count += 1
            else:
                warning_count += 1
                pass
    logger.info(f"Generated {count} target output files.")
    if warning_count > 0:
        logger.warning(f"Encountered {warning_count} issues while processing mapping specification.")
    return target_files

# --- Define Output Patterns ---
RAW_BAM_OUT_PATTERN = "{outdir}/mapping/{tax_id}/{genome_basename}/mapping.bam"
SORTED_BAM_OUT_PATTERN = "{outdir}/mapping/{tax_id}/{genome_basename}/mapping.sorted.bam"
DEDUP_BAM_OUT_PATTERN = "{outdir}/mapping/{tax_id}/{genome_basename}/mapping_dedup.bam"
DEDUP_BAM_INDEX_OUT_PATTERN = "{outdir}/mapping/{tax_id}/{genome_basename}/mapping_dedup.bam.bai"
STATS_OUT_PATTERN = "{outdir}/mapping/{tax_id}/{genome_basename}/mapping_stats.txt"

# --- Define Log Patterns ---
LOG_BBMAP_INDEX_PATTERN = "{logdir}/bbmap_index/{genome_basename}.log"
LOG_MAP_PATTERN = "{logdir}/bbmap_map/{tax_id}_{genome_basename}.log"
LOG_SORT_PATTERN = "{logdir}/samtools_sort/{tax_id}_{genome_basename}.log"
LOG_MARKDUP_PATTERN = "{logdir}/sambamba_markdup/{tax_id}_{genome_basename}.log"

# --- Define Functions for Rule All ---
# Targets are indices of sorted bam files (last step)
def get_all_target_indices():
     return get_all_target_outputs(mapping_spec_data, DEDUP_BAM_INDEX_OUT_PATTERN)


# ==================================================
#      SNAKEMAKE RULES
# ==================================================

# --- Rule: Index Reference Genome using BBMap ---
rule bbmap_index:
    input:
        fasta = get_actual_fasta_path
    output:
        idx_dir = directory(str(INDEX_DIR_BBMAP_P / "{genome_basename}"))
    params:
        build_id = 1,
        mem = config.get("BBMAP_INDEX_MEM", "20g")
    log:
        path = LOG_BBMAP_INDEX_PATTERN.format(
            logdir=LOG_DIR_P, genome_basename="{genome_basename}"
        )
    conda:
        ".." / ENVS_DIR_P / "map.yaml"
    threads: config.get("BBMAP_INDEX_THREADS", 4)
    resources:
        mem_mb=config.get("BBMAP_INDEX_MEM_MB", 22000),
        runtime=config.get("BBMAP_INDEX_RUNTIME", "1h"),
        cpus_per_task=config.get("BBMAP_INDEX_THREADS", 4)
    shell:
        "mkdir -p $(dirname {log.path}) && "
        "/usr/bin/time -v "
        "bbmap.sh "
            "ref={input.fasta} "
            "path={output.idx_dir} "
            "build={params.build_id} "
            "-Xmx{params.mem} "
            "threads={threads} "
            "pigz=t unpigz=t "
            "overwrite=t "
            f"{BBMAP_JNI_FLAG} "  # Adds 'usejni=t' or 'usejni=f'
            "> {log.path} 2>&1"


# --- Rule: Map Reads to Reference using BBMap ---
rule bbmap_map_reads:
    input:
        idx_dir = rules.bbmap_index.output.idx_dir,
        r1 = lambda wildcards: str(INPUT_DIR_P / config["FASTQ_FILE_PATTERN"].format(tax_id=wildcards.tax_id, direction="1")),
        r2 = lambda wildcards: str(INPUT_DIR_P / config["FASTQ_FILE_PATTERN"].format(tax_id=wildcards.tax_id, direction="2"))
    output:
        bam = temp(RAW_BAM_OUT_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )),
        stats = STATS_OUT_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    params:
        build_id = 1,
        mem = config.get("BBMAP_MAP_MEM", "6g"),
        minid = config["BBMAP_MINID"],
        ambig = config["BBMAP_AMBIGUOUS"],
        pairedonly = config["BBMAP_PAIREDONLY"]
    log:
        path = LOG_MAP_PATTERN.format(
            logdir=LOG_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    conda:
        ".." / ENVS_DIR_P / "map.yaml"
    threads: config.get("BBMAP_MAP_THREADS", 4)
    resources:
        mem_mb=config.get("BBMAP_MAP_MEM_MB", 8000),
        runtime=config.get("BBMAP_MAP_RUNTIME", "1h"),
        cpus_per_task=config.get("BBMAP_MAP_THREADS", 4)
    shell:
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
            "ambiguous={params.ambig} "
            "pairedonly={params.pairedonly} "
            "-Xmx{params.mem} "
            "threads={threads} "
            "pigz=t unpigz=t "
            "overwrite=t "
            f"{BBMAP_JNI_FLAG} "  # Adds 'usejni=t' or 'usejni=f'
        "> {log.path} 2>&1"


# --- Rule: Sort BAM file using Samtools ---
rule samtools_sort:
    input:
        bam = RAW_BAM_OUT_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    output:
        bam = SORTED_BAM_OUT_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    params:
        mem_per_thread = config.get("SAMTOOLS_SORT_MEM_PER_THREAD", "2G"), # Memory per thread for sorting
    log:
        path = LOG_SORT_PATTERN.format(
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
        "mkdir -p $(dirname {log.path}) && "
        "/usr/bin/time -v "
        "samtools sort "
        "-@ {threads} "                 # Number of sorting threads
        "-m {params.mem_per_thread} "   # Memory per thread
        "-o {output.bam} "              # Output file
        "{input.bam} "                  # Input file
        "2> {log.path}"                 # Log stderr


# --- Rule: Mark/Remove Duplicates using Sambamba ---
rule mark_duplicates_sambamba:
    input:
        # Input is the coordinate-sorted BAM from samtools_sort
        sorted_bam = SORTED_BAM_OUT_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    output:
        # Output is the BAM file with duplicates removed
        dedup_bam = DEDUP_BAM_OUT_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        ),
        index = DEDUP_BAM_INDEX_OUT_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    params:
        # Define a temporary directory relative to the output
        tmpdir = lambda wildcards, output: Path(output.dedup_bam).parent / "tmp_sambamba_markdup"
    log:
        path = LOG_MARKDUP_PATTERN.format(
            logdir=LOG_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    conda:
         ".." / ENVS_DIR_P / "map.yaml"
    threads: config.get("SAMBAMBA_MARKDUP_THREADS", 4)
    resources:
        mem_mb=config.get("SAMBAMBA_MARKDUP_MEM_MB", 8000),
        runtime=config.get("SAMBAMBA_MARKDUP_RUNTIME", "2h"),
        cpus_per_task=config.get("SAMBAMBA_MARKDUP_THREADS", 4)
    shell:
        "mkdir -p $(dirname {log.path}) && "
        "mkdir -p {params.tmpdir} && "
        "/usr/bin/time -v "
        "sambamba markdup -t {threads} "
        "--tmpdir={params.tmpdir} "
        "--remove-duplicates "
        "{input.sorted_bam} "
        "{output.dedup_bam} "
        "> {log.path} 2>&1 && "
        "rm -rf {params.tmpdir} "