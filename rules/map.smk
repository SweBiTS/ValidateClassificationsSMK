import yaml
import sys
from pathlib import Path
import logging

# Logger specific to this module - inherits config from root logger
logger = logging.getLogger(__name__)

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


# --- Helper Function: Find Actual Genome File Path (.fasta or .fa) ---
def get_actual_fasta_path(wildcards):
    """Finds the full path to the genome file, checking .fasta and .fa extensions."""
    basename = wildcards.genome_basename
    fasta_path = GENOMES_DIR_P / f"{basename}.fasta"
    fa_path = GENOMES_DIR_P / f"{basename}.fa"

    if fasta_path.exists():
        return str(fasta_path)
    elif fa_path.exists():
        return str(fa_path)
    else:
        raise FileNotFoundError(
            f"Genome file for basename '{basename}' not found with .fasta or .fa extension "
            f"in directory '{GENOMES_DIR_P}'"
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
BAM_INDEX_OUT_PATTERN = "{outdir}/mapping/{tax_id}/{genome_basename}/mapping.sorted.bam.bai"
STATS_OUT_PATTERN = "{outdir}/mapping/{tax_id}/{genome_basename}/mapping_stats.txt"
LOG_MAP_PATTERN = "{logdir}/bbmap_map/{tax_id}_{genome_basename}.log"
LOG_SORT_PATTERN = "{logdir}/samtools_sort/{tax_id}_{genome_basename}.log"
LOG_INDEX_PATTERN = "{logdir}/samtools_index/{tax_id}_{genome_basename}.log"

# --- Define Functions for Rule All ---
# Targets are indices of sorted bam files (last step)
def get_all_target_indices():
     return get_all_target_outputs(mapping_spec_data, BAM_INDEX_OUT_PATTERN)


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
        path = str(LOG_DIR_P / "bbmap_index" / "{genome_basename}.log")
    conda:
        ".." / ENVS_DIR_P / "map.yaml"
    threads: config.get("BBMAP_INDEX_THREADS", 8)
    shell:
        "bbmap.sh "
            "ref={input.fasta} "
            "path={output.idx_dir} "
            "build={params.build_id} "
            "-Xmx{params.mem} "
            "threads={threads} "
            "pigz=t unpigz=t "
            "overwrite=t "
            "usejni=t "
            "> {log} 2>&1"


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
        mem = config.get("BBMAP_MAP_MEM", "20g"),
        minid = config["BBMAP_MINID"],
        ambig = config["BBMAP_AMBIGUOUS"],
        pairedonly = config["BBMAP_PAIREDONLY"],
        extra_stats = lambda wildcards, output: f"statsfile={output.stats}"
    log:
        path = LOG_MAP_PATTERN.format(
            logdir=LOG_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    conda:
        ".." / ENVS_DIR_P / "map.yaml"
    threads: config.get("BBMAP_MAP_THREADS", 8)
    shell:
        # Explicitly define the path to the jni directory using the job's conda env prefix
        # NOTE: Double-check the relative path 'opt/bbmap-39.01-1/jni' if bbmap version changes or install differs.
        #"export JNI_DIR={job.conda_env_prefix}/opt/bbmap-39.01-1/jni && "
        # Prepend this directory to the LD_LIBRARY_PATH
        #"export LD_LIBRARY_PATH=$JNI_DIR:$LD_LIBRARY_PATH && "
        # Add optional debug output to SLURM log (remove later)
        #"echo 'DEBUG: JNI_DIR=$JNI_DIR' >&2 && "
        #"echo 'DEBUG: LD_LIBRARY_PATH=$LD_LIBRARY_PATH' >&2 && "
        "bbmap.sh "
            "in1={input.r1} "
            "in2={input.r2} "
            "path={input.idx_dir} "
            "build={params.build_id} "
            "out={output.bam} "
            "{params.extra_stats} "
            "minid={params.minid} "
            "ambiguous={params.ambig} "
            "pairedonly={params.pairedonly} "
            "-Xmx{params.mem} "
            "threads={threads} "
            "pigz=t unpigz=t "
            "overwrite=t "
            "usejni=t "
            "> {log} 2>&1"


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
        extra = "" # Placeholder for extra samtools sort options if needed
    log:
        path = LOG_SORT_PATTERN.format(
            logdir=LOG_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    conda:
        ".." / ENVS_DIR_P / "map.yaml"
    threads: config.get("SAMTOOLS_SORT_THREADS", 4)
    shell:
        "samtools sort "
        "-@ {threads} "                 # Number of sorting threads
        "-m {params.mem_per_thread} "   # Memory per thread
        "-o {output.bam} "              # Output file
        "{params.extra} "               # Extra options
        "{input.bam} "                  # Input file
        "2> {log}"                      # Log stderr


# --- Rule: Index sorted BAM file using Samtools ---
rule samtools_index:
    input:
        bam = SORTED_BAM_OUT_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    output:
        index = BAM_INDEX_OUT_PATTERN.format(
            outdir=OUTPUT_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    log:
        path = LOG_INDEX_PATTERN.format(
            logdir=LOG_DIR_P, tax_id="{tax_id}", genome_basename="{genome_basename}"
        )
    conda:
        ".." / ENVS_DIR_P / "map.yaml"
    threads: config.get("SAMTOOLS_INDEX_THREADS", 1) # Indexing often doesn't scale well with threads
    shell:
        "samtools index "
        "-@ {threads} "  # Number of threads
        "{input.bam} "   # Input file (index is created alongside as {input.bam}.bai)
        "2> {log}"       # Log stderr