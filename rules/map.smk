from pathlib import Path


# ========================== #
# --- FILE PATH PATTERNS --- #
# ========================== #

# --- output patterns --- #
MAPPING_OUT_DIR = OUTPUT_DIR / "mapping"
RAW_BAM_OUT_PATTERN = str(MAPPING_OUT_DIR / "{tax_id}/{genome_basename}/mapping.bam")
SORTED_BAM_OUT_PATTERN = str(MAPPING_OUT_DIR / "{tax_id}/{genome_basename}/mapping.sorted.bam")
DEDUP_BAM_OUT_PATTERN = str(MAPPING_OUT_DIR / "{tax_id}/{genome_basename}/mapping_dedup.bam")
DEDUP_BAM_INDEX_OUT_PATTERN = str(MAPPING_OUT_DIR / "{tax_id}/{genome_basename}/mapping_dedup.bam.bai")
STATS_OUT_PATTERN = str(MAPPING_OUT_DIR / "{tax_id}/{genome_basename}/mapping_stats.txt")
COVERAGE_OUTPUT_PATTERN = str(MAPPING_OUT_DIR / "{tax_id}/{genome_basename}/coverage.tsv")

# --- log patterns --- #
MAPPING_LOG_DIR = LOG_DIR / "mapping"
LOG_BBMAP_INDEX_PATTERN = str(MAPPING_LOG_DIR / "{genome_basename}__bbmap_index.log")
LOG_MAP_PATTERN = str(MAPPING_LOG_DIR / "{tax_id}_{genome_basename}__bbmap_map.log")
LOG_SORT_PATTERN = str(MAPPING_LOG_DIR / "{tax_id}_{genome_basename}__samtools_sort.log")
LOG_MARKDUP_PATTERN = str(MAPPING_LOG_DIR / "{tax_id}_{genome_basename}__sambamba_markdup.log")
LOG_CALC_COV_PATTERN = str(MAPPING_LOG_DIR / "{tax_id}_{genome_basename}__coverage.log")

# ==================== #
# --- HELPER FUNCS --- #
# ==================== #

# LÄGG TILL /usr/bin/time COMMAND AT BEGINNING OF SHELL COMMANDS (AND IN SCRIPTS)
# f"{BBMAP_JNI_FLAG} "  # Adds 'usejni=t' or 'usejni=f' <-- för att lägga inn i shell commands
# Ändra output format så att man har en tabell, och ändra output file
# då kan vi lätt samla ihop allt i slutet

# --- Find Actual Genome File Path ---
def get_actual_fasta_path(wildcards):
    """Finds the full path to the genome file, checking common extensions."""
    basename = wildcards.genome_basename
    genomes_dir = GENOMES_DIR
    
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

def get_all_coverage_files(wildcards):
    """Gets list of all per-sample coverage files. For collection rule at the end + checkpoint"""
    return get_target_outputs(MAPPING_SPEC_DATA, COVERAGE_OUTPUT_PATTERN)

# ============= #
# --- RULES --- #
# ============= #

localrules: collect_mapping_branch

# --- Index Reference Genome using BBMap --- #
rule bbmap_index:
    input:
        fasta = get_actual_fasta_path
    output:
        idx_dir = directory(BBMAP_INDEX_DIR / "{genome_basename}")
    params:
        build_id = 1,
        mem = config.get("BBMAP_INDEX_MEM", "20g")
    log:
        path = LOG_BBMAP_INDEX_PATTERN.format(
            genome_basename="{genome_basename}")
    conda:
        ".." / ENVS_DIR / "map.yaml"
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

# --- Map Reads to Reference using BBMap --- #
rule bbmap_map_reads:
    input:
        idx_dir = rules.bbmap_index.output.idx_dir,
        r1 = lambda wildcards: str(
            INPUT_DIR / config["FASTQ_FILE_PATTERN"].format(
                tax_id=wildcards.tax_id, direction="1")),
        r2 = lambda wildcards: str(
            INPUT_DIR / config["FASTQ_FILE_PATTERN"].format(
                tax_id=wildcards.tax_id, direction="2"))
    output:
        bam = RAW_BAM_OUT_PATTERN.format(
            tax_id="{tax_id}", genome_basename="{genome_basename}"),
        stats = STATS_OUT_PATTERN.format(
            tax_id="{tax_id}", genome_basename="{genome_basename}")
    params:
        build_id = 1,
        mem = config.get("BBMAP_MAP_MEM", "6g"),
        minid = config["BBMAP_MINID"],
        ambig = config["BBMAP_AMBIGUOUS"],
        pairedonly = config["BBMAP_PAIREDONLY"]
    log:
        path = LOG_MAP_PATTERN.format(
            tax_id="{tax_id}", genome_basename="{genome_basename}")
    conda:
        ".." / ENVS_DIR / "map.yaml"
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

# --- Sort BAM file using Samtools --- #
rule samtools_sort:
    input:
        bam = RAW_BAM_OUT_PATTERN.format(
            tax_id="{tax_id}", genome_basename="{genome_basename}")
    output:
        bam = SORTED_BAM_OUT_PATTERN.format(
            tax_id="{tax_id}", genome_basename="{genome_basename}")
    params:
        # Memory per thread for sorting
        mem_per_thread = config.get("SAMTOOLS_SORT_MEM_PER_THREAD", "2G")
    log:
        path = LOG_SORT_PATTERN.format(
            tax_id="{tax_id}", genome_basename="{genome_basename}")
    conda:
        ".." / ENVS_DIR / "map.yaml"
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

# --- Mark/Remove Duplicates using Sambamba --- #
rule mark_duplicates_sambamba:
    input:
        # Input is the coordinate-sorted BAM from samtools_sort
        sorted_bam = SORTED_BAM_OUT_PATTERN.format(
            tax_id="{tax_id}", genome_basename="{genome_basename}")
    output:
        # Output is the BAM file with duplicates removed
        dedup_bam = DEDUP_BAM_OUT_PATTERN.format(
            tax_id="{tax_id}", genome_basename="{genome_basename}"),
        index = DEDUP_BAM_INDEX_OUT_PATTERN.format(
            tax_id="{tax_id}", genome_basename="{genome_basename}")
    params:
        # Define a temporary directory relative to the output
        tmpdir = lambda wildcards, output: Path(output.dedup_bam).parent / "tmp_sambamba_markdup"
    log:
        path = LOG_MARKDUP_PATTERN.format(
            tax_id="{tax_id}", genome_basename="{genome_basename}")
    conda:
         ".." / ENVS_DIR / "map.yaml"
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

# --- Calculate Real Read Coverage using Samtools --- #
rule calculate_coverage:
    input:
        # Use the final deduplicated BAM
        bam = DEDUP_BAM_OUT_PATTERN.format(
            tax_id="{tax_id}", genome_basename="{genome_basename}")
    output:
        # File containing coverage
        cov_file = COVERAGE_OUTPUT_PATTERN.format(
            tax_id="{tax_id}", genome_basename="{genome_basename}")
    log:
        path = LOG_CALC_COV_PATTERN.format(
            tax_id="{tax_id}", genome_basename="{genome_basename}")
    conda:
         ".." / ENVS_DIR / "map.yaml"
    threads: config.get("CALC_REAL_COV_THREADS", 1)
    resources:
        runtime=config.get("CALC_REAL_COV_RUNTIME", "15m"),
        mem_mb=config.get("CALC_REAL_COV_MEM_MB", 1000),
        cpus_per_task=config.get("CALC_REAL_COV_THREADS", 1)
    shell:
        # Use a samtools/awk pipe to calculate coverage
        r"""
        mkdir -p $(dirname {output.cov_file}) &&
        mkdir -p $(dirname {log.path}) &&
        {{
            coverage=$(samtools depth -a {input.bam} | awk '{{sum+=$3; count++}} END {{if (count > 0) print sum/count; else print 0}}') &&
            echo -e {wildcards.tax_id}\\t{wildcards.genome_basename}\\t${{coverage}} > {output.cov_file} ;
        }} 2> {log.path}
        """

# --- Collect Coverage Results --- #
rule collect_mapping_branch:
    input:
        get_all_coverage_files
    output:
        MAPPING_BRANCH_TARGET
    shell:
        "touch {output}"