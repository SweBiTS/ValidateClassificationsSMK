from pathlib import Path

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
            genome_basename="{genome_basename}")
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
        r1 = lambda wildcards: str(
            INPUT_DIR_P / config["FASTQ_FILE_PATTERN"].format(
                tax_id=wildcards.tax_id, direction="1")),
        r2 = lambda wildcards: str(
            INPUT_DIR_P / config["FASTQ_FILE_PATTERN"].format(
                tax_id=wildcards.tax_id, direction="2"))
    output:
        bam = temp(RAW_BAM_OUT_PATTERN.format(
            tax_id="{tax_id}", genome_basename="{genome_basename}")),
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
            tax_id="{tax_id}", genome_basename="{genome_basename}")
    output:
        bam = temp(SORTED_BAM_OUT_PATTERN.format(
            tax_id="{tax_id}", genome_basename="{genome_basename}"))
    params:
        # Memory per thread for sorting
        mem_per_thread = config.get("SAMTOOLS_SORT_MEM_PER_THREAD", "2G")
    log:
        path = LOG_SORT_PATTERN.format(
            tax_id="{tax_id}", genome_basename="{genome_basename}")
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