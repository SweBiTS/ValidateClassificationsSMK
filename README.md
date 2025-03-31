# ValidateClassificationsSMK

Map taxon-specific reads (e.g., previously classified by Kraken 2) to specified reference genomes using BBMap and Samtools.

This pipeline was created to run on a cluster with the SLURM Workload Manager and Conda package manager. Therefore this README is tailored towards that kind of environment.

Currently this is tailored to work on the PDC cluster Dardel. You just need to update the SLURM details in `slurm/config.yaml` if you want to use this pipeline on another cluster.

## Workflow Overview

The pipeline performs the following main steps for each taxonomic identifier (TaxID) and associated reference genome(s) defined in the `supporting_files/mapping_specification.yaml` file:

1.  **Index Reference Genome:** Creates a BBMap index for the specified reference genome FASTA file. The indexes are stored in `output/{FOLDER}`, where `FOLDER` is the folder specified by the `INDEX_DIR_BBMAP` variable in the `config.yaml` file (default: `output/bbmap_indices`).
2.  **Map Reads:** Maps the corresponding input FASTQ reads (paired-end) against the indexed reference genome using `bbmap.sh`. Generates a BAM file and mapping statistics.
3.  **Sort BAM:** Sorts the resulting BAM file by coordinate using `samtools sort`. The initial BAM file from mapping is marked as temporary and removed automatically upon successful sorting to save disk space.
4.  **Index BAM:** Creates a standard BAM index (`.bai`) for the sorted BAM file using `samtools index`. The final outputs are the sorted, indexed BAM files and the mapping statistics files.

## Usage

### Before you run the pipeline

#### 1. Snakemake Installation (v9+ Recommended)
This pipeline utilizes features available in Snakemake v8+ (and tested with v9+). Ensure you have a compatible version installed, along with the SLURM executor plugin and PyYAML.

You can create a dedicated conda environment for Snakemake:
```bash
conda create -c conda-forge -c bioconda -n snakemake snakemake snakemake-executor-plugin-slurm pyyaml
```
*(Note: `pyyaml` is needed for parsing the `mapping_specification.yaml` file within the Snakefile).*

Activate the environment before running the pipeline:
```bash
conda activate snakemake
```

#### 2. Reference Genomes
Place the reference genome files (in FASTA format, with `.fasta` or `.fa` extension) that you intend to map reads against into the directory specified by `GENOMES_DIR` in `config.yaml` (default: `supporting_files/genomes/`).

#### 3. Mapping Specification File (`mapping_specification.yaml`)
This crucial YAML file tells the pipeline *which* reads (identified by TaxID) should be mapped to *which* reference genome(s).
* Configure the path to this file using `MAPPING_SPECIFICATION` in `config.yaml` (default: `supporting_files/mapping_specification.yaml`).
* **Format:** A YAML dictionary.
    * **Keys:** Taxonomic IDs (as numbers or strings) corresponding to your read sets.
    * **Values:** A Python *list* containing one or more reference genome *filenames* (these filenames must exist within your `GENOMES_DIR`).

```yaml
# Example mapping_specification.yaml
33189:
  - pt_128_001_nuclear_20220627.fasta # Map reads for taxid 33189 to this genome
9606: # Example: Map Human reads to two different genome versions/files
  - hg38_main_chroms.fasta
  - hg19_main_chroms.fasta
10090:
  - mm10_primary_assembly.fa
```

#### 4. Input FASTQ Files
Place the FASTQ files containing the classified reads into the `input/fastq` directory. Symlinks are acceptable.
* The pipeline expects paired-end reads (`_R1` and `_R2`).
* The filenames **must** match the pattern defined by `FASTQ_FILE_PATTERN` in `config.yaml`. The default pattern is `taxID-{tax_id}_R{direction}.fq.gz`.
* The `{tax_id}` portion extracted from the filename **must** match a key present in your `mapping_specification.yaml`.
* These FASTQ files are typically the output of an upstream process, such as the ReadExtractorSMK workflow that groups reads by taxonomic classification.

#### 5. Workflow & Cluster Configuration
1.  **Workflow Configuration (`config.yaml`):** Review the main `config.yaml` file. Check and adjust all paths (`INPUT_DIR`, `OUTPUT_DIR`, `LOG_DIR`, `GENOMES_DIR`, `INDEX_DIR_BBMAP`, `ENVS_DIR`), file patterns (`FASTQ_FILE_PATTERN`), the mapping specification path (`MAPPING_SPECIFICATION`), and tool parameters (`BBMAP_*`).
2.  **SLURM Configuration (`slurm/config.yaml`):** Review the SLURM configuration file carefully.
    * Adjust the maximum concurrent `jobs`.
    * Set appropriate `default-resources` (runtime, memory, CPUs) as a baseline.
    * Define specific resources under `set-resources` for rules like `bbmap_index`, `bbmap_map_reads`, and `samtools_sort` based on your expected data sizes and cluster capabilities. **Resource requirements likely need tuning to your use case.**
    * **Crucially, update `slurm_account` and `slurm_partition` to match your specific SLURM cluster setup.**

### Run the pipeline

1.  Activate your Snakemake conda environment:
    ```bash
    conda activate snakemake
    ```
2.  Navigate to the pipeline's root directory (where the main `Snakefile` is located).
3.  **(Recommended)** Perform a dry run first to check the execution plan and detect configuration issues:
    ```bash
    snakemake -n --executor slurm --profile slurm --use-conda --cores 1
    ```
    *(Note: `-n` is short for `--dry-run`)*

4.  Execute the pipeline:
    ```bash
    snakemake --executor slurm --profile slurm --use-conda --cores 1
    ```

    * `--executor slurm`: Use the SLURM executor plugin.
    * `--profile slurm`: Load cluster configuration from the `slurm` directory/profile. Adjust the profile name/path if needed. (Alternatively, use `--cluster-config` and `--cluster` flags if not using profiles).
    * `--use-conda`: Create/use conda environments defined in the rules.
    * `--cores 1`: Cores for the main Snakemake process (usually 1 is fine). Job cores are defined in `slurm/config.yaml`.
    * *(Hint: For long-running workflows on remote systems, consider using terminal multiplexers like `screen` or `tmux` to keep the main Snakemake process running even if your connection is interrupted).*
