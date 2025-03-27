# -*- coding: utf-8 -*-

import yaml

"""
Map reads to specified reference genomes.
"""

# For improved speed, add 'usejni=t' to the command line of bbmap tools which support the use of the compiled jni C code.

# Parameters from the configuration that we'll use
mapping_spec = config['mapping_spec']               # File with tax IDs and the reference genomes to map them against
fastq_file_pattern = config["fastq_file_pattern"]   # The pattern of the fastq files

# Load mapping specification
def load_mapping_spec(mapping_spec):
    with open(mapping_spec) as f:
        return yaml.safe_load(f)

mapping_spec = load_mapping_spec(mapping_spec)

# Define input FASTQ files (assuming a consistent naming scheme)
tax_ids = mapping_spec.keys()

def fastq_files(wildcards):
    return [fastq_file_pattern.format(tax_id=wildcards.tax_id, direction="1"), fastq_file_pattern.format(tax_id=wildcards.tax_id, direction="2")]

def index_dir(ref):
    return os.path.splitext(os.path.basename(ref))[0]

# Rule to build BBMap indexes
rule build_index:
    input:
        fasta=GENOMESDIR/"{ref}"
    output:
        index=lambda wildcards: INDEXDIR + "/" + index_dir(wildcards.ref)
    shell:
        "bbmap.sh ref={input.fasta} path={output.index}"

# Rule to map reads using BBMap
rule map_reads:
    input:
        reads=fastq_files,
        index=lambda wildcards: INDEXDIR + "/" + index_dir(wildcards.ref)
    output:
        bam=OUTDIR/"mapped/{tax_id}_{ref}.bam"
    shell:
        "bbmap.sh in={input.reads[0]} in2={input.reads[1]} ref={input.index} out={output.bam}"

# Generate rules dynamically for each TAX_ID and reference pair
rule all:
    input:
        expand("mapped/{{tax_id}}_{{ref}}.bam", tax_id=tax_ids, ref=[ref for tax_id in tax_ids for ref in mapping_spec[tax_id]])


# # Create a mapping dict between tax IDs and reference genomes to map against
# def read_mapping_spec(file):
#     with open(file) as f:
#         return yaml.safe_load(f)

# # Get the mapping dict
# taxid2ref = read_mapping_spec(mapping_spec)

# # Fail-safe: Exit if the mapping spec file was empty
# if not taxid2ref:
#     raise ValueError(f"Nothing found in {mapping_spec}. Check the file contents.")

# # The fastq output files containing taxonomic ID-specific reads
# extracted_reads = expand(
#     str(OUTDIR/"mode_{mode}/{tax_id}/{sample}_taxID-{tax_id}_{direction}.fq.gz"),
#     mode=mode,
#     sample=SAMPLES,
#     tax_id=tax_ids,
#     direction=['R1', 'R2'])
# all_outputs.extend(extracted_reads)

# rule make_index:
#     conda:
#         "../envs/map.yaml"
#     input:
#         reference = 

# rule getAll_taxID_readIDs:
#     conda:
#         "../envs/extract_reads.yaml"
#     input:
#         kraken2 = lambda wildcards: INPUTDIR/k2_classification_file.format(sample=wildcards.sample)
#     output:
#         temp(OUTDIR/"mode_{mode}/{sample}_allTaxIDs_readIDs.txt")
#     params:
#         names = names_file,
#         nodes = nodes_file,
#         extract_mode = mode,
#         tax_id_file = tax_ids_file
#     log:
#         LOGDIR / "getAll_taxID_readIDs_{sample}_{mode}.log"
#     shell:
#         """
#         /usr/bin/time -v scripts/selmeout.py \
#             --mode {params.extract_mode} \
#             --nodes {params.nodes} \
#             --names {params.names} \
#             --output {output} \
#             --input {input.kraken2} \
#             --tax_id_file {params.tax_id_file} \
#             2>&1 | tee {log}
#         """

# rule extract_taxID_reads:
#     conda:
#         "../envs/extract_reads.yaml"
#     input:
#         R1 = lambda wildcards: INPUTDIR/fastq_file_pattern.format(sample=wildcards.sample, direction='R1'),
#         R2 = lambda wildcards: INPUTDIR/fastq_file_pattern.format(sample=wildcards.sample, direction='R2'),
#         allReadIDs = OUTDIR/"mode_{mode}/{sample}_allTaxIDs_readIDs.txt"
#     output:
#         taxon_readIDs = temp(OUTDIR/"mode_{mode}/{tax_id}/{sample}_taxID-{tax_id}_readIDs.txt"),
#         R1 = OUTDIR/"mode_{mode}/{tax_id}/{sample}_taxID-{tax_id}_R1.fq.gz",
#         R2 = OUTDIR/"mode_{mode}/{tax_id}/{sample}_taxID-{tax_id}_R2.fq.gz"
#     log:
#         LOGDIR / "extract_taxID_reads_{sample}_taxID-{tax_id}_{mode}.log"
#     shell:
#         """
#         awk -v taxid={wildcards.tax_id} -F '\t' '$2 == taxid {{print $3}}' {input.allReadIDs} > {output.taxon_readIDs}
#         /usr/bin/time -v sh -c 'seqtk subseq {input.R1} {output.taxon_readIDs} | gzip - > {output.R1}' 2>> {log}
#         /usr/bin/time -v sh -c 'seqtk subseq {input.R2} {output.taxon_readIDs} | gzip - > {output.R2}' 2>> {log}
#         """
