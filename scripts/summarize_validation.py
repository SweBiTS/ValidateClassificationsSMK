import sys
import logging
import pandas as pd
import pysam
from pathlib import Path
from collections import defaultdict

# --- Snakemake objects ---
try:
    snk_input = snakemake.input
    snk_output = snakemake.output
    snk_log = snakemake.log
    snk_wildcards = snakemake.wildcards
except AttributeError as e:
    raise AttributeError(f"Missing expected Snakemake object attribute: {e}")

# --- File Paths ---
overlapping_bam_path = snk_input.overlapping_bam
classifiable_bed_path = snk_input.classifiable_bed
summary_tsv_path = snk_output.summary_tsv
log_file = snk_log.path

# --- Setup Logging ---
Path(log_file).parent.mkdir(parents=True, exist_ok=True)
logging.basicConfig(filename=log_file, filemode='w', level=logging.INFO,
                    format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

logging.info("--- Snakemake Parameters ---")
logging.info(f"Overlapping Reads BAM: {overlapping_bam_path}")
logging.info(f"Classifiable Regions BED: {classifiable_bed_path}")
logging.info(f"Output Summary TSV: {summary_tsv_path}")
logging.info(f"Log File: {log_file}")
logging.info(f"Wildcards: {snk_wildcards}")
logging.info("----------------------------")

# --- 1. Calculate Length of Classifiable Region per Contig ---
logging.info(f"Calculating classifiable region length from {classifiable_bed_path}")
contig_lengths = defaultdict(int)
try:
    bed_df = pd.read_csv(classifiable_bed_path, sep='\t', header=None,
                         usecols=[0, 1, 2], names=['chrom', 'start', 'end'],
                         dtype={'chrom': str, 'start': int, 'end': int})

    if not bed_df.empty:
        bed_df['length'] = bed_df['end'] - bed_df['start']
        contig_lengths = bed_df.groupby('chrom')['length'].sum().to_dict()
    logging.info(f"Calculated lengths for {len(contig_lengths)} contigs from BED.")

except FileNotFoundError:
    logging.error(f"Classifiable regions BED file not found: {classifiable_bed_path}")
    sys.exit(1)
except pd.errors.EmptyDataError:
     logging.warning(f"Classifiable regions BED file is empty: {classifiable_bed_path}")
except Exception as e:
    logging.error(f"Error processing classifiable regions BED {classifiable_bed_path}: {e}", exc_info=True)
    sys.exit(1)

# --- 2. Count UNIQUE Overlapping Read PAIRS per Contig ---
logging.info(f"Counting unique read pairs in {overlapping_bam_path}")
contig_read_name_sets = defaultdict(set)
all_contigs_in_bam = set()
contig_read_counts = defaultdict(int)
try:
    with pysam.AlignmentFile(overlapping_bam_path, "rb") as bamfile:
        # Get the definitive list of contigs from the BAM header
        if bamfile.header.references:
             all_contigs_in_bam = set(bamfile.header.references)
        # Initialize sets for all contigs found in BAM 
        for contig_name in all_contigs_in_bam:
             contig_read_name_sets[contig_name] = set()

        # Iterate through alignments to populate sets
        processed_alignments = 0
        for alignment in bamfile.fetch(until_eof=True):
            if not alignment.is_unmapped and alignment.query_name:
                contig_name = alignment.reference_name
                if contig_name:
                    contig_read_name_sets[contig_name].add(alignment.query_name)
                    processed_alignments += 1

        # Calculate final counts from set sizes
        for contig, read_names in contig_read_name_sets.items():
            contig_read_counts[contig] = len(read_names)
        
        total_unique_reads = sum(contig_read_counts.values())
        logging.info(f"Processed {processed_alignments} alignments, representing {total_unique_reads} unique read pairs overlapping classifiable regions across {len(contig_read_counts)} contigs.")

except FileNotFoundError:
    logging.error(f"Overlapping reads BAM file not found: {overlapping_bam_path}")
    sys.exit(1)
except ValueError as e:
     logging.error(f"Error opening or reading BAM file {overlapping_bam_path}: {e}", exc_info=True)
     sys.exit(1)
except Exception as e:
    logging.error(f"Error processing overlapping reads BAM {overlapping_bam_path}: {e}", exc_info=True)
    sys.exit(1)

# --- 3. Combine and Write Summary ---
logging.info(f"Combining results and writing to {summary_tsv_path}")
try:
    if not all_contigs_in_bam:
         logging.warning("No contigs found in BAM header. Output table will be empty.")

    summary_data = []
    for contig in sorted(list(all_contigs_in_bam)):
        summary_data.append({
            "Contig": contig,
            "NumReadsOverlapClassifiable": contig_read_counts.get(contig, 0),
            "LengthClassifiableRegion": contig_lengths.get(contig, 0)
        })

    summary_df = pd.DataFrame(summary_data)
    summary_df = summary_df[[
        "Contig", "NumReadsOverlapClassifiable", "LengthClassifiableRegion"
    ]]

    # Ensure output directory exists
    Path(summary_tsv_path).parent.mkdir(parents=True, exist_ok=True)
    summary_df.to_csv(summary_tsv_path, sep='\t', index=False, header=True)
    logging.info(f"Successfully wrote summary TSV: {summary_tsv_path}")

except Exception as e:
    logging.error(f"Failed to combine results or write summary TSV: {e}", exc_info=True)
    sys.exit(1)

logging.info("Summary script finished successfully.")