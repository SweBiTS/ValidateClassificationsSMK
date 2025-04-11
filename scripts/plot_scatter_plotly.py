# scripts/plot_scatter_plotly.py

import sys
import logging
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path
from scipy import stats
import warnings

# --- Snakemake objects ---
try:
    snk_input = snakemake.input
    snk_output = snakemake.output
    snk_log = snakemake.log
    snk_wildcards = snakemake.wildcards
except AttributeError as e:
    raise AttributeError(f"Missing expected Snakemake object attribute: {e}")

# --- File Paths ---
summary_tsv_path = snk_input.summary_tsv
plot_html_path = snk_output.plot_html
log_file = snk_log.path

# --- Setup Logging ---
Path(log_file).parent.mkdir(parents=True, exist_ok=True)
logging.basicConfig(filename=log_file, filemode='w', level=logging.INFO,
                    format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

logging.info("--- Snakemake Parameters ---")
logging.info(f"Input TSV: {summary_tsv_path}")
logging.info(f"Output Plot HTML: {plot_html_path}")
logging.info(f"Log File: {log_file}")
logging.info(f"Wildcards: {snk_wildcards}")
logging.info("----------------------------")

# --- Plotting ---
try:
    logging.info(f"Reading data from {summary_tsv_path}")
    try:
        df = pd.read_csv(summary_tsv_path, sep='\t')
    except pd.errors.EmptyDataError:
        logging.error(f"Input TSV file is empty: {summary_tsv_path}")
        sys.exit(1)

    required_columns = {'LengthClassifiableRegion', 'NumReadsOverlapClassifiable', 'Contig'}
    if df.empty or not required_columns.issubset(df.columns):
        logging.error(f"Input TSV is empty or missing required columns {required_columns}. Cannot generate plot.")
        sys.exit(1)

    logging.info(f"Generating interactive scatter plot for {len(df)} contigs.")

    # --- Calculate Correlations ---
    pearson_text = "Pearson r: "
    spearman_text = "Spearman ρ: "
    # Filter out rows where length is 0 or NaN for meaningful correlation
    df_filt = df.dropna(subset=['LengthClassifiableRegion', 'NumReadsOverlapClassifiable'])
    df_filt = df_filt[df_filt['LengthClassifiableRegion'] > 0]

    if len(df_filt) > 1:
        # Calculate Pearson
        try:
            pearson_r, p_val_p = stats.pearsonr(df_filt['LengthClassifiableRegion'], df_filt['NumReadsOverlapClassifiable'])
            if pd.isna(pearson_r):
                pearson_text += "N/A (NaN)"
                logging.warning("Pearson correlation result is NaN.")
            else:
                pearson_text += f'{pearson_r:.3f} (p={p_val_p:.2g})'
        except ValueError as e: # Catch errors like zero variance
            logging.warning(f"Could not calculate Pearson correlation: {e}")
            pearson_text += "N/A (Error)"

        # Calculate Spearman
        try:
            spearman_r, p_val_s = stats.spearmanr(df_filt['LengthClassifiableRegion'], df_filt['NumReadsOverlapClassifiable'])
            if pd.isna(spearman_r):
                spearman_text += "N/A (NaN)"
                logging.warning("Spearman correlation result is NaN.")
            else:
                spearman_text += f'{spearman_r:.3f} (p={p_val_s:.2g})'
        except ValueError as e: # Catch errors like zero variance
            logging.warning(f"Could not calculate Spearman correlation: {e}")
            spearman_text += "N/A (Error)"

        logging.info(f"Correlation Results: {pearson_text} | {spearman_text}")

    else: # Handle case with too few points for correlation
        pearson_text += "N/A (too few points)"
        spearman_text += "N/A (too few points)"
        logging.info("Correlation not calculated (too few points with length > 0).")

    # --- Create Plot ---
    plot_title = f"Reads vs Classifiable Length<br>TaxID: {snk_wildcards.tax_id} | Genome: {snk_wildcards.genome_basename}<br>{pearson_text} | {spearman_text}"

    fig = px.scatter(df, # Use original dataframe for plotting all points
                     x="LengthClassifiableRegion",
                     y="NumReadsOverlapClassifiable",
                     title=plot_title,
                     labels={
                         "LengthClassifiableRegion": "Length of Classifiable Region (bp)",
                         "NumReadsOverlapClassifiable": "Number of Overlapping Real Reads (Pairs)"
                     },
                     hover_data=['Contig']
                    )

    fig.update_layout(title_font_size=14, title_x=0.5)
    fig.update_traces(marker=dict(size=8, opacity=0.7), selector=dict(mode='markers'))

    # --- Add Buttons for Axis Scale Toggle ---
    fig.update_layout(
        updatemenus=[
            dict(
                type="buttons",
                direction="left",
                buttons=list([
                    dict(args=[{"xaxis.type": "linear", "yaxis.type": "linear"}], label="Linear Scale", method="relayout"),
                    dict(args=[{"xaxis.type": "log", "yaxis.type": "linear"}], label="Log X", method="relayout"),
                    dict(args=[{"xaxis.type": "linear", "yaxis.type": "log"}], label="Log Y", method="relayout"),
                    dict(args=[{"xaxis.type": "log", "yaxis.type": "log"}], label="Log X & Y", method="relayout")
                ]),
                pad={"r": 10, "t": 10}, showactive=True,
                x=0.1, xanchor="left", y=1.1, yanchor="top" # Position buttons
            ),
        ]
    )

    # --- Save Plot ---
    Path(plot_html_path).parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(plot_html_path)
    logging.info(f"Interactive plot saved to {plot_html_path}")

except FileNotFoundError:
    logging.error(f"Input TSV file not found: {summary_tsv_path}")
    sys.exit(1)
except Exception as e:
    logging.error(f"An error occurred during plotting: {e}", exc_info=True)
    sys.exit(1)

logging.info("Plotting script finished successfully.")