# Transcriptome and GRN Complexity Analysis

**Author:** Xoel Mato Blanco

This folder contains analyses related to transcriptome and Gene Regulatory Network (GRN) complexity.

## Analysis Steps

1.  **Data Preparation (`0_get_counts.ipynb`):**
    *   Load main AnnData object (`.h5ad`).
    *   Locate and load metadata from CellOracle output files.
    *   Extract and save count data corresponding to the metadata for specific conditions/samples into CSV files.

2.  **Complexity Calculation and Plotting (`1_transcriptome_and_grn_complexity.ipynb`):**
    *   **Setup:** Source configuration, set up plotting parameters and colors, detect cores.
    *   **Load TFs:** Read a list of transcription factors.
    *   **GRN Complexity:**
        *   Load pre-computed network scores.
        *   Analyze and plot network scores (e.g., degree, centrality) for all genes, source genes, and TF genes.
        *   Analyze and plot the distribution of network roles and modules.
    *   **Transcriptomic Complexity (scRNA-seq):**
        *   Load count matrices and metadata for different samples.
        *   Calculate per-cell complexity metrics:
            *   Number of Genes Expressed
            *   Shannon Entropy
            *   Simpson's Diversity Index
        *   Calculate per-gene complexity metric within cell types:
            *   Coefficient of Variation (CV)
        *   Plot comparisons of these scRNA-seq complexity metrics across samples and cell types for all genes and TF genes separately.

## Scripts / Notebooks

*   `0_get_counts.ipynb`: Prepares and saves count data subsets based on CellOracle metadata.
*   `1_transcriptome_and_grn_complexity.ipynb`: Calculates and visualizes GRN and transcriptomic complexity metrics.

## Related Figures

*   Supplementary Figure 12
