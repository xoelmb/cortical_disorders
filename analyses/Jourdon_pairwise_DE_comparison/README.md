# Pairwise Differential Expression Comparison in Jourdon et al. Data

**Author:** Xoel Mato Blanco

This folder contains analyses comparing gene expression variability within and between ASD (Autism Spectrum Disorder) and control donor groups using pairwise differential expression (DE) analysis on the Jourdon et al. 2023 dataset, focusing on Radial Glia (RG) cells.

## Analysis Steps

1.  **Setup & Data Loading (`1_pairwise_de.ipynb`, `2_expression_plots.ipynb`):**
    *   Load necessary R libraries (`Seurat`, `dplyr`, `ggplot2`, etc.).
    *   Define input/output paths and helper functions for plotting.
    *   Load the main Seurat object from Jourdon et al. 2023.
    *   Load supplementary data (TF lists, previous study DEGs).

2.  **RG Subset & Preprocessing (`1_pairwise_de.ipynb`, `2_expression_plots.ipynb`):**
    *   Subset the Seurat object to retain only Radial Glia (RG) cells.
    *   Perform standard preprocessing on the RG subset: SCTransform normalization, PCA, UMAP, neighbor finding, and clustering.
    *   Generate initial quality control plots: bar plots of cell counts per donor/condition, UMAP visualizations colored by metadata.

3.  **Pairwise Differential Expression (`1_pairwise_de.ipynb`):**
    *   Define all unique pairwise comparisons between individual donors (ASD vs ASD, Control vs Control, ASD vs Control).
    *   For each pair, run Seurat's `FindMarkers` function on the normalized RG data.
    *   Combine the results from all pairwise comparisons into a single table and save it.

4.  **Probability & Ratio Analysis (`1_pairwise_de.ipynb`):**
    *   Calculate the probability (frequency) of each gene being significantly differentially expressed (DEG) within different comparison types (e.g., within ASD donors, within Control donors, between ASD and Control donors).
    *   Employ a Leave-One-Donor-Out (LODO) approach to estimate these probabilities more robustly.
    *   Calculate ratios of these probabilities (e.g., P(DEG|AvA) / P(DEG|CvC)) to identify genes with higher variability in one group versus the other.
    *   Rank genes based on these probability ratios.

5.  **Visualization (`1_pairwise_de.ipynb`, `2_expression_plots.ipynb`):**
    *   Generate initial QC plots (cell counts, UMAPs) (`2_expression_plots.ipynb`).
    *   Generate various plots to visualize the probability and ratio analysis results (`1_pairwise_de.ipynb`):
        *   Dot plots showing DEG probabilities across comparison types.
        *   Scatter plots comparing probabilities (e.g., P(DEG|AvA) vs P(DEG|CvC)).
        *   Violin plots and histograms showing the distribution of probability ratios.
        *   Ranked plots highlighting top genes based on probability ratios.
    *   Generate expression plots (Violin plots, Dot plots) for selected top-ranked genes to visualize their expression patterns across donors and conditions (`1_pairwise_de.ipynb`).

## Scripts / Notebooks

*   `1_pairwise_de.ipynb`: Performs the core pairwise DE analysis, probability calculations, ratio analysis, and generates related visualizations.
*   `2_expression_plots.ipynb`: Performs initial data loading, subsetting, preprocessing, and generates preliminary QC plots (cell counts, UMAPs) and some expression plots.

## Related Figures

*   Supplementary Figure 22