# Disease Enrichment in CellOracle GRNs

**Author:** Xoel Mato Blanco

This folder contains analyses to test for the enrichment of disease-associated genes within the target gene sets of transcription factors (TFs) identified by CellOracle Gene Regulatory Networks (GRNs).

## Analysis Steps

1.  **Extract Highly Variable Genes (HVGs) (`0_get_lists_of_hvgs.ipynb`):**
    *   Load processed AnnData objects for different dataset subsets (e.g., Gliogenesis, NeuralPCW20).
    *   Extract the list of Highly Variable Genes (HVGs) identified in each subset.
    *   Save the HVG lists to text files. These lists serve as the background gene universe for enrichment tests.

2.  **Perform Enrichment Tests and Plotting (`1_enrichment_tests.ipynb`):**
    *   **Setup:** Load libraries, define parameters (datasets, cell types, diseases), load disease gene lists, define plotting palettes and helper functions.
    *   **Fisher's Exact Test:**
        *   Define functions to perform Fisher's exact test (`fisher.xy`, `fisher.xy.multiple`) comparing TF target lists from CellOracle GRNs against disease gene lists, using the corresponding HVG list as the background universe.
        *   Iterate through datasets and cell types, load GRNs, filter edges, and run enrichment tests for each TF against each disease list.
        *   Save raw enrichment results.
    *   **Process Results:**
        *   Load raw results.
        *   Perform multiple hypothesis correction (FDR) at different levels (per dataset, per dataset+disease, per dataset+disease+celltype).
        *   Calculate statistics on significant enrichments.
        *   Format results for supplementary tables and plotting.
    *   **Generate Plots:** Create various visualizations of the enrichment results, including:
        *   Line plots showing enrichment across cell types.
        *   Violin plots showing enrichment distributions.
        *   Heatmaps showing clustered enrichment patterns.
        *   Barcode plots showing top enrichments.
        *   Manhattan plots comparing top enrichments across datasets.
        *   Rank line plots detailing enrichments for specific conditions.

## Scripts / Notebooks

*   `0_get_lists_of_hvgs.ipynb`: Extracts and saves HVG lists for different dataset subsets.
*   `1_enrichment_tests.ipynb`: Performs Fisher's exact tests for disease gene enrichment in GRNs and generates plots.

## Related Figures

*    Supplementary Figure 14
