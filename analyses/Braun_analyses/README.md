# Analysis of Braun et al. 2023 Fetal Brain Data

**Author:** Xoel Mato Blanco

This folder contains analyses performed on the human fetal brain single-cell RNA-seq dataset from Braun et al. 2023 (likely referring to the Linnarsson lab data often associated with Braun). The analyses include data downloading, re-annotation, disease gene expression analysis, comparison with other datasets, and cell type enrichment analysis using EWCE.

## Analysis Steps

1.  **Data Download and Preparation (`0_download_data.ipynb`):**
    *   Downloads the primary HDF5 dataset (`HumanFetalBrainPool.h5`) and potentially an H5AD file.
    *   Extracts expression data (sparse matrix) and metadata (cell info, gene info, embeddings, factors) from the HDF5 file.
    *   Saves extracted metadata tables to CSV files.
    *   Loads updated cell annotations from an external Excel file (`cluster_lineage_annotation.summary.xlsx`, generated in step 2) and merges them with the cell metadata.
    *   Subsets the data based on brain subregion and age, saving expression matrices (MTX format) and metadata (CSV) for each subset into separate folders.

2.  **Neuronal Lineage Re-annotation (`1_neuronal_lineage_annotation.ipynb`):**
    *   Loads supplementary table S2 from the original publication and external annotation info.
    *   Classifies cell clusters as Excitatory or Inhibitory based on gene markers (e.g., VGLUT, GAD2) and existing annotations.
    *   Checks for consistency between gene-based and annotation-based classifications.
    *   Saves an intermediate table (`cluster_lineage_annotation.xlsx`).
    *   Loads a manually reviewed/commented version of the intermediate table (`cluster_lineage_annotation.commented.xlsx`).
    *   Resolves inconsistencies based on manual comments and predefined rules.
    *   Assigns final cell type classifications (`LongName`, `ShortName`, `InEx` status).
    *   Saves the final resolved annotation summary (`cluster_lineage_annotation.summary.xlsx`), which is used as input in step 1.

3.  **Disease Gene Expression Analysis (`2_disease_expression.ipynb`):**
    *   Reloads the HDF5 data.
    *   Loads disease information (`DiseaseInfo.csv`) and gene-disease associations (`parsed_lists_with_nicola.union.csv`).
    *   Subsets the data to select only Radial Glia (RG) cells from specific regions.
    *   Calculates expression statistics (total counts, percentage of expressing cells, average expression) for genes within the RG subset.
    *   For each disease gene list, calculates the number and proportion of associated genes that are expressed in the RG subset, using different expression thresholds (>=1 count, >=5 counts, >=5% cells).
    *   Saves these disease expression summary statistics to separate CSV files.

4.  **Comparison with In Vitro Data (`3_disease_expression_comparison.ipynb`):**
    *   Loads graphical objects (`graphical.rda`), disease gene lists, *in vitro* expression data (log2 RPKM), and the *in vivo* disease expression summaries generated in step 3.
    *   Calculates summary statistics for disease gene expression in the *in vitro* data.
    *   Combines *in vivo* (Braun/Linnarsson RG) and *in vitro* summary statistics.
    *   Calculates additional metrics (e.g., ratio of expressed genes to present genes).
    *   Generates and saves comparative plots (bar plots, scatter plots) visualizing disease gene expression across the two datasets.

5.  **EWCE Preparation (`4_prepare_ewce.ipynb`):**
    *   Loads gene metadata.
    *   Defines helper functions for data normalization (`sctransform`) and Cell Type Data (CTD) generation (`EWCE::generate_celltype_data`).
    *   Iterates through the subregion/age subfolders created in step 1.
    *   For each subfolder, loads the expression matrix and cell metadata.
    *   Applies SCTransform normalization to the expression data.
    *   Generates a CTD object using the normalized expression and cell type annotations (`ShortName`).
    *   Saves the generated CTD object (`ctd_*.rda`) back into the corresponding subfolder.

6.  **EWCE Analysis Execution (`5_run_ewce.ipynb`):**
    *   Lists all CTD files generated in step 5.
    *   Loads the disease gene lists used in step 3.
    *   Sets the number of bootstrap repetitions (e.g., 10000).
    *   Iterates through each CTD file.
    *   For each CTD, runs `EWCE::bootstrap_enrichment_test` for every disease gene list against the annotation level(s) present in the CTD.
    *   Aggregates the enrichment results.
    *   Saves the aggregated EWCE results (e.g., `results.csv`, `hit.cells.csv`) into the same subfolder as the input CTD file.

7.  **EWCE Result Visualization (`6_plot_ewce.ipynb`):**
    *   Loads all EWCE result files (`results.csv`) generated in step 6.
    *   Loads disease metadata (`DiseaseInfo.csv`) and cell type/region/age color palettes.
    *   Merges EWCE results with disease metadata.
    *   Calculates corrected q-values and significance flags.
    *   Defines functions for creating Manhattan plots and heatmaps.
    *   Generates various plots visualizing EWCE enrichment results across diseases, cell types, subregions, and ages (Manhattan plots, heatmaps).
    *   Saves selected plots to PDF files.
    *   Loads cell metadata (`CellInfo.annotated.csv`) and generates bar plots showing cell type proportions across subregions and ages.

## Scripts / Notebooks

*   `0_download_data.ipynb`: Data download, extraction, and subsetting.
*   `1_neuronal_lineage_annotation.ipynb`: Cell type re-annotation and consistency checks.
*   `2_disease_expression.ipynb`: Disease gene expression calculation in Radial Glia.
*   `3_disease_expression_comparison.ipynb`: Comparison of disease gene expression between *in vivo* and *in vitro* data.
*   `4_prepare_ewce.ipynb`: Preparation of Cell Type Data (CTD) files for EWCE.
*   `5_run_ewce.ipynb`: Execution of EWCE enrichment tests.
*   `6_plot_ewce.ipynb`: Visualization of EWCE results and cell proportions.

## Related Figures

*   Figure 1
*   Supplementary Figure 1
