# ASD iPSC analyses

## Bulk RNA-seq
_by Gabriel Santpere Baro & Xoel Mato Blanco_

### Prepare count matrix
- `1_CountsBulkRNA.R`

### Differential Expression Analysis
- `2_BulkLinesDEX.ipynb`

## Single-cell RNA-seq

### Annotation
_by Shaojie Ma_
- See README.md in folder

### Analyses
_by Xoel Mato Blanco_

#### Data preparation
- `0.reintegrateNoHighMito.ipynb`
- `0.Common.ipynb`
- `1.DataPreprocess.ipynb`

#### Differential expression analyses
- `2.FindMarkers.ipynb`
- `3.Pseudobulk.ipynb`
- `4.FilterMarkers.ipynb`

#### Proportion tests
- `5.ProportionTests.ipynb`

#### Plots
- `6.Plots.Vulcanos.ipynb`
- `7.Plots.Dotplot.ipynb`
- `8.Plots.VennDiagrams.Donors.ipynb`
- `9.Plots.HeatmapDE.ipynb`
- `10.Plots.densities_scanpy.ipynb`

#### Comparison with Jourdan et al. 2023
_by Alex Jourdon_
- `11.AJ_MatoBlanco_etal_code.R`