# CellOracle analyses
_by Xoel Mato Blanco_

## Noack et al. 2021
Related to Supp. Figure 15  

### Prepare data + run Cicero
- `1_DataCodePrep.ipynb`

### CellOracle pipeline
- `2_CellOracle.ipynb`


## Trevino et al. 2021
Reanalysis considering **ASD genes** under `Trevino_dataset_ASD_genes` folder (related to Fig. 5.)    

### Prepare subsets
- `0_PrepareTrevinoData.ipynb`

### CellOracle pipeline
- `1_Gliogenesis.ipynb`
- `1_NeuralPCW20.ipynb`
- `1_NeuralPCW21.ipynb`
- `1_NeuralPCW24.ipynb`

### Process results
- `2_ResultCompilation.ipynb`
- `5_NeMO_export.ipynb`

### Plot results
Related to Figure 4, Supp. Figure 10, 11, 12, 13 and 14  
- `0_Rplots.R`
- `3_GenePresence.ipynb`
- `3_PerturbationHeatmap.v5.ipynb`
- `3_ScoreComparisons`

### Comparison with mouse results
Related to Supp. Figure 15  
- `4_NoackComparisons`
