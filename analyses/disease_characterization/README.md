# Disease characterization
_by Gabriel Santpere Baro_ 

## CoGAPS patterns  
Related to Figure 1, Supp. Fig. 1, 2 and 5  
- `1_Run_CoGAPS.R` _by Carlo Colantuoni_
- `2_PrepareCoGaps.R` 

## Run MAGMA
```{bash}
# Example for Cogaps_top5.txt
for i in $(ls GWAS_MAGMA/RESULTS/*raw)
do
name=$(basename $i)
magma --gene-results ${i} --set-annot GENESETS/Cogaps_top5.txt col=2,1 --out GENESETS/ENRICH_COGAPS/${name}
done
```

## Visualize MAGMA and CoGAPS by disease
Related to Figure 1, Supp. Fig. 1, 2 and 3  
- `3_PrepareSFARI.R`
- `4_PlotEnrichmentCombined.R`

## Early to late neuronal differentiation and maturation trends
Related to Figure 1, 2 and Supp. Fig. 3, 5 and 7  
- `5_Slopes_V3.R` 

## Visualization of distribution of disease genes in CHD8/FMRP targets
Related to Supp. Figure 4  
- `6_CHD8_FMRP_targets_V3.R`
- `7_Summaries.R`

## Heatmaps of disease genes
Related to Figure 1 and Supp. Fig. 3  
- `8_Intersection_V13.R` 

## Distribution of dorsoventral bias among disease genes and patterning center marker genes
Related to Figure 2 and Supp. Fig. 3, 5 and 7  
- `9_PatterningCenters.R`
- `10_CompareFC_DESEQ_RPKM.R`
- `11_DV_bias_in_monkey_PCs.ipynb` _by Xoel Mato Blanco_
