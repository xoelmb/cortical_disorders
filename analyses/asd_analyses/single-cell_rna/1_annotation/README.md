---

title: QC of the ASD iPSC data
author: Shaojie Ma
date: Mar 06, 2024

---


### Source functions
sc.fun.R
preprocess.fun.R


### Preliminary processing
```bash
## run doublet prediction
Rscript db.scrublet.R

## Read the count and convert to Seurat objects
## add cell cycle scores & downsample data to avoid quality bias between control and ASD
## removed doublets (only a few cells)
Rscript counts.hvg.extract.R

## Integrate the control and ASD samples (with cell cycle correction)
Rscript inte.cycle-correct.seurat.v2.R
```


### Integrate with macaque data
```bash
## prepare macaque data for integration
Rscript mac.inte.predata.v2.R

## integrate with macaque NSC & neurons
Rscript mac.inte.full.R

## integrate with only macaque NSC
Rscript nsc.inte.R

## For the above integrations, visualize the results in a better way
Rscript vis.integration.R

## For the integration with macaque NSCs, further subseting the cells for better visualization
Rscript nsc.inte.sub.R

## Label transfer of NSC idnetites
## version-1: just cortical NSCs (no patterning centers). This version was finally adopted. 
Rscript nsc.label.transfer.R
## version-2: cortical NSCs+patterning center+mesenchymal
Rscript nsc-pat-mes.label.transfer.R


## Since cluster 10 contains some cells resembling FGF17/8 organizer cells, I further subclustered this cluster to specifically label these cells. 
Rscript pat.sub.R
```

### Summarize the results
```bash
## summarize the results annotate cells
## plot the sankey plot (seurat-clusters--versus--predicted macaque NSC identities)
## plot custom markers (verifying the cluster identities)
Rscript anno.combine.R
```



