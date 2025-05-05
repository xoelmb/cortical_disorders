This repository contains the code used for the paper **_Early Developmental Origins of Cortical Disorders Modeled in Human Neural Stem Cells_**.

Each folder under `analyses` contains the code used for the different analyses performed in the paper. The README.md file in each folder provides a brief description of the analyses performed in that folder.

The main analysis categories include:

* **disease_characterization:** Characterization of disease gene expression and enrichment in _in vitro_ NSC differentiation time course.
* **rcistarget:** Identification of potential regulators using RcisTarget based on co-expression modules.
* **annotate_progenitors:** Annotation of progenitor cell types.
* **celloracle:** Gene Regulatory Network inference using CellOracle.
* **asd_analyses:** Specific analyses related to Autism Spectrum Disorder datasets.
* **Braun_analyses:** Analysis of the Braun et al. 2023 fetal brain dataset, including re-annotation, disease gene expression, and EWCE enrichment.
* **Complexity:** Analysis of transcriptome and GRN complexity.
* **Disease_enrichment_in_CO_GRN:** Enrichment testing of disease genes within CellOracle GRNs.
* **Jourdon_pairwise_DE_comparison:** Pairwise differential expression analysis on the Jourdon et al. 2023 dataset to assess expression variability.

### Citation

**Early Developmental Origins of Cortical Disorders Modeled in Human Neural Stem Cells**  
Xoel Mato-Blanco, Suel-Kee Kim, Alexandre Jourdon, Shaojie Ma, Andrew T.N. Tebbenkamp, Fuchen Liu, Alvaro Duque, Flora M. Vaccarino, Nenad Sestan, Carlo Colantuoni, Pasko Rakic, Gabriel Santpere, Nicola Micali  
bioRxiv 2024.06.14.598925; doi: <https://doi.org/10.1101/2024.06.14.598925>

```bibtex
@article {Mato-Blanco2024.06.14.598925,
    author = {Mato-Blanco, Xoel and Kim, Suel-Kee and Jourdon, Alexandre and Ma, Shaojie and Tebbenkamp, Andrew T.N. and Liu, Fuchen and Duque, Alvaro and Vaccarino, Flora M. and Sestan, Nenad and Colantuoni, Carlo and Rakic, Pasko and Santpere, Gabriel and Micali, Nicola},
    title = {Early Developmental Origins of Cortical Disorders Modeled in Human Neural Stem Cells},
    elocation-id = {2024.06.14.598925},
    year = {2024},
    doi = {10.1101/2024.06.14.598925},
    publisher = {Cold Spring Harbor Laboratory},
    abstract = {The implications of the early phases of human telencephalic development, involving neural stem cells (NSCs), in the etiology of cortical disorders remain elusive. Here, we explored the expression dynamics of cortical and neuropsychiatric disorder-associated genes in datasets generated from human NSCs across telencephalic fate transitions in vitro and in vivo. We identified risk genes expressed in brain organizers and sequential gene regulatory networks across corticogenesis revealing disease-specific critical phases, when NSCs are more vulnerable to gene dysfunctions, and converging signaling across multiple diseases. Moreover, we simulated the impact of risk transcription factor (TF) depletions on different neural cell types spanning the developing human neocortex and observed a spatiotemporal-dependent effect for each perturbation. Finally, single-cell transcriptomics of newly generated autism-affected patient-derived NSCs in vitro revealed recurrent alterations of TFs orchestrating brain patterning and NSC lineage commitment. This work opens new perspectives to explore human brain dysfunctions at the early phases of development.One-sentence summary The temporal analysis of gene regulatory networks in human neural stem cells reveals multiple early critical phases associated with cortical disorders and neuropsychiatric traits.},
    URL = {https://www.biorxiv.org/content/early/2024/06/14/2024.06.14.598925},
    eprint = {https://www.biorxiv.org/content/early/2024/06/14/2024.06.14.598925.full.pdf},
    journal = {bioRxiv}
}
```
