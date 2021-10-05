# BMN characterization
This repository is a summary of codes used in  [_Deep deconvolution of the hematopoietic stem cell regulatory microenvironment reveals a high degree of specialization and conservation between mouse and human_](https://www.biorxiv.org/content/10.1101/2021.07.17.452614v1).

### Preprocessing and Integration
The cells were filtered and integrated using Seurat

### Clustering


### Bootstrapping

To measure the robustness of the clusters, we customized an existing bootstrapping-based approach based on random forest [[REF1]](#ref1). See details below.

###### Summary of the bootstrapping strategy:

<img src="Figs/Bootstrapping_summary.png" alt="Bootstrapping stratgey outline" width="100%" align ="middle" />




### Reference
 <span id="ref1"> 1. Tasic, B. et al. Adult mouse cortical cell taxonomy revealed by single cell transcriptomics. Nat. Neurosci. 19, 335–346 (2016).</span>
