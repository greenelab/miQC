# miQC
Flexible, probablistic metrics for quality control of scRNA-seq data

## Description
Single-cell RNA-sequencing (scRNA-seq) has made it possible to profile gene expression in tissues at high resolution.
An important preprocessing step prior to performing downstream analyses is to identify and remove cells with poor or degraded sample quality using quality control (QC) metrics.
Two widely used QC metrics to identify a ‘low-quality’ cell are (i) if the cell includes a high proportion of reads that map to mitochondrial DNA encoded genes (mtDNA) and (ii) if a small number of genes are detected.
miQC is data-driven QC metric that jointly models both the proportion of reads mapping to mtDNA and the number of detected genes with mixture models in a probabilistic framework to predict which cells are low-quality in a given dataset.

## Installation
This will hopefully be on bioconductor soon. Until then, the best way to install is:
`remotes::install_github("greenelab/miQC", build_vignettes = TRUE)`
