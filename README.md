# Data analysis methods for the manuscript on translation efficiency

## by Cottrell, Szczesny and Djuranovic, submitted

### Introduction

This is repository accompanying the manuscript entitled "Translation efficiency is a determinant of the magnitude of miRNA-mediated repression" by Kyle A Cottrell, Pawel Szczesny, Sergej Djuranovic.

### Manuscript methods

The accession number for ribosome profiling and RNA-seq data used in this study is GSE22004. Fold change of RPF and RNA-seq was calculated as described in Guo et al., 2010. Translation efficiency (TE) was calculated using RPF and RNA-seq rpkM from mock transfection, TE = (rpkMRPF/rpkMRNA). We obtained transcript, CDS and 3’UTR length for human genes from Ensembl using BioMart. mRNA half-lives were obtained from Tani et al., 2012. miR-155 or miR-1 targets were predicted using TargetScan. The tRNA adapative index for each gene was calculated using CodonR (https://github.com/dbgoodman/ecre_cds_analysis/tree/master/codonR). For this analysis the CDS of all human or Drosophila genes was obtained from the UCSC Table Browser and the tRNA gene table for human or Drosophila was obtained from the GtRNAdb.

### Repository data

This repository contains all the data required to run the analyses, plus the main R script. The script requires the number of libraries to be installed:

- ggplot2 (plotting)
- extrafont (Windows only, font for graphs)
- Biostrings (Bioconductor)
- biomaRt (Bioconductor)
- seqinr (Bioconductor)
- plyr
- devtools

 Running the script will produce the results of statistical tests plus the number of graphs used in the manuscript. While the plots are made using Windows-only font, the latter parts are explicitly using paths and tools available under linux operating system (such as Vienna package).