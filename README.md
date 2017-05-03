# Data analysis methods for the manuscript on translation efficiency

####  by Cottrell, Szczesny and Djuranovic, submitted

### Introduction 

The accession number for ribosome profiling and RNA-seq data used in this study is GSE22004. Translation efficiency (TE) was calculated using RPF and RNA-seq rpkM from mock transfection, TE = (rpkM RPF/rpkM RNA ). We obtained transcript, CDS and 3’UTR length
for human genes from Ensembl using BioMart. The tRNA adaptive index for each gene was calculated using CodonR
(https://github.com/dbgoodman/ecre_cds_analysis/tree/master/codonR). 

For this analysis the CDS of all human or Drosophila genes was obtained from the UCSC Table Browser and the tRNA gene table for human or Drosophila was obtained from the GtRNAdb. Analysis of ribosome profiling and RNA-seq data was performed in R 3.2.4 using packages listed below. Folding energies of RNAs were calculated with RNAfold, the part of Vienna package. 

### Repository data

This repository contains all the data required to run the analyses, plus the main R script. The script requires the number of libraries to be installed:

 - ggplot2 (plotting)
 - extrafont (Windows only, font for graphs)
 - Biostrings (Bioconductor)
 - biomaRt (Bioconductor)
 - seqinr (Bioconductor)
 - plyr
 - devtools
 
 Running the script will produce the results of statistical tests plus the number of graphs used in the manuscript. 