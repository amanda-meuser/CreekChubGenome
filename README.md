# A high-quality reference genome for the common creek chub, *Semotilus atromaculatus*

R scripts posted here were used for creating plots in the aforementioned manuscript (now available on bioRxiv: https://doi.org/10.1101/2023.07.14.549000). These scripts will also be posted on DataDryad, upon publication of the manuscript. The reference genome itself can be found on NCBI, under accession number: (coming soon...).

## Description of scripts

`circos_script.R` uses DAGchainer results from CoGe's SynMap analysis (https://genomevolution.org/coge/SynMap.pl) to create a Circos plot that illustrates syntenic matches between two genomes. 

`fishtree_script.R` can be used to create a phylogeny for fish species, indicated in the code. It draws on a database of phylogenetic relationships between species, found on the Fish Tree of Life (https://fishtreeoflife.org/). It was originally implemented on RStudio.

`genome_summary_AM.R` and `rawdata_summary_AM.R` contain both shell and R code and can be used to compute summary statistics and create summary plots for fasta files of sequencing data, including both assembeled genome contigs and raw sequence data. These were initially implemented on command line R. 