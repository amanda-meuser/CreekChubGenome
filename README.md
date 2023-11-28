# A high-quality reference genome for the common creek chub, *Semotilus atromaculatus*

R scripts posted here were used for creating plots in the aforementioned manuscript (now available on bioRxiv: https://doi.org/10.1101/2023.07.14.549000). These scripts will also be posted on DataDryad, upon publication of the manuscript. The reference genome itself can be found on NCBI, under accession number: (coming soon...).

## Description of scripts

`buscov5_script.txt` runs BUSCO v5.2.2's ray-finned fishes database on the creek chub genome. BUSCO assesses genome assembly quality by quantifying benchmarking universal single-copy orthologs.

`circos_script.conf` Circos configuration file to create a ciruclar plot that illustrates zebrafish-creek chub syntenic matches. 

`circos_script.R` edits DAGchainer results from CoGe's SynMap analysis (https://genomevolution.org/coge/SynMap.pl) to be compatible with Circos.

`fishtree_script.R` can be used to create a phylogeny for fish species, indicated in the code. It draws on a database of phylogenetic relationships between species, found on the Fish Tree of Life (https://fishtreeoflife.org/). It was originally implemented on RStudio.

`genome_summary_AM.R` and `rawdata_summary_AM.R` contain both shell and R code and can be used to compute summary statistics and create summary plots for fasta files of sequencing data, including both assembeled genome contigs and raw sequence data. These were initially implemented on command line R. 

`lets_get_kraken.txt` can be followed to create and assess genome contamination against custom Kraken2 databases. 

`repeatmodeler_script.txt` can be followed to assess repetitive regions in the creek chub genome (requires a large amount of computational resources).
