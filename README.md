# A high-quality reference genome for the common creek chub, *Semotilus atromaculatus*

R scripts posted here were used for creating plots in the aforementioned manuscript (now available on bioRxiv: https://doi.org/10.1101/2023.07.14.549000). The reference genome itself can be found on NCBI, under GenBank accession number GCA_031834385.1, and raw data from Pac-Bio long-read whole genome sequencing can be found on the Sequence Read Archive, under accession number SRX21895170. 

## Description of scripts

`buscov5_script.txt` runs BUSCO v5.2.2's ray-finned fishes database on the creek chub genome. BUSCO assesses genome assembly quality by quantifying benchmarking universal single-copy orthologs.

`circos_script.conf` Circos configuration file to create a ciruclar plot that illustrates zebrafish-creek chub syntenic matches. 

`circos_script.R` edits DAGchainer results from CoGe's SynMap analysis (https://genomevolution.org/coge/SynMap.pl) to be compatible with Circos.

`circos_tips.txt` tips for successfully running Circos.

`fishtree_script.R` can be used to create a phylogeny for fish species, indicated in the code. It draws on a database of phylogenetic relationships between species, found on the Fish Tree of Life (https://fishtreeoflife.org/). It was originally implemented on RStudio.

`genome_summary_AM.R` and `rawdata_summary_AM.R` contain both shell and R code and can be used to compute summary statistics and create summary plots for fasta files of sequencing data, including both assembeled genome contigs and raw sequence data. These were initially implemented on command line R. 

`lets_get_kraken.txt` can be followed to create and assess genome contamination against custom Kraken2 databases. 

`repeatmodeler_script.txt` can be followed to assess repetitive regions in the creek chub genome (requires a large amount of computational resources).

## DNA Barcoding files

`AMO22-0800-redo-Fwd_f12.ab1` and `AMO22-0800-redo-Fwd_f12.seq` both contain the forward sequence that was created to check the species identification of the reference genome specimen via DNA barcoding. 

`BOLD_finaltree.pdf` shows the placement of the barcode sequence among other sequences in the BOLD data base, while `IDEngine_Results_Summary.xls` shows the call of species with percent certainty. 

## Circos files 

`65989_66042.genomic-CDS.last.new.tdd10.cs0.filtered.dag.all.go_D20_g10_A5.aligncoords.gcoords_edited_R.txt` CoGe DAGChainer syntenic matches between creek chub and fathead minnow. 

`65989_66058.genomic-CDS.last.new.tdd10.cs0.filtered.dag.all.go_D20_g10_A5.aligncoords_edited_R.txt` CoGe DAGChainer syntenic matches between creek chub and zebrafish. 

`karyotype.creekchub.top50.txt`, `karyotype.fathead.txt`, and `karyotype.zebrafish.txt` are properly formatted karyotype files for Circos.

