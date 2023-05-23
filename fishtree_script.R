# script for creating phylogeny with Fishtree (https://fishtreeoflife.org/)
# Amanda Meuser -- May 2023

# install packages
install.packages("fishtree")
install.packages("ape")

# load packages
library(fishtree)
library(ape)

# load desired species into data object
phy <- fishtree_phylogeny(species = c("Danio rerio", "Semotilus atromaculatus", "Pimephales promelas"))
phy

# create an object to colour the species names
names <- factor(c("in", "out", "out"))
mycol <- c("firebrick2", "black")[names]

#plot phylogeny
plot(phy, use.edge.length = T, tip.color = mycol)
axisPhylo()



