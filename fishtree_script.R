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


# save as a pdf
pdf("fishtree_phylogeny.pdf", width = 5, height = 4)
plot(phy, use.edge.length = T, tip.color = mycol)
dev.off()

#-----------------------------------------------------------
# add more fish

# load desired species into data object
phy_more <- fishtree_phylogeny(species = c("Danio rerio", "Semotilus atromaculatus", "Pimephales promelas", "Carassius auratus", "Oncorhynchus mykiss", "Gasterosteus aculeatus", "Ctenopharyngodon idella"))
phy_more

# create an object to colour the species names
names_more <- factor(c("out", "out", "in", "in", "out", "out", "in"))
mycol_more <- c("firebrick2", "black")[names_more]

#plot phylogeny
plot(phy_more, use.edge.length = T, tip.color = mycol_more)


# save as a pdf
pdf("fishtree_phylogeny_morefish.pdf", width = 5, height = 4)
plot(phy_more, use.edge.length = T, tip.color = mycol_more)
dev.off()

