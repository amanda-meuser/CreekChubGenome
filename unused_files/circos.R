#circos plot for synteny between Semotilus atromaculatus and Danio rerio - AP Mar 2023
#resources: 
#https://jokergoo.github.io/circlize_book/book/genomic-introduction.html
#https://www.ncbi.nlm.nih.gov/assembly/GCF_000002035.6

# LIBRARIES #####
library(circlize)
library(tidyverse)

# SYNTENY #####

# Managing synteny results from CoGe (DAGChainer results)
# Remove headers of hits (start with #, 3rd and 4th columns contain chromosome-contig names)
# grep '^#' 23058_65136.CDS-genomic.last.new.tdd10.cs0.filtered.dag.all.go_D20_g10_A5.aligncoords | cut -f 3,4 > zebra_chub_synteny.txt 

# Edit headers so all we have is the contig/chromosome number (will be label specific)
synteny <- read.table("zebra_chub_synteny.txt", col.names = c("zebra", "chub"))
synteny$zebra <- str_split_i(synteny$zebra, "_", -1)
synteny$chub <- str_split_i(synteny$chub, "g", -1)
synteny$chub <- str_split_i(synteny$chub, "l", 1)
synteny <- synteny %>% 
  mutate(chub = str_sub(chub, 4, -1)) #remove zeros from chub ID

# Remove non-whole number contig names (in zebrafish thats M and a couple scaffolds)
synteny <- synteny[!is.na(as.numeric(synteny$zebra)), ] 

# Output to csv
write.csv(synteny, "synteny_chr_numbers.csv", row.names = FALSE)

# Re-upload csv made at bottom of script (gets rid of hanging zeros in creek chub contig names)
synteny <- read.csv("synteny_chr_numbers.csv", col.names = c("zebra_count", "chub_count"))

# Remove syntenic links with fewer than 15 hits 
synteny$elim <- paste(synteny$zebra_count, synteny$chub_count) # combine count columns
synteny <- synteny[synteny$elim %in% names(which(table(synteny$elim) > 15)), ] # remove rows where values in combine column repeat <15 times
synteny <- subset(synteny, select = -c(elim)) #remove combine column

# Add the chromosome name columns (for creating links)
synteny$zebra_chr = paste0("zebra_chr", synteny[, 1])
synteny$chub_chr = paste0("chub_chr", synteny[, 2])

synteny <- synteny[order(synteny$zebra_count),]

# Make an object of creek chub contigs that have links to the zebrafish chromosomes (many don't)
unique_chub <- unique(synteny[,4])
# Same object as above but without the "chub_chr" string, used for adding colour to the circos
unique_labels <- unique(synteny[,2])

# Add colour values to synteny
for (i in 1:nrow(synteny)) {
  if (synteny$zebra_count[i] %% 2 == 0) {
    synteny$colour[i] <- "lightblue"
  }
  else {
    synteny$colour[i] <- "brown"
  }
}

# CIRCOS #####

# Read in scaffold length files
zebra_chromInfo = read.table("zebrafish_scafflength.txt", col.names = c("chr","end"))
creek_chromInfo = read.table("scafflengths_p_ctg2.txt", col.names = c("chr", "end"))

# Add "start" column populated with zeros and move columns around 
creek_chromInfo$start <- 0
creek_chromInfo <- creek_chromInfo %>% relocate(end, .after=start)
zebra_chromInfo$start <- 0
zebra_chromInfo <- zebra_chromInfo %>% relocate(end, .after=start)

# Add animal name to chr values
zebra_chromInfo[ ,1] = paste0("zebra_", zebra_chromInfo[, 1])
creek_chromInfo[ ,1] = paste0("chub_", creek_chromInfo[, 1])

# Only keep creek chub contigs that are syntenic (have to make it new dataframe)
creek_chromInfo2 <- creek_chromInfo[creek_chromInfo$chr %in% unique_chub, ]

# Maintain order of 
creek_chromInfo2 <- creek_chromInfo2 %>%
  arrange(factor(chr, rev(unique_chub)), desc(start), desc(end))

chromInfo = rbind(zebra_chromInfo, creek_chromInfo2)
circos.clear()

#pdf("circos_prototype.pdf", height=6, width=6)
circos.par(gap.after = c(rep(1, nrow(zebra_chromInfo)-1), 5, rep(1, nrow(creek_chromInfo2)-1), 5))
circos.par()
circos.genomicInitialize(chromInfo, plotType = NULL)
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(2), 
              gsub(".*chr", "", CELL_META$sector.index), cex = 0.6, niceFacing = F)
}, track.height = mm_h(1), cell.padding = c(0, 0, 0, 0), bg.border = NA)
highlight.chromosome(paste0("zebra_chr", c(1:nrow(zebra_chromInfo))), 
                     col = "brown", track.index = 1)
highlight.chromosome(paste0("chub_chr", unique_labels), 
                     col = "lightblue", track.index = 1)
circos.track(ylim = c(0, 1))
#text(-0.9, -0.9, "Zebrafish")
#text(0.9, 0.9, "Creek chub")

# links between human and mouse genomes
zebra_mid = data.frame(
  chr = synteny$zebra_chr,
  mid = round((zebra_chromInfo[synteny$zebra_count, 2] + zebra_chromInfo[synteny$zebra_count, 3])/2)
)
chub_mid = data.frame(
  chr = synteny$chub_chr,
  mid = round((creek_chromInfo[synteny$chub_count, 2] + creek_chromInfo[synteny$chub_count, 3])/2)
)
circos.genomicLink(zebra_mid, chub_mid, col = synteny$colour)

#dev.off()
circos.clear()
