#circos plot for synteny between Semotilus atromaculatus and Danio rerio
#https://jokergoo.github.io/circlize_book/book/genomic-introduction.html
#https://www.ncbi.nlm.nih.gov/assembly/GCF_000002035.6

# grep '^#' 23058_65136.CDS-genomic.last.new.tdd10.cs0.filtered.dag.all.go_D20_g10_A5.aligncoords | cut -f 3,4 > zebra_chub_synteny.txt 

#install.packages("tidyverse")
library(circlize)
library(tidyverse)

synteny2 <- read.csv("synteny4.csv", col.names = c("zebra_count", "chub_count","colour"))

synteny2$zebra_chr = paste0("zebra_chr", synteny2[, 1])
synteny2$chub_chr = paste0("chub_chr", synteny2[, 2])

zebra_mid2 <- synteny2[,1]
chub_mid2 <- synteny2[,2]
zebra_chr <- synteny2[,4]
chub_chr <- synteny2[,5]

unique_chub <- unique(synteny2[,5])
unique_blue <- unique(synteny2[,2])

zebra_chromInfo = read.table("zebrafish_scafflength.txt", col.names = c("chr","end"))
creek_chromInfo = read.table("scafflengths_p_ctg2.txt", col.names = c("chr", "end"))
creek_chromInfo$start <- 0
creek_chromInfo <- creek_chromInfo %>% relocate(end, .after=start)
zebra_chromInfo$start <- 0
zebra_chromInfo <- zebra_chromInfo %>% relocate(end, .after=start)

zebra_chromInfo[ ,1] = paste0("zebra_", zebra_chromInfo[, 1])
creek_chromInfo[ ,1] = paste0("chub_", creek_chromInfo[, 1])
creek_chromInfo2 <- subset(creek_chromInfo, chr %in% unique_chub)

creek_chromInfo2 <- creek_chromInfo2 %>%
  arrange(factor(chr, rev(unique_chub)), desc(start), desc(end))
#write_csv(creek_chromInfo, "creekchromInfo.csv")
#creek_chromInfo <- read.csv("creekchromInfo.csv", col.names = c("chr", "start", "end"))

chromInfo = rbind(zebra_chromInfo, creek_chromInfo2)
head(chromInfo)
circos.clear()

pdf("circos_prototype.pdf", height=6, width=6)
circos.par(gap.after = c(rep(1, 24), 5, rep(1, 35), 5))
circos.par()
circos.genomicInitialize(chromInfo, plotType = NULL)
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(2), 
              gsub(".*chr", "", CELL_META$sector.index), cex = 0.6, niceFacing = F)
}, track.height = mm_h(1), cell.padding = c(0, 0, 0, 0), bg.border = NA)
highlight.chromosome(paste0("zebra_chr", c(1:25)), 
                     col = "brown", track.index = 1)
highlight.chromosome(paste0("chub_chr", unique_blue), 
                     col = "lightblue", track.index = 1)
circos.track(ylim = c(0, 1))
text(-0.9, -0.9, "Zebrafish")
text(0.9, 0.9, "Creek chub")

# links between human and mouse genomes
zebra_mid = data.frame(
  chr = zebra_chr,
  mid = round((zebra_chromInfo[zebra_mid2, 2] + zebra_chromInfo[zebra_mid2, 3])/2)
)
chub_mid = data.frame(
  chr = chub_chr,
  mid = round((creek_chromInfo[chub_mid2, 2] + creek_chromInfo[chub_mid2, 3])/2)
)
circos.genomicLink(zebra_mid, chub_mid, col = synteny2$colour)

dev.off()
circos.clear()

# Managing synteny results from CoGe (DAGChainer results)
# Remove headers of hits (start with #, 3rd and 4th columns contain chromosome-contig connection)
# grep '^#' 23058_65136.CDS-genomic.last.new.tdd10.cs0.filtered.dag.all.go_D20_g10_A5.aligncoords | cut -f 3,4 > zebra_chub_synteny.txt 

# Edit headers 
synteny <- read.table("zebra_chub_synteny.txt", col.names = c("zebra", "chub"))
synteny$zebra <- str_split_i(synteny$zebra, "_", -1)
synteny$chub <- str_split_i(synteny$chub, "g", -1)
synteny$chub <- str_split_i(synteny$chub, "l", 1)
synteny <- synteny %>% 
  mutate(chub = str_sub(chub, 4, -1))

synteny <- synteny[!is.na(as.numeric(synteny$zebra)), ] 
write.csv(synteny, "synteny2.csv", row.names = FALSE)

