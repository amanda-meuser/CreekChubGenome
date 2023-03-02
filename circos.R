#circos plot for synteny between Semotilus atromaculatus and Danio rerio
#https://jokergoo.github.io/circlize_book/book/genomic-introduction.html
#https://www.ncbi.nlm.nih.gov/assembly/GCF_000002035.6

# grep '^#' 23058_65136.CDS-genomic.last.new.tdd10.cs0.filtered.dag.all.go_D20_g10_A5.aligncoords | cut -f 3,4 > zebra_chub_synteny.txt 

#install.packages("tidyverse")
library(circlize)
library(tidyverse)

human_chromInfo = read.chromInfo(species = "hg19")$df
mouse_chromInfo = read.chromInfo(species = "mm10")$df
creek_chromInfo = read.chromInfo("scafflengths_p_ctg.txt")
creek_chromInfo2 = read.table("scafflengths_p_ctg2.txt", col.names = c("chr", "end"))
creek_chromInfo2$start <- 0
creek_chromInfo <- creek_chromInfo2 %>% relocate(end, .after=start)

human_chromInfo[ ,1] = paste0("human_", human_chromInfo[, 1])
mouse_chromInfo[ ,1] = paste0("mouse_", mouse_chromInfo[, 1])
creek_chromInfo[ ,1] = paste0("chub_", creek_chromInfo[, 1])
chromInfo = rbind(human_chromInfo, creek_chromInfo)
chromInfo2 = rbind(human_chromInfo, mouse_chromInfo)
head(chromInfo)

circos.par(gap.degree=1, gap.after = c(rep(1, 23), 5, rep(1, 127), 5))
circos.par()
circos.genomicInitialize(chromInfo, plotType = NULL)
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(2), 
              gsub(".*chr", "", CELL_META$sector.index), cex = 0.6, niceFacing = F)
}, track.height = mm_h(1), cell.padding = c(0, 0, 0, 0), bg.border = NA)
highlight.chromosome(paste0("human_chr", c(1:22, "X", "Y")), 
                     col = "lightblue", track.index = 1)
highlight.chromosome(paste0("chub_chr", c(1:127)), 
                     col = "lightgreen", track.index = 1)
circos.track(ylim = c(0, 1))


# links between human and mouse genomes
human_mid = data.frame(
  chr = paste0("human_chr", 1:19),
  mid = round((human_chromInfo[1:19, 2] + human_chromInfo[1:19, 3])/2)
)
numbers <- c(1,23,34,45,36,29,110,123,79,45,68,37,127,100,50,78,17,18,19)

mouse_mid = data.frame(
  chr = paste0("chub_chr", numbers),
  mid = round((creek_chromInfo[numbers, 2] + creek_chromInfo[numbers, 3])/2)
)
circos.genomicLink(human_mid, mouse_mid, col = rand_color(19))
circos.clear()
text(-0.9, -0.8, "Human\ngenome")
text(0.9, 0.8, "Mouse\ngenome")

synteny <- read.table("zebra_chub_synteny.txt")
synteny$V1 <- str_split_i(synteny$V1, "_", -1)
synteny$V2 <- str_split_i(synteny$V2, "g", -1)
synteny$V2 <- str_split_i(synteny$V2, "l", 1)
synteny <- synteny %>% 
  mutate(V2 = str_sub(V2, 4, -1))

for (i in nrow(synteny$V2)) {
  if (startsWith(synteny$V2[i], "00")) {
    synteny$V2[i] = str_split_i(synteny$V2[i], "00", -1)
  }
}
startsWith(synteny$V2[1], "00")
