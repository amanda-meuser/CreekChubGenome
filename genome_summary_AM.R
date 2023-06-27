## Script for assessing genome assembly
## Originally written by Liz Mandeville, modified by Amanda Meuser


####################################################################################################################################
# For HiFiasm assembly
####################################################################################################################################
## Unix:
## grep "^>" *.fa -c
## creekchub_assembly_hifiasm_nov2022.bp.hap1.p_ctg.fa:568
## creekchub_assembly_hifiasm_nov2022.bp.hap2.p_ctg.fa:380
## creekchub_assembly_hifiasm_nov2022.bp.p_ctg.fa:239
## creekchub_assembly_hifiasm_nov2022.bp.p_utg.fa:13966

#try to grep "n", n's mean scaffold no n's mean contig

## grep "N" *.fa -c (just looking for lines with N)
## creekchub_assembly_hifiasm_nov2022.bp.hap1.p_ctg.fa:0
## creekchub_assembly_hifiasm_nov2022.bp.hap2.p_ctg.fa:0
## creekchub_assembly_hifiasm_nov2022.bp.p_ctg.fa:0
## creekchub_assembly_hifiasm_nov2022.bp.p_utg.fa:0

## grep -P '^(?=.*>)(?=.*N)' *.fa -c (looking for lines that have both > and N)
## creekchub_assembly_hifiasm_nov2022.bp.hap1.p_ctg.fa:0
## creekchub_assembly_hifiasm_nov2022.bp.hap2.p_ctg.fa:0
## creekchub_assembly_hifiasm_nov2022.bp.p_ctg.fa:0
## creekchub_assembly_hifiasm_nov2022.bp.p_utg.fa:0

# so that means no scaffolds, just contigs. scaffolds are made up of contigs plus gaps that are filled with N's

# create scaffold (contig?) lengths file 
# awk '/^>/ {if (seqlen){print seqlen}; printf $0"\t";seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' creekchub_assembly_hifiasm_nov2022.bp.p_ctg.fa > scafflengths_p_ctg.txt


install.packages("devtools")
devtools::install_github('A-BN/fastaUtils')
devtools::install_github("karthik/wesanderson")

library(fastaUtils)
library(wesanderson)

col <- wes_palette("Darjeeling1")

# use the fastaUtils package to calculates L50/90 and N50/90
fastanalyze(fasta = 'creekchub_assembly_hifiasm_nov2022.bp.p_ctg.fa', metrics = F, plot = F, verbose = F)

scaff <- read.table("scafflengths_p_ctg.txt", sep="\t")
head(scaff)
max(scaff$V2)
mean(scaff$V2)
median(scaff$V2)

# sort, starting at largest value and decreasing
scaff_ordered <- scaff[order(-scaff$V2),]
# length of 50th contig 
scaff_ordered[50,2]
# length of 25th contig 
scaff_ordered[25,2]


# order the cumulative lengths by largest to smallest contig
accum <- cumsum(scaff_ordered$V2)
# total genome length is total with all contigs aka final value in accum
total_len <- accum[239]
# percent in top 25 contigs:
percent50 <- (accum[50] / total_len)*100
percent50 # 95.13459%
# percent in top 50 contigs:
percent25 <- (accum[25] / total_len)*100
percent25 # 74.67159%


## Plot sorted scaffold length
pdf("p_ctg_sortedscaffold.pdf")
plot(sort(scaff$V2, decreasing=T), xlab="scaffold number", ylab="scaffold length")
dev.off()

## plot accumulation of genome length
pdf("p_ctg_accumulation.pdf")
plot(cumsum(sort(as.numeric(scaff$V2), decreasing=T)), xlab="scaffold number", ylab="cumulative genome length")
dev.off()


####################################################################################################################################
# For IPA assembly
####################################################################################################################################

## Unix: 
## cd ../creekchub_nov2022/
## grep "^>" *.fa -c
## final_purged_haplotigs.fasta:7904
## final_purged_primary.fasta:873

## Creating scaffold lengths file for both the primary and haplotig fasta files
## awk '/^>/ {if (seqlen){print seqlen}; printf $0"\t";seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' final_purged_haplotigs.fasta > scafflengths_haplotigs.txt
## awk '/^>/ {if (seqlen){print seqlen}; printf $0"\t";seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' final_purged_primary.fasta > scafflengths_primary.txt


# for primary assembly
fastanalyze(fasta = '../creekchub_IPA_nov2022/final_purged_primary.fasta', metrics = F, plot = F, verbose = F)

scaff2 <- read.table("../creekchub_IPA_nov2022/scafflengths_primary.txt", sep="\t")
head(scaff2)
max(scaff2$V2)  #23528990
mean(scaff2$V2) #1257076
median(scaff2$V2) #206534
length(scaff2$V2) #873

# sort, starting at largest value and decreasing
scaff_ordered2 <- scaff2[order(-scaff2$V2),]
# length of 50th contig 
scaff_ordered2[50,2]
# length of 25th contig 
scaff_ordered2[25,2]

# order the cumulative lengths by largest to smallest contig
accum2 <- cumsum(scaff_ordered2$V2)
# total genome length is total with all contigs aka final value in accum
total_len2 <- accum2[873]
# percent in top 25 contigs:
percent50.2 <- (accum2[50] / total_len2)*100
percent50.2 # 50.77341%
# percent in top 50 contigs:
percent25.2 <- (accum2[25] / total_len2)*100
percent25.2 # 32.48558%



## Plot sorted scaffold length
pdf("primary_sortedscaffold.pdf")
plot(sort(scaff2$V2, decreasing=T), xlab="scaffold number", ylab="scaffold length")
dev.off()

## plot accumulation of genome length
pdf("primary_accumulation.pdf")
plot(cumsum(sort(as.numeric(scaff2$V2), decreasing=T)), xlab="scaffold number", ylab="cumulative genome length")
dev.off()


# for haplotigs
fastanalyze(fasta = '../creekchub_IPA_nov2022/final_purged_haplotigs.fasta', metrics = F, plot = F, verbose = F)

scaff3 <- read.table("../creekchub_IPA_nov2022/scafflengths_haplotigs.txt", sep="\t")
head(scaff3)
max(scaff3$V2)
mean(scaff3$V2)
median(scaff3$V2)

# sort, starting at largest value and decreasing
scaff_ordered3 <- scaff3[order(-scaff3$V2),]
# length of 50th contig 
scaff_ordered3[50,2]
# length of 25th contig 
scaff_ordered3[25,2]


## Plot sorted scaffold length
pdf("haplotigs_sortedscaffold.pdf")
plot(sort(scaff3$V2, decreasing=T), xlab="scaffold number", ylab="scaffold length")
dev.off()

## plot accumulation of genome length
pdf("haplotigs_accumulation.pdf")
plot(cumsum(sort(as.numeric(scaff3$V2), decreasing=T)), xlab="scaffold number", ylab="cumulative genome length")
dev.off()


####################################################################################################################################
# for HiFiasm vs IPA plot
####################################################################################################################################

# function for adding labels to multi-panel figure, written by Gregory Garner (bitbucket.org/ggg121/r_figure_letter/src/master/)
put.fig.letter <- function(label, location="topleft", x=NULL, y=NULL, 
                           offset=c(0, 0), ...) {
  if(length(label) > 1) {
    warning("length(label) > 1, using label[1]")
  }
  if(is.null(x) | is.null(y)) {
    coords <- switch(location,
                     topleft = c(0.015,0.98),
                     topcenter = c(0.5525,0.98),
                     topright = c(0.985, 0.98),
                     bottomleft = c(0.015, 0.02), 
                     bottomcenter = c(0.5525, 0.02), 
                     bottomright = c(0.985, 0.02),
                     c(0.015, 0.98) )
  } else {
    coords <- c(x,y)
  }
  this.x <- grconvertX(coords[1] + offset[1], from="nfc", to="user")
  this.y <- grconvertY(coords[2] + offset[2], from="nfc", to="user")
  text(labels=label[1], x=this.x, y=this.y, xpd=T, ...)
}



pdf("Hifiasm_IPA_comparison_labeled.pdf")
par(mfrow=c(2,2))
plot(scaff_ordered$V2, main="HiFiasm Assembly", xlab="", ylab="Contig length", ylim=c(0,60000000), xlim=c(0,900), col = col[2])
put.fig.letter(label="(a)", location="topleft", font=2, offset=c(0.1, -0.05))
plot(scaff_ordered2$V2, main="IPA Assembly", xlab="", ylab="", ylim=c(0,60000000), xlim=c(0,900), col = col[3])
put.fig.letter(label="(b)", location="topleft", font=2, offset=c(0.1, -0.05))
plot(cumsum(scaff_ordered$V2), xlab="Contig number", ylab="Cumulative genome length", ylim=c(0,1200000000), xlim=c(0,900), col = col[2])
put.fig.letter(label="(c)", location="topleft", font=2, offset=c(0.1, -0.05))
plot(cumsum(scaff_ordered2$V2), xlab="Contig number", ylab="", ylim=c(0,1200000000), xlim=c(0,900), col = col[3])
put.fig.letter(label="(d)", location="topleft", font=2, offset=c(0.1, -0.05))
dev.off()



####################################################################################################################################
# Liz's original code
####################################################################################################################################


# awk '/^>/ {if (seqlen){print seqlen}; printf $0"\t";seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' filtered.asm.cns.fa > scafflengths_filteredasmcns.txt

scaff2 <- read.table("scafflengths_filteredasmcns.txt", sep="\t")
head(scaff2)
mean(scaff2$V2)
median(scaff2$V2)

## Plot sorted scaffold length
pdf("filteredasmcns_sortedscaffold.pdf")
plot(sort(scaff2$V2, decreasing=T), xlab="scaffold number", ylab="scaffold length")
dev.off()

## plot accumulation of genome length
pdf("filteredasmcns_accumulation.pdf")
plot(cumsum(sort(as.numeric(scaff2$V2), decreasing=T)), xlab="scaffold number", ylab="cumulative genome length")
dev.off()
