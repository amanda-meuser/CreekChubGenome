## Originally written by Liz Mandeville, modified by Amanda Meuser

## Unix
## grep "^>" *.fa -c
## creekchub_assembly_hifiasm_nov2022.bp.hap1.p_ctg.fa:568
## creekchub_assembly_hifiasm_nov2022.bp.hap2.p_ctg.fa:380
## creekchub_assembly_hifiasm_nov2022.bp.p_ctg.fa:239
## creekchub_assembly_hifiasm_nov2022.bp.p_utg.fa:13966

try to grep "n", n's mean scaffold no n's mean contig


## awk '/^>/ {if (seqlen){print seqlen}; printf $0"\t";seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' creekchub_assembly_hifiasm_nov2022.bp.p_ctg.fa > scafflengths_p_ctg.txt

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("CNEr")

library(CNEr)

N50("creekchub_assembly_hifiasm_nov2022.bp.p_ctg.fa")
N90("creekchub_assembly_hifiasm_nov2022.bp.p_ctg.fa")

scaff <- read.table("scafflengths_p_ctg.txt", sep="\t")
head(scaff)
max(scaff$V2)
mean(scaff$V2)
median(scaff$V2)


## Plot sorted scaffold length
pdf("p_ctg_sortedscaffold.pdf")
plot(sort(scaff$V2, decreasing=T), xlab="scaffold number", ylab="scaffold length")
dev.off()

## plot accumulation of genome length
pdf("p_ctg_accumulation.pdf")
plot(cumsum(sort(as.numeric(scaff$V2), decreasing=T)), xlab="scaffold number", ylab="cumulative genome length")
dev.off()

####################################################################################################################

## awk '/^>/ {if (seqlen){print seqlen}; printf $0"\t";seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' filtered.asm.cns.fa > scafflengths_filteredasmcns.txt

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
