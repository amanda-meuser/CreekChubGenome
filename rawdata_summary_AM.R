## Script for assessing raw sequencing data
## Modified from the genome_summary_AM.R script by Amanda Meuser


## Written in Unix
## Count the number of reads
# grep "^>" *.fasta -c
# m64086e_221030_061849.hifi_reads.fastq:2040497
# m64190e_221106_234940.hifi_reads.fastq:2273297

## Count the read lengths and punt into files (our sequencing was done in two lanes and came back in two files/libraries)
# awk '/^>/ {if (seqlen){print seqlen}; printf $0"\t";seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' m64086e_221030_061849.hifi_reads.fasta > read_lengths.m64086e_221030_061849.txt
# awk '/^>/ {if (seqlen){print seqlen}; printf $0"\t";seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' m64190e_221106_234940.hifi_reads.fasta > read_lengths.m64190e_221106_234940.txt



## Written in R

scaff <- read.table("read_lengths.m64086e_221030_061849.txt", sep="\t")
head(scaff)
max(scaff$V2)
mean(scaff$V2) #15321.91
median(scaff$V2) #14456


scaff2 <- read.table("read_lengths.m64190e_221106_234940.txt", sep="\t")
head(scaff2)
max(scaff2$V2)
mean(scaff2$V2) #17378.37
median(scaff2$V2) #16804

# Plot sorted read lengths
pdf("read_lengths_CC_raw_data.pdf")
plot(sort(scaff$V2, decreasing=T), xlab="read number", ylab="read length", main="m64086e_221030_061849")
plot(sort(scaff2$V2, decreasing=T), xlab="read number", ylab="read length", main="m64190e_221106_234940")
dev.off()

# combine both libraries into one dataframe
scaff_total <- rbind(scaff, scaff2)
mean(scaff_total$V2) #16405.63
dim(scaff_total) #4313794


# Calculate coverage
genome_length <- 1099322644 # pulled from coge and/or the genome_summary_AM.R script
num_reads <- 4313794 # total number of reads, found by combining both files of reads and getting dimensions 
read_length <- 16405.63 # mean length of all reads

coverage <- (read_length*num_reads)/genome_length
coverage #64.37647