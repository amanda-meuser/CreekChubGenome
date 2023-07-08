## File management for circos plots, written by Amy Pitura May 2023

## Edit the DAGchainer results for circos (should work with any CoGe DAGchainer result) #########

# exclude lines with #
# grep -v '^#' 65989_66042.genomic-CDS.last.new.tdd10.cs0.filtered.dag.all.go_D20_g10_A5.aligncoords.gcoords.txt > 65989_66042.genomic-CDS.last.new.tdd10.cs0.filtered.dag.all.go_D20_g10_A5.aligncoords.gcoords_first_edit.txt

# replace all tabs and || with a space
# sed -i -e $'s/\t/ /g;s/||/ /g' 65989_66042.genomic-CDS.last.new.tdd10.cs0.filtered.dag.all.go_D20_g10_A5.aligncoords.gcoords_first_edit.txt 

# keep columns with relevant data 
# cut -f 2,3,4,11,12,13 -d" " 65989_66042.genomic-CDS.last.new.tdd10.cs0.filtered.dag.all.go_D20_g10_A5.aligncoords.gcoords_first_edit.txt > 65989_66042.genomic-CDS.last.new.tdd10.cs0.filtered.dag.all.go_D20_g10_A5.aligncoords.gcoords_second_edit.txt

## Adding link colours #########

## edited DAGchainer results 
dfCircos <- read.table("65989_66042.genomic-CDS.last.new.tdd10.cs0.filtered.dag.all.go_D20_g10_A5.aligncoords.gcoords_second_edit.txt", header=F)

## adding colours 
for(i in 1:length(dfCircos[,1])){
  if(dfCircos$V1[i] == "000008l" || dfCircos$V1[i] == "000023l" || dfCircos$V1[i] == "000038l" || dfCircos$V1[i] == "000047l"){
    dfCircos$V7[i] <- "color=chr1"
  }else if(dfCircos$V1[i] == "000044l") {
    dfCircos$V7[i] <- "color=chr2"
  }else if(dfCircos$V1[i] == "000005l" || dfCircos$V1[i] == "000031l" || dfCircos$V1[i] == "000050l") {
    dfCircos$V7[i] <- "color=chr3"
  }else if(dfCircos$V1[i] == "000030l" || dfCircos$V1[i] == "000035l" || dfCircos$V1[i] == "000049l") {
    dfCircos$V7[i] <- "color=chr5"
  }else if(dfCircos$V1[i] == "000007l" || dfCircos$V1[i] == "000016l") {
    dfCircos$V7[i] <- "color=chr6"
  }else if(dfCircos$V1[i] == "000004l" || dfCircos$V1[i] == "000009l") {
    dfCircos$V7[i] <- "color=chr7"
  }else if(dfCircos$V1[i] == "000042l") {
    dfCircos$V7[i] <- "color=chr8"
  }else if(dfCircos$V1[i] == "000055l") {
    dfCircos$V7[i] <- "color=chr9"
  }else if(dfCircos$V1[i] == "000003l" || dfCircos$V1[i] == "000036l") {
    dfCircos$V7[i] <- "color=chr10"
  }else if(dfCircos$V1[i] == "000020l") {
    dfCircos$V7[i] <- "color=chr11"
  }else if(dfCircos$V1[i] == "000012l" || dfCircos$V1[i] == "000014l") {
    dfCircos$V7[i] <- "color=chr12"
  }else if(dfCircos$V1[i] == "000025l") {
    dfCircos$V7[i] <- "color=chr13"
  }else if(dfCircos$V1[i] == "000052l") {
    dfCircos$V7[i] <- "color=chr14"
  }else if(dfCircos$V1[i] == "000015l") {
    dfCircos$V7[i] <- "color=chr15"
  }else if(dfCircos$V1[i] == "000019l" || dfCircos$V1[i] == "000041l" || dfCircos$V1[i] == "000057l") {
    dfCircos$V7[i] <- "color=chr16"
  }else if(dfCircos$V1[i] == "000028l" || dfCircos$V1[i] == "000032l") {
    dfCircos$V7[i] <- "color=chr17"
  }else if(dfCircos$V1[i] == "000002l" || dfCircos$V1[i] == "000078l" || dfCircos$V1[i] == "000046l") {
    dfCircos$V7[i] <- "color=chr18"
  }else if(dfCircos$V1[i] == "000022l") {
    dfCircos$V7[i] <- "color=chr19"
  }else if(dfCircos$V1[i] == "000062l" || dfCircos$V1[i] == "000010l") {
    dfCircos$V7[i] <- "color=chr20"
  }else if(dfCircos$V1[i] == "000001l" || dfCircos$V1[i] == "000045l" || dfCircos$V1[i] == "000053l") {
    dfCircos$V7[i] <- "color=chr21"
  }else if(dfCircos$V1[i] == "000021l" || dfCircos$V1[i] == "000013l") {
    dfCircos$V7[i] <- "color=chr22"
  }else if(dfCircos$V1[i] == "000026l" || dfCircos$V1[i] == "000006l") {
    dfCircos$V7[i] <- "color=chr23"
  }else if(dfCircos$V1[i] == "000085l") {
    dfCircos$V7[i] <- "color=chr24"
  }else if(dfCircos$V1[i] == "000082l" || dfCircos$V1[i] == "000027l" || dfCircos$V1[i] == "000018l") {
    dfCircos$V7[i] <- "color=chr25"
  }
}

## Output
write.table(dfCircos, "65989_66042.genomic-CDS.last.new.tdd10.cs0.filtered.dag.all.go_D20_g10_A5.aligncoords.gcoords_third_edit.txt", sep = " ", quote = F, row.names = F, col.names = F) 
