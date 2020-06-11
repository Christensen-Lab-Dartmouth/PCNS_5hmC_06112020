#GREAT Submission
#Embryonal v NonTumor
#Code by Nasim
#02/10/20

# ==========================================================================================
# Step 1: Initialization
# ==========================================================================================
# Set working directory
setwd("/Users/nasimazizgolshani/Dropbox/christensen/PCNS_Analysis/PCNS")

#clear workspace
rm(list=ls())
# Load Required Packages
library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(genomation)

# -------------------------- Step 2. Load data  --------------------------
EmbVsNon <- read.csv("./Plots/0.052020_wPCNS0601X/SubtypeEWAS/BetasEmbVSNon_results0520.csv", 
                     stringsAsFactors = F, row.names = 1)

#Subset to only those that have a negative fold change (differentially hypo-hydroxymethylated)
hypohydroxyCpG <- EmbVsNon[EmbVsNon$logFC < 0, ]
#subset to those with adjusted adjusted p value
SigDMRs <- hypohydroxyCpG[hypohydroxyCpG$adj.P.Val <0.02,] 

#Generate new variables for MAPINFO. A probe lists one location, to get genomic ranges we need two locations
SigDMRs$pos_start = SigDMRs$pos
SigDMRs$pos_end <- SigDMRs$pos + 1

# -------------------------- Step 3. Make GRange Object of DMRs--------------------------

# BED files are tab-delimited files with one line for each genomic region.
#Required fields are chromosome, starting and ending positions. 
#can include CpG name (Name column) so that you have region ID, not required but helfpul
#Therefore we will subset our annotation file for those regions
SigDMRsSub <- SigDMRs[,c(1,4,54, 55)]

#Create a 'GRanges' object from the Illumina annotation file 
DMR_gr <- makeGRangesFromDataFrame(SigDMRsSub, keep.extra.columns=TRUE, ignore.strand=TRUE, seqnames.field = "chr", start.field="pos_start", end.field="pos_end")
#Print the GRange object to see what it looks like
DMR_gr


# Convert GenomicRange objects to BED files
DMR_bed <-  data.frame(cbind(chrom=as.vector(seqnames(DMR_gr)),start=start(DMR_gr),end=end(DMR_gr), Region=DMR_gr$Name))
write.table(DMR_bed, file="./Plots/0.052020_wPCNS0601X/SubtypeEWAS/EmbVSNonGREATsub0520.bed", quote=F, sep="\t", row.names=F, col.names=F)



