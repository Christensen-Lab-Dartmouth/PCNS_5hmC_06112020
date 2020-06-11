#GREAT Submission
#Code by Nasim
#02/05/20
#For background on GREAT http://bejerano.stanford.edu/papers/GREAT.pdf
#https://davetang.org/muse/2014/05/16/genomic-regions-enrichment-of-annotations-tool/
#Briefly, will create a BED file with chormosomal coordinates and range of coordinates for top 5% most hydroxymethylated CpGs
#Will then submit this BED file on http://great.stanford.edu/public/html/
#Select hg19 


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
TumVsNon <- read.csv("./Script/15.EWAS/DMRsFinal/BetasTumVSNon_results020520.csv", 
                     stringsAsFactors = F, row.names = 1)

#Subset to only those that have a negative fold change (differentially hypo-hydroxymethylated)
hypohydroxyCpG <- TumVsNon[TumVsNon$logFC < 0, ]
#subset to those with adjusted adjusted p va.ue < 0.1
SigDMRs <- hypohydroxyCpG[hypohydroxyCpG$adj.P.Val <= 0.1,] 

#Generate new variables for MAPINFO. A probe lists one location, to get genomic ranges we need two locations
SigDMRs$pos_start = SigDMRs$pos
SigDMRs$pos_end <- SigDMRs$pos + 1

# -------------------------- Step 3. Make GRange Object of hypohydroxymethylated DMRs (<FDR 0.1)--------------------------

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
write.table(DMR_bed, file="./Script/15.EWAS/DMRsFinal/TumVNonDMR020520.bed", quote=F, sep="\t", row.names=F, col.names=F)


# -------------------------- Step 4. Make GRange Object of differentially hyperhydroxymethylated regions --------------------------
#Subset to only those that have a positive fold change (differentially hyper-hydroxymethylated)
hyperhydroxyCpG <- TumVsNon[TumVsNon$logFC > 0, ] #5915
#subset to those with nominal p value < 0.2 or those that have high fold change in beta value
SigDMRs <- hyperhydroxyCpG[hyperhydroxyCpG$P.Value <0.2,] 

#Generate new variables for MAPINFO. A probe lists one location, to get genomic ranges we need two locations
SigDMRs$pos_start = SigDMRs$pos
SigDMRs$pos_end <- SigDMRs$pos + 1

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
write.table(DMR_bed, file="./Script/15.EWAS/DMRsFinal/TumVNonHYPERDMR020720.bed", quote=F, sep="\t", row.names=F, col.names=F)

# -------------------------- Step 5. Make GRange Object of Top 5 CpGs--------------------------
#I ultimately did not use this
#GREAT analysis is not made to be subseted to a part of the genome!
# Load Top 5% most hyperhydroxymethylated CpGs [taken from unique individuals, no repeats]
load('./Files/01162020_top5_5hmc_unique.RData')

### Load Illumina annotation files via Bioconductor data pac kaage:
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19);
# Get Annotation
annotCpGs <- as.data.frame(getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))
#Generate new variables for MAPINFO. A probe lists one location, to get genomic ranges we need two locations
annotCpGs$pos_start = annotCpGs$pos
annotCpGs$pos_end <- annotCpGs$pos + 1

#Limit probes to top 5% CpGs
Top5annotCpGs <- annotCpGs[match(top5unique$ID,annotCpGs$Name), ];
#double check
all(Top5annotCpGs$Name %in% top5unique$ID)
identical(Top5annotCpGs$Name, top5unique$ID)

#Subset annotation to chr, pos start and pos end 
Top5Sub <- Top5annotCpGs[,c(1,4, 47, 48)]

#before proceeding check that all the fields that are in DMR subset are identically coded in top5_gr
#cannot subset great universe without this! 
Top5SubDMR <- Top5Sub[SigDMRsSub$Name,]
identical(Top5SubDMR$Name,SigDMRsSub$Name )
identical(Top5SubDMR$chr,SigDMRsSub$chr )
identical(Top5SubDMR$pos_start,SigDMRsSub$pos_start )
identical(Top5SubDMR$pos_end,SigDMRsSub$pos_end )
rm(Top5SubDMR)

#Create a 'GRanges' object from the Illumina annotation file 
top5_gr <- makeGRangesFromDataFrame(Top5Sub, keep.extra.columns=TRUE, ignore.strand=TRUE, seqnames.field = "chr", start.field="pos_start", end.field="pos_end")
#Print the GRange object to see what it looks like
top5_gr

# Convert GenomicRange objects to BED files
top5_bed <-  data.frame(cbind(chrom=as.vector(seqnames(top5_gr)),start=start(top5_gr),end=end(top5_gr), Region=top5_gr$Name))
write.table(top5_bed, file="./Script/15.EWAS/DMRsFinal/Top5GR020520.bed", quote=F, sep="\t", row.names=F, col.names=F)





