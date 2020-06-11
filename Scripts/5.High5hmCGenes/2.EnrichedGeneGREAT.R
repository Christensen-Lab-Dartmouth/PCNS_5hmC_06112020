#Pediatric CNS tumors and 5hmC
#GREAT analysis of genes enriched with high 5hmC CpGs
#Create a bed file of CpG coordinates in genes identified previously
#02/28/20
#Nasim 
#Christensen lab 
#LifeIsDope

# ==========================================================================================
# Step 1:  Initialization
# ==========================================================================================
#load relevant libraries
library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(genomation)


#Load epic Annotation dataframe with Gene Names
#For code for creating this see scripts for creating labeled volcano plots (EWAS studies)
load("./Files/epicAnnotationwGeneName.RData")
#load list of enriched genes
enrichedGenes <- read.csv("./Plots/00.052020_wPCNS0602x/10.EnrichedGenesinTop5hydroxy.csv", stringsAsFactors = F)

#Filter list of genes to those that have a proportion > 0.1
#This proportion represents the fraction of all CpGs in that gene that fall in the high 5hmC category
subEgenes <- enrichedGenes %>% filter(Proportion > 0.1)
dim(subEgenes)

#Filter epic Annotation to only those genes
epicEnrichedGenes <- epicAnnot %>% filter(Gene %in% subEgenes$Gene)
dim(epicEnrichedGenes)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Step 2: Make GRange object for all CpGs in enriched Genes# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Generate new variables for MAPINFO. A probe lists one location, to get genomic ranges we need two locations
epicEnrichedGenes$pos_start <- epicEnrichedGenes$pos
epicEnrichedGenes$pos_end <- epicEnrichedGenes$pos + 1

# BED files are tab-delimited files with one line for each genomic region.
#Required fields are chromosome, starting and ending positions. 
#can include CpG name (Name column) so that you have region ID, not required but helfpul
#Therefore we will subset our annotation file for those regions
epicGrange <- epicEnrichedGenes[,c(1,4,96, 97)]

#Create a 'GRanges' object from the Illumina annotation file 
epicEnriched_gr <- makeGRangesFromDataFrame(epicGrange, keep.extra.columns=TRUE, ignore.strand=TRUE, seqnames.field = "chr", start.field="pos_start", end.field="pos_end")
#Print the GRange object to see what it looks like
epicEnriched_gr


# Convert GenomicRange objects to BED files
epic_bed <-  data.frame(cbind(chrom=as.vector(seqnames(epicEnriched_gr)),start=start(epicEnriched_gr),end=end(epicEnriched_gr), Region=epicEnriched_gr$Name))
head(epic_bed)
write.table(epic_bed, file="./Plots/00.052020_wPCNS0602x/10.EpicEnrichedGenes.bed", quote=F, sep="\t", row.names=F, col.names=F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Step 3: Make GRange object for high 5hmC CpGs in enriched Genes# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Now subset to only genes in the top 5
load("./Files/01162020_top5_5hmc_unique.RData")
top5_enriched <- epicAnnot[epicAnnot$Name %in% top5unique$ID, ]
#subset even further to only those with "enriched genes"
Top5EnrichedGenes <- top5_enriched %>% filter(Gene %in% enrichedGenes$Gene)
dim(Top5EnrichedGenes)

#Generate new variables for MAPINFO. A probe lists one location, to get genomic ranges we need two locations
Top5EnrichedGenes$pos_start <- Top5EnrichedGenes$pos
Top5EnrichedGenes$pos_end <- Top5EnrichedGenes$pos + 1

# BED files are tab-delimited files with one line for each genomic region.
#Required fields are chromosome, starting and ending positions. 
#can include CpG name (Name column) so that you have region ID, not required but helfpul
#Therefore we will subset our annotation file for those regions
Top5Grange <- Top5EnrichedGenes[,c(1,4,96, 97)]

#Create a 'GRanges' object from the Illumina annotation file 
Top5Enriched_gr <- makeGRangesFromDataFrame(Top5Grange, keep.extra.columns=TRUE, ignore.strand=TRUE, seqnames.field = "chr", start.field="pos_start", end.field="pos_end")
#Print the GRange object to see what it looks like
Top5Enriched_gr


# Convert GenomicRange objects to BED files
Top5_bed <-  data.frame(cbind(chrom=as.vector(seqnames(Top5Enriched_gr)),start=start(Top5Enriched_gr),end=end(Top5Enriched_gr), Region=Top5Enriched_gr$Name))
head(Top5_bed)
write.table(Top5_bed, file="./Plots/00.052020_wPCNS0602x/10.Top5EnrichedGenes.bed", quote=F, sep="\t", row.names=F, col.names=F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Step 4: Make GRange object for all 5hmC CpGs in enriched Genes# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#make background bed file
#need to drop "bad probes" (sex and SNP)
#match w "Good probes file"
load("./Files/01142020_PCNS_MethOxy_GoodProbes.RData")
epicAnnot <- epicAnnot[row.names(epicAnnot) %in% row.names(PCNS_5hmc), ]
dim(epicAnnot)

#Generate new variables for MAPINFO. A probe lists one location, to get genomic ranges we need two locations
epicAnnot$pos_start <- epicAnnot$pos
epicAnnot$pos_end <- epicAnnot$pos + 1

# BED files are tab-delimited files with one line for each genomic region.
#Required fields are chromosome, starting and ending positions. 
#can include CpG name (Name column) so that you have region ID, not required but helfpul
#Therefore we will subset our annotation file for those regions
epicGrange <- epicAnnot[,c(1,4,96, 97)]

#Create a 'GRanges' object from the Illumina annotation file 
epicGrange_gr <- makeGRangesFromDataFrame(epicGrange, keep.extra.columns=TRUE, ignore.strand=TRUE, seqnames.field = "chr", start.field="pos_start", end.field="pos_end")
#Print the GRange object to see what it looks like
epicGrange_gr


# Convert GenomicRange objects to BED files
epic_bed <-  data.frame(cbind(chrom=as.vector(seqnames(epicGrange_gr)),start=start(epicGrange_gr),end=end(epicGrange_gr), Region=epicGrange_gr$Name))
head(epic_bed)

#write.table(epic_bed, file="./Script/5.GenesEnrichedwTop5CpGs/Files/epicBackground.bed", quote=F, sep="\t", row.names=F, col.names=F)
write.table(epic_bed, file="./Plots/00.052020_wPCNS0602x/10.epicBackground.bed", quote=F, sep="\t", row.names=F, col.names=F)


