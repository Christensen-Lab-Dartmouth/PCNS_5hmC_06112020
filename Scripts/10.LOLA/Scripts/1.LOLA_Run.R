
#title: "PCNS - LOLA Generate Genomic Context Objects"
#author: "Nasim"
#date: "02/07/20"


#Step 1: Initialization
library(genomation)
library(GenomicRanges)
library(tidyverse)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(LOLA)
library(missMethyl)

setwd("/Users/nasimazizgolshani/Dropbox/christensen/PCNS_Analysis/PCNS/")

#Step 2: Load Data

#To run LOLA you will need to ensure that you have downloaded the LOLA region Databases: 
#http://databio.org/regiondb. Here I have downloaded LOLACoreCaches_180412.tgz which was uploaded to the site 01-May-2018 20:47. 
#You will also have to ensure that you have simpleCache package installed.
#Direct link to tar files: http://big.databio.org/regiondb/
# Load Data
TumVsNon <- read.csv("./Script/15.EWAS/DMRsFinal/BetasTumVSNon_results020520.csv", 
                     stringsAsFactors = F, row.names = 1)

# Wrangle CpGs
#Load Top 5% bed file with Genomic regions of 37k CpGs w highest 5hmC
Top5 <- read.delim("./Script/15.EWAS/DMRsFinal/Top5GR020520.bed", header=F)
colnames(Top5) <- c("chr", "pos_start", "pos_end", "Name")
#Load bed file with DHMR
DHMR <- read.delim("./Script/15.EWAS/DMRsFinal/TumVNonDMR020520.bed", header=F)
colnames(DHMR) <- c("chr", "pos_start", "pos_end", "Name")
# Define Genomic Regions
# Base comparison - aka "Universe"
GR_Universe <- makeGRangesFromDataFrame(Top5, keep.extra.columns=TRUE, ignore.strand=TRUE, seqnames.field = "chr", start.field = "pos_start", end.field = "pos_end")

# Significant CpGs
GR_Sig_CpGs <-  makeGRangesFromDataFrame(DHMR, keep.extra.columns=TRUE, ignore.strand=TRUE, seqnames.field = "chr", start.field = "pos_start", end.field = "pos_end")


# Step 3: Run LOLA
## Load Database and Query Regions

regionDB = loadRegionDB(dbLocation = "./Script/16.LOLA/Files/LOLACore/hg19")

## Run LOLA
LOLA_Sig_CpG = runLOLA(GR_Sig_CpGs, GR_Universe, regionDB)
#Notes on interpretting results
#See About page on LOLA web: http://lolaweb.databio.org/?key=APUOB8EV4FQ7WXD
# userSet and dbSet: index into their respective input region sets.
# pvalueLog: -log10(pvalue) from the fisher's exact result
# oddsRatio: result from the fisher's exact test
# support: number of regions in userSet overlapping databaseSet
# rnkPV, rnkOR, rnkSup: rank in this table of p-value, oddsRatio, and Support respectively.
# maxRnk, meanRnk: max and mean of the 3 previous ranks, providing a combined ranking system.
# b, c, d: 3 other values completing the 2x2 contingency table (with support).
# The remaining columns describe the dbSet for the row.

#Save Output
write_csv(LOLA_Sig_CpG, "./Script/16.LOLA/Files/PCNS_LOLA_Output.csv")
save(LOLA_Sig_CpG, file = "./Script/16.LOLA/Files/PCNS_LOLA_Output.RData")



