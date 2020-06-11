#Establish subset of probes for analysis
#These will exclude sex probes and SNP associated probes
#Source of this code is MM Field Cancerization/DC RefFree Preprocess
#Author: Nasim Azizgolshani

library(tidyverse)
library(stringr)
library(ggpubr)

setwd("/Users/nasimazizgolshani/Dropbox/christensen/PCNS_Analysis/PCNS")
# Load annotation from Zhou et al 2016 for SNP annotation
ZhouGRCh38 <- read.delim('./Script/7.RefFree/ScriptBSonly/Files/EPIC.hg38.manifest.tsv')
#Load Processed Data and Beta Values
load("/Users/nasimazizgolshani/Dropbox/christensen/PCNS_Analysis/PCNS/Files/12.04.Tandem_OxBS_Betas/12062019_OxyBSFunnorm_PCNS.RData")
PCNS_5hmc <- MethOxy[, ,3] #starting #CpG sites is 865859
PCNS_5hmc <- as.data.frame(PCNS_5hmc)

library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19);
annot.850kb3 <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19);
sexProbes <- as.character(annot.850kb3$Name[annot.850kb3$chr %in% c("chrX", "chrY")]); #probes to REMOVE to get autosomal only

PCNS_5hmc <- PCNS_5hmc[!row.names(PCNS_5hmc) %in% sexProbes, ]; 
dim(PCNS_5hmc) #846232

# Filter SNP associated and cross-hybridizing probes (Zhou et. al. 2016 annotation)
cpgRemove <- ZhouGRCh38 %>% filter(MASK_general == TRUE) #105454 probes
cpgRemove <- cpgRemove[cpgRemove$probeID %in% row.names(PCNS_5hmc), ] #58199
PCNS_5hmc <- PCNS_5hmc[!row.names(PCNS_5hmc) %in%  cpgRemove$probeID,] #Probes left 743461

# Indicate number of probes removed
paste(formatC(nrow(cpgRemove),big.mark = ','),'SNP associated or cross-hybridizing  CpGs removed')

#Now Repeat for 5mc
PCNS_5mc <- MethOxy[, , 2]
PCNS_5mc <- as.data.frame(PCNS_5mc)

PCNS_5mc <- PCNS_5mc[!row.names(PCNS_5mc) %in% sexProbes, ]; 
dim(PCNS_5mc)

# Filter SNP associated and cross-hybridizing probes (Zhou et. al. 2016 annotation)
cpgRemove <- ZhouGRCh38 %>% filter(MASK_general == TRUE) #64144 probes
cpgRemove <- cpgRemove[cpgRemove$probeID %in%  row.names(PCNS_5mc), ] #43171
PCNS_5mc <- PCNS_5mc[!row.names(PCNS_5mc)%in%  cpgRemove$probeID, ]
# Indicate number of probes removed
paste(formatC(nrow(cpgRemove),big.mark = ','),'SNP associated or cross-hybridizing  CpGs removed')

#For 5C
PCNS_5c <- MethOxy[, , 1]
PCNS_5c <- as.data.frame(PCNS_5c)

PCNS_5c <- PCNS_5c[!row.names(PCNS_5c) %in% sexProbes, ]; 
dim(PCNS_5c)

# Filter SNP associated and cross-hybridizing probes (Zhou et. al. 2016 annotation)
cpgRemove <- ZhouGRCh38 %>% filter(MASK_general == TRUE) #64144 probes
cpgRemove <- cpgRemove[cpgRemove$probeID %in%  row.names(PCNS_5c), ] #43171
PCNS_5c <- PCNS_5c[!row.names(PCNS_5c)%in%  cpgRemove$probeID, ]
# Indicate number of probes removed
paste(formatC(nrow(cpgRemove),big.mark = ','),'SNP associated or cross-hybridizing  CpGs removed')

nrow(PCNS_5hmc) #743461
nrow(PCNS_5mc) #743461
nrow(PCNS_5c) #743461
save(list=c("PCNS_5hmc", "PCNS_5mc", "PCNS_5c"), file="/Users/nasimazizgolshani/Dropbox/christensen/PCNS_Analysis/PCNS/Files/01142020_PCNS_MethOxy_GoodProbes.RData", compress=TRUE)
