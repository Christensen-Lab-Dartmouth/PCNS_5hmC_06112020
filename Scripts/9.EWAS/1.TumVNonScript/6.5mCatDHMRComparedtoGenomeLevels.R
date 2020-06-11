#PCNS tumor and 5hmC
#Examine average 5mC at DHMRs as compared to 5mC in the rest of the genome
#Nasim
#02/07/20

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Step 1: Initialization~~~~~~~~~~~~~~~~~~~~~~~~~~~~~````
setwd("/Users/nasimazizgolshani/Dropbox/christensen/PCNS_Analysis/PCNS")

rm(list = ls())
library(ggplot2)
library(tidyverse)
library(wesanderson)
library(ggpubr)
library(ggsci)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Step 2: Load Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~````

# methylation data 
load("./Files/01142020_PCNS_MethOxy_GoodProbes.RData")

#Load covariates
targets <- read.csv("./Files/Metadata_complete.csv", stringsAsFactors = FALSE)

#Since metadata (targets) has duplicates for BS and OxBS samples, want to subset for just ArrayID used in my betas file
rownames(targets) <- targets$ArrayID
targets <- targets[colnames(PCNS_5mc),]

#Change column names in dataframes to patient ID from ArrayID
colnames(PCNS_5mc) <- targets$Label.ID.[match(names(PCNS_5mc),targets$ArrayID)]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Step 3: Find Mean 5mC at DHMRs  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~````
#Load the DHMR list
TumVsNon <- read.csv("./Script/15.EWAS/DMRsFinal/BetasTumVSNon_results020520.csv", 
                     stringsAsFactors = F, row.names = 1)
#subset median values to DHMRs (hypo-hydroxymethylated and adj p value < 0.1)
SigDMRs <- TumVsNon[TumVsNon$logFC < 0, ]
SigDMRs <- SigDMRs[SigDMRs$adj.P.Val <0.1, ]

#Subset 5mC data to SigDMRs
#I am only interested in differentially hypo-hydroxymethylated CpGs (due to clinical relevance)
DHMR_5mc <- PCNS_5mc[row.names(PCNS_5mc) %in% row.names(SigDMRs),]
dim(DHMR_5mc)

# Determine the average level of 5mC per CpG across all subjects 
#### mC
DHMRmethly_sum <- apply(DHMR_5mc, 2, sum, na.rm=TRUE)
DHMRavg_5mC <- (DHMRmethly_sum/nrow(DHMR_5mc))
DHMRavg_5mC <- as.data.frame(DHMRavg_5mC)
DHMRavg_5mC <- rownames_to_column(DHMRavg_5mC, var = "ID")

#Change column names
colnames(DHMRavg_5mC) <- c("ID", "avg_5mC")
#Add Tumor Type column to this
DHMRavg_5mC$TumorType <- targets$Tumor_Type[match(targets$Label.ID.,DHMRavg_5mC$ID)]
#Add Column indicating context
DHMRavg_5mC$Context <- "DHMR"

#subset to unique samples
serials <- c("PCNS0302X","PCNS0602X", "PCNS0603X")

multiregional <- c("PCNS1001A", "PCNS0301B", "PCNS0301C", "PCNS0601B", "PCNS0901B", "PCNS1901C")

DHMRavg_5mC <- DHMRavg_5mC[!DHMRavg_5mC$ID %in% c(multiregional, serials), ]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Step 4: Find Genomic Mean 5mC  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~````
# First Exclude DHMRs from list of Cpgs
PCNS_5mc <- PCNS_5mc[!row.names(PCNS_5mc) %in% row.names(DHMR_5mc), ]
#sanity check
#sum of DHMR rows and now modified 5mC data frame should equal the nrows in 5hmC dataframe
identical(nrow(PCNS_5mc) + nrow(DHMR_5mc), nrow(PCNS_5hmc))

#Determine the average level of 5mC per CpG across all subjects 
#### mC
methly_sum <- apply(PCNS_5mc, 2, sum, na.rm=TRUE)
avg_5mC <- (methly_sum/nrow(PCNS_5mc))
avg_5mC <- as.data.frame(avg_5mC)
avg_5mC <- rownames_to_column(avg_5mC, var = "ID")

#Add Tumor Type column to this
avg_5mC$TumorType <- targets$Tumor_Type[match(targets$Label.ID.,avg_5mC$ID)]
#Add Column indicating context
avg_5mC$Context <- "Genomic"

#subset to unique samples
avg_5mC <- avg_5mC[!avg_5mC$ID %in% c(multiregional, serials), ]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Step 4: Pairwise Comparison between 5mC means  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~````
#Combine the dataframes
Combo5mCMeans <- rbind(DHMRavg_5mC,avg_5mC)

#Drop aberrant control 
Combo5mCMeans <- filter(Combo5mCMeans, ID != "PCNS3311X")

#Clean up Tumor Type
Combo5mCMeans <- Combo5mCMeans %>% mutate(TumorType = ifelse(TumorType == "GLIOMA", "Glioma",
                                                             ifelse(TumorType=="EMBRYONAL", "Embryonal", TumorType)))

#Pairwise comparison between 5mC means at DHMRs and genomic 5mC levels
compare_means(avg_5mC ~ Context, data = Combo5mCMeans, 
              group.by = "TumorType", paired = TRUE)

# Box plot facetted by "TumorType"
p <- ggpaired(Combo5mCMeans, x = "Context", y = "avg_5mC",
              color = "TumorType", palette = "npg", 
              line.color = "gray", line.size = 0.4,
              facet.by = "TumorType", short.panel.labs = FALSE,
             # title="Genomic and DHMR Methylation",
              ylab="Mean 5mC Beta Value",
              xlab= "Tumor Type"
              )

# Use only p.format as label. Remove method name.
p + stat_compare_means(label = "p.format", paired = TRUE)


