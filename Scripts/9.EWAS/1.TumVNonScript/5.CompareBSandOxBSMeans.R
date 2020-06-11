#PCNS tumor and 5hmC
#Examine average 5mC and 5hmC at DHMRs
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
#Load covariates
targets <- read.csv("./Files/Metadata_complete.csv", stringsAsFactors = FALSE)

#Since metadata (targets) has duplicates for BS and OxBS samples, want to subset for just ArrayID used in my betas file
rownames(targets) <- targets$ArrayID

# methylation data 
load("./Files/01142020_PCNS_MethOxy_GoodProbes.RData")
targets <- targets[colnames(PCNS_5mc),]

#Change column names in dataframes to patient ID from ArrayID
colnames(PCNS_5mc) <- targets$Label.ID.[match(names(PCNS_5mc),targets$ArrayID)]

#load data with only unique samples
load("./Files/021520UniqueSamples.RData")

#subset to only unique samples
Unique5mc <- PCNS_5mc[, colnames(PCNS_5mc) %in% colnames(All5hmc)]

#remove other datasets to avoid confusion
rm(All5hmc, PCNS_5c, PCNS_5hmc, PCNS_5mc)

#Load the DHMR list
TumVsNon <- read.csv("./Script/15.EWAS/DMRsFinal/BetasTumVSNon_results020520.csv", 
                     stringsAsFactors = F, row.names = 1)

#subset median values to DHMRs (hypo-hydroxymethylated and adj p value < 0.1)
SigDMRs <- TumVsNon[TumVsNon$logFC < 0, ]
SigDMRs <- SigDMRs[SigDMRs$adj.P.Val <0.1, ]

#Subset 5mC data to SigDMRs
#I am only interested in differentially hypo-hydroxymethylated CpGs (due to clinical relevance)
Unique5mc <- Unique5mc[row.names(Unique5mc) %in% row.names(SigDMRs),]
dim(Unique5mc)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Step 3: Calculate Averages for 5mC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~````

# visualize total 5hmC/5mC content across samples 
# Determine the average level of 5mC per CpG across all subjects 
#### mC
methly_sum <- apply(Unique5mc, 2, sum, na.rm=TRUE)
avg_5mC <- (methly_sum/nrow(Unique5mc))
avg_5mC <- as.data.frame(avg_5mC)
avg_5mC <- rownames_to_column(avg_5mC, var = "ID")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Step 4:Subset by Tumor Type ~~~~~~~~~~~~~~~~~~~~~~~~~~~~````
#Start subsetting for tumor types and controls
GliomaSamples <- targets %>% filter(Tumor_Type=="Glioma")
EpendymomaSamples <- targets %>% filter(Tumor_Type=="Ependymoma")
EmbryonalSamples <- targets %>% filter(Tumor_Type=="Embryonal")
NonTumorSamples <- targets %>% filter(Tumor_Type=="Non-Tumor")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Step 5: Plotting Means~~~~~~~~~~~~~~~~~~~~~~~~~~~~````

#Box Plots w P Values
OxBSmethly_avg <- as.data.frame(methly_sum) %>% 
  rownames_to_column(var ="Sample") %>%
  mutate(Average5mC= methly_sum/nrow(Unique5mc)) %>%
  mutate(Type = ifelse(Sample %in% GliomaSamples$Label.ID., "Glioma", 
                       ifelse(Sample %in% EpendymomaSamples$Label.ID., "Ependymoma",
                              ifelse(Sample %in% EmbryonalSamples$Label.ID., "Embryonal",
                                     ifelse(Sample %in% NonTumorSamples$Label.ID., "Non-Tumor",NA)
                              )
                       )
  )) 

OxBSmethly_avg$Type <- as.factor(OxBSmethly_avg$Type)
compare_means(Average5mC ~ Type, data = OxBSmethly_avg, method = "t.test", paired = FALSE, ref.group = "Non-Tumor")
my_comparisons <- list( c("Non-Tumor", "Glioma"), c("Non-Tumor", "Embryonal"), c("Non-Tumor", "Ependymoma") )


#TumorCols = c(`Non-Tumor` = "seagreen", Embryonal= "palevioletred1", Ependymoma="darkgoldenrod1", Glioma= "steelblue4") 
p <- ggplot(OxBSmethly_avg, aes(x = Type, y = Average5mC, colour=Type, 
                            group=paste(Type))) +
  geom_boxplot()+
  geom_point(position=position_jitterdodge(jitter.width = 0.3, dodge.width = 1.0), alpha=0.5) +
  ggtitle("Average 5-mC Beta Value of Samples at \n Differentially Hypo-Hydroxymethylated Regions \nFrom OxBS Derived Beta Values") +
  theme_classic()+
  theme(
    text = element_text(size=16),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.text = element_text(size = 18))+
  labs(x = "Tumor Type",
       y = "Mean 5mC Beta Value") +
  scale_color_npg()+
  scale_y_continuous(limits=c(0.3, 1))

OxBSplot <- p + stat_compare_means(comparisons = my_comparisons)

print(OxBSplot)



##########################################################################################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Now Repeat for BS treated files which gives both 5hmC and 5mC ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##########################################################################################################################################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Step 2: Load Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~````

# methylation data 
load("./Script/7.RefFree/ScriptBSonly/Files/BSbetas011420.RData")

#Load the DHMR list
TumVsNon <- read.csv("./Script/15.EWAS/DMRsFinal/BetasTumVSNon_results020520.csv", 
                     stringsAsFactors = F, row.names = 1)
#subset median values to DHMRs (hypo-hydroxymethylated and adj p value < 0.1)
SigDMRs <- TumVsNon[TumVsNon$logFC < 0, ]
SigDMRs <- SigDMRs[SigDMRs$adj.P.Val <0.1, ]

#Subset 5mC data to SigDMRs
#I am only interested in differentially hypo-hydroxymethylated CpGs (due to clinical relevance)
BSbetas_pcns <- BSbetas_pcns[row.names(BSbetas_pcns) %in% row.names(SigDMRs),]
dim(BSbetas_pcns)

#Load covariates
targets <- read.csv("./Files/Metadata_complete.csv", stringsAsFactors = FALSE)

#Since metadata (targets) has duplicates for BS and OxBS samples, want to subset for just ArrayID used in my betas file
rownames(targets) <- targets$ArrayID
targets <- targets[row.names(targets) %in% colnames(BSbetas_pcns), ]

#Change column names in dataframes to patient ID from ArrayID
colnames(BSbetas_pcns) <- targets$Label.ID.[match(colnames(BSbetas_pcns),targets$ArrayID)]
BSbetas_pcns <- BSbetas_pcns[ ,targets$Label.ID.]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Step 3: Calculate Averages for 5mC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~````

# visualize total 5hmC/5mC content across samples 
# Determine the average level of 5mC per CpG across all subjects 
#### mC
methly_sum <- apply(BSbetas_pcns, 2, sum, na.rm=TRUE)
avg_5mC <- (methly_sum/nrow(BSbetas_pcns))
avg_5mC <- as.data.frame(avg_5mC)
avg_5mC <- rownames_to_column(avg_5mC, var = "ID")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Step 4: Plotting Means~~~~~~~~~~~~~~~~~~~~~~~~~~~~````

#Box Plots w P Values
BSmethly_avg <- as.data.frame(methly_sum) %>% 
  rownames_to_column(var ="Sample") %>%
  mutate(Average5mC= methly_sum/nrow(BSbetas_pcns)) %>%
  mutate(Type = ifelse(Sample %in% GliomaSamples$Label.ID., "Glioma", 
                       ifelse(Sample %in% EpendymomaSamples$Label.ID., "Ependymoma",
                              ifelse(Sample %in% EmbryonalSamples$Label.ID., "Embryonal",
                                     ifelse(Sample %in% NonTumorSamples$Label.ID., "Non-Tumor",NA)
                              )
                       )
  )) %>%
  filter(Sample != "PCNS3311X")

BSmethly_avg$Type <- as.factor(BSmethly_avg$Type)
compare_means(Average5mC ~ Type, data = BSmethly_avg, method = "t.test", paired = FALSE, ref.group = "Non-Tumor")
my_comparisons <- list( c("Non-Tumor", "Glioma"), c("Non-Tumor", "Embryonal"), c("Non-Tumor", "Ependymoma") )



p <- ggplot(BSmethly_avg, aes(x = Type, y = Average5mC, colour=Type, 
                            group=paste(Type))) +
  geom_boxplot()+
  geom_point(position=position_jitterdodge(jitter.width = 0.3, dodge.width = 1.0), alpha=0.5) +
  ggtitle("Average 5-mC Beta Value of Samples at \nDifferentially Hypo-Hydroxymethylated Regions \nFrom BS Derived Beta Values") +
  theme_classic()+
  theme(
    text = element_text(size=16),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 17, hjust = 0.5),
    legend.text = element_text(size = 16))+
  labs(x = "Tumor Type",
       y = "Mean 5mC Beta Value") +
  scale_color_npg()+
  scale_y_continuous(limits=c(0.3, 1))

BSplot <- p + stat_compare_means(comparisons = my_comparisons)
print(BSplot)



##########################################################################################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Combine Plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##########################################################################################################################################################################################
library(gridExtra)
grid.arrange(OxBSplot, BSplot, ncol=2)
