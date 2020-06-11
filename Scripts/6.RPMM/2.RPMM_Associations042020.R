#Pediatric CNS study
#Compare 5hmC means in two RPMM clusters
#Or total 5hmC (tbd)
#Nasim
#03/01/20

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Step 1 Initialization#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("/Users/nasimazizgolshani/Dropbox/christensen/PCNS_Analysis/PCNS")

rm(list = ls())
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(epitools) #OR calculation
library(grid)
library(gridExtra)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Step 2: Load and reformat Data #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# methylation data 
load("./Files/01162020_unique_tumor_betas.RData")

#Since metadata (targets) has duplicates for BS and OxBS samples, want to subset for just ArrayID used in my betas file
#load targets
targets <- read.csv("./Files/Metadata_complete.csv", stringsAsFactors = FALSE, row.names = 1)

#remove targets repeat samples
targetsRPMM <- targets %>% 
  filter(Condition=="BS",
         Tumor_Type != "Non-Tumor") %>%
  select(Label.ID., DX.AGE, GENDER,location_class, Tumor_Type, Diagnosis_Clean, WHO_GRADE, 
         Recurrence_Date, YRS_2_Recurrence, X1ST_SX_DATE, DATE_DEATH,
         TIME.BTWN.RELAPS..YRS.,Duration.of.Clinical.Follow.Up) %>%  
  mutate(Tumor_Type= ifelse(Tumor_Type=="EMBRYONAL", "Embryonal",
                            ifelse(Tumor_Type=="GLIOMA", "Glioma", Tumor_Type))) %>%
  mutate(Status = ifelse(DATE_DEATH == "Alive", 1, 0))

#load RPMM Clusters and add them to annotation
RPMM <- read.csv("./Script/9.RPMM/RPMM_TumorMosVar_Script/TumorsRPMM.csv", stringsAsFactors = FALSE)

#remove repeat samples from targets
#by matching samples to RPMM samples, will remove repeats as I ran RPMM on unique samples only!
targetsRPMM <- targetsRPMM[targetsRPMM$Label.ID. %in% RPMM$Sample_Name,]


#Add RPMM Clusters 
RPMM <- RPMM %>% mutate(Label.ID. = Sample_Name) %>%
  select(Label.ID., RPMMClusters, RPMMSampleOrder)
targetsRPMM <- merge(targetsRPMM, RPMM, by= "Label.ID.")
#27 unique tumors = 27 rows
dim(targetsRPMM)
rm(RPMM)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Step 3: Determine Average/Total 5hmC #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# visualize total 5hmC/5mC content across samples 
# Determine the average level of 5hmC/5mC per CpG across all subjects 
#### hmC
hydroxy_sum <- apply(tumor_5hmc_unique, 2, sum, na.rm=TRUE)

#Combine the dataframes into one with a column annotating tumor type
hydroxy_avg <- as.data.frame(hydroxy_sum) %>% 
  rownames_to_column(var ="Label.ID.") %>%
  mutate(Average5hmC= hydroxy_sum/nrow(tumor_5hmc_unique))

AvgWCovar <- merge(hydroxy_avg, targetsRPMM, by="Label.ID.")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Step 4: Determine Associations  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Test 1: Is RPMM Cluster Associated with Mean 5hmC?

#first turn RPMM Clusters into factors, then group data by levels
AvgWCovar$RPMMClusters <- as.factor(AvgWCovar$RPMMClusters) 
AvgWCovar <- AvgWCovar %>% group_by(RPMMClusters)

#Test Mean 5hmC and cluster membership
kruskal.test(Average5hmC  ~ RPMMClusters, data = AvgWCovar)
#Kruskal-Wallis chi-squared = 19.074, df = 1, p-value = 1.258e-05

#Test Cluster membership and tumor type
AvgWCovar$Tumor_Type <- as.factor(AvgWCovar$Tumor_Type)
AvgWCovar <- AvgWCovar %>% group_by(Tumor_Type) %>%
  mutate(is.rR = ifelse(RPMMClusters == "rR", 1, 0))


#Test 2: Is RPMM Cluster Associated with Tumor Type?
## Create OR table for different tumor types
OR_Table <- as.data.frame(matrix(ncol=5, nrow=0))
colnames(OR_Table) <- c("proportion", "OR", "lower.CI", "upper.CI", "pVal")

TumorLevels <- c("Embryonal", "Ependymoma", "Glioma")
TumorLabels <- c("Embryonal", "Ependymoma", "Glioma")

#warning will arise due to low sample sizes!
for (i in TumorLevels) {
  # Calculate proportion of sites in that context
  prop <- table(AvgWCovar$is.rR,
               AvgWCovar$Tumor_Type == i)['1','TRUE']/sum(AvgWCovar$is.rR)
  # Calculate odds ratio
  tempOR <- oddsratio.fisher(table(AvgWCovar$is.rR,AvgWCovar$Tumor_Type == i))
  # Add values to output
  OR_Table[i,] <- c(prop,tempOR$measure[2,],tempOR$p.value[2,'fisher.exact'])
}

OR_Table <- cbind(c('Tumor Type', TumorLabels),
                       c('Odds Ratio (95% CI)', paste(format(round(OR_Table$OR,2), nsmall= 2),' (',format(round(OR_Table$lower.CI,2), nsmall=2), '-',
                                                      format(round(OR_Table$upper.CI,2), nsmall=2),')', sep= '')),
                       c('P value', formatC(OR_Table$pVal, digits=2, format= 'E')))
#Results
# [,1]         [,2]                  [,3]      
# [1,] "Tumor Type" "Odds Ratio (95% CI)" "P value" 
# [2,] "Embryonal"  "0.50 (0.04- 4.05)"   "6.62E-01"
# [3,] "Ependymoma" "2.39 (0.31-21.36)"   "3.91E-01"
# [4,] "Glioma"     "0.84 (0.14- 4.98)"   "1.00E+00"
#Test 3: Is RPMM Cluster Associated with Tumor Location?
AvgWCovar$isSubtentorial <- as.factor(ifelse(AvgWCovar$location_class == "subtentorial",1,0))
tempOR <- oddsratio.fisher(table(AvgWCovar$isSubtentorial, AvgWCovar$is.rR))
#Result = not significant

#Test 4: Is RPMM Cluster Associated with Age?
kruskal.test(DX.AGE  ~ RPMMClusters, data = AvgWCovar)

#Kruskal-Wallis chi-squared = 0.95267, df = 1, p-value = 0.329


#Check to see if it correlates with grade
AvgWCovar <- AvgWCovar %>% mutate(isHighGrade = ifelse(WHO_GRADE == 3, 1, 
                                                       ifelse(WHO_GRADE ==4, 1,0))) 
AvgWCovar$isHighGrade <- factor(AvgWCovar$isHighGrade, levels = c(1,0), labels = c(1,0))

#Test 5: Is RPMM Cluster Associated with Grade?
kruskal.test(isHighGrade  ~ RPMMClusters, data = AvgWCovar)
#Kruskal-Wallis chi-squared = 2.7251, df = 1, p-value = 0.09878
