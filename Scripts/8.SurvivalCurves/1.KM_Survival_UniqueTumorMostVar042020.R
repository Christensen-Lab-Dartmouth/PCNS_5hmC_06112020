#Survival Analysis https://www.datacamp.com/community/tutorials/survival-analysis-R
#non-parametric statistic is not based on the assumption of an underlying probability 
#distribution, which makes sense since survival data has a skewed distribution
#S(t) #the survival probability at time t is given by S(t) = p.1 * p.2 * … * p.t 
#with p.1 being the proportion of all patients surviving past the first time point, 
#use the log-rank test to compare survival curves of two groups

#More on the theory:
#https://math.ucsd.edu/~rxu/math284/slect5.pdf
#The Cox Proportional Hazards model is a linear model for the log of the hazard ratio
#On interpretting Cox Prop Hazard 
#http://www.sthda.com/english/wiki/cox-proportional-hazards-model

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Step 1: Initialization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load required packages
rm(list=ls())
library(survival)
library(survminer)
library(grid)
library(dplyr)
library(tidyverse)
library(lubridate)
library(gridExtra)

setwd("/Users/nasimazizgolshani/Dropbox/christensen/PCNS_Analysis/PCNS/")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Step 2: Load Data and Reformat Covariates ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#load targets
targets <- read.csv("./Files/Metadata_complete.csv", stringsAsFactors = FALSE)

#remove targets repeat samples

targetsRPMM <- targets %>% 
  filter(Condition=="BS",
         Tumor_Type != "Non-Tumor") %>%
  select(Label.ID., DX.AGE, GENDER,location_class, Tumor_Type, DIAGNOSIS, WHO_GRADE, 
         Recurrence_Date, YRS_2_Recurrence, X1ST_SX_DATE, DATE_DEATH,
         TIME.BTWN.RELAPS..YRS.,Duration.of.Clinical.Follow.Up) %>%  
  mutate(Tumor_Type= ifelse(Tumor_Type=="EMBRYONAL", "Embryonal",
                            ifelse(Tumor_Type=="GLIOMA", "Glioma", Tumor_Type))) %>%
  mutate(Status = ifelse(DATE_DEATH == "Alive", 1, 0))

#should have 36 rows here (no controls, repeats still present)
dim(targetsRPMM)

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

#Indicate if censored or not
#data is going to be “censored” after the last time point at which 
#we know for sure that our patient did not experience death
#1 indicates censored, 0 is not
targetsRPMM <- targetsRPMM %>% mutate(CensorStat = ifelse(DATE_DEATH == "Alive", 1, 0))
#remove glioblastoma case because this sample was taken at organ harvest/no clinical data was available
targetsRPMMnogbm <- targetsRPMM[-1,]


#Recast sex as Factor 
targetsRPMMnogbm <- targetsRPMMnogbm %>%  mutate(Sex = ifelse(GENDER == "F", 1, 0),
                                                     Sex = factor(Sex, levels = c(1,0), labels = c("F","M")))

#check
class(targetsRPMMnogbm$Sex)

### Enforce data types of existing columns and clean up names for future plots
targetsRPMMnogbm$Age <- as.numeric(targetsRPMMnogbm$DX.AGE); 

#Cast RPMM cluster as factor
targetsRPMMnogbm$is.rL <- as.factor(as.integer(targetsRPMMnogbm$RPMMClusters == "rL"))

#Tumor Type as factor
targetsRPMMnogbm <- targetsRPMMnogbm %>% mutate(Tumor = ifelse(Tumor_Type == "Glioma",1,
                                                               ifelse(Tumor_Type=="Ependymoma",2,3)),
                                                Tumor = factor(Tumor, levels = c(1:3), labels = c("Glioma", "Ependymoma", "Embryonal")))

## Preview covariates
str(targetsRPMMnogbm);                                       


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Step 3: KM Survival Curve ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#create a survival object
surv_object <- Surv(time = targetsRPMMnogbm$Duration.of.Clinical.Follow.Up, event = targetsRPMMnogbm$CensorStat==0)
surv_object 

#stratify the curve by RPMM Cluster
fit1 <- survfit(surv_object ~ is.rL, data = targetsRPMMnogbm)
summary(fit1)

ggsurvplot(fit1, data = targetsRPMMnogbm,
           title = "Kaplan-Meier Curve of Survival by High and Low 5hmC Cluster",
           font.title = c(14),
           pval = TRUE, pval.method = TRUE,
           legend.title = "RPMM Cluster", 
           #be careful with labeling, this will just label the first line as rR 
           #and the second as rL, make sure rR corresponds to rR==0 by running this code 
           #without the legend labs line first then adding it
          # legend.labs = c("rR", "rL"),
          legend.labs = c("high 5hmC", "low 5hmC"),
           palette = "npg",
           xlab="Time in Years",
           risk.table = TRUE,        
         # conf.int=TRUE,
           risk.table.col = "strata",
          )


# Fit a Cox proportional hazards model
fit.coxph <- coxph(surv_object ~ Age  + Sex +Tumor + is.rL, 
                   data = targetsRPMMnogbm)

ggforest(fit.coxph, data = targetsRPMMnogbm) 

summary(fit.coxph)

#Ran sensitivity analysis with Grade, loses most significance when accounting for it 05/2020
#targetsRPMMnogbm$WHO_GRADE <- factor(targetsRPMMnogbm$WHO_GRADE)

#look at representation of tumor type in RPMM classes
TumorDistinRPMM <- targetsRPMM %>% select(RPMMClusters, Tumor_Type)
TumorDistTable <- round(prop.table(table(TumorDistinRPMM),1), 2) *100
class(TumorDistTable)
Table_Plot <- tableGrob(TumorDistTable)
grid.newpage()
grid.draw(Table_Plot)

#Look at RPMM Cluster and Time to Recurrence
#adjust years to recurrence variable to either years to recurrence or if no recurrence, 
#duration of clinical follow up
targetsRPMMnogbm <- targetsRPMMnogbm %>% 
  mutate(RecurrenceAdj = ifelse(is.na(YRS_2_Recurrence), 
                                 Duration.of.Clinical.Follow.Up, YRS_2_Recurrence)) %>%
#readjust Censor to indicate if patient had recurrence or was censored
 mutate(CensorStat = ifelse(is.na(YRS_2_Recurrence), 1, 0))


surv_object_recurr <- Surv(time = targetsRPMMnogbm$RecurrenceAdj, event = targetsRPMMnogbm$CensorStat==0)
surv_object_recurr 

fit2 <- survfit(surv_object_recurr ~ is.rL, data = targetsRPMMnogbm)


ggsurvplot(fit2, data = targetsRPMMnogbm,
           title = "Kaplan-Meier Curve of Recurrence by High and Low 5hmC Cluster",
           font.title = c(13),
           pval = TRUE, pval.method = TRUE,
           legend.title = "RPMM Cluster",    
           #Same thing, double check that rR and rL labels match!!! 
           #First run without this line
           #legend.labs = c("rR", "rL"),
           legend.labs = c("high 5hmC", "low 5hmC"),
           palette = "npg",
           xlab="Time in Years",
           ylab = "Recurrence probability",
           risk.table = TRUE,        
           risk.table.col = "strata")


# Fit a Cox proportional hazards model
fit.coxph.survival <- coxph(surv_object_recurr ~ Age + Sex +Tumor + is.rL, 
                            data = targetsRPMMnogbm)
ggforest(fit.coxph.survival, data = targetsRPMMnogbm)

summary(fit.coxph.survival)



