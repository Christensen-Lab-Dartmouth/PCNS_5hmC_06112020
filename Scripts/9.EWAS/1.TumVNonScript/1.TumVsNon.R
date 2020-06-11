# ==========================================================================================
# Pediatric CNS Tumors Hydroxymethylation
# EWAS - Limma  Differential Analysis  
#Tumor vs NonTumor 
# Code By: Nasim
# January 27, 2020
# ==========================================================================================

# ==========================================================================================
# Initialization
# ==========================================================================================
# Set working directory
setwd("/Users/nasimazizgolshani/Dropbox/christensen/PCNS_Analysis/PCNS")

# Load packages
rm(list=ls());
library(matrixStats);
library(limma);
library(qvalue);
library(ggplot2);
library(tidyverse)

# -------------------------- Step 2. Load data & Covariates --------------------------
#This data set is already subsetted to the top 37k CPGs with the highest mean 5hmC beta values
#It  includes only unique tumor and non-tumor samples (no repeats)
load("./Files/UniqueSamplesTop5CpGs.RData")
UniqueSamplesTop5CpGs <- as.matrix(UniqueSamplesTop5CpGs);
print( mode(UniqueSamplesTop5CpGs) ) #check to see if it's numeric

# Load covariate data
targets <- read.csv("./Files/Metadata_complete.csv", stringsAsFactors = FALSE, row.names = 1)

## Enforce data types of existing columns:
targets$GENDER <- as.factor(targets$GENDER);
targets$DX.AGE <- as.numeric(targets$DX.AGE); 

## Feature engineering, i.e. create new columns:
targets$isMale <- as.factor(as.integer(targets$GENDER == "M"));
targets$isTumor <- as.factor(as.integer(targets$Tumor_Type != "Non-Tumor"));

## Preview covariates
str(targets); 

#drop aberrant non-tumor sample
targets <- targets[-40,]

#add ID to rownames of targets
row.names(targets) <- targets$Label.ID.

# -------------------------- Step 3. Data Preparation for differential Analysis --------------------------
## Important: Before "set row names to match Beta Values column name"
## First, mutually subset data and M values:
targets <- targets[rownames(targets) %in% colnames(UniqueSamplesTop5CpGs), ];
Betas <- UniqueSamplesTop5CpGs[ , colnames(UniqueSamplesTop5CpGs) %in% rownames(targets)];

## Second (VERY IMPORTANT), check if the column names of `M.values`` **matches** with the row names of targets
## If so, print out a message to say "You're ready to move ahead"
##If not, execute a **match**
if(identical(row.names(targets),colnames(Betas))) {
  print("Column names of M values & row names of covariates match! OK to proceed");
} else {
  print("The column names & row names didn't match!"); #returns an error message 
  print("Proceed to matching...");
  targets <- targets[match(colnames(M.Values), row.names(targets)), ];
}


# --------------------------Step 5. LIMMA: Tumor vs. Nontumor --------------------------
## **Triple-check**: If rownames of covariates & column names of M-values don't match, the `stopifnot` function will throw an error message
## If matching is a success, nothing will happen
stopifnot( identical(row.names(targets), colnames(Betas)) );


## Preview your covariates before building the linear model
table(targets$isMale);
table(targets$isTumor);

hist(targets$DX.AGE);
qplot(as.numeric(Betas), geom="histogram", xlab="5hmC Beta values", ylab="Count")


##Create Design Matrix:
design_istumor <- model.matrix( ~ isTumor + isMale + DX.AGE, data=targets);
View(design_istumor); #check

#Perform differential analysis:
fit_istumor <- lmFit(Betas, design=design_istumor);
fit_istumor <- eBayes(fit_istumor);


### Compile a list of CpGs with Illumina annotation
### Load Illumina annotation files via Bioconductor data pac kaage:
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19);
# Get Annotation
annotCpGs <- as.data.frame(getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))
#Subset to CpGs in top 5%
ann850KSub <- annotCpGs[match(rownames(Betas),annotCpGs$Name), ];
#Double check that CpGs are the same!
all(ann850KSub$Name %in% row.names(Betas))
identical(ann850KSub$Name,row.names(Betas))

#Retrieve linear model output for the main variable (tumor vs. non-tumor):
DMP_istumor <- topTable(
  fit_istumor,
  number = Inf,
  coef = "isTumor1", #if unsure, go back to Design Matrix
  genelist = ann850KSub, 
  adjust.method = "fdr",
  sort.by = "p"
);
DMP_istumor$negLog10FDR <- -log10(DMP_istumor$adj.P.Val);

#save a copy of your results to file (edit path as necessary):
write.csv(DMP_istumor, file="./Script/15.EWAS/DMRsFinal/BetasTumVSNon_results020520.csv", quote=FALSE);



