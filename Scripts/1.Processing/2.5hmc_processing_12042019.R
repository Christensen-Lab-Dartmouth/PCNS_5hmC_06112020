
# Assembling & Visualizing 3D C-5mC-5hmC array of Pediatric CNS tumors
# Script author: Nasim
# Start date: 12/04/2019
# Last update: 12/06/2019`
# Notes:
# -- See "Christensen-Lab-Dartmouth/Normal-breast-5hmC" by github/owenwilkins for detail
# -- Source is David Chen's 5hmc Data Processing on github
# With 32 Linux CPUs, the entire procedure took 26:34:33 hh:mm:ss



rm(list=ls());

#Should run meta_Setup file to properly format metadata before starting this code
#I manually created the OX_BS, Label.ID and ArrayID columns in excel
#Pair column was created in meta setup 
meta <- read.csv("Metadata_complete.csv", head=TRUE, strip.white=TRUE, stringsAsFactors=FALSE)


#This code works better in R script vs markdown also this is more translateable to linux/HPC

## Set to source file directory
## If you are using command line (e.g. Linux Cluster), this line will NOT be run, so `cd <dir-of-this-R-script>` first
#if(interactive()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path));

## Load universal helper functions & packages:
#this must be an R script. It just doesn't work in markdown
source("2.PCNS_Utils.R");


## Set constants:
if(Sys.info()[1] == "Darwin") {
  print("Running locally...");
  IDAT_DIR <- "/Users/nasimazizgolshani/Dropbox/christensen/PCNS_Analysis/PCNS/Files/IDAT";
  EXPORT_PATH <- "/Users/nasimazizgolshani/Dropbox/christensen/PCNS_Analysis/PCNS/Plots/11072019_OxyBSFunnorm_PCNS.RData";
}  else if(Sys.info()[1] == "Linux") {
  print("Running on Linux Cluster");
  IDAT_DIR <- "/dartfs-hpc/rc/home/0/f002tf0/PCNS/Input/IDAT";
  EXPORT_PATH <- "/dartfs-hpc/rc/home/0/f002tf0/PCNS/Output/11072019_OxyBSFunnorm_PCNS.RData";
}


## Create R.list consisting of BS and OxBS arrays:
#setwd("/Users/nasimazizgolshani/Dropbox/christensen/PCNS_Analysis/PCNS/Files/IDAT")
listOfArrays <- prepareListOfArrays(path=IDAT_DIR, meta=meta)
arraysBS <- listOfArrays[["arraysBS"]];
arraysOxBS <- listOfArrays[["arraysOxBS"]];

## Compute signals for methylated & unmethylated channels:
#what DC does in this function is preprocessfunnorm as well as drops SNPS and technical probes
#minfi can access illumina annotation directly sidestepping requirement to load annotation files manually
#See utils R script for additional annotation
listOfSignals <- computeMethAndUnmeth(
  arraysBS = arraysBS,
  arraysOxBS = arraysOxBS, 
  dropTechProbes = FALSE,
  dropSnpProbes = FALSE
); #needs unpacking

methBS <- listOfSignals[["methBS"]];
methOxBS <- listOfSignals[["methOxBS"]];
unmethBS <- listOfSignals[["unmethBS"]];
unmethOxBS <- listOfSignals[["unmethOxBS"]];

## Aggregate signals by methylated vs. unmethylated status:
listOfAggSigs <- computeAggregatedSignals(
  methBS,
  methOxBS,
  idBS = meta$ArrayID[meta$OX_BS==0],
  idOxBS = meta$ArrayID[meta$OX_BS==1],
  rowsubset = NULL
); #needs unpacking

signalBS <- listOfAggSigs[["signalBS"]];
signalOxBS <- listOfAggSigs[["signalOxBS"]];

## Compute beta-values:
listOfBetas <- computeBetasByChemTx(
  methBS, 
  methOxBS,
  signalBS, 
  signalOxBS,
  idBS = meta$ArrayID[meta$OX_BS==0],
  idOxBS = meta$ArrayID[meta$OX_BS==1], 
  rowsubset = NULL
); #needs unpacking

betaBS <- listOfBetas[["betaBS"]];
betaOxBS <- listOfBetas[["betaOxBS"]];

## Fit maximum-likelhiood model:
MethOxy <- processSignals(
  betaBS,
  betaOxBS, 
  signalBS,
  signalOxBS
);

## Export:
save(list=c("MethOxy", "meta"), file="/dartfs-hpc/rc/home/0/f002tf0/PCNS/Output/11072019_OxyBSFunnorm_PCNS.RData", compress=TRUE);
sessionInfo();
