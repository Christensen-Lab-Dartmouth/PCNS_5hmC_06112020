# Recursively Partitioned Mixture Model (RPMM)
# Date: 01/10/2020
# Notes:

rm(list=ls());
library(matrixStats);
library(RPMM);
library(seriation)
library(doParallel); registerDoParallel(detectCores() - 1);
setwd("/Users/nasimazizgolshani/Dropbox/christensen/PCNS_Analysis/PCNS/")

#set up targets info
targets <- read.csv("./Files/Metadata_complete.csv", stringsAsFactors = FALSE)
load("./Files/01162020_unique_tumor_betas.RData");

OUT_FILE_NAME <-"/Users/nasimazizgolshani/Dropbox/christensen/PCNS_Analysis/PCNS/Plots/00.052020_wPCNS0602x/EWAS/TumorsRPMM.csv"

#Convert to matrix
PCNSmatrix <- as.matrix(tumor_5hmc_unique)

#Select and save unique tumors top 10k most variable CpGs
sele <- order(rowVars(PCNSmatrix), decreasing=TRUE)[1:10000]; 
TumorTopVar <- tumor_5hmc_unique[sele, ]
#save original dataframe
save(TumorTopVar, file="/Users/nasimazizgolshani/Dropbox/christensen/PCNS_Analysis/PCNS/Script/9.RPMM/RPMM_TumorMosVar_Script/Tumor_MostVar042020.RData")

#Select and save unique tumors top 10k most variable CpGs
sele <- order(rowVars(PCNSmatrix), decreasing=TRUE)[1:10000]; 
mostVar <- PCNSmatrix[sele, ]
class(mostVar)


#Create function to select most variable CpGs to submit to RPMM
load_Yinv <- function(n.CpG=10000) {
  #'@description Load, subset, and transform data for RPMM
  sele <- order(rowVars(PCNSmatrix), decreasing=TRUE)[1:n.CpG]; 
  Y_inv <- t(PCNSmatrix[sele, ]); 
  assign("Y_inv", Y_inv, envir=.GlobalEnv);
}

getRPMMClustLabels <- function(rpmmObject, Y_inv=NULL) {
  #'@description Extracts RPMM hard cluster labels
  #'@param rpmmObject RPMM object
  #'@param Y_inv Optional. Input matrix for RPMM computation. If provided, the sample name will be updated
  
  hardLabels <- blcTreeLeafClasses(rpmmObject);
  hardLabels <- as.data.frame(hardLabels);
  colnames(hardLabels) <- "RPMMClusters";
  if(! is.null(Y_inv)) rownames(hardLabels) <- hardLabels$Sample_Name <- rownames(Y_inv);
  return(hardLabels);
}

getRPMMSampOrder <- function(rpmmClusters, Y_inv) {
  #'@description Retrieves sample orders foR heat map visualization
  #'@describeIn Hinoue et al. online tutorial
  #'@param rpmmClusters data.frame with row.names = sample ID & 1 column named RPMM with cluster assignments
  #'@param Y_inv Input matrix for RPMM computation
  
  sampOrder <- c();
  for(r in names(table(rpmmClusters$RPMM))) {
    samps <- rownames(rpmmClusters)[rpmmClusters$RPMM == r];
    clu <- t(Y_inv[samps, ]);
    s_i <- seriation::seriate(clu, margin=2);
    so_i <- seriation::get_order(s_i);
    sampOrder <- c(sampOrder, samps[so_i]);
  }
  sampOrder <- data.frame(
    Sample_Name = sampOrder,
    RPMMSampleOrder = 1:length(sampOrder)
  );
  return(sampOrder);
}

main <- function() {
  ## RPMM with max_level of 2:
  load_Yinv();
  sup_rpmm <- blcTree(Y_inv, verbose=1, maxlevel=2-1); 
  print(sup_rpmm);
  
  ## Extract RPMM cluster labels:
  rpmmClusters <- getRPMMClustLabels(sup_rpmm, Y_inv);
  
  ## Retrieve RPMM sample order:
  sampOrders <- getRPMMSampOrder(rpmmClusters, Y_inv);  
  
  ## Export:
  dat_sup <- merge(rpmmClusters, sampOrders, by="Sample_Name");
  write.csv(dat_sup, file=OUT_FILE_NAME, row.names=FALSE, quote=FALSE);
}

main(); 

