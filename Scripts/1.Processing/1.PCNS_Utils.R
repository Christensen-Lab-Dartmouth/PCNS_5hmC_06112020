# R Helper Methods for Preprocessing & Extracting 5hmC Data
# Script maintainer: Nasim
# Date: (ongoing)
# Reference: David Chen Brestmilk Utils

library(OxyBS);
library(minfi);
if(c("doParallel") %in% installed.packages()) {
  print("Loading optional packages...");
  library(doParallel);
  registerDoParallel(detectCores() - 1);
}

prepareListOfArrays <- function(path, meta) {
  #'@description Initializes R lists of two sub-lists, one for all *bisulfite (BS)* the other *oxidative bisulfite (OxBS)*
  #'@param path Path of IDAT fi
  #'@param meta R data.frame that must have columns `(int) OX_BS`, `(string) Label.ID.`, `(string) ArrayID`, `(int) Pair`
  #'@source Christensen-Lab-Dartmouth/Breastmilk Pilot Project by David Chen
  
  # Find numeric files/Chips:
  chips <- list.files(path=IDAT_DIR, pattern="idat");
  chips <- gsub("\\_.*", "", chips);
  
  ## Initialize R.list associated with IDAT file:
  arrayList <- list();
  for(chip in chips) {
    ff <- grep("Grn.idat$", list.files(pattern=chip), value=TRUE);
    arrayList[[chip]] <- sort(gsub("[_]Grn.idat", "", ff));
  }
  
  ## Select arrays with BS conversion:
  arraysBS <- intersect(
    unlist(arrayList),
    meta$ArrayID[meta$OX_BS == 0]
  );
  
  ## Select arrays with oxBS conversion:
  arraysOxBS <- intersect(
    unlist(arrayList),
    meta$ArrayID[meta$OX_BS == 1]
  );
  
  return(list(arraysBS=arraysBS, arraysOxBS=arraysOxBS));
}




constrainSignals <- function(mat) {
  #'@description Helper function to handle non-finite values by resetting to finite, max values
  #'@param mat Matrix of signals (non-zero)
  stopifnot(all(mat >= 0)); #checkpoint
  infMask <- mat == Inf;
  zeroMask <- mat == 0;
  
  if(sum(infMask) > 0) {
    noInf <- as.numeric(mat);
    noInf <- noInf[is.finite(noInf)];
    mat[infMask] <- max(noInf);
  }
  
  if(sum(zeroMask) > 0) {
    mat[zeroMask] <- min(noInf);
  }
  
  return(mat);
} 


constrainBetas <- function(mat) {
  #'@description Helper function to handle atypical beta-values by constraining strictly between 0 & 1
  #'@param mat Matrix of beta values
  mat[mat == 0] <- 0.00001;
  mat[mat == 1] <- 0.99999;
  return(mat);
}


preprocessFunnormRedGreen <- function (rgSet, nPCs=2, sex=NULL, verbose=TRUE) {
  #'@description  Andy rewrote the FunNorm function so that he could pull out the functional-normalization ...
  #'@description  signals (as opposed to betas and M-values). Later Andy will use these combined signals to ...
  #'@description  estimate 5mC and 5hmC using a maximum-likelihood approach
  #'@author A. Houseman
  #'@source Christensen-Lab-Dartmouth/Normal-breast-5hmC by owenwilkins
  #'@param rgSet minfi-extracted object
  #'@param nPCs,sex Default params 
  
  # minfi:::.isRG(rgSet) #deprecated
  rgSet <- updateObject(rgSet); 
  if(verbose) cat("[preprocessFunnorm] Mapping to genome\n");
  
  gmSet <- mapToGenome(rgSet); 
  subverbose <- max(as.integer(verbose) - 1L, 0); 
  
  if(verbose) cat("[preprocessFunnorm] Quantile extraction\n");
  extractedData <- minfi:::.extractFromRGSet450k(rgSet); 
  
  if(is.null(sex)) {
   gmSet <- addSex(gmSet, getSex(gmSet, cutoff = -3));
    sex <- rep(1L, length(gmSet$predictedSex)); 
   sex[gmSet$predictedSex == "F"] <- 2L; 
  }
  rm(rgSet); #frees up memory; not strictly required
  
  if (verbose) cat("[preprocessFunnorm] Normalization\n");
   CN <- getCN(gmSet); #not used
  
  return(minfi:::.normalizeFunnorm450k(
    object = gmSet, 
    extractedData = extractedData,
    sex = sex, 
    nPCs = nPCs, 
    verbose = subverbose
  ));
}


computeMethAndUnmeth <- function(arraysBS, arraysOxBS,
                                 dropTechProbes=TRUE,
                                 dropSnpProbes=TRUE) {
  ## Use minfi to extract methylated vs. unmethylated signals:
  datListBS <- minfi::read.metharray(arraysBS, verbose=TRUE);
  datListOxBS <- minfi::read.metharray(arraysOxBS,verbose=TRUE); 
  
  ## Funnorm-preprocess using A. Houseman wrapper method:
  rgBS <- preprocessFunnormRedGreen(datListBS);
  rgOxBS <- preprocessFunnormRedGreen(datListOxBS); 
  
  #this would be a good point to do some QC etc
  #can look at Meghan's code
  
  
  ## Optional procedures with callbacks:
  if(dropTechProbes) {
    print("Dropping Technical Probes...");
    #Here you are dropping the rownames in the annotation... 
    #that have either "rs" (the control 59 probes in EPIC) and the "ch" CpH probes according to the names.
    #CpH is any non-CpG (A, T, C)
    rgBS <- dropMethylationLoci(rgBS);
    rgOxBS <- dropMethylationLoci(rgOxBS);
  } else {
    print("Technical probes NOT removed! Are you sure?");
  }
  
  if(dropSnpProbes) {
    print("Dropping Technical Probes with default MAF...");
    #This is based on the Illumina annotation file with the maf for the probe, CpG or SBE-single base extension)
    rgBS <- dropLociWithSnps(rgBS);
    rgOxBS <- dropLociWithSnps(rgOxBS);
  } else {
    print("SNP probes NOT removed! Are you sure?");
  }
  
  ## Extract beta-values & apply constraints:
  methBS <- constrainBetas(assay(rgBS, "Meth"));
  methOxBS <- constrainBetas(assay(rgOxBS, "Meth"));
  unmethBS <- constrainBetas(assay(rgBS, "Unmeth")); 
  unmethOxBS <- constrainBetas(assay(rgOxBS, "Unmeth")); 
  
  return(list(
    ## Signals for methylated channels:
    methBS = methBS,
    methOxBS = methOxBS,
    
    ## Signals for unmethylated channels:
    unmethBS = unmethBS,
    unmethOxBS = unmethOxBS
  ));
}

#Defines idBS and idOxBS in processing
computeAggregatedSignals <- function(methBS,methOxBS, idBS,idOxBS, rowsubset=NULL) {
  #'@description Companion function to `computeBetas`
  if(is.null(rowsubset)) rowsubset <- 1:dim(methBS)[1];
  
  signalBS <- methBS[rowsubset,idBS] + unmethBS[rowsubset,idBS];
  signalOxBS <- methOxBS[rowsubset,idOxBS] + unmethOxBS[rowsubset,idOxBS];
  
  signalBS <- constrainSignals(signalBS);
  signalOxBS <- constrainSignals(signalOxBS);
  
  return(list(signalBS=signalBS, signalOxBS=signalOxBS));
}


computeBetasByChemTx <- function(methBS,methOxBS, signalBS,signalOxBS, idBS,idOxBS, rowsubset=NULL) {
  #'@description Computes beta values for BS and OxBS-treated specimens
  #'@source Christensen-Lab-Dartmouth/Normal-breast-5hmC
  #'@source OxyBS R package
  if(is.null(rowsubset)) rowsubset <- 1:dim(methBS)[1];
  
  betaBS <- methBS[rowsubset,idBS] / signalBS[,idBS];
  betaOxBS <- methOxBS[rowsubset,idOxBS] / signalOxBS[,idOxBS];
  return(list(betaBS=betaBS, betaOxBS=betaOxBS));
}


processSignals <- function(betaBS, betaOxBS, signalBS, signalOxBS) {
  #'@description Processes minfi-extracted tandem-array signals
  #'@describeIn Christensen-Lab-Dartmouth/Normal-breast-5hmC;
  #'@describeIn OxyBS R pakage documentation
  
  ## Checkpoints:
  stopifnot(identical(rownames(betaBS), rownames(signalBS)));
  stopifnot(identical(rownames(betaOxBS), rownames(signalOxBS)));
  stopifnot(identical(colnames(betaBS), colnames(signalBS)));
  stopifnot(identical(colnames(betaOxBS), colnames(signalOxBS)));
  
  nCpGs <- nrow(betaBS);
  nSpecimens <- ncol(betaBS);
  
  ## Initialize C/mC/hmC data container:
  MethOxy <- array(NA, dim=c(nCpGs, nSpecimens, 3));
  dimnames(MethOxy) <- list(
    rownames(betaBS)[1:nCpGs],
    colnames(betaBS)[1:nSpecimens],
    c("C","5mC","5hmC")
  );
  
  ## Fit Maximum Likelihood model & update data container:
  print("************* Running maximum-likelihood estimation of beta-values... *************");
  for(i in 1:nSpecimens){
    print(paste0("Fitting unique sample ", i, " of ", nSpecimens));
    MethOxy[, i, ] <- fitOxBS(betaBS[,i], betaOxBS[,i], signalBS[,i], signalOxBS[,i]);
  }
  
  ## Check:
  print("Check below to see whether first two dimensions sum to 1:");
  print(table(apply(MethOxy,1:2,sum)));
  
  return(MethOxy);
}

plotVarDistr <- function(dat, pltTheme=NULL, pltTitle=NULL) {
  #'@description Plots the variance distribution of a matrix
  vars <- matrixStats::rowVars(dat);
  vars <- sort(vars, decreasing=TRUE);
  
  plt <- qplot(seq_along(vars), vars, size=I(0.1)) +
    labs(x="Index", y="Variance") +
    scale_y_continuous(limits=c(0,0.16));
  
  if(! is.null(pltTheme)) plt <- plt + pltTheme;
  if(! is.null(pltTitle))  plt <- plt + ggtitle(pltTitle);
  
  return(plt);
}

#I adjusted this to get KJ's top cpg sites for GBM instead
getJohnsonsCpGs <- function(path="/Users/nasimazizgolshani/Dropbox/christensen/Glioblastoma_5hmc_Learning/DavidCode/david_5hmC/johnson_supplemental_data.csv") {
  #'@description Loads Wilkins et al. 2018 Supp Data 1 resaved as CSV without legend
  #'@return A list of 3,876 CpGs
  JohnsonSd1 <- read.csv(path, header=TRUE, stringsAsFactors=FALSE);
  JohnsonCpGs <- JohnsonSd1$Illumina.CpG.ID;
  assign("JohnsonsCpGs", JohnsonCpGs, envir=.GlobalEnv);
}

getKMostVariableCpGs <- function(dat, k=3876) {
  #'@description Get K most variable CpGs
  #'@param dat Matrix of beta-values with rows=CpGs, columns=samples
  #'@param k Number of most variable. Default to 3876
  print(paste("Retrieving", k, "most variable CpGs..."));
  cpgs <- rownames(dat); 
  #rowVars gives Variance estimates for each row/column in a matrix 
  vars <-  matrixStats::rowVars(dat, na.rm=TRUE); 
  cpgs <- cpgs[order(vars, decreasing=TRUE)]; 
  return(cpgs[1:k]);
}

