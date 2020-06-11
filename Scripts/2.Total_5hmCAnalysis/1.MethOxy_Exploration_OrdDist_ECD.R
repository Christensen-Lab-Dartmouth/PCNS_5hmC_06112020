#Pediatric CNS and 5hmC
#Exploration of Data
#Examine ECD/Distribution of means in 5hmC/5mC tumors vs ctls
#Nasim Azizgolshani

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#Initialization#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list=ls())
#Set working directory
setwd("/Users/nasimazizgolshani/Dropbox/christensen/PCNS_Analysis/PCNS")

#Load the processed Betas w/o sex chromosomes and cross reactive probes
load("./Files/01142020_PCNS_MethOxy_GoodProbes.RData")
#load useful functions in utils 
source("./Script/1.Processing/2.PCNS_Utils.R")
#load metadata
targets <- read.csv("./Files/Metadata_complete.csv")

#Load necessary libraries
library(tidyverse)
library(ggsci)
library(lattice)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#Step 2: Plot ECD of Tumors#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Change column names in dataframes to patient ID from ArrayID
colnames(PCNS_5hmc) <- targets$Label.ID.[match(names(PCNS_5hmc),targets$ArrayID)]
colnames(PCNS_5mc) <- targets$Label.ID.[match(names(PCNS_5mc),targets$ArrayID)]


#remove repeat samples from individuals
#Identify serial samples and multiregional samples (patients with more than one sample)
serials <- c("PCNS0302X", "PCNS0603X")

multiregional <- c("PCNS1001A", "PCNS0301B", "PCNS0301C", "PCNS0601B","PCNS0601X", "PCNS0901B", "PCNS1901C")


#subset out these samples from analysis
PCNS_5hmc_noserials <- PCNS_5hmc[,!colnames(PCNS_5hmc) %in% serials]
PCNS_5hmc_unique <- PCNS_5hmc_noserials[,!colnames(PCNS_5hmc_noserials) %in% multiregional]
#drop aberrant control
PCNS_5hmc_unique <- PCNS_5hmc_unique[,-31]

#repeat for mC
PCNS_5mc_noserials <- PCNS_5mc[,!colnames(PCNS_5mc) %in% serials]
PCNS_5mc_unique <- PCNS_5mc_noserials[,!colnames(PCNS_5mc_noserials) %in% multiregional]
PCNS_5mc_unique <- PCNS_5mc_unique[,-31]

#Need to Subset for Tumors Only
TumorTargets <- targets %>% filter(Tumor_Type != "Non-Tumor")
PCNS_5hmc_Tumor <- PCNS_5hmc_unique[colnames(PCNS_5hmc_unique) %in% TumorTargets$Label.ID.]
PCNS_5mc_Tumor <- PCNS_5mc_unique[colnames(PCNS_5mc_unique) %in% TumorTargets$Label.ID.]

Non_Tumor <- targets %>% filter(Tumor_Type == "Non-Tumor")
PCNS_5hmc_Ctl <- PCNS_5hmc_unique[colnames(PCNS_5hmc_unique) %in% Non_Tumor$Label.ID.]
PCNS_5mc_Ctl <- PCNS_5mc_unique[colnames(PCNS_5mc_unique) %in% Non_Tumor$Label.ID.]

#to avoid mix-ups can remove original dataframes
#rm(PCNS_5c, PCNS_5hmc,PCNS_5mc)

#ECD by Median
plotEcdfByMedian <- function(PCNS_5hmc_Tumor, PCNS_5mc_Tumor, cpgSubset=NULL) {
  #'@description Plot cumulative densities of 5mC & 5hmC based on *median* ECDF
  #'@source Christensen-Lab-Dartmouth/Normal-breast-5hmC
  if(is.null(cpgSubset)) {
    sel <- 1:dim(PCNS_5hmc_Tumor)[1];
  } else {
    sel <- cpgSubset;
  }
  
  mdbeta1 <- apply(PCNS_5hmc_Tumor, 1, median, na.rm=TRUE); #5hmC
  mdbeta2 <- apply(PCNS_5mc_Tumor, 1, median, na.rm=TRUE); #5mC
  
  ecdfMed1 <- ecdf(mdbeta1);
  ecdfMed2 <- ecdf(mdbeta2);
  
  plot(ecdfMed2, main="Cumulative Density 5(h)mC in Tumor Samples", xlab="Median average beta",ylab="Cumulative proportion", col="blue", cex.main= 1.6 ,cex.lab=1.3, cex.axis=1.15, las=1, bty="L"); 
  plot(ecdfMed1, add=TRUE, col="red");
  legend("bottomright", c("5mC","5hmC"), lwd=3, col=c("blue","red"), bty = "n", cex = 1.6);
}


plotEcdfByMedian(PCNS_5hmc_Tumor, PCNS_5mc_Tumor, cpgSubset=NULL)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#Step 3: Plot Cummulative Density correlation between 5hmC and 5mC #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#I reload the original files (prior to being subset to non-sex/non-SNP prones) 
#because the format of this function requires it to be in matrix form

load("./Files/12.04.Tandem_OxBS_Betas/12062019_OxyBSFunnorm_PCNS.RData")

#subset metadata to unique tumors
TumorUnique <- TumorTargets %>% filter(Label.ID. %in% colnames(PCNS_5hmc_Tumor))
#subset array data to unique tumors
MethOxyTumors <- MethOxy[,colnames(MethOxy) %in% TumorUnique$ArrayID, ]
#subset to "good probes"
MethOxyTumors <- MethOxyTumors[row.names(MethOxyTumors) %in% row.names(PCNS_5hmc),,]
dim(MethOxyTumors) 

#Write function to determine and plot correlation
plotCorsBeta <- function(MethOxy, cpgSubset=NULL) {
  #'@description Plot correlation between 5hmC and mC
  #'@source Christensen-Lab-Dartmouth/Normal-breast-5hmC
  if(is.null(cpgSubset)) {
    sel <- 1:dim(MethOxy)[1];
  } else {
    sel <- cpgSubset;
  }
  
  corbeta12 <- apply(MethOxy[sel,,2:3], 1, cor, use="pair", method="pearson")[2, ]; 
  corsbeta12 <- apply(MethOxy[sel,,2:3], 1, cor, use="pair", method="spearman")[2, ];
  ecdfCorrP <- ecdf(corbeta12);
  ecdfCorrS <- ecdf(corsbeta12);
  
  plot(ecdfCorrP, main="Cumulative Density Correlation", xlab="Correlation coefficient", ylab="Cumulative proportion", col="black", cex.main= 1.6 ,cex.lab=1.3, cex.axis=1.15, las=1, bty="L"); 
  plot(ecdfCorrS, add=TRUE, col="purple"); 
  legend("topleft", c("Pearson", "Spearman"), lwd=3, col=c("black","purple"), bty = "n", cex = 1.3);
  abline(v=0, lty=3);
}

plotCorsBeta(MethOxyTumors, cpgSubset=NULL)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Step 4: Ordered Distribution #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Find Medians 
md.beta1 <- apply(PCNS_5hmc_Tumor, 1, median, na.rm=TRUE); #5hmC

#Ordered distribution by median 
md_5hmc_ordered <- md.beta1[order(md.beta1, decreasing=FALSE)]

#Find percentiles
md_5hmc_ordered <- as.data.frame(md_5hmc_ordered)  
colnames(md_5hmc_ordered) <- "Ordered_BetaValue"

md_5hmc_ordered <- md_5hmc_ordered %>% mutate(PCT = ntile(Ordered_BetaValue, 100))

#find minimum beta value of 95 percentile
PCT95 <- md_5hmc_ordered %>% filter(PCT > 95)
# produce ordered distribution plot w/ lines to indicate beta value signals of 5, 10 & 15%
ppi =300
png("./Plots/1.Plots_from_Initial_Probe_selection/3.OrderedDist/ordered_5hmC_distributionNoRptsTumors.png", width=5*ppi, height=5*ppi, res=ppi)
plot(md_5hmc_ordered$PCT, md_5hmc_ordered$Ordered_BetaValue, pch=16, bty='l', las = 1, ylab = "Median 5hmC beta-value", xlab = "Percentile Rank of CpG loci", cex.lab = 1.3, cex.axis = 1.15)
#abline(a = 0.05, b = 0, col = "purple", lwd = 2, lty = "dotdash")
abline(a = 0.093, b = 0, col = "blue", lwd = 2, lty = "dotdash")
#abline(a = 0.15, b = 0, col = "red", lwd = 2, lty = "dotdash")
abline(v=95, lwd =2, col = "gray60", lty = "dotdash")
dev.off()



# produce ordered distribution plot, highlighting the high 5hmC CpGs (top 1% median)
# png("/Users/nasimazizgolshani/Dropbox/christensen/PCNS_Analysis/PCNS/Plots/ECDF/ordered_5hmC_distribution_high_5hmC.png", width=5*ppi, height=5*ppi, res=ppi)
# plot(mean_5hmc_ordered, pch=16, xaxt="n", bty='l', las = 1, ylab = "Median 5hmC beta-value", xlab = "Ordered CpG loci", cex.lab = 1.3, cex.axis = 1.15)
# points(x=CpG_begin:CpG_end, mean_5hmc_ordered[CpG_begin:CpG_end], col='red', pch=16)
# abline(a = 0.141132, b = 0, col = "green", lwd = 2, lty = "dotdash")
# dev.off()
# 

#Stored Code:
#To plot by mean
plotEcdfByMean <- function(PCNS_5hmc_Tumor, PCNS_5mc_Tumor, cpgSubset=NULL) {
  #'@description Plot cumulative densities of 5mC & 5hmC based on *mean* ECDF
  #'@source Christensen-Lab-Dartmouth/Normal-breast-5hmC
  if(is.null(cpgSubset)) {
    sel <- 1:dim(PCNS_5hmc_Tumor)[1];
  } else {
    sel <- cpgSubset;
  }
  
  mean.beta1 <- apply(PCNS_5hmc_Tumor, 1, mean, na.rm=TRUE); #5hmC
  mean.beta2 <- apply(PCNS_5mc_Tumor, 1, mean, na.rm=TRUE); #5mC
  
  ecdfMean1 <- ecdf(mean.beta1);
  ecdfMean2 <- ecdf(mean.beta2);
  
  plot(ecdfMean2, main="Cumulative Density 5(h)mC", xlab="Mean average beta",ylab="Cumulative proportion", col="blue", cex.lab=1.3, cex.axis=1.15, las=1); 
  plot(ecdfMean1, add=TRUE, col="red");
  legend("bottomright", c("5mC","5hmC"), lwd=3, col=c("blue","red"));
}


p2 <- plotEcdfByMean(PCNS_5hmc_Tumor, PCNS_5mc_Tumor, cpgSubset = NULL)
# Define grid layout to locate plots and print each graph
library(lattice)
pushViewport(viewport(layout = grid.layout(1, 2)))
print(p1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))



