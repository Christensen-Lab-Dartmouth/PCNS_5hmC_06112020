# ==========================================================================================
# Pediatric CNS Tumors Hydroxymethylation
# EWAS - Visualizing Limma  Differential Analysis 
#Tumor vs NonTumor 
# Code By: Nasim
# February 04, 2020
# ==========================================================================================


# ==========================================================================================
# Step 1: Initialization
# ==========================================================================================
# Set working directory
setwd("/Users/nasimazizgolshani/Dropbox/christensen/PCNS_Analysis/PCNS")

# Load packages
rm(list=ls());
library(ggrepel);
library(matrixStats);
library(doParallel); registerDoParallel(detectCores() - 1);
library(ggplot2);

# -------------------------- Step 2. Load data & Covariates --------------------------
TumVsNon <- read.csv("./Script/15.EWAS/DMRsFinal/BetasTumVSNon_results020520.csv", stringsAsFactors = F)

# -------------------------- Step 3. Set Up Volcano Plot --------------------------

plt_istumor <- TumVsNon
plt_istumor$negLog10NominalP <- -log10(plt_istumor$P.Value);

## Gene names:
DMP.genes <- as.character(plt_istumor$UCSC_RefGene_Name); 
DMP.genes <- strsplit(DMP.genes, split=";");
plt_istumor$Gene <- rep(NA, length(DMP.genes) ); 
for(k in 1:length(DMP.genes)){
  DMP.genes[[k]] <- unique(DMP.genes[[k]]); 
  if(length(DMP.genes[[k]]) > 1) {
    plt_istumor$Gene[k] <- paste(DMP.genes[[k]], collapse="; "); 
  } else if(length(DMP.genes[[k]]) == 1) {
    plt_istumor$Gene[k] <- DMP.genes[[k]]; 
  } else if(length(DMP.genes[[k]]) == 0) {
    plt_istumor$Gene[k] <- NA;
  }
}

## Plotting parameters (Customize anyway you want):
#PTHRESH <- 0.01
#FDR_Thresh <- 0.05
COLORS <- c("gray","red","royalblue");
VALUES <- c(0.5, .8, .8); #scale
THEME <- theme_classic() + 
  theme(axis.text=element_text(size=13,color="black"), 
        axis.title=element_text(size=13,color="black"),
        title=element_text(size=15, color="black", face="bold"),
        legend.position="right", legend.title=element_blank(),legend.text=element_text(size=10,color="black"),
        strip.text.x=element_text(size=12,colour="black",face="bold"));

## Setting up labels and limits for colors of volcano
plt_istumor$dir <- ifelse(plt_istumor$logFC > 0, "Hyper-hydroxymethylated", "Hypo-hydroxymethylated");
plt_istumor$dir[plt_istumor$adj.P.Val > 0.1] <- "";

#Thresholds for Gene names
plt_istumor$Label <- plt_istumor$Gene; #initialize
#plt_istumor$Label[(plt_istumor$adj.P.Val >= 0.053)] <- NA; #higher threshold for gene labeling
plt_istumor$Label[(plt_istumor$logFC <= 0.166 & plt_istumor$logFC >= -0.28)] <- NA; #higher threshold for gene labeling
plt_istumor$Label <- ifelse(plt_istumor$adj.P.Val <0.05, plt_istumor$Gene, plt_istumor$Label)

## Data visualization: without limiting 
volcano_plot <- ggplot(plt_istumor, aes(x=logFC, y=negLog10NominalP, color=dir)) +
  geom_point(aes(size=dir, alpha=dir)) +
  scale_color_manual(values=COLORS) +
  scale_size_manual(values=VALUES) +
  scale_alpha_manual(values=VALUES) +
  geom_text_repel(aes(label=Label), color="black", size=4) +
  geom_hline(yintercept=2.72, linetype="dashed") +
  labs(x="log2FC in 5hmC Beta-value", y="-log10 P Value") +
  THEME;

#Create a function to format decimals of x-axis
fmt_dcimals <- function(decimals=0){
  # return a function responpsible for formatting the 
  # axis labels with a given number of decimals 
  function(x) as.character(round(x,decimals))
}

## Output a plot but constraining the x & y-axis limits:
final <- volcano_plot + 
  scale_x_continuous(limits=c(-0.6,0.4), breaks=seq(-0.6,0.4,by=0.1), labels = fmt_dcimals(2)) +
  scale_y_continuous(limits=c(0,7.5),breaks=seq(0,7,by=.5) ) +
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Tumor vs. Non-Tumor 5hmC ");

setwd("/Users/nasimazizgolshani/Dropbox/christensen/PCNS_Analysis/PCNS/Script/15.EWAS/SubtypeToNonDMR/ManuscriptEdition0520")

ggsave(filename = "TumVNon.png", final,
       width = 7, height = 7, dpi = 300, units = "in", device='png')

