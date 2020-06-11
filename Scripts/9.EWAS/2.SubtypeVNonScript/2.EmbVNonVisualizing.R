# ==========================================================================================
# Pediatric CNS Tumors Hydroxymethylation
# EWAS - Visualizing Limma  Differential Analysis 
#Visualizing Embryonal vs. Non-Tumor Comparisons
# Code By: Nasim
# Febraury 10, 2020
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
library(ggrepel);
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19);


# -------------------------- Step 2. Load data & Add Annotation to DHMR --------------------------
#Load Data of comparisons
EmbVNon <- read.csv("./Plots/0.052020_wPCNS0601X/SubtypeEWAS/BetasEmbVSNon_results0520.csv", stringsAsFactors = F, row.names =1)

#quantify loss of 5hmC
qLoss <- filter(EmbVNon, adj.P.Val < 0.1)
sum(qLoss$logFC < 0)
sum(qLoss$adj.P.Val < 0.02)
# -------------------------- Step 3. Visualize Q-Value --------------------------

COLORS <- c("gray","red","royalblue");
VALUES <- c(0.5, .8, 0.8); #scale

THEME <- theme_classic() + 
  theme(axis.text=element_text(size=8,color="black"), 
        axis.title=element_text(size=12,color="black"),
        title=element_text(size=12, color="black", face="bold"),
        legend.position="right", legend.title=element_blank(),legend.text=element_text(size=10,color="black"),
        strip.text.x=element_text(size=12,colour="black",face="bold"));

## If you are to use a new script, then run:
plt_istumor <- EmbVNon;
plt_istumor$dir <- ifelse(plt_istumor$logFC > 0, "Hyper-hydroxymethylated", "Hypo-hydroxymethylated");
plt_istumor$dir[plt_istumor$adj.P.Val > 0.1] <- "(Not Significant)";

## Data visualization: without limiting 
volcano_plot <- ggplot(plt_istumor, aes(x=logFC, y=negLog10FDR, color=dir)) +
  geom_point(aes(size=dir, alpha=dir)) +
  scale_color_manual(values=COLORS) +
  scale_size_manual(values=VALUES) +
  scale_alpha_manual(values=VALUES) +
  geom_hline(yintercept=-log10(0.1), linetype="dashed") +
  labs(x="log2FC in 5hmC Beta-value", y="-log10 FDR") +
  THEME;

## Output a plot but constraining the x & y-axis limits:
volcano_plot + 
#  scale_x_continuous(limits=c(-0.5,0.4), breaks=seq(-0.4,0.4,by=.1)) +
 # scale_y_continuous(limits=c(0,2)) +
  ggtitle("Embryonal Tumor vs. Non-Tumor 5hmC ");


# -------------------------- Step 4. Data Visualization of unadjusted P. Values --------------------------
plt_istumor2 <- EmbVNon
plt_istumor2$negLog10NominalP <- -log10(plt_istumor$P.Value);

## Plotting parameters (Customize anyway you want):
COLORS <- c("gray","red","royalblue");
VALUES <- c(0.5, .8, .8); #scale
THEME <- theme_classic() + 
  theme(axis.text=element_text(size=9,color="black"), 
        axis.title=element_text(size=12,color="black"),
        title=element_text(size=15, color="black", face="bold"),
        legend.position="right", legend.title=element_blank(),legend.text=element_text(size=10,color="black"),
        strip.text.x=element_text(size=12,colour="black",face="bold"));

## If you are to use a new script, then run:
plt_istumor2$dir <- ifelse(plt_istumor2$logFC > 0, "Hyper-hydroxymethylated", "Hypo-hydroxymethylated");
plt_istumor2$dir[plt_istumor2$adj.P.Val > 0.1] <- "(Not Significant)";

## Data visualization: without limiting 
volcano_plot <- ggplot(plt_istumor2, aes(x=logFC, y=negLog10NominalP, color=dir)) +
  geom_point(aes(size=dir, alpha=dir)) +
  scale_color_manual(values=COLORS) +
  scale_size_manual(values=VALUES) +
  scale_alpha_manual(values=VALUES) +
  #  geom_hline(yintercept=-log10(0.1), linetype="dashed") +
  labs(x="log2FC in 5hmC Beta-value", y="-log10 P Value") +
  THEME;


## Output a plot but constraining the x & y-axis limits:
volcano_plot + 
  #scale_x_continuous(limits=c(-0.5,0.4), breaks=seq(-0.4,0.4,by=.1)) +
 # scale_y_continuous(limits=c(0,6)) +
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Embryonal vs. Non-Tumor 5hmC ");


# -------------------------- Step 5. Set Up Labeled Volcano Plot --------------------------

plt_istumor <- EmbVNon
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

COLORS <- c("gray","red","royalblue");
VALUES <- c(0.5, .8, .8); #scale
THEME <- theme_classic() + 
  theme(axis.text=element_text(size=13,color="black"), 
        axis.title=element_text(size=13,color="black"),
        title=element_text(size=15, color="black", face="bold"),
        legend.position="none",
        #legend.position="right", legend.title=element_blank(),legend.text=element_text(size=10,color="black"),
        strip.text.x=element_text(size=12,colour="black",face="bold"));

## Setting up labels and limits for colors of volcano
plt_istumor$dir <- ifelse(plt_istumor$logFC > 0, "Hyper-hydroxymethylated", "Hypo-hydroxymethylated");
plt_istumor$dir[plt_istumor$adj.P.Val > 0.1] <- "";

#Thresholds for Gene names
plt_istumor$Label <- plt_istumor$Gene; #initialize
plt_istumor$Label[(plt_istumor$adj.P.Val >= 0.01384207)] <- NA; #higher threshold for gene labeling
plt_istumor$Label[(plt_istumor$logFC >= -0.3)] <- NA; #higher threshold for gene labeling

#plt_istumor$Label[(plt_istumor$logFC <= 0.3 & plt_istumor$logFC >= -0.4)] <- NA; #higher threshold for gene labeling
#plt_istumor$Label <- ifelse(plt_istumor$adj.P.Val < 0.0061, plt_istumor$Gene, plt_istumor$Label)

## Data visualization: without limiting 
volcano_plot <- ggplot(plt_istumor, aes(x=logFC, y=negLog10NominalP, color=dir)) +
  geom_point(aes(size=dir, alpha=dir)) +
  scale_color_manual(values=COLORS) +
  scale_size_manual(values=VALUES) +
  scale_alpha_manual(values=VALUES) +
  geom_text_repel(aes(label=Label), color="black", size=4) +
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
  scale_x_continuous(limits=c(-0.6,0.4), breaks=seq(-0.6,0.4,by=0.2), labels = fmt_dcimals(2)) +
  scale_y_continuous(limits=c(0,7.5),breaks=seq(0,7.5,by=.5) ) +
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Embryonal vs. Non-Tumor 5hmC ");

print(final)
setwd("/Users/nasimazizgolshani/Dropbox/christensen/PCNS_Analysis/PCNS/Plots/0.052020_wPCNS0601X/SubtypeEWAS")
ggsave(filename = "EmbVNon.png", final,
       width = 9, height = 7, dpi = 300, units = "in", device='png')
