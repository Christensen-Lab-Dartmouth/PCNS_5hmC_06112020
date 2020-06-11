#title: Genomic Context of DHMRs
#author: "nasim"
#date: "02/07/2020"
#Aims: Examine relation of Tumor DHMRs with respect to CpG Islands
#Find Enrichment at CpG Islands and Gene Regulatory Regions

setwd("/Users/nasimazizgolshani/Dropbox/christensen/PCNS_Analysis/PCNS")

# ==========================================================================================
# Step 1: Initialization
# ==========================================================================================
rm(list=ls())
# Load Required Packages
library(reshape)
library(ggpubr)
library(epitools) #OR calculation
library(forestplot)
library(tidyverse)
library(GenomicRanges)
library(genomation)
library(scales)
library(grid)
library(gridExtra)

# -------------------------- Step 2. Load data  --------------------------
TumVsNon <- read.csv("./Script/15.EWAS/DMRsFinal/BetasTumVSNon_results020520.csv", 
                     stringsAsFactors = F, row.names = 1)


# Load Top 5% most hyperhydroxymethylated CpGs
load('./Files/01162020_top5_5hmc_unique.RData')

# Load EPIC annotation file created in previous file
load('./Files/brainSpecificAnnotationComplete.RData')


#~~~~~~~~~~~~~~~~~~~~~~~~Step 3: Modify Annotation Files~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Subset the universe (aka epic annotation) to only top 5% CpGs as DMRs were extracted from top 5% list not 850k
epicAnnotation <- epicAnnotation[row.names(epicAnnotation) %in% top5unique$ID,]

# Add Phantom4  Annotation (This entails information on promoter structures/reg transcription factors/ miRNAs)
epicAnnotation$fantom4 <- ifelse(grepl("CpG", epicAnnotation$Phantom4_Enhancers), 1,0)

# Add additional annotation information
epicAnnotation$TFBS <- ifelse(grepl("chr", epicAnnotation$TFBS_NAME),1,0)

# Add results to annotation file (HDHMR=hypo-hydroxymethylated region)
HDHMR <- TumVsNon[TumVsNon$adj.P.Val <0.1, ]
epicAnnotation$HDHMR <- ifelse(row.names(epicAnnotation) %in% row.names(HDHMR),1,0)
#check
sum(epicAnnotation$HDHMR) #726 (should equal nrow(HDHMR))


# ===================================================================================
# Step 4: Calculate Enrichment for Hyperhydroxymethylated Sites 
# ===================================================================================
#Look at Enrichment with Relation to CpG Islands
#first set relation to island as factor 
epicAnnotation$Relation_to_Island <- as.factor(epicAnnotation$Relation_to_Island)

# Create OR table for different distances to CpG Islands
OR_Table_island <- as.data.frame(matrix(ncol=5, nrow=0))
colnames(OR_Table_island) <- c("proportion", "OR", "lower.CI", "upper.CI", "pVal")
islandLevels <- c("N_Shelf", "N_Shore", "Island", "S_Shore", "S_Shelf", "OpenSea")
islandLabels <- c("North Shelf", "North Shore", "CpG Island", "South Shore", "South Shelf", "Open Sea")

for (i in islandLevels) {
  # Calculate proportion of sites in that context
  prop <- table(epicAnnotation$HDHMR,
                epicAnnotation$Relation_to_Island == i)['1','TRUE']/sum(epicAnnotation$HDHMR)
  # Calculate odds ratio
  tempOR <- oddsratio.fisher(table(epicAnnotation$HDHMR,epicAnnotation$Relation_to_Island == i))
  # Add values to output
  OR_Table_island[i,] <- c(prop,tempOR$measure[2,],tempOR$p.value[2,'fisher.exact'])
}

orText_island <- cbind(c('Relation to CpG Island', islandLabels),
                       c('Odds Ratio (95% CI)', paste(format(round(OR_Table_island$OR,2), nsmall= 2),' (',format(round(OR_Table_island$lower.CI,2), nsmall=2), '-',
                                                      format(round(OR_Table_island$upper.CI,2), nsmall=2),')', sep= '')),
                       c('P value', formatC(OR_Table_island$pVal, digits=2, format= 'E')))


#Now repeat for Regulator Regions

# Enrichment for Genomic Context Only
OR_Table_RR <- as.data.frame(matrix(ncol= 4, nrow= 0))
colnames(OR_Table_RR) <- c("OR", "lower.CI", "upper.CI", "pVal")
contextList <- c("enhancer","TFBS", "utr5",  "promoters", "exon","intron", "utr3", "intergenic", "dhs")
OR_Table_RR.names <- c("Enhancer","TFBS","UTR-5",  "Promoter", "Exon","Intron", "UTR-3", "Intergenic", "DHS")
for(i in contextList){
  MH_table <- table(epicAnnotation$HDHMR, epicAnnotation[,i], epicAnnotation$Type)
  MH_test <- mantelhaen.test(MH_table, exact= T)
  OR_Table_RR[i, ] <- c(MH_test$estimate, MH_test$conf.int, MH_test$p.value)
}


orText_RR <- cbind(c('Regulatory Region', OR_Table_RR.names),
                   c('Odds Ratio (95% CI)', paste(format(round(OR_Table_RR$OR,2), nsmall= 2),' (',format(round(OR_Table_RR$lower.CI,2), nsmall=2), '-',
                                                  format(round(OR_Table_RR$upper.CI,2), nsmall=2),')', sep= '')),
                   c('P value', formatC(OR_Table_RR$pVal, digits=2, format= 'E')))

#Now repeat for super-enhancers

#Repeat for Super Enhancers
# Enrichment for Super Enhancers 
OR_Table_Super <- as.data.frame(matrix(ncol= 4, nrow= 0))
colnames(OR_Table_Super) <- c("OR", "lower.CI", "upper.CI", "pVal")
contextList <- c("astro_SuperEnhancer", "bag_SuperEnhancer", "bac_SuperEnhancer",
                 "bcg_SuperEnhancer","hip_SuperEnhancer", "templobe_SuperEnhancer","frontlobe_SuperEnhancer")
OR_Table_S.names <- c("Astrocytes", "Anterior Gyrus", "Anterior Caudate",
                      "Cingulate Gyrus","Hippocampus", "Temporal Lobe","Frontal Lobe")
for(i in contextList){
  MH_table <- table(epicAnnotation$HDHMR, epicAnnotation[,i], epicAnnotation$Type)
  MH_test <- mantelhaen.test(MH_table, exact= T)
  OR_Table_Super[i, ] <- c(MH_test$estimate, MH_test$conf.int, MH_test$p.value)
}

orText_super <- cbind(c('Super Enhancer ', OR_Table_S.names),
                      c('Odds Ratio (95% CI)', paste(format(round(OR_Table_Super$OR,2), nsmall= 2),' (',format(round(OR_Table_Super$lower.CI,2), nsmall=2), '-',
                                                     format(round(OR_Table_Super$upper.CI,2), nsmall=2),')', sep= '')),
                      c('P value', formatC(OR_Table_Super$pVal, digits=2, format= 'E')))



# ===================================================================================
# Step 4: Create 3 Forest Plots in a Panel
# ===================================================================================
#
#Plot of Enrichment with respect ot CpG Island
forestplot(labeltext= orText_island, graph.pos = 2,
           new_page = TRUE,
           mean= c(NA, OR_Table_island$OR), 
           lower= c(NA, OR_Table_island$lower.CI), 
           upper= c(NA, OR_Table_island$upper.CI), 
           graphwidth= unit(100, 'mm'),
           colgap= unit(0.05, 'mm'),
           xlab= "Odds Ratio",
           txt_gp=fpTxtGp(label=gpar(cex=0.90),
                          ticks=gpar(cex=1.1),
                          xlab=gpar(cex = 0.9),
                          legend = gpar(cex = .9)),
           lwd.zero = gpar(lwd= 2),
           lwd.xaxis = gpar(lwd= 2),
           lwd.ci = gpar(lwd=2),
           zero= 1,
           ci.vertices= TRUE,
           is.summary=c(TRUE, rep(FALSE,9)),
           xlog= T,
           col= fpColors(box= "black", line= "black", zero= "black"),
           boxsize = 0.20,
           cex=2,
           fn.ci_norm = fpDrawNormalCI, pch= 2)
plot.1 <- grid.grab()

#Plot of Enrichment of Regulatory Region Enrichment
forestplot(labeltext= orText_RR, graph.pos = 2,
           new_page = TRUE,
           mean= c(NA, OR_Table_RR$OR), 
           lower= c(NA, OR_Table_RR$lower.CI), 
           upper= c(NA, OR_Table_RR$upper.CI), 
           graphwidth= unit(100, 'mm'),
           colgap= unit(0.05, 'mm'),
           xlab= "Odds Ratio",
           txt_gp=fpTxtGp(label=gpar(cex=1.0),
                          ticks=gpar(cex=1.1),
                          xlab=gpar(cex = 1.0),
                          legend = gpar(cex = 1.0)),
           lwd.zero = gpar(lwd= 2),
           lwd.xaxis = gpar(lwd= 2),
           lwd.ci = gpar(lwd=2),
           zero= 1,
           ci.vertices= TRUE,
           is.summary=c(TRUE, rep(FALSE,9)),
           xlog= T,
           col= fpColors(box= "black", line= "black", zero= "black"),
           boxsize = 0.20,
           cex=2,
           xticks = c(0.5, 1,2.5),
           fn.ci_norm = fpDrawNormalCI, pch= 2)
plot.2 <- grid.grab()

#Plot of Enrichment of in Super Enhancers
forestplot(labeltext= orText_super, graph.pos = 2,
           new_page = TRUE,
           mean= c(NA, OR_Table_Super$OR), 
           lower= c(NA, OR_Table_Super$lower.CI), 
           upper= c(NA, OR_Table_Super$upper.CI), 
           graphwidth= unit(100, 'mm'),
           colgap= unit(0.05, 'mm'),
           xlab= "Odds Ratio",
           txt_gp=fpTxtGp(label=gpar(cex=1.0),
                          ticks=gpar(cex=1.1),
                          xlab=gpar(cex = 1.0),
                          legend = gpar(cex = 1.0)),
           lwd.zero = gpar(lwd= 2),
           lwd.xaxis = gpar(lwd= 2),
           lwd.ci = gpar(lwd=2),
           zero= 1,
           ci.vertices= TRUE,
           is.summary=c(TRUE, rep(FALSE,9)),
           xlog= T,
           col= fpColors(box= "black", line= "black", zero= "black"),
           boxsize = 0.20,
           cex=2,
           xticks = c(0.5, 1,2.5),
           fn.ci_norm = fpDrawNormalCI, pch= 2)
plot.3 <- grid.grab()

grid.newpage()
grid.arrange(plot.1, plot.2, plot.3, ncol=1)

grid.arrange(arrangeGrob(plot.2,plot.3, ncol=1, nrow=2),
             arrangeGrob(plot.1, ncol=1, nrow=1, widths = 2), heights=c(2.5,1))    
