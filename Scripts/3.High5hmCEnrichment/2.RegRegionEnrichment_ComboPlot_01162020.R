#title: "1.Genomic_Context_Analysis Tumor Samples"
#author: "nasim"
#date: "01/16/2020"
#Aims: Examine relation of Tumor Top 5% most hydroxymethylated regions with respect to CpG Islands
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

#Load 5mC data mean/median data
load("./Files/01162020_unique_hydroxy_methyl_ss.RData")

# Load Top 5% most hyperhydroxymethylated CpGs
load('./Files/01162020_top5_5hmc_unique.RData')

top5 <- top5unique
rm(top5unique)
# Load EPIC annotation file created in previous file
load('./Files/brainSpecificAnnotation.RData')


#~~~~~~~~~~~~~~~~~~~~~~~~Step 2: Modify Annotation Files~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Take out probes in sex chromosomes and cross hybridizing probes
#shortcut to match CpGs in my df of median 5hmc CpG sites since I had previously dropped those probes
epicAnnotation <- epicAnnotation[row.names(epicAnnotation) %in% hydroxy_unique_ss$ID,]

#Add Enhancer/Promoter Information using Sarah ML's code
# Change Phantom5 Enhancer Annotation
epicAnnotation$enhancer <- ifelse(grepl("chr", epicAnnotation$Phantom5_Enhancers), 1,0)

# Change Promoter Annotation
epicAnnotation$rfg_promoter <- ifelse(grepl("Promoter", epicAnnotation$Regulatory_Feature_Group), 1,0)

# Add additional annotation information
epicAnnotation$GeneBody <- ifelse(grepl("Body", epicAnnotation$UCSC_RefGene_Group),1,0)
epicAnnotation$TSS200 <- ifelse(grepl("TSS200", epicAnnotation$UCSC_RefGene_Group),1,0)
epicAnnotation$TSS1500 <- ifelse(grepl("TSS1500", epicAnnotation$UCSC_RefGene_Group),1,0)
epicAnnotation$Exon1 <- ifelse(grepl("1stExon", epicAnnotation$UCSC_RefGene_Group),1,0)
epicAnnotation$utr3 <- ifelse(grepl("3'UTR", epicAnnotation$UCSC_RefGene_Group),1,0)
epicAnnotation$utr5 <- ifelse(grepl("5'UTR", epicAnnotation$UCSC_RefGene_Group),1,0)
epicAnnotation$ExnBnd <- ifelse(grepl("ExonBnd", epicAnnotation$UCSC_RefGene_Group),1,0)
epicAnnotation$dhs <- ifelse(grepl("chr", epicAnnotation$DNase_Hypersensitivity_NAME),1,0)
epicAnnotation$TFBS <- ifelse(grepl("chr", epicAnnotation$TFBS_NAME),1,0)

#Keep the original complete version of annotation file
epicAnnotationFull <- epicAnnotation

# Subset annotation data to CpGs included in analysis
epicAnnotation <- epicAnnotationFull[row.names(epicAnnotationFull) %in% top5$ID,]

# Add results to annotation file
epicAnnotationFull$tumor_hyper <- ifelse(row.names(epicAnnotationFull) %in% row.names(epicAnnotation),1,0)

#sanity check, make sure sum of tumor_hyper column equals no. rows in top5 (37173)
sum(epicAnnotationFull$tumor_hyper)

# ===================================================================================
# Step 3: Calculate Enrichment for Hyperhydroxymethylated Sites 
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
  prop <- table(epicAnnotationFull$tumor_hyper,
                epicAnnotationFull$Relation_to_Island == i)['1','TRUE']/sum(epicAnnotationFull$tumor_hyper)
  # Calculate odds ratio
  tempOR <- oddsratio.fisher(table(epicAnnotationFull$tumor_hyper,epicAnnotationFull$Relation_to_Island == i))
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
  MH_table <- table(epicAnnotationFull$tumor_hyper, epicAnnotationFull[,i], epicAnnotationFull$Type)
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
  MH_table <- table(epicAnnotationFull$tumor_hyper, epicAnnotationFull[,i], epicAnnotationFull$Type)
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
             
             

#SAVE CPG LIST OF LOCI IN ISLANDS FOR POTENTIAL REVIEWERS
CpGIslands <- as.data.frame(epicAnnotation) %>% filter(Relation_to_Island == "Island")
write.csv(CpGIslands, file = "/Users/nasimazizgolshani/Dropbox/christensen/PCNS_Analysis/pcns_5hmC_Manuscript/05.FinalDraft_062020/AnticipateReviewers/CpGIslands281.csv")


