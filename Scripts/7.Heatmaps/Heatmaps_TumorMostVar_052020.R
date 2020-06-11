# Pediatric CNS Tumors Methylation Heat Maps
#Based on RPMM clusters 
#Top 10,000 most variable CpGs in unique tumors
# Script author: Nasim
# Date: 04/20/20
# Notes: Source https://github.com/Christensen-Lab-Dartmouth/CF_Epigenetics/blob/master/src/01b_Various_heatmaps.R

#~~~~~~~~~~~~~~~~~~~~~~~~~#Step 1 Initialization~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list=ls());
library(ggplot2);
library(matrixStats);
library(RColorBrewer)
library(tidyverse)
library(pheatmap);
library(reshape2);
library(ggsci)
library("scales")
library(doParallel); registerDoParallel(detectCores() - 1);
setwd("/Users/nasimazizgolshani/Dropbox/christensen/PCNS_Analysis/PCNS/");


#~~~~~~~~~~~~~~~~~~~~~~~~~#Step 2 Load and Format Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#load all unique tumors
#load("./Files/021520UniqueSamples.RData")
load("./Files/01142020_PCNS_MethOxy_GoodProbes.RData")
#load targets
targets <- read.csv("./Files/Metadata_complete.csv", stringsAsFactors = FALSE, row.names = 1)
#remove OxBS duplicates from targets
serials <- c("PCNS0302X", "PCNS0602X", "PCNS0603X")

multiregional <- c("PCNS1001A", "PCNS0301B", "PCNS0301C", "PCNS0601B", "PCNS0901B", "PCNS1901C")

targets <- targets[!targets$Label.ID. %in% c(multiregional, serials), ]

targets <- targets%>% filter(OX_BS == 0)

dim(targets)
#drop the aberrant control 
targets <- targets[-31, ]

#subset 5hmC data 
All5hmc <- PCNS_5hmc[, colnames(PCNS_5hmc) %in% targets$ArrayID]
colnames(All5hmc) <- targets$Label.ID.[match(names(All5hmc),targets$ArrayID)]

#~~~~~~~~~~~~~~~~~~~~~~~~~#Step 3 Load RPMM Cluster data  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#load RPMM Clusters and add them to annotation
RPMM <- read.csv("./Script/9.RPMM/RPMM_TumorMosVar_Script/TumorsRPMM.csv", stringsAsFactors = FALSE)
RPMM2 <- mutate(RPMM, Label.ID. = Sample_Name)%>%
  mutate(Cluster = ifelse(RPMMClusters == "rL", "Low_5hmC",
                             ifelse(RPMMClusters == "rR", "High_5hmC", NA)))
RPMM2 <- RPMM2[,-c(1,2)]

#need to add ctls to RPMM data as these were not included in my RPMM clustering
#Create a dataframe with the values all as NA then add them to the tumor RPMM cluster dataframe
ctl <- targets %>% filter(Tumor_Type == "Non-Tumor")
ctlRPMM <- as.data.frame(
  Label.ID. <- ctl$Label.ID.)

colnames(ctlRPMM) <- "Label.ID."

ctlRPMM <- ctlRPMM %>% mutate(Cluster = NA) %>%
  mutate(RPMMSampleOrder = NA)
head(ctlRPMM)

RPMM3 <- rbind(ctlRPMM,RPMM2)

#clean up targets dataframe
  #select for only columns necessary for heatmap
targets <- targets  %>% select(Label.ID., location_class, Tumor_Type, WHO_GRADE, 
                             YRS_2_Recurrence, DATE_DEATH) %>%  
  #Clean up Tumor Type data
                      mutate(Tumor_Type= ifelse(Tumor_Type=="EMBRYONAL", "Embryonal",
                                                ifelse(Tumor_Type=="GLIOMA", "Glioma", 
                                                       ifelse(Tumor_Type == "Non-Tumor", "Control", Tumor_Type)))) %>%
  #create a column with mortality
                      mutate(Status = ifelse(DATE_DEATH == "Alive", "Alive", "Dead"))

#Add RPMM Cluster data to targets df
targets <- merge(targets, RPMM3, "Label.ID.")



#-------------------------------------Step 5: Heat map: unsupervised clustering-------------------------------------
loadEPICannot <- function(forHeatmap=TRUE) {
  #'@description Load EPIC annotation for heatmap
  #'@param forHeatmap If TRUE, a data.frame specific for `pheatmap` will be returned. Defaults to TRUE
  require(IlluminaHumanMethylationEPICanno.ilm10b4.hg19);
  annot.850kb3 <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19);
  if(! forHeatmap) return(annot.850kb3);
  
  row_annot <- data.frame(
    row.names = annot.850kb3$Name, 
    #Island = ifelse(annot.850kb3$Relation_to_Island=="Island", "Yes", "No"),
    Promoter = ifelse(grepl("TSS",annot.850kb3$UCSC_RefGene_Group), "Yes", "No"),
    Enhancer = ifelse(annot.850kb3$X450k_Enhancer=="TRUE" | annot.850kb3$Phantom4_Enhancers != "" | annot.850kb3$Phantom5_Enhancers != "", "Yes", "No"),
    TFBS = ifelse(grepl("chr",annot.850kb3$TFBS_NAME), "Yes", "No")
    );
  return(row_annot);
}

row_annot <- loadEPICannot();

targets$Cluster2 <- ifelse(targets$Cluster == "Low_5hmC", "Low 5hmC", "High 5hmC")

heat_annot <- data.frame(
  row.names = targets$Label.ID.,
  Tumor = targets$Tumor_Type,
  Grade = targets$WHO_GRADE,
  #Location = targets$location_class,
  RPMM = targets$Cluster2,
  Status = targets$Status
  
);

shadesOfGrey <- colorRampPalette(c("grey0", "grey100"))
shadesOfGrey(5)

#help select colors as necessary
# show_col(pal_npg("nrc")(10))
 pal_npg("nrc")(10)
# brewer.pal(n = 4, name = "Dark2")
 ann_colors <- list(
      Tumor = c(Embryonal="#E64B35FF", Ependymoma="#4DBBD5FF", 
                     Glioma = "#00A087FF", Control="#3C5488FF"),
      Grade = c('1'= "#BFBFBF", '2' ="#7F7F7F", '3' = "#3F3F3F", '4'= "black"),
    #Location = c(supratentorial= "#F39B7FFF", subtentorial="#B09C85FF"),
     Status = c(Alive="lightgray", Dead="black"),
     RPMM = c(`High 5hmC`="lightgray", `Low 5hmC`="black"),
     Promoter = c(Yes="black", No="lightgray"),
     Enhancer = c(Yes="black", No="lightgray"),
     TFBS = c(Yes="black", No="lightgray"))

#select most DHMRs
 #Load DHMR analysis
 load("./Script/9.RPMM/RPMM_TumorMosVar_Script/Tumor_MostVar042020.RData")


 DHMR5hmc <- All5hmc[row.names(All5hmc) %in% row.names(TumorTopVar), ]
 dim(DHMR5hmc)
 head(DHMR5hmc)
 
mat <- as.matrix(DHMR5hmc)

sampOrders <- targets %>% arrange(RPMMSampleOrder) %>%
  pull(Label.ID.)

HEAT_COLORS <- colorRampPalette(c("yellow","black","blue"))(2300); 

png("/Users/nasimazizgolshani/Dropbox/christensen/PCNS_Analysis/PCNS/Plots/0.052020_wPCNS0601X/6.Tumor_MostVar_052020.jpeg", res=300, units="in", height=8.27, width=11.69);
pheatmap(
  mat[ , sampOrders],
  show_rownames = FALSE,
  show_colnames = FALSE, #samples
  cluster_cols = FALSE,
  annotation_col = heat_annot,
  annotation_row = row_annot,
  color = HEAT_COLORS,
  annotation_colors = ann_colors,
  gaps_col = 27,
  border_color = NA, 
  fontsize = 13)
dev.off();

