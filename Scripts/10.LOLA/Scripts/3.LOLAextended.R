
#title: "PCNS - LOLA Generate Genomic Context Objects"
#Here I use the extended files provided by LOLA (roadmaps and jaspar motifs)
#author: "Nasim"
#date: "02/10/20"
#Step 1: Initialization
rm(list=ls())
library(genomation)
library(GenomicRanges)
library(tidyverse)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(LOLA)
library(missMethyl)

setwd("/Users/nasimazizgolshani/Dropbox/christensen/PCNS_Analysis/PCNS/")

#Step 2: Load Data

#To run LOLA you will need to ensure that you have downloaded the LOLA region Databases: 
#http://databio.org/regiondb. Here I have downloaded LOLACoreCaches_180412.tgz which was uploaded to the site 01-May-2018 20:47. 
#You will also have to ensure that you have simpleCache package installed.
#Direct link to tar files: http://big.databio.org/regiondb/
# Load Data
TumVsNon <- read.csv("./Script/15.EWAS/DMRsFinal/BetasTumVSNon_results020520.csv", 
                     stringsAsFactors = F, row.names = 1)

# Wrangle CpGs
#Load Top 5% bed file with Genomic regions of 37k CpGs w highest 5hmC
Top5 <- read.delim("./Script/15.EWAS/DMRsFinal/Top5GR020520.bed", header=F)
colnames(Top5) <- c("chr", "pos_start", "pos_end", "Name")
#Load bed file with DHMR
DHMR <- read.delim("./Script/15.EWAS/DMRsFinal/TumVNonDMR020520.bed", header=F)
colnames(DHMR) <- c("chr", "pos_start", "pos_end", "Name")
# Define Genomic Regions
# Base comparison - aka "Universe"
GR_Universe <- makeGRangesFromDataFrame(Top5, keep.extra.columns=TRUE, ignore.strand=TRUE, seqnames.field = "chr", start.field = "pos_start", end.field = "pos_end")

# Significant CpGs
GR_Sig_CpGs <-  makeGRangesFromDataFrame(DHMR, keep.extra.columns=TRUE, ignore.strand=TRUE, seqnames.field = "chr", start.field = "pos_start", end.field = "pos_end")


# Step 3: Run LOLA
## Load Database and Query Regions

regionDB = loadRegionDB(dbLocation = "./Script/16.LOLA/Files/LOLAExt/hg19")

## Run LOLA
LOLA_Sig_CpG = runLOLA(GR_Sig_CpGs, GR_Universe, regionDB)

#separate data by collection 
roadmaps_epi <- LOLA_Sig_CpG %>% filter(collection=="roadmap_epigenomics")

#Separate cell line name (E0X in filename column corresponds to a cell line)
LOLAroadmaps <- separate(roadmaps_epi, col="filename", into="CellLine", sep="-")
nrow(LOLAroadmaps)

#Subset to the following cell lines: E068 (brain anterior caudate), E069 (Brain cingulate gyrus), E072, E071, E074 E067, E073, E070, E082, E081

#All CNS related cell lines
LOLAroadmaps <- LOLAroadmaps %>% filter(CellLine %in% c("E067", "E068","E069","E070", "E071","E072", "E073",
                                                      "E074", "E081", "E082", "E007", "E009", "E010"))

#Remove NAs
LOLAroadmaps <- LOLAroadmaps[!is.na(LOLAroadmaps$antibody), ]
#Designate sites as activating or repressive
#Some from textbks/wikipedia others from publications
#https://www.nature.com/articles/s41586-019-1534-3
LOLAroadmaps2 <- LOLAroadmaps %>% filter(CellLine %in% c("E069","E081", "E082",  "E009")) %>%
  mutate(CellLineDesc = ifelse(CellLine=="E069", "Adult Anterior Cingulate Gyrus",
                               ifelse(CellLine=="E081", "Male Fetal Brain",
                                      ifelse(CellLine=="E082", "Female Fetal Brain",
                                             ifelse(CellLine=="E009", "Neuronal Progenitor Cells", NA)))))


ppi=300
png("./Script/16.LOLA/Plots/RoadmapsRepressActivating.png", width=11*ppi, height=8*ppi, res=ppi)
ggplot(LOLAroadmaps2, aes(x = antibody, y = -log10(qValue), size = oddsRatio, fill = CellLineDesc),
       group_by(Effect))+
  geom_point(alpha = 0.5, shape = 21) +
  guides(size=guide_legend(override.aes = list(fill="black", alpha=1)),
         fill = guide_legend(override.aes = list(size = 5)))+
  labs(x = "", y = expression("-log"[10]*"(Q-value)"),
       size = "Odds ratio", fill = "Cell type") +
  scale_size(range = c(1, 12)) +  
  scale_x_discrete(limits=c("H3K4me3","H3K36me3","H3K9me3","H3K27me3","H3K4me1","H3K9ac","H3K27ac", "H2A.Z"))+
  geom_hline(yintercept = 2, linetype="dashed") +
  theme(legend.key.size = unit(1, "cm"),
        panel.grid.major = element_line(colour = "#d3d3d3"), 
        panel.grid.minor = element_blank(), 
       # panel.border = element_rect(fill = NA, colour = "black", size = 0.6, linetype = "solid"), 
        panel.background = element_blank(),
        axis.text.x=element_text(angle = 75, colour="black", size = 12, hjust = 1), 
        axis.text.y=element_text(colour="black", size = 12), 
        #axis.title.x=element_text(colour="black", size = 12, margin = margin(t = 30, r = 0, b = 0, l = 0)), 
        axis.title.y=element_text(colour="black", size = 12), 
        legend.key = element_blank(), 
        legend.text = element_text(colour="black", size = 12), 
        legend.title = element_text(colour="black", size = 12))
dev.off()

#separate data by collection and q value <0.0000001
roadmaps_epi <- LOLA_Sig_CpG %>% filter(collection=="roadmap_epigenomics",
                                        qValue <0.0000001)

#Repeat without subsetting to brain cell lines

