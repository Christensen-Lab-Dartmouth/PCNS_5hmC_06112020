
#Pediatric CNS Tumors
#Visualizing LOLA results of Tumor V Non Tumor DHMR
#Nasim
#02/12/20

#~~~~~~~~~~~~~~~~~~~~~~~~~`Step 1 : Initialization/ Load Data~~~~~~~~~~~~~~~~~~~~````
rm(list=ls())
library(tidyverse)
library(LOLA)
setwd("/Users/nasimazizgolshani/Dropbox/christensen/PCNS_Analysis/PCNS")
load("./Script/16.LOLA/Files/PCNS_LOLA_Output.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~Step 2: Restructure LOLA Results~~~~~~~~~~~~~~~~~~~````
#First look at all results prior to subsetting to cells of interest!
#These will include cancer cell lines such as breast and prostate
#Exploring these since Cistrome Epi and Cistrome do not have data for CNS cell lines
# look through results by collection 
cistrome_epi <- LOLA_Sig_CpG %>% filter(collection=="cistrome_epigenome")
encode <- LOLA_Sig_CpG %>% filter(collection=="encode_tfbs")
cistrome <-  LOLA_Sig_CpG %>% filter(collection=="cistrome_cistrome")
shef <- LOLA_Sig_CpG %>% filter(collection=="sheffield_dnase")
codex <- LOLA_Sig_CpG %>% filter(collection=="codex")
ucsc <- LOLA_Sig_CpG %>% filter(collection=="ucsc_features")

#Order by meanRank
cistrome_epi <- cistrome_epi[order(cistrome_epi$meanRnk, decreasing=FALSE),]
encode <- encode[order(encode$meanRnk, decreasing=FALSE),]
cistrome <- cistrome[order(cistrome$meanRnk, decreasing=FALSE),]
codex <- codex[order(codex$meanRnk, decreasing=FALSE),]

# subset to those w/ qvalue <0.05
cistrome_epi <- cistrome_epi[cistrome_epi$qValue < 0.05,]
encode <- encode[encode$qValue < 0.05,]
cistrome <- cistrome[cistrome$qValue < 0.05,]
codex <- codex[codex$qValue < 0.05,]

#Look at Cistrome data 
ggplot(cistrome, aes(x = antibody, y = -log10(qValue), size = oddsRatio, fill = cellType))+
  geom_point(alpha = 0.5, shape = 21) +
  guides(size=guide_legend(override.aes = list(fill="black", alpha=1)))+
  labs(x = "Transcription factor binding sites (cistrome database)", y = expression("-log"[10]*"(Q-value)"),
       size = "logOdds ratio", fill = "Cell type") +
  scale_size(range = c(1, 12)) 

#Look at Cistrome Epigenetics data 
  ggplot(cistrome_epi, aes(x = antibody, y = -log10(qValue), size = oddsRatio, fill = cellType))+
    geom_point(alpha = 0.5, shape = 21) +
    guides(size=guide_legend(override.aes = list(fill="black", alpha=1)))+
    labs(x = "Cistrome Histone Modifications", y = expression("-log"[10]*"(Q-value)"),
         size = "logOdds ratio", fill = "Cell type") +
    scale_size(range = c(1, 12)) 
  
  #Look at Codex  data 
  ggplot(codex, aes(x = antibody, y = -log10(qValue), size = oddsRatio, fill = cellType))+
    geom_point(alpha = 0.5, shape = 21) +
    guides(size=guide_legend(override.aes = list(fill="black", alpha=1)))+
    labs(x = "Codex TFBS", y = expression("-log"[10]*"(Q-value)"),
         size = "logOdds ratio", fill = "Cell type") +
    scale_size(range = c(1, 12)) 
  
#~~~~~~~~~~~~~~~~~~~~~~~~~Step 3: Subset to cell lines of interest and visualize~~~~~~~~~~~~~~~~~~````
# Cell_Types_Of_Interest = c("U87", "SK-N-SH", "SK-N-MC", 
#                            "Embryonic Stem Cell", "PFSK-1", "HCPEpiC", 
#                           "NH-A", "Umbilical Cord Blood Stem and Progenitor Cells", 
#                            "HBMEC","HRPEpiC","HA-sp","BE2_C",
#                            "Gliobla","WERI-Rb-1","SK-N-SH_RA",
#                            "H1hesc", "H1-hESC")

Cell_Types_Of_Interest = c("U87", "Embryonic Stem Cell",   
                           "NH-A", "Umbilical Cord Blood Stem and Progenitor Cells", 
                           "HBMEC","HA-sp","Gliobla","WERI-Rb-1",
                           "H1hesc", "H1-hESC"
                           )
Lola_CellsofInterest <- LOLA_Sig_CpG %>% filter(cellType %in% Cell_Types_Of_Interest)
Lola_CellsofInterest_ENCODE <- Lola_CellsofInterest %>% filter(collection=="encode_tfbs",
                                                               qValue != 1) 
EncodeLOLA <- separate(Lola_CellsofInterest_ENCODE, col="antibody", into="ab2")
 EncodeLOLA$cellTypeDesc <- ifelse(EncodeLOLA$cellType=="Gliobla", "Glioblastoma",
                                          ifelse(EncodeLOLA$cellType=="HA-sp","Spinal Astrocyte",
                                                 # ifelse(EncodeLOLA$cellType=="HCPEpiC", "HCPEpiC: Human Choroid Plexus Epithelial Cells",
                                                 #        ifelse(EncodeLOLA$cellType=="HRPEpiC", "HRPEpiC: Human Retinal Pigment Cells",
                                                              ifelse(EncodeLOLA$cellType=="WERI-Rb-1", "Retinoblastoma", 
                                                                      ifelse(EncodeLOLA$cellType=="H1-hESC", "Embryonic Stem Cells",
                                                                             ifelse(EncodeLOLA$cellType=="HBMEC","Cerebral Endothelium",
                                                                                    ifelse(EncodeLOLA$cellType=="NH-A", "Astrocyte", EncodeLOLA$cellType))))))

ppi=300
png("./Script/16.LOLA/Plots/lola_encode_TFBS060820.png", width=11*ppi, height=8*ppi, res=ppi)
#fill = c("chartreuse3", "lightskyblue1", "sienna1")
ggplot(EncodeLOLA, aes(x = ab2, y = -log10(qValue), size = oddsRatio, fill = cellTypeDesc))+
  geom_point(alpha = 0.5, shape = 21) +
  guides(size=guide_legend(override.aes = list(fill="black", alpha=1)),
         fill = guide_legend(override.aes = list(size = 5)))+
  labs(x = "Transcription factor binding sites (Encode database)", y = expression("-log"[10]*"(Q-value)"),
       size = "Odds ratio", fill = "Cell type") +
  scale_size(range = c(1, 15)) +
  #scale_fill_manual(values = fill) + 
  geom_hline(yintercept = 1.3, linetype="dashed") +
  theme(legend.key.size = unit(1, "cm"),
        panel.grid.major = element_line(colour = "#d3d3d3"), 
        panel.grid.minor = element_blank(), 
       # panel.border = element_rect(fill = NA, colour = "black", size = 0.6, linetype = "solid"), 
        panel.background = element_blank(),
        axis.text.x=element_text(angle = 75, colour="black", size = 12, hjust = 1), 
        axis.text.y=element_text(colour="black", size = 12), 
        axis.title.x=element_text(colour="black", size = 12), 
        axis.title.y=element_text(colour="black", size = 12), 
        legend.key = element_blank(), 
        legend.text = element_text(colour="black", size = 12), 
        legend.title = element_text(colour="black", size = 12))


dev.off()



##### Codex Plots
Lola_CellsofInterest_Codex <- Lola_CellsofInterest %>% filter(collection=="codex",
                                                               qValue !=1) 
