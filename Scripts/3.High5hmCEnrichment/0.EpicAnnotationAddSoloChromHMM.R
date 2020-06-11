# 
#   title: "EPIC_annotation_file"
# author: "nasim"
# date: "12/2/2019"

  ############################################################################################
# PCNS Tumors 
# Genomic Context Analysis - Build Annotation File for EPIC Array
# Code source MM Field Cancerization 

############################################################################################


# ==========================================================================================
# Initialization
# ==========================================================================================

# Load required packages

library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) # for Illumina annotation data
library(data.table) # for reading in results files
library(GenomicRanges) # for generating genomic ranges objects
library(genomation) # for assigning promorter, intron, exon annotation
library(tidyverse)
library(rtracklayer)

  
  # ==========================================================================================
  # Generate Annotation File
  # ==========================================================================================
  # Load illumina annotation file
  epicAnnotation <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  # Generate GRanges object for annotation file
  annotSub <- epicAnnotation[,c('Name','chr','pos')]
  epic.gr <- makeGRangesFromDataFrame(annotSub,keep.extra.columns = T,ignore.strand = T,seqnames.field = 'chr',start.field = 'pos',end.field = 'pos')
  rm(annotSub)
  rm(epicAnnotation)
  
setwd("/Users/nasimazizgolshani/Dropbox/christensen/PCNS_Analysis/PCNS/Script/3.top5_5hmc_enrichment")
load("brainSpecificAnnotation.RData")


# ==========================================================================================
# Add Solo-WCGW Loci Annotation Data
# ==========================================================================================
#Can Downlaod files here: https://zwdzwd.github.io/pmd
#Choose Illumina Epic Files

# Load annotation file for probes in highly methylated domains (HMD)
solo_hmd_probes <- read.table(file = './Script/3.top5_5hmc_enrichment/Downloaded_Annot_Files/EPIC.comHMD.probes.tsv', sep = '\t', header = F)
# Annotate EPIC annotation file
epicAnnotation$soloWCGW_HMD <- ifelse(row.names(epicAnnotation) %in% solo_hmd_probes[,4],1,0)

# Load annotation file for probes in partially methylated domains (PMD)
solo_pmd_probes <- read.table(file = './Script/3.top5_5hmc_enrichment/Downloaded_Annot_Files/EPIC.comPMD.probes.tsv', sep = '\t', header = F)
# Annotate EPIC annotation file
epicAnnotation$soloWCGW_PMD <- ifelse(row.names(epicAnnotation) %in% solo_pmd_probes[,4],1,0)

# Clean workspace
rm(solo_hmd_probes,solo_pmd_probes)

# ==========================================================================================
# Add Annotation for Marks Available in Roadmap ChIP-seq Data 
# Data for E073 Dorsolateral Prefrontal Cortex
# Adapted from Code by Meghan Muse
# ==========================================================================================

#Note on consolidated ChIP-seq data:
#Signal processing engine of MACSv2.0.10 peak caller to generate genome-wide signal coverage tracks
#Negative log10 of the Poisson p-value of ChIP-seq or DNase counts relative to expected background counts
#These signal confidence scores provides a measure of statistical significan of the observed enrichment.
# get ChIP-seq data for available chromatin marks: save paths to data online
#These files are very heavy and take time to download 
#I downloaded the files in parallel on polaris saved them and then in a separate time added them to the annotation file

E073_paths <- list("https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E073-H3K4me1.pval.signal.bigwig", 
                    "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E073-H3K4me3.pval.signal.bigwig",
                    "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E073-H3K9ac.pval.signal.bigwig",
                    "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E073-H3K9me3.pval.signal.bigwig",
                    "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E073-H3K27ac.pval.signal.bigwig",
                    "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E073-H3K27me3.pval.signal.bigwig",
                    "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E073-H3K36me3.pval.signal.bigwig")
 
 
# # import BIGWIG files from consolidated epigenomes from ROADMAP website 
# E073_data <- lapply(E073_paths, import.bw)
# # save data so we don't have to download again 
# save(E073_data, file = "roadmap_prefrontalCTX_ChIP_seq.Rdata")
load("roadmap_prefrontalCTX_ChIP_seq.Rdata")


# Function for generating binary annotation row from Granges object (from code by Owen Wilkins)
add_annot <- function(mark_object,epic.gr,epicAnnotation){
  # index for sites w/ -log10P > 2 as roadmap advise this provides good signal/noise separation 
  mark_bs <- mark_object[score(mark_object)>=2, ]
  # find the overlapping regions between the high confidence peaks and the CpGs on the 450k array 
  overlaps <- findOverlaps(mark_bs, epic.gr)
  # get the indicies of the overlapping sites in the 450K annotation file 
  indicies <- subjectHits(overlaps)
  # add dummy variable to annotation file for CpGs overlapping with this mark 
  mark_key <- rep(0, nrow(epicAnnotation))
  mark_key[indicies] <- 1
  mark_key
}

# Generate binary annotations for EPIC annotation file
epicAnnotation$Prefrontal_Ctx_H3K4me1 <- add_annot(E073_data[[1]],epic.gr,epicAnnotation)
epicAnnotation$Prefrontal_Ctx_H3K4me3 <- add_annot(E073_data[[2]],epic.gr,epicAnnotation)
epicAnnotation$Prefrontal_Ctx_H3K9ac <- add_annot(E073_data[[3]],epic.gr,epicAnnotation)
epicAnnotation$Prefrontal_Ctx_H3K9me3 <- add_annot(E073_data[[4]],epic.gr,epicAnnotation)
epicAnnotation$Prefrontal_Ctx_H3K27ac <- add_annot(E073_data[[5]],epic.gr,epicAnnotation)
epicAnnotation$Prefrontal_Ctx_H3K27me3 <- add_annot(E073_data[[6]],epic.gr,epicAnnotation)
epicAnnotation$Prefrontal_Ctx_H3K36me3 <- add_annot(E073_data[[7]],epic.gr,epicAnnotation)

# Double check annotation by looking at tables
#apply(epicAnnotation[,paste('Prefrontal_Ctx_',names(E073_data),sep = '')],2,table)

# Clean workspace
rm(E073_data)

#Repeat for Inferior Temporal Lobe
E072_paths <- list("https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E072-H3K4me1.pval.signal.bigwig", 
                   "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E072-H3K4me3.pval.signal.bigwig",
                   "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E072-H3K9ac.pval.signal.bigwig",
                   "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E072-H3K9me3.pval.signal.bigwig",
                   "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E072-H3K27ac.pval.signal.bigwig",
                   "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E072-H3K27me3.pval.signal.bigwig",
                   "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E072-H3K36me3.pval.signal.bigwig")


# import BIGWIG files from consolidated epigenomes from ROADMAP website 
E072_data <- lapply(E072_paths, import.bw)
# save data so we don't have to download again 
save(E072_data, file = "roadmap_TempLobe_ChIP_seq.Rdata")
load("roadmap_TempLobe_ChIP_seq.Rdata")

# Generate binary annotations for EPIC annotation file
epicAnnotation$TempLobe_H3K4me1 <- add_annot(E072_data[[1]],epic.gr,epicAnnotation)
epicAnnotation$TempLobe_H3K4me3 <- add_annot(E072_data[[2]],epic.gr,epicAnnotation)
epicAnnotation$TempLobe_H3K9ac <- add_annot(E072_data[[3]],epic.gr,epicAnnotation)
epicAnnotation$TempLobe_H3K9me3 <- add_annot(E072_data[[4]],epic.gr,epicAnnotation)
epicAnnotation$TempLobe_H3K27ac <- add_annot(E072_data[[5]],epic.gr,epicAnnotation)
epicAnnotation$TempLobe_H3K27me3 <- add_annot(E072_data[[6]],epic.gr,epicAnnotation)
epicAnnotation$TempLobe_H3K36me3 <- add_annot(E072_data[[7]],epic.gr,epicAnnotation)

# Double check annotation by looking at tables
#apply(epicAnnotation[,paste('tempL_',names(E072_data),sep = '')],2,table)

# Clean workspace
rm(E072_data)


#Fetal Female Brain
E082_paths <- list("https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E082-H3K4me1.pval.signal.bigwig", 
                   "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E082-H3K4me3.pval.signal.bigwig",
                   "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E082-H3K9me3.pval.signal.bigwig",
                   "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E082-H3K27me3.pval.signal.bigwig",
                   "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E082-H3K36me3.pval.signal.bigwig")


# import BIGWIG files from consolidated epigenomes from ROADMAP website 
E082_data <- lapply(E082_paths, import.bw)
# save data so we don't have to download again 
save(E082_data, file = "roadmap_FetalFB_ChIP_seq.Rdata")
load("roadmap_FetalFB_ChIP_seq.Rdata")

# Generate binary annotations for EPIC annotation file
epicAnnotation$FetalFB_H3K4me1 <- add_annot(E082_data[[1]],epic.gr,epicAnnotation)
epicAnnotation$FetalFB_H3K4me3 <- add_annot(E082_data[[2]],epic.gr,epicAnnotation)
epicAnnotation$FetalFB_H3K9me3 <- add_annot(E082_data[[3]],epic.gr,epicAnnotation)
epicAnnotation$FetalFB_H3K27me3 <- add_annot(E082_data[[4]],epic.gr,epicAnnotation)
epicAnnotation$FetalFB_H3K36me3 <- add_annot(E082_data[[5]],epic.gr,epicAnnotation)

# Double check annotation by looking at tables
#apply(epicAnnotation[,paste('FetalFB_',names(E082_data),sep = '')],2,table)

# Clean workspace
rm(E082_data)

#Fetal Male Brain
E081_paths <- list("https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E081-H3K4me1.pval.signal.bigwig", 
                   "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E081-H3K4me3.pval.signal.bigwig",
                   "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E081-H3K9me3.pval.signal.bigwig",
                   "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E081-H3K27me3.pval.signal.bigwig",
                   "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E081-H3K36me3.pval.signal.bigwig")


# import BIGWIG files from consolidated epigenomes from ROADMAP website 
E081_data <- lapply(E081_paths, import.bw)
# save data so we don't have to download again 
save(E081_data, file = "roadmap_FetalMB_ChIP_seq.Rdata")
load("roadmap_FetalMB_ChIP_seq.Rdata")

# Generate binary annotations for EPIC annotation file
epicAnnotation$FetalMB_H3K4me1 <- add_annot(E081_data[[1]],epic.gr,epicAnnotation)
epicAnnotation$FetalMB_H3K4me3 <- add_annot(E081_data[[2]],epic.gr,epicAnnotation)
epicAnnotation$FetalMB_H3K9me3 <- add_annot(E081_data[[3]],epic.gr,epicAnnotation)
epicAnnotation$FetalMB_H3K27me3 <- add_annot(E081_data[[4]],epic.gr,epicAnnotation)
epicAnnotation$FetalMB_H3K36me3 <- add_annot(E081_data[[5]],epic.gr,epicAnnotation)

# Double check annotation by looking at tables
#apply(epicAnnotation[,paste('FetalMB_',names(E081_data),sep = '')],2,table)

# Clean workspace
rm(E081_data)


# ==========================================================================================
# Save Updated Annotation File
# ==========================================================================================
# Check that appropriate annotation values are included
colnames(epicAnnotation)
# Save file
save(epicAnnotation,file = '/dartfs-hpc/rc/home/0/f002tf0/PCNS/Script/Annotation/brainSpecificAnnotationComplete.RData')

#Reload saved file 
load('./Script/3.top5_5hmc_enrichment/Downloaded_Annot_Files/brainSpecificAnnotationComplete.RData')

#Add other annotation including UTR regions

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


# Add annotation collapsing shores and shelves for island context
epicAnnotation$islandAdjust <- ifelse(epicAnnotation$Relation_to_Island %in% c('N_Shelf','S_Shelf'),'Shelf',
                                      ifelse(epicAnnotation$Relation_to_Island %in% c('N_Shore','S_Shore'),'Shore',epicAnnotation$Relation_to_Island))

#Limit CpGs to those we already selected for dropping cross hybridizing and sex probes
#do this by matching another list of our CpGs created earlier
#Load 5mC data mean/median data
load("./Files/12062019_unique_hydroxy_methyl_ss.RData")
#Take out probes in sex chromosomes and cross hybridizing probes
#shortcut to match CpGs in my df of median 5hmc CpG sites since I had previously dropped those probes
epicAnnotation <- epicAnnotation[row.names(epicAnnotation) %in% hydroxy_unique_ss$ID,]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Clean up Illumina Annotation of Gene Names
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract gene name annotation
IlluminaGenes  <- epicAnnotation[row.names(epicAnnotation),'UCSC_RefGene_Name']
# Divide results table into with and without gene annotation
Illumina_withGene <- epicAnnotation[IlluminaGenes != '',]
Illumina_withoutGene <- epicAnnotation[IlluminaGenes == '',]

# Remove CpGs not annotated to a gene
IlluminaGenes <-  IlluminaGenes[IlluminaGenes != '']
length(IlluminaGenes)

# Take only first listed gene name for each CpG annotation
TempGeneList <- lapply(IlluminaGenes,FUN = function(x) {unlist(strsplit(x,';'))[1]})
# Convert to list
IlluminaGenes <- unlist(TempGeneList)
length(TempGeneList)

# Add gene annotation column to results table
Illumina_withGene$gene <- TempGeneList
Illumina_withoutGene$gene <- NA

epicAnnotation$ID <- row.names(epicAnnotation)
# Merge results tables
annot <- as.data.frame(rbind(Illumina_withGene,Illumina_withoutGene))
epicAnnotation$gene <- annot[epicAnnotation$ID, "gene"]

#remove unnecessary columns
epicAnnotation <- epicAnnotation[,-c(95)]

# Save file
save(epicAnnotation,file = './Script/3.top5_5hmc_enrichment/Downloaded_Annot_Files/brainSpecificAnnotationComplete12202019.RData')


