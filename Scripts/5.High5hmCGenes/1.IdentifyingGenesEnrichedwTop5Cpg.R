
# ==========================================================================================
# Initialization
# ==========================================================================================
#Identify the genes most enriched with CpGs in my top 5% 

# Load Required Packages
library(reshape)
library(ggpubr)
library(tidyverse)
library(GenomicRanges)
library(genomation)
library(RColorBrewer)
setwd("/Users/nasimazizgolshani/Dropbox/christensen/PCNS_Analysis/PCNS")

# Load Top 5% most hyperhydroxymethylated CpGs
load('./Files/01162020_top5_5hmc_unique.RData')
top5 <- top5unique
rm(top5unique)

#Load 5mC data mean/median data
load("./Files/01162020_unique_hydroxy_methyl_ss.RData")
#merge 5mc and 5hmc data
methylboth <- merge(hydroxy_unique_ss,methyl_unique_ss, by="ID")
#isolate data to top5 %5hmc CpGs
methylboth <- methylboth[methylboth$ID %in% top5$ID,]
# Load EPIC annotation file created in previous file
#load('./Script/3.top5_5hmc_enrichment/Files/brainSpecificAnnotation.RData')
load("./Files/brainSpecificAnnotation.RData")

# Add results to annotation file
epicAnnotation$tumor_hyper <- ifelse(row.names(epicAnnotation) %in% top5$ID,1,0)

#sanity check, make sure sum of tumor_hyper column equals no. rows in top5 (39401)
sum(epicAnnotation$tumor_hyper)

#subset epicAnnotation CpGs to those we have previously selected (remove sex/snp probes)
epicAnnotation <- epicAnnotation[row.names(epicAnnotation) %in% hydroxy_unique_ss$ID,]
# ==========================================================================================
# Extract Gene Names for Most Hydroxymethylated CpGs
# ==========================================================================================

# Add chromosomal location
methylboth$chr <-  epicAnnotation[methylboth$ID,'chr']
methylboth$pos <-  epicAnnotation[methylboth$ID,'pos']
# Add island context annotation
methylboth$Island <- epicAnnotation[methylboth$ID,'Relation_to_Island']

# Extract gene name annotation
methylboth_genes <- epicAnnotation[methylboth$ID,'UCSC_RefGene_Name']
# Divide results table into with and without gene annotation
methylboth_withGene <- methylboth[methylboth_genes != '',]
methylboth_withoutGene <- methylboth[methylboth_genes == '',]

# Remove CpGs not annotated to a gene
methylboth_genes <-  methylboth_genes[methylboth_genes != '']
length(methylboth_genes)

# Take only first listed gene name for each CpG annotation
TempGeneList <- lapply(methylboth_genes,FUN = function(x) {unlist(strsplit(x,';'))[1]})
# Convert to list
methylboth_genes <- unlist(TempGeneList)
length(TempGeneList)

# Add gene annotation column to results table
methylboth_withGene$gene <- TempGeneList
methylboth_withoutGene$gene <- NA

# Remove duplicate genes
methylboth_genes <- unique(methylboth_genes)
length(methylboth_genes)
# 28348 highly hydroxymethylated CpGs annotate to 9601 unique genes

# Merge results tables
methylboth <- as.data.frame(rbind(methylboth_withGene,methylboth_withoutGene))
#Clean workspace
rm(hydroxy_unique_ss,methyl_unique_ss)

methylboth <- methylboth %>% arrange(-Median_5hmC)

#Calculte unmethylated cytosine in top 5% 
methylboth <-  mutate(methylboth, Median_C = 1 - (Median_5hmC + Median_5mC))
#sanity check, make sure the sum of the cytosines equal 1
methylboth <-  mutate(methylboth, Sum_C = Median_C + Median_5hmC + Median_5mC)
sum(methylboth$Sum_C) #37173 which means every row is equal to 1
#remove Sum_C
methylboth <- methylboth[,-15]

#Just take with Gene
#Calculte unmethylated cytosine in top 5% 
methylboth_withGene <-  mutate(methylboth_withGene, Median_C = 1 - (Median_5hmC + Median_5mC))
#sanity check, make sure the sum of the cytosines equal 1
methylboth_withGene <-  mutate(methylboth_withGene, Sum_C = Median_C + Median_5hmC + Median_5mC)
sum(methylboth_withGene$Sum_C) #39401 which means every row is equal to 1
#remove Sum_C
methylboth_withGene <- methylboth_withGene[,-14]

#Get a count for each gene
CountGenes <- table(unlist(methylboth_withGene$gene))
CountGenes <- as.data.frame(CountGenes)
#Arrange in descending order the number of times a gene comes up
#this is for genes that have multiple CpGs looking for which genes are enriched for high 5hmC CpG sites
CountGenes <- CountGenes %>% arrange(-Freq)
#write.csv(CountGenes, file="./Script/5.GenesEnrichedwTop5CpGs/Files/CountGenes.csv")

#Play around with finding proportion of CpGs of each gene in annotation file that are in my top 5%
# Extract gene name annotation
All_genes <- epicAnnotation[,'UCSC_RefGene_Name']
# Divide results table into with and without gene annotation
All_genenames <- epicAnnotation[All_genes != '',]
All_nogenes <- epicAnnotation[All_genes == '',]
# Remove CpGs not annotated to a gene
All_genes <-  All_genes[All_genes != '']
length(All_genes)

# Take only first listed gene name for each CpG annotation
TempGeneList <- lapply(All_genes,FUN = function(x) {unlist(strsplit(x,';'))[1]})
# Convert to list
All_genes <- unlist(TempGeneList)
length(TempGeneList)
# Add gene annotation column to results table
All_genenames$gene <- TempGeneList
All_genenames$gene <- as.character(All_genenames$gene)
All_nogenes$gene <- NA
All_nogenes$gene <- as.character(All_nogenes$gene)

#Get Gene Counts for all the sites
geneAnnotation <- as.data.frame(rbind(All_genenames,All_nogenes))
CountGenesAll <- table(unlist(geneAnnotation$gene))
CountGenesAll <- as.data.frame(CountGenesAll)
colnames(CountGenesAll) <- c("Gene", "Freq_All_Probes")
colnames(CountGenes) <- c("Gene", "Freq_Top5_Probes")

#Subset Counts for All the Genes in the annotation file for only those in my top 5%
AllGenesSubset <- CountGenesAll[CountGenesAll$Gene %in% CountGenes$Gene,]
#merge the two
BothGeneCounts <- merge(AllGenesSubset, CountGenes, by="Gene")
#Create a column with proportion of Genes
BothGeneCounts <- BothGeneCounts  %>% mutate(Proportion = Freq_Top5_Probes/Freq_All_Probes)

#Because genes with only a couple of occurences in the entire annotation can skew the proportions, 
#create a subset of genes with occurences >10 in our top 5% list
#Then order the list based on the proportion of all available CpGs for each gene that are in the Top 5% list
# A way to measure for enrichment
EnrichedGenes <-BothGeneCounts %>% filter(Freq_Top5_Probes >= 10)
EnrichedGenes <- EnrichedGenes %>% arrange(Proportion)
#write.csv(EnrichedGenes, file= "./Script/5.GenesEnrichedwTop5CpGs/Files/EnrichedGenesinTop5hydroxy.csv")
write.csv(EnrichedGenes, file= "./Plots/00.052020_wPCNS0602x/10.EnrichedGenesinTop5hydroxy.csv")


#save(methylboth_withGene,BothGeneCounts, file = "./Script/5.GenesEnrichedwTop5CpGs/Files/CountGenes.RData")
save(geneAnnotation,file = "./Script/5.GenesEnrichedwTop5CpGs/Files/epicAnnotationeGene.RData" )


