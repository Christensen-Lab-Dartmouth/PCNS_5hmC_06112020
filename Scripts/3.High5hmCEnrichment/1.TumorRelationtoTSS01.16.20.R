#title: "1.Genomic_Context_Analysis Tumor Samples"
#author: "nasim"
#date: "12/13/2019"
#Aims: Examine relation of Tumor Top 5% most hydroxymethylated regions with respect to TSS and CpG Islands



setwd('/Users/nasimazizgolshani/Dropbox/christensen/PCNS_Analysis/PCNS')


# ==========================================================================================
# Initialization
# ==========================================================================================

# Load Required Packages
library(reshape)
library(ggpubr)
library(epitools) #OR calculation
library(forestplot)
library(tidyverse)
library(GenomicRanges)
library(genomation)
# Set WD for  folder


#Load 5mC data mean/median data
load("./Files/01162020_unique_hydroxy_methyl_ss.RData")

# Load Top 5% most hyperhydroxymethylated CpGs
load('./Files/01162020_top5_5hmc_unique.RData')
top5 <- top5unique
rm(top5unique)
# Load EPIC annotation file created in previous file
load('./Files/brainSpecificAnnotation.RData')
# Add annotation collapsing shores and shelves for island context
epicAnnotation$islandAdjust <- ifelse(epicAnnotation$Relation_to_Island %in% c('N_Shelf','S_Shelf'),'Shelf',
                                      ifelse(epicAnnotation$Relation_to_Island %in% c('N_Shore','S_Shore'),'Shore',epicAnnotation$Relation_to_Island))

#Generate new variables for MAPINFO. A probe lists one location, to get genomic ranges we need two locations
epicAnnotation$pos_start = epicAnnotation$pos
epicAnnotation$pos_end <- epicAnnotation$pos + 1

#Take out probes in sex chromosomes and cross hybridizing probes
epicAnnotation <- epicAnnotation[row.names(epicAnnotation) %in% hydroxy_unique_ss$ID,]

#Keep the original complete version of annotation file
epicAnnotationFull <- epicAnnotation

# Subset annotation data to CpGs included in analysis
epicAnnotation <- epicAnnotation[row.names(epicAnnotation) %in% top5$ID,]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make Grange object for good probes from EPIC array
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#column numbers may differ depending on which annotation file you use and what you've added to it
#all you need is chr/pos/position end (position start is the same as pos so optional)
annotsub <- epicAnnotation[,c(1,2, 59,60)]

#Create a 'GRanges' object from the Illumina annotation file 
top5_gr <- makeGRangesFromDataFrame(annotsub, keep.extra.columns=TRUE, ignore.strand=TRUE, seqnames.field = "chr", start.field="pos_start", end.field="pos_end")
#Print the GRange object to see what it looks like
top5_gr

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calculate distance and orientation of 5hmC relative to nearest canonical TSS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load UCSC reference genome annotation for hg38
#up.flank is how far up-stream from TSS want to detect promoter boundaries and downstream for down.flank
transcFeat <- readTranscriptFeatures('./Script/4.MedianBetasatTSS/Files/UCSC_hg19_refGene.bed',up.flank = 2000, down.flank = 2000)

# get indices of high 5hmC CpGs in annotation file that precede transcription start sites (TSSs)
tss_precede <- precede(top5_gr, transcFeat$TSSes)
# do the same for those that follow TSSs
tss_follow <- follow(top5_gr, transcFeat$TSSes)

# identify the nearest TSS for each of the high 5hmC CpGs 
nearest_tss <- nearest(top5_gr, transcFeat$TSSes)
# index for these TSSs only
transcFeat_tss_sub <- transcFeat$TSSes[nearest_tss]

# calculate distance to nearest TSS for each high 5hmC CpG
dist_to_TSS <- distance(top5_gr, transcFeat_tss_sub)
# divide by 1000 to convert to kilobases
dist_to_TSS <- dist_to_TSS/1000


# identify which of CpGs are up- or down-stream of the nearest 
upstream_indices <- which(nearest_tss==tss_precede)
downstream_indices <- which(nearest_tss!=tss_precede)
dist_to_TSS_upstream <- dist_to_TSS[upstream_indices]
dist_to_TSS_downstream <- dist_to_TSS[downstream_indices]

# visualize distributions of each
hist(dist_to_TSS_upstream, 20)
hist(dist_to_TSS_downstream, 20)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calculate and visualize proportions of high 5hmC in various distance ranges of nearest TSS 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# write function that will calculate proportions of high 5hmC CpGs within specified 'bins' for genomic range
# (assuming distance to nearest TSS, up and downstre am, have been determined)
calc_prop <- function(dist_up, dist_down, bins){
  p_up <- c()
  p_down <- c()
  for(i in 1:(length(bins)-1)){
    t1 <- table(dist_up>bins[i] & dist_up<=bins[i+1])
    if(length(t1)>1){
      p_up[i] <- t1[[2]]/length(top5_gr)*100
    } else {
      p_up[i] <- 0
    }
  }
  for(i in 1:(length(bins)-1)){
    t2 <- table(dist_down>bins[i] & dist_down<=bins[i+1])
    if(length(t2)>1){ 
      p_down[i] <- t2[[2]]/length(top5_gr)*100
    } else {
      p_down[i] <- 0
    }
  }
  t_up_end <- table(dist_up>bins[length(bins)])
  if(length(t_up_end)>1){p_up_end <- t_up_end[[2]]/length(top5_gr)*100}
  else{p_up_end <- 0}
  
  t_down_end <- table(dist_down>bins[length(bins)])
  if(length(t_down_end)>1){p_down_end <- t_down_end[[2]]/length(top5_gr)*100}
  else{p_down_end <- 0}
  
  p_combo <- c(p_up_end, rev(p_up), p_down, p_down_end)
  p_combo
}

# caclulate & plot proportions of high 5hmC CpGs within 0-5, 5-50, 50-500, or >500kb of the nearest TSS 
props <- calc_prop(dist_to_TSS_upstream, dist_to_TSS_downstream, c(0, 5, 50, 500))
names(props) <- c("< -500", "-500 to 50", "-50 to -5", "-5 to 0", "0 to 5", "5 to 50", "50 to 500", "> 500")
ppi=300
png("./Plots/4.Genomic_Context/1.Relation_to_TSS/high_5hmC_CpGs_TSS_relation_barplot.png", width=6*ppi, height=5*ppi, res=ppi)
mp <- barplot(props, col = "steelblue3", las = 1, ylim = c(0, 30),
              ylab = "Proportion of high 5hmC CpGs (%)", xlab = "", 
              cex.lab = 1.3, cex.axis = 1.15, xaxt = "n")
labels <- c("< -500", "-500 to 50", "-50 to -5", "-5 to 0", "0 to 5", "5 to 50", "50 to 500", "> 500")
text(mp, par("usr")[3], labels = labels, srt = 45, adj = c(1.1,1.1), xpd = T, cex=1.3)
mtext("Distance to TSS (kb)", side=1, line=4, cex = 1.15)
dev.off()

# caclulate & plot proportions of high 5hmC CpGs +/-5kb of nearest TSS, giving proportions for every 1kb window 
props <- calc_prop(dist_to_TSS_upstream, dist_to_TSS_downstream, seq(0, 5, 1))
# drop 1st & last element as these relate to the proportion of CpGs outside of this region 
props <- props[-c(1, length(props))]
names(props) <- c("-5 to -4", "-4 to -3", "-3 to -2", "-2 to -1", "-1 to 0", "0 to 1", "1 to 2", "2 to 3", "3 to 4", "4 to 5")
png("./Plots/4.Genomic_Context/1.Relation_to_TSS/high_5hmC_CpGs_TSS_relation_10bp.png", width=6*ppi, height=5*ppi, res=ppi)
mp <- barplot(props, col = "steelblue3", las = 1, ylim = c(0, 16),
              ylab = "Proportion of high 5hmC CpGs (%)", xlab = "", 
              cex.lab = 1.3, cex.axis = 1.15, xaxt = "n")
text(mp, par("usr")[3], labels = names(props), srt = 45, adj = c(1.1,1.1), xpd = T, cex=1.3)
mtext("Distance to TSS (kb)", side=1, line=4, cex = 1.15)
dev.off()

#X Axis Annotation labels
neglist <- -49:0
neglist2 <- -50:-1
negcompletelist <- paste(neglist2, sep=" to ", neglist)
poslist <- 0:49
poslist2 <- 1:50
poscompletelist <- paste(poslist, sep=" to ", poslist2)
label_list <- c(negcompletelist,poscompletelist)

# caclulate & plot proportions of high 5hmC CpGs +/-5kb of nearest TSS, giving proportions for every 1kb window 
props <- calc_prop(dist_to_TSS_upstream, dist_to_TSS_downstream, seq(0, 50, 1))
# drop 1st & last element as these relate to the proportion of CpGs outside of this region 
props <- props[-c(1, length(props))]
names(props) <- label_list
png("./Plots/00.052020_wPCNS0602x/11.NA_high_5hmC_CpGs_TSS_relation_50bp.png", width=15*ppi, height=5*ppi, res=ppi)
mp <- barplot(props, col = "steelblue3", las = 1, ylim = c(0, 10),
              ylab = "Proportion of high 5hmC CpGs (%)", xlab = "", 
              cex.lab = 1.3, cex.axis = 1.15, xaxt = "n")
text(mp, par("usr")[3], labels = names(props), srt = 60, adj = c(1.1,1.1), xpd = T, cex=0.6)
mtext("Distance to TSS (kb)", side=1, line=4, cex = 1.15)
dev.off()

# x <- ggplot(prop.df, aes(x=row.names(prop.df)))+
#   geom_bar()+
#     geom_text(aes(label=row.names(prop.df)), vjust=-0.3, size=3.5)+
#   theme_minimal()

#X Axis Annotation labels
neglist <- seq(-490,0,by=10)
neglist2 <- seq(-500,-10, by=10)
negcompletelist <- paste(neglist2, sep=" to ", neglist)
poslist <- seq(0,490, by=10)
poslist2 <- seq(10,500,by=10)
poscompletelist <- paste(poslist, sep=" to ", poslist2)
label_list <- c(negcompletelist,poscompletelist)

# caclulate & plot proportions of high 5hmC CpGs +/-5kb of nearest TSS, giving proportions for every 1kb window 
props <- calc_prop(dist_to_TSS_upstream, dist_to_TSS_downstream, seq(0, 500, 10))
# drop 1st & last element as these relate to the proportion of CpGs outside of this region 
props <- props[-c(1, length(props))]
names(props) <- label_list
png("./Plots/00.052020_wPCNS0602x/11.NA_high_5hmC_CpGs_TSS_relation_500bp.png", width=15*ppi, height=5*ppi, res=ppi)
mp <- barplot(props, col = "steelblue3", las = 1, ylim = c(0, 21),
              ylab = "Proportion of high 5hmC CpGs (%)", xlab = "", 
              cex.lab = 1.3, cex.axis = 1.15, xaxt = "n")
text(mp, par("usr")[3], labels = names(props), srt = 60, adj = c(1.1,1.1), xpd = T, cex=0.4)
mtext("Distance to TSS (kb)", side=1, line=4, cex = 1.15)
dev.off()

