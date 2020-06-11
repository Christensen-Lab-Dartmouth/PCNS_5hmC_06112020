##Find Mean 5(h)mc% at distance from TSS
##Nasim Azizgolshani
##12/11/2019
#Examine porportion  level of cytosine modifications (5mC and 5hmC) of tumor vs non-tumor 
#in relation to distance from TSS 
#this is in our subset of top 5% 5hmC
# ==========================================================================================
# Step 1: Initialization
# ==========================================================================================
rm(list = ls())
# Load Required Packages
library(tidyverse)
library(GenomicRanges)
library(genomation)
library(RColorBrewer)
library(wesanderson)
library(ggsci)

# Set WD for  folder
setwd('/Users/nasimazizgolshani/Dropbox/christensen/PCNS_Analysis/PCNS')

# ==========================================================================================
# Step 2: Load Necessary Files
# ==========================================================================================

#Load all Median values for both 5mC and 5hmc for all tumors for each CpG site
#these median values are from unique tumors only 
load("./Files/01162020_unique_hydroxy_methyl_ss.RData")

# Load EPIC annotation file created in section 3
load("./Files/brainSpecificAnnotation.RData")

#Take out probes in epic annotation with sex chromosomes and cross hybridizing probes 
#by matching CpGs in my medians file
epicAnnotation <- epicAnnotation[row.names(epicAnnotation) %in% hydroxy_unique_ss$ID,] #743461

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 3: make Grange object for good probes from EPIC array
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Generate new variables for MAPINFO. A probe lists one location, to get genomic ranges we need two locations
epicAnnotation$pos_start = epicAnnotation$pos
epicAnnotation$pos_end <- epicAnnotation$pos + 1
#subset annotation to have only chr, pos, pos start and pos end
annotsub <- epicAnnotation[,c(1,2, 58, 59)]

#Create a 'GRanges' object from the Illumina annotation file 
annot_gr <- makeGRangesFromDataFrame(annotsub, keep.extra.columns=TRUE, ignore.strand=TRUE, seqnames.field = "chr", start.field="pos_start", end.field="pos_end")
#Print the GRange object to see what it looks like
annot_gr

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 4: calculate distance and orientation of 5hmC relative to nearest canonical TSS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load UCSC reference genome annotation for hg19
#up.flank is how far up-stream from TSS want to detect promoter boundaries and downstream for down.flank
transcFeat <- readTranscriptFeatures('./Script/4.MedianBetasatTSS/Files/UCSC_hg19_refGene.bed',up.flank = 2000, down.flank = 2000)

# get indices of CpGs in annotation file that precede transcription start sites (TSSs)
tss_precede <- precede(annot_gr, transcFeat$TSSes)
# do the same for those that follow TSSs
tss_follow <- follow(annot_gr, transcFeat$TSSes)

# identify the nearest TSS for each of the  CpGs 
nearest_tss <- nearest(annot_gr, transcFeat$TSSes)

# index for these TSSs only
transcFeat_tss_sub <- transcFeat$TSSes[nearest_tss]

# calculate distance to nearest TSS for each 5hmC CpG
dist_to_TSS <- distance(annot_gr, transcFeat_tss_sub)
# divide by 1000 to convert to kilobases
#dist_to_TSS <- dist_to_TSS/1000
dist_to_TSS <- round(dist_to_TSS, -1)
dist_to_TSS2 <-as.data.frame(dist_to_TSS)

#Add indices and distance to TSS to annotation file
annot_df <- as.data.frame(annot_gr)
annot_df$dist_to_TSS <- dist_to_TSS
annot_df$TSS_indices <- row.names(dist_to_TSS2)
#Add column describing direction of TSS and conditionally make upstream distances negative
annot_df$direction <- ifelse(nearest_tss==tss_precede,"upstream", 
                             ifelse(nearest_tss!=tss_precede, "downstream", NA))
annot_df$dist_w_dir <- ifelse(annot_df$direction=="upstream", paste(-annot_df$dist_to_TSS),annot_df$dist_to_TSS)
rm(dist_to_TSS2)
#save(annot_df,file="./Script/3.Top5_Enrichment/Downloaded_Annot_Files/dist_to_TSS.RData")
#combine summary stats for 5mc and 5hmc (2 separate data frames into one)
ss_unique_both <- merge(hydroxy_unique_ss, methyl_unique_ss, by ="ID")
rm(hydroxy_unique_ss)
rm(methyl_unique_ss)
#subset for only medians for simplicity
median_5mc_5hmc <- ss_unique_both[,c(1,2,5)]

#Must do the same now for control samples 
load("./Files/01162020_ctl_ss.RData")
median_ctl <- ctl_ss[,c(1,2,5)]

#merge ctl and tumors
all_medians <- merge(median_ctl,median_5mc_5hmc, by="ID")
colnames(all_medians) <- c("ID", "Ctl_5hmc", "Ctl_5mc", "Tumor_5hmc","Tumor_5mc")

#subset annot_df to only distance and CpG ID
annot_df_sub <- annot_df[,c(9,10)]
annot_df_sub <- rownames_to_column(annot_df_sub,"ID")

#combine annotation dataframe with TSS distance with data frame with all sample
matched_dist_medians <- merge(annot_df_sub, all_medians, "ID")
matched_dist_medians$dist_w_dir <- as.integer(matched_dist_medians$dist_w_dir)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 5: Calculate Means of Tumor vs Control  and Plot 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#First filter to less than 2000bp away from TSS and get medians of median beta value 
#for each distance from TSS
TSS_2kb <- matched_dist_medians %>% filter(dist_w_dir <=2000) %>% 
  filter(dist_w_dir >= (-2000)) %>%
  group_by(dist_w_dir) %>%
  summarize(M2_5hmc_Ctl= median(Ctl_5hmc), 
            M2_5mc_Ctl = median(Ctl_5mc),
            M2_5hmc_Tumor = median(Tumor_5hmc),
            M2_5mc_Tumor = median(Tumor_5mc))

#Collapse data frame 
G_DF = TSS_2kb %>%
  gather(key = "Measure", value = "Beta_value", -dist_w_dir) %>%
  #Rename column names so that legend will be cleaner
  mutate(Measure = ifelse(Measure=="M2_5hmc_Ctl", "5hmC Control", 
                          ifelse(Measure=="M2_5mc_Ctl", "5mC Control",
                                 ifelse(Measure=="M2_5hmc_Tumor", "All Tumors 5hmC",
                                        ifelse(Measure=="M2_5mc_Tumor", "All Tumors 5mC", NA ))) ))
p1 <- ggplot(G_DF, aes(x = dist_w_dir, y = Beta_value, colour = Measure)) +
  geom_line() +
  xlab("Distance from TSS (bp)") +
  ylab("Median of Median 5(h)mc") +
  theme_classic()+
  #scale_color_manual(values = wes_palette("Zissou1", n=4))
  scale_color_npg()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Repeat just with resolution at 1kb
matched_dist_medians <- matched_dist_medians %>% 
  mutate(dist_w_dir_kb = dist_w_dir / 1000) %>%
  mutate(dist_w_dir_kb = round(dist_w_dir_kb,1))


#First filter to less than 2000bp away from TSS
TSS_2kb <- matched_dist_medians %>% filter(dist_w_dir_kb <=2.0) %>% 
  filter(dist_w_dir_kb >= (-2.0)) %>%
  group_by(dist_w_dir_kb) %>%
  summarize(M2_5hmc_Ctl= median(Ctl_5hmc), 
            M2_5mc_Ctl = median(Ctl_5mc),
            M2_5hmc_Tumor = median(Tumor_5hmc),
            M2_5mc_Tumor = median(Tumor_5mc))

G_DF = TSS_2kb %>%
  gather(key = "Measure", value = "Beta_value", -dist_w_dir_kb) %>%
  #Rename column names so that legend will be cleaner
  mutate(Measure = ifelse(Measure=="M2_5hmc_Ctl", "5hmC Control", 
                          ifelse(Measure=="M2_5mc_Ctl", "5mC Control",
                                 ifelse(Measure=="M2_5hmc_Tumor", "5hmC Tumor",
                                        ifelse(Measure=="M2_5mc_Tumor", "5mC Tumor", NA ))) )) %>%
  mutate(Measure = factor(Measure, levels = c("5hmC Tumor", "5hmC Control", "5mC Tumor", "5mC Control")))
  

#Set colors for lines
library("scales")
show_col(pal_npg("nrc")(10))

ann_colors <- list(c("#E64B35FF","#4DBBD5FF", "#F39B7FFF", "#3C5488FF"))


ggplot(G_DF, aes(x = dist_w_dir_kb, y = Beta_value, colour = Measure)) +
  geom_line() +
  xlab("Distance from TSS (kilo-bp)") +
  ylab("Median of Median 5(h)mc") +
  theme_classic()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.title=element_text(size=12),
        legend.text=element_text(size=11))+
  scale_x_continuous(limits=c(-2, 2),breaks=c(-2.0, -1.5, -1.0, -0.5, 0 , 0.5, 1.0, 1.5, 2.0))+
  scale_color_manual(values = c("#E64B35FF","#4DBBD5FF", "#F39B7FFF", "#3C5488FF"))
 # geom_vline(xintercept= 0.8)+
  #scale_color_npg()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Intermediate Range
#First filter to less than 10,000bp away from TSS
TSS_2kb <- matched_dist_medians %>% filter(dist_w_dir_kb <=10.0) %>% 
  filter(dist_w_dir_kb >= (-10.0)) %>%
  group_by(dist_w_dir_kb) %>%
  summarize(M2_5hmc_Ctl= median(Ctl_5hmc), 
            M2_5mc_Ctl = median(Ctl_5mc),
            M2_5hmc_Tumor = median(Tumor_5hmc),
            M2_5mc_Tumor = median(Tumor_5mc))

G_DF = TSS_2kb %>%
  gather(key = "Measure", value = "Beta_value", -dist_w_dir_kb) %>%
  #Rename column names so that legend will be cleaner
  mutate(Measure = ifelse(Measure=="M2_5hmc_Ctl", "5hmC Control", 
                          ifelse(Measure=="M2_5mc_Ctl", "5mC Control",
                                 ifelse(Measure=="M2_5hmc_Tumor", "5hmC Tumor",
                                        ifelse(Measure=="M2_5mc_Tumor", "5mC Tumor", NA ))) )) %>%
  mutate(Measure = factor(Measure, levels = c("5mC Tumor", "5mC Control", "5hmC Control", "5hmC Tumor")))


ggplot(G_DF, aes(x = dist_w_dir_kb, y = Beta_value, colour = Measure)) +
  geom_line(size =1.2) +
  xlab("Distance from TSS (kilo-bp)") +
  ylab("Median 5(h)mc") +
  theme_classic()+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        legend.title=element_text(size=20),
        legend.text=element_text(size=18))+
  scale_x_continuous(limits=c(-10, 10), breaks= seq(-10.0, 10.0, by =1.0))+
  scale_color_manual(values = c("#E64B35FF","#4DBBD5FF", "#F39B7FFF", "#3C5488FF"))

  #scale_color_manual(values = wes_palette("Zissou1", n=4))
  # geom_vline(xintercept= 0.8)+
  #scale_color_npg()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Longer Range
matched_dist_medians$round_dist <- round(matched_dist_medians$dist_w_dir_kb, -1)
TSS_long_range <- matched_dist_medians %>% group_by(round_dist) %>%
  filter(round_dist < 500)%>%
  filter(round_dist > (-500)) %>%
  summarize(M2_5hmc_Ctl= median(Ctl_5hmc), 
            M2_5mc_Ctl = median(Ctl_5mc),
            M2_5hmc_Tumor = median(Tumor_5hmc),
            M2_5mc_Tumor = median(Tumor_5mc))

G_DF = TSS_long_range %>%
  gather(key = "Measure", value = "Beta_value", -round_dist) %>%
  #Rename column names so that legend will be cleaner
  mutate(Measure = ifelse(Measure=="M2_5hmc_Ctl", "5hmC Control", 
                          ifelse(Measure=="M2_5mc_Ctl", "5mC Control",
                                 ifelse(Measure=="M2_5hmc_Tumor", "5hmC Tumor",
                                        ifelse(Measure=="M2_5mc_Tumor", "5mC Tumor", NA ))) )) %>%
  mutate(Measure = factor(Measure, levels = c("5hmC Tumor", "5hmC Control", "5mC Tumor", "5mC Control")))



#plot
ggplot(G_DF, aes(x = round_dist, y = Beta_value, colour = Measure)) +
  geom_line() +
  xlab("Distance from TSS (kilo-bp)") +
  ylab("Median of Median 5(h)mc") +
  theme_classic()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.title=element_text(size=12),
        legend.text=element_text(size=11))+
  #scale_color_manual(values = wes_palette("Zissou1", n=4))
  scale_color_npg()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 6: Calculate Means of Tumor Subtypes vs Control  and Plot 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Look to see if tumor subtype behaves similarly
#load unique tumor beta values
load("./Files/01162020_unique_tumor_betas.RData")
#load metadata
targets <- read.csv("./Files/Metadata_complete.csv", stringsAsFactors = FALSE)

#subset for tumors only (can use a shortcut since data.topvar is already exclusive to tumors)

#Need to Subset for Tumors Only
GliomaSamples <- targets %>% filter(Tumor_Type=="Glioma")
EpendymomaSamples <- targets %>% filter(Tumor_Type=="Ependymoma")
EmbryonalSamples <- targets %>% filter(Tumor_Type=="Embryonal")

Glioma5hmc <- tumor_5hmc_unique[,colnames(tumor_5hmc_unique)%in%GliomaSamples$Label.ID.]
Ependymoma5hmc <- tumor_5hmc_unique[,colnames(tumor_5hmc_unique)%in%EpendymomaSamples$Label.ID.]
Embryonal5hmc <- tumor_5hmc_unique[,colnames(tumor_5hmc_unique)%in%EmbryonalSamples$Label.ID.]

Glioma5mc <- tumor_5mc_unique[,colnames(tumor_5mc_unique)%in%GliomaSamples$Label.ID.]
Ependymoma5mc <- tumor_5mc_unique[,colnames(tumor_5mc_unique)%in%EpendymomaSamples$Label.ID.]
Embryonal5mc <- tumor_5mc_unique[,colnames(tumor_5mc_unique)%in%EmbryonalSamples$Label.ID.]

median_5hmc_glioma <- apply(Glioma5hmc, 1, median, na.rm=TRUE)
median_5hmc_epend <- apply(Ependymoma5hmc, 1, median, na.rm=TRUE)
median_5hmc_emb <- apply(Embryonal5hmc, 1, median, na.rm=TRUE)

median_5mc_glioma <- apply(Glioma5mc, 1, median, na.rm=TRUE)
median_5mc_epend <- apply(Ependymoma5mc, 1, median, na.rm=TRUE)
median_5mc_emb <- apply(Embryonal5mc, 1, median, na.rm=TRUE)

subtype_medians <- as.data.frame(cbind(median_5hmc_glioma,median_5hmc_epend,
                                       median_5hmc_emb,median_5mc_glioma, median_5mc_epend,median_5mc_emb))

subtype_medians <- rownames_to_column(subtype_medians, "ID")
subtype_medians <- merge(subtype_medians,median_ctl, by ="ID" )
matched_subtype_medians <- merge(annot_df_sub, subtype_medians, "ID")
matched_subtype_medians$dist_w_dir <- as.integer(matched_subtype_medians$dist_w_dir)

#First filter to less than 2000bp away from TSS
TSS_2kb <- matched_subtype_medians %>% filter(dist_w_dir <= 2000) %>%
  filter(dist_w_dir >= (-2000)) %>%
  group_by(dist_w_dir) %>%
  summarize(M2_5hmc_Gl= median(median_5hmc_glioma), 
            M2_5hmc_Ep = median(median_5hmc_epend),
            M2_5hmc_Emb = median(median_5hmc_emb),
            M2_5mc_Gl= median(median_5mc_glioma), 
            M2_5mc_Ep = median(median_5mc_epend),
            M2_5mc_Emb = median(median_5mc_emb),
            M2_5hmc_Ctl = median(Ctl_Median_5hmC),
            M2_5mc_Ctl = median(Ctl_Median_5mC),)

#Rename column names so that legend will be cleaner
G_DF = TSS_2kb %>%
  gather(key = "Measure", value = "Beta_value", -dist_w_dir) %>%
  mutate(Measure = ifelse(Measure=="M2_5hmc_Gl", "5hmC Glioma", 
                          ifelse(Measure=="M2_5hmc_Ep", "5hmC Ependymoma",
                                 ifelse(Measure=="M2_5hmc_Emb", "5hmC Embryonal",
                                        ifelse(Measure== "M2_5hmc_Ctl", "5hmC Control",
                                               ifelse(Measure=="M2_5mc_Gl", "5mC Glioma",
                                                      ifelse(Measure== "M2_5mc_Ep", "5mC Ependymoma",
                                                             ifelse(Measure== "M2_5mc_Emb", "5mC Embryonal",
                                                                    ifelse(Measure== "M2_5hmc_Ctl", "5hmC Control",
                                                                           ifelse(Measure== "M2_5mc_Ctl", "5mC Control",NA)
                                                                                 )
                                                                           )
                                                                    )
                                                             )
                                                      )
                                        )
)))


#plot
ggplot(G_DF, aes(x = dist_w_dir, y = Beta_value, colour = Measure)) +
  geom_line() +
  xlab("Distance from TSS (kilo-bp)") +
  ylab("Median of Median 5(h)mc") +
  theme_classic()+
  #scale_color_brewer(palette ="Spectral")
  scale_color_npg()


#Repeat just with resolution at 5kb
matched_subtype_medians <- matched_subtype_medians %>% 
  mutate(dist_w_dir_kb = dist_w_dir / 1000) %>%
  mutate(dist_w_dir_kb = round(dist_w_dir_kb,1))


#First filter to less than 2000bp away from TSS
TSS_2kb <- matched_subtype_medians %>% filter(dist_w_dir_kb <=2.0) %>% 
  filter(dist_w_dir_kb >= (-2.0)) %>%
  group_by(dist_w_dir_kb) %>%
  summarize(M2_5hmc_Gl= median(median_5hmc_glioma), 
            M2_5hmc_Ep = median(median_5hmc_epend),
            M2_5hmc_Emb = median(median_5hmc_emb),
            M2_5mc_Gl= median(median_5mc_glioma), 
            M2_5mc_Ep = median(median_5mc_epend),
            M2_5mc_Emb = median(median_5mc_emb),
            M2_5hmc_Ctl = median(Ctl_Median_5hmC),
            M2_5mc_Ctl = median(Ctl_Median_5mC),)

G_DF = TSS_2kb %>%
  #select(dist_w_dir, contains("median"), contains("Median")) %>%
  gather(key = "Measure", value = "Beta_value", -dist_w_dir_kb)%>%
  mutate(Measure = ifelse(Measure=="M2_5hmc_Gl", "5hmC Glioma", 
                          ifelse(Measure=="M2_5hmc_Ep", "5hmC Ependymoma",
                                 ifelse(Measure=="M2_5hmc_Emb", "5hmC Embryonal",
                                        ifelse(Measure== "M2_5hmc_Ctl", "5hmC Control",
                                               ifelse(Measure=="M2_5mc_Gl", "5mC Glioma",
                                                      ifelse(Measure== "M2_5mc_Ep", "5mC Ependymoma",
                                                             ifelse(Measure== "M2_5mc_Emb", "5mC Embryonal",
                                                                    ifelse(Measure== "M2_5hmc_Ctl", "5hmC Control",
                                                                           ifelse(Measure== "M2_5mc_Ctl", "5mC Control",NA)
                                                                    )
                                                             )
                                                      )
                                               )
                                        )
                                 )
                          )))

ggplot(G_DF, aes(x = dist_w_dir_kb, y = Beta_value, colour = Measure)) +
  geom_line() +
  xlab("Distance from TSS (kilo-bp)") +
  ylab("Median of Median 5(h)mc") +
  theme_classic()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.title=element_text(size=12),
        legend.text=element_text(size=11))+
  #scale_color_brewer(palette ="Paired")
  scale_color_npg()



#Longer Range
matched_subtype_medians$round_dist <- round(matched_subtype_medians$dist_w_dir_kb, -1)
TSS_long_range <- matched_subtype_medians %>% group_by(round_dist) %>%
  filter(round_dist < 500)%>%
  filter(round_dist > (-500)) %>%
  summarize(M2_5hmc_Gl= median(median_5hmc_glioma), 
            M2_5hmc_Ep = median(median_5hmc_epend),
            M2_5hmc_Emb = median(median_5hmc_emb),
            M2_5mc_Gl= median(median_5mc_glioma), 
            M2_5mc_Ep = median(median_5mc_epend),
            M2_5mc_Emb = median(median_5mc_emb),
            M2_5hmc_Ctl = median(Ctl_Median_5hmC),
            M2_5mc_Ctl = median(Ctl_Median_5mC),)

G_DF = TSS_long_range %>%
  #select(dist_w_dir, contains("median"), contains("Median")) %>%
  gather(key = "Measure", value = "Beta_value", -round_dist)%>%
  mutate(Measure = ifelse(Measure=="M2_5hmc_Gl", "5hmC Glioma", 
                          ifelse(Measure=="M2_5hmc_Ep", "5hmC Ependymoma",
                                 ifelse(Measure=="M2_5hmc_Emb", "5hmC Embryonal",
                                        ifelse(Measure== "M2_5hmc_Ctl", "5hmC Control",
                                               ifelse(Measure=="M2_5mc_Gl", "5mC Glioma",
                                                      ifelse(Measure== "M2_5mc_Ep", "5mC Ependymoma",
                                                             ifelse(Measure== "M2_5mc_Emb", "5mC Embryonal",
                                                                    ifelse(Measure== "M2_5hmc_Ctl", "5hmC Control",
                                                                           ifelse(Measure== "M2_5mc_Ctl", "5mC Control",NA)
                                                                    )
                                                             )
                                                      )
                                               )
                                        )
                                 )
                          )))

#plot

ggplot(G_DF, aes(x = round_dist, y = Beta_value, colour = Measure)) +
  geom_line() +
  xlab("Distance from TSS (kilo-bp)") +
  ylab("Median of Median 5(h)mc") +
  theme_classic()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.title=element_text(size=12),
        legend.text=element_text(size=11))+
  scale_color_npg()

 # scale_color_brewer(palette ="Paired")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Intermediate Range
#First filter to less than 10,000bp away from TSS
TSS_10kb <- matched_subtype_medians %>% filter(dist_w_dir_kb <=10.0) %>% 
  filter(dist_w_dir_kb >= (-10.0)) %>%
  group_by(dist_w_dir_kb) %>%
  summarize(M2_5hmc_Gl= median(median_5hmc_glioma), 
            M2_5hmc_Ep = median(median_5hmc_epend),
            M2_5hmc_Emb = median(median_5hmc_emb),
            M2_5mc_Gl= median(median_5mc_glioma), 
            M2_5mc_Ep = median(median_5mc_epend),
            M2_5mc_Emb = median(median_5mc_emb),
            M2_5hmc_Ctl = median(Ctl_Median_5hmC),
            M2_5mc_Ctl = median(Ctl_Median_5mC),)

G_DF = TSS_10kb %>%
  gather(key = "Measure", value = "Beta_value", -dist_w_dir_kb) %>%
  #Rename column names so that legend will be cleaner
  mutate(Measure = ifelse(Measure=="M2_5hmc_Gl", "5hmC Glioma", 
                          ifelse(Measure=="M2_5hmc_Ep", "5hmC Ependymoma",
                                 ifelse(Measure=="M2_5hmc_Emb", "5hmC Embryonal",
                                        ifelse(Measure== "M2_5hmc_Ctl", "5hmC Control",
                                               ifelse(Measure=="M2_5mc_Gl", "5mC Glioma",
                                                      ifelse(Measure== "M2_5mc_Ep", "5mC Ependymoma",
                                                             ifelse(Measure== "M2_5mc_Emb", "5mC Embryonal",
                                                                    ifelse(Measure== "M2_5hmc_Ctl", "5hmC Control",
                                                                           ifelse(Measure== "M2_5mc_Ctl", "5mC Control",NA)
                                                                    )
                                                             )
                                                      )
                                               )
                                        )
                                 )
                          )))

ggplot(G_DF, aes(x = dist_w_dir_kb, y = Beta_value, colour = Measure)) +
  geom_line(size =1.0) +
  xlab("Distance from TSS (kilo-bp)") +
  ylab("Median 5(h)mc") +
  theme_classic()+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        legend.title=element_text(size=20),
        legend.text=element_text(size=18))+
  scale_x_continuous(limits=c(-10, 10), breaks= seq(-10.0, 10.0, by =1.0))+
  scale_color_manual(values = c("5hmC Glioma" = "#E64B35FF","5hmC Ependymoma" = "#4DBBD5FF",
                                "5hmC Embryonal"= "#F39B7FFF", "5hmC Control" = "#3C5488FF",
                                "5mC Glioma" = "#F39B7FFF","5mC Ependymoma" = "#8491B4FF",
                                "5mC Embryonal"= "#91D1C2FF", "5mC Control" = "#7E6148FF"))

