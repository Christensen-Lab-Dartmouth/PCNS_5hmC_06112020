


Prepare4JitterIsland <- function(Array5hmc, Array5mc, Array5c, GeneAnnot){
  #Add indications of treatment which will be important when we merge the three
  Array5hmc$tx <- "5hmC"
  Array5mc$tx <- "5mC"
  Array5c$tx <- "5C"

  #Need to turn row.names to column so that CpG ID will be preserved
  Array5hmc <- rownames_to_column(Array5hmc, "ID")
  Array5mc <- rownames_to_column(Array5mc, "ID")
  Array5c <- rownames_to_column(Array5c, "ID")
  
  #subset 3 arrays to only CpGs with Gene you're interested in
  Gene5hmc <- Array5hmc[Array5hmc$ID %in% GeneAnnot$ID,];
  Gene5mc <- Array5mc[Array5mc$ID %in% GeneAnnot$ID,];
  Gene5c <- Array5c[Array5c$ID %in% GeneAnnot$ID,];
  #merge annotation and array
  Gene5hmc <- merge(Gene5hmc, GeneAnnot, by="ID");
  Gene5mc <- merge(Gene5mc, GeneAnnot, by="ID");
  Gene5c <- merge(Gene5c, GeneAnnot, by="ID");
  #limit to only  veta values, tx, Relation to Island Columns
  Gene5hmc_gather <- pivot_longer(Gene5hmc,col=grep("PCNS", colnames(Gene5hmc)),names_to = "Patient", values_to = "Beta"); 
  Gene5mc_gather <- pivot_longer(Gene5mc,col=grep("PCNS", colnames(Gene5mc)),names_to = "Patient", values_to = "Beta"); 
  Gene5c_gather <- pivot_longer(Gene5c,col=grep("PCNS", colnames(Gene5c)),names_to = "Patient", values_to = "Beta"); 
  #combine the three dataframes
  Gene_allCs <- rbind(Gene5hmc_gather,Gene5mc_gather,Gene5c_gather)
  #Subset
  Gene_allCs <- Gene_allCs[,c(1,2,3,26,27)]
  return(Gene_allCs)
}

#Write function to streamline processing step

IlluminaArrayPrep <- function(Array5hmc, targets, top5, tx){
  #select for columns you're interested in (ctl vs tumor vs subtype)
  Ctl_5hmc <- Array5hmc[,colnames(Array5hmc) %in% targets$ArrayID];
  #change colnames
  colnames(Ctl_5hmc) <-targets$Label.ID.[match(names(Ctl_5hmc),targets$ArrayID)];
  #select for only CpGs in the Tumor's Top 5% Most Hydroxymethylated Regions
  Ctl_5hmc <- Ctl_5hmc[row.names(Ctl_5hmc) %in% top5$ID,]
  #Add indications of treatment which will be important when we merge the three
  Ctl_5hmc$tx <- paste("5",tx, sep="")
  # Need to turn row.names to column so that CpG ID will be preserved
  Ctl_5hmc <- rownames_to_column(Ctl_5hmc, "ID")
  return(Ctl_5hmc)
}

#This one is better for the unique tumor file as those already have column names changed
IlluminaArrayPrepUnique <- function(Array, top5, tx){
  #select for only CpGs in the Tumor's Top 5% Most Hydroxymethylated Regions
  Array <- Array[row.names(Array) %in% top5$ID,]
  #Add indications of treatment which will be important when we merge the three
  Array$tx <- paste("5",tx, sep="")
  # Need to turn row.names to column so that CpG ID will be preserved
  Array <- rownames_to_column(Array, "ID")
  return(Array)
}


#Write function to streamline processing step
IlluminaArrayPrepNoTop <- function(Array, targets, tx){
  #select for columns you're interested in (ctl vs tumor vs subtype)
  Ctl <- Array[,colnames(Array) %in% targets$ArrayID];
  #change colnames
  colnames(Ctl) <-targets$Label.ID.[match(names(Ctl),targets$ArrayID)];
  #Add indications of treatment which will be important when we merge the three
  Ctl$tx <- paste("5",tx, sep="")
  return(Ctl)
}


#This function does not add tx type (5hmc/5mc/5c), must specify that before
ArrayPrepforTSSVisualz <- function(Array5hmc,Array5mc, Array5c, GeneAnnot, TSSAnnot){
  #rownames to ID
  Array5hmc <- rownames_to_column(Array5hmc, "ID");
  Array5mc <- rownames_to_column(Array5mc, "ID");
  Array5c <- rownames_to_column(Array5c, "ID");
  #select for CpGs related to the Gene you're interested in 
  Gene5hmc <- Array5hmc[Array5hmc$ID %in% GeneAnnot$ID,];
  Gene5mc <- Array5mc[Array5mc$ID %in% GeneAnnot$ID,];
  Gene5c <- Array5c[Array5c$ID %in% GeneAnnot$ID,];
  #Add annot data
  Gene5hmc$dist_to_tss <- TSSAnnot[Gene5hmc$ID, "dist_w_dir"];
  Gene5mc$dist_to_tss <- TSSAnnot[Gene5mc$ID, "dist_w_dir"];
  Gene5c$dist_to_tss <- TSSAnnot[Gene5c$ID, "dist_w_dir"];
  #Add tx
  Gene5hmc$tx <- "5hmC"
  Gene5mc$tx <- "5mC"
  Gene5c$tx <- "5c"
  #Collapse
  hydroxy_gather <- pivot_longer(Gene5hmc, col=grep("PCNS", colnames(Gene5hmc)), names_to = "Patient", values_to = "Beta");
  methyl_gather <- pivot_longer(Gene5mc, col=grep("PCNS", colnames(Gene5mc)), names_to = "Patient", values_to = "Beta");
  unmethyl_gather <- pivot_longer(Gene5c, col=grep("PCNS", colnames(Gene5c)), names_to = "Patient", values_to = "Beta");
  #Combine all arrays
  Allbetas <- rbind(unmethyl_gather,methyl_gather,hydroxy_gather);
  return(Allbetas)}



TumorSubtypeStratify <- function(Array, targets, Subtype){
  #Need to Subset for Tumors Only
  GliomaTargets <- targets %>% filter(Tumor_Type=="GLIOMA");
  EpendymomaTargets <- targets %>% filter(Tumor_Type=="Ependymoma");
  EmbryonalTargets <- targets %>% filter(Tumor_Type=="EMBRYONAL");
  
  Glioma <- Array[,colnames(Array)%in%GliomaTargets$Label.ID.];
  Ependymoma <- Array[,colnames(Array)%in%EpendymomaTargets$Label.ID.];
  Embryonal <- Array[,colnames(Array)%in%EmbryonalTargets$Label.ID.];
  
  if (Subtype == "Glioma") {
        return(Glioma)  } 
  else {
    if(Subtype =="Ependymoma"){
      return(Ependymoma)
    }
      else{
        return(Embryonal)
      }
  }
}
