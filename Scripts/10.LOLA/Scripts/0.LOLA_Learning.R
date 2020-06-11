#LOLA
#subset to appropriate cell lines
#universe is top5%
#collection of interest -> filter
#Load LOLA
#Can figure out which databases may be most useful to you using this
#http://lolaweb.databio.org/
#Go to website, laod universe and region list. Go to tables and it will give you most relevant databases with info on them
library(LOLA)
#Set path where LOLA region files are located
dbPath = system.file("extdata", "hg19", package="LOLA")
#Load the region files
regionDB = loadRegionDB(dbPath)


regionDB = loadRegionDB("LOLACore/hg19")
#Annotation
names(regionDB)
#Dictinoary: dbLocation is where your region files are
#collectionAnno is what you are interestedin looking at like TFBS or DHS and where the fiels are
#regionGRL is your list of region (granges file)

#Load data and universe
#Both need to be GRANGES objects
data("sample_input", package="LOLA") # load userSets
data("sample_universe", package="LOLA") # load userUniverse

#How you run the test
# runLOLA tests for pairwise overlap between each user set and each
# region set in regionDB. It then uses a Fisher’s exact test to assess 
# significance of the overlap. The results are a data.table with several columns:
locResults = runLOLA(userSets, userUniverse, regionDB, cores=1)

setwd("/Users/nasimazizgolshani/Dropbox/christensen/PCNS_Analysis/PCNS/Script/16.LOLA")
library(simpleCache) 
