suppressPackageStartupMessages({ # load packages quietly
  library(devtools)
  library(phyloseq)
  library(ggplot2)
  library(vegan)
  library(ggpubr)
  #library(scales)
  library(grid)
  library(ape)
  library(plyr)
  library(dplyr)
  library(viridis)
  library(ggbiplot)
  library(readxl)
  library(metagenomeSeq)
  library(DESeq2)
  library(dplyr)
  library(magrittr)
  library(MASS)
  library(dendextend)
  library(tidyr)
  library(reshape)
  library(reshape2)
  library(wesanderson)
  library(nationalparkcolors)
  library(shades)
  library(ALDEx2)
  library(rstatix)
  library(decontam)
  library(ggvegan)
  library(microbiome)
  library(pairwiseAdonis)
  library(corrplot)
  library(SpiecEasi) # Network analysis for sparse compositional data
  library(network)
  library(intergraph)
  library(ggnet)
  library(igraph)
})

#### Load Global Env to Import Count/ASV Tables ####
load("data/SSeawater_Data_Ready.Rdata") # save global env to Rdata file
#load("data/env.driver.data.Rdata") # save global env to Rdata file

bac.dat.all[1:4,1:4]
bac.ASV_table[1:4,1:4]
bac.ASV_table[(nrow(bac.ASV_table)-4):(nrow(bac.ASV_table)),(ncol(bac.ASV_table)-4):(ncol(bac.ASV_table))] # last 4 rows & cols

head(metadata)
head(meta_scaled)
#rownames(metadata)<-metadata$SampleID

# choose to either use DO % or DO in mg/L, then drop one column
meta_scaled<-subset(meta_scaled, select=-c(DO_mgL))
head(meta_scaled)

#### Create Centered Log-Ratio Table from ASV table ####
b.clr<-decostand(bac.ASV_table[,-1],method = "clr", pseudocount = 1) #CLR transformation
b.clr[1:4,1:4]
# df must have rownames are SampleIDs, columns are ASV IDs for vegan functions below\

})

#### Load Global Env to Import Count/ASV Tables ####
load("data/SSeawater_Data_Ready.Rdata") # save global env to Rdata file
#load("data/env.driver.data.Rdata") # save global env to Rdata file

bac.dat.all[1:4,1:4]
bac.ASV_table[1:4,1:4]
bac.ASV_table[(nrow(bac.ASV_table)-4):(nrow(bac.ASV_table)),(ncol(bac.ASV_table)-4):(ncol(bac.ASV_table))] # last 4 rows & cols

head(metadata)
head(meta_scaled)
#rownames(metadata)<-metadata$SampleID

# choose to either use DO % or DO in mg/L, then drop one column
meta_scaled<-subset(meta_scaled, select=-c(DO_mgL))
head(meta_scaled)

#### Create Centered Log-Ratio Table from ASV table ####
b.clr<-decostand(bac.ASV_table[,-1],method = "clr", pseudocount = 1) #CLR transformation
b.clr[1:4,1:4]
# df must have rownames are SampleIDs, columns are ASV IDs for vegan functions below\
