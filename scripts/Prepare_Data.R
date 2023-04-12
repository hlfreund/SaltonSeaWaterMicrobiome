library(phyloseq)
library(ggplot2)
library (vegan)
library(ggpubr)
library(scales)
library(grid)
library(ape)
library(plyr)
library(dplyr)
library(viridis)
library(readxl)
library(metagenomeSeq)
library(DESeq2)
library(dplyr)
library(magrittr)
library(MASS)
library("dendextend")
library("tidyr")
library("viridis")
library("reshape")
library(reshape2)
library(wesanderson)
library(nationalparkcolors)
library(shades)

## * NOTE FOR THIS SCRIPT
# this script should be used to prepare data for input into other scripts. Allows for the dropping of Unknowns, Chloroplast, and Mitochondrial reads
# from data before running through anaylses.

#### Import metadata ####
metadata<-as.data.frame(read.csv("data/Mouse_Seqs/mapzym.csv"), header=TRUE)
head(metadata)
metadata<-na.omit(metadata)

metadata$Description<-gsub("\xca", "", metadata$Description) ## gsub is global sub (does not just remove first instance of pattern, but multiple)
rownames(metadata)<-metadata$SampleID
head(metadata)

metadata$lungTissue<-sub("yes", "lung", metadata$lungTissue)
metadata$lungTissue<-sub("no", "BALF", metadata$lungTissue)
head(metadata)

# create colors for later figures
colorset = melt(c(Alt="#540b0e",Con="#4c956c",Silica="#98c1d9"))
colorset$Exposed<-rownames(colorset)
colnames(colorset)[which(names(colorset) == "value")] <- "color"
colorset

metadata<-merge(metadata, colorset, by="Exposed")
head(metadata)
metadata$color <- as.character(metadata$color)
rownames(metadata)<-metadata$SampleID

#### Import and reformat the data ####
# bacteria first
bac.ASV_counts<-data.frame(readRDS("data/Mouse_Seqs/Results_16S/16S_ASVs_Counts_dada2_6.20.2021_Robject.rds", refhook = NULL))
dim(bac.ASV_counts)
head(bac.ASV_counts)
rownames(bac.ASV_counts)<-bac.ASV_counts$ASV_ID

# fungi next
its2.ASV_counts<-data.frame(readRDS("data/Mouse_Seqs/Results_ITS2/ITS2_ASVs_Counts_dada2_6.20.2021_Robject.rds", refhook = NULL))
dim(its2.ASV_counts)
head(its2.ASV_counts)
rownames(its2.ASV_counts)<-its2.ASV_counts$ASV_ID

### ^^^ *** Note to self - before you run stats, make sure your table is in a SAMPLE x ASV/SPECIES orientation!!

#### Import ASV to taxa sheet for ASV identification ####

# bacteria first
bac.ASV_tax<-data.frame(readRDS("data/Mouse_Seqs/Results_16S/ASVs_Taxonomy_dada2_6.16.2021_Robject.rds", refhook = NULL))
head(bac.ASV_tax)
bac.ASV_tax$Genus<-gsub("\\[(.*)\\]", "\\1", bac.ASV_tax$Genus) ## drop brackets around Eubacterium while not losing the string inside brackets

bac.ASV_tax[is.na(bac.ASV_tax)]<- "Unknown"
bac.ASV_tax$Species<-gsub("Unknown", "unknown", bac.ASV_tax$Species) ## drop brackets around Eubacterium while not losing the string inside brackets

head(bac.ASV_tax)
class(bac.ASV_tax)
bac.ASV_tax$ASV_ID<-rownames(bac.ASV_tax)
head(bac.ASV_tax)

# fungi next
its2.ASV_tax<-data.frame(readRDS("data/Mouse_Seqs/Results_ITS2/ITS2_ASVs_Taxonomy_dada2_6.16.2021_Robject.rds", refhook = NULL))
head(its2.ASV_tax)
its2.ASV_tax$ASV_ID<-rownames(its2.ASV_tax)
its2.ASV_tax<-as.data.frame(lapply(its2.ASV_tax, function(x) gsub("[a-z]__", "", x)))

its2.ASV_tax[is.na(its2.ASV_tax)]<- "Unknown"
its2.ASV_tax$Species<-gsub("Unknown", "unknown", its2.ASV_tax$Species) ## drop brackets around Eubacterium while not losing the string inside brackets

head(its2.ASV_tax)
rownames(its2.ASV_tax)<-its2.ASV_tax$ASV_ID
head(its2.ASV_tax)

#### Data Formatting and Transformation ####

# bacteria first
bac.ASV_dat<-merge(bac.ASV_counts,bac.ASV_tax, by="ASV_ID")
head(bac.ASV_dat)

bac.ASV_dat<-subset(bac.ASV_dat, Kingdom!="Unknown") ## keep only bacteria and archaean -- drop Unknowns
bac.ASV_dat<-subset(bac.ASV_dat, Phylum!="Unknown") ## keep only bacteria and archaean -- drop Unknowns=
head(bac.ASV_dat)
bac.ASV_dat<-subset(bac.ASV_dat, Class!="Chloroplast") ## keep only bacteria -- exclude Chloroplast sequences
bac.ASV_dat<-subset(bac.ASV_dat, Order!="Chloroplast") ## keep only bacteria -- exclude Chloroplast sequences
bac.ASV_dat<-subset(bac.ASV_dat, Family!="Mitochondria") ## keep only bacteria -- exclude Chloroplast sequences

'Chloroplast' %in% bac.ASV_dat # check if Chloroplast counts are still in df, should be false because they've been removed
'Mitochondria' %in% bac.ASV_dat # check if Chloroplast counts are still in df, should be false because they've been removed

head(bac.ASV_dat)
rownames(bac.ASV_dat)<-bac.ASV_dat$ASV_ID
head(bac.ASV_dat)

"Undetermined" %in% bac.ASV_dat

#### Create Sample x Species (ASV) table from counts ####
b.dat.m<-melt(bac.ASV_dat)
head(b.dat.m)
colnames(b.dat.m)[which(names(b.dat.m) == "variable")] <- "SampleID"
colnames(b.dat.m)[which(names(b.dat.m) == "value")] <- "Count"

bac.ASV_table<-as.data.frame(dcast(b.dat.m, SampleID~ASV_ID, value.var="Count", fun.aggregate=sum)) ###
head(bac.ASV_table)
rownames(bac.ASV_table)<-bac.ASV_table$SampleID
bac.ASV_table<-subset(bac.ASV_table, select=-c(SampleID))
head(bac.ASV_table)

# *** if you need a Species x Samples table, then do the following:
# bac.counts<-as.data.frame(dcast(b.dat.m, ASV_ID~SampleID, value.var="Count", fun.aggregate=sum)) ###
# rownames(bac.counts)<-bac.counts$ASV_ID
# bac.counts<-subset(bac.counts, select=-c(ASV_ID))
# head(bac.counts)

# fungi next
its2.ASV_dat<-merge(its2.ASV_counts,its2.ASV_tax, by="ASV_ID")
head(its2.ASV_dat)

its2.ASV_dat<-subset(its2.ASV_dat, Kingdom!="Unknown") ## keep only its2teria and archaean -- drop Unknowns
its2.ASV_dat<-subset(its2.ASV_dat, Phylum!="Unknown") ## keep only its2teria and archaean -- drop Unknowns=
head(its2.ASV_dat)
its2.ASV_dat<-subset(its2.ASV_dat, Class!="Chloroplast") ## keep only its2teria -- exclude Chloroplast sequences
its2.ASV_dat<-subset(its2.ASV_dat, Order!="Chloroplast") ## keep only its2teria -- exclude Chloroplast sequences
its2.ASV_dat<-subset(its2.ASV_dat, Family!="Mitochondria") ## keep only its2teria -- exclude Chloroplast sequences

'Chloroplast' %in% its2.ASV_dat # check if Chloroplast counts are still in df, should be false because they've been removed
'Mitochondria' %in% its2.ASV_dat # check if Chloroplast counts are still in df, should be false because they've been removed

head(its2.ASV_dat)
rownames(its2.ASV_dat)<-its2.ASV_dat$ASV_ID
head(its2.ASV_dat)

"Undetermined" %in% its2.ASV_dat

## Create Sample x Species (ASV) table from counts
f.dat.m<-melt(its2.ASV_dat)
head(f.dat.m)
colnames(f.dat.m)[which(names(f.dat.m) == "variable")] <- "SampleID"
colnames(f.dat.m)[which(names(f.dat.m) == "value")] <- "Count"

its2.ASV_table<-as.data.frame(dcast(f.dat.m, SampleID~ASV_ID, value.var="Count", fun.aggregate=sum)) ###
head(its2.ASV_table)
rownames(its2.ASV_table)<-its2.ASV_table$SampleID
its2.ASV_table<-subset(its2.ASV_table, select=-c(SampleID))
head(its2.ASV_table)

# *** if you need a Species x Samples table, then do the following:
# its2.counts<-as.data.frame(dcast(f.dat.m, ASV_ID~SampleID, value.var="Count", fun.aggregate=sum)) ###
# rownames(its2.counts)<-its2.counts$ASV_ID
# its2.counts<-subset(its2.counts, select=-c(ASV_ID))
# head(its2.counts)

#### Reorder metadata to have same rows as ASV tables ####
# **  this indexing method will only work if the two dfs have the same # of rows AND the same row names!

metadata=metadata[rownames(its2.ASV_table),]
bac.ASV_table=bac.ASV_table[rownames(its2.ASV_table),]
