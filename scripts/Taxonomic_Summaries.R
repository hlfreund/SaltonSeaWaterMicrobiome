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

getwd()

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

metadata$exposed_y.n<-factor(metadata$exposed_y.n, levels=c("Y", "N"))

# create colors for later figures
colorset = melt(c(Alt="#540b0e",Con="#4c956c",Silica="#98c1d9"))
colorset$Exposed<-rownames(colorset)
colnames(colorset)[which(names(colorset) == "value")] <- "color"
colorset

metadata<-merge(metadata, colorset, by="Exposed")
head(metadata)
metadata$color <- as.character(metadata$color)
rownames(metadata)<-metadata$SampleID

exp_col <- melt(c('Y'="#e63946", 'N'="#004ACE"))
exp_col$exposed_y.n<-rownames(exp_col)
names(exp_col)[which(names(exp_col) == "value")] <- "color2"
exp_col

metadata<-merge(metadata, exp_col, by="exposed_y.n")
head(metadata)
metadata$color2 <- as.character(metadata$color2)
rownames(metadata)<-metadata$SampleID
head(metadata)

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

## Create Sample x Species (ASV) table from counts
b.dat.m<-melt(bac.ASV_dat)
head(b.dat.m)
colnames(b.dat.m)[which(names(b.dat.m) == "variable")] <- "SampleID"
colnames(b.dat.m)[which(names(b.dat.m) == "value")] <- "Count"

bac.ASV_table<-as.data.frame(dcast(b.dat.m, SampleID~ASV_ID, value.var="Count", fun.aggregate=sum)) ###
head(bac.ASV_table)
rownames(bac.ASV_table)<-bac.ASV_table$SampleID
bac.ASV_table<-subset(bac.ASV_table, select=-c(SampleID))
head(bac.ASV_table)

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

## Reorder metadata to have same rows as ASV tables
metadata=metadata[rownames(its2.ASV_table),]
# ** ^ this indexing method will only work if the two dfs have the same # of rows AND the same row names!
bac.ASV_table=bac.ASV_table[rownames(its2.ASV_table),]

#### 16S: Relative Abundance by Taxa Level ####
head(b.dat.m)
### * below we use the dcast() function to "cast" the data into a wide format based on given elements (column names), taking sum of "Count"
### * decostand(df, method="total) is the function (with argument total) used to get relative abundance of OTU table
## phylum ....
b.phyla_counts <- as.data.frame(dcast(b.dat.m, SampleID~Phylum, value.var="Count", fun.aggregate=sum)) ###
head(b.phyla_counts) # counts by phyla per sample
rownames(b.phyla_counts)<-b.phyla_counts$SampleID
b.phyla_counts$SampleID<-NULL
b.phyla_RelAb<-data.frame(decostand(b.phyla_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(b.phyla_RelAb) # sanity check to make sure the transformation worked!
b.phyla_RelAb$SampleID<-rownames(b.phyla_RelAb)
head(b.phyla_RelAb)
b.phyla_m<-melt(b.phyla_RelAb)

head(b.phyla_m)
colnames(b.phyla_m)[which(names(b.phyla_m) == "variable")] <- "Phylum"
colnames(b.phyla_m)[which(names(b.phyla_m) == "value")] <- "Count"
head(b.phyla_m) ## relative abundance based on sum of counts by phyla!
b.phyla_m<-subset(b.phyla_m, SampleID!="Undetermined_S0") ## keep only bacteria -- exclude Chloroplast sequences

b.phyla_RA_meta<-merge(b.phyla_m,metadata, by="SampleID")
head(b.phyla_RA_meta) ## relative abundance based on sum of counts by phyla!

## class...

b.class_counts <- as.data.frame(dcast(b.dat.m, SampleID~Class, value.var="Count", fun.aggregate=sum)) ###
head(b.class_counts) # counts by class per sample
rownames(b.class_counts)<-b.class_counts$SampleID
b.class_counts$SampleID<-NULL
b.class_RelAb<-data.frame(decostand(b.class_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(b.class_RelAb) # sanity check to make sure the transformation worked!
b.class_RelAb$SampleID<-rownames(b.class_RelAb)
head(b.class_RelAb)
b.class_m<-melt(b.class_RelAb)

head(b.class_m)
colnames(b.class_m)[which(names(b.class_m) == "variable")] <- "Class"
colnames(b.class_m)[which(names(b.class_m) == "value")] <- "Count"
head(b.class_m) ## relative abundance based on sum of counts by class!
b.class_m<-subset(b.class_m, SampleID!="Undetermined_S0") ## keep only bacteria -- exclude Chloroplast sequences

b.class_RA_meta<-merge(b.class_m,metadata, by="SampleID")
head(b.class_RA_meta) ## relative abundance based on sum of counts by class!

## order ...

b.order_counts <- as.data.frame(dcast(b.dat.m, SampleID~Order, value.var="Count", fun.aggregate=sum)) ###
head(b.order_counts) # counts by order per sample
rownames(b.order_counts)<-b.order_counts$SampleID
b.order_counts$SampleID<-NULL
b.order_RelAb<-data.frame(decostand(b.order_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(b.order_RelAb) # sanity check to make sure the transformation worked!
b.order_RelAb$SampleID<-rownames(b.order_RelAb)
head(b.order_RelAb)
b.order_m<-melt(b.order_RelAb)

head(b.order_m)
colnames(b.order_m)[which(names(b.order_m) == "variable")] <- "Order"
colnames(b.order_m)[which(names(b.order_m) == "value")] <- "Count"
head(b.order_m) ## relative abundance based on sum of counts by order!
b.order_m<-subset(b.order_m, SampleID!="Undetermined_S0") ## keep only bacteria -- exclude Chloroplast sequences

b.order_RA_meta<-merge(b.order_m,metadata, by="SampleID")
head(b.order_RA_meta) ## relative abundance based on sum of counts by order!

## family ...

b.family_counts <- as.data.frame(dcast(b.dat.m, SampleID~Family, value.var="Count", fun.aggregate=sum)) ###
head(b.family_counts) # counts by family per sample
rownames(b.family_counts)<-b.family_counts$SampleID
b.family_counts$SampleID<-NULL
b.family_RelAb<-data.frame(decostand(b.family_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(b.family_RelAb) # sanity check to make sure the transformation worked!
b.family_RelAb$SampleID<-rownames(b.family_RelAb)
head(b.family_RelAb)
b.family_m<-melt(b.family_RelAb)

head(b.family_m)
colnames(b.family_m)[which(names(b.family_m) == "variable")] <- "Family"
colnames(b.family_m)[which(names(b.family_m) == "value")] <- "Count"
head(b.family_m) ## relative abundance based on sum of counts by family!
b.family_m<-subset(b.family_m, SampleID!="Undetermined_S0") ## keep only bacteria -- exclude Chloroplast sequences

b.family_RA_meta<-merge(b.family_m,metadata, by="SampleID")
head(b.family_RA_meta) ## relative abundance based on sum of counts by family!

## genus ...

b.genus_counts <- as.data.frame(dcast(b.dat.m, SampleID~Genus, value.var="Count", fun.aggregate=sum)) ###
head(b.genus_counts) # counts by genus per sample
rownames(b.genus_counts)<-b.genus_counts$SampleID
b.genus_counts$SampleID<-NULL
b.genus_RelAb<-data.frame(decostand(b.genus_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(b.genus_RelAb) # sanity check to make sure the transformation worked!
b.genus_RelAb$SampleID<-rownames(b.genus_RelAb)
head(b.genus_RelAb)
b.genus_m<-melt(b.genus_RelAb)

head(b.genus_m)
colnames(b.genus_m)[which(names(b.genus_m) == "variable")] <- "Genus"
colnames(b.genus_m)[which(names(b.genus_m) == "value")] <- "Count"
head(b.genus_m) ## relative abundance based on sum of counts by genus!
b.genus_m<-subset(b.genus_m, SampleID!="Undetermined_S0") ## keep only bacteria -- exclude Chloroplast sequences

b.genus_RA_meta<-merge(b.genus_m,metadata, by="SampleID")
head(b.genus_RA_meta) ## relative abundance based on sum of counts by genus!

## species ...

b.species_counts <- as.data.frame(dcast(b.dat.m, SampleID~Genus+Species, value.var="Count", fun.aggregate=sum)) ###
head(b.species_counts) # counts by species per sample
rownames(b.species_counts)<-b.species_counts$SampleID
b.species_counts$SampleID<-NULL
b.species_RelAb<-data.frame(decostand(b.species_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(b.species_RelAb) # sanity check to make sure the transformation worked!
b.species_RelAb$SampleID<-rownames(b.species_RelAb)
head(b.species_RelAb)
b.species_m<-melt(b.species_RelAb)

head(b.species_m)
colnames(b.species_m)[which(names(b.species_m) == "variable")] <- "Genus_Species"
colnames(b.species_m)[which(names(b.species_m) == "value")] <- "Count"
b.species_m$Genus_Species<-gsub("_", " ", b.species_m$Genus_Species) ## gsub is global sub (does not just remove first instance of pattern, but multiple)
head(b.species_m) ## relative abundance based on sum of counts by species!
b.species_m<-subset(b.species_m, SampleID!="Undetermined_S0") ## keep only bacteria -- exclude Chloroplast sequences

b.species_RA_meta<-merge(b.species_m,metadata, by="SampleID")
head(b.species_RA_meta) ## relative abundance based on sum of counts by species!

#### ITS2: Relative Abundance by Taxa Level ####
head(f.dat.m)
### * below we use the dcast() function to "cast" the data into a wide format based on given elements (column names), taking sum of "Count"
### * decostand(df, method="total) is the function (with argument total) used to get relative abundance of OTU table

## phylum ....
f.phyla_counts <- as.data.frame(dcast(f.dat.m, SampleID~Phylum, value.var="Count", fun.aggregate=sum)) ###
head(f.phyla_counts) # counts by phyla per sample
rownames(f.phyla_counts)<-f.phyla_counts$SampleID
f.phyla_counts$SampleID<-NULL
f.phyla_RelAb<-data.frame(decostand(f.phyla_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(f.phyla_RelAb) # sanity check to make sure the transformation worked!
f.phyla_RelAb$SampleID<-rownames(f.phyla_RelAb)
head(f.phyla_RelAb)
f.phyla_m<-melt(f.phyla_RelAb)

head(f.phyla_m)
colnames(f.phyla_m)[which(names(f.phyla_m) == "variable")] <- "Phylum"
colnames(f.phyla_m)[which(names(f.phyla_m) == "value")] <- "Count"
head(f.phyla_m) ## relative abundance based on sum of counts by phyla!

f.phyla_RA_meta<-merge(f.phyla_m,metadata, by="SampleID")
head(f.phyla_RA_meta) ## relative abundance based on sum of counts by phyla!

## class...

f.class_counts <- as.data.frame(dcast(f.dat.m, SampleID~Class, value.var="Count", fun.aggregate=sum)) ###
head(f.class_counts) # counts by class per sample
rownames(f.class_counts)<-f.class_counts$SampleID
f.class_counts$SampleID<-NULL
f.class_RelAb<-data.frame(decostand(f.class_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(f.class_RelAb) # sanity check to make sure the transformation worked!
f.class_RelAb$SampleID<-rownames(f.class_RelAb)
head(f.class_RelAb)
f.class_m<-melt(f.class_RelAb)

head(f.class_m)
colnames(f.class_m)[which(names(f.class_m) == "variable")] <- "Class"
colnames(f.class_m)[which(names(f.class_m) == "value")] <- "Count"
head(f.class_m) ## relative abundance based on sum of counts by class!

f.class_RA_meta<-merge(f.class_m,metadata, by="SampleID")
head(f.class_RA_meta) ## relative abundance based on sum of counts by class!

## order ...

f.order_counts <- as.data.frame(dcast(f.dat.m, SampleID~Order, value.var="Count", fun.aggregate=sum)) ###
head(f.order_counts) # counts by order per sample
rownames(f.order_counts)<-f.order_counts$SampleID
f.order_counts$SampleID<-NULL
f.order_RelAb<-data.frame(decostand(f.order_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(f.order_RelAb) # sanity check to make sure the transformation worked!
f.order_RelAb$SampleID<-rownames(f.order_RelAb)
head(f.order_RelAb)
f.order_m<-melt(f.order_RelAb)

head(f.order_m)
colnames(f.order_m)[which(names(f.order_m) == "variable")] <- "Order"
colnames(f.order_m)[which(names(f.order_m) == "value")] <- "Count"
head(f.order_m) ## relative abundance based on sum of counts by order!

f.order_RA_meta<-merge(f.order_m,metadata, by="SampleID")
head(f.order_RA_meta) ## relative abundance based on sum of counts by order!

## family ...

f.family_counts <- as.data.frame(dcast(f.dat.m, SampleID~Family, value.var="Count", fun.aggregate=sum)) ###
head(f.family_counts) # counts by family per sample
rownames(f.family_counts)<-f.family_counts$SampleID
f.family_counts$SampleID<-NULL
f.family_RelAb<-data.frame(decostand(f.family_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(f.family_RelAb) # sanity check to make sure the transformation worked!
f.family_RelAb$SampleID<-rownames(f.family_RelAb)
head(f.family_RelAb)
f.family_m<-melt(f.family_RelAb)

head(f.family_m)
colnames(f.family_m)[which(names(f.family_m) == "variable")] <- "Family"
colnames(f.family_m)[which(names(f.family_m) == "value")] <- "Count"
head(f.family_m) ## relative abundance based on sum of counts by family!

f.family_RA_meta<-merge(f.family_m,metadata, by="SampleID")
head(f.family_RA_meta) ## relative abundance based on sum of counts by family!

## genus ...

f.genus_counts <- as.data.frame(dcast(f.dat.m, SampleID~Genus, value.var="Count", fun.aggregate=sum)) ###
head(f.genus_counts) # counts by genus per sample
rownames(f.genus_counts)<-f.genus_counts$SampleID
f.genus_counts$SampleID<-NULL
f.genus_RelAb<-data.frame(decostand(f.genus_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(f.genus_RelAb) # sanity check to make sure the transformation worked!
f.genus_RelAb$SampleID<-rownames(f.genus_RelAb)
head(f.genus_RelAb)
f.genus_m<-melt(f.genus_RelAb)

head(f.genus_m)
colnames(f.genus_m)[which(names(f.genus_m) == "variable")] <- "Genus"
colnames(f.genus_m)[which(names(f.genus_m) == "value")] <- "Count"
head(f.genus_m) ## relative abundance based on sum of counts by genus!

f.genus_RA_meta<-merge(f.genus_m,metadata, by="SampleID")
head(f.genus_RA_meta) ## relative abundance based on sum of counts by genus!

## species ...

f.species_counts <- as.data.frame(dcast(f.dat.m, SampleID~Genus+Species, value.var="Count", fun.aggregate=sum)) ###
head(f.species_counts) # counts by species per sample
rownames(f.species_counts)<-f.species_counts$SampleID
f.species_counts$SampleID<-NULL
f.species_RelAb<-data.frame(decostand(f.species_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(f.species_RelAb) # sanity check to make sure the transformation worked!
f.species_RelAb$SampleID<-rownames(f.species_RelAb)
head(f.species_RelAb)
f.species_m<-melt(f.species_RelAb)

head(f.species_m)
colnames(f.species_m)[which(names(f.species_m) == "variable")] <- "Genus_Species"
colnames(f.species_m)[which(names(f.species_m) == "value")] <- "Count"
f.species_m$Genus_Species<-gsub("_", " ", f.species_m$Genus_Species) ## gsub is global sub (does not just remove first instance of pattern, but multiple)

head(f.species_m) ## relative abundance based on sum of counts by species!

f.species_RA_meta<-merge(f.species_m,metadata, by="SampleID")
head(f.species_RA_meta) ## relative abundance based on sum of counts by species!

#### 16S: Taxonomic Summaries by Exposure Material ####

head(b.phyla_m)

#major_taxa_proportions_tab <- apply(major_taxa_counts_tab, 2, function(x) x/sum(x)*100)
b.phy.all<-merge(b.phyla_m, metadata, by="SampleID")
head(b.phy.all)

b.phy_1<-subset(b.phyla_m,c(Count)>(1/100)) ## DROP BACTERIAL PHYLA that are less than 1% abundant!!!!!!1\
b.phy1.all<-merge(b.phy_1, metadata, by="SampleID")
head(b.phy1.all)

# 16S phyla

b1<-ggplot(b.phy.all, aes(Phylum, Count)) +
  geom_jitter(aes(color=factor(Exposed)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Exposure Material", labels=c("Alt"="Alternaria", "Con"="Control", "Silica"="Silica"),
                     values=unique(b.phy.all$color[order(b.phy.all$Exposed)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Bacterial Phyla", y="Relative Abundance", title="Bacteria/Archaea & Exposure Material")
ggsave(b1,filename = "figures/16S_taxa.summary_phyla_exp_6.20.21.pdf", width=10, height=8, dpi=600)

b2<-ggplot(b.phy1.all, aes(Phylum, Count)) +
  geom_jitter(aes(color=factor(Exposed)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Exposure Material", labels=c("Alt"="Alternaria", "Con"="Control", "Silica"="Silica"),
                     values=unique(b.phy1.all$color[order(b.phy1.all$Exposed)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Bacterial Phyla", y="Relative Abundance", title="Bacteria/Archaea & Exposure Material", subtitle="Includes Phyla > 1% Relative Abundance")
ggsave(b2,filename = "figures/16S_taxa.summary_phyla.1perc_exp_6.20.21.pdf", width=10, height=8, dpi=600)


# Class
head(b.class_m)

b.cls.dat<-merge(b.class_m, metadata, by="SampleID")
head(b.cls.dat)

b.cls_1<-subset(b.cls.dat,c(Count)>(1/100)) ## DROP BACTERIAL clsLA that are less than 1% abundant!!!!!!1\
head(b.cls_1)

# 16S Class

b3<-ggplot(b.cls.dat, aes(Class, Count)) +
  geom_jitter(aes(color=factor(Exposed)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Exposure Material", labels=c("Alt"="Alternaria", "Con"="Control", "Silica"="Silica"),
                     values=unique(b.cls.dat$color[order(b.cls.dat$Exposed)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Classes", y="Relative Abundance", title="Bacteria/Archaea & Exposure Material")
ggsave(b3,filename = "figures/16S_taxa.summary_class_exp_6.22.21.pdf", width=13, height=8, dpi=600)

b4<-ggplot(b.cls_1, aes(Class, Count)) +
  geom_jitter(aes(color=factor(Exposed)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Exposure Material", labels=c("Alt"="Alternaria", "Con"="Control", "Silica"="Silica"),
                     values=unique(b.cls_1$color[order(b.cls_1$Exposed)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Classes", y="Relative Abundance", title="Bacteria/Archaea & Exposure Material", subtitle="Includes Classes > 1% Relative Abundance")
ggsave(b4,filename = "figures/16S_taxa.summary_class.1perc_exp_6.22.21.pdf", width=10, height=8, dpi=600)

# 16S Genera
head(b.genus_m)

b.gen.dat<-merge(b.genus_m, metadata, by="SampleID")
head(b.gen.dat)
b.gen_5<-subset(b.gen.dat, c(Count)>(5/100)) ## DROP BACTERIAL GENUS that are less than 1% abundant!!!!!!1\
b.gen_10<-subset(b.gen.dat, c(Count)>(10/100)) ## DROP BACTERIAL GENUS that are less than 10% abundant!!!!!!1

b5<-ggplot(b.gen_5, aes(Genus, Count)) +
  geom_jitter(aes(color=factor(Exposed)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Exposure Material", labels=c("Alt"="Alternaria", "Con"="Control", "Silica"="Silica"),
                     values=unique(b.gen_5$color[order(b.gen_5$Exposed)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Bacterial Genera", y="Relative Abundance", title="Bacteria/Archaea & Exposure Material", subtitle="Includes Genera > 5% Relative Abundance")
ggsave(b5,filename = "figures/16S_taxa.summary_genus.5perc_exp_6.20.21.pdf", width=10, height=8, dpi=600)


# Grouped boxplots, no jitter points: do the following commented out block
b5a<-ggplot(b.gen_5, aes(Genus, Count), fill=Exposed) +
  scale_fill_manual(name ="Exposure Material",labels=c("Alt"="Alternaria", "Con"="Control", "Silica"="Silica"),
                    values=unique(b.gen_5$color[order(b.gen_5$Exposed)])) +
  geom_boxplot(aes(fill=factor(Exposed))) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Bacteria/Archaea & Exposure Material", subtitle="Includes Genera > 5% Relative Abundance")
ggsave(b5a,filename = "figures/16S_taxa.summary_Genus.5perc.2_exp_6.22.21.pdf", width=12, height=8, dpi=600)

b6<-ggplot(b.gen_10, aes(Genus, Count)) +
  geom_jitter(aes(color=factor(Exposed)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Exposure Material", labels=c("Alt"="Alternaria", "Con"="Control", "Silica"="Silica"),
                     values=unique(b.gen_10$color[order(b.gen_10$Exposed)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Bacterial Genera", y="Relative Abundance", title="Bacteria/Archaea & Exposure Material", subtitle="Includes Genera > 10% Relative Abundance")
ggsave(b6,filename = "figures/16S_taxa.summary_genus.10perc_exp_6.20.21.pdf", width=10, height=8, dpi=600)

# Species
head(b.species_m) ## relative abundance based on sum of counts by specla!
#bac.spec_m$Genus_species<-gsub('_', ' ', bac.spec_m$Genus_species)
b.species_m$Genus_Species<-gsub('Unknown', 'unknown', b.species_m$Genus_Species)
b.species_m$Genus_Species<-gsub('unknown unknown', 'Unknown', b.species_m$Genus_Species)

head(b.species_m)

b.spec.dat<-merge(b.species_m, metadata, by="SampleID")
head(b.spec.dat)

b.spec_5<-subset(b.spec.dat, c(Count)>(5/100)) ## DROP BACTERIAL specLA that are less than 1% abundant!!!!!!1\
head(b.spec_5)

b.spec_10<-subset(b.spec.dat, c(Count)>(10/100)) ## DROP BACTERIAL specLA that are less than 1% abundant!!!!!!1\
head(b.spec_10)

# 16S Species

b7<-ggplot(b.spec.dat, aes(Genus_Species, Count)) +
  geom_jitter(aes(color=factor(Exposed)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Exposure Material", labels=c("Alt"="Alternaria", "Con"="Control", "Silica"="Silica"),
                     values=unique(b.spec.dat$color[order(b.spec.dat$Exposed)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Species", y="Relative Abundance", title="Bacteria/Archaea & Exposure Material")
ggsave(b7,filename = "figures/16S_taxa.summary_Species_exp_6.22.21.pdf", width=45, height=10, dpi=600)

b8<-ggplot(b.spec_5, aes(Genus_Species, Count)) +
  geom_jitter(aes(color=factor(Exposed)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Exposure Material", labels=c("Alt"="Alternaria", "Con"="Control", "Silica"="Silica"),
                     values=unique(b.spec_5$color[order(b.spec_5$Exposed)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Species", y="Relative Abundance", title="Bacteria/Archaea & Exposure Material", subtitle="Includes Species > 5% Relative Abundance")
ggsave(b8,filename = "figures/16S_taxa.summary_Species.5perc_exp_6.22.21.pdf", width=10, height=10, dpi=600)

# Grouped boxplots, no jitter points: do the following commented out block
b8.1<-ggplot(b.spec_5, aes(Genus_Species, Count), fill=Elevation) +
  scale_fill_manual(name ="Exposure Material", labels=c("Alt"="Alternaria", "Con"="Control", "Silica"="Silica"),
                    values=unique(b.spec_5$color[order(b.spec_5$Exposed)])) +
  geom_boxplot(aes(fill=factor(Exposed))) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Species", y="Relative Abundance", title="Bacteria/Archaea & Exposure Material", subtitle="Includes Species > 5% Relative Abundance")
ggsave(b8.1,filename = "figures/16S_taxa.summary_Species.5perc.2_exp_6.22.21.pdf", width=12, height=10, dpi=600)

b9<-ggplot(b.spec_10, aes(Genus_Species, Count)) +
  geom_jitter(aes(color=factor(Exposed)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Exposure Material", labels=c("Alt"="Alternaria", "Con"="Control", "Silica"="Silica"),
                     values=unique(b.spec_10$color[order(b.spec_10$Exposed)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Species", y="Relative Abundance", title="Bacteria/Archaea & Elevation", subtitle="Includes Species > 10% Relative Abundance")
ggsave(b9,filename = "figures/16S_taxa.summary_Species.10perc_exp_6.22.21.pdf", width=10, height=10, dpi=600)

b10<-ggplot(b.spec_10, aes(Genus_Species, Count), fill=Exposed) +
  scale_fill_manual(name ="Exposure Material", labels=c("Alt"="Alternaria", "Con"="Control", "Silica"="Silica"),
                    values=unique(b.spec_10$color[order(b.spec_10$Exposed)])) +
  geom_boxplot(aes(fill=factor(Exposed))) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Species", y="Relative Abundance", title="Bacteria/Archaea & Exposure Material", subtitle="Includes Species > 10% Relative Abundance")
ggsave(b10,filename = "figures/16S_taxa.summary_Species.10perc.2_exp_6.22.21.pdf", width=10, height=10, dpi=600)



#### ITS2: Taxonomic Summaries by Exposure Material ####

head(f.phyla_m)

#major_taxa_proportions_tab <- apply(major_taxa_counts_tab, 2, function(x) x/sum(x)*100)
f.phy.all<-merge(f.phyla_m, metadata, by="SampleID")
head(f.phy.all)

f.phy_1<-subset(f.phyla_m,c(Count)>(1/100)) ## DROP Fungal PHYLA that are less than 1% abundant!!!!!!1\
f.phy1.all<-merge(f.phy_1, metadata, by="SampleID")
head(f.phy1.all)

# ITS2 phyla

f1<-ggplot(f.phy.all, aes(Phylum, Count)) +
  geom_jitter(aes(color=factor(Exposed)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Exposure Material", labels=c("Alt"="Alternaria", "Con"="Control", "Silica"="Silica"),
                     values=unique(f.phy.all$color[order(f.phy.all$Exposed)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Fungal Phyla", y="Relative Abundance", title="Fungi & Exposure Material")
ggsave(f1,filename = "figures/ITS2_taxa.summary_phyla_exp_6.20.21.pdf", width=10, height=8, dpi=600)

f2<-ggplot(f.phy1.all, aes(Phylum, Count)) +
  geom_jitter(aes(color=factor(Exposed)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Exposure Material", labels=c("Alt"="Alternaria", "Con"="Control", "Silica"="Silica"),
                     values=unique(f.phy1.all$color[order(f.phy1.all$Exposed)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Fungal Phyla", y="Relative Abundance", title="Fungi & Exposure Material", subtitle="Includes Phyla > 1% Relative Abundance")
ggsave(f2,filename = "figures/ITS2_taxa.summary_phyla.1perc_exp_6.20.21.pdf", width=10, height=8, dpi=600)


# Class
head(f.class_m)

f.cls.dat<-merge(f.class_m, metadata, by="SampleID")
head(f.cls.dat)

f.cls_1<-subset(f.cls.dat,c(Count)>(1/100)) ## DROP Fungal clsLA that are less than 1% abundant!!!!!!1\
head(f.cls_1)

# ITS2 Class

f3<-ggplot(f.cls.dat, aes(Class, Count)) +
  geom_jitter(aes(color=factor(Exposed)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Exposure Material", labels=c("Alt"="Alternaria", "Con"="Control", "Silica"="Silica"),
                     values=unique(f.cls.dat$color[order(f.cls.dat$Exposed)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Classes", y="Relative Abundance", title="Fungi & Exposure Material")
ggsave(f3,filename = "figures/ITS2_taxa.summary_class_exp_6.22.21.pdf", width=13, height=8, dpi=600)

f4<-ggplot(f.cls_1, aes(Class, Count)) +
  geom_jitter(aes(color=factor(Exposed)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Exposure Material", labels=c("Alt"="Alternaria", "Con"="Control", "Silica"="Silica"),
                     values=unique(f.cls_1$color[order(f.cls_1$Exposed)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Classes", y="Relative Abundance", title="Fungi & Exposure Material", subtitle="Includes Classes > 1% Relative Abundance")
ggsave(f4,filename = "figures/ITS2_taxa.summary_class.1perc_exp_6.22.21.pdf", width=10, height=8, dpi=600)

# ITS2 Genera
head(f.genus_m)

f.gen.dat<-merge(f.genus_m, metadata, by="SampleID")
head(f.gen.dat)
f.gen_5<-subset(f.gen.dat, c(Count)>(5/100)) ## DROP Fungal GENUS that are less than 1% abundant!!!!!!1\
f.gen_10<-subset(f.gen.dat, c(Count)>(10/100)) ## DROP Fungal GENUS that are less than 10% abundant!!!!!!1

f5<-ggplot(f.gen_5, aes(Genus, Count)) +
  geom_jitter(aes(color=factor(Exposed)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Exposure Material", labels=c("Alt"="Alternaria", "Con"="Control", "Silica"="Silica"),
                     values=unique(f.gen_5$color[order(f.gen_5$Exposed)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Fungal Genera", y="Relative Abundance", title="Fungi & Exposure Material", subtitle="Includes Genera > 5% Relative Abundance")
ggsave(f5,filename = "figures/ITS2_taxa.summary_genus.5perc_exp_6.20.21.pdf", width=10, height=8, dpi=600)


# Grouped boxplots, no jitter points: do the following commented out block
f5a<-ggplot(f.gen_5, aes(Genus, Count), fill=Exposed) +
  scale_fill_manual(name ="Exposure Material",labels=c("Alt"="Alternaria", "Con"="Control", "Silica"="Silica"),
                    values=unique(f.gen_5$color[order(f.gen_5$Exposed)])) +
  geom_boxplot(aes(fill=factor(Exposed))) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Fungi & Exposure Material", subtitle="Includes Genera > 5% Relative Abundance")
ggsave(f5a,filename = "figures/ITS2_taxa.summary_Genus.5perc.2_exp_6.22.21.pdf", width=12, height=8, dpi=600)

f6<-ggplot(f.gen_10, aes(Genus, Count)) +
  geom_jitter(aes(color=factor(Exposed)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Exposure Material", labels=c("Alt"="Alternaria", "Con"="Control", "Silica"="Silica"),
                     values=unique(f.gen_10$color[order(f.gen_10$Exposed)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Fungal Genera", y="Relative Abundance", title="Fungi & Exposure Material", subtitle="Includes Genera > 10% Relative Abundance")
ggsave(f6,filename = "figures/ITS2_taxa.summary_genus.10perc_exp_6.20.21.pdf", width=10, height=8, dpi=600)

# Species
head(f.species_m) ## relative abundance based on sum of counts by specla!
#bac.spec_m$Genus_species<-gsub('_', ' ', bac.spec_m$Genus_species)
f.species_m$Genus_Species<-gsub('Unknown', 'unknown', f.species_m$Genus_Species)
f.species_m$Genus_Species<-gsub('unknown unknown', 'Unknown', f.species_m$Genus_Species)

head(f.species_m)

f.spec.dat<-merge(f.species_m, metadata, by="SampleID")
head(f.spec.dat)

f.spec_5<-subset(f.spec.dat, c(Count)>(5/100)) ## DROP Fungal specLA that are less than 1% abundant!!!!!!1\
head(f.spec_5)

f.spec_10<-subset(f.spec.dat, c(Count)>(10/100)) ## DROP Fungal specLA that are less than 1% abundant!!!!!!1\
head(f.spec_10)

# ITS2 Species

f7<-ggplot(f.spec.dat, aes(Genus_Species, Count)) +
  geom_jitter(aes(color=factor(Exposed)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Exposure Material", labels=c("Alt"="Alternaria", "Con"="Control", "Silica"="Silica"),
                     values=unique(f.spec.dat$color[order(f.spec.dat$Exposed)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Species", y="Relative Abundance", title="Fungi & Exposure Material")
ggsave(f7,filename = "figures/ITS2_taxa.summary_Species_exp_6.22.21.pdf", width=25, height=10, dpi=600)

f8<-ggplot(f.spec_5, aes(Genus_Species, Count)) +
  geom_jitter(aes(color=factor(Exposed)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Exposure Material", labels=c("Alt"="Alternaria", "Con"="Control", "Silica"="Silica"),
                     values=unique(f.spec_5$color[order(f.spec_5$Exposed)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Species", y="Relative Abundance", title="Fungi & Exposure Material", subtitle="Includes Species > 5% Relative Abundance")
ggsave(f8,filename = "figures/ITS2_taxa.summary_Species.5perc_exp_6.22.21.pdf", width=10, height=10, dpi=600)

# Grouped boxplots, no jitter points: do the following commented out block
f8.1<-ggplot(f.spec_5, aes(Genus_Species, Count), fill=Elevation) +
  scale_fill_manual(name ="Exposure Material", labels=c("Alt"="Alternaria", "Con"="Control", "Silica"="Silica"),
                    values=unique(f.spec_5$color[order(f.spec_5$Exposed)])) +
  geom_boxplot(aes(fill=factor(Exposed))) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Species", y="Relative Abundance", title="Fungi & Exposure Material", subtitle="Includes Species > 5% Relative Abundance")
ggsave(f8.1,filename = "figures/ITS2_taxa.summary_Species.5perc.2_exp_6.22.21.pdf", width=12, height=10, dpi=600)

f9<-ggplot(f.spec_10, aes(Genus_Species, Count)) +
  geom_jitter(aes(color=factor(Exposed)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Exposure Material", labels=c("Alt"="Alternaria", "Con"="Control", "Silica"="Silica"),
                     values=unique(f.spec_10$color[order(f.spec_10$Exposed)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Species", y="Relative Abundance", title="Fungi & Elevation", subtitle="Includes Species > 10% Relative Abundance")
ggsave(f9,filename = "figures/ITS2_taxa.summary_Species.10perc_exp_6.22.21.pdf", width=10, height=10, dpi=600)

f10<-ggplot(f.spec_10, aes(Genus_Species, Count), fill=Exposed) +
  scale_fill_manual(name ="Exposure Material", labels=c("Alt"="Alternaria", "Con"="Control", "Silica"="Silica"),
                    values=unique(f.spec_10$color[order(f.spec_10$Exposed)])) +
  geom_boxplot(aes(fill=factor(Exposed))) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Species", y="Relative Abundance", title="Fungi & Exposure Material", subtitle="Includes Species > 10% Relative Abundance")
ggsave(f10,filename = "figures/ITS2_taxa.summary_Species.10perc.2_exp_6.22.21.pdf", width=10, height=10, dpi=600)



#### 16S: Taxonomic Summaries by exposed_y.n: Yes or No ####

head(b.phyla_m)

#major_taxa_proportions_tab <- apply(major_taxa_counts_tab, 2, function(x) x/sum(x)*100)
b.phy.all<-merge(b.phyla_m, metadata, by="SampleID")
head(b.phy.all)

b.phy_1<-subset(b.phyla_m,c(Count)>(1/100)) ## DROP BACTERIAL PHYLA that are less than 1% abundant!!!!!!1\
b.phy1.all<-merge(b.phy_1, metadata, by="SampleID")
head(b.phy1.all)

# 16S phyla

b1<-ggplot(b.phy.all, aes(Phylum, Count)) +
  geom_jitter(aes(color=factor(exposed_y.n)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Exposed Condition", labels=c("Y"="Yes", "N"="No"),
                     values=unique(b.phy.all$color2[order(b.phy.all$exposed_y.n)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Bacterial Phyla", y="Relative Abundance", title="Bacteria/Archaea & Exposed Condition")
ggsave(b1,filename = "figures/16S_taxa.summary_phyla_exp.Y.N_6.20.21.pdf", width=10, height=8, dpi=600)

b2<-ggplot(b.phy1.all, aes(Phylum, Count)) +
  geom_jitter(aes(color=factor(exposed_y.n)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Exposed Condition", labels=c("Y"="Yes", "N"="No"),
                     values=unique(b.phy1.all$color2[order(b.phy1.all$exposed_y.n)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Bacterial Phyla", y="Relative Abundance", title="Bacteria/Archaea & Exposed Condition", subtitle="Includes Phyla > 1% Relative Abundance")
ggsave(b2,filename = "figures/16S_taxa.summary_phyla.1perc_exp.Y.N_6.20.21.pdf", width=10, height=8, dpi=600)


# Class
head(b.class_m)

b.cls.dat<-merge(b.class_m, metadata, by="SampleID")
head(b.cls.dat)

b.cls_1<-subset(b.cls.dat,c(Count)>(1/100)) ## DROP BACTERIAL clsLA that are less than 1% abundant!!!!!!1\
head(b.cls_1)

# 16S Class

b3<-ggplot(b.cls.dat, aes(Class, Count)) +
  geom_jitter(aes(color=factor(exposed_y.n)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Exposed Condition", labels=c("Y"="Yes", "N"="No"),
                     values=unique(b.cls.dat$color2[order(b.cls.dat$exposed_y.n)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Classes", y="Relative Abundance", title="Bacteria/Archaea & Exposed Condition")
ggsave(b3,filename = "figures/16S_taxa.summary_class_exp.Y.N_6.22.21.pdf", width=13, height=8, dpi=600)

b4<-ggplot(b.cls_1, aes(Class, Count)) +
  geom_jitter(aes(color=factor(exposed_y.n)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Exposed Condition", labels=c("Y"="Yes", "N"="No"),
                     values=unique(b.cls_1$color2[order(b.cls_1$exposed_y.n)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Classes", y="Relative Abundance", title="Bacteria/Archaea & Exposed Condition", subtitle="Includes Classes > 1% Relative Abundance")
ggsave(b4,filename = "figures/16S_taxa.summary_class.1perc_exp.Y.N_6.22.21.pdf", width=10, height=8, dpi=600)

# 16S Genera
head(b.genus_m)

b.gen.dat<-merge(b.genus_m, metadata, by="SampleID")
head(b.gen.dat)
b.gen_5<-subset(b.gen.dat, c(Count)>(5/100)) ## DROP BACTERIAL GENUS that are less than 1% abundant!!!!!!1\
b.gen_10<-subset(b.gen.dat, c(Count)>(10/100)) ## DROP BACTERIAL GENUS that are less than 10% abundant!!!!!!1

b5<-ggplot(b.gen_5, aes(Genus, Count)) +
  geom_jitter(aes(color=factor(exposed_y.n)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Exposed Condition", labels=c("Y"="Yes", "N"="No"),
                     values=unique(b.gen_5$color2[order(b.gen_5$exposed_y.n)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Bacterial Genera", y="Relative Abundance", title="Bacteria/Archaea & Exposed Condition", subtitle="Includes Genera > 5% Relative Abundance")
ggsave(b5,filename = "figures/16S_taxa.summary_genus.5perc_exp.Y.N_6.20.21.pdf", width=10, height=8, dpi=600)


# Grouped boxplots, no jitter points: do the following commented out block
b5a<-ggplot(b.gen_5, aes(Genus, Count), fill=exposed_y.n) +
  scale_fill_manual(name ="Exposed Condition",labels=c("Y"="Yes", "N"="No"),
                    values=unique(b.gen_5$color2[order(b.gen_5$exposed_y.n)])) +
  geom_boxplot(aes(fill=factor(exposed_y.n))) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Bacteria/Archaea & Exposed Condition", subtitle="Includes Genera > 5% Relative Abundance")
ggsave(b5a,filename = "figures/16S_taxa.summary_Genus.5perc.2_exp.Y.N_6.22.21.pdf", width=12, height=8, dpi=600)

b6<-ggplot(b.gen_10, aes(Genus, Count)) +
  geom_jitter(aes(color=factor(exposed_y.n)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Exposed Condition", labels=c("Y"="Yes", "N"="No"),
                     values=unique(b.gen_10$color2[order(b.gen_10$exposed_y.n)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Bacterial Genera", y="Relative Abundance", title="Bacteria/Archaea & Exposed Condition", subtitle="Includes Genera > 10% Relative Abundance")
ggsave(b6,filename = "figures/16S_taxa.summary_genus.10perc_exp.Y.N_6.20.21.pdf", width=10, height=8, dpi=600)

# Species
head(b.species_m) ## relative abundance based on sum of counts by specla!
#bac.spec_m$Genus_species<-gsub('_', ' ', bac.spec_m$Genus_species)
b.species_m$Genus_Species<-gsub('Unknown', 'unknown', b.species_m$Genus_Species)
b.species_m$Genus_Species<-gsub('unknown unknown', 'Unknown', b.species_m$Genus_Species)

head(b.species_m)

b.spec.dat<-merge(b.species_m, metadata, by="SampleID")
head(b.spec.dat)

b.spec_5<-subset(b.spec.dat, c(Count)>(5/100)) ## DROP BACTERIAL specLA that are less than 1% abundant!!!!!!1\
head(b.spec_5)

b.spec_10<-subset(b.spec.dat, c(Count)>(10/100)) ## DROP BACTERIAL specLA that are less than 1% abundant!!!!!!1\
head(b.spec_10)

# 16S Species

b7<-ggplot(b.spec.dat, aes(Genus_Species, Count)) +
  geom_jitter(aes(color=factor(exposed_y.n)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Exposed Condition", labels=c("Y"="Yes", "N"="No"),
                     values=unique(b.spec.dat$color2[order(b.spec.dat$exposed_y.n)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Species", y="Relative Abundance", title="Bacteria/Archaea & Exposed Condition")
ggsave(b7,filename = "figures/16S_taxa.summary_Species_exp.Y.N_6.22.21.pdf", width=45, height=10, dpi=600)

b8<-ggplot(b.spec_5, aes(Genus_Species, Count)) +
  geom_jitter(aes(color=factor(exposed_y.n)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Exposed Condition", labels=c("Y"="Yes", "N"="No"),
                     values=unique(b.spec_5$color2[order(b.spec_5$exposed_y.n)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Species", y="Relative Abundance", title="Bacteria/Archaea & Exposed Condition", subtitle="Includes Species > 5% Relative Abundance")
ggsave(b8,filename = "figures/16S_taxa.summary_Species.5perc_exp.Y.N_6.22.21.pdf", width=10, height=10, dpi=600)

# Grouped boxplots, no jitter points: do the following commented out block
b8.1<-ggplot(b.spec_5, aes(Genus_Species, Count), fill=Elevation) +
  scale_fill_manual(name ="Exposed Condition", labels=c("Y"="Yes", "N"="No"),
                    values=unique(b.spec_5$color2[order(b.spec_5$exposed_y.n)])) +
  geom_boxplot(aes(fill=factor(exposed_y.n))) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Species", y="Relative Abundance", title="Bacteria/Archaea & Exposed Condition", subtitle="Includes Species > 5% Relative Abundance")
ggsave(b8.1,filename = "figures/16S_taxa.summary_Species.5perc.2_exp.Y.N_6.22.21.pdf", width=12, height=10, dpi=600)

b9<-ggplot(b.spec_10, aes(Genus_Species, Count)) +
  geom_jitter(aes(color=factor(exposed_y.n)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Exposed Condition", labels=c("Y"="Yes", "N"="No"),
                     values=unique(b.spec_10$color2[order(b.spec_10$exposed_y.n)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Species", y="Relative Abundance", title="Bacteria/Archaea & Elevation", subtitle="Includes Species > 10% Relative Abundance")
ggsave(b9,filename = "figures/16S_taxa.summary_Species.10perc_exp.Y.N_6.22.21.pdf", width=10, height=10, dpi=600)

b10<-ggplot(b.spec_10, aes(Genus_Species, Count), fill=exposed_y.n) +
  scale_fill_manual(name ="Exposed Condition", labels=c("Y"="Yes", "N"="No"),
                    values=unique(b.spec_10$color2[order(b.spec_10$exposed_y.n)])) +
  geom_boxplot(aes(fill=factor(exposed_y.n))) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Species", y="Relative Abundance", title="Bacteria/Archaea & Exposed Condition", subtitle="Includes Species > 10% Relative Abundance")
ggsave(b10,filename = "figures/16S_taxa.summary_Species.10perc.2_exp.Y.N_6.22.21.pdf", width=10, height=10, dpi=600)



#### ITS2: Taxonomic Summaries by exposed_y.n: Yes or No ####

head(f.phyla_m)

#major_taxa_proportions_tab <- apply(major_taxa_counts_tab, 2, function(x) x/sum(x)*100)
f.phy.all<-merge(f.phyla_m, metadata, by="SampleID")
head(f.phy.all)

f.phy_1<-subset(f.phyla_m,c(Count)>(1/100)) ## DROP Fungal PHYLA that are less than 1% abundant!!!!!!1\
f.phy1.all<-merge(f.phy_1, metadata, by="SampleID")
head(f.phy1.all)

# ITS2 phyla

f1<-ggplot(f.phy.all, aes(Phylum, Count)) +
  geom_jitter(aes(color=factor(exposed_y.n)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Exposed Condition", labels=c("Y"="Yes", "N"="No"),
                     values=unique(f.phy.all$color2[order(f.phy.all$exposed_y.n)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Fungal Phyla", y="Relative Abundance", title="Fungi & Exposed Condition")
ggsave(f1,filename = "figures/ITS2_taxa.summary_phyla_exp.Y.N_6.20.21.pdf", width=10, height=8, dpi=600)

f2<-ggplot(f.phy1.all, aes(Phylum, Count)) +
  geom_jitter(aes(color=factor(exposed_y.n)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Exposed Condition", labels=c("Y"="Yes", "N"="No"),
                     values=unique(f.phy1.all$color2[order(f.phy1.all$exposed_y.n)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Fungal Phyla", y="Relative Abundance", title="Fungi & Exposed Condition", subtitle="Includes Phyla > 1% Relative Abundance")
ggsave(f2,filename = "figures/ITS2_taxa.summary_phyla.1perc_exp.Y.N_6.20.21.pdf", width=10, height=8, dpi=600)


# Class
head(f.class_m)

f.cls.dat<-merge(f.class_m, metadata, by="SampleID")
head(f.cls.dat)

f.cls_1<-subset(f.cls.dat,c(Count)>(1/100)) ## DROP Fungal clsLA that are less than 1% abundant!!!!!!1\
head(f.cls_1)

# ITS2 Class

f3<-ggplot(f.cls.dat, aes(Class, Count)) +
  geom_jitter(aes(color=factor(exposed_y.n)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Exposed Condition", labels=c("Y"="Yes", "N"="No"),
                     values=unique(f.cls.dat$color2[order(f.cls.dat$exposed_y.n)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Classes", y="Relative Abundance", title="Fungi & Exposed Condition")
ggsave(f3,filename = "figures/ITS2_taxa.summary_class_exp.Y.N_6.22.21.pdf", width=13, height=8, dpi=600)

f4<-ggplot(f.cls_1, aes(Class, Count)) +
  geom_jitter(aes(color=factor(exposed_y.n)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Exposed Condition", labels=c("Y"="Yes", "N"="No"),
                     values=unique(f.cls_1$color2[order(f.cls_1$exposed_y.n)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Classes", y="Relative Abundance", title="Fungi & Exposed Condition", subtitle="Includes Classes > 1% Relative Abundance")
ggsave(f4,filename = "figures/ITS2_taxa.summary_class.1perc_exp.Y.N_6.22.21.pdf", width=10, height=8, dpi=600)

# ITS2 Genera
head(f.genus_m)

f.gen.dat<-merge(f.genus_m, metadata, by="SampleID")
head(f.gen.dat)
f.gen_5<-subset(f.gen.dat, c(Count)>(5/100)) ## DROP Fungal GENUS that are less than 1% abundant!!!!!!1\
f.gen_10<-subset(f.gen.dat, c(Count)>(10/100)) ## DROP Fungal GENUS that are less than 10% abundant!!!!!!1

f5<-ggplot(f.gen_5, aes(Genus, Count)) +
  geom_jitter(aes(color=factor(exposed_y.n)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Exposed Condition", labels=c("Y"="Yes", "N"="No"),
                     values=unique(f.gen_5$color2[order(f.gen_5$exposed_y.n)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Fungal Genera", y="Relative Abundance", title="Fungi & Exposed Condition", subtitle="Includes Genera > 5% Relative Abundance")
ggsave(f5,filename = "figures/ITS2_taxa.summary_genus.5perc_exp.Y.N_6.20.21.pdf", width=10, height=8, dpi=600)


# Grouped boxplots, no jitter points: do the following commented out block
f5a<-ggplot(f.gen_5, aes(Genus, Count), fill=exposed_y.n) +
  scale_fill_manual(name ="Exposed Condition",labels=c("Y"="Yes", "N"="No"),
                    values=unique(f.gen_5$color2[order(f.gen_5$exposed_y.n)])) +
  geom_boxplot(aes(fill=factor(exposed_y.n))) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Fungi & Exposed Condition", subtitle="Includes Genera > 5% Relative Abundance")
ggsave(f5a,filename = "figures/ITS2_taxa.summary_Genus.5perc.2_exp.Y.N_6.22.21.pdf", width=12, height=8, dpi=600)

f6<-ggplot(f.gen_10, aes(Genus, Count)) +
  geom_jitter(aes(color=factor(exposed_y.n)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Exposed Condition", labels=c("Y"="Yes", "N"="No"),
                     values=unique(f.gen_10$color2[order(f.gen_10$exposed_y.n)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Fungal Genera", y="Relative Abundance", title="Fungi & Exposed Condition", subtitle="Includes Genera > 10% Relative Abundance")
ggsave(f6,filename = "figures/ITS2_taxa.summary_genus.10perc_exp.Y.N_6.20.21.pdf", width=10, height=8, dpi=600)

# Species
head(f.species_m) ## relative abundance based on sum of counts by specla!
#bac.spec_m$Genus_species<-gsub('_', ' ', bac.spec_m$Genus_species)
f.species_m$Genus_Species<-gsub('Unknown', 'unknown', f.species_m$Genus_Species)
f.species_m$Genus_Species<-gsub('unknown unknown', 'Unknown', f.species_m$Genus_Species)

head(f.species_m)

f.spec.dat<-merge(f.species_m, metadata, by="SampleID")
head(f.spec.dat)

f.spec_5<-subset(f.spec.dat, c(Count)>(5/100)) ## DROP Fungal specLA that are less than 1% abundant!!!!!!1\
head(f.spec_5)

f.spec_10<-subset(f.spec.dat, c(Count)>(10/100)) ## DROP Fungal specLA that are less than 1% abundant!!!!!!1\
head(f.spec_10)

# ITS2 Species

f7<-ggplot(f.spec.dat, aes(Genus_Species, Count)) +
  geom_jitter(aes(color=factor(exposed_y.n)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Exposed Condition", labels=c("Y"="Yes", "N"="No"),
                     values=unique(f.spec.dat$color2[order(f.spec.dat$exposed_y.n)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Species", y="Relative Abundance", title="Fungi & Exposed Condition")
ggsave(f7,filename = "figures/ITS2_taxa.summary_Species_exp.Y.N_6.22.21.pdf", width=45, height=10, dpi=600)

f8<-ggplot(f.spec_5, aes(Genus_Species, Count)) +
  geom_jitter(aes(color=factor(exposed_y.n)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Exposed Condition", labels=c("Y"="Yes", "N"="No"),
                     values=unique(f.spec_5$color2[order(f.spec_5$exposed_y.n)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Species", y="Relative Abundance", title="Fungi & Exposed Condition", subtitle="Includes Species > 5% Relative Abundance")
ggsave(f8,filename = "figures/ITS2_taxa.summary_Species.5perc_exp.Y.N_6.22.21.pdf", width=10, height=10, dpi=600)

# Grouped boxplots, no jitter points: do the following commented out block
f8.1<-ggplot(f.spec_5, aes(Genus_Species, Count), fill=Elevation) +
  scale_fill_manual(name ="Exposed Condition", labels=c("Y"="Yes", "N"="No"),
                    values=unique(f.spec_5$color2[order(f.spec_5$exposed_y.n)])) +
  geom_boxplot(aes(fill=factor(exposed_y.n))) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Species", y="Relative Abundance", title="Fungi & Exposed Condition", subtitle="Includes Species > 5% Relative Abundance")
ggsave(f8.1,filename = "figures/ITS2_taxa.summary_Species.5perc.2_exp.Y.N_6.22.21.pdf", width=12, height=10, dpi=600)

f9<-ggplot(f.spec_10, aes(Genus_Species, Count)) +
  geom_jitter(aes(color=factor(exposed_y.n)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Exposed Condition", labels=c("Y"="Yes", "N"="No"),
                     values=unique(f.spec_10$color2[order(f.spec_10$exposed_y.n)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Species", y="Relative Abundance", title="Fungi & Elevation", subtitle="Includes Species > 10% Relative Abundance")
ggsave(f9,filename = "figures/ITS2_taxa.summary_Species.10perc_exp.Y.N_6.22.21.pdf", width=10, height=10, dpi=600)

f10<-ggplot(f.spec_10, aes(Genus_Species, Count), fill=exposed_y.n) +
  scale_fill_manual(name ="Exposed Condition", labels=c("Y"="Yes", "N"="No"),
                    values=unique(f.spec_10$color2[order(f.spec_10$exposed_y.n)])) +
  geom_boxplot(aes(fill=factor(exposed_y.n))) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Species", y="Relative Abundance", title="Fungi & Exposed Condition", subtitle="Includes Species > 10% Relative Abundance")
ggsave(f10,filename = "figures/ITS2_taxa.summary_Species.10perc.2_exp.Y.N_6.22.21.pdf", width=10, height=10, dpi=600)


