#### Set WD & Load Libraries ####
getwd() # use setwd("path/to/files") if you are not in the right directory
setwd("/Volumes/HLF_SSD/Aronson_Lab_Data/Salton_Sea/SaltonSeaWater")
suppressPackageStartupMessages({ # load packages quietly
  library(phyloseq)
  library(ggplot2)
  library(vegan)
  library(ggpubr)
  #library(scales)
  library(grid)
  library(ape)
  library(data.table)
  library(apeglm)
  library(plyr)
  library(dplyr)
  library(viridis)
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
  library(devtools)
  library(decontam)
})

load("data/Metagenomes/Analysis/mgm_analysis.Rdata") # load Rdata to global env

#save.image("data/Metagenomes/Analysis/mgm_analysis.Rdata") # save global env to Rdata file

## Notes:
# code & info came from two links:
## https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#quick-start
## https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/07_practical.pdf
## https://www.reneshbedre.com/blog/deseq2.html
## https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html

#### Import MGM Read Counts & Taxonomic Annotation Data ####
mgm_fxns.cov<-na.omit(data.frame(read.table(file = 'data/Metagenomes/Analysis/SSW_Samples_NoBins_Gene_Coverages_2.19.23.txt', sep='\t',header = TRUE)))
dim(mgm_fxns.cov)
head(mgm_fxns.cov)
mgm_fxn.counts_table<-dcast(mgm_fxns.cov, SampleID~KO_ID, value.var="ReadsPerGene")
mgm_fxn.counts_table[1:4,1:4] # sanity check
mgm_fxn.counts_table$SampleID<-gsub("_",".", mgm_fxn.counts_table$SampleID)
rownames(mgm_fxn.counts_table)<-mgm_fxn.counts_table$SampleID
mgm_fxn.counts_table[1:4,1:4] # sanity check

bin_fxns.cov<-na.omit(data.frame(read.table(file = 'data/Metagenomes/Analysis/SSW_Bins_Gene_Fxns_ReadCounts_2.19.23.txt', sep='\t',header = TRUE)))
dim(bin_fxns.cov)
head(bin_fxns.cov)
names(bin_fxns.cov)[names(bin_fxns.cov) == "SampleID"] <- "Bin_ID" #rename columns
names(bin_fxns.cov)[names(bin_fxns.cov) == "OriginalSampleID"] <- "SampleID"
bin_fxn.counts_table<-dcast(bin_fxns.cov, SampleID+Bin_ID~KO_ID, value.var="ReadsPerGene")
bin_fxn.counts_table[1:4,1:4] # sanity check
bin_fxn.counts_table$SampleID<-gsub("_",".", bin_fxn.counts_table$SampleID)
bin_fxn.counts_table$Bin_ID<-gsub("_",".", bin_fxn.counts_table$Bin_ID)
rownames(bin_fxn.counts_table)<-bin_fxn.counts_table$Bin_ID
bin_fxn.counts_table[1:4,1:4] # sanity check

# Remove unwanted samples
#remove_samples<-c("SS.OV.10m.seawater.0621", "SS.OV.2m.seawater.0621", "SS.OV.5m.seawater.0621")
#bac.ASV_counts<-bac.ASV_counts[,!(colnames(bac.ASV_counts) %in% remove_samples)]
#colnames(bac.ASV_counts)
#dim(bac.ASV_counts)

## Import MGM taxonomic data
mag_tax<-data.frame(read.table(file = 'data/Metagenomes/Analysis/SSW_MAGs_TaxoAnnotation_2.7.23.tsv', sep='\t',header = TRUE,fill=TRUE))
head(mag_tax)

mag_tax[mag_tax==""]<-"Unknown" # replace blank cells with Unknown label
mag_tax[is.na(mag_tax)]<- "Unknown" # turn all NAs into "Unkowns"
mag_tax$Species<-gsub("Unknown", "unknown", mag_tax$Species) # change uppercase Unkonwn to lowercase unknown for unknown species classification
head(mag_tax)

save.image("data/Metagenomes/Analysis/mgm_analysis.Rdata")

#### Update Metadata ####
# upload UNSCALED geochem data from Lyons lab
# ^^ scale env variable data in respective scripts

mgm_meta<-as.data.frame(read_xlsx("data/Metagenomes/Analysis/SaltonSeawater_Lyons_Aronson_Metadata_All.xlsx", sheet="Metagenomes_Metadata"))
head(mgm_meta)

# create factor levels for certain groups
unique(mgm_meta$Depth_m)
mgm_meta$Depth_m<-factor(mgm_meta$Depth_m, levels=c("0","5","10"))

unique(mgm_meta$SampleMonth)
mgm_meta$SampleMonth<-factor(mgm_meta$SampleMonth, levels=c("June","August","December","April"))

mgm_meta$SampDate<-interaction(mgm_meta$SampleMonth, mgm_meta$SampleYear)
head(mgm_meta)
mgm_meta$SampDate<-factor(mgm_meta$SampDate, levels=c("June.2021", "August.2021", "December.2021", "April.2022"))
head(mgm_meta)

# Create color palettes to be used for gradient colors

cold2warm1<-get_palette(paste0("#",c("BF1B0B","66ADE5","252A52")),k=3)
names(cold2warm1) <- levels(mgm_meta$Depth_m)
cold2warm1

fair_cols <- paste0("#",c("252A52", "66ADE5", "FFC465","BF1B0B"))
names(fair_cols) <- letters[1:4]
fair_ramp <- scales::colour_ramp(fair_cols)
fair_sat <- saturation(fair_ramp, 1)

# Add colors for specific variables

colorset1 = melt(c(June.2021="#36ab57",August.2021="#ff6f00",December.2021="#26547c",April.2022="#32cbff"))

colorset1$SampDate<-rownames(colorset1)
colorset1
colnames(colorset1)[which(names(colorset1) == "value")] <- "SampDate_Color"
colorset1

mgm_meta<-merge(mgm_meta, colorset1, by="SampDate")
head(mgm_meta)
mgm_meta$SampDate_Color <- as.character(mgm_meta$SampDate_Color)
rownames(mgm_meta)<-mgm_meta$SampleID
# save.image("data/SSW_analysis.Rdata")

#### Scale Chem Data in Metadata ####
head(mgm_meta)
meta_scaled<-subset(mgm_meta, select=-c(Salinity_ppt)) # drop salinity

head(meta_scaled)
meta_scaled[,8:14]<-as.data.frame(scale(meta_scaled[,8:14], center=FALSE, scale=TRUE)) #not centering before scaling
head(meta_scaled)

save.image("data/Metagenomes/Analysis/mgm_analysis.Rdata")

#### Prepare Count Data for Normalization w/ DESeq2 ####
# make sure count data & mgm_meta are in the same order
mgm_meta=mgm_meta[rownames(mgm_fxn.counts_table),] ## reorder mgm_meta to have same rows as original OTU table

mgm_fxn.counts_t.table<-as.data.frame(t(mgm_fxn.counts_table[,-1]))
mgm_fxn.counts_t.table[1:4,1:4]
#mgm_counts_matrix<-as.matrix(mgm_fxn.counts_table[,-1]) # convert count table into matrix

#### Median-Ratio Transformation ####
mgm_fxn.counts_table[1:4,1:4] # sanity check

# first get log transformation of count data, where rows are KO IDs

mgm_log.counts<-log(t(mgm_fxn.counts_table[,-1]))
mgm_log.counts[1:4,]

# get rowmeans of log data & add that column do the log data df
pseudo_ref<-data.frame(Pseudo_Reference=exp(rowMeans(mgm_log.counts)))
head(pseudo_ref)
rownames(pseudo_ref)
#mgm_log.counts<-cbind(mgm_log.counts, pseudo_ref)
#head(mgm_log.counts)

# remove genes with -Inf as average from counts table & reference df

pseudo_ref.filt<-subset(pseudo_ref, Pseudo_Reference!="-Inf" & pseudo_ref > 0)
dim(pseudo_ref.filt)
dim(pseudo_ref)

mgm.fxn.logcounts.t_filt<-as.data.frame(mgm_log.counts[which(pseudo_ref!="-Inf" & pseudo_ref > 0),])
dim(mgm.fxn.logcounts.t_filt)
head(mgm.fxn.logcounts.t_filt)

mgm.fxn.counts.t_filt<-as.data.frame(t(mgm_fxn.counts_table[,-1])[which(pseudo_ref!="-Inf" & pseudo_ref > 0),])
dim(mgm.fxn.counts.t_filt)
head(mgm.fxn.counts.t_filt)

# Divide feature counts by pseudo-reference

rownames(pseudo_ref.filt)
rownames(mgm.fxn.counts.t_filt)

## check data frame order by rownames before merging dataframes
all(rownames(pseudo_ref.filt) %in% rownames(mgm.fxn.counts.t_filt))

mgm_counts_ref_filt<-cbind(mgm.fxn.counts.t_filt, pseudo_ref.filt)
head(mgm_counts_ref_filt)

mgm_fxn_ratio.data<-mgm_counts_ref_filt[,-ncol(mgm_counts_ref_filt)]/mgm_counts_ref_filt[,ncol(mgm_counts_ref_filt)]
head(mgm_fxn_ratio.data)

# Find median of ratios --> these are your scaling factors
scaling_factors<-apply(mgm_fxn_ratio.data, 2, median)

# Compare scaling factors here to DESeq2
scaling_factors
sizeFactors(mgm_dds)

# Divide original feature counts by scaling factors
manually_normalized_counts = sweep(mgm.fxn.counts.t_filt, 2, t(scaling_factors), "/")
head(manually_normalized_counts)
dim(manually_normalized_counts)

head(mgm_fxn_counts.norm)
dim(mgm_fxn_counts.norm)
# first calculate geometric mean of gene counts
mgm.fxn.geomeans[1:4,]
rownames(mgm.fxn.geomeans)

test<-t(mgm_fxn.counts_table[,-1])/mgm.fxn.geomeans
