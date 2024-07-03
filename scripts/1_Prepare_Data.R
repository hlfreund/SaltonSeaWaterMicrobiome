#### Set WD & Load Libraries ####
getwd() # use setwd("path/to/files") if you are not in the right directory
#setwd("/Volumes/HLF_SSD/Aronson_Lab_Data/Salton_Sea/SaltonSeaWater")
suppressPackageStartupMessages({ # load packages quietly
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

#load("data/SSW_analysis.Rdata") # load Rdata to global env
#save.image("data/SSW_analysis.Rdata") # save global env to Rdata file

#### Import and Prepare Data for Analyses ####

## NOTES ABOUT DATA:
# eukaryotic counts have been removed; contaminant counts removed with decontam(); all non-Lyons-samples excluded from data set before importing into this script
# zero counts & singletons also removed

## Import ALL env plate bacterial ASV count data
bac.ASV_all<-data.frame(readRDS("data/SaltonSeawater_16S_AllData_Robject.rds", refhook = NULL))
dim(bac.ASV_all) ## has count info, ASV IDs, taxa info, and original metadata
head(bac.ASV_all)

# drop data from June 2021
bac.ASV_all<-bac.ASV_all[bac.ASV_all$SampleMonth!="June",]
"June" %in% bac.ASV_all$SampleMonth # confirm there are no Junes left

# drop technical replicates from total samps
bac.ASV_all<-bac.ASV_all[!grepl("0.2a", bac.ASV_all$SampleID),]
unique(bac.ASV_all$SampleID)

# Clean up SampleIDs
bac.ASV_all$SampleID<-gsub("m\\..*","m",bac.ASV_all$SampleID)
bac.ASV_all$SampleID<-gsub("11m","10.5m",bac.ASV_all$SampleID)

# Fix depth measurements
bac.ASV_all$Depth_m<-gsub("11","10.5",bac.ASV_all$Depth_m)

# Create ASV table
bac.ASV_table <- as.data.frame(dcast(bac.ASV_all, SampleID~ASV_ID, value.var="Count", fun.aggregate=sum)) ###
bac.ASV_table[,1:4] # counts by asvs per sample
rownames(bac.ASV_table)<-bac.ASV_table$SampleID
bac.ASV_table[1:4,1:4]

rowSums(bac.ASV_table[,-1]) # total ASVs post decontamination

# Create taxa table
bac.ASV_taxa<-as.data.frame(unique(subset(bac.ASV_all, select=c(ASV_ID,Kingdom, Phylum, Class, Order, Family, Genus, Species))))
bac.ASV_taxa[1:4,1:4]

#### Import & Update Metadata ####
# upload geochem data from Lyons lab

chem_meta<-as.data.frame(read_xlsx("data/SaltonSeawater_Lyons_Aronson_Metadata_All.xlsx", sheet="Stratification_Parameters_S"))
head(chem_meta)
chem_meta[1:4,1:4] # sanity check for gsub

# drop data from June 2021
chem_meta<-chem_meta[chem_meta$SampleMonth!="June",]
"June" %in% chem_meta$SampleMonth # confirm there are no Junes left

# drop technical replicates from total samps
chem_meta<-chem_meta[!grepl("0.2a", chem_meta$SampleID),]
unique(chem_meta$SampleID)

# Clean up SampleIDs
chem_meta$SampleID<-gsub("m\\..*","m",chem_meta$SampleID)
chem_meta$SampleID

# create color variable(s) to identify variables by colors
## color for sample type
meta1<-unique(subset(bac.ASV_all, select=c(SampleID,Sample_Type,SampleMonth,SampleYear,Depth_m,SampleSource)))
head(meta1)
dim(meta1)
rownames(meta1)<-meta1$SampleID

metadata<-merge(meta1, chem_meta, by=c("SampleID","Sample_Type","SampleMonth","SampleYear","Depth_m","SampleSource"))
head(metadata)

# create factor levels for certain groups
metadata$Depth_m<-factor(metadata$Depth_m, levels=c("0","3","4","5","7","9","10","10.5"))

metadata$SampDate<-interaction(metadata$SampleMonth,metadata$SampleYear)
head(metadata)
metadata$SampDate<-factor(metadata$SampDate, levels=c("August.2021","December.2021","April.2022"))

unique(metadata$SampleMonth)
metadata$SampleMonth<-factor(metadata$SampleMonth, levels=c("August","December","April"))
unique(metadata$SampleMonth) # sanity check

# Create color palettes to be used for gradient colors

cold2warm1<-get_palette(paste0("#",c("252A52", "66ADE5", "FFC465","BF1B0B")),k=10)
names(cold2warm1) <- levels(metadata$Depth_m)

fair_cols <- paste0("#",c("252A52", "66ADE5", "FFC465","BF1B0B"))
names(fair_cols) <- letters[1:4]
fair_ramp <- scales::colour_ramp(fair_cols)
fair_sat <- saturation(fair_ramp, 1)

# Add colors for specific variables

colorset1 = melt(c(August.2021="#ef781c",December.2021="#03045e",April.2022="#059c3f"))

colorset1$SampDate<-rownames(colorset1)
colnames(colorset1)[which(names(colorset1) == "value")] <- "SampDate_Color"
colorset1

metadata<-merge(metadata, colorset1, by="SampDate")
head(metadata)
metadata$SampDate_Color <- as.character(metadata$SampDate_Color)

rownames(metadata)<-metadata$SampleID
# save.image("data/SSW_analysis.Rdata")

#### Scale Environmental Metadata ####
head(metadata)
meta_scaled<-metadata
meta_scaled[,8:17]<-scale(meta_scaled[,8:17],center=TRUE,scale=TRUE) # only scale chem env data
head(meta_scaled)

### Merge Metadata & Count Data Together ####
# create updated file that contains all ASV data and metadata for samples & variables of interest
# ** contains UNSCALED chem data
bac.dat.all<-merge(bac.ASV_all, metadata, by=c("SampleID","Sample_Type","SampleMonth","SampleYear", "Depth_m", "SampleSource"))
#rownames(bac.dat.all)<-bac.dat.all$SampleID

# drop unnecessary variables that we do not care about
bac.dat.all<-subset(bac.dat.all, select=-c(Exposure_Duration, Exposure_Type, Deployment, ExtractionMethod, LysisType, Sample_Color))
head(bac.dat.all)
dim(bac.dat.all)

#### Compare Counts by Sample ####

raw.tot.counts<-ggplot(data=bac.ASV_table, aes(x=bac.ASV_table$SampleID, y=rowSums(bac.ASV_table[,-1]))) +
  geom_bar(stat="identity",colour="black",fill="dodgerblue")+theme_classic()+theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(title="Total ASVs per Sample",subtitle="Based on Raw ASV Counts")+ylab("Total ASVs")+xlab("SampleID")

ggsave(raw.tot.counts,filename = "figures/SSW_16S_Total_ASVs_per_Sample_barplot.png", width=15, height=12, dpi=600)

### Export Global Env for Other Scripts ####
save.image("data/SSeawater_Data_Ready.Rdata")
# ^ includes all data combined in object bac.dat.all, ASV table (samples are rows, ASVs are columns), metadata, and an ASV count table (where ASVs are rows, not columns)
# data from June 2021 has been dropped because we do not have HS/SO4 data for that timepoint
# Version Information
sessionInfo()
