
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

#load("data/Env_Seqs_All/env.seq_analysis.Rdata") # load Rdata to global env
#save.image("data/Env_Seqs_All/env.seq_analysis.Rdata") # save global env to Rdata file

#### Import and Prepare Data for Analyses ####

## Import ALL env plate bacterial ASV count data
bac.ASV_counts<-data.frame(readRDS("data/DADA2_Results_9.14.22/16S.V3V4_ASVs_Counts_dada2_Robject_copy.rds", refhook = NULL))
dim(bac.ASV_counts)
bac.ASV_counts[,1:4]

#colnames(bac.ASV_counts)<-gsub("_S[0-9]+", "", colnames(bac.ASV_counts)) # shorten sample names to match sample names in metadata file
#head(bac.ASV_counts)
colnames(bac.ASV_counts)<-gsub("_",".", colnames(bac.ASV_counts))
bac.ASV_counts$ASV_ID<-rownames(bac.ASV_counts)
head(bac.ASV_counts)
colnames(bac.ASV_counts)
#bac.ASV_counts<-bac.ASV_counts[, !duplicated(colnames(bac.ASV_counts))] # remove col duplicates
dim(bac.ASV_counts)

# Remove unwanted samples aka old Salton Seawater samples sequenced by Zymo
#remove_samples<-c("SS.OV.10m.seawater.0621", "SS.OV.2m.seawater.0621", "SS.OV.5m.seawater.0621")
#bac.ASV_counts<-bac.ASV_counts[,!(colnames(bac.ASV_counts) %in% remove_samples)]
#colnames(bac.ASV_counts)
#dim(bac.ASV_counts)

## Import ASV taxonomic data
bac.ASV_tax<-data.frame(readRDS("data/DADA2_Results_9.14.22/16S.V3V4_ASVs_Taxonomy_dada2_Robject_copy.rds", refhook = NULL))
head(bac.ASV_tax)

bac.ASV_tax[is.na(bac.ASV_tax)]<- "Unknown" # turn all NAs into "Unkowns"
bac.ASV_tax$Species<-gsub("Unknown", "unknown", bac.ASV_tax$Species) # change uppercase Unkonwn to lowercase unknown for unknown species classification
head(bac.ASV_tax)
bac.ASV_tax$ASV_ID<-rownames(bac.ASV_tax) # create ASV ID column to use for merging data frames
head(bac.ASV_tax)

#save.image("data/Env_Seqs_All/env.seq_analysis.Rdata")
#### Import metadata ####

# Import all metadata
metadata<-as.data.frame(read_excel("data/DADA2_Results_9.14.22/all_env_data.xlsx", sheet="AllOriginal_Env_Data"), header=TRUE)
head(metadata)
metadata$SampleID<-gsub("_",".", metadata$SampleID)

#metadata$SampleID<-gsub("(.*\\..*)\\..*","\\1", metadata$SampleID)
#metadata<-na.omit(metadata) # drop NAs from metadata
rownames(metadata)<-metadata$SampleID
#metadata<-subset(metadata, Project=="SaltonSea")
head(metadata)
#metadata<-subset(metadata, select=-c(Project))
#head(metadata)


#### Identify & Remove Contaminants ####
ControlDF<-metadata[metadata$Sample_Type=="Control",] # pull out samples that are controls

vector_for_decontam<-metadata$SampleOrControl # use for decontam package
# ^ tells us which are controls aka TRUE vs which are not aka FALSE

bac.ASV_counts[,-length(bac.ASV_counts)] <- as.data.frame(sapply(bac.ASV_counts[,-length(bac.ASV_counts)], as.numeric)) #convert data frame to numeric
bac.ASV_c2<-t(bac.ASV_counts[,-length(bac.ASV_counts)]) # transpose so that rows are Samples and columns are ASVs
contam_df <- isContaminant(bac.ASV_c2, neg=vector_for_decontam)

table(contam_df$contaminant) # identify contaminants aka TRUE: 1420

contam_asvs <- (contam_df[contam_df$contaminant == TRUE, ]) # pull out ASV IDs for contaminating ASVs

bac.ASV_tax[row.names(bac.ASV_tax) %in% row.names(contam_asvs),] # see which taxa are contaminants

## Create new files that EXCLUDE contaminants!!!

# making new fasta file
#contam_indices <- which(asv_fasta %in% paste0(">", contam_asvs))
#dont_want <- sort(c(contam_indices, contam_indices + 1))
#asv_fasta_no_contam <- asv_fasta[- dont_want]

# making new count table
bac.ASV_counts_no.contam <- bac.ASV_counts[!row.names(bac.ASV_counts) %in% row.names(contam_asvs), ] # drop ASVs found in contam_asvs
head(bac.ASV_counts_no.contam)

# making new taxonomy table
bac.ASV_tax.no.contam <- bac.ASV_tax[!row.names(bac.ASV_tax) %in% row.names(contam_asvs), ] # drop ASVs found in contam_asvs
head(bac.ASV_tax.no.contam)

# Remove ASVs found in Controls from samples (in addition to contaminants previously ID'd)

Control_counts<-bac.ASV_counts_no.contam[,colnames(bac.ASV_counts_no.contam) %in% ControlDF$SampleID] # see which taxa are contaminants
Control_counts
Control_counts<-Control_counts[which(rowSums(Control_counts) > 0),] # drop ASVs that don't appear in Controls
dim(Control_counts)
head(Control_counts)

bac.ASV_counts_CLEAN<-bac.ASV_counts_no.contam[!bac.ASV_counts_no.contam$ASV_ID %in% row.names(Control_counts),!colnames(bac.ASV_counts_no.contam) %in% colnames(Control_counts)]
bac.ASV_taxa_CLEAN<-bac.ASV_tax.no.contam[!bac.ASV_tax.no.contam$ASV_ID %in% row.names(Control_counts),]

# sanity check
colnames(bac.ASV_counts_CLEAN) # check for control sample IDs

## and now writing them out to files
#write(asv_fasta_no_contam, "ASVs-no-contam.fa")
write.table(bac.ASV_counts_CLEAN, "data/EnvMiSeq_F22_16S.V3V4_ASVs_Counts_NoContam.tsv",
            sep="\t", quote=F, col.names=NA)
saveRDS(bac.ASV_counts_CLEAN, file = "data/EnvMiSeq_F22_16S.V3V4_ASVs_Counts_NoContam_Robject.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

write.table(bac.ASV_taxa_CLEAN, "data/EnvMiSeq_F22_16S.V3V4_ASVs_Taxa_NoContam.tsv",
            sep="\t", quote=F, col.names=NA)
saveRDS(bac.ASV_taxa_CLEAN, file = "data/EnvMiSeq_F22_16S.V3V4_ASVs_Taxa_NoContam_Robject.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

#### Separate Sequences by Project ####
DustDF<-metadata[metadata$SampleType=="Dust",] # pull out samples that are dust
dim(DustDF)
SoilDF<-metadata[metadata$SampleType=="Soil",] # pull out samples that are soil
LungDF<-metadata[metadata$SampleType=="Lung",] # pull out samples that are lung

# find DUST samples in ASV table that are in metadata
bac.ASV_counts_CLEAN[,colnames(bac.ASV_counts_CLEAN) %in% row.names(DustDF)]
dim(bac.ASV_counts_CLEAN[,colnames(bac.ASV_counts_CLEAN) %in% row.names(DustDF)])

dust_b.ASV<-bac.ASV_counts_CLEAN[,colnames(bac.ASV_counts_CLEAN) %in% row.names(DustDF)] # find DUST samples in ASV table that are in metadata
colnames(dust_b.ASV) %in% row.names(DustDF) # sanity check that this worked

soil_b.ASV<-bac.ASV_counts_CLEAN[,colnames(bac.ASV_counts_CLEAN) %in% row.names(SoilDf)] # find DUST samples in ASV table that are in metadata
colnames(soil_b.ASV) %in% row.names(DustDF) # sanity check that this worked

lung_b.ASV<-bac.ASV_counts_CLEAN[,colnames(bac.ASV_counts_CLEAN) %in% row.names(LungDF)] # find DUST samples in ASV table that are in metadata
colnames(lung_b.ASV) %in% row.names(DustDF) # sanity check that this worked

#### Update Metadata ####
# create color variable(s) to identify variables by colors
## color for sample type
# unique(metadata$Sample_Type)
# metadata$Sample_Type<-factor(metadata$Sample_Type, levels=c("Seawater", "Soil", "Dust","Playa","Fecal", "Lung","Control"))
#
# #colorset1 = melt(c(Seawater="#1f547b",Soil="#c44536",Dust="#432818",Playa="#d00000",Fecal="#66615f",Lung="#47126b",Control="#b13d1e"))
# colorset1 = melt(c(Dust="#432818",Fecal="#66615f",Lung="#47126b",Control="#b13d1e"))
#
# colorset1$Sample_Type<-rownames(colorset1)
# colnames(colorset1)[which(names(colorset1) == "value")] <- "Sample_Color"
# colorset1

metadata<-merge(metadata, colorset1, by="Sample_Type")
head(metadata)
metadata$Sample_Color <- as.character(metadata$Sample_Color)
rownames(metadata)<-metadata$SampleID
head(metadata)
metadata$Depth_m<-factor(metadata$Depth_m, levels=c("0","2","3","4","5","7","8","9","10","11"))

# unique(metadata$SampleMonth)
# metadata$SampleMonth<-factor(metadata$SampleMonth, levels=c("Januarhy","February","March","April","May","June","August","September", "October","November","December","NA"))
#
# cold2warm1<-get_palette(paste0("#",c("252A52", "66ADE5", "FFC465","BF1B0B")),k=10)
# names(cold2warm1) <- levels(metadata$Depth_m)
#
# fair_cols <- paste0("#",c("252A52", "66ADE5", "FFC465","BF1B0B"))
# names(fair_cols) <- letters[1:4]
# fair_ramp <- scales::colour_ramp(fair_cols)
# fair_sat <- saturation(fair_ramp, 1)

### Transform Data ####
# first we merge the ASV count object and the ASV taxonomy object together by column called "ASV_ID"
## then we need to melt the separate ASV & taxonomy so that we can rbind multiple data sets

# Drop singletons & zero count ASVs
dim(bac.ASV_counts_CLEAN)
bac.ASV_counts_CLEAN<-bac.ASV_counts_CLEAN[which(rowSums(bac.ASV_counts_CLEAN[,-length(bac.ASV_counts_CLEAN)]) > 0),]
dim(bac.ASV_counts_CLEAN)

# merge CLEAN aka contaminants/controls removed count & taxa tables
bac.ASV_all<-merge(bac.ASV_counts_CLEAN,bac.ASV_taxa_CLEAN, by="ASV_ID")
head(bac.ASV_all)
dim(bac.ASV_all)
bac.ASV_all<-bac.ASV_all[, !duplicated(colnames(bac.ASV_all))] # remove col duplicates
dim(bac.ASV_all)
#bac.ASV_dat<-left_join(bac.ASV_tax.no.contam,bac.ASV_dat2, by=c("ASV_ID","Kingdom","Phylum","Class","Order","Family","Genus","Species"),all=T) # all=T to keep all rows from the dataframes, not drop rows that are missing from one DF or the other

bac.dat<-melt(bac.ASV_all)
head(bac.dat)
colnames(bac.dat)[which(names(bac.dat) == "variable")] <- "SampleID"
colnames(bac.dat)[which(names(bac.dat) == "value")] <- "Count"

# Drop all Zero counts & singletons ASVs
dim(bac.dat)
bac.dat<-bac.dat[which(bac.dat$Count > 0),]
dim(bac.dat)

#bac.ASV_dat<-bac.ASV_dat[, !duplicated(colnames(bac.ASV_dat))] # remove col duplicates
#dim(bac.ASV_dat)

#bac.asv.melt<-melt(bac.ASV_dat, id.vars=c("ASV_ID","Kingdom", "Phylum", "Class", "Order", "Family", "Genus","Species"))
#colnames(bac.asv.melt)[which(names(bac.asv.melt) == "variable")] <- "SampleID"
#colnames(bac.asv.melt)[which(names(bac.asv.melt) == "value")] <- "Count"
#head(bac.asv.melt)

# Drop unknowns and eukaryotic hits
bac.dat<-subset(bac.dat, Kingdom!="Unknown") ## drop Unknowns from Kingdom
bac.dat<-subset(bac.dat, Phylum!="Unknown") ## drop Unknowns from Phylum
head(bac.dat)
dim(bac.dat)

# Create ASV count file that is filtered of eukaryotic taxa - for later use
bac.dat.with.euks<-bac.dat
colnames(bac.dat.with.euks)

# Drop chloroplast & mitochondria seqs
bac.dat<-subset(bac.dat, Class!="Chloroplast") ## exclude Chloroplast sequences
bac.dat<-subset(bac.dat, Order!="Chloroplast") ## exclude Chloroplast sequences
bac.dat<-subset(bac.dat, Family!="Mitochondria") ## exclude Mitochondrial sequences just in case

'Chloroplast' %in% bac.dat # check if Chloroplast counts are still in df, should be false because they've been removed
'Mitochondria' %in% bac.dat # check if Mitochondria counts are still in df, should be false because they've been removed
'Undetermined' %in% bac.dat # check if undetermined taxa in data frame
#NA %in% bac.ASV_dat

head(bac.dat)
#b.dat.m<-melt(bac.ASV_dat, id.vars=c("ASV_ID","Kingdom", "Phylum", "Class", "Order", "Family", "Genus","Species"))
#head(b.dat.m)
#colnames(b.dat.m)[which(names(b.dat.m) == "variable")] <- "SampleID"
#colnames(b.dat.m)[which(names(b.dat.m) == "value")] <- "Count"
#head(b.dat.m)

# Create filtered ASV table (no contaminants & no ASVs found in controls!)
bac.ASV_table<-base::as.data.frame(dcast(bac.dat, SampleID~ASV_ID, value.var="Count", fun.aggregate=sum)) ###
head(bac.ASV_table)
bac.ASV_table[duplicated(rownames(bac.ASV_table))]
rownames(bac.ASV_table)<-bac.ASV_table$SampleID
#bac.ASV_table<-subset(bac.ASV_table, select=-c(SampleID))
rownames(bac.ASV_table)

# double check dimensions of metadata and ASV table
dim(metadata)
dim(bac.ASV_table)
# double check that the rownames exist + match
rownames(metadata)
rownames(bac.ASV_table)

# Find rows in metadata that are not in combined bacterial asv tables
setdiff(rownames(metadata), rownames(bac.ASV_table)) # check rows in metadata not in bac.ASV_table
setdiff(rownames(bac.ASV_table), rownames(metadata)) # check rows in bac.ASV_table not in metadata

bac.ASV_meta<-merge(bac.ASV_table,metadata, by.x="SampleID", by.y="SampleID") # drop samples from Sierra project
dim(bac.ASV_meta)
dim(metadata) # the difference between dimensions is that the original metadata contains info for the controls
dim(bac.ASV_table)

rownames(bac.ASV_meta)<-bac.ASV_meta$SampleID

# recreate ASV table (excluding samples we don't have metadata or counts for)
colnames(bac.ASV_meta); rownames(bac.ASV_meta)
head(bac.ASV_meta)
bac.ASV_table<-subset(bac.ASV_meta, select=-c(Sample_Type,SampleMonth, SampleYear, Sample_Color, Depth_m, SampleSource, Exposure_Duration, Exposure_Type,Deployment,SampleOrControl,ExtractionMethod,LysisType))
rownames(bac.ASV_table)<-bac.ASV_table$SampleID
head(bac.ASV_table)

# triple check dimensions of metadata and ASV table
dim(metadata)
dim(bac.ASV_table)

# Find rows in metadata that are not in combined bacterial asv tables
setdiff(rownames(metadata), rownames(bac.ASV_table)) # check rows from metadata not in bac.ASV_table
# the difference between dimensions is that the original metadata contains info for the controls
setdiff(rownames(bac.ASV_table), rownames(metadata)) # check bac.ASV_table from metadata not in rows

# reorder metadata based off of ASV table
metadata=metadata[rownames(bac.ASV_table),] ## will drop rows that are not shared by both dataframes!
# here we are reordering our metadata by rows, using the rownames from our ASV table as a guide
# this indexing method will only work if the two dfs have the same # of rows AND the same row names!

# sanity check to see if this indexing step worked
head(metadata)
dim(metadata)

head(rownames(bac.ASV_table))
dim(bac.ASV_table)

## subset data by sample type
duplicated(rownames(bac.ASV_meta))

rownames(bac.ASV_meta)
head(bac.ASV_meta)
#rownames(bac.meta.melt)<-bac.meta.melt$SampleID

# Remove ASVs from taxa table
head(bac.ASV_taxa_CLEAN)
bac.tax<-bac.ASV_taxa_CLEAN[,-length(bac.ASV_taxa_CLEAN)]
head(bac.tax)

# create super metadata + taxa + counts data frame
bac.dat.meta<-merge(bac.dat,metadata,by="SampleID")

sw.asv.meta<-subset(bac.dat.meta, Sample_Type=="Seawater")
d.asv.meta<-subset(bac.dat.meta, Sample_Type=="Dust")
s.asv.meta<-subset(bac.dat.meta, Sample_Type=="Soil")
f.asv.meta<-subset(bac.dat.meta, Sample_Type=="Fecal")
l.asv.meta<-subset(bac.dat.meta, Sample_Type=="Lung")

#write.table(sw.asv.meta, "data/Env_Seqs_All/EnvSeqsAll_9.14.2022/SaltonSeawater_16S_AllData.tsv",
#            sep="\t", quote=F, col.names=NA)
#saveRDS(sw.asv.meta, file = "data/Env_Seqs_All/EnvSeqsAll_9.14.2022/SaltonSeawater_16S_AllData_Robject.rds", ascii = FALSE, version = NULL,
#        compress = TRUE, refhook = NULL)

#### Prep Dataframe for Relative Abundance ####
head(t(bac.ASV_table[,-1]))

bac.clean.counts<-as.data.frame(t(bac.ASV_table[,-1])) # all singletons, zeros, controls, contaminants dropped
head(bac.clean.counts)
bac.clean.counts$ASV_ID<-rownames(bac.clean.counts)

#bac.tax<-subset(bac.ASV_tax.no.contam, select=c("ASV_ID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
bac.all<-merge(bac.clean.counts, bac.ASV_taxa_CLEAN, by="ASV_ID")
head(bac.all)
bac_melt<-melt(bac.all, id.vars = c("ASV_ID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
head(bac_melt)
names(bac_melt)[which(names(bac_melt) == "variable")] <- "SampleID"
names(bac_melt)[which(names(bac_melt) == "value")] <- "Counts"
head(bac_melt)

head(metadata)

all_bac<-merge(bac_melt, metadata, by = "SampleID")
head(all_bac) # contains metadata, ASV counts, and taxonomic IDs for ASVs


#### Save Global Env for Import into Other Scripts ####

save.image("data/Env_Seqs_All/Env_Seq_PreppedData.Rdata") # save global env to Rdata file

## ^ this will be loaded into other scripts for downstream analyses
