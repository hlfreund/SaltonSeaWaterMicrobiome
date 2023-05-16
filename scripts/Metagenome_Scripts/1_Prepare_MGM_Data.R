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
  library(data.table)
  library(ape)
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
  library(fitdistrplus)
  library(logspline)
  library(shades)
  library(ALDEx2)
  library(rstatix)
  library(devtools)
  library(decontam)
})

#load("data/Metagenomes/Analysis/mgm_analysis.Rdata") # load Rdata to global env

#### Import Custom Functions ####

counts_to_binary <- function(dataFrame){
  new_m <- matrix(nrow=dim(dataFrame)[1],ncol = dim(dataFrame)[2]) # create new matrix w/ same rows and cols as input dataframe
  ## dim(df)[1] gives you first dimensions (x aka rows), dim(df)[2] gives you second dimensions (y aka columns)

  for( currentRow in 1:nrow(dataFrame)){ # for every row
    for( currentCol in 1:ncol(dataFrame)){ # for every column

      if ( is.na(dataFrame[currentRow, currentCol]) & is.numeric(dataFrame[currentRow, currentCol])){ # if both row and col (specifies each cell) are NA, change val to 0
        new_m[currentRow, currentCol] = 0
        # is.numeric(df[currentRow,currentCol]) is to confirm each cell contains a numeric element
      } else if( is.numeric(dataFrame[currentRow, currentCol]) & dataFrame[currentRow, currentCol] > 0){ # if both row and col (specifies each cell) are > 0, change val to 1
        new_m[currentRow, currentCol] = 1
      } else if ( is.numeric(dataFrame[currentRow, currentCol]) & dataFrame[currentRow, currentCol] == 0){ # if both row and col (specifies each cell) == 0 , change val to 0
        new_m[currentRow, currentCol] = 0
      } else if ( is.character(dataFrame[currentRow, currentCol])){ # if both row and col (specifies each cell) == 0 , change val to 0
        new_m[currentRow, currentCol] = dataFrame[currentRow, currentCol]
      }
    }
  }
  new_df <- as.data.frame(new_m) #turns matrix into dataframe
  names(new_df) <- names(dataFrame) #names rows & cols of new dataframe to be same as row names and col names from input dataframe
  rownames(new_df) <- rownames(dataFrame)
  #  new_df2=new_df[,order(ncol(new_df):1)]
  new_df2=new_df[rownames(dataFrame),colnames(dataFrame)]
  return(new_df2) # ensures only output is the new dataframe
}

# function below is to convert large data.table NA to 0; https://stackoverflow.com/questions/7235657/fastest-way-to-replace-nas-in-a-large-data-table

remove_na <- function(x){
  dm <- as.matrix(x)
  dm[is.na(dm)] <- 0
  #dm<-as.matrix(sapply(dm[,-1], as.numeric))
  new_df <- as.data.frame(dm) #turns matrix into dataframe
  if ("SampleID" %in% colnames(new_df)){
    new_df[,!names(new_df) %in% c("SampleID")]<-as.data.frame(lapply(new_df[,!names(new_df) %in% c("SampleID")],as.numeric))
  } else{
    new_df<-as.data.frame(lapply(new_df,as.numeric))
  }
  #new_df[,-1]<-lapply(new_df[,-1],as.numeric)
  names(new_df) <- names(x) #names rows & cols of new dataframe to be same as row names and col names from input dataframe
  rownames(new_df) <- rownames(x)
  colnames(new_df) <- colnames(x)
  #  new_df2=new_df[,order(ncol(new_df):1)]
  return(new_df) # ensures only output is the new dataframe
  #data.table(dm)
}

remove_na_bins <- function(x){
  dm <- as.matrix(x)
  dm[is.na(dm)] <- 0
  #dm<-as.matrix(sapply(dm[,-1], as.numeric))
  new_df <- as.data.frame(dm) #turns matrix into dataframe
  if ("SampleID" %in% colnames(new_df) && "Bin_ID" %in% colnames(new_df)){
    new_df[,!names(new_df) %in% c("SampleID","Bin_ID")]<-as.data.frame(lapply(new_df[,!names(new_df) %in% c("SampleID","Bin_ID")],as.numeric))
  } else{
    new_df<-as.data.frame(lapply(new_df,as.numeric))
  }
  #new_df[,-1]<-lapply(new_df[,-1],as.numeric)
  names(new_df) <- names(x) #names rows & cols of new dataframe to be same as row names and col names from input dataframe
  rownames(new_df) <- rownames(x)
  colnames(new_df) <- colnames(x)
  #  new_df2=new_df[,order(ncol(new_df):1)]
  return(new_df) # ensures only output is the new dataframe
  #data.table(dm)
}

##save.image("data/Metagenomes/Analysis/mgm_analysis.Rdata") # save global env to Rdata file

## Notes:
# code & info came from :
## https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#quick-start
## https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/07_practical.pdf
## https://www.reneshbedre.com/blog/deseq2.html
## https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html

#### Import Contig Coverage Data ####
mgm_fxns.cov<-fread(file = 'data/Metagenomes/Analysis/SSW_Samples_NoBins_Gene_Coverages_5.3.23.txt', sep='\t',header = TRUE)
dim(mgm_fxns.cov)
head(mgm_fxns.cov)

# update SampleID names so they match w/ metadata
mgm_fxns.cov$SampleID<-gsub("_",".", mgm_fxns.cov$SampleID)
head(mgm_fxns.cov)
mgm_fxns.cov$SampleID<-gsub("m\\..*","m",mgm_fxns.cov$SampleID)
mgm_fxns.cov$PlotID<-gsub("^SSW.","",mgm_fxns.cov$SampleID)
head(mgm_fxns.cov)

# Drop June 2021 samples
unique(mgm_fxns.cov$SampleID)
mgm_fxns.cov<-mgm_fxns.cov[!grepl("6.15.21", mgm_fxns.cov$SampleID),]
unique(mgm_fxns.cov$SampleID) #sanity check

# count occurrences of all traits across all mgms
n_occur <- data.frame(table(mgm_fxns.cov$KO_Function))
n_occur[n_occur$Freq > 1,]

# Divide gene counts by gene length to account for sample differences in assembly
## we do this because we did not co-assemble contigs, so each KO assignments across samples may not come from genes with the same length
## need to account for gene length here first (since KO assignments can come from genes of different lengths)
mgm_fxns.cov$CovPerGene<-mgm_fxns.cov$ReadsPerGene/mgm_fxns.cov$GeneLength
head(mgm_fxns.cov)


# create list of GeneIDs & KO IDs
gene_KOs<-unique(data.frame(Gene_ID=mgm_fxns.cov$Gene_ID, KO_ID=mgm_fxns.cov$KO_ID))
dim(gene_KOs)

# create Sample ID x Gene ID count table, using reads per gene that were divided by gene length
mgm_fxn.cov_table1<-as.data.frame(dcast(mgm_fxns.cov, SampleID~Gene_ID, value.var="CovPerGene"))
mgm_fxn.cov_table1[,1:4] # sanity check
rownames(mgm_fxn.cov_table1)<-mgm_fxn.cov_table1$SampleID

# convert all NAs to 0; will take a while
mgm_fxn.cov_table<-remove_na(mgm_fxn.cov_table1)
mgm_fxn.cov_table[,1:6] # sanity check

mgm.fxn.nolows<-mgm_fxn.cov_table[,which(colSums(mgm_fxn.cov_table[,-1])>=5)] # remove functions with less than 15 total counts across mgms
#mgm_fxn.binary<-counts_to_binary(mgm_fxn.cov_table[,-1]) # custom function to convert all counts to binary (presence/absence)
mgm.fxn.nolows[1:4,1:4] # sanity check

#### Import Bin Coverage Data ####
bin_fxns.cov<-fread(file = 'data/Metagenomes/Analysis/SSW_Bins_Gene_Fxns_ReadCounts_2.19.23.txt', sep='\t',header = TRUE)
dim(bin_fxns.cov)
head(bin_fxns.cov)
names(bin_fxns.cov)[names(bin_fxns.cov) == "SampleID"] <- "Bin_ID" #rename columns
names(bin_fxns.cov)[names(bin_fxns.cov) == "OriginalSampleID"] <- "SampleID"
head(bin_fxns.cov)
bin_fxns.cov$Bin_ID<-gsub("_",".", bin_fxns.cov$Bin_ID)
bin_fxns.cov$SampleID<-gsub("_",".", bin_fxns.cov$SampleID)
head(bin_fxns.cov)
bin_fxns.cov$Bin_ID<-gsub("m\\..*\\.bin","m.bin",bin_fxns.cov$Bin_ID) # removes rep number between m and .bin
bin_fxns.cov$PlotBin<-gsub("^SSW.","",bin_fxns.cov$Bin_ID)
bin_fxns.cov$SampleID<-gsub("m\\..*","m",bin_fxns.cov$SampleID) # removes rep number between m and .
bin_fxns.cov$PlotID<-gsub("^SSW.","",bin_fxns.cov$SampleID)

head(bin_fxns.cov)

# Drop June 2021 samples
unique(bin_fxns.cov$SampleID)
bin_fxns.cov<-bin_fxns.cov[!grepl("6.15.21", bin_fxns.cov$SampleID),]
unique(bin_fxns.cov$SampleID) #sanity check

# Divide gene counts by gene length to account for sample differences in assembly
## we do this because we did not co-assemble contigs, so each KO assignments across samples may not come from genes with the same length
## need to account for gene length here first (since KO assignments can come from genes of different lengths)
bin_fxns.cov$CovPerGene<-bin_fxns.cov$ReadsPerGene/bin_fxns.cov$GeneLength
head(bin_fxns.cov)

# create list of GeneIDs & KO IDs
bin.gene_KOs<-unique(data.frame(Gene_ID=bin_fxns.cov$Gene_ID, KO_ID=bin_fxns.cov$KO_ID))
dim(bin.gene_KOs)

bin_fxn.counts_table1<-as.data.frame(dcast(bin_fxns.cov, SampleID+Bin_ID~Gene_ID, value.var="CovPerGene"))
bin_fxn.counts_table<-remove_na_bins(bin_fxn.counts_table1) # need to check if this works later - 5/16/23
bin_fxn.counts_table[1:10,1:10] # sanity check
rownames(bin_fxn.counts_table)<-bin_fxn.counts_table$Bin_ID
bin_fxn.counts_table[1:200,1:20] # sanity check

# Remove unwanted samples
#remove_samples<-c("SS.OV.10m.seawater.0621", "SS.OV.2m.seawater.0621", "SS.OV.5m.seawater.0621")
#bac.ASV_counts<-bac.ASV_counts[,!(colnames(bac.ASV_counts) %in% remove_samples)]
#colnames(bac.ASV_counts)
#dim(bac.ASV_counts)

#### Import Taxonomic Annotation Data ####
mag_tax<-fread(file = 'data/Metagenomes/Analysis/SSW_MAGs_TaxoAnnotation_2.7.23.tsv', sep='\t',header = TRUE,fill=TRUE)
head(mag_tax)
mag_tax$Bin_ID<-gsub("_",".", mag_tax$Bin_ID)
mag_tax$Bin_ID<-gsub("m\\..*\\.bin","m.bin",mag_tax$Bin_ID) # removes rep number between m and .bin
mag_tax$PlotBin<-gsub("^SSW.","",mag_tax$Bin_ID)

head(mag_tax)

# Drop June 2021 samples
unique(mag_tax$Bin_ID)
mag_tax<-mag_tax[!grepl("6.15.21", mag_tax$Bin_ID),]
unique(mag_tax$Bin_ID) #sanity check

mag_tax[mag_tax==""]<-"Unknown" # replace blank cells with Unknown label
mag_tax[is.na(mag_tax)]<- "Unknown" # turn all NAs into "Unkowns"
mag_tax$Species<-gsub("Unknown", "unknown", mag_tax$Species) # change uppercase Unknown to lowercase unknown for unknown species classification
head(mag_tax)

#### Import Contig & Bin Mapped Read Counts ####
mapped_reads<-as.data.frame(read_xlsx("data/Metagenomes/Analysis/Total_Contig_Bin_Reads.xlsx", sheet="Total_Bin_Reads"))
head(mapped_reads)
mapped_reads$PlotBin<-gsub("^SSW.","",mapped_reads$Bin_ID)
mapped_reads$PlotID<-gsub("^SSW.","",mapped_reads$SampleID)

head(mapped_reads)
mapped_all<-merge(mapped_reads,mag_tax,by=c("PlotBin","Bin_ID"))
head(mapped_all)
mapped_all$RelAb_Map<-mapped_all$Bin_Mapped_Reads/mapped_all$Total_Mapped_Reads
head(mapped_all)

# Drop June 2021 samples
unique(mapped_all$SampleID)
mapped_all<-mapped_all[!grepl("6.15.21", mapped_all$SampleID),]
unique(mapped_all$SampleID) #sanity check

#save.image("data/Metagenomes/Analysis/mgm_analysis.Rdata")

#### Update Metadata ####
# upload UNSCALED geochem data from Lyons lab
# ^^ scale env variable data in respective scripts

mgm_meta<-as.data.frame(read_xlsx("data/Metagenomes/Analysis/SSW_Lyons_Aronson_Metadata_MGM_All.xlsx", sheet="Metagenomes_Metadata"))
head(mgm_meta)
mgm_meta$SampleID<-gsub("m\\..*","m",mgm_meta$SampleID) # remove rep # after m
mgm_meta$PlotID<-gsub("^SSW.","",mgm_meta$SampleID)
head(mgm_meta)

# drop data from June 2021
mgm_meta<-mgm_meta[mgm_meta$SampleMonth!="June",]
"June" %in% mgm_meta$SampleMonth # confirm there are no Junes left
unique(mgm_meta$SampleMonth) # another sanity check to be safe

# create factor levels for certain groups
unique(mgm_meta$Depth_m)
mgm_meta$Depth_m<-factor(mgm_meta$Depth_m, levels=c("0","5","10"))

unique(mgm_meta$SampleMonth)
mgm_meta$SampleMonth<-factor(mgm_meta$SampleMonth, levels=c("August","December","April"))

mgm_meta$SampDate<-interaction(mgm_meta$SampleMonth, mgm_meta$SampleYear)
head(mgm_meta)
mgm_meta$SampDate<-factor(mgm_meta$SampDate, levels=c("August.2021", "December.2021", "April.2022"))
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

colorset1 = melt(c(August.2021="#ff6f00",December.2021="#26547c",April.2022="#32cbff"))

colorset1$SampDate<-rownames(colorset1)
colorset1
colnames(colorset1)[which(names(colorset1) == "value")] <- "SampDate_Color"
colorset1

mgm_meta<-merge(mgm_meta, colorset1, by="SampDate")
head(mgm_meta)
mgm_meta$SampDate_Color <- as.character(mgm_meta$SampDate_Color)
rownames(mgm_meta)<-mgm_meta$SampleID
# #save.image("data/SSW_analysis.Rdata")

#### Scale Chem Data in Metadata ####
head(mgm_meta)
#meta_scaled<-subset(mgm_meta, select=-c(Salinity_ppt)) # drop salinity, for now
meta_scaled<-mgm_meta

head(meta_scaled)
meta_scaled[,8:17]<-as.data.frame(scale(meta_scaled[,8:17], center=TRUE, scale=TRUE)) #centering before scaling
meta_scaled

# create numeric version of depth column for later analyses
meta_scaled$Depth.num<-as.numeric(as.character(meta_scaled$Depth_m))

#save.image("data/Metagenomes/Analysis/mgm_analysis.Rdata")

#### Import Gene Info from KEGG ####

sulf.kegg<-read.table("data/Metagenomes/Analysis/Sulfur_KOs_KEGG.txt", header = TRUE, sep = "\t", dec = ".")
arsen.kegg<-read.table("data/Metagenomes/Analysis/Arsenic_KOs_KEGG.txt", header = TRUE, sep = "\t", dec = ".")
heatshock.kegg<-read.table("data/Metagenomes/Analysis/HeatShock_KOs_KEGG.txt", header = TRUE, sep = "\t", dec = ".")

#### Check Gene Distribution in MGMs ####
mgm_fxn.cov_table[,1:4] # sanity check
mgm_counts_t<-as.data.frame(t(mgm_fxn.cov_table[,-1]))
class(mgm_counts_t$SSW.12.22.21.5m)
mgm_counts_t <- as.data.frame(sapply(mgm_counts_t, as.numeric)) # convert integer to numeric across df
class(mgm_counts_t$SSW.12.22.21.5m) # sanity check

#descdist(mgm_counts_t$SSW.4.13.22.5m, discrete = TRUE)
#descdist(mgm_counts_t$SSW.12.22.21.5m, discrete = TRUE)
#descdist(mgm_counts_t$SSW.8.24.21.0m, discrete = TRUE)

# check distribution of gene # x gene counts with two samples
ggplot(mgm_counts_t) +
  geom_histogram(aes(x = SSW.4.13.22.5m), stat = "bin", bins = 200) +
  xlab("Raw gene counts") +
  ylab("Number of genes")

ggplot(mgm_counts_t) +
  geom_histogram(aes(x = SSW.12.22.21.0m), stat = "bin", bins = 200) +
  xlab("Raw gene counts") +
  ylab("Number of genes")
## data is not normally distributed

# histogram of gene counts
hist(mgm_counts_t$SSW.12.22.21.5m, col="blue")
# visualize Q-Q plot for
qqnorm(mgm_counts_t$SSW.12.22.21.5m, pch = 1, frame = FALSE)
qqline(mgm_counts_t$SSW.12.22.21.5m, col = "red", lwd = 2)

# histogram of gene counts
hist(mgm_counts_t$SSW.4.13.22.5m, col="blue")
# visualize Q-Q plot
qqnorm(mgm_counts_t$SSW.4.13.22.5m, pch = 1, frame = FALSE)
qqline(mgm_counts_t$SSW.4.13.22.5m, col = "red", lwd = 2)

# compare mean counts vs variance of the counts across samples
mean_counts <- apply(mgm_counts_t, 1, mean)        #The second argument '1' of 'apply' function indicates the function being applied to rows. Use '2' if applied to columns
variance_counts <- apply(mgm_counts_t, 1, var)
df <- data.frame(mean_counts, variance_counts)

ggplot(df) +
  geom_point(aes(x=mean_counts, y=variance_counts)) +
  scale_y_log10(limits = c(1,1e9)) +
  scale_x_log10(limits = c(1,1e9)) +
  geom_abline(intercept = 0, slope = 1, color="red")

#save.image("data/Metagenomes/Analysis/mgm_analysis.Rdata")

#### Prepare Feature Count Data for Normalization w/ DESeq2 ####
# make sure count data & mgm_meta are in the same order
head(mgm_fxns.cov)

# scale function coverage by 100 so that you can round to create DESeq2 object
## DESeq2 requires whole integers
mgm_scale.cov_table<-mgm_fxn.cov_table
mgm_fxn.cov_table[,2:6]*100 # sanity check that this method will work
mgm_scale.cov_table<-mgm_fxn.cov_table[,-1]*100

mgm_scale.cov_table[,1:10]

mgm_scale.cov_table$SampleID<-rownames(mgm_scale.cov_table)

rownames(mgm_scale.cov_table) # sanity check
rownames(mgm_meta)

mgm_meta=mgm_meta[rownames(mgm_scale.cov_table),] ## reorder mgm_meta to have same rows as original OTU table

# convert all NAs to 0; will take a while
mgm_scale.cov_table[1:6,1:6]

#convert data frame to numeric
#mgm_scale.cov_table[,-1] <- as.data.frame(sapply(mgm_scale.cov_table[,-1], as.numeric))

scale_cov_matrix<-as.matrix(mgm_scale.cov_table[,!names(mgm_scale.cov_table) %in% c("SampleID")]) # convert count table into matrix
scale_cov_matrix2<-t(scale_cov_matrix) # transpose matrix so KO_IDs are rows, samples are columns
# ^ will be used in DESeq2 functions
rownames(mgm_meta) %in% colnames(scale_cov_matrix2) # check if rownames in metadata (SampleID) match column names in count data
dim(mgm_meta)
dim(scale_cov_matrix2)

# create the DESeq DataSet object for DGE analysis
# DESeq2 needs whole # data, so need raw read counts, NOT coverage for these tables...questionable
# mgm_dds has the counts of reads that mapped to annotated genes (aka output from featureCounts)
mgm_dds<-DESeqDataSetFromMatrix(countData=round(scale_cov_matrix2), colData = mgm_meta, design = ~ 1)

# design = ~ 1 means no design
head(counts(mgm_dds))
colSums(counts(mgm_dds)) %>% barplot

# add rownames from KO IDs to DESeq Data set object
featureData <- data.frame(KO_ID=rownames(mgm_dds)) # create object of KO IDs
mcols(mgm_dds) <- DataFrame(mcols(mgm_dds), featureData) # add KO IDs to rownames of DeSEQ dataset object
mcols(mgm_dds)

# drop KO IDs with lowest counts
keep <- rowSums(counts(mgm_dds)) >= 5
mgm_dds <- mgm_dds[keep,]

# Estimate size factor - aka normalization factor, divide all read counts by each size factor per sample
mgm_dds <- estimateSizeFactors(mgm_dds,type="ratio")
## To calculate size factor in DESeq2:
# calculates geometric mean of each gene in each sample and uses this as a pseudoreference
# calculates ratio of each sample by dividing each gene count by it's pseudoreference in each sample
# The median value of all ratios for a given sample is taken as the normalization factor (size factor)
# The differentially expressed genes should not affect the median value
# median of ratios method makes the assumption that not ALL genes are differentially expressed; therefore, the normalization factors should account for sequencing depth and RNA composition of the sample
## (large outlier genes will not represent the median ratio values)
# more here: https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html

sizeFactors(mgm_dds)

plot(sizeFactors(mgm_dds), colSums(counts(mgm_dds)))

# Does sequencing depth influence normalization?
par(mfrow=c(1,2)) # to plot the two box plots next to each other
boxplot(log2(counts(mgm_dds)), notch=TRUE,
        main = "Non-normalized read counts\n(log-transformed)",
        ylab="read counts")
boxplot(log2(counts(mgm_dds, normalize= TRUE)), notch=TRUE,
        main = "Size-factor-normalized read counts\n(log-transformed)",
        ylab="read counts")
dev.off()

#save.image("data/Metagenomes/Analysis/mgm_analysis.Rdata")

#### Median-Ratio Normalized - Gene Counts ####
#
sizeFactors(mgm_dds) # will be used to normalize counts
## to normalize counts, we divide each raw count value in a given sample by that sample’s normalization factor (aka size factor) to generate normalized count values.
### This is performed for all count values (every gene in every sample)
mgm_fxn_mr <- counts(mgm_dds, normalized=TRUE) ## median-ratio transformation is how DESeq2 Normalizes counts!
#write.table(mgm_fxn_counts.norm, file="data/Metagenomes/Analysis/MGM_NoBins_MedianRatio_GeneCounts_2.27.23.txt", sep="\t", quote=F, col.names=T)

mgm.mr<-as.data.frame(t(mgm_fxn_mr))
mgm.mr$SampleID<-rownames(mgm.mr)
mgm.mr[1:4,1:4]

#### Variance Stabilizing Transformation - Gene Counts ####

# code below is also in previous section - if you ran previous section, skip to VST step
mgm_scale.cov_table[1:6,1:6]

scale_cov_matrix<-as.matrix(mgm_scale.cov_table[,!names(mgm_scale.cov_table) %in% c("SampleID")]) # convert count table into matrix
scale_cov_matrix2<-t(scale_cov_matrix) # transpose matrix so KO_IDs are rows, samples are columns
# ^ will be used in DESeq2 functions
rownames(mgm_meta) %in% colnames(scale_cov_matrix2) # check if rownames in metadata (SampleID) match column names in count data
# sanity check
dim(mgm_meta)
dim(scale_cov_matrix2)

# variance stabilizing transformation
mgm_fxn_vst <- varianceStabilizingTransformation(mgm_dds, blind = TRUE, fitType = "parametric")
assay(mgm_fxn_vst) #see output of VST

mgm.vst<-assay(mgm_fxn_vst)
#total_fxn_vst<-colSums(mgm_fxn_vst)

#### Centered Log Ratio Transformation - Gene Counts ####
mgm_fxn.cov_table[,1:4] # sanity check
# ^ table contains gene coverage, Sample IDs as rows & genes as columns

# df must have rownames are SampleIDs, columns are ASV IDs for vegan functions below\
mgm.clr<-decostand(mgm_fxn.cov_table[,-1],method = "clr", pseudocount = 1) #CLR transformation
mgm.clr[1:4,1:4]

# check rownames of CLR transformed ASV data & metadata
rownames(mgm.clr) %in% rownames(meta_scaled)

#### Copies Per Million Transformation - Gene Counts ####
mgm_fxn.cov_table[1:4,1:4] # sanity check
mgm_fxn.cov_t.table<-as.data.frame(t(as.data.frame(mgm_fxn.cov_table[,-1])))
mgm_fxn.cov_t.table[1:4,1:4]

mgm_fxn_cpm<-(mgm_fxn.cov_t.table/colSums(mgm_fxn.cov_t.table))*10^6
mgm_fxn_cpm[1:4,1:4]
colSums(mgm_fxn_cpm)

#write.table(mgm_fxn_cpm, file="data/Metagenomes/Analysis/MGM_NoBins_CopiesPerMillion_GeneCounts_2.27.23.txt", sep="\t", quote=F, col.names=T)

#### Compare Sequencing Depth Across Samples ####
mgm.mr[1:4,(ncol(mgm.mr)-4):(ncol(mgm.mr))]
total_mr_counts<-rowSums(mgm.mr[,-(ncol(mgm.mr))])

mgm_fxn_cpm[1:4,1:4]
total_cpm_counts<-colSums(mgm_fxn_cpm)

mgm.clr[1:4,1:4]
total_clr_counts<-rowSums(mgm.clr)

total_counts<-rowSums(mgm_fxn.counts_table[,-1])
total_counts %>% barplot

total_vst_counts<-colSums(mgm.vst)

total_counts %>% barplot

par(mfrow=c(1,3)) # to plot the three box plots next to each other (1 row, 3 columns)
total_counts %>% barplot(main = "Total Counts per Sample")
total_mr_counts  %>% barplot(main = "Total Median-Ratio Transformed Counts per Sample")
total_vst_counts %>% barplot(main = "Total Variance-Stabilized Transformed Counts per Sample")
#total_cpm_counts  %>% barplot(main = "Total Copies per Million (CPM) per Sample")

dev.off()

### Pull out traits of interest ####
# create unique list of KO ID and functions
# check for duplicates to make sure each KO_ID has a unique function assignment
head(mgm_fxns.cov)
NA %in% mgm_fxns.cov$CovPerGene # just to ensure there are no NAs for genes in this df

ko_fxns1<-mgm_fxns.cov[which(mgm_fxns.cov$KO_Function %in% unique(mgm_fxns.cov$KO_Function)),] # first subset out data based on unique KO functions
ko_fxns2<-subset(ko_fxns1, select=c("Gene_ID","KO_ID","KO_Function"))

n_occur <- data.frame(table(ko_fxns2$KO_ID)) # see how many duplicates there are of KO IDs, compare duplicates
n_occur[n_occur$Freq > 1,] # what traits appear more than once?

ko_ID<-unique(data.frame(KO_ID=ko_fxns2$KO_ID)) # get a list of unique KO IDs in data

ko_fxns<-as.data.frame(unique(ko_fxns2[which(ko_fxns2$KO_ID %in% ko_ID$KO_ID),])) # use unique KO ID list to subset out KO function data
head(ko_fxns)

## pull out functions of interest
#sulfur.fxns<-as.data.frame(ko_fxns[grep("sulf|thio", ko_fxns$KO_Function), ]) # pull out sulfur functions
sulfur.fxns<-ko_fxns[which(ko_fxns$KO_ID %in% sulf.kegg$KO_ID),]
osmo.fxns<-as.data.frame(ko_fxns[grep("osmo", ko_fxns$KO_Function), ]) # pull out osmoprotectant functions
arsenic.fxns<-as.data.frame(ko_fxns[grep("arsen", ko_fxns$KO_Function), ]) # pull out arsenic functions
heatshock.fxns<-as.data.frame(ko_fxns[grep("heat shock", ko_fxns$KO_Function), ]) # pull out heat shock functions
UV.res.fxns<-as.data.frame(ko_fxns[grep("^UV ", ko_fxns$KO_Function), ]) # pull out UV related functions
endospore.fxns<-as.data.frame(ko_fxns[grep("^spore ", ko_fxns$KO_Function), ]) # pull out endospore related functions

### Merge Metadata & Count Data & Taxa Data Together ####
mgm.clr[1:4,1:4]
head(mgm_meta)
head(ko_fxns)
head(mapped_all)

# melt data with normalized feature counts to merge with all traits & traits of interest
mgm.clr$SampleID<-rownames(mgm.clr)
mgm_clr_melt<-melt(mgm.clr, by="SampleID")
head(mgm_clr_melt)
colnames(mgm_clr_melt)[which(names(mgm_clr_melt) == "variable")] <- "Gene_ID"
colnames(mgm_clr_melt)[which(names(mgm_clr_melt) == "value")] <- "CovPerGene"

mgm.clr.all<-as.data.frame(merge(ko_fxns, mgm_clr_melt, by=c("Gene_ID"),allow.cartesian = TRUE))
head(mgm.clr.all)

NA %in% mgm.clr.all$KO_ID

mgm.clr.all<-mgm.clr.all[!is.na(mgm.clr.all$KO_ID),] # only looking at functions we have KO IDs for
head(mgm.clr.all)
#NA %in% mgm.clr.all

mgm.clr.sulf<-merge(sulfur.fxns, mgm_clr_melt, by=c("Gene_ID"),allow.cartesian = TRUE)
head(mgm.clr.sulf)

mgm.clr.ars<-merge(arsenic.fxns, mgm_clr_melt, by=c("Gene_ID"),allow.cartesian = TRUE)
head(mgm.clr.ars)

mgm.clr.osmo<-merge(osmo.fxns, mgm_clr_melt, by=c("Gene_ID"),allow.cartesian = TRUE)
head(mgm.clr.osmo)

# merge mapped read counts & mgm taxonomy data to metadata
mapped_meta<-merge(mapped_all, mgm_meta, by=c("SampleID"))
head(mapped_meta)

#save.image("data/Metagenomes/Analysis/mgm_analysis.Rdata")

### Export Global Env for Other Scripts ####
#save.image("data/Metagenomes/Analysis/mgm_analysis.Rdata")
# ^ includes all data combined in object bac.dat.all, ASV table (samples are rows, ASVs are columns), mgm_meta, and an ASV count table (where ASVs are rows, not columns)
# Version Information
sessionInfo()