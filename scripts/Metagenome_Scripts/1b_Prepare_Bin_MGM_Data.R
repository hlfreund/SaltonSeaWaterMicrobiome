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

##save.image("data/Metagenomes/Analysis/mgm_MAG_analysis.Rdata") # save global env to Rdata file

## Notes:
# code & info came from :
## https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#quick-start
## https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/07_practical.pdf
## https://www.reneshbedre.com/blog/deseq2.html
## https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html

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

# create list of bins & samples
bin_list<-unique(data.frame(Bin_ID=bin_fxns.cov$Bin_ID,SampleID=bin_fxns.cov$SampleID))

# Divide gene counts by gene length to account for sample differences in assembly
## we do this because we did not co-assemble contigs, so each KO assignments across samples may not come from genes with the same length
## need to account for gene length here first (since KO assignments can come from genes of different lengths)
bin_fxns.cov$CovPerGene<-bin_fxns.cov$ReadsPerGene/bin_fxns.cov$GeneLength
head(bin_fxns.cov)

# create list of GeneIDs & KO IDs
bin.gene_KOs<-unique(data.frame(Gene_ID=bin_fxns.cov$Gene_ID, KO_ID=bin_fxns.cov$KO_ID))
dim(bin.gene_KOs)
bin.gene_KOs<-na.omit(bin.gene_KOs) # drop genes that do not have KO ID
head(bin.gene_KOs)

bin_fxn.cov_table1<-as.data.frame(dcast(bin_fxns.cov, SampleID+Bin_ID~Gene_ID, value.var="CovPerGene"))
bin_fxn.cov_table<-remove_na_bins(bin_fxn.cov_table1)
bin_fxn.cov_table[1:10,1:10] # sanity check
rownames(bin_fxn.cov_table)<-bin_fxn.cov_table$Bin_ID
bin_fxn.cov_table[,1:20] # sanity check

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

#save.image("data/Metagenomes/Analysis/mgm_MAG_analysis.Rdata")

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

colorset1 = melt(c(August.2021="#ef781c",December.2021="#03045e",April.2022="#059c3f"))

colorset1$SampDate<-rownames(colorset1)
colorset1
colnames(colorset1)[which(names(colorset1) == "value")] <- "SampDate_Color"
colorset1

mgm_meta<-merge(mgm_meta, colorset1, by="SampDate")
head(mgm_meta)
mgm_meta$SampDate_Color <- as.character(mgm_meta$SampDate_Color)
rownames(mgm_meta)<-mgm_meta$SampleID
# #save.image("data/Metagenomes/Analysis/mgm_MAG_analysis.Rdata")

# create color palette for binary heat maps of function data
binary.cols<-c("1"="red","0"="aliceblue")

#### Scale Chem Data in Metadata ####
head(mgm_meta)
#meta_scaled<-subset(mgm_meta, select=-c(Salinity_ppt)) # drop salinity, for now
meta_scaled<-mgm_meta

head(meta_scaled)
meta_scaled[,8:17]<-as.data.frame(scale(meta_scaled[,8:17], center=TRUE, scale=TRUE)) #centering before scaling
meta_scaled

# create numeric version of depth column for later analyses
meta_scaled$Depth.num<-as.numeric(as.character(meta_scaled$Depth_m))

#save.image("data/Metagenomes/Analysis/mgm_MAG_analysis.Rdata")

#### Merge Overall Metadata w/ Bin Data ####

bin_meta<-merge(mgm_meta,bin_list,by="SampleID")
head(bin_meta)

rownames(bin_meta)<-bin_meta$Bin_ID

bin_meta_scaled<-merge(meta_scaled,bin_list,by="SampleID")
rownames(bin_meta_scaled)<-bin_meta$Bin_ID
bin_meta_scaled$PlotBin<-gsub("^SSW.","",bin_meta_scaled$Bin_ID)

#### Import Gene Info from KEGG ####

all_goi.kegg<-read.table("data/Metagenomes/Analysis/Genes_of_Interest_All_KOs_Pathway_Cycle_KEGG.txt", header = TRUE, sep = "\t", dec = ".")
carb.kegg<-read.table("data/Metagenomes/Analysis/CarbonFixation_KOs_KEGG.txt", header = TRUE, sep = "\t", dec = ".")
sulf.kegg<-read.table("data/Metagenomes/Analysis/Sulfur_KOs_KEGG.txt", header = TRUE, sep = "\t", dec = ".")
nitro.kegg<-read.table("data/Metagenomes/Analysis/N_KOs_Pathway_Cycle_KEGG.txt", header = TRUE, sep = "\t", dec = ".")
arsen.kegg<-read.table("data/Metagenomes/Analysis/Arsenic_KOs_KEGG.txt", header = TRUE, sep = "\t", dec = ".")
heatshock.kegg<-read.table("data/Metagenomes/Analysis/HeatShock_KOs_KEGG.txt", header = TRUE, sep = "\t", dec = ".")
osmo.kegg<-read.table("data/Metagenomes/Analysis/Osmoprot_KOs_KEGG.txt", header = TRUE, sep = "\t", dec = ".")
selen.kegg<-read.table("data/Metagenomes/Analysis/Selenium_KOs_KEGG.txt", header = TRUE, sep = "\t", dec = ".")
metal.re.kegg<-read.table("data/Metagenomes/Analysis/MetalResistance_KOs_KEGG.txt", header = TRUE, sep = "\t", dec = ".")
photo.kegg<-read.table("data/Metagenomes/Analysis/Photo_KO_KEGG.txt", header = TRUE, sep = "\t", dec = ".")
aero.kegg<-read.table("data/Metagenomes/Analysis/OxidativePhosphorylation_KO_KEGG.txt", header = TRUE, sep = "\t", dec = ".")

#### Check Gene Distribution in MAGs ####
bin_fxn.cov_table[,1:4] # sanity check
dim(bin_fxn.cov_table)
FALSE %in% colSums(bin_fxn.cov_table[,-c(1:2)])>0 # T/F are there any genes that have 0 counts in dataframe

bin_counts_t<-as.data.frame(t(bin_fxn.cov_table[,-c(1:2)]))
dim(bin_counts_t)
class(bin_counts_t$SSW.4.13.22.10m.bin.8)
#bin_counts_t <- as.data.frame(sapply(bin_counts_t, as.numeric)) # convert integer to numeric across df
class(bin_counts_t$SSW.12.22.21.10m.bin.17) # sanity check

#descdist(bin_counts_t$SSW.4.13.22.10m.bin.8, discrete = TRUE)
#descdist(bin_counts_t$SSW.12.22.21.10m.bin.17, discrete = TRUE)
#descdist(bin_counts_t$SSW.8.24.21.0m, discrete = TRUE)

# check distribution of gene # x gene counts with two samples
ggplot(bin_counts_t) +
  geom_histogram(aes(x = SSW.4.13.22.10m.bin.8), stat = "bin", bins = 200) +
  xlab("Gene coverage") +
  ylab("Number of genes")

ggplot(bin_counts_t) +
  geom_histogram(aes(x = SSW.8.24.21.0m.bin.15), stat = "bin", bins = 200) +
  xlab("Gene coverage") +
  ylab("Number of genes")
## data is not normally distributed

# histogram of gene coverage (excluding genes with 0 coverage)
hist(bin_counts_t$SSW.12.22.21.10m.bin.17[bin_counts_t$SSW.12.22.21.10m.bin.17>0], col="blue")
# visualize Q-Q plot for
qqnorm(bin_counts_t$SSW.12.22.21.10m.bin.17[bin_counts_t$SSW.12.22.21.10m.bin.17>0], pch = 1, frame = FALSE)
qqline(bin_counts_t$SSW.12.22.21.10m.bin.17[bin_counts_t$SSW.12.22.21.10m.bin.17>0], col = "red", lwd = 2)

# histogram of gene coverage  (excluding genes with 0 coverage)
hist(bin_counts_t$SSW.4.13.22.10m.bin.8[bin_counts_t$SSW.4.13.22.10m.bin.8>0], col="blue")
# visualize Q-Q plot
qqnorm(bin_counts_t$SSW.4.13.22.10m.bin.8[bin_counts_t$SSW.4.13.22.10m.bin.8>0], pch = 1, frame = FALSE)
qqline(bin_counts_t$SSW.4.13.22.10m.bin.8[bin_counts_t$SSW.4.13.22.10m.bin.8>0], col = "red", lwd = 2)

save.image("data/Metagenomes/Analysis/mgm_MAG_analysis.Rdata")

#### Sum Coverage by KO ID in MAGs Before Transformations ####
bin_fxns.cov[1:4,]
bin_fxns.cov_noNA<-as.data.frame(bin_fxns.cov[!is.na(bin_fxns.cov$KO_ID),]) # drop KOs with NA as assignment

bin.ko.cov.sum_table<-as.data.frame(dcast(bin_fxns.cov_noNA, SampleID+Bin_ID~KO_ID, value.var="CovPerGene", fun.aggregate=sum)) ###
bin.ko.cov.sum_table[1:4,1:4]
rownames(bin.ko.cov.sum_table)<-bin.ko.cov.sum_table$Bin_ID
bin.ko.cov.sum_table[1:4,1:5]

# check rownames of summed CLR transformed feature coverage data & metadata
bin.ko.cov.sum_table$SampleID %in% rownames(bin_meta_scaled)

bin.ko.cov.sum_table[1:4,1:4] # contains the sum of coverages per gene per KO -- featureCounts was normalized by gene length across samples first

#### Prepare Bin Feature Count Data for Normalization w/ DESeq2 ####
# make sure count data & bin_meta are in the same order
head(bin.ko.cov.sum_table)

# scale function coverage by 100 so that you can round to create DESeq2 object
## DESeq2 requires whole integers
bin_scale.sum.cov_table<-bin.ko.cov.sum_table
bin_fxn.cov_table[,3:7]*100 # sanity check that this method will work; first 2 columns are SampleID & Bin_ID
bin_scale.sum.cov_table[,-c(1:2)]<-bin_scale.sum.cov_table[,-c(1:2)]*100

bin_scale.sum.cov_table[,1:10]

rownames(bin_scale.sum.cov_table) # sanity check
rownames(bin_meta)

rownames(bin_meta) %in% rownames(bin_scale.sum.cov_table)

bin_meta=bin_meta[rownames(bin_scale.sum.cov_table),] ## reorder bin_meta to have same rows as original OTU table

# convert all NAs to 0; will take a while
bin_scale.sum.cov_table[1:6,1:6]

#convert data frame to numeric
#bin_scale.sum.cov_table[,-c(1:2)] <- as.data.frame(sapply(bin_scale.sum.cov_table[,-c(1:2)], as.numeric))

scale.sum.cov.bin_matrix<-data.matrix(bin_scale.sum.cov_table[,!names(bin_scale.sum.cov_table) %in% c("SampleID","Bin_ID")]) # convert count table into matrix
scale.sum.cov.bin_matrix[1:4,1:4]
scale.sum.cov.bin_mat.pseudo<-scale.sum.cov.bin_matrix # create matrix with pseudocounts, replacing 0 as 1
scale.sum.cov.bin_mat.pseudo[scale.sum.cov.bin_mat.pseudo == 0] <- 1 # convert 0 to 1 as pseudocount

scale.sum.cov.bin_mat.pseudo[1:4,1:4]

scale.sum.cov.bin_mat.pseudo2<-t(scale.sum.cov.bin_mat.pseudo) # transpose matrix so KO_IDs are rows, samples are columns
# ^ will be used in DESeq2 functions
rownames(bin_meta) %in% colnames(scale.sum.cov.bin_mat.pseudo2) # check if rownames in metadata (SampleID) match column names in count data
dim(bin_meta)
dim(scale.sum.cov.bin_mat.pseudo2)

# create the DESeq DataSet object for DGE analysis
# DESeq2 needs whole # data, so need raw read counts, NOT coverage for these tables...questionable
# bin_dds has the scaled coverage that was calculated by dividing reads from featureCounts by gene lengths
bin_dds<-DESeqDataSetFromMatrix(countData=round(scale.sum.cov.bin_mat.pseudo2), colData = bin_meta, design = ~ 1)

# design = ~ 1 means no design
head(counts(bin_dds))
colSums(counts(bin_dds)) %>% barplot

# add rownames from KO IDs to DESeq Data set object
BinfeatureData <- data.frame(KO_ID=rownames(bin_dds)) # create object of KO IDs
mcols(bin_dds) <- DataFrame(mcols(bin_dds), BinfeatureData) # add KO IDs to rownames of DeSEQ dataset object
mcols(bin_dds)

# drop KO IDs with lowest counts
keep <- rowSums(counts(bin_dds)) >= 10
bin_dds <- bin_dds[keep,]

counts(bin_dds)

# Estimate size factor - aka normalization factor, divide all read counts by each size factor per sample
bin_dds <- estimateSizeFactors(bin_dds,type="ratio")
## To calculate size factor in DESeq2:
# calculates geometric mean of each gene in each sample and uses this as a pseudoreference
# calculates ratio of each sample by dividing each gene count by it's pseudoreference in each sample
# The median value of all ratios for a given sample is taken as the normalization factor (size factor)
# The differentially expressed genes should not affect the median value
# median of ratios method makes the assumption that not ALL genes are differentially expressed; therefore, the normalization factors should account for sequencing depth and RNA composition of the sample
## (large outlier genes will not represent the median ratio values)
# more here: https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html

sizeFactors(bin_dds)

plot(sizeFactors(bin_dds), colSums(counts(bin_dds)))

# Does sequencing depth influence normalization?
par(mfrow=c(1,2)) # to plot the two box plots next to each other
boxplot(log2(counts(bin_dds)), notch=TRUE,
        main = "Non-normalized read counts\n(log-transformed)",
        ylab="read counts")
boxplot(log2(counts(bin_dds, normalize= TRUE)), notch=TRUE,
        main = "Size-factor-normalized read counts\n(log-transformed)",
        ylab="read counts")
dev.off()

#save.image("data/Metagenomes/Analysis/mgm_MAG_analysis.Rdata")

#### Median-Ratio Normalized - Gene Counts in MAGs ####
#
sizeFactors(bin_dds) # will be used to normalize counts
## to normalize counts, we divide each raw count value in a given sample by that sample’s normalization factor (aka size factor) to generate normalized count values.
### This is performed for all count values (every gene in every sample)
bin_fxn_mr <- counts(bin_dds, normalized=TRUE) ## median-ratio transformation is how DESeq2 Normalizes counts!
#write.table(bin_fxn_counts.norm, file="data/Metagenomes/Analysis/bin_NoBins_MedianRatio_GeneCounts_2.27.23.txt", sep="\t", quote=F, col.names=T)

bin.mr<-as.data.frame(t(bin_fxn_mr))
bin.mr$SampleID<-rownames(bin.mr)
bin.mr[1:4,1:4]

#### Variance Stabilizing Transformation - Gene Counts in MAGs ####

# you should be able to use matrix or DESeq2 object for this next function, but matrix was not working?
# variance stabilizing transformation
bin_fxn_vst <- varianceStabilizingTransformation(bin_dds, blind = TRUE, fitType = "parametric")
assay(bin_fxn_vst) #see output of VST

bin.vst<-assay(bin_fxn_vst)
#total_fxn_vst<-colSums(bin_fxn_vst)

#### Centered Log Ratio Transformation - Gene Counts in MAGs ####
bin.ko.cov.sum_table[1:4,1:4]
# ^ table contains gene coverage, Sample IDs as rows & genes as columns

# df must have rownames are SampleIDs, columns are ASV IDs for vegan functions below\
bin.clr<-decostand(bin.ko.cov.sum_table[,-c(1:2)],method = "clr", pseudocount = 1) #CLR transformation
bin.clr[1:4,1:4]

# check rownames of CLR transformed ASV data & metadata
rownames(bin.clr) %in% rownames(bin_meta)

#### Robust Centered Log Ratio Transformation - Gene in MAGs ####
bin.ko.cov.sum_table[1:4,1:4]
# ^ table contains gene coverage, Sample IDs as rows & genes as columns

# df must have rownames are SampleIDs, columns are ASV IDs for vegan functions below\
bin.Rclr<-decostand(bin.ko.cov.sum_table[,-c(1:2)],method = "rclr") #CLR transformation
bin.Rclr[1:4,1:4]
# NOTE: Robust CLR just excludes 0s and performs CLR transformation without pseudocount
# robust clr ("rclr") is similar to regular clr (see above) but allows data that contains zeroes.
# This method does not use pseudocounts, unlike the standard clr. Robust clr divides the values by geometric mean of the observed features; zero values are kept as zeroes, and not taken into account.
#In high dimensional data, the geometric mean of rclr is a good approximation of the true geometric mean

# check rownames of CLR transformed ASV data & metadata
rownames(bin.Rclr) %in% rownames(bin_meta)

#### Copies Per Million Transformation - Gene Counts in MAGs ####
bin.ko.cov.sum_table[1:4,1:4] # sanity check
bin_fxn.sum.cov_t.table<-as.data.frame(t(as.data.frame(bin.ko.cov.sum_table[,-c(1:2)])))
bin_fxn.sum.cov_t.table[1:4,1:4]

bin_fxn_cpm<-(bin_fxn.sum.cov_t.table/colSums(bin_fxn.sum.cov_t.table))*10^6
bin_fxn_cpm[1:4,1:4]
colSums(bin_fxn_cpm)

#write.table(bin_fxn_cpm, file="data/Metagenomes/Analysis/bin_NoBins_CopiesPerMillion_GeneCounts_2.27.23.txt", sep="\t", quote=F, col.names=T)

#### Create Presence/Absence Table of Functions in MAGs ####
bin.ko.cov.sum_table[,-c(1:2)][1:4,1:4]

bin_fxn.binary<-counts_to_binary(bin.ko.cov.sum_table[,-c(1:2)]) # custom function to convert all counts to binary (presence/absence)
# sanity check that function worked below
bin_fxn.binary[1:5,1:5]
bin.ko.cov.sum_table[,-c(1:2)][1:5,1:5]

#### Create CLR-transformed Coveraged Table w/ NAs for Absent Functions ####

# are NA table and bin.clr in the same order?
rownames(bin_fxn.binary)
rownames(bin.clr)

colnames(bin_fxn.binary)
colnames(bin.clr)
colnames(bin_fxn.binary)[which(colnames(bin_fxn.binary) %in% colnames(bin.clr) == FALSE)] # bin_fxn.binary does not have SampleID column
colnames(bin.clr)[which(colnames(bin.clr) %in% colnames(bin_fxn.binary) == FALSE)] # bin.clr does not have SampleID column

identical(bin.clr,bin_fxn.binary) # SampleID columns are in different places in these two data frames...

# then create logical table saying which functions are NA in this ko.cov.na table AND have a negative coverage value in bin.clr
# TRUE means they are absent, FALSE means they are present
NA.fxns <- (bin_fxn.binary==0)

# create data frame that will be CLR transformed sum coverages + NA values
bin.clr.na<-bin.clr
bin.clr.na[NA.fxns] <- NA
#bin.clr.na[NA.fxns == TRUE]<- NA

# did this work?
bin.clr[1:4,1:4]
bin.clr.na[1:4,1:4]
bin_fxn.binary[1:4,1:4]

# create sample ID column
bin.clr.na$Bin_ID<-rownames(bin.clr.na)
# * use bin.clr.na for heatmaps of functions and coverage!
# to get rid of sampleID column later for bin.clr.na, use the following code
## bin.clr.na[,!names(bin.clr.na) %in% c("SampleID")]



### Pull out traits of interest ####
# create unique list of KO ID and functions
# check for duplicates to make sure each KO_ID has a unique function assignment
head(bin_fxns.cov)
NA %in% bin_fxns.cov$CovPerGene # just to ensure there are no NAs for genes in this df

bin.ko_fxns1<-unique(bin_fxns.cov[,1:3]) # first subset out data based on unique KO functions

n_occur <- data.frame(table(bin.ko_fxns1$KO_ID)) # see how many duplicates there are of KO IDs, compare duplicates
n_occur[n_occur$Freq > 1,] # what traits appear more than once?

bin.ko_ID<-unique(data.frame(KO_ID=bin.ko_fxns1$KO_ID)) # get a list of unique KO IDs in data

bin.ko_fxns<-as.data.frame(bin.ko_fxns1[!is.na(bin.ko_fxns1$KO_ID),])  # use unique KO ID list to subset out KO function data
head(bin.ko_fxns)

## pull out functions of interest
carb.fxns.bins<-bin.ko_fxns[which(bin.ko_fxns$KO_ID %in% carb.kegg$KO_ID),]
All_GOI.fxns.bins<-bin.ko_fxns[which(bin.ko_fxns$KO_ID %in% all_goi.kegg$KO_ID),]
sulfur.fxns.bins<-bin.ko_fxns[which(bin.ko_fxns$KO_ID %in% sulf.kegg$KO_ID),]
nitro.fxns.bins<-bin.ko_fxns[which(bin.ko_fxns$KO_ID %in% nitro.kegg$KO_ID),]
osmo.fxns.bins<-bin.ko_fxns[which(bin.ko_fxns$KO_ID %in% osmo.kegg$KO_ID),]
selen.fxns.bins<-bin.ko_fxns[which(bin.ko_fxns$KO_ID %in% selen.kegg$KO_ID),]
arsen.fxns.bins<-bin.ko_fxns[which(bin.ko_fxns$KO_ID %in% arsen.kegg$KO_ID),]
HS.fxns.bins<-bin.ko_fxns[which(bin.ko_fxns$KO_ID %in% heatshock.kegg$KO_ID),]
metal.fxns.bins<-bin.ko_fxns[which(bin.ko_fxns$KO_ID %in% metal.re.kegg$KO_ID),]
photo.fxn.bins<-bin.ko_fxns[which(bin.ko_fxns$KO_ID %in% photo.kegg$KO_ID),]
aero.fxn.bins<-bin.ko_fxns[which(bin.ko_fxns$KO_ID %in% aero.kegg$KO_ID),]

### Merge Metadata & Contig Count Data & Taxa Data Together ####
bin.clr[1:4,1:4]
head(bin_meta)
head(bin.ko_fxns)
head(mapped_all)

# melt data with normalized feature counts to merge with all traits & traits of interest
bin.clr$Bin_ID<-rownames(bin.clr)
bin_clr_melt<-melt(bin.clr, by="Bin_ID")
head(bin_clr_melt)
colnames(bin_clr_melt)[which(names(bin_clr_melt) == "variable")] <- "KO_ID"
colnames(bin_clr_melt)[which(names(bin_clr_melt) == "value")] <- "SumCovPerKO"
bin_clr_melt[1:4,]

bin.clr.all<-as.data.frame(merge(bin.ko_fxns, bin_clr_melt, by=c("KO_ID"),allow.cartesian = TRUE))
head(bin.clr.all)

NA %in% bin.clr.all$KO_ID

bin.clr.all<-bin.clr.all[!is.na(bin.clr.all$KO_ID),] # only looking at functions we have KO IDs for
head(bin.clr.all)
#NA %in% bin.clr.all

bin_clr_melt[1:4,] #sanity check

bin.clr.sulf<-merge(sulfur.fxns.bins, bin_clr_melt, by=c("KO_ID"),allow.cartesian = TRUE)
head(bin.clr.sulf)

bin.clr.ars<-merge(arsen.fxns.bins, bin_clr_melt, by=c("KO_ID"),allow.cartesian = TRUE)
head(bin.clr.ars)

bin.clr.osmo<-merge(osmo.fxns.bins, bin_clr_melt, by=c("KO_ID"),allow.cartesian = TRUE)
head(bin.clr.osmo)

bin.clr.HS<-merge(HS.fxns.bins, bin_clr_melt, by=c("KO_ID"),allow.cartesian = TRUE)
head(bin.clr.HS)

bin.clr.metal<-merge(metal.fxns.bins, bin_clr_melt, by=c("KO_ID"),allow.cartesian = TRUE)
head(bin.clr.metal)

# merge mapped read counts & bin taxonomy data to metadata
mapped_meta<-merge(mapped_all, bin_meta, by=c("SampleID"))
head(mapped_meta)

#save.image("data/Metagenomes/Analysis/mgm_MAG_analysis.Rdata")

### Export Global Env for Other Scripts ####
save.image("data/Metagenomes/Analysis/mgm_MAG_analysis.Rdata")
# ^ includes all data combined in object bac.dat.all, ASV table (samples are rows, ASVs are columns), mgm_meta, and an ASV count table (where ASVs are rows, not columns)
# Version Information
sessionInfo()