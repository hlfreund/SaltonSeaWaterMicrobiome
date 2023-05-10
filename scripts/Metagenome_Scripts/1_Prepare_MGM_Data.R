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

dt.to.na <- function(dataFrame){
  new_m <- matrix(nrow=dim(dataFrame)[1],ncol = dim(dataFrame)[2]) # create new matrix w/ same rows and cols as input dataframe
  ## dim(df)[1] gives you first dimensions (x aka rows), dim(df)[2] gives you second dimensions (y aka columns)

  for( currentRow in 1:nrow(dataFrame)){ # for every row
    for( currentCol in 1:ncol(dataFrame)){ # for every column

      if ( is.na(dataFrame[currentRow, currentCol]) & is.numeric(dataFrame[currentRow, currentCol])){ # if both row and col (specifies each cell) are NA, change val to 0
        new_m[currentRow, currentCol] = 0
        # is.numeric(df[currentRow,currentCol]) is to confirm each cell contains a numeric element
      } else if( is.numeric(dataFrame[currentRow, currentCol]) & dataFrame[currentRow, currentCol] > 0){ # if both row and col (specifies each cell) are > 0, change val to 1
        next
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

##save.image("data/Metagenomes/Analysis/mgm_analysis.Rdata") # save global env to Rdata file

## Notes:
# code & info came from :
## https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#quick-start
## https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/07_practical.pdf
## https://www.reneshbedre.com/blog/deseq2.html
## https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html

#### Import MGM Read Counts & Taxonomic Annotation Data ####
mgm_fxns.cov<-fread(file = 'data/Metagenomes/Analysis/SSW_Samples_NoBins_Gene_Coverages_5.3.23.txt', sep='\t',header = TRUE)
dim(mgm_fxns.cov)
head(mgm_fxns.cov)

# count occurrences of all traits across all mgms
n_occur <- data.frame(table(mgm_fxns.cov$KO_Function))
n_occur[n_occur$Freq > 1,]

# update SampleID names so they match w/ metadata
mgm_fxns.cov$SampleID<-gsub("_",".", mgm_fxns.cov$SampleID)
head(mgm_fxns.cov)
mgm_fxns.cov$SampleID<-gsub("m\\..*","m",mgm_fxns.cov$SampleID)
mgm_fxns.cov$PlotID<-gsub("^SSW.","",mgm_fxns.cov$SampleID)
head(mgm_fxns.cov)

# Divide gene counts by gene length to account for sample differences in assembly
## we do this because we did not co-assemble contigs, so each KO assignments across samples may not come from genes with the same length
## need to account for gene length here first (since KO assignments can come from genes of different lengths)
mgm_fxns.cov$CovPerGene<-mgm_fxns.cov$ReadsPerGene/mgm_fxns.cov$GeneLength
head(mgm_fxns.cov)

# create list of GeneIDs & KO IDs
gene_KOs<-unique(data.frame(Gene_ID=mgm_fxns.cov$Gene_ID, KO_ID=mgm_fxns.cov$KO_ID))
dim(gene_KOs)

# create Sample ID x Gene ID count table, using reads per gene that were divided by gene length
mgm_fxn.cov_table<-dcast(mgm_fxns.cov, SampleID~Gene_ID, value.var="CovPerGene")
mgm_fxn.cov_table[,1:4] # sanity check
rownames(mgm_fxn.cov_table)<-mgm_fxn.cov_table$SampleID

# convert all NAs to 0; will take a while
mgm_fxn.cov_table[,-1][is.na(mgm_fxn.cov_table[,-1])] <- 0
mgm_fxn.cov_table[1:6,1:6] # sanity check

mgm.fxn.nolows<-mgm_fxn.cov_table[,which(colSums(mgm_fxn.cov_table[,-1])>=5)] # remove functions with less than 15 total counts across mgms
#mgm_fxn.binary<-counts_to_binary(mgm_fxn.cov_table[,-1]) # custom function to convert all counts to binary (presence/absence)
mgm.fxn.nolows[1:4,1:4] # sanity check

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

# Divide gene counts by gene length to account for sample differences in assembly
## we do this because we did not co-assemble contigs, so each KO assignments across samples may not come from genes with the same length
## need to account for gene length here first (since KO assignments can come from genes of different lengths)
bin_fxns.cov$CovPerGene<-bin_fxns.cov$ReadsPerGene/bin_fxns.cov$GeneLength
head(bin_fxns.cov)

# create list of GeneIDs & KO IDs
bin.gene_KOs<-unique(data.frame(Gene_ID=bin_fxns.cov$Gene_ID, KO_ID=bin_fxns.cov$KO_ID))
dim(bin.gene_KOs)

bin_fxn.counts_table<-dcast(bin_fxns.cov, SampleID+Bin_ID~Gene_ID, value.var="CovPerGene")
bin_fxn.counts_table[1:4,1:4] # sanity check
rownames(bin_fxn.counts_table)<-bin_fxn.counts_table$Bin_ID
bin_fxn.counts_table[1:4,1:4] # sanity check

# Remove unwanted samples
#remove_samples<-c("SS.OV.10m.seawater.0621", "SS.OV.2m.seawater.0621", "SS.OV.5m.seawater.0621")
#bac.ASV_counts<-bac.ASV_counts[,!(colnames(bac.ASV_counts) %in% remove_samples)]
#colnames(bac.ASV_counts)
#dim(bac.ASV_counts)

## Import MGM taxonomic data
mag_tax<-fread(file = 'data/Metagenomes/Analysis/SSW_MAGs_TaxoAnnotation_2.7.23.tsv', sep='\t',header = TRUE,fill=TRUE)
head(mag_tax)
mag_tax$Bin_ID<-gsub("_",".", mag_tax$Bin_ID)
mag_tax$Bin_ID<-gsub("m\\..*\\.bin","m.bin",mag_tax$Bin_ID) # removes rep number between m and .bin
mag_tax$PlotBin<-gsub("^SSW.","",mag_tax$Bin_ID)

head(mag_tax)

mag_tax[mag_tax==""]<-"Unknown" # replace blank cells with Unknown label
mag_tax[is.na(mag_tax)]<- "Unknown" # turn all NAs into "Unkowns"
mag_tax$Species<-gsub("Unknown", "unknown", mag_tax$Species) # change uppercase Unknown to lowercase unknown for unknown species classification
head(mag_tax)

## Import Contig & Bin Mapped Read Counts
mapped_reads<-as.data.frame(read_xlsx("data/Metagenomes/Analysis/Total_Contig_Bin_Reads.xlsx", sheet="Total_Bin_Reads"))
head(mapped_reads)
mapped_reads$PlotBin<-gsub("^SSW.","",mapped_reads$Bin_ID)
mapped_reads$PlotID<-gsub("^SSW.","",mapped_reads$SampleID)

head(mapped_reads)
mapped_all<-merge(mapped_reads,mag_tax,by=c("PlotBin","Bin_ID"))
head(mapped_all)
mapped_all$RelAb_Map<-mapped_all$Bin_Mapped_Reads/mapped_all$Total_Mapped_Reads
head(mapped_all)

#save.image("data/Metagenomes/Analysis/mgm_analysis.Rdata")

#### Update Metadata ####
# upload UNSCALED geochem data from Lyons lab
# ^^ scale env variable data in respective scripts

mgm_meta<-as.data.frame(read_xlsx("data/Metagenomes/Analysis/SSW_Lyons_Aronson_Metadata_MGM_All.xlsx", sheet="Metagenomes_Metadata"))
head(mgm_meta)
mgm_meta$SampleID<-gsub("m\\..*","m",mgm_meta$SampleID) # remove rep # after m
mgm_meta$PlotID<-gsub("^SSW.","",mgm_meta$SampleID)
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
# #save.image("data/SSW_analysis.Rdata")

#### Scale Chem Data in Metadata ####
head(mgm_meta)
#meta_scaled<-subset(mgm_meta, select=-c(Salinity_ppt)) # drop salinity, for now
meta_scaled<-mgm_meta

head(meta_scaled)
meta_scaled[,8:17]<-as.data.frame(scale(meta_scaled[,8:17], center=TRUE, scale=TRUE)) #centering before scaling
meta_scaled

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

# scale function coverage by 10 so that you can round to create DESeq2 object
## DESeq2 requires whole integers
mgm_fxns.cov$ScaleCov<-mgm_fxns.cov$CovPerGene*10 # scale coverages by 10

mgm_scale.cov_table<-dcast(mgm_fxns.cov, SampleID~Gene_ID, value.var="ScaleCov")

rownames(mgm_scale.cov_table)<-mgm_scale.cov_table$SampleID
rownames(mgm_scale.cov_table) # sanity check
rownames(mgm_meta)

mgm_meta=mgm_meta[rownames(mgm_scale.cov_table),] ## reorder mgm_meta to have same rows as original OTU table

# convert all NAs to 0; will take a while
mgm_scale.cov_table[1:6,1:6]
mgm_scale.cov_table[,-1][is.na(mgm_scale.cov_table[,-1])] <- 0
mgm_scale.cov_table[1:6,1:6]

#convert data frame to numeric
#mgm_scale.cov_table[,-1] <- as.data.frame(sapply(mgm_scale.cov_table[,-1], as.numeric))

scale_cov_matrix<-as.matrix(mgm_scale.cov_table[,-1]) # convert count table into matrix
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

scale_cov_matrix<-as.matrix(mgm_scale.cov_table[,-1]) # convert count table into matrix
scale_cov_matrix2<-t(scale_cov_matrix) # transpose matrix so KO_IDs are rows, samples are columns
# ^ will be used in DESeq2 functions
rownames(mgm_meta) %in% colnames(scale_cov_matrix2) # check if rownames in metadata (SampleID) match column names in count data
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
"NA" %in% mgm_fxns.cov$CovPerGene # just to ensure there are no NAs for genes in this df

ko_fxns1<-mgm_fxns.cov[which(mgm_fxns.cov$KO_Function %in% unique(mgm_fxns.cov$KO_Function)),] # first subset out data based on unique KO functions
ko_fxns2<-subset(ko_fxns1, select=c("Gene_ID","KO_ID","KO_Function"))

n_occur <- data.frame(table(ko_fxns2$KO_ID)) # see how many duplicates there are of KO IDs, compare duplicates
n_occur[n_occur$Freq > 1,] # what traits appear more than once?

ko_ID<-unique(data.frame(KO_ID=ko_fxns2$KO_ID)) # get a list of unique KO IDs in data

ko_fxns<-unique(ko_fxns2[which(ko_fxns2$KO_ID %in% ko_ID$KO_ID),]) # use unique KO ID list to subset out KO function data
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

mgm.clr.all<-merge(ko_fxns, mgm_clr_melt, by=c("Gene_ID"),allow.cartesian = TRUE)
head(mgm.clr.all)
mgm.clr.all<-subset(mgm.clr.all, KO_Function!="NA") # only looking at functions we have KO IDs for
head(mgm.clr.all)

mgm.clr.sulf<-merge(sulfur.fxns, mgm_clr_melt, by=c("Gene_ID"),allow.cartesian = TRUE)
head(mgm.clr.sulf)

mgm.clr.ars<-merge(arsenic.fxns, mgm_clr_melt, by=c("Gene_ID"),allow.cartesian = TRUE)
head(mgm.clr.ars)

# merge mapped read counts & mgm taxonomy data to metadata
mapped_meta<-merge(mapped_all, mgm_meta, by=c("SampleID"))
head(mapped_meta)

#save.image("data/Metagenomes/Analysis/mgm_analysis.Rdata")

#### Clustering by Features Across Samples ####

# using CLR data first
mgm.clr[1:4,1:4]

# calculate our Euclidean distance matrix using CLR data
mgm.euc_dist1 <- dist(mgm.clr, method = "euclidean")

# creating our hierarcical clustering dendrogram
mgm.euc_clust1 <- hclust(mgm.euc_dist1, method="ward.D2")

# let's make it a little nicer...
mgm.euc_dend <- as.dendrogram(mgm.euc_clust1, hang=0.2)
mgm.dend_cols <- as.character(mgm_meta$SampDate_Color[order.dendrogram(mgm.euc_dend)])
labels_colors(mgm.euc_dend) <- mgm.dend_cols

plot(mgm.euc_dend, ylab="CLR Euclidean Distance",cex = 0.5) + title(main = "Functional Trait Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 3)
legend("topright",legend = c("June 2021","August 2021","December 2021","April 2022"),cex=.8,col = c( "#36ab57","#ff6f00","#26547c","#32cbff"),pch = 15, bty = "n")

# Control is dark blue ("#218380"), #Alternaria is light blue ("#73d2de")
dev.off()

#### Functional Beta Diversity - MRT data ####
# MR = median-ratio transformation
mgm.mr[1:4,1:4] # sample IDs are rows, genes are columns
mgm_fxn.counts_table[1:4,1:4] # sanity check

# check rownames of MR & Mr transformed feature count data & metadata
rownames(mgm.mr) %in% rownames(meta_scaled)

## PCOA with VST transformed data first
# calculate our Euclidean distance matrix using VST data
mgm.euc_dist.mr <- dist(mgm.mr, method = "euclidean")

# creating our hierarcical clustering dendrogram
mgm.euc.mr_clust <- hclust(mgm.euc_dist.mr, method="ward.D2")

# let's make it a little nicer...
mgm.euc.mr_dend <- as.dendrogram(mgm.euc.mr_clust, hang=0.2)
mgm.dend_cols <- as.character(mgm_meta$SampDate_Color[order.dendrogram(mgm.euc.mr_dend)])
labels_colors(mgm.euc.mr_dend) <- mgm.dend_cols

plot(mgm.euc.mr_dend, ylab="Median-Ratio, Euclidean Distance",cex = 0.5) + title(main = "Bacteria/Archaea Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 3)
legend("topright",legend = c("June 2021","August 2021","December 2021","April 2022"),cex=.8,col = c( "#36ab57","#ff6f00","#26547c","#32cbff"),pch = 15, bty = "n")
# Control is dark blue ("#218380"), #Alternaria is light blue ("#73d2de")
dev.off()

# let's use our Euclidean distance matrix from before
mgm.pcoa.mr <- pcoa(mgm.euc_dist.mr) # pcoa of euclidean distance matrix = PCA of euclidean distance matrix
##save.image("data/ssw_mr.euc.dist1_3.7.23.Rdata")

# The proportion of variances explained is in its element values$Relative_eig
mgm.pcoa.mr$values

# extract principal coordinates
mgm.pcoa.mr.vectors<-data.frame(mgm.pcoa.mr$vectors)
mgm.pcoa.mr.vectors$SampleID<-rownames(mgm.pcoa.mr$vectors)

# merge pcoa coordinates w/ metadata
mgm.pcoa.mr.meta<-merge(mgm.pcoa.mr.vectors, mgm_meta, by.x="SampleID", by.y="SampleID")
mgm.pcoa.mr.meta$SampleMonth
mgm.pcoa.mr.meta$SampDate

head(mgm.pcoa.mr.meta)

head(mgm.pcoa.mr$values) # pull out Relative (Relative_eig) variation % to add to axes labels

# create PCoA ggplot fig
pcoa3<-ggplot(mgm.pcoa.mr.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(SampDate)), size=4)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Seawater",subtitle="Using Median-Ratio Transformed Feature Data",xlab="PC1 [41.14%]", ylab="PC2 [9.04%]",color="Sample Type")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Type",values=unique(mgm.pcoa.mr.meta$SampDate_Color[order(mgm.pcoa.mr.meta$SampDate)]),labels=c("June.2021"="June 2021","August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("PC1 [33.04%]") + ylab("PC2 [29.24%]")

ggsave(pcoa3,filename = "figures/MGM_Figs/SSW_MGM_pcoa_MR_sampdate.png", width=12, height=10, dpi=600)

# sample month shape, depth color
pcoa4<-ggplot(mgm.pcoa.mr.meta, aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampleMonth), size=5)+theme_bw()+
  labs(title="PCoA: Metagenome Functions in Salton Seawater",subtitle="Using Median-Ratio Transformed Feature Data",xlab="PC1", ylab="PC2",color="Depth (m)")+
  theme_classic()+ theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),axis.text = element_text(size=12),axis.text.x = element_text(vjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17))+
  scale_color_continuous(low="blue3",high="red",trans = 'reverse') + scale_shape_discrete(labels=c("June 2021","August 2021","December 2021","April 2022"),name="Sample Date") +
  xlab("PC1 [33.04%]") + ylab("PC2 [29.24%]")

ggsave(pcoa4,filename = "figures/MGM_Figs/SSW_MGM_pcoa_MR.traits_depth.png", width=12, height=10, dpi=600)

## betadisper to look at within group variance

# first by sampling date
mgm.disper1<-betadisper(mgm.euc_dist.mr, mgm_meta$SampDate)
mgm.disper1

permutest(mgm.disper1, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
#Pairwise comparisons:
#  (Observed p-value below diagonal, permuted p-value above diagonal)
#               June.2021 August.2021 December.2021 April.2022
#June.2021                   0.69800       0.26000      0.221
#August.2021     0.62084                   0.36600      0.425
#December.2021   0.25263     0.33869                    0.747
#April.2022      0.23545     0.37980       0.69090

anova(mgm.disper1) # p = 0.5253 --> accept the Null H, spatial medians are NOT significantly difference across sample dates

TukeyHSD(mgm.disper1) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other

#                             diff        lwr       upr     p adj
#August.2021-June.2021      33.096688 -101.9651 168.15849 0.8474848
#December.2021-June.2021   -18.757013 -153.8188 116.30479 0.9655375
#April.2022-June.2021      -13.429156 -148.4910 121.63265 0.9866629
#December.2021-August.2021 -51.853701 -172.6567  68.94925 0.5260566
#April.2022-August.2021    -46.525845 -167.3288  74.27711 0.6046202
#April.2022-December.2021    5.327857 -115.4751 126.13081 0.9987814

# Visualize dispersions
png('figures/MGM_Figs/SSW_MGM_pcoa_MR_betadispersion_sampledate.png',width = 700, height = 600, res=100)
plot(mgm.disper1,main = "Centroids and Dispersion (Median-Ratio Data)", col=colorset1$SampDate_Color)
dev.off()

png('figures/MGM_Figs/SSW_MGM_boxplot_MR_centroid_distance_sampledate.png',width = 700, height = 600, res=100)
boxplot(mgm.disper1,xlab="Sample Collection Date", main = "Distance to Centroid by Category (Median-Ratio Data)", sub="Euclidean Distance of Median-Ratio Transformed Data", col=colorset1$SampDate_Color)
dev.off()

# What about between sampling depths?
mgm.disper2<-betadisper(mgm.euc_dist.mr, mgm_meta$Depth_m)
mgm.disper2

permutest(mgm.disper2, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons

anova(mgm.disper2) # p = 0.3914 --> accept the Null H, spatial medians are NOT significantly difference across sample dates

TukeyHSD(mgm.disper2) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other
#         diff       lwr      upr     p adj
#5-0  67.213614  -84.49614 218.9234 0.4511753
#10-0 69.573085  -82.13667 221.2828 0.4287061
#10-5  2.359471 -138.09647 142.8154 0.9987307

colfunc <- colorRampPalette(c("red", "blue"))
colfunc(3)

# Visualize dispersions
png('figures/MGM_Figs/ssw_mgm_pcoa_MR_betadispersion_depth.png',width = 700, height = 600, res=100)
plot(mgm.disper2,main = "Centroids and Dispersion (Median-Ratio Data)", col=colfunc(3))
dev.off()

png('figures/MGM_Figs/ssw_mgm_boxplot_MR_centroid_distance_depth.png',width = 700, height = 600, res=100)
boxplot(mgm.disper2,xlab="Sample Collection Depth", main = "Distance to Centroid by Category (Median-Ratio Data)", sub="Euclidean Distance of Median-Ratio Transformed Data", col=colfunc(3))
dev.off()


#### Functional Beta Diversity - VST data ####
mgm.vst[1:4,1:4] # sample IDs are rows, genes are columns
mgm_fxn.counts_table[1:4,1:4] # sanity check

# check rownames of VST & VST transformed feature count data & metadata
colnames(mgm.vst) %in% rownames(meta_scaled)

## PCOA with VST transformed data first
# calculate our Euclidean distance matrix using VST data
mgm.euc.vst_dist <- dist(t(mgm.vst), method = "euclidean")

# creating our hierarcical clustering dendrogram
mgm.euc.vst_clust <- hclust(mgm.euc.vst_dist, method="ward.D2")

# let's make it a little nicer...
mgm.euc.vst_dend <- as.dendrogram(mgm.euc.vst_clust, hang=0.2)
mgm.dend_cols <- as.character(mgm_meta$SampDate_Color[order.dendrogram(mgm.euc.vst_dend)])
labels_colors(mgm.euc.vst_dend) <- mgm.dend_cols

plot(mgm.euc.vst_dend, ylab="VST Euclidean Distance",cex = 0.5) + title(main = "Bacteria/Archaea Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 3)
legend("topright",legend = c("June 2021","August 2021","December 2021","April 2022"),cex=.8,col = c( "#26547c","#36ab57","#32cbff","#ff6f00"),pch = 15, bty = "n")
# Control is dark blue ("#218380"), #Alternaria is light blue ("#73d2de")
dev.off()

# let's use our euc.vstlidean distance matrix from before
mgm.pcoa.vst <- pcoa(mgm.euc.vst_dist) # pcoa of euc.vstlidean distance matrix = PCA of euc.vstlidean distance matrix
##save.image("data/ssw_vst.euc.vst.dist1_3.7.23.Rdata")

# The proportion of variances explained is in its element values$Relative_eig
mgm.pcoa.vst$values

# extract principal coordinates
mgm.pcoa.vst.vectors<-data.frame(mgm.pcoa.vst$vectors)
mgm.pcoa.vst.vectors$SampleID<-rownames(mgm.pcoa.vst$vectors)

# merge pcoa coordinates w/ metadata
mgm.pcoa.vst.meta<-merge(mgm.pcoa.vst.vectors, mgm_meta, by.x="SampleID", by.y="SampleID")
mgm.pcoa.vst.meta$SampleMonth
mgm.pcoa.vst.meta$SampDate

head(mgm.pcoa.vst.meta)

mgm.pcoa.vst$values # pull out Relative (Relative_eig) variation % to add to axes labels

# create PCoA ggplot fig
pcoa3<-ggplot(mgm.pcoa.vst.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(SampDate)), size=4)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Seawater",subtitle="Using Variance Stabilization Transformed Feature Data",xlab="PC1 [41.14%]", ylab="PC2 [9.04%]",color="Sample Type")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Type",values=unique(mgm.pcoa.vst.meta$SampDate_Color[order(mgm.pcoa.vst.meta$SampDate)]),labels=c("June.2021"="June 2021","August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("PC1 [25.06%]") + ylab("PC2 [21.98%]")

ggsave(pcoa3,filename = "figures/MGM_Figs/SSW_MGM_pcoa_VST_sampdate.png", width=12, height=10, dpi=600)

# sample month shape, depth color
pcoa4<-ggplot(mgm.pcoa.vst.meta, aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampleMonth), size=5)+theme_bw()+
  labs(title="PCoA: Metagenome Functions in Salton Seawater",subtitle="Using Variance Stabilization Transformed Feature Data",xlab="PC1", ylab="PC2",color="Depth (m)")+
  theme_classic()+ theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),axis.text = element_text(size=12),axis.text.x = element_text(vjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17))+
  scale_color_continuous(low="blue3",high="red",trans = 'reverse') + scale_shape_discrete(labels=c("June 2021","August 2021","December 2021","April 2022"),name="Sample Date") +
  xlab("PC1 [25.06%]") + ylab("PC2 [21.98%]")

ggsave(pcoa4,filename = "figures/MGM_Figs/SSW_MGM_pcoa_VST.traits_depth.png", width=12, height=10, dpi=600)

## betadisper to look at within group variance

# first by sampling date
mgm.disper3<-betadisper(mgm.euc.vst_dist, mgm_meta$SampDate)
mgm.disper3

permutest(mgm.disper3, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
#Pairwise comparisons:
#  (Observed p-value below diagonal, permuted p-value above diagonal)
#               June.2021 August.2021 December.2021 April.2022
#June.2021                  0.347000      0.083000      0.037
#August.2021    0.306448                  0.638000      0.793
#December.2021  0.072562    0.635576                    0.283
#April.2022     0.039442    0.818560      0.281216

anova(mgm.disper3) # p = 0.2825 --> accept the Null H, spatial medians are NOT significantly difference across sample dates

TukeyHSD(mgm.disper3) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other

#                             diff        lwr       upr     p adj
#August.2021-June.2021      0.42756584 -0.4879045 1.3430362 0.4622918
#December.2021-June.2021    0.59012975 -0.3253406 1.5056001 0.2312457
#April.2022-June.2021       0.35879411 -0.5566762 1.2742644 0.5923352
#December.2021-August.2021  0.16256391 -0.6562577 0.9813855 0.9097581
#April.2022-August.2021    -0.06877173 -0.8875933 0.7500498 0.9918298
#April.2022-December.2021  -0.23133564 -1.0501572 0.5874859 0.7879737

# Visualize dispersions
png('figures/MGM_Figs/SSW_MGM_pcoa_vst_betadispersion_sampledate.png',width = 700, height = 600, res=100)
plot(mgm.disper3,main = "Centroids and Dispersion (VST Data)", col=colorset1$SampDate_Color)
dev.off()

png('figures/MGM_Figs/SSW_MGM_boxplot_vst_centroid_distance_sampledate.png',width = 700, height = 600, res=100)
boxplot(mgm.disper3,xlab="Sample Collection Date", main = "Distance to Centroid by Category (VST Data)", sub="Euclidean Distance of VST Data", col=colorset1$SampDate_Color)
dev.off()

# What about between sampling depths?
mgm.disper4<-betadisper(mgm.euc.vst_dist, mgm_meta$Depth_m)
mgm.disper4

permutest(mgm.disper4, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons

anova(mgm.disper4) # p = 0.5946 --> accept the Null H, spatial medians are NOT significantly difference across sample dates

TukeyHSD(mgm.disper4) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other
#         diff       lwr      upr     p adj
#5-0  0.4686185 -1.335515 2.272752 0.7464998
#10-0 0.6553526 -1.148780 2.459486 0.5754756
#10-5 0.1867341 -1.483569 1.857037 0.9456750

colfunc <- colorRampPalette(c("red", "blue"))
colfunc(3)

# Visualize dispersions
png('figures/MGM_Figs/ssw_mgm_pcoa_VST_betadispersion_depth.png',width = 700, height = 600, res=100)
plot(mgm.disper2,main = "Centroids and Dispersion (VST Data)", col=colfunc(3))
dev.off()

png('figures/MGM_Figs/ssw_mgm_boxplot_VST_centroid_distance_depth.png',width = 700, height = 600, res=100)
boxplot(mgm.disper2,xlab="Sample Collection Depth", main = "Distance to Centroid by Category (VST Data)", sub="Euclidean Distance of VST Data", col=colfunc(3))
dev.off()

#### Functional Beta Diversity - CLR data ####
mgm.clr[1:4,1:4] # sample IDs are rows, genes are columns
mgm_fxn.counts_table[1:4,1:4] # sanity check

# check rownames of CLR & VST transformed feature count data & metadata
rownames(mgm.clr) %in% rownames(meta_scaled)

## PCOA with CLR transformed data first
# calculate our Euclidean distance matrix using CLR data
mgm.euc.clr_dist <- dist(mgm.clr, method = "euclidean")

# creating our hierarcical clustering dendrogram
mgm.euc.clr_clust <- hclust(mgm.euc.clr_dist, method="ward.D2")

# let's make it a little nicer...
mgm.euc.clr_dend <- as.dendrogram(mgm.euc.clr_clust, hang=0.2)
mgm.dend_cols <- as.character(mgm_meta$SampDate_Color[order.dendrogram(mgm.euc.clr_dend)])
labels_colors(mgm.euc.clr_dend) <- mgm.dend_cols

plot(mgm.euc.clr_dend, ylab="CLR Euclidean Distance",cex = 0.5) + title(main = "Bacteria/Archaea Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 3)
legend("topright",legend = c("June 2021","August 2021","December 2021","April 2022"),cex=.8,col = c( "#26547c","#36ab57","#32cbff","#ff6f00"),pch = 15, bty = "n")
# Control is dark blue ("#218380"), #Alternaria is light blue ("#73d2de")
dev.off()

# let's use our Euclidean distance matrix from before
mgm.pcoa.clr <- pcoa(mgm.euc.clr_dist) # pcoa of euclidean distance matrix = PCA of euclidean distance matrix
##save.image("data/ssw_clr.euc.dist1_3.7.23.Rdata")

# The proportion of variances explained is in its element values$Relative_eig
mgm.pcoa.clr$values

# extract principal coordinates
mgm.pcoa.clr.vectors<-data.frame(mgm.pcoa.clr$vectors)
mgm.pcoa.clr.vectors$SampleID<-rownames(mgm.pcoa.clr$vectors)

# merge pcoa coordinates w/ metadata
mgm.pcoa.clr.meta<-merge(mgm.pcoa.clr.vectors, mgm_meta, by.x="SampleID", by.y="SampleID")
mgm.pcoa.clr.meta$SampleMonth
mgm.pcoa.clr.meta$SampDate

head(mgm.pcoa.clr.meta)

mgm.pcoa.clr$values # pull out Relative (Relative_eig) variation % to add to axes labels

# create PCoA ggplot fig
pcoa5<-ggplot(mgm.pcoa.clr.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(SampDate)), size=4)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Seawater",subtitle="Using Centered-Log Ratio Feature Data",xlab="PC1 [41.14%]", ylab="PC2 [9.04%]",color="Sample Type")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Type",values=unique(mgm.pcoa.clr.meta$SampDate_Color[order(mgm.pcoa.clr.meta$SampDate)]),labels=c("June.2021"="June 2021","August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("PC1 [23.56%]") + ylab("PC2 [20.45%]")

ggsave(pcoa5,filename = "figures/MGM_Figs/SSW_MGM_pcoa_CLR_sampdate.png", width=12, height=10, dpi=600)

# sample month shape, depth color
pcoa6<-ggplot(mgm.pcoa.clr.meta, aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampleMonth), size=5)+theme_bw()+
  labs(title="PCoA: Metagenome Functions in Salton Seawater",subtitle="Using Centered-Log Ratio Feature Data",xlab="PC1", ylab="PC2",color="Depth (m)")+
  theme_classic()+ theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),axis.text = element_text(size=12),axis.text.x = element_text(vjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17))+
  scale_color_continuous(low="blue3",high="red",trans = 'reverse') + scale_shape_discrete(labels=c("June 2021","August 2021","December 2021","April 2022"),name="Sample Date") +
  xlab("PC1 [23.56%]") + ylab("PC2 [20.45%]")

ggsave(pcoa6,filename = "figures/MGM_Figs/SSW_MGM_pcoa_CLR.traits_depth.png", width=12, height=10, dpi=600)

## betadisper to look at within group variance

# first by sampling date
mgm.disper5<-betadisper(mgm.euc.clr_dist, mgm_meta$SampDate)
mgm.disper5

permutest(mgm.disper5, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
#Pairwise comparisons:
#  (Observed p-value below diagonal, permuted p-value above diagonal)
#               June.2021 August.2021 December.2021 April.2022
#June.2021                 0.13000000    0.03700000      0.010
#August.2021   0.12533884                0.50300000      0.189
#December.2021 0.00835232  0.56548064                    0.038
#April.2022    0.00043402  0.20006152    0.00656846

anova(mgm.disper5) # p = 0.02656 --> reject the Null H, spatial medians ARE significantly difference across sample dates

TukeyHSD(mgm.disper5) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other

#                             diff        lwr       upr     p adj
#August.2021-June.2021      3.611143 -0.4208285 7.6431151 0.0784173 .
#December.2021-June.2021    4.511716  0.4797445 8.5436881 0.0302728 *
#April.2022-June.2021       1.576242 -2.4557301 5.6082134 0.5941200
#December.2021-August.2021  0.900573 -2.7057322 4.5068782 0.8404583
#April.2022-August.2021    -2.034902 -5.6412069 1.5714035 0.3208184
#April.2022-December.2021  -2.935475 -6.5417798 0.6708306 0.1118431

# Visualize dispersions
png('figures/MGM_Figs/SSW_MGM_pcoa_CLR_betadispersion_sampledate.png',width = 700, height = 600, res=100)
plot(mgm.disper5,main = "Centroids and Dispersion based on Aitchison Distance (CLR Data)", col=colorset1$SampDate_Color)
dev.off()

png('figures/MGM_Figs/SSW_MGM_boxplot_CLR_centroid_distance_sampledate.png',width = 700, height = 600, res=100)
boxplot(mgm.disper5,xlab="Sample Collection Date", main = "Distance to Centroid by Category (CLR Data)", sub="Based on Aitchison Distance", col=colorset1$SampDate_Color)
dev.off()

# What about between sampling depths?
mgm.disper6<-betadisper(mgm.euc.clr_dist, mgm_meta$Depth_m)
mgm.disper6

permutest(mgm.disper6, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons

anova(mgm.disper6) # p = 0.6044 --> accept the Null H, spatial medians are NOT significantly difference across sample dates

TukeyHSD(mgm.disper6) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other
#         diff       lwr      upr     p adj
#5-0   2.3886761 -4.507223 9.284575 0.6031569
#10-0  1.9710200 -4.924879 8.866919 0.7037811
#10-5 -0.4176561 -6.802018 5.966706 0.9809665

colfunc <- colorRampPalette(c("red", "blue"))
colfunc(3)

# Visualize dispersions
png('figures/MGM_Figs/ssw_mgm_pcoa_CLR_betadispersion_depth.png',width = 700, height = 600, res=100)
plot(mgm.disper6,main = "Centroids and Dispersion based on Aitchison Distance (CLR Data)", col=colfunc(3))
dev.off()

png('figures/MGM_Figs/ssw_mgm_boxplot_CLR_centroid_distance_depth.png',width = 700, height = 600, res=100)
boxplot(mgm.disper6,xlab="Sample Collection Depth", main = "Distance to Centroid by Category (CLR Data)", sub="Based on Aitchison Distance", col=colfunc(3))
dev.off()

#### Traits of Interest - Heat Maps ####
## heatmaps of traits of interest
test<-mgm.fxn.nolows[,which(colnames(mgm.fxn.nolows) %in% sulfur.fxns$KO_ID)]
test2<-test[,which(colSums(test)>=200)]

heatmap(as.matrix(test), scale = "none")

test3<-mgm.fxn.nolows[,which(colnames(mgm.fxn.nolows) %in% arsenic.fxns$KO_ID)]

heatmap(as.matrix(test3), scale = "none")

#### Relative Abundance of MGM Bins Across Samples, Depths, Seasons ####

head(mapped_meta)
# pull out total mapped reads per sample, results from BWA-mem
total_mapped_samp<-unique(data.frame(SampleID=mapped_meta$SampleID,Total_Mapped_Reads=mapped_meta$Total_Mapped_Reads))

# calculate total mapped reads per bin, and recalculate them to find total mapped reads per taxa
## this is because some bins were mapped to the same genus & species, so I am curious about these specific taxa's relative abundance in the mapped reads
gen.spec_mapped<-dcast(mapped_meta, SampleID~Genus+Species, fun.aggregate = sum, value.var="Bin_Mapped_Reads")
rownames(gen.spec_mapped)<-gen.spec_mapped$SampleID
head(gen.spec_mapped)
gen.spec_mapped2<-merge(gen.spec_mapped,total_mapped_samp, by="SampleID")
head(gen.spec_mapped2)
rownames(gen.spec_mapped2)<-gen.spec_mapped2$SampleID

# divide total mapped reads per genus&species by total mapped reads per sample to get relative abundance of reads for specific taxa across mgms
gen.sp_mapped_RA<-gen.spec_mapped2[,-c(1,ncol(gen.spec_mapped2))]/gen.spec_mapped2[,ncol(gen.spec_mapped2)]
head(gen.sp_mapped_RA)
png('figures/MGM_Figs/SSW_MGM_Genus_RelAb_Reads_4.8.23.png',width = 2200, height = 2200, res=100)
heatmap(as.matrix(t(gen.sp_mapped_RA)), scale = "none")
dev.off()

# melt down relativized ASV counts to merge with metadata
gen.sp_mapped_RA$SampleID<-rownames(gen.sp_mapped_RA)
gs_map_RA.m<-melt(gen.sp_mapped_RA)

head(gs_map_RA.m)
colnames(gs_map_RA.m)[which(names(gs_map_RA.m) == "variable")] <- "Genus_species"
colnames(gs_map_RA.m)[which(names(gs_map_RA.m) == "value")] <- "RelAb"
head(gs_map_RA.m) ## relative abundance based on sum of counts by class!

gs_map_RA.meta<-merge(gs_map_RA.m,meta_scaled, by="SampleID")
head(gs_map_RA.meta) ## relative abundance based on sum of counts by class!

# find the midpoint of RelAb
max(gs_map_RA.meta$RelAb)
min(gs_map_RA.meta$RelAb)

gs_map_RA.meta<-gs_map_RA.meta[with(gs_map_RA.meta, order(gs_map_RA.meta$SampDate, gs_map_RA.meta$Depth_m)),]
gs_map_RA.meta$SampleID <- factor(gs_map_RA.meta$SampleID, levels=unique(gs_map_RA.meta$SampleID[order(gs_map_RA.meta$SampDate, gs_map_RA.meta$Depth_m)]))

bin.hm1<-ggplot(gs_map_RA.meta, aes(gs_map_RA.meta$SampleID[order(gs_map_RA.meta$SampDate, gs_map_RA.meta$Depth_m)], Genus_species, fill= RelAb)) +geom_tile()+scale_fill_gradient2(low="skyblue",mid="white",high="orange",midpoint=0.035)+
  theme_classic()+theme(axis.title.x = element_text(size=13,vjust=-0.5),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Sample ID", y="MAG Bin Assignments", title="Relative Abundance of MAG Bins by Sample",fill="Read Relative Abundance")+scale_x_discrete(expand = c(0,0))

ggsave(bin.hm1,filename = "figures/MGM_Figs/Heatmap_MGM_Bins_RelAb_bySample_4.10.23.png", width=12, height=10, dpi=600)


tsum1<-ggplot(mapped_meta, aes(Genus, RelAb_Map)) +
  geom_jitter(aes(color=as.numeric(as.character(Depth_m))), size=2, width=0.15, height=0) +
  scale_colour_gradient(low="red",high="blue",guide = guide_colourbar(reverse = TRUE)) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Genera of MAG Bins", y="Read Relative Abundance", title="MAG Bins (Genera) by Depth",color="Depth (m)")

ggsave(tsum1,filename = "figures/MGM_Figs/TaxaSummary_MGM_Bins_RelAb_Depth_4.10.23.png", width=12, height=10, dpi=600)

tsum2<-ggplot(mapped_meta, aes(Genus, RelAb_Map)) +
  geom_jitter(aes(color=as.numeric(as.character(Depth_m)),shape=SampDate), size=2, width=0.15, height=0) +
  scale_colour_gradient(low="red",high="blue",guide = guide_colourbar(reverse = TRUE)) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Genera of MAG Bins", y="Read Relative Abundance", title="MAG Bins (Genera) by Depth & Sample Date",color="Depth (m)", shape="Sampling Date")

ggsave(tsum2,filename = "figures/MGM_Figs/TaxaSummary_MGM_Bins_RelAb_Depth_Date_4.10.23.png", width=12, height=10, dpi=600)

#### Taxonomic Beta Diversity - Read Relative Abundance data ####
mgm.clr[1:4,1:4] # sample IDs are rows, genes are columns
mgm_fxn.counts_table[1:4,1:4] # sanity check

# check rownames of CLR & VST transformed feature count data & metadata
rownames(mgm.clr) %in% rownames(meta_scaled)

## PCOA with CLR transformed data first
# calculate our Euclidean distance matrix using CLR data
mgm.euc_dist.clr <- dist(mgm.clr, method = "euclidean")

# creating our hierarcical clustering dendrogram
mgm.euc_clust <- hclust(mgm.euc_dist.clr, method="ward.D2")

# let's make it a little nicer...
mgm.euc_dend <- as.dendrogram(mgm.euc_clust, hang=0.2)
mgm.dend_cols <- as.character(mgm_meta$SampDate_Color[order.dendrogram(mgm.euc_dend)])
labels_colors(mgm.euc_dend) <- mgm.dend_cols

plot(mgm.euc_dend, ylab="CLR Euclidean Distance",cex = 0.5) + title(main = "Bacteria/Archaea Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 3)
legend("topright",legend = c("June 2021","August 2021","December 2021","April 2022"),cex=.8,col = c( "#26547c","#36ab57","#32cbff","#ff6f00"),pch = 15, bty = "n")
# Control is dark blue ("#218380"), #Alternaria is light blue ("#73d2de")
dev.off()

# let's use our Euclidean distance matrix from before
mgm.pcoa.clr <- pcoa(mgm.euc_dist.clr) # pcoa of euclidean distance matrix = PCA of euclidean distance matrix
##save.image("data/ssw_clr.euc.dist1_3.7.23.Rdata")

# The proportion of variances explained is in its element values$Relative_eig
mgm.pcoa.clr$values

# extract principal coordinates
mgm.pcoa.clr.vectors<-data.frame(mgm.pcoa.clr$vectors)
mgm.pcoa.clr.vectors$SampleID<-rownames(mgm.pcoa.clr$vectors)

# merge pcoa coordinates w/ metadata
mgm.pcoa.clr.meta<-merge(mgm.pcoa.clr.vectors, mgm_meta, by.x="SampleID", by.y="SampleID")
mgm.pcoa.clr.meta$SampleMonth
mgm.pcoa.clr.meta$SampDate

head(mgm.pcoa.clr.meta)

head(mgm.pcoa.clr.meta$values) # pull out Relative (Relative_eig) variation % to add to axes labels

# create PCoA ggplot fig
pcoa1<-ggplot(mgm.pcoa.clr.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(SampDate)), size=4)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Seawater",subtitle="Using Centered-Log Ratio Feature Data",xlab="PC1 [41.14%]", ylab="PC2 [9.04%]",color="Sample Type")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Type",values=unique(mgm.pcoa.clr.meta$SampDate_Color[order(mgm.pcoa.clr.meta$SampDate)]),labels=c("June.2021"="June 2021","August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("PC1 [23.56%]") + ylab("PC2 [20.45%]")

ggsave(pcoa1,filename = "figures/MGM_Figs/SSW_MGM_pcoa_CLR_sampdate.png", width=12, height=10, dpi=600)

# sample month shape, depth color
pcoa2<-ggplot(mgm.pcoa.clr.meta, aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampleMonth), size=5)+theme_bw()+
  labs(title="PCoA: Metagenome Functions in Salton Seawater",subtitle="Using Centered-Log Ratio Feature Data",xlab="PC1", ylab="PC2",color="Depth (m)")+
  theme_classic()+ theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),axis.text = element_text(size=12),axis.text.x = element_text(vjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17))+
  scale_color_continuous(low="blue3",high="red",trans = 'reverse') + scale_shape_discrete(labels=c("June 2021","August 2021","December 2021","April 2022"),name="Sample Date") +
  xlab("PC1 [23.56%]") + ylab("PC2 [20.45%]")

ggsave(pcoa2,filename = "figures/MGM_Figs/SSW_MGM_pcoa_CLR.traits_depth.png", width=12, height=10, dpi=600)

## betadisper to look at within group variance

# first by sampling date
mgm.disper1<-betadisper(mgm.euc_dist.clr, mgm_meta$SampDate)
mgm.disper1

permutest(mgm.disper1, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
#Pairwise comparisons:
#  (Observed p-value below diagonal, permuted p-value above diagonal)
#               June.2021 August.2021 December.2021 April.2022
#June.2021                 0.13000000    0.03700000      0.010
#August.2021   0.12533884                0.50300000      0.189
#December.2021 0.00835232  0.56548064                    0.038
#April.2022    0.00043402  0.20006152    0.00656846

anova(mgm.disper1) # p = 0.02656 --> reject the Null H, spatial medians ARE significantly difference across sample dates

TukeyHSD(mgm.disper1) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other

#                             diff        lwr       upr     p adj
#August.2021-June.2021      3.611143 -0.4208285 7.6431151 0.0784173 .
#December.2021-June.2021    4.511716  0.4797445 8.5436881 0.0302728 *
#April.2022-June.2021       1.576242 -2.4557301 5.6082134 0.5941200
#December.2021-August.2021  0.900573 -2.7057322 4.5068782 0.8404583
#April.2022-August.2021    -2.034902 -5.6412069 1.5714035 0.3208184
#April.2022-December.2021  -2.935475 -6.5417798 0.6708306 0.1118431

# Visualize dispersions
png('figures/MGM_Figs/SSW_MGM_pcoa_clr_betadispersion_sampledate.png',width = 700, height = 600, res=100)
plot(mgm.disper1,main = "Centroids and Dispersion based on Aitchison Distance (CLR Data)", col=colorset1$SampDate_Color)
dev.off()

png('figures/MGM_Figs/SSW_MGM_boxplot_clr_centroid_distance_sampledate.png',width = 700, height = 600, res=100)
boxplot(mgm.disper1,xlab="Sample Collection Date", main = "Distance to Centroid by Category (CLR Data)", sub="Based on Aitchison Distance", col=colorset1$SampDate_Color)
dev.off()

# What about between sampling depths?
mgm.disper2<-betadisper(mgm.euc_dist.clr, mgm_meta$Depth_m)
mgm.disper2

permutest(mgm.disper2, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons

anova(mgm.disper2) # p = 0.6044 --> accept the Null H, spatial medians are NOT significantly difference across sample dates

TukeyHSD(mgm.disper2) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other
#         diff       lwr      upr     p adj
#5-0   2.3886761 -4.507223 9.284575 0.6031569
#10-0  1.9710200 -4.924879 8.866919 0.7037811
#10-5 -0.4176561 -6.802018 5.966706 0.9809665

colfunc <- colorRampPalette(c("red", "blue"))
colfunc(3)

# Visualize dispersions
png('figures/MGM_Figs/ssw_mgm_pcoa_clr_betadispersion_depth.png',width = 700, height = 600, res=100)
plot(mgm.disper2,main = "Centroids and Dispersion based on Aitchison Distance (CLR Data)", col=colfunc(3))
dev.off()

png('figures/MGM_Figs/ssw_mgm_boxplot_clr_centroid_distance_depth.png',width = 700, height = 600, res=100)
boxplot(mgm.disper2,xlab="Sample Collection Depth", main = "Distance to Centroid by Category (CLR Data)", sub="Based on Aitchison Distance", col=colfunc(3))
dev.off()


### Export Global Env for Other Scripts ####
#save.image("data/Metagenomes/Analysis/mgm_analysis.Rdata")
# ^ includes all data combined in object bac.dat.all, ASV table (samples are rows, ASVs are columns), mgm_meta, and an ASV count table (where ASVs are rows, not columns)
# Version Information
sessionInfo()