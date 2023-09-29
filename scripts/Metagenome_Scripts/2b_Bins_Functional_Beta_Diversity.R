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
  #library(heatmaply)
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

#### Load Data ####
#load("data/Metagenomes/Analysis/mgm_analysis.Rdata") # load Rdata to global env
load("data/Metagenomes/Analysis/mgm_MAG_analysis.Rdata")
#load("data/Metagenomes/Analysis/SSW_MAG_Bin_Fxn_BetaDiv.Rdata")
#load("data/Metagenomes/Analysis/SSW_MAG_Bin_Fxn_BetaDiv_SulfurMetab.Rdata")

head(bin_meta_scaled)
arsen.kegg[1:4,]
bin.ko.cov.sum_table[1:4,1:4]
head(bin.clr.ars)

# fixing some col names in bin_meta_scaled
#colnames(bin_meta_scaled)[which(names(bin_meta_scaled) == "DO_Percent_Local")] <- "DO_Percent_Local"
#colnames(bin_meta_scaled)[which(names(bin_meta_scaled) == "Dissolved_Organic Matter_RFU")] <- "Dissolved_OrganicMatter_RFU"

#bin_meta_scaled$SampDate<-gsub("\\."," ",bin_meta_scaled$SampDate) # drop period between month & year in SampDate col
#bin_meta_scaled$SampDate<-factor(bin_meta_scaled$SampDate,levels=c("August 2021","December 2021","April 2022"))

# Before transformations (i.e., VST, CLR, etc) were done, the following was performed
# featureCounts counted reads that mapped to genes
# Reads mapped to genes / gene length for all genes across all samples
# Gene coverage was then added together for each KO ID, since multiple genes were assigned the same KO ID
# Summed coverage per KO was then transformed via median-ratio, vst, and clr

#### Functional Beta Diversity - CLR data ####
bin.clr[1:4,1:4] # sample IDs are rows, genes are columns
bin.ko.cov.sum_table[1:4,1:4] # sanity check

# check rownames of CLR & VST transformed feature count data & metadata
rownames(bin.clr) %in% rownames(bin_meta_scaled)

## PCOA with CLR transformed data first
# calculate our Euclidean distance matrix using CLR data
bin.euc.clr_dist <- dist(bin.clr, method = "euclidean")

# creating our hierarcical clustering dendrogram
bin.euc.clr_clust <- hclust(bin.euc.clr_dist, method="ward.D2")

# let's make it a little nicer...
bin.euc.clr_dend <- as.dendrogram(bin.euc.clr_clust, hang=0.2)
bin.dend_cols <- as.character(bin_meta_scaled$SampDate_Color[order.dendrogram(bin.euc.clr_dend)])
labels_colors(bin.euc.clr_dend) <- bin.dend_cols

plot(bin.euc.clr_dend, ylab="CLR Euclidean Distance",cex = 0.5) + title(main = "Bacteria/Archaea Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 3)
legend("topright",legend = c("August 2021","December 2021","April 2022"),cex=.8,col = c("#ef781c","#03045e","#059c3f"),pch = 15, bty = "n")
# Control is dark blue ("#218380"), #Alternaria is light blue ("#73d2de")
dev.off()

# let's use our Euclidean distance matrix from before
bin.pcoa.clr <- pcoa(bin.euc.clr_dist) # pcoa of euclidean distance matrix = PCA of euclidean distance matrix

# The proportion of variances explained is in its element values$Relative_eig
bin.pcoa.clr$values

# extract principal coordinates
bin.pcoa.clr.vectors<-data.frame(bin.pcoa.clr$vectors)
bin.pcoa.clr.vectors$Bin_ID<-rownames(bin.pcoa.clr$vectors)

# merge pcoa coordinates w/ metadata
bin.pcoa.clr.meta<-merge(bin.pcoa.clr.vectors, bin_meta, by.x="Bin_ID", by.y="Bin_ID")
bin.pcoa.clr.meta$SampleMonth
bin.pcoa.clr.meta$SampDate

head(bin.pcoa.clr.meta)

head(bin.pcoa.clr$values) # pull out Relative (Relative_eig) variation % to add to axes labels

# create PCoA ggplot fig
bin.pcoa1<-ggplot(bin.pcoa.clr.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(SampDate)), size=4)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Seawater",subtitle="Using CLR Transformed, Summed Gene Coverage per KO Function",xlab="PC1 [41.14%]", ylab="PC2 [9.04%]",color="Sample Type")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Type",values=unique(bin.pcoa.clr.meta$SampDate_Color[order(bin.pcoa.clr.meta$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("PC1 [32.17%]") + ylab("PC2 [11.98%]")

#ggsave(bin.pcoa1,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/SSW_MAG_Bin_pcoa_CLR_SummedCoverage_Per_KO_sampdate.png", width=12, height=10, dpi=600)

# sample month shape, depth color
bin.pcoa2<-ggplot(bin.pcoa.clr.meta, aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(color=as.numeric(as.character(Depth_m)),shape=SampleMonth), size=5)+theme_bw()+
  labs(title="PCoA: Metagenome Functions in Salton Seawater",subtitle="Using CLR Transformed, Summed Gene Coverage per KO Function",xlab="PC1", ylab="PC2",color="Depth (m)")+
  theme_classic()+ theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),axis.text = element_text(size=12),axis.text.x = element_text(vjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17))+
  scale_color_continuous(low="blue3",high="red",trans = 'reverse') + scale_shape_discrete(labels=c("August 2021","December 2021","April 2022"),name="Sample Date") +
  xlab("PC1 [32.17%]") + ylab("PC2 [11.98%]")

ggplot(bin.pcoa.clr.meta, aes(x=Axis.1, y=Axis.2)) +theme_bw()+
  geom_text(aes(color=as.numeric(as.character(Depth_m)),label=Bin_ID), size=3) +
  labs(title="PCoA: Metagenome Functions in Salton Seawater",subtitle="Using CLR Transformed, Summed Gene Coverage per KO Function",xlab="PC1", ylab="PC2",color="Depth (m)")+
  theme_classic()+ theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),axis.text = element_text(size=12),axis.text.x = element_text(vjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17))+
  scale_color_continuous(low="blue3",high="red",trans = 'reverse') + scale_shape_discrete(labels=c("August 2021","December 2021","April 2022"),name="Sample Date") +
  xlab("PC1 [32.17%]") + ylab("PC2 [11.98%]")
#ggsave(bin.pcoa2,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/SSW_MAG_Bin_pcoa_CLR_SummedCoverage_Per_KO.traits_depth.png", width=12, height=10, dpi=600)


#### Pull Out Sulfur Metabolic Fxns from CLR data ####
## heatmaps of traits of interest

bin.clr[1:4,1:4]

# pull out sulfur functions from CLR transformed, summed coverages (summed coverage per KO)
sulf.ko.bin<-bin.clr[,which(colnames(bin.clr) %in% sulfur.fxns.bins$KO_ID)] # merge CLR data w/ S fxns found in contigs from KOFamScan
sulf.ko.bin$Bin_ID<-rownames(sulf.ko.bin)
sulf.ko.bin.melt<-melt(sulf.ko.bin, by="Bin_ID")
colnames(sulf.ko.bin.melt)[which(names(sulf.ko.bin.melt) == "variable")] <- "KO_ID"
colnames(sulf.ko.bin.melt)[which(names(sulf.ko.bin.melt) == "value")] <- "CLR_SumCovPerKO"
head(sulf.ko.bin.melt) #sanity check

clr.sulf.ko.bin<-merge(sulf.ko.bin.melt,sulf.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
head(clr.sulf.ko.bin)
colnames(clr.sulf.ko.bin)[which(names(clr.sulf.ko.bin) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
clr.cov.sum.sulf.ko.bin<-as.data.frame(dcast(clr.sulf.ko.bin, Bin_ID~KO_Function.KEGG, value.var="CLR_SumCovPerKO", fun.aggregate=sum)) ###just dcast, nothing is being added here!
rownames(clr.cov.sum.sulf.ko.bin)<-clr.cov.sum.sulf.ko.bin$Bin_ID
clr.cov.sum.sulf.ko.bin[1:4,]

# sanity check
clr.cov.sum.sulf.ko.bin$`cysH; phosphoadenosine phosphosulfate reductase [EC:1.8.4.8 1.8.4.10]`[1:4]
head(clr.sulf.ko.bin)

#### Pull Out Sulfur Metabolic Fxns from CLR data - with NAs ####
## heatmaps of traits of interest

bin.clr.na[1:4,1:4]

# pull out sulfur functions from CLR transformed, summed coverages (summed gene coverage per KO)
sulf.ko.na.bin<-bin.clr.na[,which(colnames(bin.clr.na) %in% sulfur.fxns.bins$KO_ID)] # merge CLR data w/ S fxns found in contigs from KOFamScan
dim(sulf.ko.na.bin) # 35 total bins, 14 columns aka functions
sulf.ko.na.bin<-sulf.ko.na.bin[rowSums(is.na(sulf.ko.na.bin)) != ncol(sulf.ko.na.bin), ] # drop all rows that contain ONLY NAs
# is.na() tells us which elements are NA (TRUE) or not NA (FALSE)
# rowSums(is.na(df)) tells us how many columns have NAs; in this df, there are a total of 14 columns
# if rowSums(is.na(df)) != total # of rows, aka if the total columns per row with NAs is LESS than the total # of columns per row, then keep the row
# ^ if the row contains only NAs, and rowSums(is.na(df)) is NOT 14 (total # of columns), then we keep it

dim(sulf.ko.na.bin) # dropped 20 bins that do not have any S fxns
sulf.ko.na.bin$Bin_ID<-rownames(sulf.ko.na.bin)

sulf.ko.na.bin.melt<-melt(sulf.ko.na.bin, by="Bin_ID")
colnames(sulf.ko.na.bin.melt)[which(names(sulf.ko.na.bin.melt) == "variable")] <- "KO_ID"
colnames(sulf.ko.na.bin.melt)[which(names(sulf.ko.na.bin.melt) == "value")] <- "CLR_SumCovPerKO"
head(sulf.ko.na.bin.melt) #sanity check

clr.sulf.ko.na.bin<-merge(sulf.ko.na.bin.melt,sulf.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
head(clr.sulf.ko.na.bin)
colnames(clr.sulf.ko.na.bin)[which(names(clr.sulf.ko.na.bin) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
clr.cov.sum.sulf.ko.na.bin<-as.data.frame(dcast(clr.sulf.ko.na.bin, Bin_ID~KO_Function.KEGG, value.var="CLR_SumCovPerKO", fun.aggregate=sum)) ###just dcast, nothing is being added here!
rownames(clr.cov.sum.sulf.ko.na.bin)<-clr.cov.sum.sulf.ko.na.bin$Bin_ID
clr.cov.sum.sulf.ko.na.bin[1:4,]

# sanity check
clr.cov.sum.sulf.ko.na.bin$`cysH; phosphoadenosine phosphosulfate reductase [EC:1.8.4.8 1.8.4.10]`[1:6]
head(clr.sulf.ko.na.bin)

#### Pull Out Specific Genes by S Pathways ####
bin.ko.cov.sum_table[1:4,1:4] # contains the sum of coverages per gene per KO -- featureCounts was normalized by gene length across samples first to get coverage, then summed up per KO ID

## pull out all KOs in each Pathway
assim.sulfate.red.bin<-data.frame(KO_Function.KEGG=unique(clr.sulf.ko.bin$KO_Function[which(clr.sulf.ko.bin$Pathway=="Assimilatory Sulfate Reduction")]))
dissim.sulfate.redox.bin<-data.frame(KO_Function.KEGG=unique(clr.sulf.ko.bin$KO_Function.KEGG[which(clr.sulf.ko.bin$Pathway=="Dissimilatory Sulfate Redox")]))
mult.sulf.bin<-data.frame(KO_Function.KEGG=unique(clr.sulf.ko.bin$KO_Function.KEGG[which(clr.sulf.ko.bin$Pathway=="Multiple Pathways")]))
sox.system.bin<-data.frame(KO_Function.KEGG=unique(clr.sulf.ko.bin$KO_Function.KEGG[which(clr.sulf.ko.bin$Pathway=="SOX System")]))

# pull out functions & CLR info per pathway
asSO4.ko.cov.bin<-clr.cov.sum.sulf.ko.bin[,-1][,colnames(clr.cov.sum.sulf.ko.bin[,-1]) %in% assim.sulfate.red.bin$KO_Function.KEGG] # pull out sox genes from gene list found in CLR transformed cov per KO
disSO4.ko.cov.bin<-clr.cov.sum.sulf.ko.bin[,-1][,colnames(clr.cov.sum.sulf.ko.bin[,-1]) %in% dissim.sulfate.redox.bin$KO_Function.KEGG] # pull out sox genes from gene list found in CLR transformed cov per KO
multiS.ko.cov.bin<-clr.cov.sum.sulf.ko.bin[,-1][,colnames(clr.cov.sum.sulf.ko.bin[,-1]) %in% mult.sulf.bin$KO_Function.KEGG] # pull out sox genes from gene list found in CLR transformed cov per KO
sox.ko.cov.bin<-clr.cov.sum.sulf.ko.bin[,-1][,colnames(clr.cov.sum.sulf.ko.bin[,-1]) %in% sox.system.bin$KO_Function.KEGG] # pull out sox genes from gene list found in CLR transformed cov per KO

#### Sulfur Heat Maps ####
# see max & mean of summed
max(clr.cov.sum.sulf.ko.bin[,-1])
min(as.matrix(clr.cov.sum.sulf.ko.bin[,-1]))

# first heat map of sulfur KOs
heatmap(as.matrix(clr.cov.sum.sulf.ko.bin[,-1]), scale = "none")

colSums(clr.cov.sum.sulf.ko.bin[,-1])
#clr.cov.sum.sulf.ko2 <- clr.cov.sum.sulf.ko.bin[,which(colSums(clr.cov.sum.sulf.ko.bin[,-1])>10)]

heatmap(as.matrix(clr.cov.sum.sulf.ko.bin[,-1]), scale = "none")

# prep for ggplot2 heatmap
head(clr.sulf.ko.na.bin)
clr.sulf.all.bin1<-merge(clr.sulf.ko.na.bin,bin_meta_scaled,by="Bin_ID")
clr.sulf.all.bin<-merge(clr.sulf.all.bin1,mag_tax,by=c("Bin_ID","PlotBin"))

head(clr.sulf.all.bin)
clr.sulf.all.bin$PlotBin = factor(clr.sulf.all.bin$PlotBin, levels=unique(clr.sulf.all.bin$PlotBin[order(clr.sulf.all.bin$SampDate,clr.sulf.all.bin$Depth_m)]), ordered=TRUE)
clr.sulf.all.bin$SampDate<-gsub("\\."," ",clr.sulf.all.bin$SampDate)
clr.sulf.all.bin$SampDate<-factor(clr.sulf.all.bin$SampDate, levels=c("August 2021","December 2021","April 2022"))

clr.sulf.all.bin$PathShort<-clr.sulf.all.bin$Pathway
clr.sulf.all.bin$PathShort[(clr.sulf.all.bin$PathShort) == "Dissimilatory Sulfate Redox"] <- "D.SO4 RedOx"
clr.sulf.all.bin$PathShort[(clr.sulf.all.bin$PathShort) == "Assimilatory Sulfate Reduction"] <- "A.SO4 Red"
clr.sulf.all.bin$PathShort[(clr.sulf.all.bin$PathShort) == "Multiple Pathways"] <- "Multi Paths"
clr.sulf.all.bin$PathShort[(clr.sulf.all.bin$PathShort) == "S Disproportionation"] <- "S Disprop."

clr.sulf.all.bin$Pathway<-factor(clr.sulf.all.bin$Pathway,levels=c("Assimilatory Sulfate Reduction","Dissimilatory Sulfate Redox","Multiple Pathways","SOX","S Disproportionation"))
clr.sulf.all.bin$PathShort<-factor(clr.sulf.all.bin$PathShort,levels=c("A.SO4 Red","D.SO4 RedOx","Multi Paths","SOX","S Disprop."))

clr.sulf.all.bin$PathSpecShort<-clr.sulf.all.bin$PathwaySpecific
clr.sulf.all.bin$PathSpecShort[(clr.sulf.all.bin$PathSpecShort) == "Dissimilatory Sulfate Redox"] <- "D.SO4 RedOx"
clr.sulf.all.bin$PathSpecShort[(clr.sulf.all.bin$PathSpecShort) == "Assimilatory Sulfate Reduction"] <- "A.SO4 Red"
clr.sulf.all.bin$PathSpecShort[(clr.sulf.all.bin$PathSpecShort) == "Multiple Pathways"] <- "MultiPaths"
clr.sulf.all.bin$PathSpecShort[(clr.sulf.all.bin$PathSpecShort) == "Sulfur Disproportionation"] <- "S Disprop."
clr.sulf.all.bin$PathSpecShort[(clr.sulf.all.bin$PathSpecShort) == "Sulfide Oxidation"] <- "H2S Ox"
clr.sulf.all.bin$PathSpecShort[(clr.sulf.all.bin$PathSpecShort) == "Sulfite Oxidation"] <- "SO3 Ox"
clr.sulf.all.bin$PathSpecShort[(clr.sulf.all.bin$PathSpecShort) == "Thiosulfate Oxidation"] <- "S2O3 Ox"

clr.sulf.all.bin$PathwaySpecific<-factor(clr.sulf.all.bin$PathwaySpecific,levels=c("Assimilatory Sulfate Reduction","Dissimilatory Sulfate Redox","Multiple Pathways","SOX","S Disproportionation","Sulfide Oxidation","Sulfite Oxidation","Thiosulfate Oxidation"))
clr.sulf.all.bin$PathSpecShort<-factor(clr.sulf.all.bin$PathSpecShort,levels=c("A.SO4 Red","D.SO4 RedOx","MultiPaths","S Disprop.","H2S Ox","SO3 Ox","S2O3 Ox"))

head(clr.sulf.all.bin)

# For heatmap color gradient
max(clr.sulf.all.bin$CLR_SumCovPerKO, na.rm=TRUE)
max(clr.sulf.all.bin$CLR_SumCovPerKO, na.rm=TRUE)/2
min(clr.sulf.all.bin$CLR_SumCovPerKO, na.rm=TRUE)

# Figures

sulf.hm1a<-ggplot(clr.sulf.all.bin, aes(PlotBin, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.9","0.55","0.2"),breaks=c(0.9,0.55,0.2)) + labs(title="Sulfur Metabolism in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(sulf.hm1a,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Bins_BinID_by_Function_heatmap.png", width=20, height=13, dpi=600)

sulf.hm1a2<-ggplot(clr.sulf.all.bin, aes(PlotBin, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.9","0.55","0.2"),breaks=c(0.9,0.55,0.2)) + labs(title="Sulfur Metabolism in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(SampDate~.,scales="free_y", space = "free")

ggsave(sulf.hm1a2,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Bins_BinID_by_Function_heatmap2.png", width=17, height=15, dpi=600)

sulf.hm1b<-ggplot(clr.sulf.all.bin, aes(PlotBin, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.9","0.55","0.2"),breaks=c(0.9,0.55,0.2)) + labs(title="Sulfur Metabolism in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Pathway~.,scales="free_y", space = "free")

ggsave(sulf.hm1b,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Bins_BinID_by_Function_Pathway_heatmap.png", width=20, height=15, dpi=600)

sulf.hm1b2<-ggplot(clr.sulf.all.bin, aes(PlotBin, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.9","0.55","0.2"),breaks=c(0.9,0.55,0.2)) + labs(title="Sulfur Metabolism in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~.,scales="free_y", space = "free")

ggsave(sulf.hm1b2,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Bins_BinID_by_Function_Pathway_heatmap2.png", width=17, height=15, dpi=600)

sulf.hm1b3<-ggplot(clr.sulf.all.bin, aes(Genus, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.9","0.55","0.2"),breaks=c(0.9,0.55,0.2)) + labs(title="Sulfur Metabolism in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~.,scales="free_y", space = "free")

ggsave(sulf.hm1b3,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Bins_Genus_by_Function_Pathway_heatmap.png", width=17, height=15, dpi=600)

sulf.hm1c<-ggplot(clr.sulf.all.bin, aes(PlotBin, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.9","0.55","0.2"),breaks=c(0.9,0.55,0.2)) + labs(title="Sulfur Metabolism in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathwaySpecific~SampDate, scales="free", space = "free")

ggsave(sulf.hm1c,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Bins_BinID_by_Function_SampDate_PathwaySpecific_best_heatmap.png", width=20, height=20, dpi=600)

sulf.hm1c2<-ggplot(clr.sulf.all.bin, aes(PlotBin, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.9","0.55","0.2"),breaks=c(0.9,0.55,0.2)) + labs(title="Sulfur Metabolism in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathSpecShort~SampDate, scales="free", space = "free")

ggsave(sulf.hm1c2,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Bins_BinID_by_Function_SampDate_PathwaySpecific_best_heatmap2.png", width=20, height=20, dpi=600)

sulf.hm1c3<-ggplot(clr.sulf.all.bin, aes(Genus, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.9","0.55","0.2"),breaks=c(0.9,0.55,0.2)) + labs(title="Sulfur Metabolism in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathwaySpecific~SampDate, scales="free", space = "free")

ggsave(sulf.hm1c3,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Bins_Genus_by_Function_SampDate_PathwaySpecific_best_heatmap.png", width=20, height=20, dpi=600)

sulf.hm1c4<-ggplot(clr.sulf.all.bin, aes(Genus, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.9","0.55","0.2"),breaks=c(0.9,0.55,0.2)) + labs(title="Sulfur Metabolism in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathSpecShort~SampDate, scales="free", space = "free")

ggsave(sulf.hm1c4,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Bins_Genus_by_Function_SampDate_PathwaySpecific_best_heatmap2.png", width=20, height=20, dpi=600)

sulf.hm1d2<-ggplot(clr.sulf.all.bin, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.9","0.55","0.2"),breaks=c(0.9,0.55,0.2)) + labs(title="Sulfur Metabolism in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathwaySpecific~SampDate, scales="free", space = "free")

ggsave(sulf.hm1d2,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Bins_Depth_by_Function_SampDate_PathwaySpecific_best_heatmap.png", width=20, height=20, dpi=600)

sulf.hm1d2a<-ggplot(clr.sulf.all.bin, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.9","0.55","0.2"),breaks=c(0.9,0.55,0.2)) + labs(title="Sulfur Metabolism in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathSpecShort~SampDate, scales="free", space = "free")

ggsave(sulf.hm1d2a,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Bins_Depth_by_Function_SampDate_PathwaySpecific_best_heatmap2.png", width=20, height=20, dpi=600)

# sulf.hm1c<-ggplot(clr.sulf.all.bin, aes(interaction(SampDate,Depth_m), KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.9","0.55","0.2"),breaks=c(0.9,0.55,0.2)) + labs(title="Sulfur Metabolism in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Pathway~.,scales="free_y", space = "free")
#
# ggsave(sulf.hm1c,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Bins_SampDate_Depth_by_Function_Pathway_heatmap.png", width=15, height=18, dpi=600)

sulf.hm1d<-ggplot(clr.sulf.all.bin, aes(PlotBin, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.9","0.55","0.2"),breaks=c(0.9,0.55,0.2)) + labs(title="Sulfur Metabolism in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,vjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate,scales="free_x", space = "free")

ggsave(sulf.hm1d,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Bins_Depth_by_Function_SampDate_best_heatmap.png", width=20, height=15, dpi=600)

# sulf.hm1e<-ggplot(clr.sulf.all.bin, aes(Bin_ID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.9","0.55","0.2"),breaks=c(0.9,0.55,0.2)) + labs(title="Sulfur Metabolism in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,vjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Pathway~SampDate,scales="free_x", space = "free")
#
# ggsave(sulf.hm1e,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Bins_Depth_by_Function_SampDate_best_heatmap.png", width=20, height=13, dpi=600)

# sulf.hm1f<-ggplot(clr.sulf.all.bin, aes(Pathway, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.9","0.55","0.2"),breaks=c(0.9,0.55,0.2)) + labs(title="Sulfur Metabolism in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11,face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Depth_m~SampDate,scales="free", space = "free")
#
# ggsave(sulf.hm1f,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Bins_heatmap1d.png", width=18, height=18, dpi=600)
#
# sulf.hm1g<-ggplot(clr.sulf.all.bin, aes(Pathway, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.9","0.55","0.2"),breaks=c(0.9,0.55,0.2)) + labs(title="Sulfur Metabolism in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11,face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_wrap(.~SampDate)
#
# ggsave(sulf.hm1g,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Bins_heatmap1d.png", width=18, height=18, dpi=600)

# sulf.hm1e<-ggplot(clr.sulf.all.bin[clr.sulf.all.bin$Depth_m==0,], aes(Pathway, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.9","0.55","0.2"),breaks=c(0.9,0.55,0.2)) + labs(title="Sulfur Metabolism in Salton Seawater MAGs - 0m",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11,face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)
#
# ggsave(sulf.hm1e,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/Sulfur_KOFxns_Pathways_MGMs_Bins_0m_heatmap.png", width=18, height=18, dpi=600)
#
# sulf.hm1f<-ggplot(clr.sulf.all.bin[clr.sulf.all.bin$Depth_m==5,], aes(Pathway, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.9","0.55","0.2"),breaks=c(0.9,0.55,0.2)) + labs(title="Sulfur Metabolism in Salton Seawater MAGs - 5m",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)
#
# ggsave(sulf.hm1f,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/Sulfur_KOFxns_Pathways_MGMs_Bins_5m_heatmap.png", width=18, height=18, dpi=600)
#
# sulf.hm1g<-ggplot(clr.sulf.all.bin[clr.sulf.all.bin$Depth_m==10,], aes(Pathway, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.9","0.55","0.2"),breaks=c(0.9,0.55,0.2)) + labs(title="Sulfur Metabolism in Salton Seawater MAGs - 10m",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)
#
# ggsave(sulf.hm1g,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/Sulfur_KOFxns_Pathways_MGMs_Bins_10m_heatmap.png", width=18, height=18, dpi=600)
# #
# # pull out specific S functions
# ## first, SOX
# clr.Sox<-clr.sulf.all.bin[grepl('Sox', clr.sulf.all.bin$KO_Function),] # pull out just Sox functions
# clr.Sox$SampDate = factor(clr.Sox$SampDate, levels=c("August 2021","December 2021", "April 2022"))
# clr.Sox$SampleID = factor(clr.Sox$SampleID, levels=unique(clr.Sox$SampleID[order(clr.Sox$SampDate,clr.Sox$Depth_m)]), ordered=TRUE)
#
# NA %in% clr.Sox$CLR_SumCovPerKO
#
# s.sox.hm<-ggplot(clr.Sox, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
#   geom_tile(colour="white",size=0.25)  +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8",na.value="grey50",labels=c("0.9","0.55","0.2"),breaks=c(0.9,0.55,0.2)) + labs(title="SOX Functions in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO, Grouped by Bin Assigment",fill="CLR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=15),
#         axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.border=element_blank(),panel.background = element_rect(fill = "white", colour = NA)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0)) + facet_grid(.~SampDate)
#
# ggsave(s.sox.hm,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/SSW_S_SOX_Contigs_bySampDate_Depth_heatmap.png", width=15, height=10, dpi=600)
#
# # Assimilatory sulfate reduction
# clr.as.S.redox<-clr.sulf.all.bin[grepl('1.8.7.1|1.8.1.2|2.7.7.4|1.8.4.10|1.8.4.8|2.7.1.25', clr.sulf.all.bin$KO_Function),] # pull out just assimilatory sulfate reduction functions
# clr.as.S.redox$SampDate = factor(clr.as.S.redox$SampDate, levels=c("August 2021","December 2021", "April 2022"))
# clr.as.S.redox$SampleID = factor(clr.as.S.redox$SampleID, levels=unique(clr.as.S.redox$SampleID[order(clr.as.S.redox$SampDate,clr.as.S.redox$Depth_m)]), ordered=TRUE)
#
# NA %in% clr.as.S.redox$CLR_SumCovPerKO
#
# s.R.hm1<-ggplot(clr.as.S.redox, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
#   geom_tile(colour="white",size=0.25)  +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8",na.value="grey50",labels=c("0.9","0.55","0.2"),breaks=c(0.9,0.55,0.2)) + labs(title="Assimilatory Sulfuate Reduction in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO, Grouped by Bin Assigment",fill="CLR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=15),
#         axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.border=element_blank(),panel.background = element_rect(fill = "white", colour = NA)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0)) + facet_grid(.~SampDate)
#
# ggsave(s.R.hm1,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/SSW_S_AssSO4_Reduction_Contigs_bySampDate_Depth_heatmap.png", width=15, height=10, dpi=600)
#
# # Dissimilatory sulfate reduction and oxidation
# clr.dis.S.redox<-clr.sulf.all.bin[grepl('1.8.99.2|1.8.99.5|2.7.7.4', clr.sulf.all.bin$KO_Function),] # pull out just Sox functions
# clr.dis.S.redox$SampDate = factor(clr.dis.S.redox$SampDate, levels=c("August 2021","December 2021", "April 2022"))
# clr.dis.S.redox$SampleID = factor(clr.dis.S.redox$SampleID, levels=unique(clr.dis.S.redox$SampleID[order(clr.dis.S.redox$SampDate,clr.dis.S.redox$Depth_m)]), ordered=TRUE)
#
# NA %in% clr.dis.S.redox$CLR_SumCovPerKO
#
# s.RO.hm1<-ggplot(clr.dis.S.redox, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
#   geom_tile(colour="white",size=0.25)  +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8",na.value="grey50",labels=c("0.9","0.55","0.2"),breaks=c(0.9,0.55,0.2)) + labs(title="Dissimilarity Sulfuate RedOx in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO, Grouped by Bin Assigment",fill="CLR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=15),
#         axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.border=element_blank(),panel.background = element_rect(fill = "white", colour = NA)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0)) + facet_grid(.~SampDate)
#
# ggsave(s.RO.hm1,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/SSW_S_DissSO4_RedOx_Contigs_bySampDate_Depth_heatmap.png", width=15, height=10, dpi=600)

#### Look at Specific S Gene Coverage Across Samples ####

head(clr.cov.sum.sulf.ko.na.bin) # columns are genes in this df

# merge with scaled metadata and prep for scatterplots of traits across samples
clr.sulf.trait.table.bin1<-merge(clr.cov.sum.sulf.ko.na.bin,bin_meta_scaled,by="Bin_ID")
clr.sulf.trait.table.bin<-merge(clr.sulf.trait.table.bin1,mag_tax,by=c("Bin_ID","PlotBin"))

head(clr.sulf.trait.table.bin)
clr.sulf.trait.table.bin$PlotBin = factor(clr.sulf.trait.table.bin$PlotBin, levels=unique(clr.sulf.trait.table.bin$PlotBin[order(clr.sulf.trait.table.bin$SampDate,clr.sulf.trait.table.bin$Depth_m)]), ordered=TRUE)
clr.sulf.trait.table.bin$SampDate<-gsub("\\."," ",clr.sulf.trait.table.bin$SampDate)
clr.sulf.trait.table.bin$SampDate<-factor(clr.sulf.trait.table.bin$SampDate, levels=c("August 2021","December 2021","April 2022"))

head(clr.sulf.trait.table.bin)

# Note: not looking at every S cycling gene included in this project but looking at ones that appear to have noticeable trends in heat maps

# First by Bin ID

### SOX genes
# `soxY; sulfur-oxidizing protein SoxY`
soxy.bin.fs<-ggplot(clr.sulf.trait.table.bin, aes(x=PlotBin, y=`soxY; sulfur-oxidizing protein SoxY`,color=SampDate,group=SampDate)) + geom_jitter(aes(color=SampDate), size=3, width=0.15, height=0) + theme_bw()+
  labs(title="SoxY Depth of Coverage in MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table.bin$SampDate_Color[order(clr.sulf.trait.table.bin$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Bin ID") + ylab("CLR-Transformed Coverage")

ggsave(soxy.bin.fs,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/FxnScatterplots/SoxY_CLR_Coverage_BinID_fxn.sum.png", width=12, height=10, dpi=600)

# `soxA; L-cysteine S-thiosulfotransferase [EC:2.8.5.2]`
soxa.bin.fs<-ggplot(clr.sulf.trait.table.bin, aes(x=PlotBin, y=`soxA; L-cysteine S-thiosulfotransferase [EC:2.8.5.2]`,color=SampDate,group=SampDate)) + geom_jitter(aes(color=SampDate), size=3, width=0.15, height=0) + theme_bw()+
  labs(title="SoxA Depth of Coverage in MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table.bin$SampDate_Color[order(clr.sulf.trait.table.bin$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Bin ID") + ylab("CLR-Transformed Coverage")

ggsave(soxa.bin.fs,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/FxnScatterplots/SoxA_CLR_Coverage_BinID_fxn.sum.png", width=12, height=10, dpi=600)

#`soxD; S-disulfanyl-L-cysteine oxidoreductase SoxD [EC:1.8.2.6]`
soxd.bin.fs<-ggplot(clr.sulf.trait.table.bin, aes(x=PlotBin, y=`soxD; S-disulfanyl-L-cysteine oxidoreductase SoxD [EC:1.8.2.6]`,color=SampDate,group=SampDate)) + geom_jitter(aes(color=SampDate), size=3, width=0.15, height=0) + theme_bw()+
  labs(title="SoxD Depth of Coverage in MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table.bin$SampDate_Color[order(clr.sulf.trait.table.bin$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Bin ID") + ylab("CLR-Transformed Coverage")

ggsave(soxd.bin.fs,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/FxnScatterplots/SoxD_CLR_Coverage_BinID_fxn.sum.png", width=12, height=10, dpi=600)

# `soxB; S-sulfosulfanyl-L-cysteine sulfohydrolase [EC:3.1.6.20]`
soxb.bin.fs<-ggplot(clr.sulf.trait.table.bin, aes(x=PlotBin, y=`soxB; S-sulfosulfanyl-L-cysteine sulfohydrolase [EC:3.1.6.20]`,color=SampDate,group=SampDate)) + geom_jitter(aes(color=SampDate), size=3, width=0.15, height=0) + theme_bw()+
  labs(title="SoxB Depth of Coverage in MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table.bin$SampDate_Color[order(clr.sulf.trait.table.bin$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Bin ID") + ylab("CLR-Transformed Coverage")

ggsave(soxb.bin.fs,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/FxnScatterplots/SoxB_CLR_Coverage_BinID_fxn.sum.png", width=12, height=10, dpi=600)

#### H2S --> S oxidation genes
# `sqr; sulfide:quinone oxidoreductase [EC:1.8.5.4]`
sqr.bin.fs<-ggplot(clr.sulf.trait.table.bin, aes(x=PlotBin, y=`sqr; sulfide:quinone oxidoreductase [EC:1.8.5.4]`,color=SampDate,group=SampDate)) + geom_jitter(aes(color=SampDate), size=3, width=0.15, height=0) + theme_bw()+
  labs(title="sqr Depth of Coverage in MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table.bin$SampDate_Color[order(clr.sulf.trait.table.bin$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Bin ID") + ylab("CLR-Transformed Coverage")

ggsave(sqr.bin.fs,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/FxnScatterplots/sqr_CLR_Coverage_BinID_fxn.sum.png", width=12, height=10, dpi=600)

# `fccB; sulfide dehydrogenase [flavocytochrome c] flavoprotein chain [EC:1.8.2.3]`
fccB.bin.fs<-ggplot(clr.sulf.trait.table.bin, aes(x=PlotBin, y=`fccB; sulfide dehydrogenase [flavocytochrome c] flavoprotein chain [EC:1.8.2.3]`,color=SampDate,group=SampDate)) + geom_jitter(aes(color=SampDate), size=3, width=0.15, height=0) + theme_bw()+
  labs(title="fccB Depth of Coverage in MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table.bin$SampDate_Color[order(clr.sulf.trait.table.bin$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Bin ID") + ylab("CLR-Transformed Coverage")

ggsave(fccB.bin.fs,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/FxnScatterplots/fccB_CLR_Coverage_BinID_fxn.sum.png", width=12, height=10, dpi=600)

#### Dissimilatory SO4 RedOx genes
# `aprA; adenylylsulfate reductase, subunit A [EC:1.8.99.2]`
# aprA.bin.fs<-ggplot(clr.sulf.trait.table.bin, aes(x=PlotBin, y=`aprA; adenylylsulfate reductase, subunit A [EC:1.8.99.2]`,color=SampDate,group=SampDate)) + geom_jitter(aes(color=SampDate), size=3, width=0.15, height=0) + theme_bw()+
#   labs(title="aprA Depth of Coverage in MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
#   guides(shape = guide_legend(override.aes = list(size = 5)))+
#   scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table.bin$SampDate_Color[order(clr.sulf.trait.table.bin$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
#   xlab("Bin ID") + ylab("CLR-Transformed Coverage")
#
# ggsave(aprA.bin.fs,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/FxnScatterplots/aprA_CLR_Coverage_BinID_fxn.sum.png", width=12, height=10, dpi=600)

# `dsrB; dissimilatory sulfite reductase beta subunit [EC:1.8.99.5]`
dsrB.bin.fs<-ggplot(clr.sulf.trait.table.bin, aes(x=PlotBin, y=`dsrB; dissimilatory sulfite reductase beta subunit [EC:1.8.99.5]`,color=SampDate,group=SampDate)) + geom_jitter(aes(color=SampDate), size=3, width=0.15, height=0) + theme_bw()+
  labs(title="dsrB Depth of Coverage in MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table.bin$SampDate_Color[order(clr.sulf.trait.table.bin$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Bin ID") + ylab("CLR-Transformed Coverage")

ggsave(dsrB.bin.fs,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/FxnScatterplots/dsrB_CLR_Coverage_BinID_fxn.sum.png", width=12, height=10, dpi=600)

#### Assimilatory SO4 Reduction

# `cysNC; bifunctional enzyme CysN/CysC [EC:2.7.7.4 2.7.1.25]`
cysNC.bin.fs<-ggplot(clr.sulf.trait.table.bin, aes(x=PlotBin, y=`cysNC; bifunctional enzyme CysN/CysC [EC:2.7.7.4 2.7.1.25]`,color=SampDate,group=SampDate)) + geom_jitter(aes(color=SampDate), size=3, width=0.15, height=0) + theme_bw()+
  labs(title="CysNC Depth of Coverage in MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table.bin$SampDate_Color[order(clr.sulf.trait.table.bin$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Bin ID") + ylab("CLR-Transformed Coverage")

ggsave(cysNC.bin.fs,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/FxnScatterplots/CysNC_CLR_Coverage_BinID_fxn.sum.png", width=12, height=10, dpi=600)

# `cysD; sulfate adenylyltransferase subunit 2 [EC:2.7.7.4]`
cysD.bin.fs<-ggplot(clr.sulf.trait.table.bin, aes(x=PlotBin, y=`cysD; sulfate adenylyltransferase subunit 2 [EC:2.7.7.4]`,color=SampDate,group=SampDate)) + geom_jitter(aes(color=SampDate), size=3, width=0.15, height=0) + theme_bw()+
  labs(title="CysD Depth of Coverage in MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table.bin$SampDate_Color[order(clr.sulf.trait.table.bin$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Bin ID") + ylab("CLR-Transformed Coverage")

ggsave(cysD.bin.fs,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/FxnScatterplots/CysD_CLR_Coverage_BinID_fxn.sum.png", width=12, height=10, dpi=600)

# `cysH; phosphoadenosine phosphosulfate reductase [EC:1.8.4.8 1.8.4.10]`
cysH.bin.fs<-ggplot(clr.sulf.trait.table.bin, aes(x=PlotBin, y=`cysH; phosphoadenosine phosphosulfate reductase [EC:1.8.4.8 1.8.4.10]`,color=SampDate,group=SampDate)) + geom_jitter(aes(color=SampDate), size=3, width=0.15, height=0) + theme_bw()+
  labs(title="CysH Depth of Coverage in MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table.bin$SampDate_Color[order(clr.sulf.trait.table.bin$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Bin ID") + ylab("CLR-Transformed Coverage")

ggsave(cysH.bin.fs,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/FxnScatterplots/CysH_CLR_Coverage_BinID_fxn.sum.png", width=12, height=10, dpi=600)

# `sir; sulfite reductase (ferredoxin) [EC:1.8.7.1]`
# sir.bin.fs<-ggplot(clr.sulf.trait.table.bin, aes(x=PlotBin, y=`sir; sulfite reductase (ferredoxin) [EC:1.8.7.1]`,color=SampDate,group=SampDate)) + geom_jitter(aes(color=SampDate), size=3, width=0.15, height=0) + theme_bw()+
#   labs(title="sir Depth of Coverage in MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
#   guides(shape = guide_legend(override.aes = list(size = 5)))+
#   scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table.bin$SampDate_Color[order(clr.sulf.trait.table.bin$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
#   xlab("Bin ID") + ylab("CLR-Transformed Coverage")
#
# ggsave(sir.bin.fs,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/FxnScatterplots/sir_CLR_Coverage_BinID_fxn.sum.png", width=12, height=10, dpi=600)

#### S Disproportionation
#
# # `phsA,psrA; thiosulfate reductase / polysulfide reductase chain A [EC:1.8.5.5]`
# phsA.psrA.bin.fs<-ggplot(clr.sulf.trait.table.bin, aes(x=PlotBin, y=`phsA,psrA; thiosulfate reductase / polysulfide reductase chain A [EC:1.8.5.5]`,color=SampDate,group=SampDate)) + geom_jitter(aes(color=SampDate), size=3, width=0.15, height=0) + theme_bw()+
#   labs(title="phsA Depth of Coverage in MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
#   guides(shape = guide_legend(override.aes = list(size = 5)))+
#   scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table.bin$SampDate_Color[order(clr.sulf.trait.table.bin$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
#   xlab("Bin ID") + ylab("CLR-Transformed Coverage")
#
# ggsave(phsA.psrA.bin.fs,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/FxnScatterplots/phsA.psrA_CLR_Coverage_BinID_fxn.sum.png", width=12, height=10, dpi=600)


## By Genus

# Note: not looking at every S cycling gene included in this project but looking at ones that appear to have noticeable trends in heat maps
### SOX genes
# `soxY; sulfur-oxidizing protein SoxY`
soxy.bin.fs1<-ggplot(clr.sulf.trait.table.bin, aes(x=Genus, y=`soxY; sulfur-oxidizing protein SoxY`,color=SampDate,group=SampDate)) + geom_jitter(aes(color=SampDate), size=3, width=0.15, height=0) + theme_bw()+
  labs(title="SoxY Depth of Coverage in MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table.bin$SampDate_Color[order(clr.sulf.trait.table.bin$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Genus") + ylab("CLR-Transformed Coverage")

ggsave(soxy.bin.fs1,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/FxnScatterplots/SoxY_CLR_Coverage_Genus_fxn.sum.png", width=12, height=10, dpi=600)

# `soxA; L-cysteine S-thiosulfotransferase [EC:2.8.5.2]`
soxa.bin.fs1<-ggplot(clr.sulf.trait.table.bin, aes(x=Genus, y=`soxA; L-cysteine S-thiosulfotransferase [EC:2.8.5.2]`,color=SampDate,group=SampDate)) + geom_jitter(aes(color=SampDate), size=3, width=0.15, height=0) + theme_bw()+
  labs(title="SoxA Depth of Coverage in MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table.bin$SampDate_Color[order(clr.sulf.trait.table.bin$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Genus") + ylab("CLR-Transformed Coverage")

ggsave(soxa.bin.fs1,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/FxnScatterplots/SoxA_CLR_Coverage_Genus_fxn.sum.png", width=12, height=10, dpi=600)

#`soxD; S-disulfanyl-L-cysteine oxidoreductase SoxD [EC:1.8.2.6]`
soxd.bin.fs1<-ggplot(clr.sulf.trait.table.bin, aes(x=Genus, y=`soxD; S-disulfanyl-L-cysteine oxidoreductase SoxD [EC:1.8.2.6]`,color=SampDate,group=SampDate)) + geom_jitter(aes(color=SampDate), size=3, width=0.15, height=0) + theme_bw()+
  labs(title="SoxD Depth of Coverage in MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table.bin$SampDate_Color[order(clr.sulf.trait.table.bin$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Genus") + ylab("CLR-Transformed Coverage")

ggsave(soxd.bin.fs1,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/FxnScatterplots/SoxD_CLR_Coverage_Genus_fxn.sum.png", width=12, height=10, dpi=600)

# `soxB; S-sulfosulfanyl-L-cysteine sulfohydrolase [EC:3.1.6.20]`
soxb.bin.fs1<-ggplot(clr.sulf.trait.table.bin, aes(x=Genus, y=`soxB; S-sulfosulfanyl-L-cysteine sulfohydrolase [EC:3.1.6.20]`,color=SampDate,group=SampDate)) + geom_jitter(aes(color=SampDate), size=3, width=0.15, height=0) + theme_bw()+
  labs(title="SoxB Depth of Coverage in MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table.bin$SampDate_Color[order(clr.sulf.trait.table.bin$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Genus") + ylab("CLR-Transformed Coverage")

ggsave(soxb.bin.fs1,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/FxnScatterplots/SoxB_CLR_Coverage_Genus_fxn.sum.png", width=12, height=10, dpi=600)

#### H2S --> S oxidation genes
# `sqr; sulfide:quinone oxidoreductase [EC:1.8.5.4]`
sqr.bin.fs1<-ggplot(clr.sulf.trait.table.bin, aes(x=Genus, y=`sqr; sulfide:quinone oxidoreductase [EC:1.8.5.4]`,color=SampDate,group=SampDate)) + geom_jitter(aes(color=SampDate), size=3, width=0.15, height=0) + theme_bw()+
  labs(title="sqr Depth of Coverage in MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table.bin$SampDate_Color[order(clr.sulf.trait.table.bin$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Genus") + ylab("CLR-Transformed Coverage")

ggsave(sqr.bin.fs1,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/FxnScatterplots/sqr_CLR_Coverage_Genus_fxn.sum.png", width=12, height=10, dpi=600)

# `fccB; sulfide dehydrogenase [flavocytochrome c] flavoprotein chain [EC:1.8.2.3]`
fccB.bin.fs1<-ggplot(clr.sulf.trait.table.bin, aes(x=Genus, y=`fccB; sulfide dehydrogenase [flavocytochrome c] flavoprotein chain [EC:1.8.2.3]`,color=SampDate,group=SampDate)) + geom_jitter(aes(color=SampDate), size=3, width=0.15, height=0) + theme_bw()+
  labs(title="fccB Depth of Coverage in MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table.bin$SampDate_Color[order(clr.sulf.trait.table.bin$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Genus") + ylab("CLR-Transformed Coverage")

ggsave(fccB.bin.fs1,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/FxnScatterplots/fccB_CLR_Coverage_Genus_fxn.sum.png", width=12, height=10, dpi=600)

#### Dissimilatory SO4 RedOx genes
# `aprA; adenylylsulfate reductase, subunit A [EC:1.8.99.2]`
# aprA.bin.fs1<-ggplot(clr.sulf.trait.table.bin, aes(x=Genus, y=`aprA; adenylylsulfate reductase, subunit A [EC:1.8.99.2]`,color=SampDate,group=SampDate)) + geom_jitter(aes(color=SampDate), size=3, width=0.15, height=0) + theme_bw()+
#   labs(title="aprA Depth of Coverage in MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
#   guides(shape = guide_legend(override.aes = list(size = 5)))+
#   scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table.bin$SampDate_Color[order(clr.sulf.trait.table.bin$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
#   xlab("Genus") + ylab("CLR-Transformed Coverage")
#
# ggsave(aprA.bin.fs1,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/FxnScatterplots/aprA_CLR_Coverage_Genus_fxn.sum.png", width=12, height=10, dpi=600)

# `dsrB; dissimilatory sulfite reductase beta subunit [EC:1.8.99.5]`
dsrB.bin.fs1<-ggplot(clr.sulf.trait.table.bin, aes(x=Genus, y=`dsrB; dissimilatory sulfite reductase beta subunit [EC:1.8.99.5]`,color=SampDate,group=SampDate)) + geom_jitter(aes(color=SampDate), size=3, width=0.15, height=0) + theme_bw()+
  labs(title="dsrB Depth of Coverage in MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table.bin$SampDate_Color[order(clr.sulf.trait.table.bin$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Genus") + ylab("CLR-Transformed Coverage")

ggsave(dsrB.bin.fs1,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/FxnScatterplots/dsrB_CLR_Coverage_Genus_fxn.sum.png", width=12, height=10, dpi=600)

#### Assimilatory SO4 Reduction

# `cysNC; bifunctional enzyme CysN/CysC [EC:2.7.7.4 2.7.1.25]`
cysNC.bin.fs1<-ggplot(clr.sulf.trait.table.bin, aes(x=Genus, y=`cysNC; bifunctional enzyme CysN/CysC [EC:2.7.7.4 2.7.1.25]`,color=SampDate,group=SampDate)) + geom_jitter(aes(color=SampDate), size=3, width=0.15, height=0) + theme_bw()+
  labs(title="CysNC Depth of Coverage in MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table.bin$SampDate_Color[order(clr.sulf.trait.table.bin$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Genus") + ylab("CLR-Transformed Coverage")

ggsave(cysNC.bin.fs1,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/FxnScatterplots/CysNC_CLR_Coverage_Genus_fxn.sum.png", width=12, height=10, dpi=600)

# `cysD; sulfate adenylyltransferase subunit 2 [EC:2.7.7.4]`
cysD.bin.fs1<-ggplot(clr.sulf.trait.table.bin, aes(x=Genus, y=`cysD; sulfate adenylyltransferase subunit 2 [EC:2.7.7.4]`,color=SampDate,group=SampDate)) + geom_jitter(aes(color=SampDate), size=3, width=0.15, height=0) + theme_bw()+
  labs(title="CysD Depth of Coverage in MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table.bin$SampDate_Color[order(clr.sulf.trait.table.bin$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Genus") + ylab("CLR-Transformed Coverage")

ggsave(cysD.bin.fs1,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/FxnScatterplots/CysD_CLR_Coverage_Genus_fxn.sum.png", width=12, height=10, dpi=600)

# `cysH; phosphoadenosine phosphosulfate reductase [EC:1.8.4.8 1.8.4.10]`
cysH.bin.fs1<-ggplot(clr.sulf.trait.table.bin, aes(x=Genus, y=`cysH; phosphoadenosine phosphosulfate reductase [EC:1.8.4.8 1.8.4.10]`,color=SampDate,group=SampDate)) + geom_jitter(aes(color=SampDate), size=3, width=0.15, height=0) + theme_bw()+
  labs(title="CysH Depth of Coverage in MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table.bin$SampDate_Color[order(clr.sulf.trait.table.bin$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Genus") + ylab("CLR-Transformed Coverage")

ggsave(cysH.bin.fs1,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/FxnScatterplots/CysH_CLR_Coverage_Genus_fxn.sum.png", width=12, height=10, dpi=600)

# `sir; sulfite reductase (ferredoxin) [EC:1.8.7.1]`
# sir.bin.fs1<-ggplot(clr.sulf.trait.table.bin, aes(x=Genus, y=`sir; sulfite reductase (ferredoxin) [EC:1.8.7.1]`,color=SampDate,group=SampDate)) + geom_jitter(aes(color=SampDate), size=3, width=0.15, height=0) + theme_bw()+
#   labs(title="sir Depth of Coverage in MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
#   guides(shape = guide_legend(override.aes = list(size = 5)))+
#   scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table.bin$SampDate_Color[order(clr.sulf.trait.table.bin$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
#   xlab("Genus") + ylab("CLR-Transformed Coverage")
#
# ggsave(sir.bin.fs1,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/FxnScatterplots/sir_CLR_Coverage_Genus_fxn.sum.png", width=12, height=10, dpi=600)

#### S Disproportionation

# `phsA,psrA; thiosulfate reductase / polysulfide reductase chain A [EC:1.8.5.5]`
# phsA.psrA.bin.fs1<-ggplot(clr.sulf.trait.table.bin, aes(x=Genus, y=`phsA,psrA; thiosulfate reductase / polysulfide reductase chain A [EC:1.8.5.5]`,color=SampDate,group=SampDate)) + geom_jitter(aes(color=SampDate), size=3, width=0.15, height=0) + theme_bw()+
#   labs(title="phsA Depth of Coverage in MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
#   guides(shape = guide_legend(override.aes = list(size = 5)))+
#   scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table.bin$SampDate_Color[order(clr.sulf.trait.table.bin$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
#   xlab("Genus") + ylab("CLR-Transformed Coverage")
#
# ggsave(phsA.psrA.bin.fs1,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/FxnScatterplots/phsA.psrA_CLR_Coverage_Genus_fxn.sum.png", width=12, height=10, dpi=600)


#### Pull out Sulfur Metabolic Fxns from Binary Data ####

sulf.ko.bin.bi<-bin_fxn.binary[,which(colnames(bin_fxn.binary) %in% sulfur.fxns.bins$KO_ID)] # merge CLR data w/ S fxns found in contigs from KOFamScan
sulf.ko.bin.bi$Bin_ID<-rownames(sulf.ko.bin.bi)
sulf.ko.bin.bi.melt<-melt(sulf.ko.bin.bi, by="Bin_ID")
colnames(sulf.ko.bin.bi.melt)[which(names(sulf.ko.bin.bi.melt) == "variable")] <- "KO_ID"
colnames(sulf.ko.bin.bi.melt)[which(names(sulf.ko.bin.bi.melt) == "value")] <- "PresAb"
head(sulf.ko.bin.bi.melt) #sanity check

clr.sulf.ko.bin.bi<-merge(sulf.ko.bin.bi.melt,sulf.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
head(clr.sulf.ko.bin.bi)
colnames(clr.sulf.ko.bin.bi)[which(names(clr.sulf.ko.bin.bi) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
clr.cov.sum.sulf.ko.bin.bi<-as.data.frame(dcast(clr.sulf.ko.bin.bi, Bin_ID~KO_Function.KEGG, value.var="PresAb", fun.aggregate=sum)) ###just dcast, nothing is being added here!
rownames(clr.cov.sum.sulf.ko.bin.bi)<-clr.cov.sum.sulf.ko.bin.bi$Bin_ID
clr.cov.sum.sulf.ko.bin.bi[1:4,]

# sanity check
clr.cov.sum.sulf.ko.bin.bi$`cysH; phosphoadenosine phosphosulfate reductase [EC:1.8.4.8 1.8.4.10]`[1:4]
head(clr.sulf.ko.bin.bi)

#### Sulfur Binary Heat Maps ####
# see max & mean of summed

# prep for ggplot2 heatmap
clr.sulf.ko.bin.bi[1:4,]

clr.sulf.all.bin.bi1<-merge(clr.sulf.ko.bin.bi,bin_meta_scaled,by="Bin_ID")
clr.sulf.all.bin.bi<-merge(clr.sulf.all.bin.bi1,mag_tax,by="Bin_ID")

head(clr.sulf.ko.bin.bi)
clr.sulf.all.bin.bi$PlotBin = factor(clr.sulf.all.bin.bi$PlotBin, levels=unique(clr.sulf.all.bin.bi$PlotBin[order(clr.sulf.all.bin.bi$SampDate,clr.sulf.all.bin.bi$Depth_m)]), ordered=TRUE)
clr.sulf.all.bin.bi$SampDate<-gsub("\\."," ",clr.sulf.all.bin.bi$SampDate)
clr.sulf.all.bin.bi$SampDate<-factor(clr.sulf.all.bin.bi$SampDate, levels=c("August 2021","December 2021","April 2022"))
clr.sulf.all.bin.bi$KO_Function.KEGG = factor(clr.sulf.all.bin.bi$KO_Function.KEGG, levels=unique(clr.sulf.all.bin.bi$KO_Function.KEGG[order(clr.sulf.all.bin.bi$Pathway)]), ordered=TRUE)

clr.sulf.all.bin.bi$PathShort<-clr.sulf.all.bin.bi$Pathway
clr.sulf.all.bin.bi$PathShort[(clr.sulf.all.bin.bi$PathShort) == "Dissimilatory Sulfate Redox"] <- "D.SO4 RedOx"
clr.sulf.all.bin.bi$PathShort[(clr.sulf.all.bin.bi$PathShort) == "Assimilatory Sulfate Reduction"] <- "A.SO4 Red"
clr.sulf.all.bin.bi$PathShort[(clr.sulf.all.bin.bi$PathShort) == "Multiple Pathways"] <- "Multi Paths"
clr.sulf.all.bin.bi$PathShort[(clr.sulf.all.bin.bi$PathShort) == "S Disproportionation"] <- "S Disprop."

clr.sulf.all.bin.bi$Pathway<-factor(clr.sulf.all.bin.bi$Pathway,levels=c("Assimilatory Sulfate Reduction","Dissimilatory Sulfate Redox","Multiple Pathways","SOX","S Disproportionation"))
clr.sulf.all.bin.bi$PathShort<-factor(clr.sulf.all.bin.bi$PathShort,levels=c("A.SO4 Red","D.SO4 RedOx","Multi Paths","SOX","S Disprop."))

clr.sulf.all.bin.bi$PathSpecShort<-clr.sulf.all.bin.bi$PathwaySpecific
clr.sulf.all.bin.bi$PathSpecShort[(clr.sulf.all.bin.bi$PathSpecShort) == "Dissimilatory Sulfate Redox"] <- "1"
clr.sulf.all.bin.bi$PathSpecShort[(clr.sulf.all.bin.bi$PathSpecShort) == "Assimilatory Sulfate Reduction"] <- "2"
clr.sulf.all.bin.bi$PathSpecShort[(clr.sulf.all.bin.bi$PathSpecShort) == "Multiple Pathways"] <- "3"
clr.sulf.all.bin.bi$PathSpecShort[(clr.sulf.all.bin.bi$PathSpecShort) == "Sulfur Disproportionation"] <- "4"
clr.sulf.all.bin.bi$PathSpecShort[(clr.sulf.all.bin.bi$PathSpecShort) == "Sulfide Oxidation"] <- "5"
clr.sulf.all.bin.bi$PathSpecShort[(clr.sulf.all.bin.bi$PathSpecShort) == "Sulfite Oxidation"] <- "6"
clr.sulf.all.bin.bi$PathSpecShort[(clr.sulf.all.bin.bi$PathSpecShort) == "Thiosulfate Oxidation"] <- "7"

clr.sulf.all.bin.bi$PathwaySpecific<-factor(clr.sulf.all.bin.bi$PathwaySpecific,levels=c("Assimilatory Sulfate Reduction","Dissimilatory Sulfate Redox","Multiple Pathways","SOX","S Disproportionation","Sulfide Oxidation","Sulfite Oxidation","Thiosulfate Oxidation"))
clr.sulf.all.bin.bi$PathSpecShort<-factor(clr.sulf.all.bin.bi$PathSpecShort,levels=c("1","2","3","4","5","6","7"))

head(clr.sulf.all.bin.bi)

# Figures
sulf.bi.hm1a<-ggplot(clr.sulf.all.bin.bi, aes(PlotBin, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Sulfur Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(sulf.bi.hm1a,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/PresenceAbsence/Sulfur_KOFxns_MGMs_BinID_by_Function_Binary_heatmap.png", width=18, height=13, dpi=600)

sulf.bi.hm1b<-ggplot(clr.sulf.all.bin.bi, aes(PlotBin, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Sulfur Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~.,scales="free_y", space = "free")

ggsave(sulf.bi.hm1b,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/PresenceAbsence/Sulfur_KOFxns_MGMs_BinID_by_Function_Pathway_Binary_heatmap.png", width=17, height=15, dpi=600)

sulf.bi.hm1b2<-ggplot(clr.sulf.all.bin.bi, aes(Genus, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Sulfur Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~.,scales="free_y", space = "free")

ggsave(sulf.bi.hm1b2,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/PresenceAbsence/Sulfur_KOFxns_MGMs_Genus_by_Function_Pathway_Binary_heatmap.png", width=17, height=15, dpi=600)

sulf.bi.hm1b3<-ggplot(clr.sulf.all.bin.bi[clr.sulf.all.bin.bi$Genus=="HIMB30",], aes(PlotBin, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Sulfur Fixation in HIMB30 MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~.,scales="free_y", space = "free")

ggsave(sulf.bi.hm1b3,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/PresenceAbsence/Sulfur_KOFxns_MGMs_HIMB30_Bins_Only_by_Function_Pathway_Binary_heatmap.png", width=17, height=15, dpi=600)

sulf.bi.hm1b3<-ggplot(clr.sulf.all.bin.bi, aes(PlotBin, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Sulfur Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~.,scales="free_y", space = "free")

ggsave(sulf.bi.hm1b3,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/PresenceAbsence/Sulfur_KOFxns_MGMs_Genus_by_Function_Pathway_Binary_heatmap.png", width=17, height=15, dpi=600)

sulf.bi.hm1c2<-ggplot(clr.sulf.all.bin.bi, aes(PlotBin, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Sulfur Metabolism in Salton Seawater MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathwaySpecific~SampDate, scales="free", space = "free")

ggsave(sulf.bi.hm1c2,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/PresenceAbsence/Sulfur_KOFxns_MGMs_Bins_BinID_by_Function_SampDate_PathwaySpecific_Binary_best_heatmap.png", width=20, height=20, dpi=600)

sulf.bi.hm1c3<-ggplot(clr.sulf.all.bin.bi, aes(Genus, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Sulfur Metabolism in Salton Seawater MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathwaySpecific~SampDate, scales="free", space = "free")

ggsave(sulf.bi.hm1c3,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/PresenceAbsence/Sulfur_KOFxns_MGMs_Bins_Genus_by_Function_SampDate_PathwaySpecific_Binary_best_heatmap.png", width=20, height=20, dpi=600)

sulf.bi.hm1c4<-ggplot(clr.sulf.all.bin.bi[clr.sulf.all.bin.bi$Genus=="HIMB30",], aes(PlotBin, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Sulfur Fixation in HIMB30 MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathwaySpecific~.,scales="free_y", space = "free")

ggsave(sulf.bi.hm1c4,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/PresenceAbsence/Sulfur_KOFxns_MGMs_HIMB30_Bins_Only_by_Function_PathwaySpecific_Binary_heatmap.png", width=17, height=15, dpi=600)

sulf.bi.hm1d<-ggplot(clr.sulf.all.bin.bi, aes(PlotBin, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Sulfur Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate,scales="free_x", space = "free")

ggsave(sulf.bi.hm1d,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/PresenceAbsence/Sulfur_KOFxns_MGMs_BinID_by_Function_SampDate_Binary_best_heatmap.png", width=20, height=13, dpi=600)

sulf.bi.hm1d2<-ggplot(clr.sulf.all.bin.bi, aes(Genus, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Sulfur Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate,scales="free_x", space = "free")

ggsave(sulf.bi.hm1d2,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/PresenceAbsence/Sulfur_KOFxns_MGMs_Genus_by_Function_SampDate_Binary_best_heatmap.png", width=20, height=13, dpi=600)

sulf.bi.hm1e<-ggplot(clr.sulf.all.bin.bi, aes(PlotBin, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Sulfur Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~SampDate, scales="free", space = "free")

ggsave(sulf.bi.hm1e,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/PresenceAbsence/Sulfur_KOFxns_MGMs_BinID_by_Function_SampDate_Pathway_Binary_best_heatmap.png", width=20, height=15, dpi=600)

sulf.bi.hm1e2<-ggplot(clr.sulf.all.bin.bi, aes(Genus, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Sulfur Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~SampDate, scales="free", space = "free")

ggsave(sulf.bi.hm1e,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/PresenceAbsence/Sulfur_KOFxns_MGMs_Genus_by_Function_SampDate_Pathway_Binary_best_heatmap.png", width=20, height=15, dpi=600)

sulf.bi.hm1e2<-ggplot(clr.sulf.all.bin.bi, aes(Depth_m, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Sulfur Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~SampDate, scales="free", space = "free")

ggsave(sulf.bi.hm1e,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/PresenceAbsence/Sulfur_KOFxns_MGMs_Depth_by_Function_SampDate_Pathway_Binary_best_heatmap.png", width=20, height=15, dpi=600)

# sulf.bi.hm1e0<-ggplot(clr.sulf.all.bin.bi[clr.sulf.all.bin.bi$Depth_m==0,], aes(PathShort, KO_Function.KEGG, fill=factor(PresAb))) +
#   geom_tile(colour="black",size=0.25) +
#   scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Sulfur Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate)
#
# ggsave(sulf.bi.hm1e0,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/PresenceAbsence/Sulfur_KOFxns_Bins_Pathways_Binary_0m_heatmap.png", width=18, height=18, dpi=600)
#
# sulf.bi.hm1e5<-ggplot(clr.sulf.all.bin.bi[clr.sulf.all.bin.bi$Depth_m==5,], aes(PathShort, KO_Function.KEGG, fill=factor(PresAb))) +
#   geom_tile(colour="black",size=0.25) +
#   scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Sulfur Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)
#
# ggsave(sulf.bi.hm1e5,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/PresenceAbsence/Sulfur_KOFxns_Bins_Pathways_Binary_5m_heatmap.png", width=18, height=18, dpi=600)
#
# sulf.bi.hm1e6<-ggplot(clr.sulf.all.bin.bi[clr.sulf.all.bin.bi$Depth_m==10,], aes(PathShort, KO_Function.KEGG, fill=factor(PresAb))) +
#   geom_tile(colour="black",size=0.25) +
#   scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Sulfur Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)
#
# ggsave(sulf.bi.hm1e6,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Sulfur/PresenceAbsence/Sulfur_KOFxns_Bins_Pathways_Binary_10m_heatmap.png", width=18, height=18, dpi=600)


#### Pull Out Nitrogen Metabolic Fxns from CLR data - with NAs ####
## heatmaps of traits of interest

bin.clr.na[1:4,1:4]

# pull out nitro functions from CLR transformed, summed coverages (summed coverage per KO)
nitro.ko.na.bin<-bin.clr.na[,which(colnames(bin.clr.na) %in% nitro.fxns.bins$KO_ID)] # merge CLR data w/ S fxns found in contigs from KOFamScan
dim(nitro.ko.na.bin) # 35 total bins, 3 columns aka functions
nitro.ko.na.bin<-nitro.ko.na.bin[rowSums(is.na(nitro.ko.na.bin)) != ncol(nitro.ko.na.bin), ] # drop all rows that contain ONLY NAs
dim(nitro.ko.na.bin) # only 4 bins have these functions...
# is.na() tells us which elements are NA (TRUE) or not NA (FALSE)
# rowSums(is.na(df)) tells us how many columns have NAs; in this df, there are a total of 14 columns
# if rowSums(is.na(df)) != total # of rows, aka if the total columns per row with NAs is LESS than the total # of columns per row, then keep the row
# ^ if the row contains only NAs, and rowSums(is.na(df)) is NOT 14 (total # of columns), then we keep it

nitro.ko.na.bin$Bin_ID<-rownames(nitro.ko.na.bin)
nitro.ko.na.bin.melt<-melt(nitro.ko.na.bin, by="Bin_ID")
colnames(nitro.ko.na.bin.melt)[which(names(nitro.ko.na.bin.melt) == "variable")] <- "KO_ID"
colnames(nitro.ko.na.bin.melt)[which(names(nitro.ko.na.bin.melt) == "value")] <- "CLR_SumCovPerKO"
head(nitro.ko.na.bin.melt) #sanity check

clr.nitro.ko.na.bin<-merge(nitro.ko.na.bin.melt,nitro.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
head(clr.nitro.ko.na.bin)
colnames(clr.nitro.ko.na.bin)[which(names(clr.nitro.ko.na.bin) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
clr.cov.sum.nitro.ko.na.bin<-as.data.frame(dcast(clr.nitro.ko.na.bin, Bin_ID~KO_Function.KEGG, value.var="CLR_SumCovPerKO", fun.aggregate=sum)) ###just dcast, nothing is being added here!
rownames(clr.cov.sum.nitro.ko.na.bin)<-clr.cov.sum.nitro.ko.na.bin$Bin_ID
clr.cov.sum.nitro.ko.na.bin[1:4,]

# sanity check
clr.cov.sum.nitro.ko.na.bin$`nirK; nitrite reductase (NO-forming) [EC:1.7.2.1]`[1:4]
head(clr.nitro.ko.na.bin)

#### Nitrogen Heat Maps ####
# see max & mean of summed
max(clr.nitro.ko.na.bin$CLR_SumCovPerKO,na.rm=TRUE)
mean(clr.nitro.ko.na.bin$CLR_SumCovPerKO,na.rm=TRUE)

# prep for ggplot2 heatmap
clr.nitro.ko.na.bin[1:4,]

clr.nitro.all.bin1<-merge(clr.nitro.ko.na.bin,bin_meta_scaled,by="Bin_ID")
clr.nitro.all.bin<-merge(clr.nitro.all.bin1,mag_tax,by=c("Bin_ID","PlotBin"))

head(clr.nitro.all.bin)
clr.nitro.all.bin$PlotBin = factor(clr.nitro.all.bin$PlotBin, levels=unique(clr.nitro.all.bin$PlotBin[order(clr.nitro.all.bin$SampDate,clr.nitro.all.bin$Depth_m)]), ordered=TRUE)
clr.nitro.all.bin$SampDate<-gsub("\\."," ",clr.nitro.all.bin$SampDate)
clr.nitro.all.bin$SampDate<-factor(clr.nitro.all.bin$SampDate, levels=c("August 2021","December 2021","April 2022"))
clr.nitro.all.bin$KO_Function.KEGG = factor(clr.nitro.all.bin$KO_Function.KEGG, levels=unique(clr.nitro.all.bin$KO_Function.KEGG[order(clr.nitro.all.bin$Pathway)]), ordered=TRUE)

# create shortened name for pathways
clr.nitro.all.bin$PathShort<-clr.nitro.all.bin$Pathway
# vvv can only do this type of renaming if variables are characters, not factors
clr.nitro.all.bin$PathShort[(clr.nitro.all.bin$PathShort) == "Dissimilatory Nitrate Reduction"] <- "D. NO3 Red"
clr.nitro.all.bin$PathShort[(clr.nitro.all.bin$PathShort) == "Assimilatory Nitrate Reduction"] <- "A. NO3 Red"

# turn pathways & pathshort into factors
clr.nitro.all.bin$Pathway<-factor(clr.nitro.all.bin$Pathway,levels=c("Assimilatory Nitrate Reduction","Dissimilatory Nitrate Reduction","Multiple Pathways","Denitrification","Anammox"))
clr.nitro.all.bin$PathShort<-factor(clr.nitro.all.bin$PathShort,levels=c("A. NO3 Red","D. NO3 Red","Multiple Pathways","Denitrification","Anammox"))

head(clr.nitro.all.bin)

# for heatmap color gradient
max(clr.nitro.all.bin$CLR_SumCovPerKO,na.rm=TRUE)
max(clr.nitro.all.bin$CLR_SumCovPerKO,na.rm=TRUE)/2
min(clr.nitro.all.bin$CLR_SumCovPerKO,na.rm=TRUE)

# Figures
nitro.hm1a<-ggplot(clr.nitro.all.bin, aes(PlotBin, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.4","0.2","0.1"),breaks=c(0.4,0.2,0.1)) + labs(title="Nitrogen Metabolism in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(nitro.hm1a,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Nitrogen/Nitrogen_KOFxns_MGMs_BinID_by_Function_heatmap.png", width=18, height=13, dpi=600)

nitro.hm1b<-ggplot(clr.nitro.all.bin, aes(PlotBin, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.4","0.2","0.1"),breaks=c(0.4,0.2,0.1)) + labs(title="Nitrogen Metabolism in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~.,scales="free_y", space = "free")

ggsave(nitro.hm1b,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Nitrogen/Nitrogen_KOFxns_MGMs_BinID_by_Function_Pathway_heatmap.png", width=17, height=15, dpi=600)

# nitro.hm1c<-ggplot(clr.nitro.all.bin, aes(interaction(SampDate,Depth_m), KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.4","0.2","0.1"),breaks=c(0.4,0.2,0.1)) + labs(title="Nitrogen Metabolism in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~.,scales="free_y", space = "free")
#
# ggsave(nitro.hm1c,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Nitrogen/Nitrogen_KOFxns_MGMs_SampDate_Depth_by_Function_Pathway_heatmap.png", width=15, height=18, dpi=600)

nitro.hm1d<-ggplot(clr.nitro.all.bin, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.4","0.2","0.1"),breaks=c(0.4,0.2,0.1)) + labs(title="Nitrogen Metabolism in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate,scales="free_x", space = "free")

ggsave(nitro.hm1d,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Nitrogen/Nitrogen_KOFxns_MGMs_Depth_by_Function_SampDate_best_heatmap.png", width=20, height=13, dpi=600)

nitro.hm1e<-ggplot(clr.nitro.all.bin, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.4","0.2","0.1"),breaks=c(0.4,0.2,0.1)) + labs(title="Nitrogen Metabolism in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~SampDate, scales="free", space = "free")

ggsave(nitro.hm1e,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Nitrogen/Nitrogen_KOFxns_MGMs_Depth_by_Function_SampDate_Pathway_best_heatmap.png", width=20, height=15, dpi=600)
#
# nitro.hm1f<-ggplot(clr.nitro.all.bin, aes(PathShort, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.4","0.2","0.1"),breaks=c(0.4,0.2,0.1)) + labs(title="Nitrogen Metabolism in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11,face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Depth_m~SampDate,scales="free", space = "free")
#
# ggsave(nitro.hm1f,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Nitrogen/Nitrogen_KOFxns_MGMs_heatmap1d.png", width=18, height=18, dpi=600)
#
# nitro.hm1g<-ggplot(clr.nitro.all.bin, aes(PathShort, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.4","0.2","0.1"),breaks=c(0.4,0.2,0.1)) + labs(title="Nitrogen Metabolism in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11,face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_wrap(.~SampDate)
#
# ggsave(nitro.hm1g,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Nitrogen/Nitrogen_KOFxns_MGMs_heatmap1d.png", width=18, height=18, dpi=600)

nitro.hm1e<-ggplot(clr.nitro.all.bin[clr.nitro.all.bin$Depth_m==0,], aes(PathShort, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.4","0.2","0.1"),breaks=c(0.4,0.2,0.1)) + labs(title="Nitrogen Metabolism in Salton Seawater MAGs - 0m",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)

ggsave(nitro.hm1e,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Nitrogen/Nitrogen_KOFxns_Pathways_MGMs_0m_heatmap.png", width=18, height=18, dpi=600)

nitro.hm1f<-ggplot(clr.nitro.all.bin[clr.nitro.all.bin$Depth_m==5,], aes(PathShort, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.4","0.2","0.1"),breaks=c(0.4,0.2,0.1)) + labs(title="Nitrogen Metabolism in Salton Seawater MAGs - 5m",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)

ggsave(nitro.hm1f,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Nitrogen/Nitrogen_KOFxns_Pathways_MGMs_5m_heatmap.png", width=18, height=18, dpi=600)

nitro.hm1g<-ggplot(clr.nitro.all.bin[clr.nitro.all.bin$Depth_m==10,], aes(PathShort, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.4","0.2","0.1"),breaks=c(0.4,0.2,0.1)) + labs(title="Nitrogen Metabolism in Salton Seawater MAGs - 10m",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)

ggsave(nitro.hm1g,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Nitrogen/Nitrogen_KOFxns_Pathways_MGMs_10m_heatmap.png", width=18, height=18, dpi=600)

#### Pull out Nitrogen Metabolic Fxns from Binary Data ####

nitro.ko.bin.bi<-bin_fxn.binary[,which(colnames(bin_fxn.binary) %in% nitro.fxns.bins$KO_ID)] # merge CLR data w/ N fxns found in contigs from KOFamScan
nitro.ko.bin.bi$Bin_ID<-rownames(nitro.ko.bin.bi)
nitro.ko.bin.bi.melt<-melt(nitro.ko.bin.bi, by="Bin_ID")
colnames(nitro.ko.bin.bi.melt)[which(names(nitro.ko.bin.bi.melt) == "variable")] <- "KO_ID"
colnames(nitro.ko.bin.bi.melt)[which(names(nitro.ko.bin.bi.melt) == "value")] <- "PresAb"
head(nitro.ko.bin.bi.melt) #sanity check

clr.nitro.ko.bin.bi<-merge(nitro.ko.bin.bi.melt,nitro.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
head(clr.nitro.ko.bin.bi)
colnames(clr.nitro.ko.bin.bi)[which(names(clr.nitro.ko.bin.bi) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
clr.cov.sum.nitro.ko.bin.bi<-as.data.frame(dcast(clr.nitro.ko.bin.bi, Bin_ID~KO_Function.KEGG, value.var="PresAb", fun.aggregate=sum)) ###just dcast, nothing is being added here!
rownames(clr.cov.sum.nitro.ko.bin.bi)<-clr.cov.sum.nitro.ko.bin.bi$Bin_ID
clr.cov.sum.nitro.ko.bin.bi[1:4,]

# sanity check
clr.cov.sum.nitro.ko.bin.bi$`nosZ; nitrous-oxide reductase [EC:1.7.2.4]`[1:4]
head(clr.nitro.ko.bin.bi)

#### Nitrogen Binary Heat Maps ####
# see max & mean of summed
max(clr.cov.sum.nitro.ko.bin[,-1])
mean(as.matrix(clr.cov.sum.nitro.ko.bin[,-1]))

# first heat map of nitro KOs
heatmap(as.matrix(clr.cov.sum.nitro.ko.bin[,-1]), scale = "none")

colSums(clr.cov.sum.nitro.ko.bin[,-1])
#clr.cov.sum.nitro.ko.bin2 <- clr.cov.sum.nitro.ko.bin[,which(colSums(clr.cov.sum.nitro.ko.bin[,-1])>10)]

heatmap(as.matrix(clr.cov.sum.nitro.ko.bin[,-1]), scale = "none")

# prep for ggplot2 heatmap
clr.nitro.ko.bin.bi[1:4,]

clr.nitro.all.bin.bi1<-merge(clr.nitro.ko.bin.bi,bin_meta_scaled,by="Bin_ID")
clr.nitro.all.bin.bi<-merge(clr.nitro.all.bin.bi1,mag_tax,by=c("Bin_ID","PlotBin"))

head(clr.nitro.all.bin.bi)
clr.nitro.all.bin.bi$PlotBin = factor(clr.nitro.all.bin.bi$PlotBin, levels=unique(clr.nitro.all.bin.bi$PlotBin[order(clr.nitro.all.bin.bi$SampDate,clr.nitro.all.bin.bi$Depth_m)]), ordered=TRUE)
clr.nitro.all.bin.bi$SampDate<-gsub("\\."," ",clr.nitro.all.bin.bi$SampDate)
clr.nitro.all.bin.bi$SampDate<-factor(clr.nitro.all.bin.bi$SampDate, levels=c("August 2021","December 2021","April 2022"))
clr.nitro.all.bin.bi$KO_Function.KEGG = factor(clr.nitro.all.bin.bi$KO_Function.KEGG, levels=unique(clr.nitro.all.bin.bi$KO_Function.KEGG[order(clr.nitro.all.bin.bi$Pathway)]), ordered=TRUE)

# create shortened name for pathways
clr.nitro.all.bin.bi$PathShort<-clr.nitro.all.bin.bi$Pathway
# vvv can only do this type of renaming if variables are characters, not factors
clr.nitro.all.bin.bi$PathShort[(clr.nitro.all.bin.bi$PathShort) == "Dissimilatory Nitrate Reduction"] <- "D. NO3 Red"
clr.nitro.all.bin.bi$PathShort[(clr.nitro.all.bin.bi$PathShort) == "Assimilatory Nitrate Reduction"] <- "A. NO3 Red"

# turn pathways & pathshort into factors
clr.nitro.all.bin.bi$Pathway<-factor(clr.nitro.all.bin.bi$Pathway,levels=c("Assimilatory Nitrate Reduction","Dissimilatory Nitrate Reduction","Multiple Pathways","Denitrification","Anammox"))
clr.nitro.all.bin.bi$PathShort<-factor(clr.nitro.all.bin.bi$PathShort,levels=c("A. NO3 Red","D. NO3 Red","Multiple Pathways","Denitrification","Anammox"))

head(clr.nitro.all.bin.bi)

# Figures
nitro.bi.hm1a<-ggplot(clr.nitro.all.bin.bi, aes(PlotBin, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Nitrogen Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(nitro.bi.hm1a,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Nitrogen/PresenceAbsence/Nitrogen_KOFxns_MGMs_BinID_by_Function_Binary_heatmap.png", width=18, height=13, dpi=600)

nitro.bi.hm1b<-ggplot(clr.nitro.all.bin.bi, aes(PlotBin, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Nitrogen Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~.,scales="free_y", space = "free")

ggsave(nitro.bi.hm1b,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Nitrogen/PresenceAbsence/Nitrogen_KOFxns_MGMs_BinID_by_Function_Pathway_Binary_heatmap.png", width=17, height=15, dpi=600)

nitro.bi.hm1b2<-ggplot(clr.nitro.all.bin.bi, aes(Genus, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Nitrogen Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~.,scales="free_y", space = "free")

ggsave(nitro.bi.hm1b2,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Nitrogen/PresenceAbsence/Nitrogen_KOFxns_MGMs_Genus_by_Function_Pathway_Binary_heatmap.png", width=17, height=15, dpi=600)

nitro.bi.hm1d<-ggplot(clr.nitro.all.bin.bi, aes(PlotBin, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Nitrogen Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate,scales="free_x", space = "free")

ggsave(nitro.bi.hm1d,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Nitrogen/PresenceAbsence/Nitrogen_KOFxns_MGMs_BinID_by_Function_SampDate_Binary_best_heatmap.png", width=20, height=13, dpi=600)

nitro.bi.hm1d2<-ggplot(clr.nitro.all.bin.bi, aes(Genus, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Nitrogen Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate,scales="free_x", space = "free")

ggsave(nitro.bi.hm1d2,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Nitrogen/PresenceAbsence/Nitrogen_KOFxns_MGMs_Genus_by_Function_SampDate_Binary_best_heatmap.png", width=20, height=13, dpi=600)

nitro.bi.hm1e<-ggplot(clr.nitro.all.bin.bi, aes(PlotBin, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Nitrogen Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~SampDate, scales="free", space = "free")

ggsave(nitro.bi.hm1e,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Nitrogen/PresenceAbsence/Nitrogen_KOFxns_MGMs_BinID_by_Function_SampDate_Pathway_Binary_best_heatmap.png", width=20, height=15, dpi=600)

nitro.bi.hm1e2<-ggplot(clr.nitro.all.bin.bi, aes(Genus, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Nitrogen Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~SampDate, scales="free", space = "free")

ggsave(nitro.bi.hm1e,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Nitrogen/PresenceAbsence/Nitrogen_KOFxns_MGMs_Genus_by_Function_SampDate_Pathway_Binary_best_heatmap.png", width=20, height=15, dpi=600)

nitro.bi.hm1e2<-ggplot(clr.nitro.all.bin.bi, aes(Depth_m, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Nitrogen Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~SampDate, scales="free", space = "free")

ggsave(nitro.bi.hm1e,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Nitrogen/PresenceAbsence/Nitrogen_KOFxns_MGMs_Depth_by_Function_SampDate_Pathway_Binary_best_heatmap.png", width=20, height=15, dpi=600)
#
# nitro.bi.hm1e0<-ggplot(clr.nitro.all.bin.bi[clr.nitro.all.bin.bi$Depth_m==0,], aes(PathShort, KO_Function.KEGG, fill=factor(PresAb))) +
#   geom_tile(colour="black",size=0.25) +
#   scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Nitrogen Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate)
#
# ggsave(nitro.bi.hm1e0,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Nitrogen/PresenceAbsence/Nitrogen_KOFxns_Bins_Pathways_Binary_0m_heatmap.png", width=18, height=18, dpi=600)
#
# nitro.bi.hm1e5<-ggplot(clr.nitro.all.bin.bi[clr.nitro.all.bin.bi$Depth_m==5,], aes(PathShort, KO_Function.KEGG, fill=factor(PresAb))) +
#   geom_tile(colour="black",size=0.25) +
#   scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Nitrogen Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)
#
# ggsave(nitro.bi.hm1e5,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Nitrogen/PresenceAbsence/Nitrogen_KOFxns_Bins_Pathways_Binary_5m_heatmap.png", width=18, height=18, dpi=600)
#
# nitro.bi.hm1e6<-ggplot(clr.nitro.all.bin.bi[clr.nitro.all.bin.bi$Depth_m==10,], aes(PathShort, KO_Function.KEGG, fill=factor(PresAb))) +
#   geom_tile(colour="black",size=0.25) +
#   scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Nitrogen Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)
#
# ggsave(nitro.bi.hm1e6,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Nitrogen/PresenceAbsence/Nitrogen_KOFxns_Bins_Pathways_Binary_10m_heatmap.png", width=18, height=18, dpi=600)


#### Pull Out Carbon Metabolic Fxns from CLR data - with NAs ####
## heatmaps of traits of interest

bin.clr.na[1:4,1:4]

# pull out Carbon metabolism functions from CLR transformed, summed coverages (summed coverage per KO)
carb.ko.na.bin<-bin.clr.na[,which(colnames(bin.clr.na) %in% carb.fxns.bins$KO_ID)] # merge CLR data w/ carbon-related fxns found in contigs from KOFamScan
dim(carb.ko.na.bin) # 35 total bins, 43 columns aka functions
carb.ko.na.bin<-carb.ko.na.bin[rowSums(is.na(carb.ko.na.bin)) != ncol(carb.ko.na.bin), ] # drop all rows that contain ONLY NAs
dim(carb.ko.na.bin) # 31 bins have the 43 functions
# is.na() tells us which elements are NA (TRUE) or not NA (FALSE)
# rowSums(is.na(df)) tells us how many columns have NAs; in this df, there are a total of 43 columns
# if rowSums(is.na(df)) != total # of rows, aka if the total columns per row with NAs is LESS than the total # of columns per row, then keep the row
# ^ if the row contains only NAs, and rowSums(is.na(df)) is NOT 43 (total # of columns), then we keep it

carb.ko.na.bin$Bin_ID<-rownames(carb.ko.na.bin)
carb.ko.na.bin.melt<-melt(carb.ko.na.bin, by="Bin_ID")
colnames(carb.ko.na.bin.melt)[which(names(carb.ko.na.bin.melt) == "variable")] <- "KO_ID"
colnames(carb.ko.na.bin.melt)[which(names(carb.ko.na.bin.melt) == "value")] <- "CLR_SumCovPerKO"
head(carb.ko.na.bin.melt) #sanity check

clr.carb.ko.na.bin<-merge(carb.ko.na.bin.melt,carb.kegg,by.x=c("KO_ID"),by.y=c("KO_ID"))
head(clr.carb.ko.na.bin)
colnames(clr.carb.ko.na.bin)[which(names(clr.carb.ko.na.bin) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
clr.cov.sum.carb.ko.na.bin<-as.data.frame(dcast(clr.carb.ko.na.bin, Bin_ID~KO_Function.KEGG, value.var="CLR_SumCovPerKO", fun.aggregate=sum)) ###
rownames(clr.cov.sum.carb.ko.na.bin)<-clr.cov.sum.carb.ko.na.bin$Bin_ID
clr.cov.sum.carb.ko.na.bin[1:4,1:4]

# sanity check
clr.cov.sum.carb.ko.na.bin$`accB, bccP; acetyl-CoA carboxylase biotin carboxyl carrier protein`[1:4]
head(clr.cov.sum.carb.ko.na.bin)

#### Carbon Heat Maps ####
# see max & mean of summed
max(clr.cov.sum.carb.ko.na.bin[,-1],na.rm = TRUE)
mean(as.matrix(clr.cov.sum.carb.ko.na.bin[,-1]),na.rm = TRUE)

colSums(clr.cov.sum.carb.ko.na.bin[,-1],,na.rm = TRUE)

# prep for ggplot2 heatmap
clr.carb.ko.na.bin[1:4,]
clr.carb.all.bin1<-merge(clr.carb.ko.na.bin,bin_meta_scaled,by="Bin_ID")
clr.carb.all.bin<-merge(clr.carb.all.bin1,mag_tax,by=c("Bin_ID","PlotBin"))

head(clr.carb.all.bin)
#clr.carb.all.bin$PlotBin = factor(clr.carb.all.bin$PlotBin, levels=unique(clr.carb.all.bin$PlotBin[order(clr.carb.all.bin$SampDate,clr.carb.all.bin$Depth_m)]), ordered=TRUE)
clr.carb.all.bin$SampDate<-gsub("\\."," ",clr.carb.all.bin$SampDate)
clr.carb.all.bin$SampDate<-factor(clr.carb.all.bin$SampDate, levels=c("August 2021","December 2021","April 2022"))

unique(clr.carb.all.bin$Pathway)
clr.carb.all.bin<-subset(clr.carb.all.bin, clr.carb.all.bin$Pathway!="Multiple Pathways")
"Multiple Pathways" %in% clr.carb.all.bin$Pathway
clr.carb.all.bin$PathShort<-clr.carb.all.bin$Pathway
clr.carb.all.bin$PathShort[(clr.carb.all.bin$PathShort) == "Reductive Tricarboxylic Acid Cycle"] <- "rTCA"
clr.carb.all.bin$PathShort[(clr.carb.all.bin$PathShort) == "3-Hydroxypropionate Bi-cycle"] <- "3HP"
clr.carb.all.bin$PathShort[(clr.carb.all.bin$PathShort) == "Reductive acetyl-CoA Pathway"] <- "RAcCoa"
clr.carb.all.bin$PathShort[(clr.carb.all.bin$PathShort) == "Calvin Cycle"] <- "CBB"
#clr.carb.all.bin$PathShort[(clr.carb.all.bin$PathShort) == "TCA Cycle"] <- "TCA"


clr.carb.all.bin$Pathway<-factor(clr.carb.all.bin$Pathway,levels=c("3-Hydroxypropionate Bi-cycle","Reductive Tricarboxylic Acid Cycle","Reductive acetyl-CoA Pathway","Calvin Cycle"))
clr.carb.all.bin$PathShort<-factor(clr.carb.all.bin$PathShort,levels=c("3HP","rTCA","RAcCoa","CBB"))

#clr.carb.all.bin$KO_Function.KEGG = factor(clr.carb.all.bin$KO_Function.KEGG, levels=unique(clr.carb.all.bin$KO_Function.KEGG[order(clr.carb.all.bin$Pathway)]), ordered=TRUE)
clr.carb.all.bin$PlotBin = factor(clr.carb.all.bin$PlotBin, levels=unique(clr.carb.all.bin$PlotBin[order(clr.carb.all.bin$SampDate,clr.carb.all.bin$Depth_m)]), ordered=TRUE)

head(clr.carb.all.bin)

# For heatmap color gradient
max(clr.carb.all.bin$CLR_SumCovPerKO,na.rm = TRUE)
max(clr.carb.all.bin$CLR_SumCovPerKO,na.rm = TRUE)/2
min(clr.carb.all.bin$CLR_SumCovPerKO,na.rm = TRUE)

# Figures
carb.hm1a<-ggplot(clr.carb.all.bin, aes(PlotBin, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.8","0.45","0.1"),breaks=c(0.8,0.45,0.1)) + labs(title="Carbon Fixation in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(carb.hm1a,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Carbon/Carbon_KOFxns_MGMs_BinID_by_Function_heatmap.png", width=18, height=16, dpi=600)

carb.hm1b<-ggplot(clr.carb.all.bin, aes(PlotBin, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.8","0.45","0.1"),breaks=c(0.8,0.45,0.1)) + labs(title="Carbon Fixation in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~.,scales="free_y", space = "free")

ggsave(carb.hm1b,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Carbon/Carbon_KOFxns_MGMs_Bins_BinID_by_Function_Pathway_heatmap.png", width=22, height=20, dpi=600)

carb.hm1b1<-ggplot(clr.carb.all.bin, aes(Genus, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.8","0.45","0.1"),breaks=c(0.8,0.45,0.1)) + labs(title="Carbon Fixation in Salton Seawater MAGs by Genus",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~.,scales="free_y", space = "free")

ggsave(carb.hm1b1,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Carbon/Carbon_KOFxns_MGMs_Bins_Genus_by_Function_Pathway_heatmap.png", width=22, height=20, dpi=600)

carb.hm1d<-ggplot(clr.carb.all.bin, aes(PlotBin, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.8","0.45","0.1"),breaks=c(0.8,0.45,0.1)) + labs(title="Carbon Fixation in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~SampDate, scales="free", space = "free")

ggsave(carb.hm1d,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Carbon/Carbon_KOFxns_MGMs_Bins_BinID_by_Function_SampDate_Pathway_best_heatmap.png", width=24, height=20, dpi=600)

carb.hm1e<-ggplot(clr.carb.all.bin, aes(Genus, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.8","0.45","0.1"),breaks=c(0.8,0.45,0.1)) + labs(title="Carbon Fixation in Salton Seawater MAGs by Genus",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~SampDate, scales="free", space = "free")

ggsave(carb.hm1e,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Carbon/Carbon_KOFxns_MGMs_Bins_Genus_by_Function_SampDate_Pathway_best_heatmap.png", width=24, height=20, dpi=600)

# carb.hm1f<-ggplot(clr.carb.all.bin, aes(Pathway, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.8","0.45","0.1"),breaks=c(0.8,0.45,0.1)) + labs(title="Carbon Fixation in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11,face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Depth_m~SampDate,scales="free", space = "free")
#
# ggsave(carb.hm1f,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Carbon/Carbon_KOFxns_MGMs_heatmap1d.png", width=18, height=18, dpi=600)
#
# carb.hm1g<-ggplot(clr.carb.all.bin, aes(Pathway, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.8","0.45","0.1"),breaks=c(0.8,0.45,0.1)) + labs(title="Carbon Fixation in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11,face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_wrap(.~SampDate)
#
# ggsave(carb.hm1g,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Carbon/Carbon_KOFxns_MGMs_heatmap1d.png", width=18, height=18, dpi=600)

carb.hm1e<-ggplot(clr.carb.all.bin[clr.carb.all.bin$Depth_m==0,], aes(Pathway, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.8","0.45","0.1"),breaks=c(0.8,0.45,0.1)) + labs(title="Carbon Fixation in Salton Seawater MAGs - 0m",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)

ggsave(carb.hm1e,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Carbon/Carbon_KOFxns_Bins_Pathways_MGMs_0m_heatmap.png", width=18, height=18, dpi=600)

carb.hm1f<-ggplot(clr.carb.all.bin[clr.carb.all.bin$Depth_m==5,], aes(Pathway, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.8","0.45","0.1"),breaks=c(0.8,0.45,0.1)) + labs(title="Carbon Fixation in Salton Seawater MAGs - 5m",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)

ggsave(carb.hm1f,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Carbon/Carbon_KOFxns_Bins_Pathways_MGMs_5m_heatmap.png", width=18, height=18, dpi=600)

carb.hm1g<-ggplot(clr.carb.all.bin[clr.carb.all.bin$Depth_m==10,], aes(Pathway, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.8","0.45","0.1"),breaks=c(0.8,0.45,0.1)) + labs(title="Carbon Fixation in Salton Seawater MAGs - 10m",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)

ggsave(carb.hm1g,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Carbon/Carbon_KOFxns_Bins_Pathways_MGMs_10m_heatmap.png", width=18, height=18, dpi=600)


#### Pull out Carbon Metabolic Fxns from Binary Data ####
carb.ko.bin.bi<-bin_fxn.binary[,which(colnames(bin_fxn.binary) %in% carb.fxns.bins$KO_ID)] # merge CLR data w/ N fxns found in contigs from KOFamScan
carb.ko.bin.bi$Bin_ID<-rownames(carb.ko.bin.bi)
carb.ko.bin.bi.melt<-melt(carb.ko.bin.bi, by="Bin_ID")
colnames(carb.ko.bin.bi.melt)[which(names(carb.ko.bin.bi.melt) == "variable")] <- "KO_ID"
colnames(carb.ko.bin.bi.melt)[which(names(carb.ko.bin.bi.melt) == "value")] <- "PresAb"
head(carb.ko.bin.bi.melt) #sanity check

clr.carb.ko.bin.bi<-merge(carb.ko.bin.bi.melt,carb.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
head(clr.carb.ko.bin.bi)
colnames(clr.carb.ko.bin.bi)[which(names(clr.carb.ko.bin.bi) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
clr.cov.sum.carb.ko.bin.bi<-as.data.frame(dcast(clr.carb.ko.bin.bi, Bin_ID~KO_Function.KEGG, value.var="PresAb", fun.aggregate=sum)) ###just dcast, nothing is being added here!
rownames(clr.cov.sum.carb.ko.bin.bi)<-clr.cov.sum.carb.ko.bin.bi$Bin_ID
clr.cov.sum.carb.ko.bin.bi[1:4,]

# sanity check
clr.cov.sum.carb.ko.bin.bi$`accB, bccP; acetyl-CoA carboxylase biotin carboxyl carrier protein`[1:4]
head(clr.cov.sum.carb.ko.bin.bi)

#### Carbon Binary Heat Maps ####
# prep for ggplot2 heatmap
clr.carb.ko.bin.bi[1:4,]
clr.carb.all.bin.bi1<-merge(clr.carb.ko.bin.bi,bin_meta_scaled,by="Bin_ID")
clr.carb.all.bin.bi<-merge(clr.carb.all.bin.bi1,mag_tax,by=c("Bin_ID","PlotBin"))

head(clr.carb.all.bin.bi)
clr.carb.all.bin.bi$PlotBin = factor(clr.carb.all.bin.bi$PlotBin, levels=unique(clr.carb.all.bin.bi$PlotBin[order(clr.carb.all.bin.bi$SampDate,clr.carb.all.bin.bi$Depth_m)]), ordered=TRUE)
clr.carb.all.bin.bi$SampDate<-gsub("\\."," ",clr.carb.all.bin.bi$SampDate)
clr.carb.all.bin.bi$SampDate<-factor(clr.carb.all.bin.bi$SampDate, levels=c("August 2021","December 2021","April 2022"))

unique(clr.carb.all.bin.bi$Pathway)
clr.carb.all.bin.bi<-subset(clr.carb.all.bin.bi, clr.carb.all.bin.bi$Pathway!="Multiple Pathways")
"Multiple Pathways" %in% clr.carb.all.bin.bi$Pathway
clr.carb.all.bin.bi$PathShort<-clr.carb.all.bin.bi$Pathway
clr.carb.all.bin.bi$PathShort[(clr.carb.all.bin.bi$PathShort) == "Reductive Tricarboxylic Acid Cycle"] <- "rTCA"
clr.carb.all.bin.bi$PathShort[(clr.carb.all.bin.bi$PathShort) == "3-Hydroxypropionate Bi-cycle"] <- "3HP"
clr.carb.all.bin.bi$PathShort[(clr.carb.all.bin.bi$PathShort) == "Reductive acetyl-CoA Pathway"] <- "RAcCoa"
clr.carb.all.bin.bi$PathShort[(clr.carb.all.bin.bi$PathShort) == "Calvin Cycle"] <- "CBB"
#clr.carb.all.bin.bi$PathShort[(clr.carb.all.bin.bi$PathShort) == "TCA Cycle"] <- "TCA"

clr.carb.all.bin.bi$Pathway<-factor(clr.carb.all.bin.bi$Pathway,levels=c("3-Hydroxypropionate Bi-cycle","Reductive Tricarboxylic Acid Cycle","Reductive acetyl-CoA Pathway","Calvin Cycle"))
clr.carb.all.bin.bi$PathShort<-factor(clr.carb.all.bin.bi$PathShort,levels=c("3HP","rTCA","RAcCoa","CBB"))

clr.carb.all.bin.bi$KO_Function.KEGG = factor(clr.carb.all.bin.bi$KO_Function.KEGG, levels=unique(clr.carb.all.bin.bi$KO_Function.KEGG[order(clr.carb.all.bin.bi$Pathway)]), ordered=TRUE)

head(clr.carb.all.bin.bi)

# Figures
carb.bi.hm1a<-ggplot(clr.carb.all.bin.bi, aes(PlotBin, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Carbon Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(carb.bi.hm1a,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Carbon/PresenceAbsence/Carbon_KOFxns_MGMs_BinID_by_Function_Binary_heatmap.png", width=18, height=13, dpi=600)

carb.bi.hm1b<-ggplot(clr.carb.all.bin.bi, aes(PlotBin, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Carbon Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~.,scales="free_y", space = "free")

ggsave(carb.bi.hm1b,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Carbon/PresenceAbsence/Carbon_KOFxns_MGMs_BinID_by_Function_Pathway_Binary_heatmap.png", width=20, height=20, dpi=600)

carb.bi.hm1b2<-ggplot(clr.carb.all.bin.bi, aes(Genus, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Carbon Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~.,scales="free_y", space = "free")

ggsave(carb.bi.hm1b2,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Carbon/PresenceAbsence/Carbon_KOFxns_MGMs_Genus_by_Function_Pathway_Binary_heatmap.png", width=20, height=20, dpi=600)

carb.bi.hm1d<-ggplot(clr.carb.all.bin.bi, aes(PlotBin, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Carbon Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate,scales="free_x", space = "free")

ggsave(carb.bi.hm1d,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Carbon/PresenceAbsence/Carbon_KOFxns_MGMs_BinID_by_Function_SampDate_Binary_best_heatmap.png", width=20, height=20, dpi=600)

carb.bi.hm1d2<-ggplot(clr.carb.all.bin.bi, aes(Genus, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Carbon Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate,scales="free_x", space = "free")

ggsave(carb.bi.hm1d2,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Carbon/PresenceAbsence/Carbon_KOFxns_MGMs_Genus_by_Function_SampDate_Binary_best_heatmap.png", width=20, height=20, dpi=600)

carb.bi.hm1e<-ggplot(clr.carb.all.bin.bi, aes(PlotBin, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Carbon Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~SampDate, scales="free", space = "free")

ggsave(carb.bi.hm1e,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Carbon/PresenceAbsence/Carbon_KOFxns_MGMs_BinID_by_Function_SampDate_Pathway_Binary_best_heatmap.png", width=22, height=20, dpi=600)

carb.bi.hm1e2<-ggplot(clr.carb.all.bin.bi, aes(Genus, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Carbon Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~SampDate, scales="free", space = "free")

ggsave(carb.bi.hm1e,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Carbon/PresenceAbsence/Carbon_KOFxns_MGMs_Genus_by_Function_SampDate_Pathway_Binary_best_heatmap.png", width=20, height=20, dpi=600)

# carb.bi.hm1e0<-ggplot(clr.carb.all.bin.bi[clr.carb.all.bin.bi$Depth_m==0,], aes(PathShort, KO_Function.KEGG, fill=factor(PresAb))) +
#   geom_tile(colour="black",size=0.25) +
#   scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Carbon Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate)
#
# ggsave(carb.bi.hm1e0,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Carbon/PresenceAbsence/Carbon_KOFxns_Bins_Pathways_Binary_0m_heatmap.png", width=18, height=18, dpi=600)
#
# carb.bi.hm1e5<-ggplot(clr.carb.all.bin.bi[clr.carb.all.bin.bi$Depth_m==5,], aes(PathShort, KO_Function.KEGG, fill=factor(PresAb))) +
#   geom_tile(colour="black",size=0.25) +
#   scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Carbon Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)
#
# ggsave(carb.bi.hm1e5,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Carbon/PresenceAbsence/Carbon_KOFxns_Bins_Pathways_Binary_5m_heatmap.png", width=18, height=18, dpi=600)
#
# carb.bi.hm1e6<-ggplot(clr.carb.all.bin.bi[clr.carb.all.bin.bi$Depth_m==10,], aes(PathShort, KO_Function.KEGG, fill=factor(PresAb))) +
#   geom_tile(colour="black",size=0.25) +
#   scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Carbon Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)
#
# ggsave(carb.bi.hm1e6,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Carbon/PresenceAbsence/Carbon_KOFxns_Bins_Pathways_Binary_10m_heatmap.png", width=18, height=18, dpi=600)

#### Pull Out Phototrophy Fxns from CLR data ####
## heatmaps of traits of interest

bin.clr.na[1:4,1:4]

# pull out Phototrophy functions from CLR transformed, summed coverages (summed coverage per KO)
photo.ko.na.bin<-bin.clr.na[,which(colnames(bin.clr.na) %in% photo.fxn.bins$KO_ID)] # merge CLR data w/ photo-related fxns found in contigs from KOFamScan
photo.ko.na.bin$Bin_ID<-rownames(photo.ko.na.bin)

photo.ko.na.bin.melt<-melt(photo.ko.na.bin, by="Bin_ID")
colnames(photo.ko.na.bin.melt)[which(names(photo.ko.na.bin.melt) == "variable")] <- "KO_ID"
colnames(photo.ko.na.bin.melt)[which(names(photo.ko.na.bin.melt) == "value")] <- "CLR_SumCovPerKO"
head(photo.ko.na.bin.melt) #sanity check

clr.photo.ko.na.bin<-merge(photo.ko.na.bin.melt,photo.kegg,by.x=c("KO_ID"),by.y=c("KO_ID"))
head(clr.photo.ko.na.bin)
unique(clr.photo.ko.na.bin$KO_Function)

colnames(clr.photo.ko.na.bin)[which(names(clr.photo.ko.na.bin) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
clr.cov.sum.photo.ko.na.bin<-as.data.frame(dcast(clr.photo.ko.na.bin, Bin_ID~KO_Function.KEGG, value.var="CLR_SumCovPerKO", fun.aggregate=sum)) ###
rownames(clr.cov.sum.photo.ko.na.bin)<-clr.cov.sum.photo.ko.na.bin$Bin_ID
clr.cov.sum.photo.ko.na.bin[1:4,1:3]

# ## Phototrophy Heat Maps ####
# see max & mean of summed
max(clr.cov.sum.photo.ko.na.bin[,-1])
mean(as.matrix(clr.cov.sum.photo.ko.na.bin[,-1]))

# # first heat map of sulfur KOs

heatmap(as.matrix(clr.cov.sum.photo.ko.na.bin[,-1]), scale = "none")

# prep for ggplot2 heatmap
clr.photo.ko.na.bin[1:4,]
clr.photo.all.bin1<-merge(clr.photo.ko.na.bin,bin_meta_scaled,by="Bin_ID")
clr.photo.all.bin<-merge(clr.photo.all.bin1,mag_tax,by=c("Bin_ID","PlotBin"))

head(clr.photo.all.bin)
clr.photo.all.bin$PlotBin = factor(clr.photo.all.bin$PlotBin, levels=unique(clr.photo.all.bin$PlotBin[order(clr.photo.all.bin$SampDate,clr.photo.all.bin$Genus)]), ordered=TRUE)
clr.photo.all.bin$SampDate<-gsub("\\."," ",clr.photo.all.bin$SampDate)
clr.photo.all.bin$SampDate<-factor(clr.photo.all.bin$SampDate, levels=c("August 2021","December 2021","April 2022"))

clr.photo.all.bin$PathShort<-clr.photo.all.bin$Pathway
clr.photo.all.bin$PathShort[(clr.photo.all.bin$PathShort) == "Proteorhodopsin"] <- "PR"
clr.photo.all.bin$PathShort[(clr.photo.all.bin$PathShort) == "Sensory Rhodopsin"] <- "SR"

clr.photo.all.bin$Pathway<-factor(clr.photo.all.bin$Pathway,levels=c("Proteorhodopsin","Sensory Rhodopsin"))
clr.photo.all.bin$PathShort<-factor(clr.photo.all.bin$PathShort,levels=c("PR","SR"))

clr.photo.all.bin$MethShort<-clr.photo.all.bin$Method
clr.photo.all.bin$MethShort[(clr.photo.all.bin$MethShort) == "Bacterial Rhodopsin"] <- "Bac Rhod"

#clr.photo.all.bin$Method<-factor(clr.photo.all.bin$Method,levels=c("Bacterial Rhodopsin","Oxygenic Photosynthesis","Anoxygenic Photosynthesis"))
#clr.photo.all.bin$MethShort<-factor(clr.photo.all.bin$MethShort,levels=c("Bac Rhod","Ox PS","AnOx PS"))

unique(clr.photo.all.bin$Phototrophy)
#clr.photo.all.bin$Phototrophy<-factor(clr.photo.all.bin$Phototrophy,levels=c("Hetero","Auto"))

clr.photo.all.bin$KO_Function.KEGG = factor(clr.photo.all.bin$KO_Function.KEGG, levels=unique(clr.photo.all.bin$KO_Function.KEGG[order(clr.photo.all.bin$Phototrophy)]), ordered=TRUE)

head(clr.photo.all.bin)

# For heatmap color gradient
max(clr.photo.all.bin$CLR_SumCovPerKO, na.rm=TRUE)
median(clr.photo.all.bin$CLR_SumCovPerKO, na.rm=TRUE)
min(clr.photo.all.bin$CLR_SumCovPerKO, na.rm=TRUE)

# Figures below
# by SampleID

photo.hm1a<-ggplot(clr.photo.all.bin, aes(PlotBin, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.23","0.22","0.21"),breaks=c(0.23,0.22,0.21)) + labs(title="Phototrophy Metabolism in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(photo.hm1a,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_BinID_by_Function_heatmap.png", width=18, height=13, dpi=600)

photo.hm1a2<-ggplot(clr.photo.all.bin, aes(PlotBin, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.23","0.22","0.21"),breaks=c(0.23,0.22,0.21)) + labs(title="Phototrophy Metabolism in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Phototrophy~.,scales="free_y", space = "free")

ggsave(photo.hm1a2,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_BinID_by_Function_Phototrophy_heatmap.png", width=17, height=15, dpi=600)

photo.hm1a3<-ggplot(clr.photo.all.bin, aes(PlotBin, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.23","0.22","0.21"),breaks=c(0.23,0.22,0.21)) + labs(title="Phototrophy Metabolism in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate,scales="free_x", space = "free")

ggsave(photo.hm1a3,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_BinID_by_Function_SampDate_best_heatmap.png", width=20, height=13, dpi=600)

photo.hm1a4<-ggplot(clr.photo.all.bin, aes(PlotBin, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.23","0.22","0.21"),breaks=c(0.23,0.22,0.21)) + labs(title="Phototrophy Metabolism in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~.,scales="free_y", space = "free")

ggsave(photo.hm1a4,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_BinID_by_Function_Phototrophy_System_heatmap2.png", width=17, height=15, dpi=600)

photo.hm1a5<-ggplot(clr.photo.all.bin, aes(PlotBin, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.23","0.22","0.21"),breaks=c(0.23,0.22,0.21)) + labs(title="Phototrophy Metabolism in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Phototrophy~SampDate, scales="free", space = "free")

ggsave(photo.hm1a5,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_BinID_by_Function_SampDate_Phototrophy_best_heatmap.png", width=20, height=15, dpi=600)

photo.hm1a6<-ggplot(clr.photo.all.bin, aes(PlotBin, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.23","0.22","0.21"),breaks=c(0.23,0.22,0.21)) + labs(title="Phototrophy Metabolism in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~SampDate, scales="free", space = "free")

ggsave(photo.hm1a6,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_BinID_by_Function_SampDate_Phototrophy_System_best_heatmap2.png", width=20, height=15, dpi=600)

# by Genus
photo.hm1b<-ggplot(clr.photo.all.bin, aes(Genus, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.23","0.22","0.21"),breaks=c(0.23,0.22,0.21)) + labs(title="Phototrophy Metabolism in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(photo.hm1b,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_Genus_by_Function_heatmap.png", width=18, height=13, dpi=600)

photo.hm1b2<-ggplot(clr.photo.all.bin, aes(Genus, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.23","0.22","0.21"),breaks=c(0.23,0.22,0.21)) + labs(title="Phototrophy Metabolism in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Phototrophy~.,scales="free_y", space = "free")

ggsave(photo.hm1b2,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_Genus_by_Function_Phototrophy_heatmap.png", width=17, height=15, dpi=600)

photo.hm1b3<-ggplot(clr.photo.all.bin, aes(Genus, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.23","0.22","0.21"),breaks=c(0.23,0.22,0.21)) + labs(title="Phototrophy Metabolism in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate,scales="free_x", space = "free")

ggsave(photo.hm1b3,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_Genus_by_Function_SampDate_best_heatmap.png", width=20, height=13, dpi=600)

photo.hm1b4<-ggplot(clr.photo.all.bin, aes(Genus, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.23","0.22","0.21"),breaks=c(0.23,0.22,0.21)) + labs(title="Phototrophy Metabolism in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~.,scales="free_y", space = "free")

ggsave(photo.hm1b4,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_Genus_by_Function_Phototrophy_System_heatmap2.png", width=17, height=15, dpi=600)

photo.hm1b5<-ggplot(clr.photo.all.bin, aes(Genus, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.23","0.22","0.21"),breaks=c(0.23,0.22,0.21)) + labs(title="Phototrophy Metabolism in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Phototrophy~SampDate, scales="free", space = "free")

ggsave(photo.hm1b5,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_Genus_by_Function_SampDate_Phototrophy_best_heatmap.png", width=20, height=15, dpi=600)

photo.hm1b6<-ggplot(clr.photo.all.bin, aes(Genus, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.23","0.22","0.21"),breaks=c(0.23,0.22,0.21)) + labs(title="Phototrophy Metabolism in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~SampDate, scales="free", space = "free")

ggsave(photo.hm1b6,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_Genus_by_Function_SampDate_Phototrophy_System_best_heatmap2.png", width=20, height=15, dpi=600)

photo.hm1e<-ggplot(clr.photo.all.bin[clr.photo.all.bin$Genus==0,], aes(Phototrophy, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.23","0.22","0.21"),breaks=c(0.23,0.22,0.21)) + labs(title="Phototrophy Metabolism in Salton Seawater MAGs - 0m",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)

ggsave(photo.hm1e,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Phototrophy/Phototrophy_KOFxns_Phototrophys_MGMs_0m_heatmap.png", width=18, height=18, dpi=600)

photo.hm1f<-ggplot(clr.photo.all.bin[clr.photo.all.bin$Genus==5,], aes(Phototrophy, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.23","0.22","0.21"),breaks=c(0.23,0.22,0.21)) + labs(title="Phototrophy Metabolism in Salton Seawater MAGs - 5m",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)

ggsave(photo.hm1f,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Phototrophy/Phototrophy_KOFxns_Phototrophys_MGMs_5m_heatmap.png", width=18, height=18, dpi=600)

photo.hm1g<-ggplot(clr.photo.all.bin[clr.photo.all.bin$Genus==10,], aes(Phototrophy, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.23","0.22","0.21"),breaks=c(0.23,0.22,0.21)) + labs(title="Phototrophy Metabolism in Salton Seawater MAGs - 10m",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)

ggsave(photo.hm1g,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Phototrophy/Phototrophy_KOFxns_Phototrophys_MGMs_10m_heatmap.png", width=18, height=18, dpi=600)

#### Pull out Phototrophy Metabolic Fxns from Binary Data ####
#photo.ko.bin.bi<-bin_fxn.binary[,which(colnames(bin_fxn.binary) %in% photo.fxn.bins$KO_ID)] # merge CLR data w/ photoon-related fxns found in contigs from KOFamScan
# pull out Phototrophy functions from CLR transformed, summed coverages (summed coverage per KO)
photo.ko.na.bin<-bin.bi.clr[,which(colnames(bin.bi.clr) %in% photo.fxn.bins$KO_ID)] # merge CLR data w/ photo-related fxns found in contigs from KOFamScan
photo.ko.na.bin$Bin_ID<-rownames(photo.ko.na.bin)

photo.ko.na.bin.bi.melt<-melt(photo.ko.na.bin, by="Bin_ID")
colnames(photo.ko.na.bin.bi.melt)[which(names(photo.ko.na.bin.bi.melt) == "variable")] <- "KO_ID"
colnames(photo.ko.na.bin.bi.melt)[which(names(photo.ko.na.bin.bi.melt) == "value")] <- "CLR_SumCovPerKO"
head(photo.ko.na.bin.bi.melt) #sanity check

clr.photo.ko.na.bin<-merge(photo.ko.na.bin.bi.melt,photo.kegg,by.x=c("KO_ID"),by.y=c("KO_ID"))
head(clr.photo.ko.na.bin)
unique(clr.photo.ko.na.bin$KO_Function)

colnames(clr.photo.ko.na.bin)[which(names(clr.photo.ko.na.bin) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
clr.cov.sum.photo.ko.na.bin<-as.data.frame(dcast(clr.photo.ko.na.bin, Bin_ID~KO_Function.KEGG, value.var="CLR_SumCovPerKO", fun.aggregate=sum)) ###
rownames(clr.cov.sum.photo.ko.na.bin)<-clr.cov.sum.photo.ko.na.bin$Bin_ID
clr.cov.sum.photo.ko.na.bin[1:4,1:3]

#### Phototrophy Binary Heat Maps ####
# prep for ggplot2 heatmap
clr.photo.ko.bin.bi[1:4,]

clr.photo.all.bin.bi1<-merge(clr.photo.ko.bin.bi,bin_meta_scaled,by="Bin_ID")
clr.photo.all.bin.bi<-merge(clr.photo.all.bin.bi1,mag_tax,by=c("Bin_ID","PlotBin"))

head(clr.photo.all.bin.bi)
clr.photo.all.bin.bi$PlotBin = factor(clr.photo.all.bin.bi$PlotBin, levels=unique(clr.photo.all.bin.bi$PlotBin[order(clr.photo.all.bin.bi$SampDate,clr.photo.all.bin.bi$Genus)]), ordered=TRUE)
clr.photo.all.bin.bi$SampDate<-gsub("\\."," ",clr.photo.all.bin.bi$SampDate)
clr.photo.all.bin.bi$SampDate<-factor(clr.photo.all.bin.bi$SampDate, levels=c("August 2021","December 2021","April 2022"))

clr.photo.all.bin.bi$PathShort<-clr.photo.all.bin.bi$Pathway
clr.photo.all.bin.bi$PathShort[(clr.photo.all.bin.bi$PathShort) == "Proteorhodopsin"] <- "PR"
clr.photo.all.bin.bi$PathShort[(clr.photo.all.bin.bi$PathShort) == "Sensory Rhodopsin"] <- "SR"

clr.photo.all.bin.bi$Pathway<-factor(clr.photo.all.bin.bi$Pathway,levels=c("Proteorhodopsin","Sensory Rhodopsin"))
clr.photo.all.bin.bi$PathShort<-factor(clr.photo.all.bin.bi$PathShort,levels=c("PR","SR"))

clr.photo.all.bin.bi$MethShort<-clr.photo.all.bin.bi$Method
clr.photo.all.bin.bi$MethShort[(clr.photo.all.bin.bi$MethShort) == "Bacterial Rhodopsin"] <- "Bac Rhod"

#clr.photo.all.bin.bi$Method<-factor(clr.photo.all.bin.bi$Method,levels=c("Bacterial Rhodopsin","Oxygenic Photosynthesis","Anoxygenic Photosynthesis"))
#clr.photo.all.bin.bi$MethShort<-factor(clr.photo.all.bin.bi$MethShort,levels=c("Bac Rhod","Ox PS","AnOx PS"))

unique(clr.photo.all.bin.bi$Phototrophy)
#clr.photo.all.bin.bi$Phototrophy<-factor(clr.photo.all.bin.bi$Phototrophy,levels=c("Hetero","Auto"))

clr.photo.all.bin.bi$KO_Function.KEGG = factor(clr.photo.all.bin.bi$KO_Function.KEGG, levels=unique(clr.photo.all.bin.bi$KO_Function.KEGG[order(clr.photo.all.bin.bi$Phototrophy)]), ordered=TRUE)

head(clr.photo.all.bin.bi)

# Figures

# By sample ID
photo.bi.hm1a<-ggplot(clr.photo.all.bin.bi, aes(PlotBin, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Phototrophy Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(photo.bi.hm1a,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Phototrophy/PresenceAbsence/Phototrophy_KOFxns_MGMs_by_Function_Binary_heatmap.png", width=18, height=13, dpi=600)

photo.bi.hm1b<-ggplot(clr.photo.all.bin.bi, aes(PlotBin, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Phototrophy Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~.,scales="free_y", space = "free")

ggsave(photo.bi.hm1b,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Phototrophy/PresenceAbsence/Phototrophy_KOFxns_MGMs_by_Function_System_Binary_heatmap.png", width=17, height=15, dpi=600)

photo.bi.hm1c<-ggplot(clr.photo.all.bin.bi, aes(PlotBin, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Phototrophy Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Phototrophy~.,scales="free_y", space = "free")

ggsave(photo.bi.hm1c,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Phototrophy/PresenceAbsence/Phototrophy_KOFxns_MGMs_by_Function_Phototrophy_Binary_heatmap.png", width=17, height=15, dpi=600)

photo.bi.hm1d<-ggplot(clr.photo.all.bin.bi, aes(PlotBin, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Phototrophy Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate,scales="free_x", space = "free")

ggsave(photo.bi.hm1d,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Phototrophy/PresenceAbsence/Phototrophy_KOFxns_MGMs_by_Function_SampDate_Binary_best_heatmap.png", width=20, height=13, dpi=600)

photo.bi.hm1e<-ggplot(clr.photo.all.bin.bi, aes(PlotBin, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Phototrophy Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~SampDate, scales="free", space = "free")

ggsave(photo.bi.hm1e,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Phototrophy/PresenceAbsence/Phototrophy_KOFxns_MGMs_by_Function_SampDate_System_Binary_best_heatmap.png", width=20, height=15, dpi=600)

photo.bi.hm1f<-ggplot(clr.photo.all.bin.bi, aes(PlotBin, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Phototrophy Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Phototrophy~SampDate, scales="free", space = "free")

ggsave(photo.bi.hm1f,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Phototrophy/PresenceAbsence/Phototrophy_KOFxns_MGMs_by_Function_SampDate_Phototrophy_Binary_best_heatmap.png", width=20, height=15, dpi=600)

# By genus

photo.bi.hm1b2<-ggplot(clr.photo.all.bin.bi, aes(Genus, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Phototrophy Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~.,scales="free_y", space = "free")

ggsave(photo.bi.hm1b2,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Phototrophy/PresenceAbsence/Phototrophy_KOFxns_MGMs_Genus_by_Function_System_Binary_heatmap.png", width=17, height=15, dpi=600)

photo.bi.hm1c2<-ggplot(clr.photo.all.bin.bi, aes(Genus, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Phototrophy Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Phototrophy~.,scales="free_y", space = "free")

ggsave(photo.bi.hm1c2,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Phototrophy/PresenceAbsence/Phototrophy_KOFxns_MGMs_Genus_by_Function_Phototrophy_Binary_heatmap.png", width=17, height=15, dpi=600)

photo.bi.hm1d2<-ggplot(clr.photo.all.bin.bi, aes(Genus, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Phototrophy Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate,scales="free_x", space = "free")

ggsave(photo.bi.hm1d2,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Phototrophy/PresenceAbsence/Phototrophy_KOFxns_MGMs_Genus_by_Function_SampDate_Binary_best_heatmap.png", width=20, height=13, dpi=600)

photo.bi.hm1e2<-ggplot(clr.photo.all.bin.bi, aes(Genus, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Phototrophy Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~SampDate, scales="free", space = "free")

ggsave(photo.bi.hm1e2,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Phototrophy/PresenceAbsence/Phototrophy_KOFxns_MGMs_Genus_by_Function_SampDate_System_Binary_best_heatmap.png", width=20, height=15, dpi=600)

photo.bi.hm1f2<-ggplot(clr.photo.all.bin.bi, aes(Genus, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Phototrophy Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Phototrophy~SampDate, scales="free", space = "free")

ggsave(photo.bi.hm1f2,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Phototrophy/PresenceAbsence/Phototrophy_KOFxns_MGMs_Genus_by_Function_SampDate_Phototrophy_Binary_best_heatmap.png", width=20, height=15, dpi=600)

# photo.bi.hm1e0<-ggplot(clr.photo.all.bin.bi[clr.photo.all.bin.bi$Genus==0,], aes(PathShort, KO_Function.KEGG, fill=factor(PresAb))) +
#   geom_tile(colour="black",size=0.25) +
#   scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Phototrophy Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate)
#
# ggsave(photo.bi.hm1e0,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Phototrophy/PresenceAbsence/Phototrophy_KOFxns_Phototrophys_Binary_0m_heatmap.png", width=18, height=18, dpi=600)
#
# photo.bi.hm1e5<-ggplot(clr.photo.all.bin.bi[clr.photo.all.bin.bi$Genus==5,], aes(PathShort, KO_Function.KEGG, fill=factor(PresAb))) +
#   geom_tile(colour="black",size=0.25) +
#   scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Phototrophy Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)
#
# ggsave(photo.bi.hm1e5,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Phototrophy/PresenceAbsence/Phototrophy_KOFxns_Phototrophys_Binary_5m_heatmap.png", width=18, height=18, dpi=600)
#
# photo.bi.hm1e6<-ggplot(clr.photo.all.bin.bi[clr.photo.all.bin.bi$Genus==10,], aes(PathShort, KO_Function.KEGG, fill=factor(PresAb))) +
#   geom_tile(colour="black",size=0.25) +
#   scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Phototrophy Fixation in Salton Seawater MAGs",fill="Presence/Absence") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)
#
# ggsave(photo.bi.hm1e6,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Phototrophy/PresenceAbsence/Phototrophy_KOFxns_Phototrophys_Binary_10m_heatmap.png", width=18, height=18, dpi=600)


#### Pull Out Aerobic Respiration Fxns from CLR data ####
## heatmaps of traits of interest

bin.clr[1:4,1:4]

# pull out Carbon metabolism functions from CLR transformed, summed coverages (summed coverage per KO)
aero.ko.bin<-bin.clr[,which(colnames(bin.clr) %in% aero.fxn.bins$KO_ID)] # merge CLR data w/ aeroon-related fxns found in contigs from KOFamScan
aero.ko.bin$Bin_ID<-rownames(aero.ko.bin)
aero.ko.bin.melt<-melt(aero.ko.bin, by="Bin_ID")
colnames(aero.ko.bin.melt)[which(names(aero.ko.bin.melt) == "variable")] <- "KO_ID"
colnames(aero.ko.bin.melt)[which(names(aero.ko.bin.melt) == "value")] <- "CLR_SumCovPerKO"
head(aero.ko.bin.melt) #sanity check

clr.aero.ko.bin<-merge(aero.ko.bin.melt,aero.kegg,by.x=c("KO_ID"),by.y=c("KO_ID"))
head(clr.aero.ko.bin)
colnames(clr.aero.ko.bin)[which(names(clr.aero.ko.bin) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
clr.cov.sum.aero.ko.bin<-as.data.frame(dcast(clr.aero.ko.bin, Bin_ID~KO_Function.KEGG, value.var="CLR_SumCovPerKO", fun.aggregate=sum)) ###
rownames(clr.cov.sum.aero.ko.bin)<-clr.cov.sum.aero.ko.bin$Bin_ID
clr.cov.sum.aero.ko.bin[1:4,1:4]

# sanity check
clr.cov.sum.aero.ko.bin$`sdhB, frdB; succinate dehydrogenase iron-sulfur subunit [EC:1.3.5.1]`[1:4]
head(clr.cov.sum.aero.ko.bin)

#### Aerobic Respiration Heat Maps ####
# see max & mean of summed
max(clr.cov.sum.aero.ko.bin[,-1])
mean(as.matrix(clr.cov.sum.aero.ko.bin[,-1]))

# first heat map of sulfur KOs
heatmap(as.matrix(clr.cov.sum.aero.ko.bin[,-1]), scale = "none")

colSums(clr.cov.sum.aero.ko.bin[,-1])

heatmap(as.matrix(clr.cov.sum.aero.ko.bin[,-1]), scale = "none")

# prep for ggplot2 heatmap
clr.aero.ko.bin[1:4,]
clr.aero.all.bin1<-merge(clr.aero.ko.bin,bin_meta_scaled,by="Bin_ID")
clr.aero.all.bin<-merge(clr.aero.all.bin1,mag_tax,by=c("Bin_ID","PlotBin"))

head(clr.aero.all.bin)
clr.aero.all.bin$PlotBin = factor(clr.aero.all.bin$PlotBin, levels=unique(clr.aero.all.bin$PlotBin[order(clr.aero.all.bin$SampDate,clr.aero.all.bin$Depth_m)]), ordered=TRUE)
clr.aero.all.bin$SampDate<-gsub("\\."," ",clr.aero.all.bin$SampDate)
clr.aero.all.bin$SampDate<-factor(clr.aero.all.bin$SampDate, levels=c("August 2021","December 2021","April 2022"))

unique(clr.aero.all.bin$Enzyme)

clr.aero.all.bin$EnzShort<-clr.aero.all.bin$Enzyme
clr.aero.all.bin$EnzShort[(clr.aero.all.bin$EnzShort) == "Cytochrome c oxidase"] <- "Cox"
clr.aero.all.bin$EnzShort[(clr.aero.all.bin$EnzShort) == "F-type ATPase"] <- "F-ATPase"
clr.aero.all.bin$EnzShort[(clr.aero.all.bin$EnzShort) == "NADH:quinone oxidoreductase"] <- "NDH-2"
clr.aero.all.bin$EnzShort[(clr.aero.all.bin$EnzShort) == "Fumarate reductase"] <- "FRD"
clr.aero.all.bin$EnzShort[(clr.aero.all.bin$EnzShort) == "Succinate dehydrogenase"] <- "SDH"

clr.aero.all.bin$Enzyme<-factor(clr.aero.all.bin$Enzyme,levels=c("Cytochrome c oxidase","F-type ATPase","NADH:quinone oxidoreductase","Fumarate reductase","Succinate dehydrogenase"))
clr.aero.all.bin$EnzShort<-factor(clr.aero.all.bin$EnzShort,levels=c("Cox","F-ATPase","NDH-2","FRD","SDH"))

clr.aero.all.bin$KO_Function.KEGG = factor(clr.aero.all.bin$KO_Function.KEGG, levels=unique(clr.aero.all.bin$KO_Function.KEGG[order(clr.aero.all.bin$Enzyme)]), ordered=TRUE)

head(clr.aero.all.bin)

# For heatmap color gradient
max(clr.aero.all.bin$CLR_SumCovPerKO)
max(clr.aero.all.bin$CLR_SumCovPerKO)/2
min(clr.aero.all.bin$CLR_SumCovPerKO)

# Figures

# First by Plot Bin
aero.hm1a<-ggplot(clr.aero.all.bin, aes(PlotBin, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.70","0.35","-0.1"),breaks=c(0.70,0.35,-0.1)) + labs(title="Aerobic Respiration in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.35),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(aero.hm1a,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Aerobic_Respiration/Aerobic_Respiration_KOFxns_MGMs_BinID_by_Function_heatmap.png", width=18, height=16, dpi=600)

aero.hm1b<-ggplot(clr.aero.all.bin, aes(PlotBin, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.70","0.35","-0.1"),breaks=c(0.70,0.35,-0.1)) + labs(title="Aerobic Respiration in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.35),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(EnzShort~.,scales="free_y", space = "free")

ggsave(aero.hm1b,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Aerobic_Respiration/Aerobic_Respiration_KOFxns_MGMs_Bins_BinID_by_Function_Enzyme_heatmap.png", width=22, height=20, dpi=600)

aero.hm1c<-ggplot(clr.aero.all.bin, aes(PlotBin, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.70","0.35","-0.1"),breaks=c(0.70,0.35,-0.1)) + labs(title="Aerobic Respiration in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.35),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(EnzShort~SampDate, scales="free", space = "free")

ggsave(aero.hm1c,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Aerobic_Respiration/Aerobic_Respiration_KOFxns_MGMs_Bins_BinID_by_Function_SampDate_Enzyme_best_heatmap.png", width=24, height=20, dpi=600)

# by Genus
aero.hm1b1<-ggplot(clr.aero.all.bin, aes(Genus, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.70","0.35","-0.1"),breaks=c(0.70,0.35,-0.1)) + labs(title="Aerobic Respiration in Salton Seawater MAGs by Genus",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.35),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(EnzShort~.,scales="free_y", space = "free")

ggsave(aero.hm1b1,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Aerobic_Respiration/Aerobic_Respiration_KOFxns_MGMs_Bins_Genus_by_Function_Enzyme_heatmap.png", width=22, height=20, dpi=600)

aero.hm1c1<-ggplot(clr.aero.all.bin, aes(Genus, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.70","0.35","-0.1"),breaks=c(0.70,0.35,-0.1)) + labs(title="Aerobic Respiration in Salton Seawater MAGs by Genus",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.35),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(EnzShort~SampDate, scales="free", space = "free")

ggsave(aero.hm1c1,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Aerobic_Respiration/Aerobic_Respiration_KOFxns_MGMs_Bins_Genus_by_Function_SampDate_Enzyme_best_heatmap.png", width=24, height=20, dpi=600)

# aero.hm1f<-ggplot(clr.aero.all.bin, aes(Enzyme, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.70","0.35","-0.1"),breaks=c(0.70,0.35,-0.1)) + labs(title="Aerobic Respiration in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.35),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11,face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Depth_m~SampDate,scales="free", space = "free")
#
# ggsave(aero.hm1f,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Aerobic_Respiration/Aerobic_Respiration_KOFxns_MGMs_heatmap1d.png", width=18, height=18, dpi=600)
#
# aero.hm1g<-ggplot(clr.aero.all.bin, aes(Enzyme, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.70","0.35","-0.1"),breaks=c(0.70,0.35,-0.1)) + labs(title="Aerobic Respiration in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.35),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11,face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_wrap(.~SampDate)
#
# ggsave(aero.hm1g,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Aerobic_Respiration/Aerobic_Respiration_KOFxns_MGMs_heatmap1d.png", width=18, height=18, dpi=600)
#
# aero.hm1e<-ggplot(clr.aero.all.bin[clr.aero.all.bin$Depth_m==0,], aes(Enzyme, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.70","0.35","-0.1"),breaks=c(0.70,0.35,-0.1)) + labs(title="Aerobic Respiration in Salton Seawater MAGs - 0m",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.35),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11,face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)
#
# ggsave(aero.hm1e,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Aerobic_Respiration/Aerobic_Respiration_KOFxns_Bins_Enzymes_MGMs_0m_heatmap.png", width=18, height=18, dpi=600)
#
# aero.hm1f<-ggplot(clr.aero.all.bin[clr.aero.all.bin$Depth_m==5,], aes(Enzyme, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.70","0.35","-0.1"),breaks=c(0.70,0.35,-0.1)) + labs(title="Aerobic Respiration in Salton Seawater MAGs - 5m",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.35),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)
#
# ggsave(aero.hm1f,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Aerobic_Respiration/Aerobic_Respiration_KOFxns_Bins_Enzymes_MGMs_5m_heatmap.png", width=18, height=18, dpi=600)
#
# aero.hm1g<-ggplot(clr.aero.all.bin[clr.aero.all.bin$Depth_m==10,], aes(Enzyme, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("0.70","0.35","-0.1"),breaks=c(0.70,0.35,-0.1)) + labs(title="Aerobic Respiration in Salton Seawater MAGs - 10m",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.35),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)
#
# ggsave(aero.hm1g,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Aerobic_Respiration/Aerobic_Respiration_KOFxns_Bins_Enzymes_MGMs_10m_heatmap.png", width=18, height=18, dpi=600)


#### Pull out Aerobic Respiration Fxns from Binary Data ####
aero.ko.bin.bi<-bin_fxn.binary[,which(colnames(bin_fxn.binary) %in% aero.fxn.bins$KO_ID)] # merge CLR data w/ N fxns found in contigs from KOFamScan
aero.ko.bin.bi$Bin_ID<-rownames(aero.ko.bin.bi)
aero.ko.bin.bi.melt<-melt(aero.ko.bin.bi, by="Bin_ID")
colnames(aero.ko.bin.bi.melt)[which(names(aero.ko.bin.bi.melt) == "variable")] <- "KO_ID"
colnames(aero.ko.bin.bi.melt)[which(names(aero.ko.bin.bi.melt) == "value")] <- "PresAb"
head(aero.ko.bin.bi.melt) #sanity check

clr.aero.ko.bin.bi<-merge(aero.ko.bin.bi.melt,aero.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
head(clr.aero.ko.bin.bi)
colnames(clr.aero.ko.bin.bi)[which(names(clr.aero.ko.bin.bi) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
clr.cov.sum.aero.ko.bin.bi<-as.data.frame(dcast(clr.aero.ko.bin.bi, Bin_ID~KO_Function.KEGG, value.var="PresAb", fun.aggregate=sum)) ###just dcast, nothing is being added here!
rownames(clr.cov.sum.aero.ko.bin.bi)<-clr.cov.sum.aero.ko.bin.bi$Bin_ID
clr.cov.sum.aero.ko.bin.bi[1:4,]

# sanity check
clr.cov.sum.aero.ko.bin.bi$`sdhB, frdB; succinate dehydrogenase iron-sulfur subunit [EC:1.3.5.1]`[1:4]
head(clr.cov.sum.aero.ko.bin.bi)

#### Aerobic Respiration Binary Heat Maps ####
# prep for ggplot2 heatmap
clr.aero.ko.bin.bi[1:4,]
clr.aero.all.bin.bi1<-merge(clr.aero.ko.bin.bi,bin_meta_scaled,by="Bin_ID")
clr.aero.all.bin.bi<-merge(clr.aero.all.bin.bi1,mag_tax,by=c("Bin_ID","PlotBin"))

head(clr.aero.all.bin.bi)
clr.aero.all.bin.bi$PlotBin = factor(clr.aero.all.bin.bi$PlotBin, levels=unique(clr.aero.all.bin.bi$PlotBin[order(clr.aero.all.bin.bi$SampDate,clr.aero.all.bin.bi$Depth_m)]), ordered=TRUE)
clr.aero.all.bin.bi$SampDate<-gsub("\\."," ",clr.aero.all.bin.bi$SampDate)
clr.aero.all.bin.bi$SampDate<-factor(clr.aero.all.bin.bi$SampDate, levels=c("August 2021","December 2021","April 2022"))

unique(clr.aero.all.bin.bi$Enzyme)
clr.aero.all.bin.bi$EnzShort<-clr.aero.all.bin.bi$Enzyme
clr.aero.all.bin.bi$EnzShort[(clr.aero.all.bin.bi$EnzShort) == "Cytochrome c oxidase"] <- "Cox"
clr.aero.all.bin.bi$EnzShort[(clr.aero.all.bin.bi$EnzShort) == "F-type ATPase"] <- "F-ATPase"
clr.aero.all.bin.bi$EnzShort[(clr.aero.all.bin.bi$EnzShort) == "NADH:quinone oxidoreductase"] <- "NDH-2"
clr.aero.all.bin.bi$EnzShort[(clr.aero.all.bin.bi$EnzShort) == "Fumarate reductase"] <- "FRD"
clr.aero.all.bin.bi$EnzShort[(clr.aero.all.bin.bi$EnzShort) == "Succinate dehydrogenase"] <- "SDH"

clr.aero.all.bin.bi$Enzyme<-factor(clr.aero.all.bin.bi$Enzyme,levels=c("Cytochrome c oxidase","F-type ATPase","NADH:quinone oxidoreductase","Fumarate reductase","Succinate dehydrogenase"))
clr.aero.all.bin.bi$EnzShort<-factor(clr.aero.all.bin.bi$EnzShort,levels=c("Cox","F-ATPase","NDH-2","FRD","SDH"))

clr.aero.all.bin.bi$KO_Function.KEGG = factor(clr.aero.all.bin.bi$KO_Function.KEGG, levels=unique(clr.aero.all.bin.bi$KO_Function.KEGG[order(clr.aero.all.bin.bi$Enzyme)]), ordered=TRUE)

head(clr.aero.all.bin.bi)

# Figures

# First by Plot Bin
aero.bi.hm1a<-ggplot(clr.aero.all.bin.bi, aes(PlotBin, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Aerobic Respiration in Salton Seawater MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(aero.bi.hm1a,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Aerobic_Respiration/PresenceAbsence/Aerobic_Respiration_KOFxns_MGMs_BinID_by_Function_Binary_heatmap.png", width=18, height=13, dpi=600)

aero.bi.hm1b<-ggplot(clr.aero.all.bin.bi, aes(PlotBin, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Aerobic Respiration in Salton Seawater MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(EnzShort~.,scales="free_y", space = "free")

ggsave(aero.bi.hm1b,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Aerobic_Respiration/PresenceAbsence/Aerobic_Respiration_KOFxns_MGMs_BinID_by_Function_Enzyme_Binary_heatmap.png", width=20, height=20, dpi=600)

aero.bi.hm1c<-ggplot(clr.aero.all.bin.bi, aes(PlotBin, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Aerobic Respiration in Salton Seawater MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate,scales="free_x", space = "free")

ggsave(aero.bi.hm1c,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Aerobic_Respiration/PresenceAbsence/Aerobic_Respiration_KOFxns_MGMs_BinID_by_Function_SampDate_Binary_best_heatmap.png", width=20, height=20, dpi=600)

aero.bi.hm1d<-ggplot(clr.aero.all.bin.bi, aes(PlotBin, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Aerobic Respiration in Salton Seawater MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(EnzShort~SampDate, scales="free", space = "free")

ggsave(aero.bi.hm1d,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Aerobic_Respiration/PresenceAbsence/Aerobic_Respiration_KOFxns_MGMs_BinID_by_Function_SampDate_Enzyme_Binary_best_heatmap.png", width=22, height=20, dpi=600)

# by Genus

aero.bi.hm1b1<-ggplot(clr.aero.all.bin.bi, aes(Genus, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Aerobic Respiration in Salton Seawater MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(EnzShort~.,scales="free_y", space = "free")

ggsave(aero.bi.hm1b1,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Aerobic_Respiration/PresenceAbsence/Aerobic_Respiration_KOFxns_MGMs_Genus_by_Function_Enzyme_Binary_heatmap.png", width=20, height=20, dpi=600)

aero.bi.hm1c1<-ggplot(clr.aero.all.bin.bi, aes(Genus, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Aerobic Respiration in Salton Seawater MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate,scales="free_x", space = "free")

ggsave(aero.bi.hm1c1,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Aerobic_Respiration/PresenceAbsence/Aerobic_Respiration_KOFxns_MGMs_Genus_by_Function_SampDate_Binary_best_heatmap.png", width=20, height=20, dpi=600)

aero.bi.hm1d1<-ggplot(clr.aero.all.bin.bi, aes(Genus, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Aerobic Respiration in Salton Seawater MAGs",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(EnzShort~SampDate, scales="free", space = "free")

ggsave(aero.bi.hm1d1,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Aerobic_Respiration/PresenceAbsence/Aerobic_Respiration_KOFxns_MGMs_Genus_by_Function_SampDate_Enzyme_Binary_best_heatmap.png", width=20, height=20, dpi=600)

# aero.bi.hm1e0<-ggplot(clr.aero.all.bin.bi[clr.aero.all.bin.bi$Depth_m==0,], aes(EnzShort, KO_Function.KEGG, fill=factor(PresAb))) +
#   geom_tile(colour="black",size=0.25) +
#   scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Aerobic Respiration in Salton Seawater MAGs",fill="Presence/Absence") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate)
#
# ggsave(aero.bi.hm1e0,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Aerobic_Respiration/PresenceAbsence/Aerobic_Respiration_KOFxns_Bins_Enzymes_Binary_0m_heatmap.png", width=18, height=18, dpi=600)
#
# aero.bi.hm1e5<-ggplot(clr.aero.all.bin.bi[clr.aero.all.bin.bi$Depth_m==5,], aes(EnzShort, KO_Function.KEGG, fill=factor(PresAb))) +
#   geom_tile(colour="black",size=0.25) +
#   scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Aerobic Respiration in Salton Seawater MAGs",fill="Presence/Absence") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)
#
# ggsave(aero.bi.hm1e5,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Aerobic_Respiration/PresenceAbsence/Aerobic_Respiration_KOFxns_Bins_Enzymes_Binary_5m_heatmap.png", width=18, height=18, dpi=600)
#
# aero.bi.hm1e6<-ggplot(clr.aero.all.bin.bi[clr.aero.all.bin.bi$Depth_m==10,], aes(EnzShort, KO_Function.KEGG, fill=factor(PresAb))) +
#   geom_tile(colour="black",size=0.25) +
#   scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Aerobic Respiration in Salton Seawater MAGs",fill="Presence/Absence") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)
#
# ggsave(aero.bi.hm1e6,filename = "figures/MGM_Figs/BinsOnly/FxnDiv/Aerobic_Respiration/PresenceAbsence/Aerobic_Respiration_KOFxns_Bins_Enzymes_Binary_10m_heatmap.png", width=18, height=18, dpi=600)

#### Homogeneity of Variance (CLR data only)- Composition by Groups ####
## betadisper to look at homogeneity of group dispersions (aka variance) when considering multiple variables
# multivariate analogue to Levene's test of homogeneity of variances
# program finds spatial median or centroid of the group, & compare distances of group to centroid/spatial median via ANOVA

#While PERMANOVA tests differences in group means (analogous to MANOVA),
## a related test called PERMDISP can be used to evaluate homogeneity of group dispersion
#(analogous to Levene's test for equal variances). The vegan function for this test is “betadisper”:
## * need a distance matrix!

bin.clr[1:4,1:4] # sample IDs are rows, genes are columns
bin.ko.cov.sum_table[1:4,1:4] # sanity check

# check rownames of CLR & VST transformed feature count data & metadata
rownames(bin_meta_scaled) %in% rownames(bin.clr) #bin.clr was used to make the distance matrix b.euc_dist

# calculate our Euclidean distance matrix using CLR data
mgm.euc.clr_dist <- dist(bin.clr, method = "euclidean")

## betadisper to look at within group variance
# first by compare dispersions by sampling date
mgm.disper5<-betadisper(mgm.euc.clr_dist, mgm_meta$SampDate)
mgm.disper5

permutest(mgm.disper5, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
#Pairwise comparisons:
#  (Observed p-value below diagonal, permuted p-value above diagonal)
#               August.2021 December.2021 April.2022
# August.2021                     0.35300      0.150
# December.2021     0.34305                    0.891
# April.2022        0.13953       0.91214

anova(mgm.disper5) # p = 0.364 --> accept the Null H, spatial medians ARE NOT significantly difference across sample dates

TukeyHSD(mgm.disper5) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other

#                                diff       lwr       upr     p adj
# December.2021-August.2021 -1.2242178 -4.159177 1.710742 0.4551913
# April.2022-August.2021    -1.3368736 -4.271833 1.598086 0.3996944
# April.2022-December.2021  -0.1126558 -3.047615 2.822304 0.9923919

# Visualize dispersions
png('figures/MGM_Figs/BinsOnly/FxnDiv/SSW_MAG_Bin_pcoa_CLR_SummedCoverage_perKO_betadispersion_sampledate.png',width = 700, height = 600, res=100)
plot(mgm.disper5,main = "Centroids and Dispersion based on Aitchison Distance (CLR Data)", col=colorset1$SampDate_Color)
dev.off()

png('figures/MGM_Figs/BinsOnly/FxnDiv/SSW_MAG_Bin_boxplot_CLR_SummedCoverage_perKO_centroid_distance_sampledate.png',width = 700, height = 600, res=100)
boxplot(mgm.disper5,xlab="Sample Collection Date", main = "Distance to Centroid by Category (CLR Data)", sub="Based on Aitchison Distance", col=colorset1$SampDate_Color)
dev.off()

# What about between sampling depths?
mgm.disper6<-betadisper(mgm.euc.clr_dist, mgm_meta$Depth_m)
mgm.disper6

permutest(mgm.disper6, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons

anova(mgm.disper6) # p = 0.7948 --> accept the Null H, spatial medians are NOT significantly difference across sample dates

TukeyHSD(mgm.disper6) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other
#         diff       lwr      upr     p adj
#5-0  0.4449734 -3.466870 4.356816 0.9357599
#10-0 0.8808174 -3.031026 4.792660 0.7772850
#10-5 0.4358440 -3.475999 4.347687 0.9382591

colfunc <- colorRampPalette(c("red", "blue"))
colfunc(3)

# Visualize dispersions
png('figures/MGM_Figs/BinsOnly/FxnDiv/ssw_MAG_Bin_pcoa_CLR_SummedCoverage_per_KO_betadispersion_depth.png',width = 700, height = 600, res=100)
plot(mgm.disper6,main = "Centroids and Dispersion based on Aitchison Distance (CLR Data)", col=colfunc(3))
dev.off()

png('figures/MGM_Figs/BinsOnly/FxnDiv/ssw_MAG_Bin_boxplot_CLR_SummedCoverage_per_KO_centroid_distance_depth.png',width = 700, height = 600, res=100)
boxplot(mgm.disper6,xlab="Sample Collection Depth", main = "Distance to Centroid by Category (CLR Data)", sub="Based on Aitchison Distance", col=colfunc(3))
dev.off()
## Significant differences in homogeneities can be tested using either parametric or permutational tests,
##and parametric post hoc contrasts can also be investigated:

#### PERMANOVAs to Env Variables Across Groups ####

## The currently preferred analysis for evaluating differences among groups is PERMANOVA.
## This analysis partitions sums of squares using dissimilarities,
##  evaluating differences in the centroids of groups in multivariate space.
##  The vegan functions “adonis” and “adonis2” are used to compute PERMANOVA in R.

help(adonis)

## can specify dataframes for analysis, or we can alternatively specify a dissimilarity matrix:

#Other advantages of using PERMANOVA are that we can test for interactions between predictor variables,
## and we can use both categorical and continuous predictor variables.
## An advantage of adonis2 is that we can test for overall model fit, setting by=NULL, or by individual terms (w/ by="terms")
## w/ distance matrices - The adonis2 tests are identical to anova.cca of dbrda. With Euclidean distances, the tests are also identical to anova.cca of rda.

# First make sure your data frames you're comparing are in the same exact order!!
rownames(bin.clr) %in% rownames(bin_meta_scaled)
bin_meta_scaled=bin_meta_scaled[rownames(bin.clr),] ## reorder metadata to match order of CLR data
perm <- with(bin_meta_scaled, how(nperm = 1000, blocks = SampDate))

pnova1<-adonis2(bin.clr ~ DO_Percent_Local*ORP_mV*Temp_DegC*Dissolved_OrganicMatter_RFU*Depth_m*Sulfate_milliM*Sulfide_microM,data=bin_meta_scaled,method = "euclidean",by="terms",permutations=perm)
pnova1
## none are significant

adonis2(bin.clr ~ DO_Percent_Local*ORP_mV*Temp_DegC*Dissolved_OrganicMatter_RFU*Depth_m*Sulfate_milliM*Sulfide_microM,data=bin_meta_scaled,method = "euclidean",by=NULL,permutations=perm)
#         Df SumOfSqs     R2    F Pr(>F)
#Model    23    34412 0.73114 1.8918 0.4825
#Residual 16    12654 0.26886
#Total    39    47066 1.00000

# remove categorical variables
pnova2<-adonis2(bin.clr ~ DO_Percent_Local*ORP_mV*Temp_DegC*Dissolved_OrganicMatter_RFU*Sulfate_milliM*Sulfide_microM,data=bin_meta_scaled,method = "euclidean",by="terms",permutations=perm)
pnova2
# nothing significant

adonis2(bin.clr ~ DO_Percent_Local*ORP_mV*Temp_DegC*Dissolved_OrganicMatter_RFU*Sulfate_milliM*Sulfide_microM,data=bin_meta_scaled,method = "euclidean",by=NULL,permutations=perm)
#         Df SumOfSqs      R2      F   Pr(>F)
#Model    23    34412 0.73114 1.8918 0.4615
#Residual 16    12654 0.26886
#Total    39    47066 1.00000

pnova3<-adonis2(bin.clr ~ DO_Percent_Local*Temp_DegC*Dissolved_OrganicMatter_RFU*Sulfate_milliM*Sulfide_microM,data=bin_meta_scaled,method = "euclidean",by="terms",permutations=perm)
pnova3
#                                   Df SumOfSqs      R2       F   Pr(>F)
#Sulfide_microM                                               1     1165 0.02474  1.4725 0.006993 **
#Temp_DegC:Dissolved_OrganicMatter_RFU                        1     1349 0.02865  1.7052 0.045954 *
#DO_Percent_Local:Sulfide_microM                              1      944 0.02006  1.1935 0.061938 .

adonis2(bin.clr ~ DO_Percent_Local*Temp_DegC*Dissolved_OrganicMatter_RFU*Sulfate_milliM*Sulfide_microM,data=bin_meta_scaled,method = "euclidean",by=NULL,permutations=perm)
#         Df SumOfSqs      R2      F   Pr(>F)
#Model    23    34412 0.73114 1.8918 0.4775
#Residual 16    12654 0.26886
#Total    39    47066 1.00000

pnova4<-adonis2(bin.clr ~ DO_Percent_Local*Temp_DegC*Dissolved_OrganicMatter_RFU*Sulfide_microM,data=bin_meta_scaled,method = "euclidean",by="terms",permutations=perm)
pnova4
#                                         Df SumOfSqs      R2       F   Pr(>F)
#Sulfide_microM                           1     1122 0.02383  1.5127 0.004995 **
#Temp_DegC:Dissolved_OrganicMatter_RFU    1     1256 0.02669  1.6944 0.052947 .

adonis2(bin.clr ~ DO_Percent_Local*Temp_DegC*Dissolved_OrganicMatter_RFU*Sulfide_microM,data=bin_meta_scaled,method = "euclidean",by=NULL,permutations=perm)
#         Df SumOfSqs      R2      F   Pr(>F)
#Model   15    29270 0.62189 2.6316 0.1339
#Residual 24    17796 0.37811
#Total    39    47066 1.00000

pnova4b<-adonis2(bin.clr ~ Dissolved_OrganicMatter_RFU*Temp_DegC*Sulfide_microM,data=bin_meta_scaled,method = "euclidean",by="terms",permutations=perm)
pnova4b
#                                   Df SumOfSqs      R2       F   Pr(>F)
#Sulfide_microM                                        1     1355 0.02880 1.7221 0.003996 **
#Dissolved_OrganicMatter_RFU:Temp_DegC                 1     3882 0.08249 4.9329 0.055944 .

pnova4c<-adonis2(bin.clr ~ Dissolved_OrganicMatter_RFU*Temp_DegC*Sulfide_microM,data=bin_meta_scaled,method = "euclidean",by="terms",permutations=perm)
pnova4c
#                                   Df SumOfSqs      R2       F   Pr(>F)
#Sulfide_microM                                        1     1355 0.02880 1.7221 0.003996 **
#Dissolved_OrganicMatter_RFU:Temp_DegC                 1     3882 0.08249 4.9329 0.047952 *

pnova5<-adonis2(bin.clr ~ ORP_mV*Dissolved_OrganicMatter_RFU*Temp_DegC*Sulfide_microM,data=bin_meta_scaled,method = "euclidean",by="terms",permutations=perm)
pnova5
#                                               Df SumOfSqs      R2       F   Pr(>F)
#ORP_mV                                         1     4239 0.09006 5.5903 0.03397 *
#Dissolved_OrganicMatter_RFU:Temp_DegC          1     1545 0.03283 2.0378 0.05295 .
#ORP_mV:Dissolved_OrganicMatter_RFU:Temp_DegC   1     1519 0.03227 2.0033 0.06993 .

adonis2(bin.clr ~ ORP_mV*Dissolved_OrganicMatter_RFU*Temp_DegC*Sulfide_microM,data=bin_meta_scaled,method = "euclidean",by=NULL,permutations=perm)
#         Df SumOfSqs      R2      F   Pr(>F)
#Model    15    28868 0.61336 2.5383 0.1748
#Residual 24    18197 0.38664
#Total    39    47066 1.00000

pnova6a<-adonis2(bin.clr ~ ORP_mV*Dissolved_OrganicMatter_RFU*Temp_DegC,data=bin_meta_scaled,method = "euclidean",by="terms",permutations=perm)
pnova6a
#                                             Df SumOfSqs      R2      F   Pr(>F)
#ORP_mV                                        1     4239 0.09006 5.7372 0.003996 **
#Dissolved_OrganicMatter_RFU                   1     5542 0.11776 7.5017 0.044955 *
#Temp_DegC                                     1     5995 0.12736 8.1137 0.093906 .
#ORP_mV:Dissolved_OrganicMatter_RFU            1     1261 0.02679 1.7069 0.267732
#ORP_mV:Temp_DegC                              1     3457 0.07345 4.6791 0.167832
#Dissolved_OrganicMatter_RFU:Temp_DegC         1     1521 0.03231 2.0584 0.040959 *
#ORP_mV:Dissolved_OrganicMatter_RFU:Temp_DegC  1     1409 0.02994 1.9075 0.059940 .

adonis2(bin.clr ~ ORP_mV*Dissolved_OrganicMatter_RFU*Temp_DegC,data=bin_meta_scaled,method = "euclidean",by=NULL,permutations=perm) #significant
#         Df SumOfSqs      R2      F  Pr(>F)
#Model     7    23424 0.49768 4.5292 0.01698 *
#Residual 32    23642 0.50232
#Total    39    47066 1.00000

pnova6b<-adonis2(bin.clr ~ ORP_mV*Dissolved_OrganicMatter_RFU,data=bin_meta_scaled,method = "euclidean",by="terms",permutations=perm)
pnova6b # only ORP is significant
adonis2(bin.clr ~ ORP_mV*Dissolved_OrganicMatter_RFU,data=bin_meta_scaled,method = "euclidean",by=NULL,permutations=perm) # significant
# ^ model explains 23.15% of R^2 aka variation

pnova6c<-adonis2(bin.clr ~ ORP_mV*Temp_DegC,data=bin_meta_scaled,method = "euclidean",by="terms",permutations=perm)
pnova6c  # only ORP is significant
adonis2(bin.clr ~ ORP_mV*Temp_DegC,data=bin_meta_scaled,method = "euclidean",by=NULL,permutations=perm) # insignificant
# ^ model explains 37.8% of R^2 aka variation

pnova6d<-adonis2(bin.clr ~ ORP_mV*Sulfide_microM,data=bin_meta_scaled,method = "euclidean",by="terms",permutations=perm)
pnova6d # only ORP is significant
adonis2(bin.clr ~ ORP_mV*Sulfide_microM,data=bin_meta_scaled,method = "euclidean",by=NULL,permutations=perm) #insignificant
# ^ model explains 16.14% of R^2 aka variation

## BEST MODEL as of 5/11/23: explains 49.77% of variation in composition, p=0.023
adonis2(bin.clr ~ ORP_mV*Temp_DegC*Dissolved_OrganicMatter_RFU,data=bin_meta_scaled,method = "euclidean",by=NULL,permutations=perm)
#         Df SumOfSqs      R2      F  Pr(>F)
#Model    7    23424 0.49768 4.5292 0.02298 *
#Residual 32    23642 0.50232
#Total    39    47066 1.00000

### SELF REMINDER FOR R^2
### Coefficient of Determination, denoted R2 or r2
### is the proportion of the variance in the dependent variable that is predictable from the independent variable(s)

### Pseudo F stat for PERMANOVA
### pseudo F-ratio: It compares the total sum of squared dissimilarities (or ranked dissimilarities) among objects belonging to different groups to that of objects belonging to the same group.
### Larger F-ratios indicate more pronounced group separation, however, the significance of this ratio is usually of more interest than its magnitude.

#### Using Shapiro-Wilk test for Normality ####

# are these functions normally distributed? Need to know for linear models in next step
shapiro.test(clr.sulf.ko.bin$CLR_SumCovPerKO) # what is the p-value?

# p > 0.05 states distribution of data are not significantly different from normal distribution
# p < 0.05 means that data is significantly different from a normal distribution
hist(clr.sulf.ko.bin$CLR_SumCovPerKO, col="blue") # with outliars

# visualize Q-Q plot for species richness
qqnorm(clr.sulf.ko.bin$CLR_SumCovPerKO, pch = 1, frame = FALSE)
qqline(clr.sulf.ko.bin$CLR_SumCovPerKO, col = "red", lwd = 2)

shapiro.test(clr.ars.ko$Bac_Species_Richness) # what is the p-value?
# p-value = 0.02873 w/ outliars
shapiro.test(clr.ars.ko2$Bac_Species_Richness) # what is the p-value? * No outliars
# p-value =  0.01219; no outliars
# p > 0.05 states distribution of data are not significantly different from normal distribution
# p < 0.05 means that data is significantly different from a normal distribution
hist(clr.ars.ko$Bac_Species_Richness, col="blue")

# visualize Q-Q plot for species richness
qqnorm(clr.ars.ko$Bac_Species_Richness, pch = 1, frame = FALSE) # with outliars
qqline(clr.ars.ko$Bac_Species_Richness, col = "red", lwd = 2)

### NOTE: sulf.ko.bin.clr.all has dropped outliers based on Shannon Diversity!

shapiro.test(bin_meta_scaled$DO_Percent_Local) # p-value = 0.3369
hist(bin_meta_scaled$DO_Percent_Local, col="blue")
# visualize Q-Q plot for species richness
qqnorm(bin_meta_scaled$DO_Percent_Local, pch = 1, frame = FALSE) # with outliars
qqline(bin_meta_scaled$DO_Percent_Local, col = "red", lwd = 2)

shapiro.test(bin_meta_scaled$ORP_mV) # p-value = 9.373e-06
hist(bin_meta_scaled$ORP_mV, col="blue")
# visualize Q-Q plot for species richness
qqnorm(bin_meta_scaled$ORP_mV, pch = 1, frame = FALSE) # with outliars
qqline(bin_meta_scaled$ORP_mV, col = "red", lwd = 2)

shapiro.test(bin_meta_scaled$Temp_DegC) # p-value = 5.39e-06
hist(bin_meta_scaled$Temp_DegC, col="blue")
# visualize Q-Q plot for species richness
qqnorm(bin_meta_scaled$Temp_DegC, pch = 1, frame = FALSE) # with outliars
qqline(bin_meta_scaled$Temp_DegC, col = "red", lwd = 2)

shapiro.test(bin_meta_scaled$Dissolved_OrganicMatter_RFU) #  p-value = 0.0003007
hist(bin_meta_scaled$Dissolved_OrganicMatter_RFU, col="blue")
# visualize Q-Q plot for species richness
qqnorm(bin_meta_scaled$Dissolved_OrganicMatter_RFU, pch = 1, frame = FALSE) # with outliars
qqline(bin_meta_scaled$Dissolved_OrganicMatter_RFU, col = "red", lwd = 2)

shapiro.test(bin_meta_scaled$Sulfate_milliM) # p-value = 0.4055
hist(bin_meta_scaled$Sulfate_milliM, col="blue")
# visualize Q-Q plot for species richness
qqnorm(bin_meta_scaled$Sulfate_milliM, pch = 1, frame = FALSE) # with outliars
qqline(bin_meta_scaled$Sulfate_milliM, col = "red", lwd = 2)

shapiro.test(bin_meta_scaled$Sulfide_microM) # p-value =  3.462e-06
hist(bin_meta_scaled$Sulfide_microM, col="blue")
# visualize Q-Q plot for species richness
qqnorm(bin_meta_scaled$Sulfide_microM, pch = 1, frame = FALSE) # with outliars
qqline(bin_meta_scaled$Sulfide_microM, col = "red", lwd = 2)

#### Prep Data for Linear Regressions ####

sulf.ko.bin.clr.all<-merge(clr.sulf.ko.bin,bin_meta_scaled,by=c("SampleID"))
ars.ko.clr.all<-merge(clr.ars.ko,bin_meta_scaled,by=c("SampleID"))
osmo.ko.clr.all<-merge(clr.osmo.ko,bin_meta_scaled,by=c("SampleID"))

#### Linear Regression Comparisons ####
head(sulf.ko.bin.clr.all)

sulf.fxn.glm.fit1<-glm(formula = CLR_SumCovPerKO ~ DO_Percent_Local, data=sulf.ko.bin.clr.all)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(sulf.fxn.glm.fit1)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)       0.0112706  0.0004664  24.166   <2e-16 ***
#DO_Percent_Local -0.0009547  0.0005092  -1.875   0.0687 .

sulf.fxn.glm.fit2<-glm(formula = CLR_SumCovPerKO ~ ORP_mV, data=sulf.ko.bin.clr.all)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(sulf.fxn.glm.fit2)

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)   0.0112343  0.0004731  23.745   <2e-16 ***
#ORP_mV      -0.0001512  0.0004968  -0.304    0.763

sulf.fxn.glm.fit3<-glm(formula = CLR_SumCovPerKO ~ Temp_DegC, family = Gamma, data=sulf.ko.bin.clr.all)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(sulf.fxn.glm.fit3)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept) 0.0113225  0.0004451  25.439   <2e-16 ***
#Temp_DegC   0.0011947  0.0004861   2.458   0.0188 *

sulf.fxn.glm.fit5<-glm(formula = CLR_SumCovPerKO ~ Dissolved_OrganicMatter_RFU, family = Gamma, data=sulf.ko.bin.clr.all)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(sulf.fxn.glm.fit5)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)                 0.0112493  0.0004708  23.896   <2e-16 ***
#Dissolved_OrganicMatter_RFU 0.0004269  0.0004659   0.916    0.365

sulf.fxn.glm.fit6<-glm(formula = CLR_SumCovPerKO ~ Sulfate_milliM, family = Gamma, data=sulf.ko.bin.clr.all)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(sulf.fxn.glm.fit6)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)     0.0112325  0.0004772  23.536   <2e-16 ***
#Sulfate_milliM -0.0002664  0.0004893  -0.545    0.589

sulf.fxn.glm.fit7<-glm(formula = CLR_SumCovPerKO ~ Sulfide_microM, family = Gamma, data=sulf.ko.bin.clr.all)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(sulf.fxn.glm.fit7)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)    0.0112379  0.0004723   23.80   <2e-16 ***
#Sulfide_microM 0.0002944  0.0005160    0.57    0.572

sulf.fxn.glm.fit8<-glm(formula = CLR_SumCovPerKO ~ as.numeric(as.character(Depth_m)), family = Gamma, data=sulf.ko.bin.clr.all)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(sulf.fxn.glm.fit8)

fit1<-aov(CLR_SumCovPerKO ~ as.factor(Depth_m), data=sulf.ko.bin.clr.all)
#pairwise.adonis(sulf.ko.bin.clr.all$CLR_SumCovPerKO, sulf.ko.bin.clr.all$Depth_m, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different

summary(fit1)
#Df           Sum Sq Mean Sq    F value   Pr(>F)
#Depth_m      7   4097   585.3   1.114   0.38
#Residuals   31  16294   525.6
Tuk1<-TukeyHSD(fit1)
Tuk1$Depth_m

#plot(CLR_SumCovPerKO ~ Depth_m, data=sulf.ko.bin.clr.all)
#abline(aov(DustComplexity ~ Elevation, data=sulf.ko.bin.clr.all))

# fit.0<-aov(DustComplexity ~ as.factor(Elevation), data=sulf.ko.bin.clr.all)
# summary(fit.0)
# TukeyHSD(fit.0)
# Levene's test with one independent variable
## Levene's tests whether variances of 2 samples are equal
## we want variances to be the same -- want NON SIGNIFICANCE!
## t test assumes that variances are the same, so Levene's test needs to be non significant
## Fligner's test is a Levene's test for data that are not normally distributed
## more here: https://www.geeksforgeeks.org/fligner-killeen-test-in-r-programming/
fligner.test(CLR_SumCovPerKO ~ Depth_m, data = sulf.ko.bin.clr.all)
# Fligner-Killeen:med chi-squared = 4.091, df = 7, p-value = 0.7692
# Which shows that the data do not deviate significantly from homogeneity.
compare_means(CLR_SumCovPerKO ~ Depth_m, data=sulf.ko.bin.clr.all, method="anova",p.adjust.method = "bonferroni") # won't take as.factor(Elevation) as input

### Export Global Env for Other Scripts ####
#save.image("data/Metagenomes/Analysis/SSW_MAG_Bin_Fxn_BetaDiv.Rdata")
# ^ includes all data combined in object bac.dat.all, ASV table (samples are rows, ASVs are columns), mgm_meta, and an ASV count table (where ASVs are rows, not columns)
# Version Information
sessionInfo()