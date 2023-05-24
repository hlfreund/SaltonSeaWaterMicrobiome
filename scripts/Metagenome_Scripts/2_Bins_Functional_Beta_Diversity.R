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
load("data/Metagenomes/Analysis/SSW_MAG_Bin_Fxn_BetaDiv.Rdata")

head(bin_meta_scaled)
arsen.kegg[1:4,]
bin.ko.cov.sum_table[1:4,1:4]
head(bin.clr.ars)

# fixing some col names in meta_scaled
#colnames(meta_scaled)[which(names(meta_scaled) == "DO_Percent_Local")] <- "DO_Percent_Local"
#colnames(meta_scaled)[which(names(meta_scaled) == "Dissolved_Organic Matter_RFU")] <- "Dissolved_OrganicMatter_RFU"

bin_meta_scaled$SampDate<-gsub("\\."," ",bin_meta_scaled$SampDate) # drop period between month & year in SampDate col
bin_meta_scaled$SampDate<-factor(bin_meta_scaled$SampDate,levels=c("August 2021","December 2021","April 2022"))

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
bin.dend_cols <- as.character(meta_scaled$SampDate_Color[order.dendrogram(bin.euc.clr_dend)])
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

#ggsave(bin.pcoa1,filename = "figures/MGM_Figs/SSW_MAG_Bin_pcoa_CLR_SummedCoverage_Per_KO_sampdate.png", width=12, height=10, dpi=600)

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
#ggsave(bin.pcoa2,filename = "figures/MGM_Figs/SSW_MAG_Bin_pcoa_CLR_SummedCoverage_Per_KO.traits_depth.png", width=12, height=10, dpi=600)

#### Traits of Interest - Heat Maps ####
## heatmaps of traits of interest

# Sulfur fxns

sulf.ko<-bin.clr[,which(colnames(bin.clr) %in% sulfur.fxns.bins$KO_ID)]
sulf.ko$Bin_ID<-rownames(sulf.ko)
sulf.ko.melt<-melt(sulf.ko, by="Bin_ID")
colnames(sulf.ko.melt)[which(names(sulf.ko.melt) == "variable")] <- "KO_ID"
colnames(sulf.ko.melt)[which(names(sulf.ko.melt) == "value")] <- "CLR_SumCovPerKO"
head(sulf.ko.melt) #sanity check

sulf.ko.bintax<-merge(sulf.ko.melt,mag_tax,by=c("Bin_ID"))
clr.sulf.ko<-merge(sulf.ko.bintax,sulfur.fxns.bins,by=c("KO_ID"))
clr.cov.sum.sulf.ko<-as.data.frame(dcast(clr.sulf.ko, Bin_ID~KO_Function, value.var="CLR_SumCovPerKO", fun.aggregate=sum)) ###
rownames(clr.cov.sum.sulf.ko)<-clr.cov.sum.sulf.ko$Bin_ID
clr.cov.sum.sulf.ko[1:4,1:4]

# see max & mean of summed
max(clr.cov.sum.sulf.ko[,-1])
mean(as.matrix(clr.cov.sum.sulf.ko[,-1]))

# first heat map of sulfur KOs
heatmap(as.matrix(clr.cov.sum.sulf.ko[,-1]), scale = "none")

colSums(clr.cov.sum.sulf.ko[,-1])

heatmap(as.matrix(clr.cov.sum.sulf.ko[,-1]), scale = "none")

# prep for ggplot2 heatmap
clr.sulf.ko<-clr.sulf.ko[clr.sulf.ko$Genus!="Unknown",]
clr.sulf.ko[1:4,]
clr.sulf.all<-merge(clr.sulf.ko,bin_meta_scaled,by="Bin_ID")
clr.sulf.all$SampDate = gsub("\\."," ",clr.sulf.all$SampDate)
clr.sulf.all$SampDate = factor(clr.sulf.all$SampDate,levels=c("August 2021","December 2021","April 2022"))
clr.sulf.all$Bin_ID = factor(clr.sulf.all$Bin_ID, levels=unique(clr.sulf.all$Bin_ID[order(clr.sulf.all$SampDate,clr.sulf.all$Depth_m)]), ordered=TRUE)

sulf.hm1<-ggplot(clr.sulf.all, aes(Bin_ID, KO_Function, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Heatmap: Sulfur Functions in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Summed Gene Coverage by KO",fill="CLR Coverage Sums Per KO") +
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),
        axis.text = element_text(size=12),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=12),plot.title = element_text(size=17),
        axis.ticks=element_line(size=0.4),panel.border=element_blank()) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(sulf.hm1,filename = "figures/MGM_Figs/SSW_Bins_SulfMetabolism_heatmap.png", width=15, height=15, dpi=600)

sulf.hm2<-ggplot(clr.sulf.all, aes(Genus, KO_Function, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Heatmap: Sulfur Functions in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Summed Gene Coverage by KO",fill="CLR Coverage Sums Per KO") +
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),
        axis.text = element_text(size=12),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=12),plot.title = element_text(size=17),
        axis.ticks=element_line(size=0.4),panel.border=element_blank()) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(sulf.hm2,filename = "figures/MGM_Figs/SSW_Bins_SulfMetabolism_by_Genus_heatmap.png", width=15, height=15, dpi=600)

sulf.hm3<-ggplot(clr.sulf.all, aes(Genus, KO_Function, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) + facet_grid(SampDate~Depth_m) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",na.value="grey50") + labs(title="Heatmap: Sulfur Functions in Salton Seawater MAGs",subtitle="Using CLR-Transformed, Summed Gene Coverage by KO, Grouped by Bin Assigment",fill="CLR Coverage Sums Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=15),
        axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.border=element_blank(),panel.background = element_rect(fill = "white", colour = NA)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0)) + facet_grid(.~SampDate)

ggsave(sulf.hm3,filename = "figures/MGM_Figs/SSW_Bins_SulfMetabolism_by_Genus_Date_Depth_heatmap.png", width=20, height=15, dpi=600)

# only look at two MAGs w/ highest coverage of Sulfur fxns
clr.sulf.some<-clr.sulf.all[clr.sulf.all$Genus %in% c("Casp-actino5","HIMB30"),]
clr.sulf.some$Bin_ID = factor(clr.sulf.some$Bin_ID, levels=unique(clr.sulf.some$Bin_ID[order(clr.sulf.some$SampDate,clr.sulf.some$Depth_m)]), ordered=TRUE)
NA %in% clr.sulf.some$CLR_SumCovPerKO

some.sulf.hm<-ggplot(clr.sulf.some, aes(Depth_m, KO_Function, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25)  +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",na.value="grey50") + labs(title="MAGs Containing Highest Coverages for Sulfur Metabolic Functions",subtitle="Using CLR-Transformed, Summed Gene Coverage by KO, Grouped by Bin Assigment",fill="CLR Coverage Sums Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=15),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.border=element_blank(),panel.background = element_rect(fill = "white", colour = NA)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0)) + facet_grid(.~SampDate)

ggsave(some.sulf.hm,filename = "figures/MGM_Figs/Bins_by_Genus_HighSulfMetabolism_bySampDate_Depth_heatmap.png", width=20, height=15, dpi=600)

some.sulf.hm1<-ggplot(clr.sulf.some, aes(Bin_ID, KO_Function, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25)  +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="MAGs Containing Highest Coverages for Sulfur Metabolic Functions",subtitle="Using CLR-Transformed, Summed Gene Coverage by KO, Grouped by Bin Assigment",fill="CLR Coverage Sums Per KO") +
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),
        axis.text = element_text(size=12),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=12),plot.title = element_text(size=17),
        axis.ticks=element_line(size=0.4),panel.border=element_blank(),panel.background = element_rect(fill = "white", colour = NA)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0)) + facet_grid(.~Genus)

ggsave(some.sulf.hm1,filename = "figures/MGM_Figs/Bins_by_Genus_HighSulfMetabolism_bySampDate_Depth_heatmap.png", width=20, height=15, dpi=600)

# pull out specific S functions
## first, SOX
clr.bin.Sox<-clr.sulf.all[grepl('Sox', clr.sulf.all$KO_Function),] # pull out just Sox functions
#clr.bin.Sox$SampDate = factor(clr.bin.Sox$SampDate, levels=c("August 2021","December 2021", "April 2022"))
clr.bin.Sox$Bin_ID = factor(clr.bin.Sox$Bin_ID, levels=unique(clr.bin.Sox$Bin_ID[order(clr.bin.Sox$SampDate,clr.bin.Sox$Depth_m)]), ordered=TRUE)

NA %in% clr.bin.Sox$CLR_SumCovPerKO

s.sox.hm<-ggplot(clr.bin.Sox, aes(Genus, KO_Function, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25)  +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",na.value="grey50") + labs(title="SOX Functions in Salton Seawater MAGs by Genus",subtitle="Using CLR-Transformed, Summed Gene Coverage by KO, Grouped by Bin Assigment",fill="CLR Coverage Sums Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=15),
        axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.border=element_blank(),panel.background = element_rect(fill = "white", colour = NA)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0)) + facet_grid(SampDate~Depth_m)

ggsave(s.sox.hm,filename = "figures/MGM_Figs/SSW_S_SOX_Bins_bySampDate_Depth_heatmap.png", width=15, height=10, dpi=600)

# Assimilatory sulfate reduction
clr.bin.as.S.redox<-clr.sulf.all[grepl('1.8.7.1|1.8.1.2|2.7.7.4|1.8.4.10|1.8.4.8|2.7.1.25', clr.sulf.all$KO_Function),] # pull out just assimilatory sulfate reduction functions
clr.bin.as.S.redox$SampDate = factor(clr.bin.as.S.redox$SampDate, levels=c("August 2021","December 2021", "April 2022"))
clr.bin.as.S.redox$Bin_ID = factor(clr.bin.as.S.redox$Bin_ID, levels=unique(clr.bin.as.S.redox$Bin_ID[order(clr.bin.as.S.redox$SampDate,clr.bin.as.S.redox$Depth_m)]), ordered=TRUE)

NA %in% clr.bin.as.S.redox$CLR_SumCovPerKO

s.R.hm1<-ggplot(clr.bin.as.S.redox, aes(Genus, KO_Function, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25)  +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",na.value="grey50") + labs(title="Assimilatory Sulfuate Reduction in Salton Seawater MAGs by Genus",subtitle="Using CLR-Transformed, Summed Gene Coverage by KO, Grouped by Bin Assigment",fill="CLR Coverage Sums Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=15),
        axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.border=element_blank(),panel.background = element_rect(fill = "white", colour = NA)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0)) + facet_grid(SampDate~Depth_m)

ggsave(s.R.hm1,filename = "figures/MGM_Figs/SSW_S_AssSO4_Reduction_Bins_bySampDate_Depth_heatmap.png", width=15, height=10, dpi=600)

# Dissimilatory sulfate reduction and oxidation
clr.bin.dis.S.redox<-clr.sulf.all[grepl('1.8.99.2|1.8.99.5|2.7.7.4', clr.sulf.all$KO_Function),] # pull out just Sox functions
clr.bin.dis.S.redox$SampDate = factor(clr.bin.dis.S.redox$SampDate, levels=c("August 2021","December 2021", "April 2022"))
clr.bin.dis.S.redox$Bin_ID = factor(clr.bin.dis.S.redox$Bin_ID, levels=unique(clr.bin.dis.S.redox$Bin_ID[order(clr.bin.dis.S.redox$SampDate,clr.bin.dis.S.redox$Depth_m)]), ordered=TRUE)

NA %in% clr.bin.dis.S.redox$CLR_SumCovPerKO

s.RO.hm1<-ggplot(clr.bin.dis.S.redox, aes(Genus, KO_Function, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25)  +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",na.value="grey50") + labs(title="Dissimilarity Sulfuate RedOx in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Summed Gene Coverage by KO, Grouped by Bin Assigment",fill="CLR Coverage Sums Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=15),
        axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.border=element_blank(),panel.background = element_rect(fill = "white", colour = NA)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0)) + facet_grid(SampDate~Depth_m)

ggsave(s.RO.hm1,filename = "figures/MGM_Figs/SSW_S_DissSO4_RedOx_Bins_bySampDate_Depth_heatmap.png", width=15, height=10, dpi=600)

## Arsenic Fxns
ars.ko<-bin.clr[,which(colnames(bin.clr) %in% arsenic.fxns$KO_ID)]
ars.ko$Bin_ID<-rownames(ars.ko)
ars.ko.melt<-melt(ars.ko, by="Bin_ID")
colnames(ars.ko.melt)[which(names(ars.ko.melt) == "variable")] <- "KO_ID"
colnames(ars.ko.melt)[which(names(ars.ko.melt) == "value")] <- "CLR_SumCovPerKO"
ars.ko.melt #sanity check

clr.ars.ko<-merge(ars.ko.melt,arsenic.fxns,by=c("KO_ID"))
clr.cov.sum.ars.ko<-as.data.frame(dcast(clr.ars.ko, Bin_ID~KO_Function, value.var="CLR_SumCovPerKO", fun.aggregate=sum)) ###
rownames(clr.cov.sum.ars.ko)<-clr.cov.sum.ars.ko$Bin_ID
clr.cov.sum.ars.ko[1:4,1:4]

# see max & mean of summed
max(clr.cov.sum.ars.ko[,-1])
mean(as.matrix(clr.cov.sum.ars.ko[,-1]))

# first heat map of arsenic KOs
heatmap(as.matrix(clr.cov.sum.ars.ko[,-1]), scale = "none")

colSums(clr.cov.sum.ars.ko[,-1])
clr.cov.sum.ars.ko_2 <- clr.cov.sum.ars.ko[,which(colSums(clr.cov.sum.ars.ko[,-1])>10)]

heatmap(as.matrix(clr.cov.sum.ars.ko_2[,-1]), scale = "none")

# prep for ggplot2 heatmap
clr.ars.ko[1:4,]
clr.ars.all<-merge(clr.ars.ko,meta_scaled,by="Bin_ID")
clr.ars.all$Bin_ID = factor(clr.ars.all$Bin_ID, levels=unique(clr.ars.all$Bin_ID[order(clr.ars.all$SampDate,clr.ars.all$Depth_m)]), ordered=TRUE)

ars.hm1<-ggplot(clr.ars.all, aes(Bin_ID, KO_Function, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Heatmap: Arsenic Functions in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Summed Gene Coverage by KO",fill="CLR Coverage Sums Per KO") +
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),
        axis.text = element_text(size=12),axis.text.x = element_text(vjust=-0.00000000001,angle=45),legend.text = element_text(size=12),plot.title = element_text(size=17),
        axis.ticks=element_line(size=0.4),panel.border=element_blank()) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(ars.hm1,filename = "figures/MGM_Figs/Arsenic_KOFxns_MGMs_heatmap1.png", width=18, height=15, dpi=600)

## Osmoprotectant Fxns
osmo.ko<-bin.clr[,which(colnames(bin.clr) %in% osmo.fxns$KO_ID)]
osmo.ko$Bin_ID<-rownames(osmo.ko)
osmo.ko.melt<-melt(osmo.ko, by="Bin_ID")
colnames(osmo.ko.melt)[which(names(osmo.ko.melt) == "variable")] <- "KO_ID"
colnames(osmo.ko.melt)[which(names(osmo.ko.melt) == "value")] <- "CLR_SumCovPerKO"
osmo.ko.melt #sanity check

clr.osmo.ko<-merge(osmo.ko.melt,osmo.fxns,by=c("KO_ID"))
clr.cov.sum.osmo.ko<-as.data.frame(dcast(clr.osmo.ko, Bin_ID~KO_Function, value.var="CLR_SumCovPerKO", fun.aggregate=sum)) ###
rownames(clr.cov.sum.osmo.ko)<-clr.cov.sum.osmo.ko$Bin_ID
clr.cov.sum.osmo.ko[1:4,1:4]

# see max & mean of summed
max(clr.cov.sum.osmo.ko[,-1])
mean(as.matrix(clr.cov.sum.osmo.ko[,-1]))

# first heat map of arsenic KOs
heatmap(as.matrix(clr.cov.sum.osmo.ko[,-1]), scale = "none")

colSums(clr.cov.sum.osmo.ko[,-1])
#clr.cov.sum.osmo.ko_2 <- clr.cov.sum.osmo.ko[,which(colSums(clr.cov.sum.osmo.ko[,-1])>10)]

heatmap(as.matrix(clr.cov.sum.osmo.ko[,-1]), scale = "none")

# prep for ggplot2 heatmap
clr.osmo.ko[1:4,]
clr.osmo.all<-merge(clr.osmo.ko,meta_scaled,by="Bin_ID")
clr.osmo.all$Bin_ID = factor(clr.osmo.all$Bin_ID, levels=unique(clr.osmo.all$Bin_ID[order(clr.osmo.all$SampDate,clr.osmo.all$Depth_m)]), ordered=TRUE)

osmo.hm1<-ggplot(clr.osmo.all, aes(Bin_ID, KO_Function, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Heatmap: Osmoprotectant Functions in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Summed Gene Coverage by KO",fill="CLR Coverage Sums Per KO") +
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),
        axis.text = element_text(size=12),axis.text.x = element_text(vjust=-0.1,angle=45),legend.text = element_text(size=12),plot.title = element_text(size=17),
        axis.ticks=element_line(size=0.4),panel.border=element_blank()) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(osmo.hm1,filename = "figures/MGM_Figs/OsmoProtectant_KOFxns_MGMs_heatmap1.png", width=18, height=12, dpi=600)

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
rownames(meta_scaled) %in% rownames(bin.clr) #bin.clr was used to make the distance matrix b.euc_dist

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
png('figures/MGM_Figs/SSW_MAG_Bin_pcoa_CLR_SummedCoverage_perKO_betadispersion_sampledate.png',width = 700, height = 600, res=100)
plot(mgm.disper5,main = "Centroids and Dispersion based on Aitchison Distance (CLR Data)", col=colorset1$SampDate_Color)
dev.off()

png('figures/MGM_Figs/SSW_MAG_Bin_boxplot_CLR_SummedCoverage_perKO_centroid_distance_sampledate.png',width = 700, height = 600, res=100)
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
png('figures/MGM_Figs/ssw_MAG_Bin_pcoa_CLR_SummedCoverage_per_KO_betadispersion_depth.png',width = 700, height = 600, res=100)
plot(mgm.disper6,main = "Centroids and Dispersion based on Aitchison Distance (CLR Data)", col=colfunc(3))
dev.off()

png('figures/MGM_Figs/ssw_MAG_Bin_boxplot_CLR_SummedCoverage_per_KO_centroid_distance_depth.png',width = 700, height = 600, res=100)
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
rownames(bin.clr) %in% rownames(meta_scaled)
meta_scaled=meta_scaled[rownames(bin.clr),] ## reorder metadata to match order of CLR data
perm <- with(meta_scaled, how(nperm = 1000, blocks = SampDate))

pnova1<-adonis2(bin.clr ~ DO_Percent_Local*ORP_mV*Temp_DegC*Dissolved_OrganicMatter_RFU*Depth_m*Sulfate_milliM*Sulfide_microM,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
pnova1
## none are significant

adonis2(bin.clr ~ DO_Percent_Local*ORP_mV*Temp_DegC*Dissolved_OrganicMatter_RFU*Depth_m*Sulfate_milliM*Sulfide_microM,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm)
#         Df SumOfSqs     R2    F Pr(>F)
#Model    23    34412 0.73114 1.8918 0.4825
#Residual 16    12654 0.26886
#Total    39    47066 1.00000

# remove categorical variables
pnova2<-adonis2(bin.clr ~ DO_Percent_Local*ORP_mV*Temp_DegC*Dissolved_OrganicMatter_RFU*Sulfate_milliM*Sulfide_microM,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
pnova2
# nothing significant

adonis2(bin.clr ~ DO_Percent_Local*ORP_mV*Temp_DegC*Dissolved_OrganicMatter_RFU*Sulfate_milliM*Sulfide_microM,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm)
#         Df SumOfSqs      R2      F   Pr(>F)
#Model    23    34412 0.73114 1.8918 0.4615
#Residual 16    12654 0.26886
#Total    39    47066 1.00000

pnova3<-adonis2(bin.clr ~ DO_Percent_Local*Temp_DegC*Dissolved_OrganicMatter_RFU*Sulfate_milliM*Sulfide_microM,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
pnova3
#                                   Df SumOfSqs      R2       F   Pr(>F)
#Sulfide_microM                                               1     1165 0.02474  1.4725 0.006993 **
#Temp_DegC:Dissolved_OrganicMatter_RFU                        1     1349 0.02865  1.7052 0.045954 *
#DO_Percent_Local:Sulfide_microM                              1      944 0.02006  1.1935 0.061938 .

adonis2(bin.clr ~ DO_Percent_Local*Temp_DegC*Dissolved_OrganicMatter_RFU*Sulfate_milliM*Sulfide_microM,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm)
#         Df SumOfSqs      R2      F   Pr(>F)
#Model    23    34412 0.73114 1.8918 0.4775
#Residual 16    12654 0.26886
#Total    39    47066 1.00000

pnova4<-adonis2(bin.clr ~ DO_Percent_Local*Temp_DegC*Dissolved_OrganicMatter_RFU*Sulfide_microM,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
pnova4
#                                         Df SumOfSqs      R2       F   Pr(>F)
#Sulfide_microM                           1     1122 0.02383  1.5127 0.004995 **
#Temp_DegC:Dissolved_OrganicMatter_RFU    1     1256 0.02669  1.6944 0.052947 .

adonis2(bin.clr ~ DO_Percent_Local*Temp_DegC*Dissolved_OrganicMatter_RFU*Sulfide_microM,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm)
#         Df SumOfSqs      R2      F   Pr(>F)
#Model   15    29270 0.62189 2.6316 0.1339
#Residual 24    17796 0.37811
#Total    39    47066 1.00000

pnova4b<-adonis2(bin.clr ~ Dissolved_OrganicMatter_RFU*Temp_DegC*Sulfide_microM,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
pnova4b
#                                   Df SumOfSqs      R2       F   Pr(>F)
#Sulfide_microM                                        1     1355 0.02880 1.7221 0.003996 **
#Dissolved_OrganicMatter_RFU:Temp_DegC                 1     3882 0.08249 4.9329 0.055944 .

pnova4c<-adonis2(bin.clr ~ Dissolved_OrganicMatter_RFU*Temp_DegC*Sulfide_microM,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
pnova4c
#                                   Df SumOfSqs      R2       F   Pr(>F)
#Sulfide_microM                                        1     1355 0.02880 1.7221 0.003996 **
#Dissolved_OrganicMatter_RFU:Temp_DegC                 1     3882 0.08249 4.9329 0.047952 *

pnova5<-adonis2(bin.clr ~ ORP_mV*Dissolved_OrganicMatter_RFU*Temp_DegC*Sulfide_microM,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
pnova5
#                                               Df SumOfSqs      R2       F   Pr(>F)
#ORP_mV                                         1     4239 0.09006 5.5903 0.03397 *
#Dissolved_OrganicMatter_RFU:Temp_DegC          1     1545 0.03283 2.0378 0.05295 .
#ORP_mV:Dissolved_OrganicMatter_RFU:Temp_DegC   1     1519 0.03227 2.0033 0.06993 .

adonis2(bin.clr ~ ORP_mV*Dissolved_OrganicMatter_RFU*Temp_DegC*Sulfide_microM,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm)
#         Df SumOfSqs      R2      F   Pr(>F)
#Model    15    28868 0.61336 2.5383 0.1748
#Residual 24    18197 0.38664
#Total    39    47066 1.00000

pnova6a<-adonis2(bin.clr ~ ORP_mV*Dissolved_OrganicMatter_RFU*Temp_DegC,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
pnova6a
#                                             Df SumOfSqs      R2      F   Pr(>F)
#ORP_mV                                        1     4239 0.09006 5.7372 0.003996 **
#Dissolved_OrganicMatter_RFU                   1     5542 0.11776 7.5017 0.044955 *
#Temp_DegC                                     1     5995 0.12736 8.1137 0.093906 .
#ORP_mV:Dissolved_OrganicMatter_RFU            1     1261 0.02679 1.7069 0.267732
#ORP_mV:Temp_DegC                              1     3457 0.07345 4.6791 0.167832
#Dissolved_OrganicMatter_RFU:Temp_DegC         1     1521 0.03231 2.0584 0.040959 *
#ORP_mV:Dissolved_OrganicMatter_RFU:Temp_DegC  1     1409 0.02994 1.9075 0.059940 .

adonis2(bin.clr ~ ORP_mV*Dissolved_OrganicMatter_RFU*Temp_DegC,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm) #significant
#         Df SumOfSqs      R2      F  Pr(>F)
#Model     7    23424 0.49768 4.5292 0.01698 *
#Residual 32    23642 0.50232
#Total    39    47066 1.00000

pnova6b<-adonis2(bin.clr ~ ORP_mV*Dissolved_OrganicMatter_RFU,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
pnova6b # only ORP is significant
adonis2(bin.clr ~ ORP_mV*Dissolved_OrganicMatter_RFU,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm) # significant
# ^ model explains 23.15% of R^2 aka variation

pnova6c<-adonis2(bin.clr ~ ORP_mV*Temp_DegC,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
pnova6c  # only ORP is significant
adonis2(bin.clr ~ ORP_mV*Temp_DegC,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm) # insignificant
# ^ model explains 37.8% of R^2 aka variation

pnova6d<-adonis2(bin.clr ~ ORP_mV*Sulfide_microM,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
pnova6d # only ORP is significant
adonis2(bin.clr ~ ORP_mV*Sulfide_microM,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm) #insignificant
# ^ model explains 16.14% of R^2 aka variation

## BEST MODEL as of 5/11/23: explains 49.77% of variation in composition, p=0.023
adonis2(bin.clr ~ ORP_mV*Temp_DegC*Dissolved_OrganicMatter_RFU,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm)
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
shapiro.test(clr.sulf.ko$CLR_SumCovPerKO) # what is the p-value?

# p > 0.05 states distribution of data are not significantly different from normal distribution
# p < 0.05 means that data is significantly different from a normal distribution
hist(clr.sulf.ko$CLR_SumCovPerKO, col="blue") # with outliars

# visualize Q-Q plot for species richness
qqnorm(clr.sulf.ko$CLR_SumCovPerKO, pch = 1, frame = FALSE)
qqline(clr.sulf.ko$CLR_SumCovPerKO, col = "red", lwd = 2)

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

### NOTE: sulf.ko.clr.all has dropped outliers based on Shannon Diversity!

shapiro.test(meta_scaled$DO_Percent_Local) # p-value = 0.3369
hist(meta_scaled$DO_Percent_Local, col="blue")
# visualize Q-Q plot for species richness
qqnorm(meta_scaled$DO_Percent_Local, pch = 1, frame = FALSE) # with outliars
qqline(meta_scaled$DO_Percent_Local, col = "red", lwd = 2)

shapiro.test(meta_scaled$ORP_mV) # p-value = 9.373e-06
hist(meta_scaled$ORP_mV, col="blue")
# visualize Q-Q plot for species richness
qqnorm(meta_scaled$ORP_mV, pch = 1, frame = FALSE) # with outliars
qqline(meta_scaled$ORP_mV, col = "red", lwd = 2)

shapiro.test(meta_scaled$Temp_DegC) # p-value = 5.39e-06
hist(meta_scaled$Temp_DegC, col="blue")
# visualize Q-Q plot for species richness
qqnorm(meta_scaled$Temp_DegC, pch = 1, frame = FALSE) # with outliars
qqline(meta_scaled$Temp_DegC, col = "red", lwd = 2)

shapiro.test(meta_scaled$Dissolved_OrganicMatter_RFU) #  p-value = 0.0003007
hist(meta_scaled$Dissolved_OrganicMatter_RFU, col="blue")
# visualize Q-Q plot for species richness
qqnorm(meta_scaled$Dissolved_OrganicMatter_RFU, pch = 1, frame = FALSE) # with outliars
qqline(meta_scaled$Dissolved_OrganicMatter_RFU, col = "red", lwd = 2)

shapiro.test(meta_scaled$Sulfate_milliM) # p-value = 0.4055
hist(meta_scaled$Sulfate_milliM, col="blue")
# visualize Q-Q plot for species richness
qqnorm(meta_scaled$Sulfate_milliM, pch = 1, frame = FALSE) # with outliars
qqline(meta_scaled$Sulfate_milliM, col = "red", lwd = 2)

shapiro.test(meta_scaled$Sulfide_microM) # p-value =  3.462e-06
hist(meta_scaled$Sulfide_microM, col="blue")
# visualize Q-Q plot for species richness
qqnorm(meta_scaled$Sulfide_microM, pch = 1, frame = FALSE) # with outliars
qqline(meta_scaled$Sulfide_microM, col = "red", lwd = 2)

#### Prep Data for Linear Regressions ####

sulf.ko.clr.all<-merge(clr.sulf.ko,meta_scaled,by=c("SampleID"))
ars.ko.clr.all<-merge(clr.ars.ko,meta_scaled,by=c("SampleID"))
osmo.ko.clr.all<-merge(clr.osmo.ko,meta_scaled,by=c("SampleID"))

#### Linear Regression Comparisons ####
head(sulf.ko.clr.all)

sulf.fxn.glm.fit1<-glm(formula = CLR_SumCovPerKO ~ DO_Percent_Local, data=sulf.ko.clr.all)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(sulf.fxn.glm.fit1)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)       0.0112706  0.0004664  24.166   <2e-16 ***
#DO_Percent_Local -0.0009547  0.0005092  -1.875   0.0687 .

sulf.fxn.glm.fit2<-glm(formula = CLR_SumCovPerKO ~ ORP_mV, data=sulf.ko.clr.all)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(sulf.fxn.glm.fit2)

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)   0.0112343  0.0004731  23.745   <2e-16 ***
#ORP_mV      -0.0001512  0.0004968  -0.304    0.763

sulf.fxn.glm.fit3<-glm(formula = CLR_SumCovPerKO ~ Temp_DegC, family = Gamma, data=sulf.ko.clr.all)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(sulf.fxn.glm.fit3)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept) 0.0113225  0.0004451  25.439   <2e-16 ***
#Temp_DegC   0.0011947  0.0004861   2.458   0.0188 *

sulf.fxn.glm.fit5<-glm(formula = CLR_SumCovPerKO ~ Dissolved_OrganicMatter_RFU, family = Gamma, data=sulf.ko.clr.all)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(sulf.fxn.glm.fit5)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)                 0.0112493  0.0004708  23.896   <2e-16 ***
#Dissolved_OrganicMatter_RFU 0.0004269  0.0004659   0.916    0.365

sulf.fxn.glm.fit6<-glm(formula = CLR_SumCovPerKO ~ Sulfate_milliM, family = Gamma, data=sulf.ko.clr.all)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(sulf.fxn.glm.fit6)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)     0.0112325  0.0004772  23.536   <2e-16 ***
#Sulfate_milliM -0.0002664  0.0004893  -0.545    0.589

sulf.fxn.glm.fit7<-glm(formula = CLR_SumCovPerKO ~ Sulfide_microM, family = Gamma, data=sulf.ko.clr.all)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(sulf.fxn.glm.fit7)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)    0.0112379  0.0004723   23.80   <2e-16 ***
#Sulfide_microM 0.0002944  0.0005160    0.57    0.572

sulf.fxn.glm.fit8<-glm(formula = CLR_SumCovPerKO ~ as.numeric(as.character(Depth_m)), family = Gamma, data=sulf.ko.clr.all)%>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(sulf.fxn.glm.fit8)

fit1<-aov(CLR_SumCovPerKO ~ as.factor(Depth_m), data=sulf.ko.clr.all)
#pairwise.adonis(sulf.ko.clr.all$CLR_SumCovPerKO, sulf.ko.clr.all$Depth_m, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different

summary(fit1)
#Df           Sum Sq Mean Sq    F value   Pr(>F)
#Depth_m      7   4097   585.3   1.114   0.38
#Residuals   31  16294   525.6
Tuk1<-TukeyHSD(fit1)
Tuk1$Depth_m

#plot(CLR_SumCovPerKO ~ Depth_m, data=sulf.ko.clr.all)
#abline(aov(DustComplexity ~ Elevation, data=sulf.ko.clr.all))

# fit.0<-aov(DustComplexity ~ as.factor(Elevation), data=sulf.ko.clr.all)
# summary(fit.0)
# TukeyHSD(fit.0)
# Levene's test with one independent variable
## Levene's tests whether variances of 2 samples are equal
## we want variances to be the same -- want NON SIGNIFICANCE!
## t test assumes that variances are the same, so Levene's test needs to be non significant
## Fligner's test is a Levene's test for data that are not normally distributed
## more here: https://www.geeksforgeeks.org/fligner-killeen-test-in-r-programming/
fligner.test(CLR_SumCovPerKO ~ Depth_m, data = sulf.ko.clr.all)
# Fligner-Killeen:med chi-squared = 4.091, df = 7, p-value = 0.7692
# Which shows that the data do not deviate significantly from homogeneity.
compare_means(CLR_SumCovPerKO ~ Depth_m, data=sulf.ko.clr.all, method="anova",p.adjust.method = "bonferroni") # won't take as.factor(Elevation) as input

### Export Global Env for Other Scripts ####
#save.image("data/Metagenomes/Analysis/SSW_MAG_Bin_Fxn_BetaDiv.Rdata")
# ^ includes all data combined in object bac.dat.all, ASV table (samples are rows, ASVs are columns), mgm_meta, and an ASV count table (where ASVs are rows, not columns)
# Version Information
sessionInfo()