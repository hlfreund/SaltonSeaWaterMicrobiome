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
load("data/Metagenomes/Analysis/SSW_mgm_analysis.Rdata") # load Rdata to global env
load("data/Metagenomes/Analysis/SSW_MGM_Fxn_BetaDiv.Rdata") # load Rdata to global env

head(meta_scaled)
arsen.fxns[1:4,1:4]
ko.cov.sum_table[1:4,1:4]
head(mgm.clr.ars)

# fixing some col names in meta_scaled
#colnames(meta_scaled)[which(names(meta_scaled) == "Dissolved_Oxygen_Percent_Local")] <- "DO_Percent_Local"
#colnames(meta_scaled)[which(names(meta_scaled) == "Dissolved_Organic Matter_RFU")] <- "Dissolved_OrganicMatter_RFU"

# Before transformations (i.e., VST, CLR, etc) were done, the following was performed
# featureCounts counted reads that mapped to genes
# Reads mapped to genes / gene length for all genes across all samples
# Gene coverage was then added together for each KO ID, since multiple genes were assigned the same KO ID
# Summed coverage per KO was then transformed via median-ratio, vst, and clr

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
legend("topright",legend = c("August 2021","December 2021","April 2022"),cex=.8,col = c( "#ef781c","#03045e","#059c3f"),pch = 15, bty = "n")

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
legend("topright",legend = c("August 2021","December 2021","April 2022"),cex=.8,col = c( "#ef781c","#03045e","#059c3f"),pch = 15, bty = "n")
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
  labs(title="PCoA: Bacteria/Archaea in Salton Seawater",subtitle="Using Median-Ratio Transformed Feature Data",xlab="PC1 [41.14%]", ylab="PC2 [9.04%]",color="Sample Date")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(mgm.pcoa.mr.meta$SampDate_Color[order(mgm.pcoa.mr.meta$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("PC1 [33.04%]") + ylab("PC2 [29.24%]")

ggsave(pcoa3,filename = "figures/MGM_Figs/SSW_MGM_pcoa_MR_sampdate.png", width=12, height=10, dpi=600)

# sample month shape, depth color
pcoa4<-ggplot(mgm.pcoa.mr.meta, aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampleMonth), size=5)+theme_bw()+
  labs(title="PCoA: Metagenome Functions in Salton Seawater",subtitle="Using Median-Ratio Transformed Feature Data",xlab="PC1", ylab="PC2",color="Depth (m)")+
  theme_classic()+ theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),axis.text = element_text(size=12),axis.text.x = element_text(vjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17))+
  scale_color_continuous(low="blue3",high="red",trans = 'reverse') + scale_shape_discrete(labels=c("August 2021","December 2021","April 2022"),name="Sample Date") +
  xlab("PC1 [33.04%]") + ylab("PC2 [29.24%]")

ggsave(pcoa4,filename = "figures/MGM_Figs/SSW_MGM_pcoa_MR.traits_depth.png", width=12, height=10, dpi=600)

## betadisper to look at homogeneity of group dispersions (aka variance) when considering multiple variables
# multivariate analogue to Levene's test of homogeneity of variances
# program finds spatial median or centroid of the group, & compare distances of group to centroid/spatial median via ANOVA

#While PERMANOVA tests differences in group means (analogous to MANOVA),
## a related test called PERMDISP can be used to evaluate homogeneity of group dispersion
#(analogous to Levene's test for equal variances). The vegan function for this test is “betadisper”:
## * need a distance matrix!

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
legend("topright",legend = c("August 2021","December 2021","April 2022"),cex=.8,col = c( "#03045e","#059c3f","#ef781c"),pch = 15, bty = "n")
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
  labs(title="PCoA: Bacteria/Archaea in Salton Seawater",subtitle="Using Variance Stabilization Transformed Feature Data",xlab="PC1 [41.14%]", ylab="PC2 [9.04%]",color="Sample Date")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(mgm.pcoa.vst.meta$SampDate_Color[order(mgm.pcoa.vst.meta$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("PC1 [25.06%]") + ylab("PC2 [21.98%]")

ggsave(pcoa3,filename = "figures/MGM_Figs/SSW_MGM_pcoa_VST_sampdate.png", width=12, height=10, dpi=600)

# sample month shape, depth color
pcoa4<-ggplot(mgm.pcoa.vst.meta, aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampleMonth), size=5)+theme_bw()+
  labs(title="PCoA: Metagenome Functions in Salton Seawater",subtitle="Using Variance Stabilization Transformed Feature Data",xlab="PC1", ylab="PC2",color="Depth (m)")+
  theme_classic()+ theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),axis.text = element_text(size=12),axis.text.x = element_text(vjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17))+
  scale_color_continuous(low="blue3",high="red",trans = 'reverse') + scale_shape_discrete(labels=c("August 2021","December 2021","April 2022"),name="Sample Date") +
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
ko.cov.sum_table[1:4,1:4] # sanity check

# check rownames of CLR & VST transformed feature count data & metadata
rownames(mgm.clr) %in% rownames(meta_scaled)

## PCOA with CLR transformed data first
# calculate our Euclidean distance matrix using CLR data
mgm.euc.clr_dist <- dist(mgm.clr, method = "euclidean")

# creating our hierarcical clustering dendrogram
mgm.euc.clr_clust <- hclust(mgm.euc.clr_dist, method="ward.D2")

# let's make it a little nicer...
mgm.euc.clr_dend <- as.dendrogram(mgm.euc.clr_clust, hang=0.2)
mgm.dend_cols <- as.character(meta_scaled$SampDate_Color[order.dendrogram(mgm.euc.clr_dend)])
labels_colors(mgm.euc.clr_dend) <- mgm.dend_cols

plot(mgm.euc.clr_dend, ylab="CLR Euclidean Distance",cex = 0.5) + title(main = "Bacteria/Archaea Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 3)
legend("topright",legend = c("August 2021","December 2021","April 2022"),cex=.8,col = c("#ef781c","#03045e","#059c3f"),pch = 15, bty = "n")
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
  labs(title="PCoA: Bacteria/Archaea in Salton Seawater",subtitle="Using CLR Transformed, Summed Gene Coverage per KO Function",xlab="PC1 [41.14%]", ylab="PC2 [9.04%]",color="Sample Date")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(mgm.pcoa.clr.meta$SampDate_Color[order(mgm.pcoa.clr.meta$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("PC1 [34.06%]") + ylab("PC2 [22.48%]")

ggsave(pcoa5,filename = "figures/MGM_Figs/SSW_MGM_pcoa_CLR_SummedCoverage_Per_KO_sampdate.png", width=12, height=10, dpi=600)

# sample month shape, depth color
pcoa6<-ggplot(mgm.pcoa.clr.meta, aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampleMonth), size=5)+theme_bw()+
  labs(title="PCoA: Metagenome Functions in Salton Seawater",subtitle="Using CLR Transformed, Summed Gene Coverage per KO Function",xlab="PC1", ylab="PC2",color="Depth (m)")+
  theme_classic()+ theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),axis.text = element_text(size=12),axis.text.x = element_text(vjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17))+
  scale_color_continuous(low="blue3",high="red",trans = 'reverse') + scale_shape_discrete(labels=c("August 2021","December 2021","April 2022"),name="Sample Date") +
  xlab("PC1 [34.06%]") + ylab("PC2 [22.48%]")

ggsave(pcoa6,filename = "figures/MGM_Figs/SSW_MGM_pcoa_CLR_SummedCoverage_Per_KO.traits_depth.png", width=12, height=10, dpi=600)

#### Traits of Interest - Heat Maps ####
## heatmaps of traits of interest

# Sulfur fxns
mgm.clr[1:4,1:4]

sulf.ko<-mgm.clr[,which(colnames(mgm.clr) %in% sulfur.fxns$KO_ID)]
sulf.ko$SampleID<-rownames(sulf.ko)
sulf.ko.melt<-melt(sulf.ko, by="SampleID")
colnames(sulf.ko.melt)[which(names(sulf.ko.melt) == "variable")] <- "KO_ID"
colnames(sulf.ko.melt)[which(names(sulf.ko.melt) == "value")] <- "CLR_SumCovPerKO"
sulf.ko.melt #sanity check

clr.sulf.ko<-merge(sulf.ko.melt,sulfur.fxns,by=c("KO_ID"))
clr.cov.sum.sulf.ko<-as.data.frame(dcast(clr.sulf.ko, SampleID~KO_Function, value.var="CLR_SumCovPerKO", fun.aggregate=sum)) ###
rownames(clr.cov.sum.sulf.ko)<-clr.cov.sum.sulf.ko$SampleID
clr.cov.sum.sulf.ko[1:4,1:4]

# see max & mean of summed
max(clr.cov.sum.sulf.ko[,-1])
mean(as.matrix(clr.cov.sum.sulf.ko[,-1]))

# first heat map of sulfur KOs
heatmap(as.matrix(clr.cov.sum.sulf.ko[,-1]), scale = "none")

colSums(clr.cov.sum.sulf.ko[,-1])
clr.cov.sum.sulf.ko_2 <- clr.cov.sum.sulf.ko[,which(colSums(clr.cov.sum.sulf.ko[,-1])>10)]

heatmap(as.matrix(clr.cov.sum.sulf.ko_2[,-1]), scale = "none")

# prep for ggplot2 heatmap
clr.sulf.ko[1:4,]
clr.sulf.all<-merge(clr.sulf.ko,meta_scaled,by="SampleID")
clr.sulf.all$SampleID = factor(clr.sulf.all$SampleID, levels=unique(clr.sulf.all$SampleID[order(clr.sulf.all$SampDate,clr.sulf.all$Depth_m)]), ordered=TRUE)
clr.sulf.all$SampDate<-gsub("\\."," ",clr.sulf.all$SampDate)
head(clr.sulf.all)

sulf.hm1a<-ggplot(clr.sulf.all, aes(SampleID, KO_Function, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Heatmap: Sulfur Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Summed Gene Coverage by KO",fill="CLR Coverage Sums Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=20),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(sulf.hm1a,filename = "figures/MGM_Figs/Sulfur_KOFxns_MGMs_heatmap1a.png", width=18, height=13, dpi=600)

sulf.hm1b<-ggplot(clr.sulf.all, aes(Depth_m, KO_Function, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Heatmap: Sulfur Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Summed Gene Coverage by KO",fill="CLR Coverage Sums Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=20),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)

ggsave(sulf.hm1b,filename = "figures/MGM_Figs/Sulfur_KOFxns_MGMs_heatmap1b.png", width=18, height=13, dpi=600)

# pull out specific S functions
## first, SOX
clr.Sox<-clr.sulf.all[grepl('Sox', clr.sulf.all$KO_Function),] # pull out just Sox functions
clr.Sox$SampDate = factor(clr.Sox$SampDate, levels=c("August 2021","December 2021", "April 2022"))
clr.Sox$SampleID = factor(clr.Sox$SampleID, levels=unique(clr.Sox$SampleID[order(clr.Sox$SampDate,clr.Sox$Depth_m)]), ordered=TRUE)

NA %in% clr.Sox$CLR_SumCovPerKO

s.sox.hm<-ggplot(clr.Sox, aes(Depth_m, KO_Function, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25)  +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",na.value="grey50") + labs(title="SOX Functions in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Summed Gene Coverage by KO, Grouped by Bin Assigment",fill="CLR Coverage Sums Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=15),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.border=element_blank(),panel.background = element_rect(fill = "white", colour = NA)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0)) + facet_grid(.~SampDate)

ggsave(s.sox.hm,filename = "figures/MGM_Figs/SSW_S_SOX_Contigs_bySampDate_Depth_heatmap.png", width=15, height=10, dpi=600)

# Assimilatory sulfate reduction
clr.as.S.redox<-clr.sulf.all[grepl('1.8.7.1|1.8.1.2|2.7.7.4|1.8.4.10|1.8.4.8|2.7.1.25', clr.sulf.all$KO_Function),] # pull out just assimilatory sulfate reduction functions
clr.as.S.redox$SampDate = factor(clr.as.S.redox$SampDate, levels=c("August 2021","December 2021", "April 2022"))
clr.as.S.redox$SampleID = factor(clr.as.S.redox$SampleID, levels=unique(clr.as.S.redox$SampleID[order(clr.as.S.redox$SampDate,clr.as.S.redox$Depth_m)]), ordered=TRUE)

NA %in% clr.as.S.redox$CLR_SumCovPerKO

s.R.hm1<-ggplot(clr.as.S.redox, aes(Depth_m, KO_Function, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25)  +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",na.value="grey50") + labs(title="Assimilatory Sulfuate Reduction in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Summed Gene Coverage by KO, Grouped by Bin Assigment",fill="CLR Coverage Sums Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=15),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.border=element_blank(),panel.background = element_rect(fill = "white", colour = NA)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0)) + facet_grid(.~SampDate)

ggsave(s.R.hm1,filename = "figures/MGM_Figs/SSW_S_AssSO4_Reduction_Contigs_bySampDate_Depth_heatmap.png", width=15, height=10, dpi=600)

# Dissimilatory sulfate reduction and oxidation
clr.dis.S.redox<-clr.sulf.all[grepl('1.8.99.2|1.8.99.5|2.7.7.4', clr.sulf.all$KO_Function),] # pull out just Sox functions
clr.dis.S.redox$SampDate = factor(clr.dis.S.redox$SampDate, levels=c("August 2021","December 2021", "April 2022"))
clr.dis.S.redox$SampleID = factor(clr.dis.S.redox$SampleID, levels=unique(clr.dis.S.redox$SampleID[order(clr.dis.S.redox$SampDate,clr.dis.S.redox$Depth_m)]), ordered=TRUE)

NA %in% clr.dis.S.redox$CLR_SumCovPerKO

s.RO.hm1<-ggplot(clr.dis.S.redox, aes(Depth_m, KO_Function, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25)  +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",na.value="grey50") + labs(title="Dissimilarity Sulfuate RedOx in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Summed Gene Coverage by KO, Grouped by Bin Assigment",fill="CLR Coverage Sums Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=15),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.border=element_blank(),panel.background = element_rect(fill = "white", colour = NA)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0)) + facet_grid(.~SampDate)

ggsave(s.RO.hm1,filename = "figures/MGM_Figs/SSW_S_DissSO4_RedOx_Contigs_bySampDate_Depth_heatmap.png", width=15, height=10, dpi=600)

## Arsenic Fxns
ars.ko<-mgm.clr[,which(colnames(mgm.clr) %in% arsen.fxns$KO_ID)]
ars.ko$SampleID<-rownames(ars.ko)
ars.ko.melt<-melt(ars.ko, by="SampleID")
colnames(ars.ko.melt)[which(names(ars.ko.melt) == "variable")] <- "KO_ID"
colnames(ars.ko.melt)[which(names(ars.ko.melt) == "value")] <- "CLR_SumCovPerKO"
ars.ko.melt #sanity check

clr.ars.ko<-merge(ars.ko.melt,arsen.fxns,by=c("KO_ID"))
clr.cov.sum.ars.ko<-as.data.frame(dcast(clr.ars.ko, SampleID~KO_Function, value.var="CLR_SumCovPerKO", fun.aggregate=sum)) ###
rownames(clr.cov.sum.ars.ko)<-clr.cov.sum.ars.ko$SampleID
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
clr.ars.all<-merge(clr.ars.ko,meta_scaled,by="SampleID")
clr.ars.all$SampDate<-gsub("\\."," ",clr.ars.all$SampDate)
clr.ars.all$SampDate = factor(clr.ars.all$SampDate, levels=c("August 2021","December 2021", "April 2022"))
clr.ars.all$SampleID = factor(clr.ars.all$SampleID, levels=unique(clr.ars.all$SampleID[order(clr.ars.all$SampDate,clr.ars.all$Depth_m)]), ordered=TRUE)

ars.hm1<-ggplot(clr.ars.all, aes(SampleID, KO_Function, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Heatmap: Arsenic Functions in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Summed Gene Coverage by KO",fill="CLR Coverage Sums Per KO") +
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),
        axis.text = element_text(size=12),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17),
        axis.ticks=element_line(size=0.4),panel.border=element_blank()) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(ars.hm1,filename = "figures/MGM_Figs/Arsenic_KOFxns_MGMs_heatmap1.png", width=18, height=15, dpi=600)

ars.hm2<-ggplot(clr.ars.all, aes(Depth_m, KO_Function, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Heatmap: Arsenic Functions in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Summed Gene Coverage by KO",fill="CLR Coverage Sums Per KO") +
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),
        axis.text = element_text(size=12),axis.text.x = element_text(),legend.text = element_text(size=12),plot.title = element_text(size=17),
        axis.ticks=element_line(size=0.4),panel.border=element_blank()) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)

ggsave(ars.hm2,filename = "figures/MGM_Figs/Arsenic_KOFxns_MGMs_heatmap_better.png", width=18, height=12, dpi=600)

## Selenium Fxns
sel.ko<-mgm.clr[,which(colnames(mgm.clr) %in% selen.fxns$KO_ID)]
sel.ko$SampleID<-rownames(sel.ko)
sel.ko.melt<-melt(sel.ko, by="SampleID")
colnames(sel.ko.melt)[which(names(sel.ko.melt) == "variable")] <- "KO_ID"
colnames(sel.ko.melt)[which(names(sel.ko.melt) == "value")] <- "CLR_SumCovPerKO"
sel.ko.melt #sanity check

clr.sel.ko<-merge(sel.ko.melt,selen.fxns,by=c("KO_ID"))
clr.cov.sum.sel.ko<-as.data.frame(dcast(clr.sel.ko, SampleID~KO_Function, value.var="CLR_SumCovPerKO", fun.aggregate=sum)) ###
rownames(clr.cov.sum.sel.ko)<-clr.cov.sum.sel.ko$SampleID
clr.cov.sum.sel.ko[1:4,1:4]

# see max & mean of summed
max(clr.cov.sum.sel.ko[,-1])
mean(as.matrix(clr.cov.sum.sel.ko[,-1]))

# first heat map of selenic KOs
heatmap(as.matrix(clr.cov.sum.sel.ko[,-1]), scale = "none")

colSums(clr.cov.sum.sel.ko[,-1])
clr.cov.sum.sel.ko_2 <- clr.cov.sum.sel.ko[,which(colSums(clr.cov.sum.sel.ko[,-1])>10)]

heatmap(as.matrix(clr.cov.sum.sel.ko_2[,-1]), scale = "none")

# prep for ggplot2 heatmap
clr.sel.ko[1:4,]
clr.sel.all<-merge(clr.sel.ko,meta_scaled,by="SampleID")
clr.sel.all$SampDate<-gsub("\\."," ",clr.sel.all$SampDate)
clr.sel.all$SampDate = factor(clr.sel.all$SampDate, levels=c("August 2021","December 2021", "April 2022"))
clr.sel.all$SampleID = factor(clr.sel.all$SampleID, levels=unique(clr.sel.all$SampleID[order(clr.sel.all$SampDate,clr.sel.all$Depth_m)]), ordered=TRUE)

sel.hm1<-ggplot(clr.sel.all, aes(SampleID, KO_Function, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Heatmap: Selenium Functions in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Summed Gene Coverage by KO",fill="CLR Coverage Sums Per KO") +
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),
        axis.text = element_text(size=12),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17),
        axis.ticks=element_line(size=0.4),panel.border=element_blank()) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(sel.hm1,filename = "figures/MGM_Figs/Selenium_KOFxns_MGMs_heatmap1.png", width=18, height=15, dpi=600)

sel.hm2<-ggplot(clr.sel.all, aes(Depth_m, KO_Function, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Heatmap: Selenium Functions in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Summed Gene Coverage by KO",fill="CLR Coverage Sums Per KO") +
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),
        axis.text = element_text(size=12),axis.text.x = element_text(),legend.text = element_text(size=12),plot.title = element_text(size=17),
        axis.ticks=element_line(size=0.4),panel.border=element_blank()) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)

ggsave(sel.hm2,filename = "figures/MGM_Figs/Selenium_KOFxns_MGMs_heatmap_better.png", width=18, height=12, dpi=600)

## Osmoprotectant Fxns
osmo.ko<-mgm.clr[,which(colnames(mgm.clr) %in% osmo.fxns$KO_ID)]
osmo.ko$SampleID<-rownames(osmo.ko)
osmo.ko.melt<-melt(osmo.ko, by="SampleID")
colnames(osmo.ko.melt)[which(names(osmo.ko.melt) == "variable")] <- "KO_ID"
colnames(osmo.ko.melt)[which(names(osmo.ko.melt) == "value")] <- "CLR_SumCovPerKO"
osmo.ko.melt #sanity check

clr.osmo.ko<-merge(osmo.ko.melt,osmo.fxns,by=c("KO_ID"))
clr.cov.sum.osmo.ko<-as.data.frame(dcast(clr.osmo.ko, SampleID~KO_Function, value.var="CLR_SumCovPerKO", fun.aggregate=sum)) ###
rownames(clr.cov.sum.osmo.ko)<-clr.cov.sum.osmo.ko$SampleID
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
clr.osmo.all<-merge(clr.osmo.ko,meta_scaled,by="SampleID")
clr.osmo.all$SampDate<-gsub("\\."," ",clr.osmo.all$SampDate)
clr.osmo.all$SampDate = factor(clr.osmo.all$SampDate, levels=c("August 2021","December 2021", "April 2022"))
clr.osmo.all$SampleID = factor(clr.osmo.all$SampleID, levels=unique(clr.osmo.all$SampleID[order(clr.osmo.all$SampDate,clr.osmo.all$Depth_m)]), ordered=TRUE)

osmo.hm1<-ggplot(clr.osmo.all, aes(SampleID, KO_Function, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Heatmap: Osmoprotectant Functions in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Summed Gene Coverage by KO",fill="CLR Coverage Sums Per KO") +
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),
        axis.text = element_text(size=12),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17),
        axis.ticks=element_line(size=0.4),panel.border=element_blank()) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(osmo.hm1,filename = "figures/MGM_Figs/OsmoProtectant_KOFxns_MGMs_heatmap1.png", width=18, height=12, dpi=600)

osmo.hm2<-ggplot(clr.osmo.all, aes(Depth_m, KO_Function, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Heatmap: Osmoprotectant & Tolerance Functions in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Summed Gene Coverage by KO",fill="CLR Coverage Sums Per KO") +
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),
        axis.text = element_text(size=12),axis.text.x = element_text(),legend.text = element_text(size=12),plot.title = element_text(size=17),
        axis.ticks=element_line(size=0.4),panel.border=element_blank()) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)

ggsave(osmo.hm2,filename = "figures/MGM_Figs/OsmoProtectant_KOFxns_MGMs_heatmap_better.png", width=18, height=12, dpi=600)

## Metal Resistance/Tolerance Fxns
met.ko<-mgm.clr[,which(colnames(mgm.clr) %in% metal.fxns$KO_ID)]
met.ko$SampleID<-rownames(met.ko)
met.ko.melt<-melt(met.ko, by="SampleID")
colnames(met.ko.melt)[which(names(met.ko.melt) == "variable")] <- "KO_ID"
colnames(met.ko.melt)[which(names(met.ko.melt) == "value")] <- "CLR_SumCovPerKO"
met.ko.melt #sanity check

clr.met.ko<-merge(met.ko.melt,metal.fxns,by=c("KO_ID"))
clr.cov.sum.met.ko<-as.data.frame(dcast(clr.met.ko, SampleID~KO_Function, value.var="CLR_SumCovPerKO", fun.aggregate=sum)) ###
rownames(clr.cov.sum.met.ko)<-clr.cov.sum.met.ko$SampleID
clr.cov.sum.met.ko[1:4,1:4]

# see max & mean of summed
max(clr.cov.sum.met.ko[,-1])
mean(as.matrix(clr.cov.sum.met.ko[,-1]))

# first heat map of metalic KOs
heatmap(as.matrix(clr.cov.sum.met.ko[,-1]), scale = "none")

colSums(clr.cov.sum.met.ko[,-1])
clr.cov.sum.met.ko_2 <- clr.cov.sum.met.ko[,which(colSums(clr.cov.sum.met.ko[,-1])>10)]

#heatmap(as.matrix(clr.cov.sum.met.ko_2[,-1]), scale = "none")

# prep for ggplot2 heatmap
clr.met.ko[1:4,]
clr.met.all<-merge(clr.met.ko,meta_scaled,by="SampleID")
clr.met.all$SampDate<-gsub("\\."," ",clr.met.all$SampDate)
clr.met.all$SampDate = factor(clr.met.all$SampDate, levels=c("August 2021","December 2021", "April 2022"))
clr.met.all$SampleID = factor(clr.met.all$SampleID, levels=unique(clr.met.all$SampleID[order(clr.met.all$SampDate,clr.met.all$Depth_m)]), ordered=TRUE)

met.hm1<-ggplot(clr.met.all, aes(SampleID, KO_Function, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Heatmap: Metal Resistance & Tolerance Functions in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Summed Gene Coverage by KO",fill="CLR Coverage Sums Per KO") +
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),
        axis.text = element_text(size=12),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17),
        axis.ticks=element_line(size=0.4),panel.border=element_blank()) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(met.hm1,filename = "figures/MGM_Figs/Metal_ResistToler_KOFxns_MGMs_heatmap1.png", width=18, height=15, dpi=600)

met.hm2<-ggplot(clr.met.all, aes(Depth_m, KO_Function, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Heatmap: Metal Resistance & Tolerance Functions in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Summed Gene Coverage by KO",fill="CLR Coverage Sums Per KO") +
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),
        axis.text = element_text(size=12),axis.text.x = element_text(),legend.text = element_text(size=12),plot.title = element_text(size=17),
        axis.ticks=element_line(size=0.4),panel.border=element_blank()) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)

ggsave(met.hm2,filename = "figures/MGM_Figs/Metal_ResistToler_KOFxns_MGMs_heatmap_better.png", width=18, height=15, dpi=600)

#### Sulfur Metabolism PCoA ####
## PCOA with CLR transformed data first
# calculate our Euclidean distance matrix using CLR data
sulf.ko<-mgm.clr[,which(colnames(mgm.clr) %in% sulf.kegg$KO_ID)]
sulf.ko[1:4,1:4]

sulf.euc.clr_dist <- dist(sulf.ko, method = "euclidean")

# creating our hierarcical clustering dendrogram
sulf.euc.clr_clust <- hclust(sulf.euc.clr_dist, method="ward.D2")

# let's make it a little nicer...
sulf.euc.clr_dend <- as.dendrogram(sulf.euc.clr_clust, hang=0.2)
sulf.dend_cols <- as.character(meta_scaled$SampDate_Color[order.dendrogram(sulf.euc.clr_dend)])
labels_colors(sulf.euc.clr_dend) <- sulf.dend_cols

plot(sulf.euc.clr_dend, ylab="CLR Euclidean Distance",cex = 0.5) + title(main = "Bacteria/Archaea Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 3)
legend("topright",legend = c("August 2021","December 2021","April 2022"),cex=.8,col = c("#ef781c","#03045e","#059c3f"),pch = 15, bty = "n")
# Control is dark blue ("#218380"), #Alternaria is light blue ("#73d2de")
dev.off()

# let's use our Euclidean distance matrix from before
sulf.pcoa.clr <- pcoa(sulf.euc.clr_dist) # pcoa of euclidean distance matrix = PCA of euclidean distance matrix
##save.image("data/ssw_clr.euc.dist1_3.7.23.Rdata")

# The proportion of variances explained is in its element values$Relative_eig
sulf.pcoa.clr$values

# extract principal coordinates
sulf.pcoa.clr.vectors<-data.frame(sulf.pcoa.clr$vectors)
sulf.pcoa.clr.vectors$SampleID<-rownames(sulf.pcoa.clr$vectors)

# merge pcoa coordinates w/ metadata
sulf.pcoa.clr.meta<-merge(sulf.pcoa.clr.vectors, mgm_meta, by.x="SampleID", by.y="SampleID")
sulf.pcoa.clr.meta$SampleMonth
sulf.pcoa.clr.meta$SampDate

head(sulf.pcoa.clr.meta)

sulf.pcoa.clr$values # pull out Relative (Relative_eig) variation % to add to axes labels

# create PCoA ggplot fig
pcoa.s1<-ggplot(sulf.pcoa.clr.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(SampDate)), size=4)+theme_bw()+
  labs(title="PCoA: Sulfur Metabolism in Salton Seawater",subtitle="Using CLR Transformed, Summed Gene Coverage per KO Function",xlab="PC1 [41.14%]", ylab="PC2 [9.04%]",color="Sample Date")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(sulf.pcoa.clr.meta$SampDate_Color[order(sulf.pcoa.clr.meta$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("PC1 [58.37%]") + ylab("PC2 [15.42%]")

ggsave(pcoa.s1,filename = "figures/MGM_Figs/SSW_SulfurOnly_pcoa_CLR_SummedCoverage_Per_KO_sampdate.png", width=12, height=10, dpi=600)

# sample month shape, depth color
pcoa.s2<-ggplot(sulf.pcoa.clr.meta, aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampleMonth), size=5)+theme_bw()+
  labs(title="PCoA: Sulfur Metabolsim in Salton Seawater",subtitle="Using CLR Transformed, Summed Gene Coverage per KO Function",xlab="PC1", ylab="PC2",color="Depth (m)")+
  theme_classic()+ theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),axis.text = element_text(size=12),axis.text.x = element_text(vjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17))+
  scale_color_continuous(low="blue3",high="red",trans = 'reverse') + scale_shape_discrete(labels=c("August 2021","December 2021","April 2022"),name="Sample Date") +
  xlab("PC1 [66.54%]") + ylab("PC2 [16.69%]")

ggsave(pcoa.s2,filename = "figures/MGM_Figs/SSW_SulfurOnly_pcoa_CLR_SummedCoverage_Per_KO.traits_depth.png", width=12, height=10, dpi=600)

#### Arsenic Metabolism PCoA ####
## PCOA with CLR transformed data first
# calculate our Euclidean distance matrix using CLR data
ars.ko<-mgm.clr[,which(colnames(mgm.clr) %in% arsen.kegg$KO_ID)]
ars.ko[1:4,1:4]

ars.euc.clr_dist <- dist(ars.ko, method = "euclidean")

# creating our hierarcical clustering dendrogram
ars.euc.clr_clust <- hclust(ars.euc.clr_dist, method="ward.D2")

# let's make it a little nicer...
ars.euc.clr_dend <- as.dendrogram(ars.euc.clr_clust, hang=0.2)
ars.dend_cols <- as.character(meta_scaled$SampDate_Color[order.dendrogram(ars.euc.clr_dend)])
labels_colors(ars.euc.clr_dend) <- ars.dend_cols

plot(ars.euc.clr_dend, ylab="CLR Euclidean Distance",cex = 0.5) + title(main = "Bacteria/Archaea Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 3)
legend("topright",legend = c("August 2021","December 2021","April 2022"),cex=.8,col = c("#ef781c","#03045e","#059c3f"),pch = 15, bty = "n")
# Control is dark blue ("#218380"), #Alternaria is light blue ("#73d2de")
dev.off()

# let's use our Euclidean distance matrix from before
ars.pcoa.clr <- pcoa(ars.euc.clr_dist) # pcoa of euclidean distance matrix = PCA of euclidean distance matrix
##save.image("data/ssw_clr.euc.dist1_3.7.23.Rdata")

# The proportion of variances explained is in its element values$Relative_eig
ars.pcoa.clr$values

# extract principal coordinates
ars.pcoa.clr.vectors<-data.frame(ars.pcoa.clr$vectors)
ars.pcoa.clr.vectors$SampleID<-rownames(ars.pcoa.clr$vectors)

# merge pcoa coordinates w/ metadata
ars.pcoa.clr.meta<-merge(ars.pcoa.clr.vectors, mgm_meta, by.x="SampleID", by.y="SampleID")
ars.pcoa.clr.meta$SampleMonth
ars.pcoa.clr.meta$SampDate

head(ars.pcoa.clr.meta)

ars.pcoa.clr$values # pull out Relative (Relative_eig) variation % to add to axes labels

# create PCoA ggplot fig
pcoa.a1<-ggplot(ars.pcoa.clr.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(SampDate)), size=4)+theme_bw()+
  labs(title="PCoA: Arsenic Metabolism in Salton Seawater",subtitle="Using CLR Transformed, Summed Gene Coverage per KO Function",xlab="PC1 [41.14%]", ylab="PC2 [9.04%]",color="Sample Date")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(ars.pcoa.clr.meta$SampDate_Color[order(ars.pcoa.clr.meta$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("PC1 [53.60%]") + ylab("PC2 [30.40%]")

ggsave(pcoa.a1,filename = "figures/MGM_Figs/SSW_ArsenicOnly_pcoa_CLR_SummedCoverage_Per_KO_sampdate.png", width=12, height=10, dpi=600)

# sample month shape, depth color
pcoa.a2<-ggplot(ars.pcoa.clr.meta, aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampleMonth), size=5)+theme_bw()+
  labs(title="PCoA: Arsenic Metabolsim in Salton Seawater",subtitle="Using CLR Transformed, Summed Gene Coverage per KO Function",xlab="PC1", ylab="PC2",color="Depth (m)")+
  theme_classic()+ theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),axis.text = element_text(size=12),axis.text.x = element_text(vjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17))+
  scale_color_continuous(low="blue3",high="red",trans = 'reverse') + scale_shape_discrete(labels=c("August 2021","December 2021","April 2022"),name="Sample Date") +
  xlab("PC1 [53.60%]") + ylab("PC2 [30.40%]")

ggsave(pcoa.a2,filename = "figures/MGM_Figs/SSW_ArsenicOnly_pcoa_CLR_SummedCoverage_Per_KO.traits_depth.png", width=12, height=10, dpi=600)


#### Metal Resistance/Tolerance PCoA ####
## PCOA with CLR transformed data first
# calculate our Euclidean distance matrix using CLR data
met.ko<-mgm.clr[,which(colnames(mgm.clr) %in% metal.fxns$KO_ID)]
met.ko[1:4,1:4]

met.euc.clr_dist <- dist(met.ko, method = "euclidean")

# creating our hierarcical clustering dendrogram
met.euc.clr_clust <- hclust(met.euc.clr_dist, method="ward.D2")

# let's make it a little nicer...
met.euc.clr_dend <- as.dendrogram(met.euc.clr_clust, hang=0.2)
met.dend_cols <- as.character(meta_scaled$SampDate_Color[order.dendrogram(met.euc.clr_dend)])
labels_colors(met.euc.clr_dend) <- met.dend_cols

plot(met.euc.clr_dend, ylab="CLR Euclidean Distance",cex = 0.5) + title(main = "Bacteria/Archaea Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 3)
legend("topright",legend = c("August 2021","December 2021","April 2022"),cex=.8,col = c("#ef781c","#03045e","#059c3f"),pch = 15, bty = "n")
# Control is dark blue ("#218380"), #Alternaria is light blue ("#73d2de")
dev.off()

# let's use our Euclidean distance matrix from before
met.pcoa.clr <- pcoa(met.euc.clr_dist) # pcoa of euclidean distance matrix = PCA of euclidean distance matrix
##save.image("data/ssw_clr.euc.dist1_3.7.23.Rdata")

# The proportion of variances explained is in its element values$Relative_eig
met.pcoa.clr$values

# extract principal coordinates
met.pcoa.clr.vectors<-data.frame(met.pcoa.clr$vectors)
met.pcoa.clr.vectors$SampleID<-rownames(met.pcoa.clr$vectors)

# merge pcoa coordinates w/ metadata
met.pcoa.clr.meta<-merge(met.pcoa.clr.vectors, mgm_meta, by.x="SampleID", by.y="SampleID")
met.pcoa.clr.meta$SampleMonth
met.pcoa.clr.meta$SampDate

head(met.pcoa.clr.meta)

met.pcoa.clr$values # pull out Relative (Relative_eig) variation % to add to axes labels

# create PCoA ggplot fig
pcoa.m1<-ggplot(met.pcoa.clr.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(SampDate)), size=4)+theme_bw()+
  labs(title="PCoA: Metal Resistance & Tolerance in Salton Seawater",subtitle="Using CLR Transformed, Summed Gene Coverage per KO Function",xlab="PC1 [41.14%]", ylab="PC2 [9.04%]",color="Sample Date")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(met.pcoa.clr.meta$SampDate_Color[order(met.pcoa.clr.meta$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("PC1 [61.21%]") + ylab("PC2 [19.15%]")

ggsave(pcoa.m1,filename = "figures/MGM_Figs/SSW_MetalResistTolerance_pcoa_CLR_SummedCoverage_Per_KO_sampdate.png", width=12, height=10, dpi=600)

# sample month shape, depth color
pcoa.m2<-ggplot(met.pcoa.clr.meta, aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampleMonth), size=5)+theme_bw()+
  labs(title="PCoA: Metal Resistance & Tolerance in Salton Seawater",subtitle="Using CLR Transformed, Summed Gene Coverage per KO Function",xlab="PC1", ylab="PC2",color="Depth (m)")+
  theme_classic()+ theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),axis.text = element_text(size=12),axis.text.x = element_text(vjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17))+
  scale_color_continuous(low="blue3",high="red",trans = 'reverse') + scale_shape_discrete(labels=c("August 2021","December 2021","April 2022"),name="Sample Date") +
  xlab("PC1 [61.21%]") + ylab("PC2 [19.15]")

ggsave(pcoa.m2,filename = "figures/MGM_Figs/SSW_MetalResistTolerance_pcoa_CLR_SummedCoverage_Per_KO.traits_depth.png", width=12, height=10, dpi=600)

#### Homogeneity of Variance (CLR data only)- Composition by Groups ####
## betadisper to look at homogeneity of group dispersions (aka variance) when considering multiple variables
# multivariate analogue to Levene's test of homogeneity of variances
# program finds spatial median or centroid of the group, & compare distances of group to centroid/spatial median via ANOVA

#While PERMANOVA tests differences in group means (analogous to MANOVA),
## a related test called PERMDISP can be used to evaluate homogeneity of group dispersion
#(analogous to Levene's test for equal variances). The vegan function for this test is “betadisper”:
## * need a distance matrix!

mgm.clr[1:4,1:4] # sample IDs are rows, genes are columns
ko.cov.sum_table[1:4,1:4] # sanity check

# check rownames of CLR & VST transformed feature count data & metadata
rownames(meta_scaled) %in% rownames(mgm.clr) #mgm.clr was used to make the distance matrix b.euc_dist

# calculate our Euclidean distance matrix using CLR data
mgm.euc.clr_dist <- dist(mgm.clr, method = "euclidean")

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
png('figures/MGM_Figs/SSW_MGM_pcoa_CLR_SummedCoverage_perKO_betadispersion_sampledate.png',width = 700, height = 600, res=100)
plot(mgm.disper5,main = "Centroids and Dispersion based on Aitchison Distance (CLR Data)", col=colorset1$SampDate_Color)
dev.off()

png('figures/MGM_Figs/SSW_MGM_boxplot_CLR_SummedCoverage_perKO_centroid_distance_sampledate.png',width = 700, height = 600, res=100)
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
png('figures/MGM_Figs/ssw_mgm_pcoa_CLR_SummedCoverage_per_KO_betadispersion_depth.png',width = 700, height = 600, res=100)
plot(mgm.disper6,main = "Centroids and Dispersion based on Aitchison Distance (CLR Data)", col=colfunc(3))
dev.off()

png('figures/MGM_Figs/ssw_mgm_boxplot_CLR_SummedCoverage_per_KO_centroid_distance_depth.png',width = 700, height = 600, res=100)
boxplot(mgm.disper6,xlab="Sample Collection Depth", main = "Distance to Centroid by Category (CLR Data)", sub="Based on Aitchison Distance", col=colfunc(3))
dev.off()
## Significant differences in homogeneities can be tested using either parametric or permutational tests,
##and parametric post hoc contrasts can also be investigated:

#### Homogeneity of Variance for Specific Fxns - Composition by Groups ####
## betadisper to look at homogeneity of group dispersions (aka variance) when considering multiple variables
# multivariate analogue to Levene's test of homogeneity of variances
# program finds spatial median or centroid of the group, & compare distances of group to centroid/spatial median via ANOVA

#While PERMANOVA tests differences in group means (analogous to MANOVA),
## a related test called PERMDISP can be used to evaluate homogeneity of group dispersion
#(analogous to Levene's test for equal variances). The vegan function for this test is “betadisper”:
## * need a distance matrix!

mgm.clr[1:4,1:4] # sample IDs are rows, genes are columns
clr.cov.sum.sulf.ko[1:4,1:4] # sanity check

# check rownames of CLR & VST transformed feature count data & metadata
rownames(meta_scaled) %in% rownames(clr.cov.sum.sulf.ko) #mgm.clr was used to make the distance matrix b.euc_dist

# calculate our Euclidean distance matrix using CLR data
sulf.euc.clr_dist <- dist(clr.cov.sum.sulf.ko[,-1], method = "euclidean")

## betadisper to look at within group variance
# first by compare dispersions by sampling date
sulf.disper1<-betadisper(sulf.euc.clr_dist, mgm_meta$SampDate)
sulf.disper1

permutest(sulf.disper1, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
#Pairwise comparisons:
#  (Observed p-value below diagonal, permuted p-value above diagonal)
#              August.2021 December.2021 April.2022
#August.2021                     0.30400      0.470
#December.2021     0.27483                    0.874
#April.2022        0.44830       0.86299

anova(sulf.disper1) # p = 0.4965 --> accept the Null H, spatial medians ARE NOT significantly difference across sample dates

TukeyHSD(sulf.disper1) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other

#                                diff       lwr       upr     p adj
# December.2021-August.2021 -41.597964 -151.7919  68.59594 0.5171715
# April.2022-August.2021    -35.900953 -146.0949  74.29295 0.6036850
# April.2022-December.2021    5.697011 -104.4969 115.89091 0.9862550

# Visualize dispersions
png('figures/MGM_Figs/SSW_MGM_pcoa_CLR_SummedCoverage_perKO_SulfurOnly_betadispersion_sampledate.png',width = 700, height = 600, res=100)
plot(sulf.disper1,main = "Centroids and Dispersion based on Aitchison Distance (Sulfur CLR Data)", col=colorset1$SampDate_Color)
dev.off()

png('figures/MGM_Figs/SSW_MGM_boxplot_CLR_SummedCoverage_perKO_SulfurOnly_centroid_distance_sampledate.png',width = 700, height = 600, res=100)
boxplot(sulf.disper1,xlab="Sample Collection Date", main = "Distance to Centroid by Category (Sulfur CLR Data)", sub="Based on Aitchison Distance", col=colorset1$SampDate_Color)
dev.off()

# What about between sampling depths?
sulf.disper2<-betadisper(sulf.euc.clr_dist, mgm_meta$Depth_m)
sulf.disper2

permutest(sulf.disper2, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons

anova(sulf.disper2) # p = 0.313 --> accept the Null H, spatial medians are NOT significantly difference across sample dates

TukeyHSD(sulf.disper2) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other
#         diff       lwr      upr     p adj
# 5-0   36.65950 -31.13395 104.45295 0.2942902
# 10-0  12.78902 -55.00443  80.58247 0.8360255
# 10-5 -23.87048 -91.66393  43.92296 0.5590104

colfunc <- colorRampPalette(c("red", "blue"))
colfunc(3)

# Visualize dispersions
png('figures/MGM_Figs/ssw_mgm_pcoa_CLR_SummedCoverage_per_KO_SulfOnly_betadispersion_depth.png',width = 700, height = 600, res=100)
plot(sulf.disper2,main = "Centroids and Dispersion based on Aitchison Distance (Sulfur CLR Data)", col=colfunc(3))
dev.off()

png('figures/MGM_Figs/ssw_mgm_boxplot_CLR_SummedCoverage_per_KO_SulfOnly_centroid_distance_depth.png',width = 700, height = 600, res=100)
boxplot(sulf.disper2,xlab="Sample Collection Depth", main = "Distance to Centroid by Category (Sulfur CLR Data)", sub="Based on Aitchison Distance", col=colfunc(3))
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

# drop SampleID column from mgm.clr
mgm.clr<-mgm.clr[,!names(mgm.clr) %in% c("SampleID")]
mgm.clr[1:4,1:4]

# First make sure your data frames you're comparing are in the same exact order!!
rownames(mgm.clr) %in% rownames(meta_scaled)

meta_scaled=meta_scaled[rownames(mgm.clr),] ## reorder metadata to match order of CLR data
perm <- with(meta_scaled, how(nperm = 1000, blocks = SampDate))

pnova1<-adonis2(mgm.clr ~ DO_Percent_Local*ORP_mV*Dissolved_OrganicMatter_RFU*Depth.num*Sulfate_milliM*Sulfide_microM,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
pnova1
## none are significant

adonis2(mgm.clr ~ DO_Percent_Local*ORP_mV*Dissolved_OrganicMatter_RFU*Depth.num*Sulfate_milliM*Sulfide_microM,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm)
#         Df SumOfSqs     R2    F Pr(>F)
#Model    23    34412 0.73114 1.8918 0.4825
#Residual 16    12654 0.26886
#Total    39    47066 1.00000

pnova2<-adonis2(mgm.clr ~ DO_Percent_Local*ORP_mV*Dissolved_OrganicMatter_RFU*Sulfate_milliM*Sulfide_microM,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
pnova2
# nothing significant

adonis2(mgm.clr ~ DO_Percent_Local*ORP_mV*Dissolved_OrganicMatter_RFU*Sulfate_milliM*Sulfide_microM,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm)

pnova3<-adonis2(mgm.clr ~ ORP_mV*Dissolved_OrganicMatter_RFU*Sulfate_milliM*Sulfide_microM,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
pnova3

adonis2(mgm.clr ~ DO_Percent_Local*Temp_DegC*Dissolved_OrganicMatter_RFU*Sulfate_milliM*Sulfide_microM,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm)

pnova4<-adonis2(mgm.clr ~ Dissolved_OrganicMatter_RFU*Sulfide_microM,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
pnova4
# nothing

adonis2(mgm.clr ~ Dissolved_OrganicMatter_RFU*Sulfide_microM,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm)


### SELF REMINDER FOR R^2
### Coefficient of Determination, denoted R2 or r2
### is the proportion of the variance in the dependent variable that is predictable from the independent variable(s)

### Pseudo F stat for PERMANOVA
### pseudo F-ratio: It compares the total sum of squared dissimilarities (or ranked dissimilarities) among objects belonging to different groups to that of objects belonging to the same group.
### Larger F-ratios indicate more pronounced group separation, however, the significance of this ratio is usually of more interest than its magnitude.

#### PERMANOVAs to Env Variables Across Groups - Specific Fxns ####

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

# drop SampleID column from mgm.clr
clr.cov.sum.sulf.ko[1:4,1:4]

# First make sure your data frames you're comparing are in the same exact order!!
rownames(clr.cov.sum.sulf.ko) %in% rownames(meta_scaled)

meta_scaled=meta_scaled[rownames(clr.cov.sum.sulf.ko),] ## reorder metadata to match order of CLR data
perm <- with(meta_scaled, how(nperm = 1000, blocks = SampDate))

s.pnov0<-adonis2(clr.cov.sum.sulf.ko[,-1] ~ Depth_m,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
s.pnov0

adonis2(clr.cov.sum.sulf.ko[,-1] ~ Depth_m,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm)

s.pnov1<-adonis2(clr.cov.sum.sulf.ko[,-1] ~ ORP_mV*Dissolved_OrganicMatter_RFU*Depth.num*Sulfate_milliM*Sulfide_microM,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
s.pnov1
## none are significant

adonis2(clr.cov.sum.sulf.ko[,-1] ~ ORP_mV*Dissolved_OrganicMatter_RFU*Depth.num*Sulfate_milliM*Sulfide_microM,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm)
#         Df SumOfSqs     R2    F Pr(>F)
#Model    23    34412 0.73114 1.8918 0.4825
#Residual 16    12654 0.26886
#Total    39    47066 1.00000

s.pnov2<-adonis2(clr.cov.sum.sulf.ko[,-1] ~ Dissolved_OrganicMatter_RFU*Sulfate_milliM*Sulfide_microM,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
s.pnov2
# nothing significant

adonis2(clr.cov.sum.sulf.ko[,-1] ~ Dissolved_OrganicMatter_RFU*Sulfate_milliM*Sulfide_microM,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm)

s.pnov3<-adonis2(clr.cov.sum.sulf.ko[,-1] ~ Sulfate_milliM*Sulfide_microM,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
s.pnov3

adonis2(clr.cov.sum.sulf.ko[,-1] ~ Sulfate_milliM*Sulfide_microM,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm)

s.pnov4<-adonis2(clr.cov.sum.sulf.ko[,-1] ~ Dissolved_OrganicMatter_RFU*Sulfate_milliM,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
s.pnov4
# nothing

adonis2(clr.cov.sum.sulf.ko[,-1] ~ Dissolved_OrganicMatter_RFU*Sulfate_milliM,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm)


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

sulf.fxn.glm.fit8<-glm(formula = CLR_SumCovPerKO ~ as.numeric(Depth_m), family = Gamma, data=sulf.ko.clr.all)%>%
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
save.image("data/Metagenomes/Analysis/SSW_MGM_Fxn_BetaDiv.Rdata")
# ^ includes all data combined in object bac.dat.all, ASV table (samples are rows, ASVs are columns), mgm_meta, and an ASV count table (where ASVs are rows, not columns)
# Version Information
sessionInfo()