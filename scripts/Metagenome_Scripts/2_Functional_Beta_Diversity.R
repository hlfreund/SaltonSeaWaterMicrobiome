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

#### Load Data & See Info About Data ####
load("data/Metagenomes/Analysis/SSW_mgm_analysis.Rdata") # load Rdata to global env
load("data/Metagenomes/Analysis/SSW_MGM_FxnBetaDiv.Rdata")

head(meta_scaled)
arsen.fxns[1:4,]
ko.cov.sum_table[1:4,1:4]
head(mgm.clr.ars)

# ABOUT THE DATA:
# Before transformations (i.e., VST, CLR, etc) were done, the following was performed:
# featureCounts counted reads that mapped to genes in contigs
# Reads mapped to genes were divided by gene length for all genes across all samples
# Gene coverage was then added together for each KO ID, since multiple genes were assigned the same KO ID
# Summed coverage per KO was then transformed via median-ratio, vst, and clr

## For pathway analyses -- after gene coverage was calculated and added together per KO ID, they were added together for each pathway
## summed coverages per KO ID, then per pathway were transformed by CLR

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
mgm_fxn.cov_table[1:4,1:4] # sanity check

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
mgm_fxn.cov_table[1:4,1:4] # sanity check

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

#### Pull Out Sulfur Metabolic Fxns from CLR data ####
## heatmaps of traits of interest

mgm.clr[1:4,1:4]

# pull out sulfur functions from CLR transformed, summed coverages (summed coverage per KO)
sulf.ko<-mgm.clr[,which(colnames(mgm.clr) %in% sulfur.fxns$KO_ID)] # merge CLR data w/ S fxns found in contigs from KOFamScan
sulf.ko$SampleID<-rownames(sulf.ko)
sulf.ko.melt<-melt(sulf.ko, by="SampleID")
colnames(sulf.ko.melt)[which(names(sulf.ko.melt) == "variable")] <- "KO_ID"
colnames(sulf.ko.melt)[which(names(sulf.ko.melt) == "value")] <- "CLR_SumCovPerKO"
head(sulf.ko.melt) #sanity check

clr.sulf.ko<-merge(sulf.ko.melt,sulf.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
head(clr.sulf.ko)
colnames(clr.sulf.ko)[which(names(clr.sulf.ko) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
clr.cov.sum.sulf.ko<-as.data.frame(dcast(clr.sulf.ko, SampleID~KO_Function.KEGG, value.var="CLR_SumCovPerKO", fun.aggregate=sum)) ###just dcast, nothing is being added here!
rownames(clr.cov.sum.sulf.ko)<-clr.cov.sum.sulf.ko$SampleID
clr.cov.sum.sulf.ko[1:4,]

# sanity check
clr.cov.sum.sulf.ko$`cysJ; sulfite reductase (NADPH) flavoprotein alpha-component [EC:1.8.1.2]`[1:4]
head(clr.sulf.ko)

#### Add Sum Coverage per KO per S Pathway - then CLR transform ####
ko.cov.sum_table[1:4,1:4] # contains the sum of coverages per gene per KO -- featureCounts was normalized by gene length across samples first to get coverage, then summed up per KO ID
ko.cov.sum.m<-melt(ko.cov.sum_table, by="SampleID")
colnames(ko.cov.sum.m)[which(names(ko.cov.sum.m) == "variable")] <- "KO_ID"
colnames(ko.cov.sum.m)[which(names(ko.cov.sum.m) == "value")] <- "SumCovPerKO"
head(ko.cov.sum.m) #sanity check

# count up gene coverage by KO, by Pathway
sulf.cov.sum.m<-merge(ko.cov.sum.m,sulf.kegg,by.x=c("KO_ID"),by.y=c("KO_ID"))
head(sulf.cov.sum.m)
sulf.path.cov.sum<-as.data.frame(dcast(sulf.cov.sum.m, SampleID~Pathway, value.var="SumCovPerKO", fun.aggregate=sum)) ###
rownames(sulf.path.cov.sum)<-sulf.path.cov.sum$SampleID
sulf.path.cov.sum[1:4,1:4]

# then CLR transform
# df must have rownames are SampleIDs, columns are ASV IDs for vegan functions below\
s.path.clr<-decostand(sulf.path.cov.sum[,-1],method = "clr", pseudocount = 1) #CLR transformation
s.path.clr[1:4,1:4]

# check rownames of CLR transformed ASV data & metadata
rownames(s.path.clr) %in% rownames(meta_scaled)

## pull out all KOs in each Pathway
unique(sulf.cov.sum.m$Pathway)
assim.sulfate.red<-data.frame(KO_Function=unique(sulf.cov.sum.m$KO_Function[which(sulf.cov.sum.m$Pathway=="Assimilatory Sulfate Reduction")]))
dissim.sulfate.redox<-data.frame(KO_Function.KEGG=unique(clr.sulf.ko$KO_Function.KEGG[which(clr.sulf.ko$Pathway=="Dissimilatory Sulfate Redox")]))
mult.sulf<-data.frame(KO_Function.KEGG=unique(clr.sulf.ko$KO_Function.KEGG[which(clr.sulf.ko$Pathway=="Multiple Pathways")]))
sox.system<-data.frame(KO_Function.KEGG=unique(clr.sulf.ko$KO_Function.KEGG[which(clr.sulf.ko$Pathway=="SOX System")]))


### Sulfur Heat Maps ####
# see max & mean of summed
max(clr.cov.sum.sulf.ko[,-1])
mean(as.matrix(clr.cov.sum.sulf.ko[,-1]))

# first heat map of sulfur KOs
heatmap(as.matrix(clr.cov.sum.sulf.ko[,-1]), scale = "none")

colSums(clr.cov.sum.sulf.ko[,-1])
#clr.cov.sum.sulf.ko2 <- clr.cov.sum.sulf.ko[,which(colSums(clr.cov.sum.sulf.ko[,-1])>10)]

heatmap(as.matrix(clr.cov.sum.sulf.ko[,-1]), scale = "none")

# prep for ggplot2 heatmap
clr.sulf.ko[1:4,]
clr.sulf.all<-merge(clr.sulf.ko,meta_scaled,by="SampleID")
head(clr.sulf.all)
clr.sulf.all$SampleID = factor(clr.sulf.all$SampleID, levels=unique(clr.sulf.all$SampleID[order(clr.sulf.all$SampDate,clr.sulf.all$Depth_m)]), ordered=TRUE)
clr.sulf.all$SampDate<-gsub("\\."," ",clr.sulf.all$SampDate)
clr.sulf.all$SampDate<-factor(clr.sulf.all$SampDate, levels=c("August 2021","December 2021","April 2022"))
clr.sulf.all$Pathway<-factor(clr.sulf.all$Pathway,levels=c("Assimilatory Sulfate Reduction","Dissimilatory Sulfate Redox","Multiple Pathways","SOX System"))
clr.sulf.all$KO_Function.KEGG = factor(clr.sulf.all$KO_Function.KEGG, levels=unique(clr.sulf.all$KO_Function.KEGG[order(clr.sulf.all$Pathway)]), ordered=TRUE)

head(clr.sulf.all)

sulf.hm1a<-ggplot(clr.sulf.all, aes(SampleID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.5","-0.5"),breaks=c(1.5,0.5,-0.5)) + labs(title="Sulfur Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(sulf.hm1a,filename = "figures/MGM_Figs/Sulfur_KOFxns_MGMs_SampID_by_Function_heatmap.png", width=18, height=13, dpi=600)

sulf.hm1b<-ggplot(clr.sulf.all, aes(SampleID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.5","-0.5"),breaks=c(1.5,0.5,-0.5)) + labs(title="Sulfur Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Pathway~.,scales="free_y", space = "free")

ggsave(sulf.hm1b,filename = "figures/MGM_Figs/Sulfur_KOFxns_MGMs_SampID_by_Function_Pathway_heatmap.png", width=17, height=15, dpi=600)

sulf.hm1c<-ggplot(clr.sulf.all, aes(interaction(SampDate,Depth_m), KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.5","-0.5"),breaks=c(1.5,0.5,-0.5)) + labs(title="Sulfur Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Pathway~.,scales="free_y", space = "free")

ggsave(sulf.hm1c,filename = "figures/MGM_Figs/Sulfur_KOFxns_MGMs_SampDate_Depth_by_Function_Pathway_heatmap.png", width=15, height=18, dpi=600)

sulf.hm1d<-ggplot(clr.sulf.all, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.5","-0.5"),breaks=c(1.5,0.5,-0.5)) + labs(title="Sulfur Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate,scales="free_x", space = "free")

ggsave(sulf.hm1d,filename = "figures/MGM_Figs/Sulfur_KOFxns_MGMs_Depth_by_Function_SampDate_best_heatmap.png", width=20, height=13, dpi=600)

sulf.hm1e<-ggplot(clr.sulf.all, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.5","-0.5"),breaks=c(1.5,0.5,-0.5)) + labs(title="Sulfur Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Pathway~SampDate, scales="free", space = "free")

ggsave(sulf.hm1e,filename = "figures/MGM_Figs/Sulfur_KOFxns_MGMs_Depth_by_Function_SampDate_Pathway_best_heatmap.png", width=20, height=15, dpi=600)
#
# sulf.hm1f<-ggplot(clr.sulf.all, aes(Pathway, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.5","-0.5"),breaks=c(1.5,0.5,-0.5)) + labs(title="Sulfur Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11,face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Depth_m~SampDate,scales="free", space = "free")
#
# ggsave(sulf.hm1f,filename = "figures/MGM_Figs/Sulfur_KOFxns_MGMs_heatmap1d.png", width=18, height=18, dpi=600)
#
# sulf.hm1g<-ggplot(clr.sulf.all, aes(Pathway, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.5","-0.5"),breaks=c(1.5,0.5,-0.5)) + labs(title="Sulfur Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11,face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_wrap(.~SampDate)
#
# ggsave(sulf.hm1g,filename = "figures/MGM_Figs/Sulfur_KOFxns_MGMs_heatmap1d.png", width=18, height=18, dpi=600)

sulf.hm1e<-ggplot(clr.sulf.all[clr.sulf.all$Depth_m==0,], aes(Pathway, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.5","-0.5"),breaks=c(1.5,0.5,-0.5)) + labs(title="Sulfur Metabolism in Salton Seawater Metagenomes - 0m",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)

ggsave(sulf.hm1e,filename = "figures/MGM_Figs/Sulfur_KOFxns_Pathways_MGMs_0m_heatmap.png", width=18, height=18, dpi=600)

sulf.hm1f<-ggplot(clr.sulf.all[clr.sulf.all$Depth_m==5,], aes(Pathway, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.5","-0.5"),breaks=c(1.5,0.5,-0.5)) + labs(title="Sulfur Metabolism in Salton Seawater Metagenomes - 5m",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)

ggsave(sulf.hm1f,filename = "figures/MGM_Figs/Sulfur_KOFxns_Pathways_MGMs_5m_heatmap.png", width=18, height=18, dpi=600)

sulf.hm1g<-ggplot(clr.sulf.all[clr.sulf.all$Depth_m==10,], aes(Pathway, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.5","-0.5"),breaks=c(1.5,0.5,-0.5)) + labs(title="Sulfur Metabolism in Salton Seawater Metagenomes - 10m",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)

ggsave(sulf.hm1g,filename = "figures/MGM_Figs/Sulfur_KOFxns_Pathways_MGMs_10m_heatmap.png", width=18, height=18, dpi=600)

# pull out specific S functions
## first, SOX
clr.Sox<-clr.sulf.all[grepl('Sox', clr.sulf.all$KO_Function),] # pull out just Sox functions
clr.Sox$SampDate = factor(clr.Sox$SampDate, levels=c("August 2021","December 2021", "April 2022"))
clr.Sox$SampleID = factor(clr.Sox$SampleID, levels=unique(clr.Sox$SampleID[order(clr.Sox$SampDate,clr.Sox$Depth_m)]), ordered=TRUE)

NA %in% clr.Sox$CLR_SumCovPerKO

s.sox.hm<-ggplot(clr.Sox, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25)  +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",na.value="grey50",labels=c("1.5","0.5","-0.5"),breaks=c(1.5,0.5,-0.5)) + labs(title="SOX Functions in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO, Grouped by Bin Assigment",fill="CLR Coverage Per KO") +
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

s.R.hm1<-ggplot(clr.as.S.redox, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25)  +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",na.value="grey50",labels=c("1.5","0.5","-0.5"),breaks=c(1.5,0.5,-0.5)) + labs(title="Assimilatory Sulfuate Reduction in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO, Grouped by Bin Assigment",fill="CLR Coverage Per KO") +
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

s.RO.hm1<-ggplot(clr.dis.S.redox, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25)  +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",na.value="grey50",labels=c("1.5","0.5","-0.5"),breaks=c(1.5,0.5,-0.5)) + labs(title="Dissimilarity Sulfuate RedOx in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO, Grouped by Bin Assigment",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=15),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.border=element_blank(),panel.background = element_rect(fill = "white", colour = NA)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0)) + facet_grid(.~SampDate)

ggsave(s.RO.hm1,filename = "figures/MGM_Figs/SSW_S_DissSO4_RedOx_Contigs_bySampDate_Depth_heatmap.png", width=15, height=10, dpi=600)

#### Pull Out DOM Metabolic Fxns from CLR data ####
## heatmaps of traits of interest

mgm.clr[1:4,1:4]

# pull out sulfur functions from CLR transformed, summed coverages (summed coverage per KO)
DOM.ko<-mgm.clr[,which(colnames(mgm.clr) %in% DOM.fxns$KO_ID)] # merge CLR data w/ DOM-related fxns found in contigs from KOFamScan
DOM.ko$SampleID<-rownames(DOM.ko)
DOM.ko.melt<-melt(DOM.ko, by="SampleID")
colnames(DOM.ko.melt)[which(names(DOM.ko.melt) == "variable")] <- "KO_ID"
colnames(DOM.ko.melt)[which(names(DOM.ko.melt) == "value")] <- "CLR_SumCovPerKO"
head(DOM.ko.melt) #sanity check

clr.DOM.ko<-merge(DOM.ko.melt,dom.kegg,by.x=c("KO_ID"),by.y=c("KO_ID"))
head(clr.DOM.ko)
colnames(clr.DOM.ko)[which(names(clr.DOM.ko) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
clr.cov.sum.DOM.ko<-as.data.frame(dcast(clr.DOM.ko, SampleID~KO_Function.KEGG, value.var="CLR_SumCovPerKO", fun.aggregate=sum)) ###
rownames(clr.cov.sum.DOM.ko)<-clr.cov.sum.DOM.ko$SampleID
clr.cov.sum.DOM.ko[1:4,1:4]

#### Add Sum Coverage per KO per DOM Pathway - then CLR transform ####
ko.cov.sum_table[1:4,1:4] # contains the sum of coverages per gene per KO -- featureCounts was normalized by gene length across samples first to get coverage, then summed up per KO ID
ko.cov.sum.m<-melt(ko.cov.sum_table, by="SampleID")
colnames(ko.cov.sum.m)[which(names(ko.cov.sum.m) == "variable")] <- "KO_ID"
colnames(ko.cov.sum.m)[which(names(ko.cov.sum.m) == "value")] <- "SumCovPerKO"
head(ko.cov.sum.m) #sanity check

# count up gene coverage by KO, by Pathway
DOM.cov.sum.m<-merge(ko.cov.sum.m,dom.kegg,by.x=c("KO_ID"),by.y=c("KO_ID"))
head(DOM.cov.sum.m)
DOM.path.cov.sum<-as.data.frame(dcast(DOM.cov.sum.m, SampleID~Pathway, value.var="SumCovPerKO", fun.aggregate=sum)) ###
rownames(DOM.path.cov.sum)<-DOM.path.cov.sum$SampleID
DOM.path.cov.sum[1:4,]

unique(dom.kegg$Pathway)
names(DOM.path.cov.sum) # not all DOM pathways considered are found in contigs

# then CLR transform
# df must have rownames are SampleIDs, columns are ASV IDs for vegan functions below\
DOM.path.clr<-decostand(DOM.path.cov.sum[,-1],method = "clr", pseudocount = 1) #CLR transformation
DOM.path.clr[1:4,]

# check rownames of CLR transformed ASV data & metadata
rownames(DOM.path.clr) %in% rownames(meta_scaled)

unique(DOM.cov.sum.m$Pathway)
unique(dom.kegg$Pathway)

head(DOM.cov.sum.m)

### DOM Heat Maps ####
# see max & mean of summed
max(clr.cov.sum.DOM.ko[,-1])
mean(as.matrix(clr.cov.sum.DOM.ko[,-1]))

# first heat map of sulfur KOs
heatmap(as.matrix(clr.cov.sum.DOM.ko[,-1]), scale = "none")

colSums(clr.cov.sum.DOM.ko[,-1])

heatmap(as.matrix(clr.cov.sum.DOM.ko[,-1]), scale = "none")

# prep for ggplot2 heatmap
clr.DOM.ko[1:4,]
clr.DOM.all<-merge(clr.DOM.ko,meta_scaled,by="SampleID")
head(clr.DOM.all)
clr.DOM.all$SampleID = factor(clr.DOM.all$SampleID, levels=unique(clr.DOM.all$SampleID[order(clr.DOM.all$SampDate,clr.DOM.all$Depth_m)]), ordered=TRUE)
clr.DOM.all$SampDate<-gsub("\\."," ",clr.DOM.all$SampDate)
clr.DOM.all$SampDate<-factor(clr.DOM.all$SampDate, levels=c("August 2021","December 2021","April 2022"))
unique(clr.DOM.all$Pathway)
unique(clr.DOM.all$Cycle)
clr.DOM.all$Cycle<-factor(clr.DOM.all$Cycle, levels=c("Sulfur Cycle","Carbon Cycle","Nitrogen Cycle"))

clr.DOM.all$PathShort<-clr.DOM.all$Pathway
# vvv can only do this type of renaming if variables are characters, not factors
clr.DOM.all$PathShort[(clr.DOM.all$PathShort) == "Assimilatory Sulfate Reduction"] <- "A. Sulfate Red"
clr.DOM.all$PathShort[(clr.DOM.all$PathShort) == "Dissimilatory Sulfate Redox"] <- "D. Sulfate Redox"
clr.DOM.all$PathShort[(clr.DOM.all$PathShort) == "SOX System"] <- "SOX"
clr.DOM.all$PathShort[(clr.DOM.all$PathShort) == "Reductive Citrate Cycle"] <- "Red. Citrate Cycle"
clr.DOM.all$PathShort[(clr.DOM.all$PathShort) == "3-Hydroxypropionate Bi-cycle"] <- "3-H BC"
clr.DOM.all$PathShort[(clr.DOM.all$PathShort) == "Phosphate acetyltransferase-acetate kinase Pathway"] <- "P.A.A.K."
clr.DOM.all$PathShort[(clr.DOM.all$PathShort) == "Assimilatory Nitrate Reduction"] <- "A. Nitrate Red"

clr.DOM.all$Pathway<-factor(clr.DOM.all$Pathway,levels=c("Assimilatory Sulfate Reduction","Dissimilatory Sulfate Redox","SOX System",
                                                         "Reductive Citrate Cycle","3-Hydroxypropionate Bi-cycle","Methanogenesis","Phosphate acetyltransferase-acetate kinase Pathway",
                                                         "Denitrification","Assimilatory Nitrate Reduction","Anammox"))

clr.DOM.all$PathShort<-factor(clr.DOM.all$PathShort,levels=c("A. Sulfate Red","D. Sulfate Redox","SOX",
                                                         "Red. Citrate Cycle","3-H BC","Methanogenesis","P.A.A.K.",
                                                         "Denitrification","A. Nitrate Red","Anammox"))
clr.DOM.all$KO_Function.KEGG = factor(clr.DOM.all$KO_Function.KEGG, levels=unique(clr.DOM.all$KO_Function.KEGG[order(clr.DOM.all$Cycle,clr.DOM.all$Pathway)]), ordered=TRUE)

head(clr.DOM.all)

DOM.hm1a<-ggplot(clr.DOM.all, aes(SampleID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.5","-0.5"),breaks=c(1.5,0.5,-0.5)) + labs(title="DOM Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(DOM.hm1a,filename = "figures/MGM_Figs/DOM_KOFxns_MGMs_SampID_by_Function_heatmap.png", width=18, height=13, dpi=600)

DOM.hm1b1<-ggplot(clr.DOM.all, aes(SampleID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.5","-0.5"),breaks=c(1.5,0.5,-0.5)) + labs(title="DOM Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~.,scales="free_y", space = "free")

ggsave(DOM.hm1b1,filename = "figures/MGM_Figs/DOM_KOFxns_MGMs_SampID_by_Function_Pathway_heatmap.png", width=17, height=22, dpi=600)

DOM.hm1b2<-ggplot(clr.DOM.all, aes(SampleID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.5","-0.5"),breaks=c(1.5,0.5,-0.5)) + labs(title="DOM Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Cycle~.,scales="free_y", space = "free")

ggsave(DOM.hm1b2,filename = "figures/MGM_Figs/DOM_KOFxns_MGMs_SampID_by_Function_Cycle_heatmap.png", width=17, height=15, dpi=600)

DOM.hm1c<-ggplot(clr.DOM.all, aes(interaction(SampDate,Depth_m), KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.5","-0.5"),breaks=c(1.5,0.5,-0.5)) + labs(title="DOM Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Pathway~.,scales="free_y", space = "free")

ggsave(DOM.hm1c,filename = "figures/MGM_Figs/DOM_KOFxns_MGMs_SampDate_Depth_by_Function_Pathway_heatmap.png", width=15, height=20, dpi=600)

DOM.hm1c1<-ggplot(clr.DOM.all, aes(interaction(SampDate,Depth_m), KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.5","-0.5"),breaks=c(1.5,0.5,-0.5)) + labs(title="DOM Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Cycle+Pathway~.,scales="free_y", space = "free")

ggsave(DOM.hm1c1,filename = "figures/MGM_Figs/DOM_KOFxns_MGMs_SampDate_Depth_by_Function_Pathway_heatmap.png", width=15, height=20, dpi=600)

DOM.hm1c2<-ggplot(clr.DOM.all[which(clr.DOM.all$CLR_SumCovPerKO>0),], aes(interaction(SampDate,Depth_m), KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.5","0"),breaks=c(1.5,0.5,0)) + labs(title="DOM Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO (CLR > 0)",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Pathway~.,scales="free_y", space = "free")

ggsave(DOM.hm1c2,filename = "figures/MGM_Figs/DOM_KOFxns_MGMs_SampDate_Depth_by_Function_Pathway_HigherCLR_heatmap.png", width=15, height=18, dpi=600)

DOM.hm1c3<-ggplot(clr.DOM.all[which(clr.DOM.all$CLR_SumCovPerKO>0),], aes(interaction(SampDate,Depth_m), KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.5","0"),breaks=c(1.5,0.5,0)) + labs(title="DOM Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO (CLR > 0)",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Cycle+Pathway~.,scales="free_y", space = "free")

ggsave(DOM.hm1c3,filename = "figures/MGM_Figs/DOM_KOFxns_MGMs_SampDate_Depth_by_Function_Pathway_Cycle_HigherCLR_heatmap.png", width=15, height=18, dpi=600)

DOM.hm1d<-ggplot(clr.DOM.all, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.5","-0.5"),breaks=c(1.5,0.5,-0.5)) + labs(title="DOM Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate,scales="free_x", space = "free")

ggsave(DOM.hm1d,filename = "figures/MGM_Figs/DOM_KOFxns_MGMs_Depth_by_Function_SampDate_best_heatmap.png", width=20, height=13, dpi=600)

DOM.hm1e<-ggplot(clr.DOM.all, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.5","-0.5"),breaks=c(1.5,0.5,-0.5)) + labs(title="DOM Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Cycle~SampDate, scales="free", space = "free")

ggsave(DOM.hm1e,filename = "figures/MGM_Figs/DOM_KOFxns_MGMs_Depth_by_Function_SampDate_Pathway_best_heatmap.png", width=15, height=20, dpi=600)

DOM.hm1e2<-ggplot(clr.DOM.all, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.5","-0.5"),breaks=c(1.5,0.5,-0.5)) + labs(title="DOM Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Cycle~SampDate, scales="free", space = "free")

ggsave(DOM.hm1e2,filename = "figures/MGM_Figs/DOM_KOFxns_MGMs_Depth_by_Function_SampDate_Cycle_best_heatmap.png", width=15, height=15, dpi=600)

#
# DOM.hm1f<-ggplot(clr.DOM.all, aes(Pathway, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.5","-0.5"),breaks=c(1.5,0.5,-0.5)) + labs(title="DOM Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11,face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Depth_m~SampDate,scales="free", space = "free")
#
# ggsave(DOM.hm1f,filename = "figures/MGM_Figs/DOM_KOFxns_MGMs_heatmap1d.png", width=18, height=18, dpi=600)
#
# DOM.hm1g<-ggplot(clr.DOM.all, aes(Pathway, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.5","-0.5"),breaks=c(1.5,0.5,-0.5)) + labs(title="DOM Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11,face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_wrap(.~SampDate)
#
# ggsave(DOM.hm1g,filename = "figures/MGM_Figs/DOM_KOFxns_MGMs_heatmap1d.png", width=18, height=18, dpi=600)

DOM.hm1e<-ggplot(clr.DOM.all[clr.DOM.all$Depth_m==0,], aes(Pathway, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.5","-0.5"),breaks=c(1.5,0.5,-0.5)) + labs(title="DOM Metabolism in Salton Seawater Metagenomes - 0m",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)

ggsave(DOM.hm1e,filename = "figures/MGM_Figs/DOM_KOFxns_Pathways_MGMs_0m_heatmap.png", width=18, height=18, dpi=600)

DOM.hm1f<-ggplot(clr.DOM.all[clr.DOM.all$Depth_m==5,], aes(Pathway, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.5","-0.5"),breaks=c(1.5,0.5,-0.5)) + labs(title="DOM Metabolism in Salton Seawater Metagenomes - 5m",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)

ggsave(DOM.hm1f,filename = "figures/MGM_Figs/DOM_KOFxns_Pathways_MGMs_5m_heatmap.png", width=18, height=18, dpi=600)

DOM.hm1g<-ggplot(clr.DOM.all[clr.DOM.all$Depth_m==10,], aes(Pathway, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.5","-0.5"),breaks=c(1.5,0.5,-0.5)) + labs(title="DOM Metabolism in Salton Seawater Metagenomes - 10m",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)

ggsave(DOM.hm1g,filename = "figures/MGM_Figs/DOM_KOFxns_Pathways_MGMs_10m_heatmap.png", width=18, height=18, dpi=600)


# Subset DOM df by pathways that are more prevalent, then rerun figs
unique(clr.DOM.all$PathShort)
clr.DOM.some<-clr.DOM.all[grepl('A. Sulfate Red|D. Sulfate Redox|SOX|3-H BC|P.A.A.K.|Anammox', clr.DOM.all$PathShort),] # pull out just assimilatory sulfate reduction functions
clr.DOM.some$PathShort
#clr.DOM.some$SampDate = factor(clr.DOM.some$SampDate, levels=c("August 2021","December 2021", "April 2022"))
#clr.DOM.some$SampleID = factor(clr.DOM.some$SampleID, levels=unique(clr.DOM.some$SampleID[order(clr.DOM.some$Depth_m,clr.DOM.some$SampDate)]), ordered=TRUE)

DOM.some.hm1<-ggplot(clr.DOM.some, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.5","-0.5"),breaks=c(1.5,0.5,-0.5)) +
  labs(title="DOM Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO",caption = "Only Includes Pathways w/ Higher CLR") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold"),plot.caption = element_text(size=14, face = "italic")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~SampDate, scales="free", space = "free")

ggsave(DOM.some.hm1,filename = "figures/MGM_Figs/DOM_KOFxns_MGMs_Depth_by_Function_SampDate_HigherCov_Pathways_Only_heatmap.png", width=15, height=15, dpi=600)

DOM.some.hm2<-ggplot(clr.DOM.some, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.5","-0.5"),breaks=c(1.5,0.5,-0.5)) +
  labs(title="DOM Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO",caption = "Only Includes Pathways w/ Higher CLR") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold"),plot.caption = element_text(size=14, face = "italic")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Cycle+PathShort~SampDate, scales="free", space = "free")

ggsave(DOM.some.hm2,filename = "figures/MGM_Figs/DOM_KOFxns_MGMs_Depth_by_Function_SampDate_HigherCov_Pathways_Only_heatmap_v2.png", width=15, height=18, dpi=600)


#### Pull Out Carbon Metabolic Fxns from CLR data ####
## heatmaps of traits of interest

mgm.clr[1:4,1:4]

# pull out Carbon metabolism functions from CLR transformed, summed coverages (summed coverage per KO)
carb.ko<-mgm.clr[,which(colnames(mgm.clr) %in% carb.fxns$KO_ID)] # merge CLR data w/ carbon-related fxns found in contigs from KOFamScan
carb.ko$SampleID<-rownames(carb.ko)
carb.ko.melt<-melt(carb.ko, by="SampleID")
colnames(carb.ko.melt)[which(names(carb.ko.melt) == "variable")] <- "KO_ID"
colnames(carb.ko.melt)[which(names(carb.ko.melt) == "value")] <- "CLR_SumCovPerKO"
head(carb.ko.melt) #sanity check

clr.carb.ko<-merge(carb.ko.melt,carb.kegg,by.x=c("KO_ID"),by.y=c("KO_ID"))
head(clr.carb.ko)
colnames(clr.carb.ko)[which(names(clr.carb.ko) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
clr.cov.sum.carb.ko<-as.data.frame(dcast(clr.carb.ko, SampleID~KO_Function.KEGG, value.var="CLR_SumCovPerKO", fun.aggregate=sum)) ###
rownames(clr.cov.sum.carb.ko)<-clr.cov.sum.carb.ko$SampleID
clr.cov.sum.carb.ko[1:4,1:4]

#### Add Sum Coverage per KO per DOM Pathway - then CLR transform ####
ko.cov.sum_table[1:4,1:4] # contains the sum of coverages per gene per KO -- featureCounts was normalized by gene length across samples first to get coverage, then summed up per KO ID
ko.cov.sum.m<-melt(ko.cov.sum_table, by="SampleID")
colnames(ko.cov.sum.m)[which(names(ko.cov.sum.m) == "variable")] <- "KO_ID"
colnames(ko.cov.sum.m)[which(names(ko.cov.sum.m) == "value")] <- "SumCovPerKO"
head(ko.cov.sum.m) #sanity check

# count up gene coverage by KO, by Pathway
carb.cov.sum.m<-merge(ko.cov.sum.m,carb.kegg,by.x=c("KO_ID"),by.y=c("KO_ID"))
head(carb.cov.sum.m)
carb.path.cov.sum<-as.data.frame(dcast(carb.cov.sum.m, SampleID~Pathway, value.var="SumCovPerKO", fun.aggregate=sum)) ###
rownames(carb.path.cov.sum)<-carb.path.cov.sum$SampleID
carb.path.cov.sum[1:4,]

unique(carb.kegg$Pathway)
names(carb.path.cov.sum) # not all DOM pathways considered are found in contigs

# then CLR transform
# df must have rownames are SampleIDs, columns are ASV IDs for vegan functions below\
carb.path.clr<-decostand(carb.path.cov.sum[,-1],method = "clr", pseudocount = 1) #CLR transformation
carb.path.clr[1:4,]

# check rownames of CLR transformed ASV data & metadata
rownames(carb.path.clr) %in% rownames(meta_scaled)

unique(carb.cov.sum.m$Pathway)
unique(carb.kegg$Pathway)

head(carb.cov.sum.m)

### Carbon Heat Maps ####
# see max & mean of summed
max(clr.cov.sum.carb.ko[,-1])
mean(as.matrix(clr.cov.sum.carb.ko[,-1]))

# first heat map of sulfur KOs
heatmap(as.matrix(clr.cov.sum.carb.ko[,-1]), scale = "none")

colSums(clr.cov.sum.carb.ko[,-1])

heatmap(as.matrix(clr.cov.sum.carb.ko[,-1]), scale = "none")

# prep for ggplot2 heatmap
clr.carb.ko[1:4,]
clr.carb.all<-merge(clr.carb.ko,meta_scaled,by="SampleID")
head(clr.carb.all)
clr.carb.all$SampleID = factor(clr.carb.all$SampleID, levels=unique(clr.carb.all$SampleID[order(clr.carb.all$SampDate,clr.carb.all$Depth_m)]), ordered=TRUE)
clr.carb.all$SampDate<-gsub("\\."," ",clr.carb.all$SampDate)
clr.carb.all$SampDate<-factor(clr.carb.all$SampDate, levels=c("August 2021","December 2021","April 2022"))
unique(clr.carb.all$Pathway)
clr.carb.all<-subset(clr.carb.all, clr.carb.all$Pathway!="Multiple Pathways")
"Multiple Pathways" %in% clr.carb.all$Pathway
clr.carb.all$Pathway<-factor(clr.carb.all$Pathway,levels=c("3-Hydroxypropionate Bi-cycle","Reductive Citrate Cycle","Phosphate acetyltransferase-acetate kinase Pathway","Reductive acetyl-CoA Pathway"))
clr.carb.all$KO_Function.KEGG = factor(clr.carb.all$KO_Function.KEGG, levels=unique(clr.carb.all$KO_Function.KEGG[order(clr.carb.all$Pathway)]), ordered=TRUE)

head(clr.carb.all)

carb.hm1a<-ggplot(clr.carb.all, aes(SampleID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.5","-0.5"),breaks=c(1.5,0.5,-0.5)) + labs(title="Carbon Fixation in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(carb.hm1a,filename = "figures/MGM_Figs/Carbon_KOFxns_MGMs_SampID_by_Function_heatmap.png", width=18, height=13, dpi=600)

carb.hm1b<-ggplot(clr.carb.all, aes(SampleID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.5","-0.5"),breaks=c(1.5,0.5,-0.5)) + labs(title="Carbon Fixation in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Pathway~.,scales="free_y", space = "free")

ggsave(carb.hm1b,filename = "figures/MGM_Figs/Carbon_KOFxns_MGMs_SampID_by_Function_Pathway_heatmap.png", width=17, height=20, dpi=600)

carb.hm1c<-ggplot(clr.carb.all, aes(interaction(SampDate,Depth_m), KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.5","-0.5"),breaks=c(1.5,0.5,-0.5)) + labs(title="Carbon Fixation in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Pathway~.,scales="free_y", space = "free")

ggsave(carb.hm1c,filename = "figures/MGM_Figs/Carbon_KOFxns_MGMs_SampDate_Depth_by_Function_Pathway_heatmap.png", width=15, height=20, dpi=600)

carb.hm1d<-ggplot(clr.carb.all, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.5","-0.5"),breaks=c(1.5,0.5,-0.5)) + labs(title="Carbon Fixation in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate,scales="free_x", space = "free")

ggsave(carb.hm1d,filename = "figures/MGM_Figs/Carbon_KOFxns_MGMs_Depth_by_Function_SampDate_best_heatmap.png", width=20, height=13, dpi=600)

carb.hm1e<-ggplot(clr.carb.all, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.5","-0.5"),breaks=c(1.5,0.5,-0.5)) + labs(title="Carbon Fixation in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Pathway~SampDate, scales="free", space = "free")

ggsave(carb.hm1e,filename = "figures/MGM_Figs/Carbon_KOFxns_MGMs_Depth_by_Function_SampDate_Pathway_best_heatmap.png", width=20, height=15, dpi=600)
#
# carb.hm1f<-ggplot(clr.carb.all, aes(Pathway, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.5","-0.5"),breaks=c(1.5,0.5,-0.5)) + labs(title="Carbon Fixation in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11,face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Depth_m~SampDate,scales="free", space = "free")
#
# ggsave(carb.hm1f,filename = "figures/MGM_Figs/Carbon_KOFxns_MGMs_heatmap1d.png", width=18, height=18, dpi=600)
#
# carb.hm1g<-ggplot(clr.carb.all, aes(Pathway, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.5","-0.5"),breaks=c(1.5,0.5,-0.5)) + labs(title="Carbon Fixation in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11,face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_wrap(.~SampDate)
#
# ggsave(carb.hm1g,filename = "figures/MGM_Figs/Carbon_KOFxns_MGMs_heatmap1d.png", width=18, height=18, dpi=600)

carb.hm1e<-ggplot(clr.carb.all[clr.carb.all$Depth_m==0,], aes(Pathway, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.5","-0.5"),breaks=c(1.5,0.5,-0.5)) + labs(title="Carbon Fixation in Salton Seawater Metagenomes - 0m",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)

ggsave(carb.hm1e,filename = "figures/MGM_Figs/Carbon_KOFxns_Pathways_MGMs_0m_heatmap.png", width=18, height=18, dpi=600)

carb.hm1f<-ggplot(clr.carb.all[clr.carb.all$Depth_m==5,], aes(Pathway, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.5","-0.5"),breaks=c(1.5,0.5,-0.5)) + labs(title="Carbon Fixation in Salton Seawater Metagenomes - 5m",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)

ggsave(carb.hm1f,filename = "figures/MGM_Figs/Carbon_KOFxns_Pathways_MGMs_5m_heatmap.png", width=18, height=18, dpi=600)

carb.hm1g<-ggplot(clr.carb.all[clr.carb.all$Depth_m==10,], aes(Pathway, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.5","-0.5"),breaks=c(1.5,0.5,-0.5)) + labs(title="Carbon Fixation in Salton Seawater Metagenomes - 10m",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)

ggsave(carb.hm1g,filename = "figures/MGM_Figs/Carbon_KOFxns_Pathways_MGMs_10m_heatmap.png", width=18, height=18, dpi=600)

#### Pull Out Arsenic Metabolic Fxns from CLR data ####
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

#### Arsenic Metabolism Heat Maps ####
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
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Arsenic Functions in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),
        axis.text = element_text(size=12),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17),
        axis.ticks=element_line(size=0.4),panel.border=element_blank()) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(ars.hm1,filename = "figures/MGM_Figs/Arsenic_KOFxns_MGMs_heatmap1.png", width=18, height=15, dpi=600)

ars.hm2<-ggplot(clr.ars.all, aes(Depth_m, KO_Function, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Arsenic Functions in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),
        axis.text = element_text(size=12),axis.text.x = element_text(),legend.text = element_text(size=12),plot.title = element_text(size=17),
        axis.ticks=element_line(size=0.4),panel.border=element_blank()) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)

ggsave(ars.hm2,filename = "figures/MGM_Figs/Arsenic_KOFxns_MGMs_heatmap_better.png", width=18, height=12, dpi=600)

#### Pull Out Selenium Metabolism Fxns from CLR data ####
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

#### Selenium Metabolism Heat Maps ####

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
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Selenium Functions in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),
        axis.text = element_text(size=12),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17),
        axis.ticks=element_line(size=0.4),panel.border=element_blank()) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(sel.hm1,filename = "figures/MGM_Figs/Selenium_KOFxns_MGMs_heatmap1.png", width=18, height=15, dpi=600)

sel.hm2<-ggplot(clr.sel.all, aes(Depth_m, KO_Function, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Selenium Functions in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),
        axis.text = element_text(size=12),axis.text.x = element_text(),legend.text = element_text(size=12),plot.title = element_text(size=17),
        axis.ticks=element_line(size=0.4),panel.border=element_blank()) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)

ggsave(sel.hm2,filename = "figures/MGM_Figs/Selenium_KOFxns_MGMs_heatmap_better.png", width=18, height=12, dpi=600)

#### Pull Out Osmoprotectant/tolerance Fxns from CLR data ####
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

#### Osmoprotectant/tolerance Functions Heat Maps ####

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
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Osmoprotectant Functions in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),
        axis.text = element_text(size=12),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17),
        axis.ticks=element_line(size=0.4),panel.border=element_blank()) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(osmo.hm1,filename = "figures/MGM_Figs/OsmoProtectant_KOFxns_MGMs_heatmap1.png", width=18, height=12, dpi=600)

osmo.hm2<-ggplot(clr.osmo.all, aes(Depth_m, KO_Function, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Osmoprotectant & Tolerance Functions in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),
        axis.text = element_text(size=12),axis.text.x = element_text(),legend.text = element_text(size=12),plot.title = element_text(size=17),
        axis.ticks=element_line(size=0.4),panel.border=element_blank()) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)

ggsave(osmo.hm2,filename = "figures/MGM_Figs/OsmoProtectant_KOFxns_MGMs_heatmap_better.png", width=18, height=12, dpi=600)

#### Pull Out Metal Resistance/Tolerance Fxns from CLR data ####
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

#### Metal Resistance/Tolerance Functions Heat Maps ####

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
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Metal Resistance & Tolerance Functions in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),
        axis.text = element_text(size=12),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17),
        axis.ticks=element_line(size=0.4),panel.border=element_blank()) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(met.hm1,filename = "figures/MGM_Figs/Metal_ResistToler_KOFxns_MGMs_heatmap1.png", width=18, height=15, dpi=600)

met.hm2<-ggplot(clr.met.all, aes(Depth_m, KO_Function, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Metal Resistance & Tolerance Functions in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),
        axis.text = element_text(size=12),axis.text.x = element_text(),legend.text = element_text(size=12),plot.title = element_text(size=17),
        axis.ticks=element_line(size=0.4),panel.border=element_blank()) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)

ggsave(met.hm2,filename = "figures/MGM_Figs/Metal_ResistToler_KOFxns_MGMs_heatmap_better.png", width=18, height=15, dpi=600)

#### Sulfur Metabolism PCoA ####
## PCOA with CLR transformed data first
# calculate our Euclidean distance matrix using CLR data

# pull out sulfur functions from CLR transformed, summed coverages (summed coverage per KO)
sulf.ko<-mgm.clr[,which(colnames(mgm.clr) %in% sulfur.fxns$KO_ID)] # merge CLR data w/ S fxns found in contigs from KOFamScan
sulf.ko$SampleID<-rownames(sulf.ko)
sulf.ko.melt<-melt(sulf.ko, by="SampleID")
colnames(sulf.ko.melt)[which(names(sulf.ko.melt) == "variable")] <- "KO_ID"
colnames(sulf.ko.melt)[which(names(sulf.ko.melt) == "value")] <- "CLR_SumCovPerKO"
head(sulf.ko.melt) #sanity check

clr.sulf.ko<-merge(sulf.ko.melt,sulf.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
head(clr.sulf.ko)
colnames(clr.sulf.ko)[which(names(clr.sulf.ko) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
clr.cov.sum.sulf.ko<-as.data.frame(dcast(clr.sulf.ko, SampleID~KO_Function.KEGG, value.var="CLR_SumCovPerKO", fun.aggregate=sum)) ###
rownames(clr.cov.sum.sulf.ko)<-clr.cov.sum.sulf.ko$SampleID
clr.cov.sum.sulf.ko[1:4,1:4]

sulf.euc.clr_dist <- dist(clr.cov.sum.sulf.ko[,-1], method = "euclidean")

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
  xlab("PC1 [67.36%]") + ylab("PC2 [16.73%]")

ggsave(pcoa.s1,filename = "figures/MGM_Figs/SSW_SulfurOnly_pcoa_CLR_SummedCoverage_Per_KO_sampdate.png", width=12, height=10, dpi=600)

# sample month shape, depth color
pcoa.s2<-ggplot(sulf.pcoa.clr.meta, aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampleMonth), size=5)+theme_bw()+
  labs(title="PCoA: Sulfur Metabolsim in Salton Seawater",subtitle="Using CLR Transformed, Summed Gene Coverage per KO Function",xlab="PC1", ylab="PC2",color="Depth (m)")+
  theme_classic()+ theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),axis.text = element_text(size=12),axis.text.x = element_text(vjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17))+
  scale_color_continuous(low="blue3",high="red",trans = 'reverse') + scale_shape_discrete(labels=c("August 2021","December 2021","April 2022"),name="Sample Date") +
  xlab("PC1 [67.36%]") + ylab("PC2 [16.73%]")

ggsave(pcoa.s2,filename = "figures/MGM_Figs/SSW_SulfurOnly_pcoa_CLR_SummedCoverage_Per_KO.traits_depth.png", width=12, height=10, dpi=600)

#### DOM Metabolism PCoA ####
## PCOA with CLR transformed data first
# calculate our Euclidean distance matrix using CLR data


mgm.clr[1:4,1:4]

# pull out sulfur functions from CLR transformed, summed coverages (summed coverage per KO)
DOM.ko<-mgm.clr[,which(colnames(mgm.clr) %in% DOM.fxns$KO_ID)] # merge CLR data w/ DOM-related fxns found in contigs from KOFamScan
DOM.ko$SampleID<-rownames(DOM.ko)
DOM.ko.melt<-melt(DOM.ko, by="SampleID")
colnames(DOM.ko.melt)[which(names(DOM.ko.melt) == "variable")] <- "KO_ID"
colnames(DOM.ko.melt)[which(names(DOM.ko.melt) == "value")] <- "CLR_SumCovPerKO"
head(DOM.ko.melt) #sanity check

clr.DOM.ko<-merge(DOM.ko.melt,dom.kegg,by.x=c("KO_ID"),by.y=c("KO_ID"))
head(clr.DOM.ko)
colnames(clr.DOM.ko)[which(names(clr.DOM.ko) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
clr.cov.sum.DOM.ko<-as.data.frame(dcast(clr.DOM.ko, SampleID~KO_Function.KEGG, value.var="CLR_SumCovPerKO", fun.aggregate=sum)) ###
rownames(clr.cov.sum.DOM.ko)<-clr.cov.sum.DOM.ko$SampleID
clr.cov.sum.DOM.ko[1:4,1:4]

DOM.euc.clr_dist <- dist(clr.cov.sum.DOM.ko[,-1], method = "euclidean")

# creating our hierarcical clustering dendrogram
DOM.euc.clr_clust <- hclust(DOM.euc.clr_dist, method="ward.D2")

# let's make it a little nicer...
DOM.euc.clr_dend <- as.dendrogram(DOM.euc.clr_clust, hang=0.2)
DOM.dend_cols <- as.character(meta_scaled$SampDate_Color[order.dendrogram(DOM.euc.clr_dend)])
labels_colors(DOM.euc.clr_dend) <- DOM.dend_cols

plot(DOM.euc.clr_dend, ylab="CLR Euclidean Distance",cex = 0.5) + title(main = "Bacteria/Archaea Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 3)
legend("topright",legend = c("August 2021","December 2021","April 2022"),cex=.8,col = c("#ef781c","#03045e","#059c3f"),pch = 15, bty = "n")
# Control is dark blue ("#218380"), #Alternaria is light blue ("#73d2de")
dev.off()

# let's use our Euclidean distance matrix from before
DOM.pcoa.clr <- pcoa(DOM.euc.clr_dist) # pcoa of euclidean distance matrix = PCA of euclidean distance matrix
##save.image("data/ssw_clr.euc.dist1_3.7.23.Rdata")

# The proportion of variances explained is in its element values$Relative_eig
DOM.pcoa.clr$values

# extract principal coordinates
DOM.pcoa.clr.vectors<-data.frame(DOM.pcoa.clr$vectors)
DOM.pcoa.clr.vectors$SampleID<-rownames(DOM.pcoa.clr$vectors)

# merge pcoa coordinates w/ metadata
DOM.pcoa.clr.meta<-merge(DOM.pcoa.clr.vectors, mgm_meta, by.x="SampleID", by.y="SampleID")
DOM.pcoa.clr.meta$SampleMonth
DOM.pcoa.clr.meta$SampDate

head(DOM.pcoa.clr.meta)

DOM.pcoa.clr$values # pull out Relative (Relative_eig) variation % to add to axes labels

# create PCoA ggplot fig
pcoa.s1<-ggplot(DOM.pcoa.clr.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(SampDate)), size=4)+theme_bw()+
  labs(title="PCoA: DOM Metabolism in Salton Seawater",subtitle="Using CLR Transformed, Summed Gene Coverage per KO Function",xlab="PC1 [41.14%]", ylab="PC2 [9.04%]",color="Sample Date")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(DOM.pcoa.clr.meta$SampDate_Color[order(DOM.pcoa.clr.meta$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("PC1 [63.86%]") + ylab("PC2 [19.14%]")

ggsave(pcoa.s1,filename = "figures/MGM_Figs/SSW_DOMOnly_pcoa_CLR_SummedCoverage_Per_KO_sampdate.png", width=12, height=10, dpi=600)

# sample month shape, depth color
pcoa.s2<-ggplot(DOM.pcoa.clr.meta, aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampleMonth), size=5)+theme_bw()+
  labs(title="PCoA: DOM Metabolsim in Salton Seawater",subtitle="Using CLR Transformed, Summed Gene Coverage per KO Function",xlab="PC1", ylab="PC2",color="Depth (m)")+
  theme_classic()+ theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),axis.text = element_text(size=12),axis.text.x = element_text(vjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17))+
  scale_color_continuous(low="blue3",high="red",trans = 'reverse') + scale_shape_discrete(labels=c("August 2021","December 2021","April 2022"),name="Sample Date") +
  xlab("PC1 [63.86%]") + ylab("PC2 [19.14%]")

ggsave(pcoa.s2,filename = "figures/MGM_Figs/SSW_DOMOnly_pcoa_CLR_SummedCoverage_Per_KO.traits_depth.png", width=12, height=10, dpi=600)


#### Carbon Fixation PCoA ####
## PCOA with CLR transformed data first
# calculate our Euclidean distance matrix using CLR data


mgm.clr[1:4,1:4]

# pull out sulfur functions from CLR transformed, summed coverages (summed coverage per KO)
carb.ko<-mgm.clr[,which(colnames(mgm.clr) %in% carb.fxns$KO_ID)] # merge CLR data w/ DOM-related fxns found in contigs from KOFamScan
carb.ko$SampleID<-rownames(carb.ko)
carb.ko.melt<-melt(carb.ko, by="SampleID")
colnames(carb.ko.melt)[which(names(carb.ko.melt) == "variable")] <- "KO_ID"
colnames(carb.ko.melt)[which(names(carb.ko.melt) == "value")] <- "CLR_SumCovPerKO"
head(carb.ko.melt) #sanity check

clr.carb.ko<-merge(carb.ko.melt,carb.kegg,by.x=c("KO_ID"),by.y=c("KO_ID"))
head(clr.carb.ko)
colnames(clr.carb.ko)[which(names(clr.carb.ko) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
clr.cov.sum.carb.ko<-as.data.frame(dcast(clr.carb.ko, SampleID~KO_Function.KEGG, value.var="CLR_SumCovPerKO", fun.aggregate=sum)) ###
rownames(clr.cov.sum.carb.ko)<-clr.cov.sum.carb.ko$SampleID
clr.cov.sum.carb.ko[1:4,1:4]

carb.euc.clr_dist <- dist(clr.cov.sum.carb.ko[,-1], method = "euclidean")

# creating our hierarcical clustering dendrogram
carb.euc.clr_clust <- hclust(carb.euc.clr_dist, method="ward.D2")

# let's make it a little nicer...
carb.euc.clr_dend <- as.dendrogram(carb.euc.clr_clust, hang=0.2)
carb.dend_cols <- as.character(meta_scaled$SampDate_Color[order.dendrogram(carb.euc.clr_dend)])
labels_colors(carb.euc.clr_dend) <- carb.dend_cols

plot(carb.euc.clr_dend, ylab="CLR Euclidean Distance",cex = 0.5) + title(main = "Bacteria/Archaea Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 3)
legend("topright",legend = c("August 2021","December 2021","April 2022"),cex=.8,col = c("#ef781c","#03045e","#059c3f"),pch = 15, bty = "n")
# Control is dark blue ("#218380"), #Alternaria is light blue ("#73d2de")
dev.off()

# let's use our Euclidean distance matrix from before
carb.pcoa.clr <- pcoa(carb.euc.clr_dist) # pcoa of euclidean distance matrix = PCA of euclidean distance matrix
##save.image("data/ssw_clr.euc.dist1_3.7.23.Rdata")

# The proportion of variances explained is in its element values$Relative_eig
carb.pcoa.clr$values

# extract principal coordinates
carb.pcoa.clr.vectors<-data.frame(carb.pcoa.clr$vectors)
carb.pcoa.clr.vectors$SampleID<-rownames(carb.pcoa.clr$vectors)

# merge pcoa coordinates w/ metadata
carb.pcoa.clr.meta<-merge(carb.pcoa.clr.vectors, mgm_meta, by.x="SampleID", by.y="SampleID")
carb.pcoa.clr.meta$SampleMonth
carb.pcoa.clr.meta$SampDate

head(carb.pcoa.clr.meta)

carb.pcoa.clr$values # pull out Relative (Relative_eig) variation % to add to axes labels

# create PCoA ggplot fig
pcoa.s1<-ggplot(carb.pcoa.clr.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(SampDate)), size=4)+theme_bw()+
  labs(title="PCoA: DOM Metabolism in Salton Seawater",subtitle="Using CLR Transformed, Summed Gene Coverage per KO Function",xlab="PC1 [41.14%]", ylab="PC2 [9.04%]",color="Sample Date")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(carb.pcoa.clr.meta$SampDate_Color[order(carb.pcoa.clr.meta$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("PC1 [29.64%]") + ylab("PC2 [22.51%]")

ggsave(pcoa.s1,filename = "figures/MGM_Figs/SSW_CarbonFixationOnly_pcoa_CLR_SummedCoverage_Per_KO_sampdate.png", width=12, height=10, dpi=600)

# sample month shape, depth color
pcoa.s2<-ggplot(carb.pcoa.clr.meta, aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampleMonth), size=5)+theme_bw()+
  labs(title="PCoA: DOM Metabolsim in Salton Seawater",subtitle="Using CLR Transformed, Summed Gene Coverage per KO Function",xlab="PC1", ylab="PC2",color="Depth (m)")+
  theme_classic()+ theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),axis.text = element_text(size=12),axis.text.x = element_text(vjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17))+
  scale_color_continuous(low="blue3",high="red",trans = 'reverse') + scale_shape_discrete(labels=c("August 2021","December 2021","April 2022"),name="Sample Date") +
  xlab("PC1 [29.64%]") + ylab("PC2 [22.51%]")

ggsave(pcoa.s2,filename = "figures/MGM_Figs/SSW_CarbonFixationOnly_pcoa_CLR_SummedCoverage_Per_KO.traits_depth.png", width=12, height=10, dpi=600)

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
mgm.disper<-betadisper(mgm.euc.clr_dist, mgm_meta$SampDate)
mgm.disper

permutest(mgm.disper, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
#Pairwise comparisons:
#  (Observed p-value below diagonal, permuted p-value above diagonal)
#               August.2021 December.2021 April.2022
# August.2021                     0.35300      0.150
# December.2021     0.34305                    0.891
# April.2022        0.13953       0.91214

anova(mgm.disper) # p = 0.364 --> accept the Null H, spatial medians ARE NOT significantly difference across sample dates

TukeyHSD(mgm.disper) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other

#                                diff       lwr       upr     p adj
# December.2021-August.2021 -1.2242178 -4.159177 1.710742 0.4551913
# April.2022-August.2021    -1.3368736 -4.271833 1.598086 0.3996944
# April.2022-December.2021  -0.1126558 -3.047615 2.822304 0.9923919

# Visualize dispersions
png('figures/MGM_Figs/SSW_MGM_pcoa_CLR_SummedCoverage_perKO_betadispersion_sampledate.png',width = 700, height = 600, res=100)
plot(mgm.disper,main = "Centroids and Dispersion based on Aitchison Distance (CLR Data)", col=colorset1$SampDate_Color)
dev.off()

png('figures/MGM_Figs/SSW_MGM_boxplot_CLR_SummedCoverage_perKO_centroid_distance_sampledate.png',width = 700, height = 600, res=100)
boxplot(mgm.disper,xlab="Sample Collection Date", main = "Distance to Centroid by Category (CLR Data)", sub="Based on Aitchison Distance", col=colorset1$SampDate_Color)
dev.off()

# What about between sampling depths?
mgm.disper2<-betadisper(mgm.euc.clr_dist, mgm_meta$Depth_m)
mgm.disper2

permutest(mgm.disper2, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons

anova(mgm.disper2) # p = 0.7948 --> accept the Null H, spatial medians are NOT significantly difference across sample dates

TukeyHSD(mgm.disper2) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other
#         diff       lwr      upr     p adj
#5-0  0.4449734 -3.466870 4.356816 0.9357599
#10-0 0.8808174 -3.031026 4.792660 0.7772850
#10-5 0.4358440 -3.475999 4.347687 0.9382591

colfunc <- colorRampPalette(c("red", "blue"))
colfunc(3)

# Visualize dispersions
png('figures/MGM_Figs/ssw_mgm_pcoa_CLR_SummedCoverage_per_KO_betadispersion_depth.png',width = 700, height = 600, res=100)
plot(mgm.disper2,main = "Centroids and Dispersion based on Aitchison Distance (CLR Data)", col=colfunc(3))
dev.off()

png('figures/MGM_Figs/ssw_mgm_boxplot_CLR_SummedCoverage_per_KO_centroid_distance_depth.png',width = 700, height = 600, res=100)
boxplot(mgm.disper2,xlab="Sample Collection Depth", main = "Distance to Centroid by Category (CLR Data)", sub="Based on Aitchison Distance", col=colfunc(3))
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
perm <- with(meta_scaled, how(nperm = 1000))

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

#### PERMANOVAs to Env Variables Across Groups - Sulfur Fxns ####

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

clr.cov.sum.sulf.ko[1:4,1:4]

# First make sure your data frames you're comparing are in the same exact order!!
rownames(clr.cov.sum.sulf.ko) %in% rownames(meta_scaled)

meta_scaled=meta_scaled[rownames(clr.cov.sum.sulf.ko),] ## reorder metadata to match order of CLR data
perm <- with(meta_scaled, how(nperm = 1000))

s.pnov0<-adonis2(clr.cov.sum.sulf.ko[,-1] ~ Depth_m,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
s.pnov0

adonis2(clr.cov.sum.sulf.ko[,-1] ~ Depth_m,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm)

s.pnov1<-adonis2(clr.cov.sum.sulf.ko[,-1] ~ ORP_mV*Dissolved_OrganicMatter_RFU*Depth.num*Sulfate_milliM*Sulfide_microM*Temp_DegC,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
s.pnov1
## none are significant

adonis2(clr.cov.sum.sulf.ko[,-1] ~ ORP_mV*Dissolved_OrganicMatter_RFU*Depth.num*Sulfate_milliM*Sulfide_microM*Temp_DegC,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm)
#         Df SumOfSqs     R2    F Pr(>F)
#Model    23    34412 0.73114 1.8918 0.4825
#Residual 16    12654 0.26886
#Total    39    47066 1.00000

s.pnov2<-adonis2(clr.cov.sum.sulf.ko[,-1] ~ ORP_mV*Dissolved_OrganicMatter_RFU*Sulfate_milliM*Sulfide_microM*Temp_DegC,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
s.pnov2
# nothing significant

adonis2(clr.cov.sum.sulf.ko[,-1] ~ Dissolved_OrganicMatter_RFU*Sulfate_milliM*Sulfide_microM,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm)

s.pnov3<-adonis2(clr.cov.sum.sulf.ko[,-1] ~ Sulfate_milliM*Sulfide_microM,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
s.pnov3

adonis2(clr.cov.sum.sulf.ko[,-1] ~ Sulfate_milliM*Sulfide_microM,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm)

s.pnov4<-adonis2(clr.cov.sum.sulf.ko[,-1] ~ Dissolved_OrganicMatter_RFU*Sulfate_milliM,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
s.pnov4
# nothing

adonis2(clr.cov.sum.sulf.ko[,-1] ~ Dissolved_OrganicMatter_RFU,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm)

s.pnov5<-adonis2(clr.cov.sum.sulf.ko[,-1] ~ Dissolved_OrganicMatter_RFU*Sulfate_milliM*Temp_DegC,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
s.pnov5

# what about by S metabolic pathway?
rownames(s.path.clr) %in% rownames(meta_scaled)

head(s.path.clr)

spath1<-adonis2(s.path.clr ~ SampDate,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
spath1

spath2<-adonis2(s.path.clr ~ SampDate*Depth_m,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
spath2


fit1<-aov(s.path.clr$`Assimilatory Sulfate Reduction` ~ meta_scaled$SampDate)
#pairwise.adonis(aug21.div$Bac_Shannon_Diversity, aug21.div$Depth_m, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different

summary(fit1)
#             Df           Sum Sq Mean Sq    F value   Pr(>F)
#meta_scaled$SampDate  2 0.8403  0.4202   40.94 0.000318 ***
#Residuals             6 0.0616  0.0103
Tuk1<-TukeyHSD(fit1)
#                               diff        lwr        upr     p adj
# December.2021-August.2021  0.01600118 -0.2378062  0.2698086 0.9796528
# April.2022-August.2021    -0.64004161 -0.8938490 -0.3862342 0.0005992
# April.2022-December.2021  -0.65604278 -0.9098502 -0.4022354 0.0005231

fit2<-aov(s.path.clr$`SOX System` ~ meta_scaled$SampDate)
#pairwise.adonis(aug21.div$Bac_Shannon_Diversity, aug21.div$Depth_m, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different

summary(fit2)
#             Df           Sum Sq Mean Sq    F value   Pr(>F)
#meta_scaled$SampDate  2 0.8403  0.4202   40.94 0.000318 ***
#Residuals             6 0.0616  0.0103
Tuk1<-TukeyHSD(fit2)

### SELF REMINDER FOR R^2
### Coefficient of Determination, denoted R2 or r2
### is the proportion of the variance in the dependent variable that is predictable from the independent variable(s)

### Pseudo F stat for PERMANOVA
### pseudo F-ratio: It compares the total sum of squared dissimilarities (or ranked dissimilarities) among objects belonging to different groups to that of objects belonging to the same group.
### Larger F-ratios indicate more pronounced group separation, however, the significance of this ratio is usually of more interest than its magnitude.


#### PERMANOVAs to Env Variables Across Groups - DOM Fxns ####

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

clr.cov.sum.DOM.ko[1:4,1:4]

# First make sure your data frames you're comparing are in the same exact order!!
rownames(clr.cov.sum.DOM.ko) %in% rownames(meta_scaled)

meta_scaled=meta_scaled[rownames(clr.cov.sum.DOM.ko),] ## reorder metadata to match order of CLR data
perm <- with(meta_scaled, how(nperm = 1000))

dom.pnov0<-adonis2(clr.cov.sum.DOM.ko[,-1] ~ Depth_m,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
dom.pnov0

adonis2(clr.cov.sum.DOM.ko[,-1] ~ Depth_m,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm)

dom.pnov1<-adonis2(clr.cov.sum.DOM.ko[,-1] ~ ORP_mV*Dissolved_OrganicMatter_RFU*Depth.num*Sulfate_milliM*Sulfide_microM,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
dom.pnov1
## none are significant

adonis2(clr.cov.sum.DOM.ko[,-1] ~ ORP_mV*Dissolved_OrganicMatter_RFU*Depth.num*Sulfate_milliM*Sulfide_microM,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm)
#         Df SumOfSqs     R2    F Pr(>F)
#Model    23    34412 0.73114 1.8918 0.4825
#Residual 16    12654 0.26886
#Total    39    47066 1.00000

dom.pnov2<-adonis2(clr.cov.sum.DOM.ko[,-1] ~ ORP_mV*Dissolved_OrganicMatter_RFU*Sulfate_milliM*Sulfide_microM,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
dom.pnov2
# nothing significant

adonis2(clr.cov.sum.DOM.ko[,-1] ~ Dissolved_OrganicMatter_RFU*Sulfate_milliM*Sulfide_microM,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm)

dom.pnov3<-adonis2(clr.cov.sum.DOM.ko[,-1] ~ ORP_mV*Sulfate_milliM,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
dom.pnov3

adonis2(clr.cov.sum.DOM.ko[,-1] ~ Sulfate_milliM*Sulfide_microM,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm)

dom.pnov4<-adonis2(clr.cov.sum.DOM.ko[,-1] ~ Dissolved_OrganicMatter_RFU*Sulfate_milliM,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
dom.pnov4
# nothing

adonis2(clr.cov.sum.DOM.ko[,-1] ~ Dissolved_OrganicMatter_RFU,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm)

dom.pnov5<-adonis2(clr.cov.sum.DOM.ko[,-1] ~ ORP_mV*Temp_DegC,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
dom.pnov5
#                   Df SumOfSqs      R2      F  Pr(>F)
# ORP_mV            1   0.6345 0.17070 3.2173 0.06019 .
# Temp_DegC         1   1.5428 0.41507 7.8231 0.04167 *
#   ORP_mV:Temp_DegC  1   0.5537 0.14895 2.8074 0.45833
# Residual          5   0.9861 0.26528
# Total             8   3.7171 1.00000


### SELF REMINDER FOR R^2
### Coefficient of Determination, denoted R2 or r2
### is the proportion of the variance in the dependent variable that is predictable from the independent variable(s)

### Pseudo F stat for PERMANOVA
### pseudo F-ratio: It compares the total sum of squared dissimilarities (or ranked dissimilarities) among objects belonging to different groups to that of objects belonging to the same group.
### Larger F-ratios indicate more pronounced group separation, however, the significance of this ratio is usually of more interest than its magnitude.



#### PERMANOVAs to Env Variables Across Groups - Carbon Fxns ####

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

clr.cov.sum.carb.ko[1:4,1:4]

# First make sure your data frames you're comparing are in the same exact order!!
rownames(clr.cov.sum.carb.ko) %in% rownames(meta_scaled)

meta_scaled=meta_scaled[rownames(clr.cov.sum.carb.ko),] ## reorder metadata to match order of CLR data
perm <- with(meta_scaled, how(nperm = 1000))

c.pnov0<-adonis2(clr.cov.sum.carb.ko[,-1] ~ Depth_m,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
c.pnov0

adonis2(clr.cov.sum.carb.ko[,-1] ~ Depth_m,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm)

c.pnov1<-adonis2(clr.cov.sum.carb.ko[,-1] ~ ORP_mV*Dissolved_OrganicMatter_RFU*Depth.num*Sulfate_milliM*Sulfide_microM*Temp_DegC,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
c.pnov1
## none are significant

adonis2(clr.cov.sum.carb.ko[,-1] ~ ORP_mV*Dissolved_OrganicMatter_RFU*Depth.num*Sulfate_milliM*Sulfide_microM,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm)

c.pnov2<-adonis2(clr.cov.sum.carb.ko[,-1] ~ ORP_mV*Dissolved_OrganicMatter_RFU*Sulfate_milliM*Sulfide_microM*Temp_DegC,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
c.pnov2
# nothing significant

adonis2(clr.cov.sum.carb.ko[,-1] ~ ORP_mV*Dissolved_OrganicMatter_RFU*Sulfate_milliM*Sulfide_microM*Temp_DegC,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm)

c.pnov3<-adonis2(clr.cov.sum.carb.ko[,-1] ~ ORP_mV*Sulfate_milliM*Temp_DegC,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
c.pnov3

adonis2(clr.cov.sum.carb.ko[,-1] ~ Sulfate_milliM*Sulfide_microM,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm)

c.pnov4<-adonis2(clr.cov.sum.carb.ko[,-1] ~ Dissolved_OrganicMatter_RFU*Sulfate_milliM,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
c.pnov4
# nothing

adonis2(clr.cov.sum.carb.ko[,-1] ~ Dissolved_OrganicMatter_RFU,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm)

c.pnov5<-adonis2(clr.cov.sum.carb.ko[,-1] ~ ORP_mV*Temp_DegC,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
c.pnov5
#                   Df SumOfSqs      R2      F  Pr(>F)
# ORP_mV            1   0.6345 0.17070 3.2173 0.06019 .
# Temp_DegC         1   1.5428 0.41507 7.8231 0.04167 *
#   ORP_mV:Temp_DegC  1   0.5537 0.14895 2.8074 0.45833
# Residual          5   0.9861 0.26528
# Total             8   3.7171 1.00000


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
# clr.cov.sum.sulf.ko[which(colnames(clr.cov.sum.sulf.ko) %in% assim.sulfate.red$KO_Function.KEGG)]
# clr.cov.sum.DOM.ko[1:4,1:4]
#
# head(clr.cov.sum.sulf.ko)
# sulf.ko.clr.all<-merge(clr.sulf.ko,meta_scaled,by=c("SampleID"))
# ars.ko.clr.all<-merge(clr.ars.ko,meta_scaled,by=c("SampleID"))
# osmo.ko.clr.all<-merge(clr.osmo.ko,meta_scaled,by=c("SampleID"))

#### Linear Regression Comparisons ####
head(sulf.pcoa.clr.meta)

step1<-step(glm(formula = Axis.1 ~ ., data=sulf.pcoa.clr.meta[,c(2,16,18:19,23:25)]))
summary(step1)
# Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)
# (Intercept)     2.41148    4.35439   0.554   0.6092
# ORP_mV          0.03297    0.01002   3.291   0.0302 *
#   Temp_DegC       0.09329    0.02300   4.057   0.0154 *
#   Sulfate_milliM -0.03820    0.01905  -2.006   0.1154
# Sulfide_microM  0.18190    0.05676   3.205   0.0327 *

step2<-step(glm(formula = Axis.1 ~ ., data=DOM.pcoa.clr.meta[,c(2,16,18:19,23:25)]))
summary(step2)
#
# Coefficients:
#                             Estimate Std. Error t value Pr(>|t|)
# (Intercept)                  6.27047    5.08754   1.233    0.306
# ORP_mV                       0.01941    0.01100   1.764    0.176
# Temp_DegC                    0.06344    0.02224   2.853    0.065 .
# Dissolved_OrganicMatter_RFU -0.08206    0.07147  -1.148    0.334
# Sulfate_milliM              -0.03692    0.01830  -2.018    0.137
# Sulfide_microM               0.10999    0.06076   1.810    0.168

step3<-step(glm(formula = Axis.1 ~ ., data=carb.pcoa.clr.meta[,c(2,16,18:19,23:25)]))
summary(step3)

# Coefficients:
#                               Estimate Std. Error t value Pr(>|t|)
# (Intercept)                 -7.22530    9.58689  -0.754   0.4930
# ORP_mV                      -0.03745    0.02677  -1.399   0.2343
# Dissolved_OrganicMatter_RFU -0.30280    0.20883  -1.450   0.2207
# Sulfate_milliM               0.10231    0.03519   2.907   0.0438 *
# Sulfide_microM              -0.23002    0.15127  -1.521   0.2030

sulf.fxn.glm.fit1<-glm(formula = Axis.1 ~ DO_Percent_Local, data=sulf.pcoa.clr.meta)%>%
  adjust_pvalue(method="bonferroni")
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(sulf.fxn.glm.fit1)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)       0.0112706  0.0004664  24.166   <2e-16 ***
#DO_Percent_Local -0.0009547  0.0005092  -1.875   0.0687 .

sulf.fxn.glm.fit2<-glm(formula = Axis.1 ~ ORP_mV, data=sulf.pcoa.clr.meta)%>%
  adjust_pvalue(method="bonferroni")

# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(sulf.fxn.glm.fit2)

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)   0.0112343  0.0004731  23.745   <2e-16 ***
#ORP_mV      -0.0001512  0.0004968  -0.304    0.763

sulf.fxn.glm.fit3<-glm(formula = Axis.1 ~ Temp_DegC, family = Gamma, data=sulf.pcoa.clr.meta)%>%
  adjust_pvalue(method="bonferroni")
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(sulf.fxn.glm.fit3)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept) 0.0113225  0.0004451  25.439   <2e-16 ***
#Temp_DegC   0.0011947  0.0004861   2.458   0.0188 *

sulf.fxn.glm.fit5<-glm(formula = Axis.1 ~ Dissolved_OrganicMatter_RFU, family = Gamma, data=sulf.pcoa.clr.meta)%>%
  adjust_pvalue(method="bonferroni")
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(sulf.fxn.glm.fit5)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)                 0.0112493  0.0004708  23.896   <2e-16 ***
#Dissolved_OrganicMatter_RFU 0.0004269  0.0004659   0.916    0.365

sulf.fxn.glm.fit6<-glm(formula = Axis.1 ~ Sulfate_milliM, family = Gamma, data=sulf.pcoa.clr.meta)%>%
  adjust_pvalue(method="bonferroni")
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(sulf.fxn.glm.fit6)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)     0.0112325  0.0004772  23.536   <2e-16 ***
#Sulfate_milliM -0.0002664  0.0004893  -0.545    0.589

sulf.fxn.glm.fit7<-glm(formula = Axis.1 ~ Sulfide_microM, family = Gamma, data=sulf.pcoa.clr.meta)%>%
  adjust_pvalue(method="bonferroni")
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(sulf.fxn.glm.fit7)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)    0.0112379  0.0004723   23.80   <2e-16 ***
#Sulfide_microM 0.0002944  0.0005160    0.57    0.572

sulf.fxn.glm.fit8<-glm(formula = Axis.1 ~ as.numeric(Depth_m), family = Gamma, data=sulf.pcoa.clr.meta)%>%
  adjust_pvalue(method="bonferroni")
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(sulf.fxn.glm.fit8)

#fit1<-aov(Axis.1 ~ as.factor(Depth_m), data=sulf.pcoa.clr.meta)
#pairwise.adonis(sulf.pcoa.clr.meta$Axis.1, sulf.pcoa.clr.meta$Depth_m, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different

summary(fit1)
#Df           Sum Sq Mean Sq    F value   Pr(>F)
#Depth_m      7   4097   585.3   1.114   0.38
#Residuals   31  16294   525.6
Tuk1<-TukeyHSD(fit1)
Tuk1$Depth_m

#plot(Axis.1 ~ Depth_m, data=sulf.pcoa.clr.meta)
#abline(aov(DustComplexity ~ Elevation, data=sulf.pcoa.clr.meta))

# fit.0<-aov(DustComplexity ~ as.factor(Elevation), data=sulf.pcoa.clr.meta)
# summary(fit.0)
# TukeyHSD(fit.0)
# Levene's test with one independent variable
## Levene's tests whether variances of 2 samples are equal
## we want variances to be the same -- want NON SIGNIFICANCE!
## t test assumes that variances are the same, so Levene's test needs to be non significant
## Fligner's test is a Levene's test for data that are not normally distributed
## more here: https://www.geeksforgeeks.org/fligner-killeen-test-in-r-programming/
#fligner.test(Axis.1 ~ Depth_m, data = sulf.pcoa.clr.meta)
# Fligner-Killeen:med chi-squared = 4.091, df = 7, p-value = 0.7692
# Which shows that the data do not deviate significantly from homogeneity.
#compare_means(Axis.1 ~ Depth_m, data=sulf.pcoa.clr.meta, method="anova",p.adjust.method = "bonferroni") # won't take as.factor(Elevation) as input

### Export Global Env for Other Scripts ####
save.image("data/Metagenomes/Analysis/SSW_MGM_Fxn_BetaDiv.Rdata")
# ^ includes all data combined in object bac.dat.all, ASV table (samples are rows, ASVs are columns), mgm_meta, and an ASV count table (where ASVs are rows, not columns)
# Version Information
sessionInfo()