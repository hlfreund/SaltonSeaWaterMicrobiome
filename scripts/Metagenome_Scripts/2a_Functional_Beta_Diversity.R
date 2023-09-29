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
  library(pairwiseAdonis)
})

#### Load Data & See Info About Data ####
load("data/Metagenomes/Analysis/SSW_mgm_analysis.Rdata") # load Rdata to global env
#load("data/Metagenomes/Analysis/SSW_MGM_FxnBetaDiv.Rdata")

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

# NOTE about CLR transformation:
## uses a pseudocount of 1 to replace 0s, which is why not all 0s are treated equally
## need to look into robustCLR, which uses CLR transformation without 0s. Need more info on this methodology...

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

ggsave(pcoa3,filename = "figures/MGM_Figs/FxnDiv/PCoAs/MedianRatioTransformation/SSW_MGM_pcoa_MR_sampdate.png", width=12, height=10, dpi=600)

# sample month shape, depth color
pcoa4<-ggplot(mgm.pcoa.mr.meta, aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampleMonth), size=5)+theme_bw()+
  labs(title="PCoA: Metagenome Functions in Salton Seawater",subtitle="Using Median-Ratio Transformed Feature Data",xlab="PC1", ylab="PC2",color="Depth (m)")+
  theme_classic()+ theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),axis.text = element_text(size=12),axis.text.x = element_text(vjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17))+
  scale_color_continuous(low="blue3",high="red",trans = 'reverse') + scale_shape_discrete(labels=c("August 2021","December 2021","April 2022"),name="Sample Date") +
  xlab("PC1 [33.04%]") + ylab("PC2 [29.24%]")

ggsave(pcoa4,filename = "figures/MGM_Figs/FxnDiv/PCoAs/MedianRatioTransformation/SSW_MGM_pcoa_MR.traits_depth.png", width=12, height=10, dpi=600)

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
png('figures/MGM_Figs/FxnDiv/PCoAs/MedianRatioTransformation/SSW_MGM_pcoa_MR_betadispersion_sampledate.png',width = 700, height = 600, res=100)
plot(mgm.disper1,main = "Centroids and Dispersion (Median-Ratio Data)", col=colorset1$SampDate_Color)
dev.off()

png('figures/MGM_Figs/FxnDiv/PCoAs/MedianRatioTransformation/SSW_MGM_boxplot_MR_centroid_distance_sampledate.png',width = 700, height = 600, res=100)
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
png('figures/MGM_Figs/FxnDiv/PCoAs/MedianRatioTransformation/ssw_mgm_pcoa_MR_betadispersion_depth.png',width = 700, height = 600, res=100)
plot(mgm.disper2,main = "Centroids and Dispersion (Median-Ratio Data)", col=colfunc(3))
dev.off()

png('figures/MGM_Figs/FxnDiv/PCoAs/MedianRatioTransformation/ssw_mgm_boxplot_MR_centroid_distance_depth.png',width = 700, height = 600, res=100)
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

ggsave(pcoa3,filename = "figures/MGM_Figs/FxnDiv/PCoAs/VarianceStabilizingTransformation/SSW_MGM_pcoa_VST_sampdate.png", width=12, height=10, dpi=600)

# sample month shape, depth color
pcoa4<-ggplot(mgm.pcoa.vst.meta, aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampleMonth), size=5)+theme_bw()+
  labs(title="PCoA: Metagenome Functions in Salton Seawater",subtitle="Using Variance Stabilization Transformed Feature Data",xlab="PC1", ylab="PC2",color="Depth (m)")+
  theme_classic()+ theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),axis.text = element_text(size=12),axis.text.x = element_text(vjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17))+
  scale_color_continuous(low="blue3",high="red",trans = 'reverse') + scale_shape_discrete(labels=c("August 2021","December 2021","April 2022"),name="Sample Date") +
  xlab("PC1 [25.06%]") + ylab("PC2 [21.98%]")

ggsave(pcoa4,filename = "figures/MGM_Figs/FxnDiv/PCoAs/VarianceStabilizingTransformation/SSW_MGM_pcoa_VST.traits_depth.png", width=12, height=10, dpi=600)

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
png('figures/MGM_Figs/FxnDiv/PCoAs/VarianceStabilizingTransformation/SSW_MGM_pcoa_vst_betadispersion_sampledate.png',width = 700, height = 600, res=100)
plot(mgm.disper3,main = "Centroids and Dispersion (VST Data)", col=colorset1$SampDate_Color)
dev.off()

png('figures/MGM_Figs/FxnDiv/PCoAs/VarianceStabilizingTransformation/SSW_MGM_boxplot_vst_centroid_distance_sampledate.png',width = 700, height = 600, res=100)
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
png('figures/MGM_Figs/FxnDiv/PCoAs/VarianceStabilizingTransformation/ssw_mgm_pcoa_VST_betadispersion_depth.png',width = 700, height = 600, res=100)
plot(mgm.disper2,main = "Centroids and Dispersion (VST Data)", col=colfunc(3))
dev.off()

png('figures/MGM_Figs/FxnDiv/PCoAs/VarianceStabilizingTransformation/ssw_mgm_boxplot_VST_centroid_distance_depth.png',width = 700, height = 600, res=100)
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

ggsave(pcoa5,filename = "figures/MGM_Figs/FxnDiv/PCoAs/CenterLogRatioTransformation/SSW_MGM_pcoa_CLR_SummedCoverage_Per_KO_sampdate.png", width=12, height=10, dpi=600)

# sample month shape, depth color
pcoa6<-ggplot(mgm.pcoa.clr.meta, aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampleMonth), size=5)+theme_bw()+
  labs(title="PCoA: Metagenome Functions in Salton Seawater",subtitle="Using CLR Transformed, Summed Gene Coverage per KO Function",xlab="PC1", ylab="PC2",color="Depth (m)")+
  theme_classic()+ theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),axis.text = element_text(size=12),axis.text.x = element_text(vjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17))+
  scale_color_continuous(low="blue3",high="red",trans = 'reverse') + scale_shape_discrete(labels=c("August 2021","December 2021","April 2022"),name="Sample Date") +
  xlab("PC1 [34.06%]") + ylab("PC2 [22.48%]")

ggsave(pcoa6,filename = "figures/MGM_Figs/FxnDiv/PCoAs/CenterLogRatioTransformation/SSW_MGM_pcoa_CLR_SummedCoverage_Per_KO.traits_depth.png", width=12, height=10, dpi=600)

#### Functional Beta Diversity - Robust CLR data ####
mgm.Rclr[1:4,1:4] # sample IDs are rows, genes are columns
ko.cov.sum_table[1:4,1:4] # sanity check

# check rownames of Robust CLR & VST transformed feature count data & metadata
rownames(mgm.Rclr) %in% rownames(meta_scaled)

## PCOA with Robust CLR transformed data first
# calculate our Euclidean distance matrix using Robust CLR data
mgm.euc.Rclr_dist <- dist(mgm.Rclr, method = "euclidean")

# creating our hierarcical clustering dendrogram
mgm.euc.Rclr_clust <- hclust(mgm.euc.Rclr_dist, method="ward.D2")

# let's make it a little nicer...
mgm.euc.Rclr_dend <- as.dendrogram(mgm.euc.Rclr_clust, hang=0.2)
mgm.dend_cols <- as.character(meta_scaled$SampDate_Color[order.dendrogram(mgm.euc.Rclr_dend)])
labels_colors(mgm.euc.Rclr_dend) <- mgm.dend_cols

plot(mgm.euc.Rclr_dend, ylab="Robust CLR Euclidean Distance",cex = 0.5) + title(main = "Bacteria/Archaea Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 3)
legend("topright",legend = c("August 2021","December 2021","April 2022"),cex=.8,col = c("#ef781c","#03045e","#059c3f"),pch = 15, bty = "n")
# Control is dark blue ("#218380"), #Alternaria is light blue ("#73d2de")
dev.off()

# let's use our Euclidean distance matrix from before
mgm.pcoa.Rclr <- pcoa(mgm.euc.Rclr_dist) # pcoa of euclidean distance matrix = PCA of euclidean distance matrix
##save.image("data/ssw_clr.euc.dist1_3.7.23.Rdata")

# The proportion of variances explained is in its element values$Relative_eig
mgm.pcoa.Rclr$values

# extract principal coordinates
mgm.pcoa.Rclr.vectors<-data.frame(mgm.pcoa.Rclr$vectors)
mgm.pcoa.Rclr.vectors$SampleID<-rownames(mgm.pcoa.Rclr$vectors)

# merge pcoa coordinates w/ metadata
mgm.pcoa.Rclr.meta<-merge(mgm.pcoa.Rclr.vectors, mgm_meta, by.x="SampleID", by.y="SampleID")
mgm.pcoa.Rclr.meta$SampleMonth
mgm.pcoa.Rclr.meta$SampDate

head(mgm.pcoa.Rclr.meta)

mgm.pcoa.Rclr$values # pull out Relative (Relative_eig) variation % to add to axes labels

# create PCoA ggplot fig
pcoa5<-ggplot(mgm.pcoa.Rclr.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(SampDate)), size=4)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Seawater",subtitle="Using Robust CLR Transformed, Summed Gene Coverage per KO Function",xlab="PC1 [41.14%]", ylab="PC2 [9.04%]",color="Sample Date")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(mgm.pcoa.Rclr.meta$SampDate_Color[order(mgm.pcoa.Rclr.meta$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("PC1 [34.00%]") + ylab("PC2 [16.63%]")

ggsave(pcoa5,filename = "figures/MGM_Figs/FxnDiv/PCoAs/RobustCenterLogRatioTransformation/SSW_MGM_pcoa_Robust CLR_SummedCoverage_Per_KO_sampdate.png", width=12, height=10, dpi=600)

# sample month shape, depth color
pcoa6<-ggplot(mgm.pcoa.Rclr.meta, aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampleMonth), size=5)+theme_bw()+
  labs(title="PCoA: Metagenome Functions in Salton Seawater",subtitle="Using Robust CLR Transformed, Summed Gene Coverage per KO Function",xlab="PC1", ylab="PC2",color="Depth (m)")+
  theme_classic()+ theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),axis.text = element_text(size=12),axis.text.x = element_text(vjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17))+
  scale_color_continuous(low="blue3",high="red",trans = 'reverse') + scale_shape_discrete(labels=c("August 2021","December 2021","April 2022"),name="Sample Date") +
  xlab("PC1 [34.00%]") + ylab("PC2 [16.63%]")

ggsave(pcoa6,filename = "figures/MGM_Figs/FxnDiv/PCoAs/RobustCenterLogRatioTransformation/SSW_MGM_pcoa_Robust CLR_SummedCoverage_Per_KO.traits_depth.png", width=12, height=10, dpi=600)


#### Pull Out Sulfur Metabolic Fxns from CLR data - NO NAs ####
## heatmaps of traits of interest

mgm.clr[1:4,1:4]

# pull out sulfur functions from CLR transformed, summed coverages (summed gene coverage per KO)
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

# Find highest (normalized) coverage fxns per sample per timepoint
aug.sulf<-clr.sulf.ko[grep("8.24.21.",clr.sulf.ko$SampleID),]
aug.sulf[which.max(aug.sulf$CLR_SumCovPerKO),] # Aug 21, 0m, 1.727594 cysD; sulfate adenylyltransferase subunit 2 [EC:2.7.7.4] Assimilatory Sulfate Reduction
aug.sulf[order(aug.sulf$CLR_SumCovPerKO,decreasing=TRUE),]

dec.sulf<-clr.sulf.ko[grep("12.22.21.",clr.sulf.ko$SampleID),]
dec.sulf[which.max(dec.sulf$CLR_SumCovPerKO),] # Dec 21, 10m, 1.788872 cysD; sulfate adenylyltransferase subunit 2 [EC:2.7.7.4] Assimilatory Sulfate Reduction
dec.sulf[order(dec.sulf$CLR_SumCovPerKO,decreasing=TRUE),]

apr.sulf<-clr.sulf.ko[grep("4.13.22.",clr.sulf.ko$SampleID),]
apr.sulf[which.max(apr.sulf$CLR_SumCovPerKO),] # Apr 22, 5m, 1.853456 sqr; sulfide:quinone oxidoreductase [EC:1.8.5.4] Sulfide Oxidation
apr.sulf[order(apr.sulf$CLR_SumCovPerKO,decreasing=TRUE),]

#### Pull Out Sulfur Metabolic Fxns from CLR data - with NAs ####
## heatmaps of traits of interest

mgm.clr.na[1:4,1:4]

# pull out sulfur functions from CLR transformed, summed coverages (summed gene coverage per KO)
sulf.ko.na<-mgm.clr.na[,which(colnames(mgm.clr.na) %in% sulfur.fxns$KO_ID)] # merge CLR data w/ S fxns found in contigs from KOFamScan
sulf.ko.na$SampleID<-rownames(sulf.ko.na)
sulf.ko.na.melt<-melt(sulf.ko.na, by="SampleID")
colnames(sulf.ko.na.melt)[which(names(sulf.ko.na.melt) == "variable")] <- "KO_ID"
colnames(sulf.ko.na.melt)[which(names(sulf.ko.na.melt) == "value")] <- "CLR_SumCovPerKO"
head(sulf.ko.na.melt) #sanity check

clr.sulf.ko.na<-merge(sulf.ko.na.melt,sulf.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
head(clr.sulf.ko.na)
colnames(clr.sulf.ko.na)[which(names(clr.sulf.ko.na) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
clr.cov.sum.sulf.ko.na<-as.data.frame(dcast(clr.sulf.ko.na, SampleID~KO_Function.KEGG, value.var="CLR_SumCovPerKO", fun.aggregate=sum)) ###just dcast, nothing is being added here!
rownames(clr.cov.sum.sulf.ko.na)<-clr.cov.sum.sulf.ko.na$SampleID
clr.cov.sum.sulf.ko.na[1:4,]

# sanity check
clr.cov.sum.sulf.ko.na$`cysJ; sulfite reductase (NADPH) flavoprotein alpha-component [EC:1.8.1.2]`[1:4]
head(clr.sulf.ko.na)

#### Pull out CLR Cov Per S Genes in S Pathways ####
ko.cov.sum_table[1:4,1:4] # contains the sum of coverages per gene per KO -- featureCounts was normalized by gene length across samples first to get coverage, then summed up per KO ID

## pull out all KOs in each Pathway
unique(clr.sulf.ko$Pathway)
head(clr.sulf.ko)
clr.cov.sum.sulf.ko[1:4,]

assim.sulfate.red<-data.frame(KO_Function.KEGG=unique(clr.sulf.ko$KO_Function[which(clr.sulf.ko$Pathway=="Assimilatory Sulfate Reduction")]))
dissim.sulfate.redox<-data.frame(KO_Function.KEGG=unique(clr.sulf.ko$KO_Function.KEGG[which(clr.sulf.ko$Pathway=="Dissimilatory Sulfate Redox")]))
mult.sulf<-data.frame(KO_Function.KEGG=unique(clr.sulf.ko$KO_Function.KEGG[which(clr.sulf.ko$Pathway=="Multiple Pathways")]))
sox.system<-data.frame(KO_Function.KEGG=unique(clr.sulf.ko$KO_Function.KEGG[which(clr.sulf.ko$Pathway=="SOX")]))

# pull out functions & CLR info per pathway
asSO4.ko.cov<-clr.cov.sum.sulf.ko[,-1][,colnames(clr.cov.sum.sulf.ko[,-1]) %in% assim.sulfate.red$KO_Function.KEGG] # pull out sox genes from gene list found in CLR transformed cov per KO
disSO4.ko.cov<-clr.cov.sum.sulf.ko[,-1][,colnames(clr.cov.sum.sulf.ko[,-1]) %in% dissim.sulfate.redox$KO_Function.KEGG] # pull out sox genes from gene list found in CLR transformed cov per KO
multiS.ko.cov<-clr.cov.sum.sulf.ko[,-1][,colnames(clr.cov.sum.sulf.ko[,-1]) %in% mult.sulf$KO_Function.KEGG] # pull out sox genes from gene list found in CLR transformed cov per KO
sox.ko.cov<-clr.cov.sum.sulf.ko[,-1][,colnames(clr.cov.sum.sulf.ko[,-1]) %in% sox.system$KO_Function.KEGG] # pull out sox genes from gene list found in CLR transformed cov per KO

# * NOTE - using summed coverage per KO, CLR transformed data --> then using info about genes in each pathway to pull out specific gene CLRs to compare pathways
# we are NOT summering coverage per KO, then per pathway, then CLR transforming. Note added 7/13/23

#### Sulfur Heat Maps ####
# see max & mean of summed
max(clr.cov.sum.sulf.ko[,-1])
min(clr.cov.sum.sulf.ko[,-1])
max(clr.cov.sum.sulf.ko[,-1])/2

# first heat map of sulfur KOs
heatmap(as.matrix(clr.cov.sum.sulf.ko[,-1]), scale = "none")

colSums(clr.cov.sum.sulf.ko[,-1])
#clr.cov.sum.sulf.ko2 <- clr.cov.sum.sulf.ko[,which(colSums(clr.cov.sum.sulf.ko[,-1])>10)]

heatmap(as.matrix(clr.cov.sum.sulf.ko[,-1]), scale = "none")

# prep for ggplot2 heatmap -- using CLR data that includes NAs so they are blocked out on heatmaps
clr.sulf.ko.na[1:4,]
clr.sulf.all<-merge(clr.sulf.ko.na,meta_scaled,by="SampleID")
head(clr.sulf.all)
clr.sulf.all$PlotID = factor(clr.sulf.all$PlotID, levels=unique(clr.sulf.all$PlotID[order(clr.sulf.all$SampDate,clr.sulf.all$Depth_m)]), ordered=TRUE)
clr.sulf.all$SampDate<-gsub("\\."," ",clr.sulf.all$SampDate)
clr.sulf.all$SampDate<-factor(clr.sulf.all$SampDate, levels=c("August 2021","December 2021","April 2022"))

clr.sulf.all$PathShort<-clr.sulf.all$Pathway
clr.sulf.all$PathShort[(clr.sulf.all$PathShort) == "Dissimilatory Sulfate Redox"] <- "D.SO4 RedOx"
clr.sulf.all$PathShort[(clr.sulf.all$PathShort) == "Assimilatory Sulfate Reduction"] <- "A.SO4 Red"
clr.sulf.all$PathShort[(clr.sulf.all$PathShort) == "Multiple Pathways"] <- "Multi Paths"
clr.sulf.all$PathShort[(clr.sulf.all$PathShort) == "S Disproportionation"] <- "S Disprop."

clr.sulf.all$Pathway<-factor(clr.sulf.all$Pathway,levels=c("Assimilatory Sulfate Reduction","Dissimilatory Sulfate Redox","Multiple Pathways","SOX","S Disproportionation"))
clr.sulf.all$PathShort<-factor(clr.sulf.all$PathShort,levels=c("A.SO4 Red","D.SO4 RedOx","Multi Paths","SOX","S Disprop."))

clr.sulf.all$PathSpecShort<-clr.sulf.all$PathwaySpecific
clr.sulf.all$PathSpecShort[(clr.sulf.all$PathSpecShort) == "Dissimilatory Sulfate Redox"] <- "D.SO4 RedOx"
clr.sulf.all$PathSpecShort[(clr.sulf.all$PathSpecShort) == "Assimilatory Sulfate Reduction"] <- "A.SO4 Red"
clr.sulf.all$PathSpecShort[(clr.sulf.all$PathSpecShort) == "Multiple Pathways"] <- "Multi"
clr.sulf.all$PathSpecShort[(clr.sulf.all$PathSpecShort) == "Sulfur Disproportionation"] <- "Dispro"
clr.sulf.all$PathSpecShort[(clr.sulf.all$PathSpecShort) == "Sulfide Oxidation"] <- "H2S Ox"
clr.sulf.all$PathSpecShort[(clr.sulf.all$PathSpecShort) == "Sulfite Oxidation"] <- "SO3 Ox"
clr.sulf.all$PathSpecShort[(clr.sulf.all$PathSpecShort) == "Thiosulfate Oxidation"] <- "S2O3 Ox"

clr.sulf.all$PathwaySpecific<-factor(clr.sulf.all$PathwaySpecific,levels=c("Assimilatory Sulfate Reduction","Dissimilatory Sulfate Redox","Multiple Pathways","SOX","Sulfur Disproportionation","Sulfide Oxidation","Sulfite Oxidation","Thiosulfate Oxidation"))
clr.sulf.all$PathSpecShort<-factor(clr.sulf.all$PathSpecShort,levels=c("A.SO4 Red","D.SO4 RedOx","Multi","Dispro","H2S Ox","SO3 Ox","S2O3 Ox"))

clr.sulf.all$KO_Function.KEGG = factor(clr.sulf.all$KO_Function.KEGG, levels=unique(clr.sulf.all$KO_Function.KEGG[order(clr.sulf.all$Pathway)]), ordered=TRUE)

head(clr.sulf.all)

# For heatmap color gradient
max(clr.sulf.all$CLR_SumCovPerKO, na.rm=TRUE)
max(clr.sulf.all$CLR_SumCovPerKO, na.rm=TRUE)/2
min(clr.sulf.all$CLR_SumCovPerKO, na.rm=TRUE)

# Figures below
# by SampleID

sulf.hm1a<-ggplot(clr.sulf.all, aes(PlotID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.6","-0.3"),breaks=c(1.5,0.6,-0.3)) + labs(title="Sulfur Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(sulf.hm1a,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_SampID_by_Function_heatmap.png", width=18, height=13, dpi=600)

sulf.hm1a2<-ggplot(clr.sulf.all, aes(PlotID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.6","-0.3"),breaks=c(1.5,0.6,-0.3)) + labs(title="Sulfur Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Pathway~.,scales="free_y", space = "free")

ggsave(sulf.hm1a2,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_SampID_by_Function_Pathway_heatmap.png", width=17, height=15, dpi=600)

sulf.hm1a3<-ggplot(clr.sulf.all, aes(PlotID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.6","-0.3"),breaks=c(1.5,0.6,-0.3)) + labs(title="Sulfur Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate,scales="free_x", space = "free")

ggsave(sulf.hm1a3,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_SampleID_by_Function_SampDate_best_heatmap.png", width=20, height=13, dpi=600)

sulf.hm1a4<-ggplot(clr.sulf.all, aes(PlotID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.6","-0.3"),breaks=c(1.5,0.6,-0.3)) + labs(title="Sulfur Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~.,scales="free_y", space = "free")

ggsave(sulf.hm1a4,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_SampID_by_Function_Pathway_heatmap2.png", width=17, height=15, dpi=600)

sulf.hm1a5<-ggplot(clr.sulf.all, aes(PlotID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.6","-0.3"),breaks=c(1.5,0.6,-0.3)) + labs(title="Sulfur Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Pathway~SampDate, scales="free", space = "free")

ggsave(sulf.hm1a5,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_SampleID_by_Function_SampDate_Pathway_heatmap.png", width=20, height=15, dpi=600)

sulf.hm1a6<-ggplot(clr.sulf.all, aes(PlotID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.6","-0.3"),breaks=c(1.5,0.6,-0.3)) + labs(title="Sulfur Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~SampDate, scales="free", space = "free")

ggsave(sulf.hm1a6,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_SampleID_by_Function_SampDate_Pathway_best_heatmap2.png", width=20, height=15, dpi=600)

sulf.hm1a7<-ggplot(clr.sulf.all, aes(PlotID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.6","-0.3"),breaks=c(1.5,0.6,-0.3)) + labs(title="Sulfur Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathwaySpecific~SampDate, scales="free", space = "free")

ggsave(sulf.hm1a7,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_SampleID_by_Function_SampDate_PathwaySpecific_best_heatmap.png", width=20, height=20, dpi=600)

sulf.hm1a8<-ggplot(clr.sulf.all, aes(PlotID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.6","-0.3"),breaks=c(1.5,0.6,-0.3)) + labs(title="Sulfur Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathSpecShort~SampDate, scales="free", space = "free")

ggsave(sulf.hm1a8,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_SampleID_by_Function_SampDate_PathwaySpecific_best_heatmap2.png", width=20, height=20, dpi=600)


# sulf.hm1c<-ggplot(clr.sulf.all, aes(interaction(SampDate,Depth_m), KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.6","-0.3"),breaks=c(1.5,0.6,-0.3)) + labs(title="Sulfur Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Pathway~.,scales="free_y", space = "free")
#
# ggsave(sulf.hm1c,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_SampDate_Depth_by_Function_Pathway_heatmap.png", width=15, height=18, dpi=600)
#
# sulf.hm1c2<-ggplot(clr.sulf.all, aes(interaction(SampDate,Depth_m), KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.6","-0.3"),breaks=c(1.5,0.6,-0.3)) + labs(title="Sulfur Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~.,scales="free_y", space = "free")
#
# ggsave(sulf.hm1c2,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_SampDate_Depth_by_Function_Pathway_heatmap2.png", width=15, height=18, dpi=600)

# by Depth
sulf.hm1b<-ggplot(clr.sulf.all, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.6","-0.3"),breaks=c(1.5,0.6,-0.3)) + labs(title="Sulfur Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(sulf.hm1b,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Depth_by_Function_heatmap.png", width=18, height=13, dpi=600)

sulf.hm1b2<-ggplot(clr.sulf.all, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.6","-0.3"),breaks=c(1.5,0.6,-0.3)) + labs(title="Sulfur Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Pathway~.,scales="free_y", space = "free")

ggsave(sulf.hm1b2,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Depth_by_Function_Pathway_heatmap.png", width=17, height=15, dpi=600)

sulf.hm1b3<-ggplot(clr.sulf.all, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.6","-0.3"),breaks=c(1.5,0.6,-0.3)) + labs(title="Sulfur Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate,scales="free_x", space = "free")

ggsave(sulf.hm1b3,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Depth_by_Function_SampDate_heatmap.png", width=20, height=13, dpi=600)

sulf.hm1b4<-ggplot(clr.sulf.all, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.6","-0.3"),breaks=c(1.5,0.6,-0.3)) + labs(title="Sulfur Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~.,scales="free_y", space = "free")

ggsave(sulf.hm1b4,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Depth_by_Function_Pathway_heatmap2.png", width=17, height=15, dpi=600)

sulf.hm1b5<-ggplot(clr.sulf.all, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.6","-0.3"),breaks=c(1.5,0.6,-0.3)) + labs(title="Sulfur Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Pathway~SampDate, scales="free", space = "free")

ggsave(sulf.hm1b5,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Depth_by_Function_SampDate_Pathway_heatmap.png", width=20, height=15, dpi=600)

sulf.hm1b6<-ggplot(clr.sulf.all, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.6","-0.3"),breaks=c(1.5,0.6,-0.3)) + labs(title="Sulfur Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~SampDate, scales="free", space = "free")

ggsave(sulf.hm1b6,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Depth_by_Function_SampDate_Pathway_best_heatmap.png", width=20, height=15, dpi=600)

sulf.hm1b7<-ggplot(clr.sulf.all, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.6","-0.3"),breaks=c(1.5,0.6,-0.3)) + labs(title="Sulfur Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathwaySpecific~SampDate, scales="free", space = "free")

ggsave(sulf.hm1b7,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Depth_by_Function_SampDate_PathwaySpecific_heatmap.png", width=20, height=20, dpi=600)

sulf.hm1b8<-ggplot(clr.sulf.all, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.6","-0.3"),breaks=c(1.5,0.6,-0.3)) + labs(title="Sulfur Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathSpecShort~SampDate, scales="free", space = "free")

ggsave(sulf.hm1b8,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Depth_by_Function_SampDate_PathwaySpecific_best_heatmap.png", width=20, height=20, dpi=600)

sulf.hm1b8a<-ggplot(clr.sulf.all, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.6","-0.3"),breaks=c(1.5,0.6,-0.3)) + labs(fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=26,face="bold"),legend.text = element_text(size=15),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),strip.text = element_text(size = 12,face="bold"),strip.text.y=element_text(face="bold",size = 17.5),strip.text.x = element_text(size = 25)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathSpecShort~SampDate, scales="free", space = "free")

ggsave(sulf.hm1b8a,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Depth_by_Function_SampDate_PathwaySpecific_best_heatmap_poster.png", width=25, height=25, dpi=600)
ggsave(sulf.hm1b8a,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_Depth_by_Function_SampDate_PathwaySpecific_best_heatmap_poster2.png", width=30, height=25, dpi=600)

# sulf.hm1e<-ggplot(clr.sulf.all[clr.sulf.all$Depth_m==0,], aes(Pathway, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.6","-0.3"),breaks=c(1.5,0.6,-0.3)) + labs(title="Sulfur Metabolism in Salton Seawater Metagenomes - 0m",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11,face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)
#
# ggsave(sulf.hm1e,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Sulfur_KOFxns_Pathways_MGMs_0m_heatmap.png", width=18, height=18, dpi=600)
#
# sulf.hm1f<-ggplot(clr.sulf.all[clr.sulf.all$Depth_m==5,], aes(Pathway, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.6","-0.3"),breaks=c(1.5,0.6,-0.3)) + labs(title="Sulfur Metabolism in Salton Seawater Metagenomes - 5m",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)
#
# ggsave(sulf.hm1f,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Sulfur_KOFxns_Pathways_MGMs_5m_heatmap.png", width=18, height=18, dpi=600)
#
# sulf.hm1g<-ggplot(clr.sulf.all[clr.sulf.all$Depth_m==10,], aes(Pathway, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.6","-0.3"),breaks=c(1.5,0.6,-0.3)) + labs(title="Sulfur Metabolism in Salton Seawater Metagenomes - 10m",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)
#
# ggsave(sulf.hm1g,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Sulfur_KOFxns_Pathways_MGMs_10m_heatmap.png", width=18, height=18, dpi=600)

#### Look at Specific S Gene Coverage Across Samples ####

# Note: Must run section "Pull Out Sulfur Metabolic Fxns from CLR data" before running this section
head(clr.cov.sum.sulf.ko.na) # columns are genes in this df

# merge with scaled metadata and prep for scatterplots of traits across samples
clr.sulf.trait.table<-merge(clr.cov.sum.sulf.ko.na,meta_scaled,by="SampleID")
head(clr.sulf.trait.table)
clr.sulf.trait.table$PlotID = factor(clr.sulf.trait.table$PlotID, levels=unique(clr.sulf.trait.table$PlotID[order(clr.sulf.trait.table$SampDate,clr.sulf.trait.table$Depth_m)]), ordered=TRUE)
clr.sulf.trait.table$SampDate<-gsub("\\."," ",clr.sulf.trait.table$SampDate)
clr.sulf.trait.table$SampDate<-factor(clr.sulf.trait.table$SampDate, levels=c("August 2021","December 2021","April 2022"))

head(clr.sulf.trait.table)

# Note: not looking at every S cycling gene included in this project but looking at ones that appear to have noticeable trends in heat maps
### SOX genes

# `soxA; L-cysteine S-thiosulfotransferase [EC:2.8.5.2]`
soxa.scat<-ggplot(clr.sulf.trait.table, aes(x=PlotID, y=`soxA; L-cysteine S-thiosulfotransferase [EC:2.8.5.2]`,color=SampDate,group=SampDate)) + geom_point(size=4) + geom_line() + theme_bw()+
  labs(title="SoxA Depth of Coverage in Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table$SampDate_Color[order(clr.sulf.trait.table$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("SampleID") + ylab("CLR-Transformed Coverage")

ggsave(soxa.scat,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Fxn_Scatterplots/SoxA_CLR_Coverage_SampleID_scatterplot.png", width=12, height=10, dpi=600)

# `soxB; S-sulfosulfanyl-L-cysteine sulfohydrolase [EC:3.1.6.20]`
soxb.scat<-ggplot(clr.sulf.trait.table, aes(x=PlotID, y=`soxB; S-sulfosulfanyl-L-cysteine sulfohydrolase [EC:3.1.6.20]`,color=SampDate,group=SampDate)) + geom_point(size=4) + geom_line() + theme_bw()+
  labs(title="SoxB Depth of Coverage in Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table$SampDate_Color[order(clr.sulf.trait.table$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("SampleID") + ylab("CLR-Transformed Coverage")

ggsave(soxb.scat,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Fxn_Scatterplots/SoxB_CLR_Coverage_SampleID_scatterplot.png", width=12, height=10, dpi=600)

# `soxX; L-cystein S-thiosulfotransferase`
soxx.scat<-ggplot(clr.sulf.trait.table, aes(x=PlotID, y=`soxX; L-cysteine S-thiosulfotransferase [EC:2.8.5.2]`,color=SampDate,group=SampDate)) + geom_point(size=4) + geom_line() + theme_bw()+
  labs(title="SoxX Depth of Coverage in Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table$SampDate_Color[order(clr.sulf.trait.table$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("SampleID") + ylab("CLR-Transformed Coverage")

ggsave(soxx.scat,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Fxn_Scatterplots/SoxX_CLR_Coverage_SampleID_scatterplot.png", width=12, height=10, dpi=600)

# `soxY; sulfur-oxidizing protein SoxY`
soxy.scat<-ggplot(clr.sulf.trait.table, aes(x=PlotID, y=`soxY; sulfur-oxidizing protein SoxY`,color=SampDate,group=SampDate)) + geom_point(size=4) + geom_line() + theme_bw()+
  labs(title="SoxY Depth of Coverage in Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table$SampDate_Color[order(clr.sulf.trait.table$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("SampleID") + ylab("CLR-Transformed Coverage")

ggsave(soxy.scat,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Fxn_Scatterplots/SoxY_CLR_Coverage_SampleID_scatterplot.png", width=12, height=10, dpi=600)

#`soxZ; sulfur-oxidizing protein SoxZ`
soxz.scat<-ggplot(clr.sulf.trait.table, aes(x=PlotID, y=`soxZ; sulfur-oxidizing protein SoxZ`,color=SampDate,group=SampDate)) + geom_point(size=4) + geom_line() + theme_bw()+
  labs(title="SoxZ Depth of Coverage in Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table$SampDate_Color[order(clr.sulf.trait.table$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("SampleID") + ylab("CLR-Transformed Coverage")

ggsave(soxz.scat,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Fxn_Scatterplots/SoxZ_CLR_Coverage_SampleID_scatterplot.png", width=12, height=10, dpi=600)

#`soxD; S-disulfanyl-L-cysteine oxidoreductase SoxD [EC:1.8.2.6]`
soxd.scat<-ggplot(clr.sulf.trait.table, aes(x=PlotID, y=`soxD; S-disulfanyl-L-cysteine oxidoreductase SoxD [EC:1.8.2.6]`,color=SampDate,group=SampDate)) + geom_point(size=4) + geom_line() + theme_bw()+
  labs(title="SoxD Depth of Coverage in Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table$SampDate_Color[order(clr.sulf.trait.table$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("SampleID") + ylab("CLR-Transformed Coverage")

ggsave(soxd.scat,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Fxn_Scatterplots/SoxD_CLR_Coverage_SampleID_scatterplot.png", width=12, height=10, dpi=600)

#`soxC; sulfane dehydrogenase subunit SoxC`
soxc.scat<-ggplot(clr.sulf.trait.table, aes(x=PlotID, y=`soxC; sulfane dehydrogenase subunit SoxC`,color=SampDate,group=SampDate)) + geom_point(size=4) + geom_line() + theme_bw()+
  labs(title="SoxC Depth of Coverage in Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table$SampDate_Color[order(clr.sulf.trait.table$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("SampleID") + ylab("CLR-Transformed Coverage")

ggsave(soxc.scat,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Fxn_Scatterplots/SoxC_CLR_Coverage_SampleID_scatterplot.png", width=12, height=10, dpi=600)

#### SO3 Oxidation genes
# `soeA; sulfite dehydrogenase (quinone) subunit SoeA [EC:1.8.5.6]`
soeA.scat<-ggplot(clr.sulf.trait.table, aes(x=PlotID, y=`soeA; sulfite dehydrogenase (quinone) subunit SoeA [EC:1.8.5.6]`,color=SampDate,group=SampDate)) + geom_point(size=4) + geom_line() + theme_bw()+
  labs(title="soeA Depth of Coverage in Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table$SampDate_Color[order(clr.sulf.trait.table$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("SampleID") + ylab("CLR-Transformed Coverage")

ggsave(soeA.scat,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Fxn_Scatterplots/soeA_CLR_Coverage_SampleID_scatterplot.png", width=12, height=10, dpi=600)

#### H2S --> S oxidation genes
# `sqr; sulfide:quinone oxidoreductase [EC:1.8.5.4]`
sqr.scat<-ggplot(clr.sulf.trait.table, aes(x=PlotID, y=`sqr; sulfide:quinone oxidoreductase [EC:1.8.5.4]`,color=SampDate,group=SampDate)) + geom_point(size=4) + geom_line() + theme_bw()+
  labs(title="sqr Depth of Coverage in Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table$SampDate_Color[order(clr.sulf.trait.table$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("SampleID") + ylab("CLR-Transformed Coverage")

ggsave(sqr.scat,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Fxn_Scatterplots/sqr_CLR_Coverage_SampleID_scatterplot.png", width=12, height=10, dpi=600)

# `fccB; sulfide dehydrogenase [flavocytochrome c] flavoprotein chain [EC:1.8.2.3]`
fccB.scat<-ggplot(clr.sulf.trait.table, aes(x=PlotID, y=`fccB; sulfide dehydrogenase [flavocytochrome c] flavoprotein chain [EC:1.8.2.3]`,color=SampDate,group=SampDate)) + geom_point(size=4) + geom_line() + theme_bw()+
  labs(title="fccB Depth of Coverage in Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table$SampDate_Color[order(clr.sulf.trait.table$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("SampleID") + ylab("CLR-Transformed Coverage")

ggsave(fccB.scat,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Fxn_Scatterplots/fccB_CLR_Coverage_SampleID_scatterplot.png", width=12, height=10, dpi=600)

#### Dissimilatory SO4 RedOx genes
# `aprA; adenylylsulfate reductase, subunit A [EC:1.8.99.2]`
aprA.scat<-ggplot(clr.sulf.trait.table, aes(x=PlotID, y=`aprA; adenylylsulfate reductase, subunit A [EC:1.8.99.2]`,color=SampDate,group=SampDate)) + geom_point(size=4) + geom_line() + theme_bw()+
  labs(title="aprA Depth of Coverage in Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table$SampDate_Color[order(clr.sulf.trait.table$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("SampleID") + ylab("CLR-Transformed Coverage")

ggsave(aprA.scat,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Fxn_Scatterplots/aprA_CLR_Coverage_SampleID_scatterplot.png", width=12, height=10, dpi=600)

# `dsrB; dissimilatory sulfite reductase beta subunit [EC:1.8.99.5]`
dsrB.scat<-ggplot(clr.sulf.trait.table, aes(x=PlotID, y=`dsrB; dissimilatory sulfite reductase beta subunit [EC:1.8.99.5]`,color=SampDate,group=SampDate)) + geom_point(size=4) + geom_line() + theme_bw()+
  labs(title="dsrB Depth of Coverage in Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table$SampDate_Color[order(clr.sulf.trait.table$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("SampleID") + ylab("CLR-Transformed Coverage")

ggsave(dsrB.scat,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Fxn_Scatterplots/dsrB_CLR_Coverage_SampleID_scatterplot.png", width=12, height=10, dpi=600)

# `dsrA; dissimilatory sulfite reductase alpha subunit [EC:1.8.99.5]`
dsrA.scat<-ggplot(clr.sulf.trait.table, aes(x=PlotID, y=`dsrA; dissimilatory sulfite reductase alpha subunit [EC:1.8.99.5]`,color=SampDate,group=SampDate)) + geom_point(size=4) + geom_line() + theme_bw()+
  labs(title="dsrA Depth of Coverage in Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table$SampDate_Color[order(clr.sulf.trait.table$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("SampleID") + ylab("CLR-Transformed Coverage")

ggsave(dsrA.scat,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Fxn_Scatterplots/dsrA_CLR_Coverage_SampleID_scatterplot.png", width=12, height=10, dpi=600)

# `aprB; adenylylsulfate reductase, subunit B [EC:1.8.99.2]`
aprB.scat<-ggplot(clr.sulf.trait.table, aes(x=PlotID, y=`aprB; adenylylsulfate reductase, subunit B [EC:1.8.99.2]`,color=SampDate,group=SampDate)) + geom_point(size=4) + geom_line() + theme_bw()+
  labs(title="aprB Depth of Coverage in Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table$SampDate_Color[order(clr.sulf.trait.table$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("SampleID") + ylab("CLR-Transformed Coverage")

ggsave(aprB.scat,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Fxn_Scatterplots/aprB_CLR_Coverage_SampleID_scatterplot.png", width=12, height=10, dpi=600)

#### Assimilatory SO4 Reduction

# `sir; sulfite reductase (ferredoxin) [EC:1.8.7.1]`
sir.scat<-ggplot(clr.sulf.trait.table, aes(x=PlotID, y=`sir; sulfite reductase (ferredoxin) [EC:1.8.7.1]`,color=SampDate,group=SampDate)) + geom_point(size=4) + geom_line() + theme_bw()+
  labs(title="sir Depth of Coverage in Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table$SampDate_Color[order(clr.sulf.trait.table$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("SampleID") + ylab("CLR-Transformed Coverage")

ggsave(sir.scat,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Fxn_Scatterplots/sir_CLR_Coverage_SampleID_scatterplot.png", width=12, height=10, dpi=600)

# `cysNC; bifunctional enzyme CysN/CysC [EC:2.7.7.4 2.7.1.25]`
cysNC.scat<-ggplot(clr.sulf.trait.table, aes(x=PlotID, y=`cysNC; bifunctional enzyme CysN/CysC [EC:2.7.7.4 2.7.1.25]`,color=SampDate,group=SampDate)) + geom_point(size=4) + geom_line() + theme_bw()+
  labs(title="CysNC Depth of Coverage in Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table$SampDate_Color[order(clr.sulf.trait.table$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("SampleID") + ylab("CLR-Transformed Coverage")

ggsave(cysNC.scat,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Fxn_Scatterplots/CysNC_CLR_Coverage_SampleID_scatterplot.png", width=12, height=10, dpi=600)

# `cysN; sulfate adenylyltransferase subunit 1 [EC:2.7.7.4]`
# CysN.scat<-ggplot(clr.sulf.trait.table, aes(x=PlotID, y=`cysN; sulfate adenylyltransferase subunit 1 [EC:2.7.7.4]`,color=SampDate,group=SampDate)) + geom_point(size=4) + geom_line() + theme_bw()+
#   labs(title="CysN Depth of Coverage in Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
#   guides(shape = guide_legend(override.aes = list(size = 5)))+
#   scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table$SampDate_Color[order(clr.sulf.trait.table$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
#   xlab("SampleID") + ylab("CLR-Transformed Coverage")
#
# ggsave(CysN.scat,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Fxn_Scatterplots/CysN_CLR_Coverage_SampleID_scatterplot.png", width=12, height=10, dpi=600)

# `cysH; phosphoadenosine phosphosulfate reductase [EC:1.8.4.8 1.8.4.10]`
cysH.scat<-ggplot(clr.sulf.trait.table, aes(x=PlotID, y=`cysH; phosphoadenosine phosphosulfate reductase [EC:1.8.4.8 1.8.4.10]`,color=SampDate,group=SampDate)) + geom_point(size=4) + geom_line() + theme_bw()+
  labs(title="CysH Depth of Coverage in Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table$SampDate_Color[order(clr.sulf.trait.table$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("SampleID") + ylab("CLR-Transformed Coverage")

ggsave(cysH.scat,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Fxn_Scatterplots/CysH_CLR_Coverage_SampleID_scatterplot.png", width=12, height=10, dpi=600)

# `cysC; adenylylsulfate kinase [EC:2.7.1.25]`
# cysC.scat<-ggplot(clr.sulf.trait.table, aes(x=PlotID, y=`cysC; adenylylsulfate kinase [EC:2.7.1.25]`,color=SampDate,group=SampDate)) + geom_point(size=4) + geom_line() + theme_bw()+
#   labs(title="CysC Depth of Coverage in Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
#   guides(shape = guide_legend(override.aes = list(size = 5)))+
#   scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table$SampDate_Color[order(clr.sulf.trait.table$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
#   xlab("SampleID") + ylab("CLR-Transformed Coverage")
#
# ggsave(cysC.scat,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Fxn_Scatterplots/CysC_CLR_Coverage_SampleID_scatterplot.png", width=12, height=10, dpi=600)

#`cysI; sulfite reductase (NADPH) hemoprotein beta-component [EC:1.8.1.2]`
cysI.scat<-ggplot(clr.sulf.trait.table, aes(x=PlotID, y=`cysI; sulfite reductase (NADPH) hemoprotein beta-component [EC:1.8.1.2]`,color=SampDate,group=SampDate)) + geom_point(size=4) + geom_line() + theme_bw()+
  labs(title="CysI Depth of Coverage in Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table$SampDate_Color[order(clr.sulf.trait.table$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("SampleID") + ylab("CLR-Transformed Coverage")

ggsave(cysI.scat,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Fxn_Scatterplots/CysI_CLR_Coverage_SampleID_scatterplot.png", width=12, height=10, dpi=600)

# `cysD; sulfate adenylyltransferase subunit 2 [EC:2.7.7.4]`
cysD.scat<-ggplot(clr.sulf.trait.table, aes(x=PlotID, y=`cysD; sulfate adenylyltransferase subunit 2 [EC:2.7.7.4]`,color=SampDate,group=SampDate)) + geom_point(size=4) + geom_line() + theme_bw()+
  labs(title="CysD Depth of Coverage in Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table$SampDate_Color[order(clr.sulf.trait.table$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("SampleID") + ylab("CLR-Transformed Coverage")

ggsave(cysD.scat,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Fxn_Scatterplots/CysD_CLR_Coverage_SampleID_scatterplot.png", width=12, height=10, dpi=600)

# `sir; sulfite reductase (ferredoxin) [EC:1.8.7.1]`
sir.scat<-ggplot(clr.sulf.trait.table, aes(x=PlotID, y=`sir; sulfite reductase (ferredoxin) [EC:1.8.7.1]`,color=SampDate,group=SampDate)) + geom_point(size=4) + geom_line() + theme_bw()+
  labs(title="sir Depth of Coverage in Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table$SampDate_Color[order(clr.sulf.trait.table$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("SampleID") + ylab("CLR-Transformed Coverage")

ggsave(sir.scat,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Fxn_Scatterplots/sir_CLR_Coverage_SampleID_scatterplot.png", width=12, height=10, dpi=600)

#### S Disproportionation

# `phsA,psrA; thiosulfate reductase / polysulfide reductase chain A [EC:1.8.5.5]`
phsA.psrA.scat<-ggplot(clr.sulf.trait.table, aes(x=PlotID, y=`phsA,psrA; thiosulfate reductase / polysulfide reductase chain A [EC:1.8.5.5]`,color=SampDate,group=SampDate)) + geom_point(size=4) + geom_line() + theme_bw()+
  labs(title="phsA/psrA Depth of Coverage in Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table$SampDate_Color[order(clr.sulf.trait.table$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("SampleID") + ylab("CLR-Transformed Coverage")

ggsave(phsA.psrA.scat,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Fxn_Scatterplots/phsA.psrA_CLR_Coverage_SampleID_scatterplot.png", width=12, height=10, dpi=600)

#### Genes in Multiple S Pathways

# `sat, met3; sulfate adenylyltransferase [EC:2.7.7.4]`
sat.scat<-ggplot(clr.sulf.trait.table, aes(x=PlotID, y=`sat, met3; sulfate adenylyltransferase [EC:2.7.7.4]`,color=SampDate,group=SampDate)) + geom_point(size=4) + geom_line() + theme_bw()+
  labs(title="sat/met3 Depth of Coverage in Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(clr.sulf.trait.table$SampDate_Color[order(clr.sulf.trait.table$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("SampleID") + ylab("CLR-Transformed Coverage")

ggsave(sat.scat,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Fxn_Scatterplots/sat_CLR_Coverage_SampleID_scatterplot.png", width=12, height=10, dpi=600)

#### Pull out Sulfur Metabolic Fxns from Binary Data ####

# mgm_fxn.binary --> ran counts_to_binary() custom function on ko.cov.sum_table to determine if KO appears at all in samples or not

sulf.ko.bi<-mgm_fxn.binary[,which(colnames(mgm_fxn.binary) %in% sulfur.fxns$KO_ID)] # merge CLR data w/ S fxns found in contigs from KOFamScan
sulf.ko.bi$SampleID<-rownames(sulf.ko.bi)
sulf.ko.bi.melt<-melt(sulf.ko.bi, by="SampleID")
colnames(sulf.ko.bi.melt)[which(names(sulf.ko.bi.melt) == "variable")] <- "KO_ID"
colnames(sulf.ko.bi.melt)[which(names(sulf.ko.bi.melt) == "value")] <- "PresAb"
head(sulf.ko.bi.melt) #sanity check

clr.sulf.ko.bi<-merge(sulf.ko.bi.melt,sulf.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
head(clr.sulf.ko.bi)
colnames(clr.sulf.ko.bi)[which(names(clr.sulf.ko.bi) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
clr.cov.sum.sulf.ko.bi<-as.data.frame(dcast(clr.sulf.ko.bi, SampleID~KO_Function.KEGG, value.var="PresAb", fun.aggregate=sum)) ###just dcast, nothing is being added here!
rownames(clr.cov.sum.sulf.ko.bi)<-clr.cov.sum.sulf.ko.bi$SampleID
clr.cov.sum.sulf.ko.bi[1:4,]

# sanity check
clr.cov.sum.sulf.ko.bi$`cysH; phosphoadenosine phosphosulfate reductase [EC:1.8.4.8 1.8.4.10]`[1:4]
head(clr.sulf.ko.bi)

#### Sulfur Binary Heat Maps ####

# prep for ggplot2 heatmap
clr.sulf.ko.bi[1:4,]

clr.sulf.all.bi<-merge(clr.sulf.ko.bi,meta_scaled,by="SampleID")

head(clr.sulf.ko.bi)
clr.sulf.all.bi$PlotID = factor(clr.sulf.all.bi$PlotID, levels=unique(clr.sulf.all.bi$PlotID[order(clr.sulf.all.bi$SampDate,clr.sulf.all.bi$Depth_m)]), ordered=TRUE)
clr.sulf.all.bi$SampDate<-gsub("\\."," ",clr.sulf.all.bi$SampDate)
clr.sulf.all.bi$SampDate<-factor(clr.sulf.all.bi$SampDate, levels=c("August 2021","December 2021","April 2022"))
clr.sulf.all.bi$KO_Function.KEGG = factor(clr.sulf.all.bi$KO_Function.KEGG, levels=unique(clr.sulf.all.bi$KO_Function.KEGG[order(clr.sulf.all.bi$Pathway)]), ordered=TRUE)

clr.sulf.all.bi$PathShort<-clr.sulf.all.bi$Pathway
clr.sulf.all.bi$PathShort[(clr.sulf.all.bi$PathShort) == "Dissimilatory Sulfate Redox"] <- "D.SO4 RedOx"
clr.sulf.all.bi$PathShort[(clr.sulf.all.bi$PathShort) == "Assimilatory Sulfate Reduction"] <- "A.SO4 Red"
clr.sulf.all.bi$PathShort[(clr.sulf.all.bi$PathShort) == "Multiple Pathways"] <- "Multi Paths"
clr.sulf.all.bi$PathShort[(clr.sulf.all.bi$PathShort) == "S Disproportionation"] <- "S Disprop."

clr.sulf.all.bi$Pathway<-factor(clr.sulf.all.bi$Pathway,levels=c("Assimilatory Sulfate Reduction","Dissimilatory Sulfate Redox","Multiple Pathways","SOX","S Disproportionation"))
clr.sulf.all.bi$PathShort<-factor(clr.sulf.all.bi$PathShort,levels=c("A.SO4 Red","D.SO4 RedOx","Multi Paths","SOX","S Disprop."))

clr.sulf.all.bi$PathSpecShort<-clr.sulf.all.bi$PathwaySpecific
clr.sulf.all.bi$PathSpecShort[(clr.sulf.all.bi$PathSpecShort) == "Dissimilatory Sulfate Redox"] <- "D.SO4 RedOx"
clr.sulf.all.bi$PathSpecShort[(clr.sulf.all.bi$PathSpecShort) == "Assimilatory Sulfate Reduction"] <- "A.SO4 Red"
clr.sulf.all.bi$PathSpecShort[(clr.sulf.all.bi$PathSpecShort) == "Multiple Pathways"] <- "Multi"
clr.sulf.all.bi$PathSpecShort[(clr.sulf.all.bi$PathSpecShort) == "Sulfur Disproportionation"] <- "S Disprop."
clr.sulf.all.bi$PathSpecShort[(clr.sulf.all.bi$PathSpecShort) == "Sulfide Oxidation"] <- "H2S Ox"
clr.sulf.all.bi$PathSpecShort[(clr.sulf.all.bi$PathSpecShort) == "Sulfite Oxidation"] <- "SO3 Ox"
clr.sulf.all.bi$PathSpecShort[(clr.sulf.all.bi$PathSpecShort) == "Thiosulfate Oxidation"] <- "S2O3 Ox"

clr.sulf.all.bi$PathwaySpecific<-factor(clr.sulf.all.bi$PathwaySpecific,levels=c("Assimilatory Sulfate Reduction","Dissimilatory Sulfate Redox","Multiple Pathways","Sulfur Disproportionation","Sulfide Oxidation","Sulfite Oxidation","Thiosulfate Oxidation"))
clr.sulf.all.bi$PathSpecShort<-factor(clr.sulf.all.bi$PathSpecShort,levels=c("A.SO4 Red","D.SO4 RedOx","Multi","S Disprop.","H2S Ox","SO3 Ox","S2O3 Ox"))

head(clr.sulf.all.bi)

# Figures
sulf.bi.hm1a<-ggplot(clr.sulf.all.bi, aes(PlotID, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Sulfur Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(sulf.bi.hm1a,filename = "figures/MGM_Figs/FxnDiv/Sulfur/PresenceAbsence/Sulfur_KOFxns_MGMs_by_Function_Binary_heatmap.png", width=18, height=13, dpi=600)

# sulf.bi.hm1b<-ggplot(clr.sulf.all.bi, aes(PlotID, KO_Function.KEGG, fill=factor(PresAb))) +
#   geom_tile(colour="black",size=0.25) +
#   scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Sulfur Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~.,scales="free_y", space = "free")
#
# ggsave(sulf.bi.hm1b,filename = "figures/MGM_Figs/FxnDiv/Sulfur/PresenceAbsence/Sulfur_KOFxns_MGMs_by_Function_Pathway_Binary_heatmap.png", width=17, height=15, dpi=600)

sulf.bi.hm1b2<-ggplot(clr.sulf.all.bi, aes(Depth_m, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Sulfur Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~.,scales="free_y", space = "free")

ggsave(sulf.bi.hm1b2,filename = "figures/MGM_Figs/FxnDiv/Sulfur/PresenceAbsence/Sulfur_KOFxns_MGMs_Depth_by_Function_Pathway_Binary_heatmap.png", width=17, height=15, dpi=600)

sulf.bi.hm1c<-ggplot(clr.sulf.all.bi, aes(PlotID, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Sulfur Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate,scales="free_x", space = "free")

ggsave(sulf.bi.hm1c,filename = "figures/MGM_Figs/FxnDiv/Sulfur/PresenceAbsence/Sulfur_KOFxns_MGMs_by_Function_SampDate_Binary_best_heatmap.png", width=20, height=13, dpi=600)

sulf.bi.hm1c2<-ggplot(clr.sulf.all.bi, aes(Depth_m, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Sulfur Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate,scales="free_x", space = "free")

ggsave(sulf.bi.hm1c2,filename = "figures/MGM_Figs/FxnDiv/Sulfur/PresenceAbsence/Sulfur_KOFxns_MGMs_Depth_by_Function_SampDate_Binary_best_heatmap.png", width=20, height=13, dpi=600)

sulf.bi.hm1d<-ggplot(clr.sulf.all.bi, aes(PlotID, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Sulfur Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~SampDate, scales="free", space = "free")

ggsave(sulf.bi.hm1d,filename = "figures/MGM_Figs/FxnDiv/Sulfur/PresenceAbsence/Sulfur_KOFxns_MGMs_by_Function_SampDate_Pathway_Binary_best_heatmap.png", width=20, height=15, dpi=600)

sulf.bi.hm1d2<-ggplot(clr.sulf.all.bi, aes(Depth_m, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Sulfur Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~SampDate, scales="free", space = "free")

ggsave(sulf.bi.hm1d2,filename = "figures/MGM_Figs/FxnDiv/Sulfur/PresenceAbsence/Sulfur_KOFxns_MGMs_Depth_by_Function_SampDate_Pathway_Binary_best_heatmap.png", width=20, height=15, dpi=600)

# sulf.bi.hm1c2<-ggplot(clr.sulf.all.bi, aes(PlotID, KO_Function.KEGG, fill=factor(PresAb))) +
#   geom_tile(colour="black",size=0.25) +
#   scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Sulfur Metabolism in Salton Seawater Metagenomes",fill="Presence/Absence") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathwaySpecific~SampDate, scales="free", space = "free")
#
# ggsave(sulf.bi.hm1c2,filename = "figures/MGM_Figs/FxnDiv/Sulfur/PresenceAbsence/Sulfur_KOFxns_MGMs_by_Function_SampDate_PathwaySpecific_Binary_best_heatmap.png", width=20, height=20, dpi=600)

sulf.bi.hm1e<-ggplot(clr.sulf.all.bi, aes(Depth_m, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Sulfur Metabolism in Salton Seawater Metagenomes",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathSpecShort~SampDate, scales="free", space = "free")

ggsave(sulf.bi.hm1e,filename = "figures/MGM_Figs/FxnDiv/Sulfur/PresenceAbsence/Sulfur_KOFxns_MGMs_Depth_by_Function_SampDate_PathwaySpecific_Binary_best_heatmap.png", width=20, height=20, dpi=600)

# sulf.bi.hm1e0<-ggplot(clr.sulf.all.bi[clr.sulf.all.bi$Depth_m==0,], aes(PathShort, KO_Function.KEGG, fill=factor(PresAb))) +
#   geom_tile(colour="black",size=0.25) +
#   scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Sulfur Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate)
#
# ggsave(sulf.bi.hm1e0,filename = "figures/MGM_Figs/FxnDiv/Sulfur/PresenceAbsence/Sulfur_KOFxns_Pathways_Binary_0m_heatmap.png", width=18, height=18, dpi=600)
#
# sulf.bi.hm1e5<-ggplot(clr.sulf.all.bi[clr.sulf.all.bi$Depth_m==5,], aes(PathShort, KO_Function.KEGG, fill=factor(PresAb))) +
#   geom_tile(colour="black",size=0.25) +
#   scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Sulfur Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)
#
# ggsave(sulf.bi.hm1e5,filename = "figures/MGM_Figs/FxnDiv/Sulfur/PresenceAbsence/Sulfur_KOFxns_Pathways_Binary_5m_heatmap.png", width=18, height=18, dpi=600)
#
# sulf.bi.hm1e6<-ggplot(clr.sulf.all.bi[clr.sulf.all.bi$Depth_m==10,], aes(PathShort, KO_Function.KEGG, fill=factor(PresAb))) +
#   geom_tile(colour="black",size=0.25) +
#   scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Sulfur Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)
#
# ggsave(sulf.bi.hm1e6,filename = "figures/MGM_Figs/FxnDiv/Sulfur/PresenceAbsence/Sulfur_KOFxns_Pathways_Binary_10m_heatmap.png", width=18, height=18, dpi=600)

#### Pull Out Nitrogen Metabolic Fxns from CLR data - NO NAs ####
## heatmaps of traits of interest

mgm.clr[1:4,1:4]

# pull out nitro functions from CLR transformed, summed coverages (summed coverage per KO)
nitro.ko<-mgm.clr[,which(colnames(mgm.clr) %in% nitro.fxns$KO_ID)] # merge CLR data w/ S fxns found in contigs from KOFamScan
nitro.ko$SampleID<-rownames(nitro.ko)
nitro.ko.melt<-melt(nitro.ko, by="SampleID")
colnames(nitro.ko.melt)[which(names(nitro.ko.melt) == "variable")] <- "KO_ID"
colnames(nitro.ko.melt)[which(names(nitro.ko.melt) == "value")] <- "CLR_SumCovPerKO"
head(nitro.ko.melt) #sanity check

clr.nitro.ko<-merge(nitro.ko.melt,nitro.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
head(clr.nitro.ko)
colnames(clr.nitro.ko)[which(names(clr.nitro.ko) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
clr.cov.sum.nitro.ko<-as.data.frame(dcast(clr.nitro.ko, SampleID~KO_Function.KEGG, value.var="CLR_SumCovPerKO", fun.aggregate=sum)) ###just dcast, nothing is being added here!
rownames(clr.cov.sum.nitro.ko)<-clr.cov.sum.nitro.ko$SampleID
clr.cov.sum.nitro.ko[1:4,]

# sanity check
clr.cov.sum.nitro.ko$`nirK; nitrite reductase (NO-forming) [EC:1.7.2.1]`[1:4]
head(clr.nitro.ko)


#### Pull Out Nitrogen Metabolic Fxns from CLR data - with NAs ####
## heatmaps of traits of interest

mgm.clr.na[1:4,1:4]

# pull out nitro functions from CLR transformed, summed coverages (summed coverage per KO)
nitro.ko.na<-mgm.clr.na[,which(colnames(mgm.clr.na) %in% nitro.fxns$KO_ID)] # merge CLR data w/ S fxns found in contigs from KOFamScan
nitro.ko.na$SampleID<-rownames(nitro.ko.na)
nitro.ko.na.melt<-melt(nitro.ko.na, by="SampleID")
colnames(nitro.ko.na.melt)[which(names(nitro.ko.na.melt) == "variable")] <- "KO_ID"
colnames(nitro.ko.na.melt)[which(names(nitro.ko.na.melt) == "value")] <- "CLR_SumCovPerKO"
head(nitro.ko.na.melt) #sanity check

clr.nitro.ko.na<-merge(nitro.ko.na.melt,nitro.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
head(clr.nitro.ko.na)
colnames(clr.nitro.ko.na)[which(names(clr.nitro.ko.na) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
clr.cov.sum.nitro.ko.na<-as.data.frame(dcast(clr.nitro.ko.na, SampleID~KO_Function.KEGG, value.var="CLR_SumCovPerKO", fun.aggregate=sum)) ###just dcast, nothing is being added here!
rownames(clr.cov.sum.nitro.ko.na)<-clr.cov.sum.nitro.ko.na$SampleID
clr.cov.sum.nitro.ko.na[1:4,]

# sanity check
clr.cov.sum.nitro.ko.na$`nirK; nitrite reductase (NO-forming) [EC:1.7.2.1]`[1:4]
head(clr.nitro.ko.na)

### Nitrogen Heat Maps ####
# see max & mean of summed
max(clr.cov.sum.nitro.ko[,-1])
mean(as.matrix(clr.cov.sum.nitro.ko[,-1]))

# first heat map of nitro KOs
heatmap(as.matrix(clr.cov.sum.nitro.ko[,-1]), scale = "none")

colSums(clr.cov.sum.nitro.ko[,-1])
#clr.cov.sum.nitro.ko2 <- clr.cov.sum.nitro.ko[,which(colSums(clr.cov.sum.nitro.ko[,-1])>10)]

heatmap(as.matrix(clr.cov.sum.nitro.ko[,-1]), scale = "none")

# prep for ggplot2 heatmap
clr.nitro.ko.na[1:4,]
clr.nitro.all<-merge(clr.nitro.ko.na,meta_scaled,by="SampleID")
head(clr.nitro.all)
clr.nitro.all$PlotID = factor(clr.nitro.all$PlotID, levels=unique(clr.nitro.all$PlotiD[order(clr.nitro.all$SampDate,clr.nitro.all$Depth_m)]), ordered=TRUE)
clr.nitro.all$SampDate<-gsub("\\."," ",clr.nitro.all$SampDate)
clr.nitro.all$SampDate<-factor(clr.nitro.all$SampDate, levels=c("August 2021","December 2021","April 2022"))
clr.nitro.all$KO_Function.KEGG = factor(clr.nitro.all$KO_Function.KEGG, levels=unique(clr.nitro.all$KO_Function.KEGG[order(clr.nitro.all$Pathway)]), ordered=TRUE)

# create shortened name for pathways
clr.nitro.all$PathShort<-clr.nitro.all$Pathway
# vvv can only do this type of renaming if variables are characters, not factors
clr.nitro.all$PathShort[(clr.nitro.all$PathShort) == "Dissimilatory Nitrate Reduction"] <- "D. NO3 Red"
clr.nitro.all$PathShort[(clr.nitro.all$PathShort) == "Assimilatory Nitrate Reduction"] <- "A. NO3 Red"

# turn pathways & pathshort into factors
clr.nitro.all$Pathway<-factor(clr.nitro.all$Pathway,levels=c("Assimilatory Nitrate Reduction","Dissimilatory Nitrate Reduction","Multiple Pathways","Denitrification","Anammox"))
clr.nitro.all$PathShort<-factor(clr.nitro.all$PathShort,levels=c("A. NO3 Red","D. NO3 Red","Multiple Pathways","Denitrification","Anammox"))

head(clr.nitro.all)

# For heatmap color gradient
max(clr.nitro.all$CLR_SumCovPerKO, na.rm=TRUE)
max(clr.nitro.all$CLR_SumCovPerKO, na.rm=TRUE)/2
min(clr.nitro.all$CLR_SumCovPerKO, na.rm=TRUE)

# Figures below
# by SampleID

nitro.hm1a<-ggplot(clr.nitro.all, aes(PlotID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.4","0.5","-0.4"),breaks=c(1.4,0.5,-0.4)) + labs(title="Nitrogen Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(nitro.hm1a,filename = "figures/MGM_Figs/FxnDiv/Nitrogen/Nitrogen_KOFxns_MGMs_SampID_by_Function_heatmap.png", width=18, height=13, dpi=600)

nitro.hm1a2<-ggplot(clr.nitro.all, aes(PlotID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.4","0.5","-0.4"),breaks=c(1.4,0.5,-0.4)) + labs(title="Nitrogen Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Pathway~.,scales="free_y", space = "free")

ggsave(nitro.hm1a2,filename = "figures/MGM_Figs/FxnDiv/Nitrogen/Nitrogen_KOFxns_MGMs_SampID_by_Function_Pathway_heatmap.png", width=17, height=15, dpi=600)

nitro.hm1a3<-ggplot(clr.nitro.all, aes(PlotID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.4","0.5","-0.4"),breaks=c(1.4,0.5,-0.4)) + labs(title="Nitrogen Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate,scales="free_x", space = "free")

ggsave(nitro.hm1a3,filename = "figures/MGM_Figs/FxnDiv/Nitrogen/Nitrogen_KOFxns_MGMs_SampleID_by_Function_SampDate_best_heatmap.png", width=20, height=13, dpi=600)

nitro.hm1a4<-ggplot(clr.nitro.all, aes(PlotID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.4","0.5","-0.4"),breaks=c(1.4,0.5,-0.4)) + labs(title="Nitrogen Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~.,scales="free_y", space = "free")

ggsave(nitro.hm1a4,filename = "figures/MGM_Figs/FxnDiv/Nitrogen/Nitrogen_KOFxns_MGMs_SampID_by_Function_Pathway_heatmap2.png", width=17, height=15, dpi=600)

nitro.hm1a5<-ggplot(clr.nitro.all, aes(PlotID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.4","0.5","-0.4"),breaks=c(1.4,0.5,-0.4)) + labs(title="Nitrogen Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Pathway~SampDate, scales="free", space = "free")

ggsave(nitro.hm1a5,filename = "figures/MGM_Figs/FxnDiv/Nitrogen/Nitrogen_KOFxns_MGMs_SampleID_by_Function_SampDate_Pathway_best_heatmap.png", width=20, height=15, dpi=600)

nitro.hm1a6<-ggplot(clr.nitro.all, aes(PlotID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.4","0.5","-0.4"),breaks=c(1.4,0.5,-0.4)) + labs(title="Nitrogen Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~SampDate, scales="free", space = "free")

ggsave(nitro.hm1a6,filename = "figures/MGM_Figs/FxnDiv/Nitrogen/Nitrogen_KOFxns_MGMs_SampleID_by_Function_SampDate_Pathway_best_heatmap2.png", width=20, height=15, dpi=600)

# by Depth
nitro.hm1b<-ggplot(clr.nitro.all, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.4","0.5","-0.4"),breaks=c(1.4,0.5,-0.4)) + labs(title="Nitrogen Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(nitro.hm1b,filename = "figures/MGM_Figs/FxnDiv/Nitrogen/Nitrogen_KOFxns_MGMs_Depth_by_Function_heatmap.png", width=18, height=13, dpi=600)

nitro.hm1b2<-ggplot(clr.nitro.all, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.4","0.5","-0.4"),breaks=c(1.4,0.5,-0.4)) + labs(title="Nitrogen Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Pathway~.,scales="free_y", space = "free")

ggsave(nitro.hm1b2,filename = "figures/MGM_Figs/FxnDiv/Nitrogen/Nitrogen_KOFxns_MGMs_Depth_by_Function_Pathway_heatmap.png", width=17, height=15, dpi=600)

nitro.hm1b3<-ggplot(clr.nitro.all, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.4","0.5","-0.4"),breaks=c(1.4,0.5,-0.4)) + labs(title="Nitrogen Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate,scales="free_x", space = "free")

ggsave(nitro.hm1b3,filename = "figures/MGM_Figs/FxnDiv/Nitrogen/Nitrogen_KOFxns_MGMs_Depth_by_Function_SampDate_best_heatmap.png", width=20, height=13, dpi=600)

nitro.hm1b4<-ggplot(clr.nitro.all, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.4","0.5","-0.4"),breaks=c(1.4,0.5,-0.4)) + labs(title="Nitrogen Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~.,scales="free_y", space = "free")

ggsave(nitro.hm1b4,filename = "figures/MGM_Figs/FxnDiv/Nitrogen/Nitrogen_KOFxns_MGMs_Depth_by_Function_Pathway_heatmap2.png", width=17, height=15, dpi=600)

nitro.hm1b5<-ggplot(clr.nitro.all, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.4","0.5","-0.4"),breaks=c(1.4,0.5,-0.4)) + labs(title="Nitrogen Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Pathway~SampDate, scales="free", space = "free")

ggsave(nitro.hm1b5,filename = "figures/MGM_Figs/FxnDiv/Nitrogen/Nitrogen_KOFxns_MGMs_Depth_by_Function_SampDate_Pathway_best_heatmap.png", width=20, height=15, dpi=600)

nitro.hm1b6<-ggplot(clr.nitro.all, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.4","0.5","-0.4"),breaks=c(1.4,0.5,-0.4)) + labs(title="Nitrogen Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~SampDate, scales="free", space = "free")

ggsave(nitro.hm1b6,filename = "figures/MGM_Figs/FxnDiv/Nitrogen/Nitrogen_KOFxns_MGMs_Depth_by_Function_SampDate_Pathway_best_heatmap2.png", width=20, height=15, dpi=600)

# nitro.hm1e<-ggplot(clr.nitro.all[clr.nitro.all$Depth_m==0,], aes(Pathway, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.4","0.5","-0.4"),breaks=c(1.4,0.5,-0.4)) + labs(title="Nitrogen Metabolism in Salton Seawater Metagenomes - 0m",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11,face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)
#
# ggsave(nitro.hm1e,filename = "figures/MGM_Figs/FxnDiv/Nitrogen/Nitrogen_KOFxns_Pathways_MGMs_0m_heatmap.png", width=18, height=18, dpi=600)
#
# nitro.hm1f<-ggplot(clr.nitro.all[clr.nitro.all$Depth_m==5,], aes(Pathway, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.4","0.5","-0.4"),breaks=c(1.4,0.5,-0.4)) + labs(title="Nitrogen Metabolism in Salton Seawater Metagenomes - 5m",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)
#
# ggsave(nitro.hm1f,filename = "figures/MGM_Figs/FxnDiv/Nitrogen/Nitrogen_KOFxns_Pathways_MGMs_5m_heatmap.png", width=18, height=18, dpi=600)
#
# nitro.hm1g<-ggplot(clr.nitro.all[clr.nitro.all$Depth_m==10,], aes(Pathway, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.4","0.5","-0.4"),breaks=c(1.4,0.5,-0.4)) + labs(title="Nitrogen Metabolism in Salton Seawater Metagenomes - 10m",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)
#
# ggsave(nitro.hm1g,filename = "figures/MGM_Figs/FxnDiv/Nitrogen/Nitrogen_KOFxns_Pathways_MGMs_10m_heatmap.png", width=18, height=18, dpi=600)

#### Pull out Nitrogen Metabolic Fxns from Binary Data ####

nitro.ko.bi<-mgm_fxn.binary[,which(colnames(mgm_fxn.binary) %in% nitro.fxns$KO_ID)] # merge CLR data w/ N fxns found in contigs from KOFamScan
nitro.ko.bi$SampleID<-rownames(nitro.ko.bi)
nitro.ko.bi.melt<-melt(nitro.ko.bi, by="SampleID")
colnames(nitro.ko.bi.melt)[which(names(nitro.ko.bi.melt) == "variable")] <- "KO_ID"
colnames(nitro.ko.bi.melt)[which(names(nitro.ko.bi.melt) == "value")] <- "PresAb"
head(nitro.ko.bi.melt) #sanity check

clr.nitro.ko.bi<-merge(nitro.ko.bi.melt,nitro.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
head(clr.nitro.ko.bi)
colnames(clr.nitro.ko.bi)[which(names(clr.nitro.ko.bi) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
clr.cov.sum.nitro.ko.bi<-as.data.frame(dcast(clr.nitro.ko.bi, SampleID~KO_Function.KEGG, value.var="PresAb", fun.aggregate=sum)) ###just dcast, nothing is being added here!
rownames(clr.cov.sum.nitro.ko.bi)<-clr.cov.sum.nitro.ko.bi$SampleID
clr.cov.sum.nitro.ko.bi[1:4,]

# sanity check
clr.cov.sum.nitro.ko.bi$`nosZ; nitrous-oxide reductase [EC:1.7.2.4]`[1:4]
head(clr.nitro.ko.bi)

#### Nitrogen Binary Heat Maps ####
# see max & mean of summed
max(clr.cov.sum.nitro.ko[,-1])
mean(as.matrix(clr.cov.sum.nitro.ko[,-1]))

# first heat map of nitro KOs
heatmap(as.matrix(clr.cov.sum.nitro.ko[,-1]), scale = "none")

colSums(clr.cov.sum.nitro.ko[,-1])
#clr.cov.sum.nitro.ko2 <- clr.cov.sum.nitro.ko[,which(colSums(clr.cov.sum.nitro.ko[,-1])>10)]

heatmap(as.matrix(clr.cov.sum.nitro.ko[,-1]), scale = "none")

# prep for ggplot2 heatmap
clr.nitro.ko.bi[1:4,]

clr.nitro.all.bi<-merge(clr.nitro.ko.bi,meta_scaled,by="SampleID")

head(clr.sulf.all)
clr.nitro.all.bi$PlotID = factor(clr.nitro.all.bi$PlotID, levels=unique(clr.nitro.all.bi$PlotID[order(clr.nitro.all.bi$SampDate,clr.nitro.all.bi$Depth_m)]), ordered=TRUE)
clr.nitro.all.bi$SampDate<-gsub("\\."," ",clr.nitro.all.bi$SampDate)
clr.nitro.all.bi$SampDate<-factor(clr.nitro.all.bi$SampDate, levels=c("August 2021","December 2021","April 2022"))
clr.nitro.all.bi$KO_Function.KEGG = factor(clr.nitro.all.bi$KO_Function.KEGG, levels=unique(clr.nitro.all.bi$KO_Function.KEGG[order(clr.nitro.all.bi$Pathway)]), ordered=TRUE)

# create shortened name for pathways
clr.nitro.all.bi$PathShort<-clr.nitro.all.bi$Pathway
# vvv can only do this type of renaming if variables are characters, not factors
clr.nitro.all.bi$PathShort[(clr.nitro.all.bi$PathShort) == "Dissimilatory Nitrate Reduction"] <- "D. NO3 Red"
clr.nitro.all.bi$PathShort[(clr.nitro.all.bi$PathShort) == "Assimilatory Nitrate Reduction"] <- "A. NO3 Red"

# turn pathways & pathshort into factors
clr.nitro.all.bi$Pathway<-factor(clr.nitro.all.bi$Pathway,levels=c("Assimilatory Nitrate Reduction","Dissimilatory Nitrate Reduction","Multiple Pathways","Denitrification","Anammox"))
clr.nitro.all.bi$PathShort<-factor(clr.nitro.all.bi$PathShort,levels=c("A. NO3 Red","D. NO3 Red","Multiple Pathways","Denitrification","Anammox"))

head(clr.nitro.all.bi)

# Figures
nitro.bi.hm1a<-ggplot(clr.nitro.all.bi, aes(PlotID, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Nitrogen Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(nitro.bi.hm1a,filename = "figures/MGM_Figs/FxnDiv/Nitrogen/PresenceAbsence/Nitrogen_KOFxns_MGMs_by_Function_Binary_heatmap.png", width=18, height=13, dpi=600)

nitro.bi.hm1b<-ggplot(clr.nitro.all.bi, aes(PlotID, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Nitrogen Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~.,scales="free_y", space = "free")

ggsave(nitro.bi.hm1b,filename = "figures/MGM_Figs/FxnDiv/Nitrogen/PresenceAbsence/Nitrogen_KOFxns_MGMs_by_Function_Pathway_Binary_heatmap.png", width=17, height=15, dpi=600)

nitro.bi.hm1b2<-ggplot(clr.nitro.all.bi, aes(Depth_m, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Nitrogen Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~.,scales="free_y", space = "free")

ggsave(nitro.bi.hm1b2,filename = "figures/MGM_Figs/FxnDiv/Nitrogen/PresenceAbsence/Nitrogen_KOFxns_MGMs_Depth_by_Function_Pathway_Binary_heatmap.png", width=17, height=15, dpi=600)

nitro.bi.hm1d<-ggplot(clr.nitro.all.bi, aes(PlotID, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Nitrogen Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate,scales="free_x", space = "free")

ggsave(nitro.bi.hm1d,filename = "figures/MGM_Figs/FxnDiv/Nitrogen/PresenceAbsence/Nitrogen_KOFxns_MGMs_by_Function_SampDate_Binary_best_heatmap.png", width=20, height=13, dpi=600)

nitro.bi.hm1d2<-ggplot(clr.nitro.all.bi, aes(Depth_m, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Nitrogen Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate,scales="free_x", space = "free")

ggsave(nitro.bi.hm1d2,filename = "figures/MGM_Figs/FxnDiv/Nitrogen/PresenceAbsence/Nitrogen_KOFxns_MGMs_Depth_by_Function_SampDate_Binary_best_heatmap.png", width=20, height=13, dpi=600)

nitro.bi.hm1e<-ggplot(clr.nitro.all.bi, aes(PlotID, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Nitrogen Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~SampDate, scales="free", space = "free")

ggsave(nitro.bi.hm1e,filename = "figures/MGM_Figs/FxnDiv/Nitrogen/PresenceAbsence/Nitrogen_KOFxns_MGMs_by_Function_SampDate_Pathway_Binary_best_heatmap.png", width=20, height=15, dpi=600)

nitro.bi.hm1e2<-ggplot(clr.nitro.all.bi, aes(Depth_m, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Nitrogen Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~SampDate, scales="free", space = "free")

ggsave(nitro.bi.hm1e,filename = "figures/MGM_Figs/FxnDiv/Nitrogen/PresenceAbsence/Nitrogen_KOFxns_MGMs_Depth_by_Function_SampDate_Pathway_Binary_best_heatmap.png", width=20, height=15, dpi=600)

nitro.bi.hm1e2<-ggplot(clr.nitro.all.bi, aes(Depth_m, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Nitrogen Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~SampDate, scales="free", space = "free")

ggsave(nitro.bi.hm1e,filename = "figures/MGM_Figs/FxnDiv/Nitrogen/PresenceAbsence/Nitrogen_KOFxns_MGMs_Depth_by_Function_SampDate_Pathway_Binary_best_heatmap.png", width=20, height=15, dpi=600)

# nitro.bi.hm1e0<-ggplot(clr.nitro.all.bi[clr.nitro.all.bi$Depth_m==0,], aes(PathShort, KO_Function.KEGG, fill=factor(PresAb))) +
#   geom_tile(colour="black",size=0.25) +
#   scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Nitrogen Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate)
#
# ggsave(nitro.bi.hm1e0,filename = "figures/MGM_Figs/FxnDiv/Nitrogen/PresenceAbsence/Nitrogen_KOFxns_Pathways_Binary_0m_heatmap.png", width=18, height=18, dpi=600)
#
# nitro.bi.hm1e5<-ggplot(clr.nitro.all.bi[clr.nitro.all.bi$Depth_m==5,], aes(PathShort, KO_Function.KEGG, fill=factor(PresAb))) +
#   geom_tile(colour="black",size=0.25) +
#   scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Nitrogen Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)
#
# ggsave(nitro.bi.hm1e5,filename = "figures/MGM_Figs/FxnDiv/Nitrogen/PresenceAbsence/Nitrogen_KOFxns_Pathways_Binary_5m_heatmap.png", width=18, height=18, dpi=600)
#
# nitro.bi.hm1e6<-ggplot(clr.nitro.all.bi[clr.nitro.all.bi$Depth_m==10,], aes(PathShort, KO_Function.KEGG, fill=factor(PresAb))) +
#   geom_tile(colour="black",size=0.25) +
#   scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Nitrogen Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)
#
# ggsave(nitro.bi.hm1e6,filename = "figures/MGM_Figs/FxnDiv/Nitrogen/PresenceAbsence/Nitrogen_KOFxns_Pathways_Binary_10m_heatmap.png", width=18, height=18, dpi=600)

#### Pull Out Carbon Metabolic Fxns from CLR data - NO NAs ####
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

#### Pull Out Carbon Metabolic Fxns from CLR data - with NAs ####
## heatmaps of traits of interest

mgm.clr.na[1:4,1:4]

# pull out Carbon metabolism functions from CLR transformed, summed coverages (summed coverage per KO)
carb.ko.na<-mgm.clr.na[,which(colnames(mgm.clr.na) %in% carb.fxns$KO_ID)] # merge CLR data w/ carbon-related fxns found in contigs from KOFamScan
carb.ko.na$SampleID<-rownames(carb.ko.na)
carb.ko.na.melt<-melt(carb.ko.na, by="SampleID")
colnames(carb.ko.na.melt)[which(names(carb.ko.na.melt) == "variable")] <- "KO_ID"
colnames(carb.ko.na.melt)[which(names(carb.ko.na.melt) == "value")] <- "CLR_SumCovPerKO"
head(carb.ko.na.melt) #sanity check

clr.carb.ko.na<-merge(carb.ko.na.melt,carb.kegg,by.x=c("KO_ID"),by.y=c("KO_ID"))
head(clr.carb.ko.na)
colnames(clr.carb.ko.na)[which(names(clr.carb.ko.na) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
clr.cov.sum.carb.ko.na<-as.data.frame(dcast(clr.carb.ko.na, SampleID~KO_Function.KEGG, value.var="CLR_SumCovPerKO", fun.aggregate=sum)) ###
rownames(clr.cov.sum.carb.ko.na)<-clr.cov.sum.carb.ko.na$SampleID
clr.cov.sum.carb.ko.na[1:4,1:4]


### Carbon Heat Maps ####
# see max & mean of summed
max(clr.cov.sum.carb.ko[,-1])
mean(as.matrix(clr.cov.sum.carb.ko[,-1]))

# first heat map of sulfur KOs
heatmap(as.matrix(clr.cov.sum.carb.ko[,-1]), scale = "none")

colSums(clr.cov.sum.carb.ko[,-1])

heatmap(as.matrix(clr.cov.sum.carb.ko[,-1]), scale = "none")

# prep for ggplot2 heatmap
clr.carb.ko.na[1:4,]
clr.carb.all<-merge(clr.carb.ko.na,meta_scaled,by="SampleID")
head(clr.carb.all)
clr.carb.all$PlotID = factor(clr.carb.all$PlotID, levels=unique(clr.carb.all$PlotID[order(clr.carb.all$SampDate,clr.carb.all$Depth_m)]), ordered=TRUE)
clr.carb.all$SampDate<-gsub("\\."," ",clr.carb.all$SampDate)
clr.carb.all$SampDate<-factor(clr.carb.all$SampDate, levels=c("August 2021","December 2021","April 2022"))

unique(clr.carb.all$Pathway)
clr.carb.all<-subset(clr.carb.all, clr.carb.all$Pathway!="Multiple Pathways")
"Multiple Pathways" %in% clr.carb.all$Pathway
clr.carb.all$PathShort<-clr.carb.all$Pathway
clr.carb.all$PathShort[(clr.carb.all$PathShort) == "Reductive Tricarboxylic Acid Cycle"] <- "rTCA"
clr.carb.all$PathShort[(clr.carb.all$PathShort) == "3-Hydroxypropionate Bi-cycle"] <- "3HP"
clr.carb.all$PathShort[(clr.carb.all$PathShort) == "Reductive acetyl-CoA Pathway"] <- "RAcCoa"
clr.carb.all$PathShort[(clr.carb.all$PathShort) == "Calvin Cycle"] <- "CBB"

clr.carb.all$Pathway<-factor(clr.carb.all$Pathway,levels=c("3-Hydroxypropionate Bi-cycle","Reductive Tricarboxylic Acid Cycle","Reductive acetyl-CoA Pathway","Calvin Cycle"))
clr.carb.all$PathShort<-factor(clr.carb.all$PathShort,levels=c("3HP","rTCA","RAcCoa","CBB"))

clr.carb.all$KO_Function.KEGG = factor(clr.carb.all$KO_Function.KEGG, levels=unique(clr.carb.all$KO_Function.KEGG[order(clr.carb.all$Pathway)]), ordered=TRUE)

head(clr.carb.all)

# For heatmap color gradient
max(clr.carb.all$CLR_SumCovPerKO, na.rm=TRUE)
max(clr.carb.all$CLR_SumCovPerKO, na.rm=TRUE)/2
min(clr.carb.all$CLR_SumCovPerKO, na.rm=TRUE)

# Figures below
# by SampleID

carb.hm1a<-ggplot(clr.carb.all, aes(PlotID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.8","0.7","-0.4"),breaks=c(1.8,0.7,-0.4)) + labs(title="Carbon Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(carb.hm1a,filename = "figures/MGM_Figs/FxnDiv/Carbon/Carbon_KOFxns_MGMs_SampID_by_Function_heatmap.png", width=18, height=13, dpi=600)

carb.hm1a2<-ggplot(clr.carb.all, aes(PlotID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.8","0.7","-0.4"),breaks=c(1.8,0.7,-0.4)) + labs(title="Carbon Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Pathway~.,scales="free_y", space = "free")

ggsave(carb.hm1a2,filename = "figures/MGM_Figs/FxnDiv/Carbon/Carbon_KOFxns_MGMs_SampID_by_Function_Pathway_heatmap.png", width=17, height=15, dpi=600)

carb.hm1a3<-ggplot(clr.carb.all, aes(PlotID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.8","0.7","-0.4"),breaks=c(1.8,0.7,-0.4)) + labs(title="Carbon Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate,scales="free_x", space = "free")

ggsave(carb.hm1a3,filename = "figures/MGM_Figs/FxnDiv/Carbon/Carbon_KOFxns_MGMs_SampleID_by_Function_SampDate_best_heatmap.png", width=20, height=13, dpi=600)

carb.hm1a4<-ggplot(clr.carb.all, aes(PlotID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.8","0.7","-0.4"),breaks=c(1.8,0.7,-0.4)) + labs(title="Carbon Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~.,scales="free_y", space = "free")

ggsave(carb.hm1a4,filename = "figures/MGM_Figs/FxnDiv/Carbon/Carbon_KOFxns_MGMs_SampID_by_Function_Pathway_heatmap2.png", width=17, height=15, dpi=600)

carb.hm1a5<-ggplot(clr.carb.all, aes(PlotID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.8","0.7","-0.4"),breaks=c(1.8,0.7,-0.4)) + labs(title="Carbon Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Pathway~SampDate, scales="free", space = "free")

ggsave(carb.hm1a5,filename = "figures/MGM_Figs/FxnDiv/Carbon/Carbon_KOFxns_MGMs_SampleID_by_Function_SampDate_Pathway_best_heatmap.png", width=20, height=15, dpi=600)

carb.hm1a6<-ggplot(clr.carb.all, aes(PlotID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.8","0.7","-0.4"),breaks=c(1.8,0.7,-0.4)) + labs(title="Carbon Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~SampDate, scales="free", space = "free")

ggsave(carb.hm1a6,filename = "figures/MGM_Figs/FxnDiv/Carbon/Carbon_KOFxns_MGMs_SampleID_by_Function_SampDate_Pathway_best_heatmap2.png", width=20, height=15, dpi=600)

# by Depth
carb.hm1b<-ggplot(clr.carb.all, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.8","0.7","-0.4"),breaks=c(1.8,0.7,-0.4)) + labs(title="Carbon Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(carb.hm1b,filename = "figures/MGM_Figs/FxnDiv/Carbon/Carbon_KOFxns_MGMs_Depth_by_Function_heatmap.png", width=18, height=13, dpi=600)

carb.hm1b2<-ggplot(clr.carb.all, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.8","0.7","-0.4"),breaks=c(1.8,0.7,-0.4)) + labs(title="Carbon Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Pathway~.,scales="free_y", space = "free")

ggsave(carb.hm1b2,filename = "figures/MGM_Figs/FxnDiv/Carbon/Carbon_KOFxns_MGMs_Depth_by_Function_Pathway_heatmap.png", width=17, height=15, dpi=600)

carb.hm1b3<-ggplot(clr.carb.all, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.8","0.7","-0.4"),breaks=c(1.8,0.7,-0.4)) + labs(title="Carbon Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate,scales="free_x", space = "free")

ggsave(carb.hm1b3,filename = "figures/MGM_Figs/FxnDiv/Carbon/Carbon_KOFxns_MGMs_Depth_by_Function_SampDate_best_heatmap.png", width=20, height=13, dpi=600)

carb.hm1b4<-ggplot(clr.carb.all, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.8","0.7","-0.4"),breaks=c(1.8,0.7,-0.4)) + labs(title="Carbon Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~.,scales="free_y", space = "free")

ggsave(carb.hm1b4,filename = "figures/MGM_Figs/FxnDiv/Carbon/Carbon_KOFxns_MGMs_Depth_by_Function_Pathway_heatmap2.png", width=17, height=15, dpi=600)

carb.hm1b5<-ggplot(clr.carb.all, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.8","0.7","-0.4"),breaks=c(1.8,0.7,-0.4)) + labs(title="Carbon Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Pathway~SampDate, scales="free", space = "free")

ggsave(carb.hm1b5,filename = "figures/MGM_Figs/FxnDiv/Carbon/Carbon_KOFxns_MGMs_Depth_by_Function_SampDate_Pathway_best_heatmap.png", width=20, height=15, dpi=600)

carb.hm1b6<-ggplot(clr.carb.all, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.8","0.7","-0.4"),breaks=c(1.8,0.7,-0.4)) + labs(title="Carbon Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~SampDate, scales="free", space = "free")

ggsave(carb.hm1b6,filename = "figures/MGM_Figs/FxnDiv/Carbon/Carbon_KOFxns_MGMs_Depth_by_Function_SampDate_Pathway_best_heatmap2.png", width=20, height=15, dpi=600)

# carb.hm1e<-ggplot(clr.carb.all[clr.carb.all$Depth_m==0,], aes(Pathway, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.8","0.7","-0.4"),breaks=c(1.8,0.7,-0.4)) + labs(title="Carbon Metabolism in Salton Seawater Metagenomes - 0m",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11,face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)
#
# ggsave(carb.hm1e,filename = "figures/MGM_Figs/FxnDiv/Carbon/Carbon_KOFxns_Pathways_MGMs_0m_heatmap.png", width=18, height=18, dpi=600)
#
# carb.hm1f<-ggplot(clr.carb.all[clr.carb.all$Depth_m==5,], aes(Pathway, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.8","0.7","-0.4"),breaks=c(1.8,0.7,-0.4)) + labs(title="Carbon Metabolism in Salton Seawater Metagenomes - 5m",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)
#
# ggsave(carb.hm1f,filename = "figures/MGM_Figs/FxnDiv/Carbon/Carbon_KOFxns_Pathways_MGMs_5m_heatmap.png", width=18, height=18, dpi=600)
#
# carb.hm1g<-ggplot(clr.carb.all[clr.carb.all$Depth_m==10,], aes(Pathway, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.8","0.7","-0.4"),breaks=c(1.8,0.7,-0.4)) + labs(title="Carbon Metabolism in Salton Seawater Metagenomes - 10m",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(hjust=1,angle=45),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)
#
# ggsave(carb.hm1g,filename = "figures/MGM_Figs/FxnDiv/Carbon/Carbon_KOFxns_Pathways_MGMs_10m_heatmap.png", width=18, height=18, dpi=600)

#### Pull out Carbon Metabolic Fxns from Binary Data ####
carb.ko.bi<-mgm_fxn.binary[,which(colnames(mgm_fxn.binary) %in% carb.fxns$KO_ID)] # merge CLR data w/ N fxns found in contigs from KOFamScan
carb.ko.bi$SampleID<-rownames(carb.ko.bi)
carb.ko.bi.melt<-melt(carb.ko.bi, by="SampleID")
colnames(carb.ko.bi.melt)[which(names(carb.ko.bi.melt) == "variable")] <- "KO_ID"
colnames(carb.ko.bi.melt)[which(names(carb.ko.bi.melt) == "value")] <- "PresAb"
head(carb.ko.bi.melt) #sanity check

clr.carb.ko.bi<-merge(carb.ko.bi.melt,carb.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
head(clr.carb.ko.bi)
colnames(clr.carb.ko.bi)[which(names(clr.carb.ko.bi) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
clr.cov.sum.carb.ko.bi<-as.data.frame(dcast(clr.carb.ko.bi, SampleID~KO_Function.KEGG, value.var="PresAb", fun.aggregate=sum)) ###just dcast, nothing is being added here!
rownames(clr.cov.sum.carb.ko.bi)<-clr.cov.sum.carb.ko.bi$SampleID
clr.cov.sum.carb.ko.bi[1:4,]

# sanity check
clr.cov.sum.carb.ko.bi$`accB, bccP; acetyl-CoA carboxylase biotin carboxyl carrier protein`[1:4]
head(clr.cov.sum.carb.ko.bi)

#### Carbon Binary Heat Maps ####
# prep for ggplot2 heatmap
clr.carb.ko.bi[1:4,]
clr.carb.all.bi<-merge(clr.carb.ko.bi,meta_scaled,by="SampleID")

head(clr.carb.all.bi)
clr.carb.all.bi$PlotID = factor(clr.carb.all.bi$PlotID, levels=unique(clr.carb.all.bi$PlotID[order(clr.carb.all.bi$SampDate,clr.carb.all.bi$Depth_m)]), ordered=TRUE)
clr.carb.all.bi$SampDate<-gsub("\\."," ",clr.carb.all.bi$SampDate)
clr.carb.all.bi$SampDate<-factor(clr.carb.all.bi$SampDate, levels=c("August 2021","December 2021","April 2022"))

unique(clr.carb.all.bi$Pathway)
clr.carb.all.bi<-subset(clr.carb.all.bi, clr.carb.all.bi$Pathway!="Multiple Pathways")
"Multiple Pathways" %in% clr.carb.all.bi$Pathway
clr.carb.all.bi$PathShort<-clr.carb.all.bi$Pathway
clr.carb.all.bi$PathShort[(clr.carb.all.bi$PathShort) == "Reductive Tricarboxylic Acid Cycle"] <- "rTCA"
clr.carb.all.bi$PathShort[(clr.carb.all.bi$PathShort) == "3-Hydroxypropionate Bi-cycle"] <- "3HP"
clr.carb.all.bi$PathShort[(clr.carb.all.bi$PathShort) == "Reductive acetyl-CoA Pathway"] <- "RAcCoa"
clr.carb.all.bi$PathShort[(clr.carb.all.bi$PathShort) == "Calvin Cycle"] <- "CBB"

clr.carb.all.bi$Pathway<-factor(clr.carb.all.bi$Pathway,levels=c("3-Hydroxypropionate Bi-cycle","Reductive Tricarboxylic Acid Cycle","Reductive acetyl-CoA Pathway","Calvin Cycle"))
clr.carb.all.bi$PathShort<-factor(clr.carb.all.bi$PathShort,levels=c("3HP","rTCA","RAcCoa","CBB"))

clr.carb.all.bi$KO_Function.KEGG = factor(clr.carb.all.bi$KO_Function.KEGG, levels=unique(clr.carb.all.bi$KO_Function.KEGG[order(clr.carb.all.bi$Pathway)]), ordered=TRUE)

head(clr.carb.all.bi)

# Figures
carb.bi.hm1a<-ggplot(clr.carb.all.bi, aes(PlotID, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Carbon Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(carb.bi.hm1a,filename = "figures/MGM_Figs/FxnDiv/Carbon/PresenceAbsence/Carbon_KOFxns_MGMs_by_Function_Binary_heatmap.png", width=18, height=13, dpi=600)

carb.bi.hm1b<-ggplot(clr.carb.all.bi, aes(PlotID, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Carbon Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~.,scales="free_y", space = "free")

ggsave(carb.bi.hm1b,filename = "figures/MGM_Figs/FxnDiv/Carbon/PresenceAbsence/Carbon_KOFxns_MGMs_by_Function_Pathway_Binary_heatmap.png", width=17, height=15, dpi=600)

carb.bi.hm1b2<-ggplot(clr.carb.all.bi, aes(Depth_m, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Carbon Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~.,scales="free_y", space = "free")

ggsave(carb.bi.hm1b2,filename = "figures/MGM_Figs/FxnDiv/Carbon/PresenceAbsence/Carbon_KOFxns_MGMs_Depth_by_Function_Pathway_Binary_heatmap.png", width=17, height=15, dpi=600)

carb.bi.hm1d<-ggplot(clr.carb.all.bi, aes(PlotID, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Carbon Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate,scales="free_x", space = "free")

ggsave(carb.bi.hm1d,filename = "figures/MGM_Figs/FxnDiv/Carbon/PresenceAbsence/Carbon_KOFxns_MGMs_by_Function_SampDate_Binary_best_heatmap.png", width=20, height=13, dpi=600)

carb.bi.hm1d2<-ggplot(clr.carb.all.bi, aes(Depth_m, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Carbon Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate,scales="free_x", space = "free")

ggsave(carb.bi.hm1d2,filename = "figures/MGM_Figs/FxnDiv/Carbon/PresenceAbsence/Carbon_KOFxns_MGMs_Depth_by_Function_SampDate_Binary_best_heatmap.png", width=20, height=13, dpi=600)

carb.bi.hm1e<-ggplot(clr.carb.all.bi, aes(PlotID, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Carbon Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~SampDate, scales="free", space = "free")

ggsave(carb.bi.hm1e,filename = "figures/MGM_Figs/FxnDiv/Carbon/PresenceAbsence/Carbon_KOFxns_MGMs_by_Function_SampDate_Pathway_Binary_best_heatmap.png", width=20, height=15, dpi=600)

carb.bi.hm1e2<-ggplot(clr.carb.all.bi, aes(Depth_m, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Carbon Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~SampDate, scales="free", space = "free")

ggsave(carb.bi.hm1e,filename = "figures/MGM_Figs/FxnDiv/Carbon/PresenceAbsence/Carbon_KOFxns_MGMs_Depth_by_Function_SampDate_Pathway_Binary_best_heatmap.png", width=20, height=15, dpi=600)

carb.bi.hm1e2<-ggplot(clr.carb.all.bi, aes(Depth_m, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Carbon Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~SampDate, scales="free", space = "free")

ggsave(carb.bi.hm1e,filename = "figures/MGM_Figs/FxnDiv/Carbon/PresenceAbsence/Carbon_KOFxns_MGMs_Depth_by_Function_SampDate_Pathway_Binary_best_heatmap.png", width=20, height=15, dpi=600)

# carb.bi.hm1e0<-ggplot(clr.carb.all.bi[clr.carb.all.bi$Depth_m==0,], aes(PathShort, KO_Function.KEGG, fill=factor(PresAb))) +
#   geom_tile(colour="black",size=0.25) +
#   scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Carbon Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate)
#
# ggsave(carb.bi.hm1e0,filename = "figures/MGM_Figs/FxnDiv/Carbon/PresenceAbsence/Carbon_KOFxns_Pathways_Binary_0m_heatmap.png", width=18, height=18, dpi=600)
#
# carb.bi.hm1e5<-ggplot(clr.carb.all.bi[clr.carb.all.bi$Depth_m==5,], aes(PathShort, KO_Function.KEGG, fill=factor(PresAb))) +
#   geom_tile(colour="black",size=0.25) +
#   scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Carbon Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)
#
# ggsave(carb.bi.hm1e5,filename = "figures/MGM_Figs/FxnDiv/Carbon/PresenceAbsence/Carbon_KOFxns_Pathways_Binary_5m_heatmap.png", width=18, height=18, dpi=600)
#
# carb.bi.hm1e6<-ggplot(clr.carb.all.bi[clr.carb.all.bi$Depth_m==10,], aes(PathShort, KO_Function.KEGG, fill=factor(PresAb))) +
#   geom_tile(colour="black",size=0.25) +
#   scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Carbon Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)
#
# ggsave(carb.bi.hm1e6,filename = "figures/MGM_Figs/FxnDiv/Carbon/PresenceAbsence/Carbon_KOFxns_Pathways_Binary_10m_heatmap.png", width=18, height=18, dpi=600)

#### Pull Out Phototrophy Fxns from CLR data - NO NAs ####
## heatmaps of traits of interest

mgm.clr[1:4,1:4]

# pull out Phototrophy functions from CLR transformed, summed coverages (summed coverage per KO)
photo.ko<-mgm.clr[,which(colnames(mgm.clr) %in% photo.fxn$KO_ID)] # merge CLR data w/ photoon-related fxns found in contigs from KOFamScan
photo.ko$SampleID<-rownames(photo.ko)
photo.ko.melt<-melt(photo.ko, by="SampleID")
colnames(photo.ko.melt)[which(names(photo.ko.melt) == "variable")] <- "KO_ID"
colnames(photo.ko.melt)[which(names(photo.ko.melt) == "value")] <- "CLR_SumCovPerKO"
head(photo.ko.melt) #sanity check

clr.photo.ko<-merge(photo.ko.melt,photo.kegg,by.x=c("KO_ID"),by.y=c("KO_ID"))
head(clr.photo.ko)
colnames(clr.photo.ko)[which(names(clr.photo.ko) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
clr.cov.sum.photo.ko<-as.data.frame(dcast(clr.photo.ko, SampleID~KO_Function.KEGG, value.var="CLR_SumCovPerKO", fun.aggregate=sum)) ###
rownames(clr.cov.sum.photo.ko)<-clr.cov.sum.photo.ko$SampleID
clr.cov.sum.photo.ko[1:4,1:4]

#### Pull Out Phototrophy Fxns from CLR data - with NAs ####
## heatmaps of traits of interest

mgm.clr.na[1:4,1:4]

# pull out Phototrophy functions from CLR transformed, summed coverages (summed coverage per KO)
photo.ko.na<-mgm.clr.na[,which(colnames(mgm.clr.na) %in% photo.fxn$KO_ID)] # merge CLR data w/ photoon-related fxns found in contigs from KOFamScan
photo.ko.na$SampleID<-rownames(photo.ko.na)
photo.ko.na.melt<-melt(photo.ko.na, by="SampleID")
colnames(photo.ko.na.melt)[which(names(photo.ko.na.melt) == "variable")] <- "KO_ID"
colnames(photo.ko.na.melt)[which(names(photo.ko.na.melt) == "value")] <- "CLR_SumCovPerKO"
head(photo.ko.na.melt) #sanity check

clr.photo.ko.na<-merge(photo.ko.na.melt,photo.kegg,by.x=c("KO_ID"),by.y=c("KO_ID"))
head(clr.photo.ko.na)
colnames(clr.photo.ko.na)[which(names(clr.photo.ko.na) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
clr.cov.sum.photo.ko.na<-as.data.frame(dcast(clr.photo.ko.na, SampleID~KO_Function.KEGG, value.var="CLR_SumCovPerKO", fun.aggregate=sum)) ###
rownames(clr.cov.sum.photo.ko.na)<-clr.cov.sum.photo.ko.na$SampleID
clr.cov.sum.photo.ko.na[1:4,1:4]

### Phototrophy Heat Maps ####
# see max & mean of summed
max(clr.cov.sum.photo.ko[,-1])
mean(as.matrix(clr.cov.sum.photo.ko[,-1]))

# first heat map of sulfur KOs
heatmap(as.matrix(clr.cov.sum.photo.ko[,-1]), scale = "none")

colSums(clr.cov.sum.photo.ko[,-1])

heatmap(as.matrix(clr.cov.sum.photo.ko[,-1]), scale = "none")

# prep for ggplot2 heatmap
clr.photo.ko.na[1:4,]
clr.photo.all<-merge(clr.photo.ko.na,meta_scaled,by="SampleID")
head(clr.photo.all)
clr.photo.all$PlotID = factor(clr.photo.all$PlotID, levels=unique(clr.photo.all$PlotID[order(clr.photo.all$SampDate,clr.photo.all$Depth_m)]), ordered=TRUE)
clr.photo.all$SampDate<-gsub("\\."," ",clr.photo.all$SampDate)
clr.photo.all$SampDate<-factor(clr.photo.all$SampDate, levels=c("August 2021","December 2021","April 2022"))

unique(clr.photo.all$Pathway)
clr.photo.all$PathShort<-clr.photo.all$Pathway
clr.photo.all$PathShort[(clr.photo.all$PathShort) == "Proteorhodopsin"] <- "PR"
#clr.photo.all$PathShort[(clr.photo.all$PathShort) == "Bacteriorhodopsin"] <- "BR"
clr.photo.all$PathShort[(clr.photo.all$PathShort) == "Sensory Rhodopsin"] <- "SR"
#clr.photo.all$PathShort[(clr.photo.all$PathShort) == "Halorhodopsin"] <- "HR"
clr.photo.all$PathShort[(clr.photo.all$PathShort) == "Photosystem II"] <- "PS II"
clr.photo.all$PathShort[(clr.photo.all$PathShort) == "Photosystem I"] <- "PS I"
clr.photo.all$PathShort[(clr.photo.all$PathShort) == "Anoxygenic Photosystem II"] <- "AnOx PS"

clr.photo.all$Pathway<-factor(clr.photo.all$Pathway,levels=c("Proteorhodopsin","Sensory Rhodopsin","Photosystem II","Photosystem I","Anoxygenic Photosystem II"))
clr.photo.all$PathShort<-factor(clr.photo.all$PathShort,levels=c("PR","SR","PS II","PS I","AnOx PS"))

clr.photo.all$MethShort<-clr.photo.all$Method
clr.photo.all$MethShort[(clr.photo.all$MethShort) == "Bacterial Rhodopsin"] <- "Bac Rhod"
clr.photo.all$MethShort[(clr.photo.all$MethShort) == "Anoxygenic Photosynthesis"] <- "AnOx PS"
clr.photo.all$MethShort[(clr.photo.all$MethShort) == "Oxygenic Photosynthesis"] <- "Ox PS"

clr.photo.all$Method<-factor(clr.photo.all$Method,levels=c("Bacterial Rhodopsin","Oxygenic Photosynthesis","Anoxygenic Photosynthesis"))
clr.photo.all$MethShort<-factor(clr.photo.all$MethShort,levels=c("Bac Rhod","Ox PS","AnOx PS"))

unique(clr.photo.all$Phototrophy)
clr.photo.all$Phototrophy<-factor(clr.photo.all$Phototrophy,levels=c("Hetero","Auto"))

clr.photo.all$KO_Function.KEGG = factor(clr.photo.all$KO_Function.KEGG, levels=unique(clr.photo.all$KO_Function.KEGG[order(clr.photo.all$Phototrophy)]), ordered=TRUE)

head(clr.photo.all)

# For heatmap color gradient
max(clr.photo.all$CLR_SumCovPerKO, na.rm=TRUE)
median(clr.photo.all$CLR_SumCovPerKO, na.rm=TRUE)
min(clr.photo.all$CLR_SumCovPerKO, na.rm=TRUE)

# Figures below
# by SampleID

photo.hm1a<-ggplot(clr.photo.all, aes(PlotID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.25","0.5","-0.3"),breaks=c(1.25,0.5,-0.3)) + labs(title="Phototrophy Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(photo.hm1a,filename = "figures/MGM_Figs/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_SampID_by_Function_heatmap.png", width=18, height=13, dpi=600)

photo.hm1b<-ggplot(clr.photo.all, aes(PlotID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.25","0.5","-0.3"),breaks=c(1.25,0.5,-0.3)) + labs(title="Phototrophy Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Phototrophy~.,scales="free_y", space = "free")

ggsave(photo.hm1b,filename = "figures/MGM_Figs/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_SampID_by_Function_Phototrophy_heatmap.png", width=17, height=15, dpi=600)

photo.hm1c<-ggplot(clr.photo.all, aes(PlotID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.25","0.5","-0.3"),breaks=c(1.25,0.5,-0.3)) + labs(title="Phototrophy Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate,scales="free_x", space = "free")

ggsave(photo.hm1c,filename = "figures/MGM_Figs/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_SampleID_by_Function_SampDate_best_heatmap.png", width=20, height=13, dpi=600)

photo.hm1d<-ggplot(clr.photo.all, aes(PlotID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.25","0.5","-0.3"),breaks=c(1.25,0.5,-0.3)) + labs(title="Phototrophy Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~.,scales="free_y", space = "free")

ggsave(photo.hm1d,filename = "figures/MGM_Figs/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_SampID_by_Function_Phototrophy_System_heatmap2.png", width=17, height=15, dpi=600)

photo.hm1e<-ggplot(clr.photo.all, aes(PlotID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.25","0.5","-0.3"),breaks=c(1.25,0.5,-0.3)) + labs(title="Phototrophy Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Phototrophy~SampDate, scales="free", space = "free")

ggsave(photo.hm1e,filename = "figures/MGM_Figs/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_SampleID_by_Function_SampDate_Phototrophy_best_heatmap.png", width=20, height=15, dpi=600)

photo.hm1f<-ggplot(clr.photo.all, aes(PlotID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.25","0.5","-0.3"),breaks=c(1.25,0.5,-0.3)) + labs(title="Phototrophy Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~SampDate, scales="free", space = "free")

ggsave(photo.hm1f,filename = "figures/MGM_Figs/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_SampleID_by_Function_SampDate_Phototrophy_System_best_heatmap2.png", width=20, height=15, dpi=600)

photo.hm1g<-ggplot(clr.photo.all, aes(PlotID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.25","0.5","-0.3"),breaks=c(1.25,0.5,-0.3)) + labs(title="Phototrophy Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(MethShort~SampDate, scales="free", space = "free")

ggsave(photo.hm1g,filename = "figures/MGM_Figs/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_SampleID_by_Function_SampDate_Phototrophy_Method_best_heatmap.png", width=20, height=15, dpi=600)

# by Depth
# photo.hm1b<-ggplot(clr.photo.all, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.25","0.5","-0.3"),breaks=c(1.25,0.5,-0.3)) + labs(title="Phototrophy Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))
#
# ggsave(photo.hm1b,filename = "figures/MGM_Figs/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_Depth_by_Function_heatmap.png", width=18, height=13, dpi=600)

photo.hm1b2<-ggplot(clr.photo.all, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.25","0.5","-0.3"),breaks=c(1.25,0.5,-0.3)) + labs(title="Phototrophy Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Phototrophy~.,scales="free_y", space = "free")

ggsave(photo.hm1b2,filename = "figures/MGM_Figs/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_Depth_by_Function_Phototrophy_heatmap.png", width=17, height=15, dpi=600)

photo.hm1c2<-ggplot(clr.photo.all, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.25","0.5","-0.3"),breaks=c(1.25,0.5,-0.3)) + labs(title="Phototrophy Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate,scales="free_x", space = "free")

ggsave(photo.hm1c2,filename = "figures/MGM_Figs/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_Depth_by_Function_SampDate_best_heatmap.png", width=20, height=13, dpi=600)

photo.hm1d2<-ggplot(clr.photo.all, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.25","0.5","-0.3"),breaks=c(1.25,0.5,-0.3)) + labs(title="Phototrophy Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~.,scales="free_y", space = "free")

ggsave(photo.hm1d2,filename = "figures/MGM_Figs/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_Depth_by_Function_Phototrophy_System_heatmap2.png", width=17, height=15, dpi=600)

photo.hm1e2<-ggplot(clr.photo.all, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.25","0.5","-0.3"),breaks=c(1.25,0.5,-0.3)) + labs(title="Phototrophy Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Phototrophy~SampDate, scales="free", space = "free")

ggsave(photo.hm1e2,filename = "figures/MGM_Figs/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_Depth_by_Function_SampDate_Phototrophy_best_heatmap.png", width=20, height=15, dpi=600)

photo.hm1f2<-ggplot(clr.photo.all, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.25","0.5","-0.3"),breaks=c(1.25,0.5,-0.3)) + labs(title="Phototrophy Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~SampDate, scales="free", space = "free")

ggsave(photo.hm1f2,filename = "figures/MGM_Figs/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_Depth_by_Function_SampDate_Phototrophy_System_best_heatmap2.png", width=20, height=15, dpi=600)

photo.hm1e2<-ggplot(clr.photo.all, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.25","0.5","-0.3"),breaks=c(1.25,0.5,-0.3)) + labs(title="Phototrophy Metabolism in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(MethShort~SampDate, scales="free", space = "free")

ggsave(photo.hm1e2,filename = "figures/MGM_Figs/FxnDiv/Phototrophy/Phototrophy_KOFxns_MGMs_Depth_by_Function_SampDate_Phototrophy_Method_best_heatmap.png", width=20, height=15, dpi=600)


#### Look at Specific Phototrophy Genes ####
# Note: Must run section "Pull Out photour Metabolic Fxns from CLR data" before running this section
head(clr.cov.sum.photo.ko.na) # columns are genes in this df

# merge with scaled metadata and prep for scatterplots of traits across samples
clr.photo.trait.table<-merge(clr.cov.sum.photo.ko.na,meta_scaled,by="SampleID")
head(clr.photo.trait.table)
clr.photo.trait.table$PlotID = factor(clr.photo.trait.table$PlotID, levels=unique(clr.photo.trait.table$PlotID[order(clr.photo.trait.table$SampDate,clr.photo.trait.table$Depth_m)]), ordered=TRUE)
clr.photo.trait.table$SampDate<-gsub("\\."," ",clr.photo.trait.table$SampDate)
clr.photo.trait.table$SampDate<-factor(clr.photo.trait.table$SampDate, levels=c("August 2021","December 2021","April 2022"))

head(clr.photo.trait.table)

# Note: not looking at every S cycling gene included in this project but looking at ones that appear to have noticeable trends in heat maps

# `blh; beta-carotene 15,15'-dioxygenase [EC:1.13.11.63]`
blh.scat<-ggplot(clr.photo.trait.table, aes(x=PlotID, y=`blh; beta-carotene 15,15'-dioxygenase [EC:1.13.11.63]`,color=SampDate,group=SampDate)) + geom_point(size=4) + geom_line() + theme_bw()+
  labs(title="Blh Depth of Coverage in Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(clr.photo.trait.table$SampDate_Color[order(clr.photo.trait.table$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("SampleID") + ylab("CLR-Transformed Coverage")

ggsave(blh.scat,filename = "figures/MGM_Figs/FxnDiv/Phototrophy/Fxn_Scatterplots/Blh_CLR_Coverage_SampleID_scatterplot.png", width=12, height=10, dpi=600)

# `sop; sensory rhodopsin`
sop.scat<-ggplot(clr.photo.trait.table, aes(x=PlotID, y=`sop; sensory rhodopsin`,color=SampDate,group=SampDate)) + geom_point(size=4) + geom_line() + theme_bw()+
  labs(title="Sop Depth of Coverage in Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(clr.photo.trait.table$SampDate_Color[order(clr.photo.trait.table$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("SampleID") + ylab("CLR-Transformed Coverage")

ggsave(sop.scat,filename = "figures/MGM_Figs/FxnDiv/Phototrophy/Fxn_Scatterplots/Sop_CLR_Coverage_SampleID_scatterplot.png", width=12, height=10, dpi=600)

# `psbA; photosystem II P680 reaction center D1 protein [EC:1.10.3.9]`
psbA.scat<-ggplot(clr.photo.trait.table, aes(x=PlotID, y=`psbA; photosystem II P680 reaction center D1 protein [EC:1.10.3.9]`,color=SampDate,group=SampDate)) + geom_point(size=4) + geom_line() + theme_bw()+
  labs(title="PsbA Depth of Coverage in Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(clr.photo.trait.table$SampDate_Color[order(clr.photo.trait.table$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("SampleID") + ylab("CLR-Transformed Coverage")

ggsave(psbA.scat,filename = "figures/MGM_Figs/FxnDiv/Phototrophy/Fxn_Scatterplots/PsbA_CLR_Coverage_SampleID_scatterplot.png", width=12, height=10, dpi=600)

#### Pull out Phototrophy Metabolic Fxns from Binary Data ####
photo.ko.bi<-mgm_fxn.binary[,which(colnames(mgm_fxn.binary) %in% photo.fxn$KO_ID)] # merge CLR data w/ N fxns found in contigs from KOFamScan
photo.ko.bi$SampleID<-rownames(photo.ko.bi)
photo.ko.bi.melt<-melt(photo.ko.bi, by="SampleID")
colnames(photo.ko.bi.melt)[which(names(photo.ko.bi.melt) == "variable")] <- "KO_ID"
colnames(photo.ko.bi.melt)[which(names(photo.ko.bi.melt) == "value")] <- "PresAb"
head(photo.ko.bi.melt) #sanity check

clr.photo.ko.bi<-merge(photo.ko.bi.melt,photo.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
head(clr.photo.ko.bi)
colnames(clr.photo.ko.bi)[which(names(clr.photo.ko.bi) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
clr.cov.sum.photo.ko.bi<-as.data.frame(dcast(clr.photo.ko.bi, SampleID~KO_Function.KEGG, value.var="PresAb", fun.aggregate=sum)) ###just dcast, nothing is being added here!
rownames(clr.cov.sum.photo.ko.bi)<-clr.cov.sum.photo.ko.bi$SampleID
clr.cov.sum.photo.ko.bi[1:4,]

# sanity check
clr.cov.sum.photo.ko.bi$`accB, bccP; acetyl-CoA photooxylase biotin photooxyl carrier protein`[1:4]
head(clr.cov.sum.photo.ko.bi)

#### Phototrophy Binary Heat Maps ####
# prep for ggplot2 heatmap
clr.photo.ko.bi[1:4,]
clr.photo.all.bi<-merge(clr.photo.ko.bi,meta_scaled,by="SampleID")

head(clr.photo.all.bi)
clr.photo.all.bi$PlotID = factor(clr.photo.all.bi$PlotID, levels=unique(clr.photo.all.bi$PlotID[order(clr.photo.all.bi$SampDate,clr.photo.all.bi$Depth_m)]), ordered=TRUE)
clr.photo.all.bi$SampDate<-gsub("\\."," ",clr.photo.all.bi$SampDate)
clr.photo.all.bi$SampDate<-factor(clr.photo.all.bi$SampDate, levels=c("August 2021","December 2021","April 2022"))

unique(clr.photo.all.bi$Pathway)
clr.photo.all.bi$PathShort<-clr.photo.all.bi$Pathway
clr.photo.all.bi$PathShort[(clr.photo.all.bi$PathShort) == "Proteorhodopsin"] <- "PR"
#clr.photo.all.bi$PathShort[(clr.photo.all.bi$PathShort) == "Bacteriorhodopsin"] <- "BR"
clr.photo.all.bi$PathShort[(clr.photo.all.bi$PathShort) == "Sensory Rhodopsin"] <- "SR"
#clr.photo.all.bi$PathShort[(clr.photo.all.bi$PathShort) == "Halorhodopsin"] <- "HR"
clr.photo.all.bi$PathShort[(clr.photo.all.bi$PathShort) == "Photosystem II"] <- "PS II"
clr.photo.all.bi$PathShort[(clr.photo.all.bi$PathShort) == "Photosystem I"] <- "PS I"
clr.photo.all.bi$PathShort[(clr.photo.all.bi$PathShort) == "Anoxygenic Photosystem II"] <- "AnOx PS"

clr.photo.all.bi$Pathway<-factor(clr.photo.all.bi$Pathway,levels=c("Proteorhodopsin","Sensory Rhodopsin","Photosystem II","Photosystem I","Anoxygenic Photosystem II"))
clr.photo.all.bi$PathShort<-factor(clr.photo.all.bi$PathShort,levels=c("PR","SR","PS II","PS I","AnOx PS"))

clr.photo.all.bi$MethShort<-clr.photo.all.bi$Method
clr.photo.all.bi$MethShort[(clr.photo.all.bi$MethShort) == "Bacterial Rhodopsin"] <- "Bac Rhod"
clr.photo.all.bi$MethShort[(clr.photo.all.bi$MethShort) == "Anoxygenic Photosynthesis"] <- "AnOx PS"
clr.photo.all.bi$MethShort[(clr.photo.all.bi$MethShort) == "Oxygenic Photosynthesis"] <- "Ox PS"

clr.photo.all.bi$Method<-factor(clr.photo.all.bi$Method,levels=c("Bacterial Rhodopsin","Oxygenic Photosynthesis","Anoxygenic Photosynthesis"))
clr.photo.all.bi$MethShort<-factor(clr.photo.all.bi$MethShort,levels=c("Bac Rhod","Ox PS","AnOx PS"))

clr.photo.all.bi$Phototrophy<-factor(clr.photo.all.bi$Phototrophy,levels=c("Hetero","Auto"))

clr.photo.all.bi$KO_Function.KEGG = factor(clr.photo.all.bi$KO_Function.KEGG, levels=unique(clr.photo.all.bi$KO_Function.KEGG[order(clr.photo.all.bi$Phototrophy)]), ordered=TRUE)

head(clr.photo.all.bi)

# Figures

# By sample ID
photo.bi.hm1a<-ggplot(clr.photo.all.bi, aes(PlotID, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Phototrophy Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(photo.bi.hm1a,filename = "figures/MGM_Figs/FxnDiv/Phototrophy/PresenceAbsence/Phototrophy_KOFxns_MGMs_by_Function_Binary_heatmap.png", width=18, height=13, dpi=600)

photo.bi.hm1b<-ggplot(clr.photo.all.bi, aes(PlotID, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Phototrophy Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~.,scales="free_y", space = "free")

ggsave(photo.bi.hm1b,filename = "figures/MGM_Figs/FxnDiv/Phototrophy/PresenceAbsence/Phototrophy_KOFxns_MGMs_by_Function_System_Binary_heatmap.png", width=17, height=15, dpi=600)

photo.bi.hm1c<-ggplot(clr.photo.all.bi, aes(PlotID, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Phototrophy Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Phototrophy~.,scales="free_y", space = "free")

ggsave(photo.bi.hm1c,filename = "figures/MGM_Figs/FxnDiv/Phototrophy/PresenceAbsence/Phototrophy_KOFxns_MGMs_by_Function_Phototrophy_Binary_heatmap.png", width=17, height=15, dpi=600)

photo.bi.hm1d<-ggplot(clr.photo.all.bi, aes(PlotID, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Phototrophy Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate,scales="free_x", space = "free")

ggsave(photo.bi.hm1d,filename = "figures/MGM_Figs/FxnDiv/Phototrophy/PresenceAbsence/Phototrophy_KOFxns_MGMs_by_Function_SampDate_Binary_best_heatmap.png", width=20, height=13, dpi=600)

photo.bi.hm1e<-ggplot(clr.photo.all.bi, aes(PlotID, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Phototrophy Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~SampDate, scales="free", space = "free")

ggsave(photo.bi.hm1e,filename = "figures/MGM_Figs/FxnDiv/Phototrophy/PresenceAbsence/Phototrophy_KOFxns_MGMs_by_Function_SampDate_System_Binary_best_heatmap.png", width=20, height=15, dpi=600)

photo.bi.hm1f<-ggplot(clr.photo.all.bi, aes(PlotID, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Phototrophy Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Phototrophy~SampDate, scales="free", space = "free")

ggsave(photo.bi.hm1f,filename = "figures/MGM_Figs/FxnDiv/Phototrophy/PresenceAbsence/Phototrophy_KOFxns_MGMs_by_Function_SampDate_Phototrophy_Binary_best_heatmap.png", width=20, height=15, dpi=600)

photo.bi.hm1g<-ggplot(clr.photo.all.bi, aes(PlotID, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Phototrophy Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(MethShort~SampDate, scales="free", space = "free")

ggsave(photo.bi.hm1g,filename = "figures/MGM_Figs/FxnDiv/Phototrophy/PresenceAbsence/Phototrophy_KOFxns_MGMs_by_Function_SampDate_Photo_Method_Binary_heatmap.png", width=20, height=15, dpi=600)

# By depth

photo.bi.hm1b2<-ggplot(clr.photo.all.bi, aes(Depth_m, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Phototrophy Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~.,scales="free_y", space = "free")

ggsave(photo.bi.hm1b2,filename = "figures/MGM_Figs/FxnDiv/Phototrophy/PresenceAbsence/Phototrophy_KOFxns_MGMs_Depth_by_Function_System_Binary_heatmap.png", width=17, height=15, dpi=600)

photo.bi.hm1c2<-ggplot(clr.photo.all.bi, aes(Depth_m, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Phototrophy Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Phototrophy~.,scales="free_y", space = "free")

ggsave(photo.bi.hm1c2,filename = "figures/MGM_Figs/FxnDiv/Phototrophy/PresenceAbsence/Phototrophy_KOFxns_MGMs_Depth_by_Function_Phototrophy_Binary_heatmap.png", width=17, height=15, dpi=600)

photo.bi.hm1d2<-ggplot(clr.photo.all.bi, aes(Depth_m, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Phototrophy Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate,scales="free_x", space = "free")

ggsave(photo.bi.hm1d2,filename = "figures/MGM_Figs/FxnDiv/Phototrophy/PresenceAbsence/Phototrophy_KOFxns_MGMs_Depth_by_Function_SampDate_Binary_best_heatmap.png", width=20, height=13, dpi=600)

photo.bi.hm1e2<-ggplot(clr.photo.all.bi, aes(Depth_m, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Phototrophy Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~SampDate, scales="free", space = "free")

ggsave(photo.bi.hm1e2,filename = "figures/MGM_Figs/FxnDiv/Phototrophy/PresenceAbsence/Phototrophy_KOFxns_MGMs_Depth_by_Function_SampDate_System_Binary_best_heatmap.png", width=20, height=15, dpi=600)

photo.bi.hm1f2<-ggplot(clr.photo.all.bi, aes(Depth_m, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Phototrophy Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Phototrophy~SampDate, scales="free", space = "free")

ggsave(photo.bi.hm1f2,filename = "figures/MGM_Figs/FxnDiv/Phototrophy/PresenceAbsence/Phototrophy_KOFxns_MGMs_Depth_by_Function_SampDate_Phototrophy_Binary_best_heatmap.png", width=20, height=15, dpi=600)

photo.bi.hm1g2<-ggplot(clr.photo.all.bi, aes(Depth_m, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Phototrophy Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(MethShort~SampDate, scales="free", space = "free")

ggsave(photo.bi.hm1g2,filename = "figures/MGM_Figs/FxnDiv/Phototrophy/PresenceAbsence/Phototrophy_KOFxns_MGMs_Depth_by_Function_SampDate_Photo_Method_Binary_best_heatmap.png", width=20, height=15, dpi=600)

# photo.bi.hm1e0<-ggplot(clr.photo.all.bi[clr.photo.all.bi$Depth_m==0,], aes(PathShort, KO_Function.KEGG, fill=factor(PresAb))) +
#   geom_tile(colour="black",size=0.25) +
#   scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Phototrophy Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate)
#
# ggsave(photo.bi.hm1e0,filename = "figures/MGM_Figs/FxnDiv/Phototrophy/PresenceAbsence/Phototrophy_KOFxns_Phototrophys_Binary_0m_heatmap.png", width=18, height=18, dpi=600)
#
# photo.bi.hm1e5<-ggplot(clr.photo.all.bi[clr.photo.all.bi$Depth_m==5,], aes(PathShort, KO_Function.KEGG, fill=factor(PresAb))) +
#   geom_tile(colour="black",size=0.25) +
#   scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Phototrophy Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)
#
# ggsave(photo.bi.hm1e5,filename = "figures/MGM_Figs/FxnDiv/Phototrophy/PresenceAbsence/Phototrophy_KOFxns_Phototrophys_Binary_5m_heatmap.png", width=18, height=18, dpi=600)
#
# photo.bi.hm1e6<-ggplot(clr.photo.all.bi[clr.photo.all.bi$Depth_m==10,], aes(PathShort, KO_Function.KEGG, fill=factor(PresAb))) +
#   geom_tile(colour="black",size=0.25) +
#   scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Phototrophy Fixation in Salton Seawater Metagenomes",fill="Presence/Absence") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)
#
# ggsave(photo.bi.hm1e6,filename = "figures/MGM_Figs/FxnDiv/Phototrophy/PresenceAbsence/Phototrophy_KOFxns_Phototrophys_Binary_10m_heatmap.png", width=18, height=18, dpi=600)

#### Pull Out Aerobic Respiration Fxns from CLR data ####
## heatmaps of traits of interest

mgm.clr.na[1:4,1:4]

# pull out Aerobic Respiration functions from CLR transformed, summed coverages (summed coverage per KO)
aero.ko<-mgm.clr.na[,which(colnames(mgm.clr.na) %in% aero.fxn$KO_ID)] # merge CLR data w/ aeroon-related fxns found in contigs from KOFamScan
aero.ko$SampleID<-rownames(aero.ko)
aero.ko.melt<-melt(aero.ko, by="SampleID")
colnames(aero.ko.melt)[which(names(aero.ko.melt) == "variable")] <- "KO_ID"
colnames(aero.ko.melt)[which(names(aero.ko.melt) == "value")] <- "CLR_SumCovPerKO"
head(aero.ko.melt) #sanity check

clr.aero.ko<-merge(aero.ko.melt,aero.kegg,by.x=c("KO_ID"),by.y=c("KO_ID"))
head(clr.aero.ko)
colnames(clr.aero.ko)[which(names(clr.aero.ko) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
clr.cov.sum.aero.ko<-as.data.frame(dcast(clr.aero.ko, SampleID~KO_Function.KEGG, value.var="CLR_SumCovPerKO", fun.aggregate=sum)) ###
rownames(clr.cov.sum.aero.ko)<-clr.cov.sum.aero.ko$SampleID
clr.cov.sum.aero.ko[1:4,1:4]

### Aerobic Respiration Heat Maps ####
# see max & mean of summed
max(clr.cov.sum.aero.ko[,-1])
mean(as.matrix(clr.cov.sum.aero.ko[,-1]))

# first heat map of sulfur KOs
heatmap(as.matrix(clr.cov.sum.aero.ko[,-1]), scale = "none")

colSums(clr.cov.sum.aero.ko[,-1])

heatmap(as.matrix(clr.cov.sum.aero.ko[,-1]), scale = "none")

# prep for ggplot2 heatmap
clr.aero.ko[1:4,]
clr.aero.all<-merge(clr.aero.ko,meta_scaled,by="SampleID")
head(clr.aero.all)
clr.aero.all$PlotID = factor(clr.aero.all$PlotID, levels=unique(clr.aero.all$PlotID[order(clr.aero.all$SampDate,clr.aero.all$Depth_m)]), ordered=TRUE)
clr.aero.all$SampDate<-gsub("\\."," ",clr.aero.all$SampDate)
clr.aero.all$SampDate<-factor(clr.aero.all$SampDate, levels=c("August 2021","December 2021","April 2022"))

clr.aero.all$EnzShort<-clr.aero.all$Enzyme
clr.aero.all$EnzShort[(clr.aero.all$EnzShort) == "Cytochrome c oxidase"] <- "Cox"
clr.aero.all$EnzShort[(clr.aero.all$EnzShort) == "F-type ATPase"] <- "F-ATPase"
clr.aero.all$EnzShort[(clr.aero.all$EnzShort) == "NADH:quinone oxidoreductase"] <- "NDH-2"
clr.aero.all$EnzShort[(clr.aero.all$EnzShort) == "Fumarate reductase"] <- "FRD"
clr.aero.all$EnzShort[(clr.aero.all$EnzShort) == "Succinate dehydrogenase"] <- "SDH"

clr.aero.all$Enzyme<-factor(clr.aero.all$Enzyme,levels=c("Cytochrome c oxidase","F-type ATPase","NADH:quinone oxidoreductase","Fumarate reductase","Succinate dehydrogenase"))
clr.aero.all$EnzShort<-factor(clr.aero.all$EnzShort,levels=c("Cox","F-ATPase","NDH-2","FRD","SDH"))

clr.aero.all$KO_Function.KEGG = factor(clr.aero.all$KO_Function.KEGG, levels=unique(clr.aero.all$KO_Function.KEGG[order(clr.aero.all$Enzyme)]), ordered=TRUE)

head(clr.aero.all)

# For heatmap color gradient
max(clr.aero.all$CLR_SumCovPerKO, na.rm=TRUE)
max(clr.aero.all$CLR_SumCovPerKO, na.rm=TRUE)/2
min(clr.aero.all$CLR_SumCovPerKO, na.rm=TRUE)

# Figures below
# by SampleID

aero.hm1a<-ggplot(clr.aero.all, aes(PlotID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.75","0.75","0","-0.5"),breaks=c(1.75,1,0,-0.5)) + labs(title="Aerobic Respiration in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(aero.hm1a,filename = "figures/MGM_Figs/FxnDiv/Aerobic_Respiration/Aerobic_Respiration_KOFxns_MGMs_SampID_by_Function_heatmap.png", width=18, height=13, dpi=600)

# aero.hm1b<-ggplot(clr.aero.all, aes(PlotID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
#   geom_tile(colour="white",size=0.15) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.75","0.75","0","-0.5"),breaks=c(1.75,0.75,0,-0.5)) + labs(title="Aerobic Respiration in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Enzyme~.,scales="free_y", space = "free")
#
# ggsave(aero.hm1b,filename = "figures/MGM_Figs/FxnDiv/Aerobic_Respiration/Aerobic_Respiration_KOFxns_MGMs_SampID_by_Function_Aerobic_Respiration_heatmap.png", width=17, height=15, dpi=600)

aero.hm1c<-ggplot(clr.aero.all, aes(PlotID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.75","0.75","0","-0.5"),breaks=c(1.75,0.75,0,-0.5)) + labs(title="Aerobic Respiration in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate,scales="free_x", space = "free")

ggsave(aero.hm1c,filename = "figures/MGM_Figs/FxnDiv/Aerobic_Respiration/Aerobic_Respiration_KOFxns_MGMs_SampleID_by_Function_SampDate_best_heatmap.png", width=20, height=13, dpi=600)

aero.hm1d<-ggplot(clr.aero.all, aes(PlotID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.75","0.75","0","-0.5"),breaks=c(1.75,0.75,0,-0.5)) + labs(title="Aerobic Respiration in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(EnzShort~.,scales="free_y", space = "free")

ggsave(aero.hm1d,filename = "figures/MGM_Figs/FxnDiv/Aerobic_Respiration/Aerobic_Respiration_KOFxns_MGMs_SampID_by_Function_Aerobic_Respiration_Enzyme_heatmap2.png", width=17, height=15, dpi=600)
#
# aero.hm1e<-ggplot(clr.aero.all, aes(PlotID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.75","0.75","0","-0.5"),breaks=c(1.75,0.75,0,-0.5)) + labs(title="Aerobic Respiration in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Enzyme~SampDate, scales="free", space = "free")
#
# ggsave(aero.hm1e,filename = "figures/MGM_Figs/FxnDiv/Aerobic_Respiration/Aerobic_Respiration_KOFxns_MGMs_SampleID_by_Function_SampDate_Aerobic_Respiration_best_heatmap.png", width=20, height=15, dpi=600)

aero.hm1f<-ggplot(clr.aero.all, aes(PlotID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.75","0.75","0","-0.5"),breaks=c(1.75,0.75,0,-0.5)) + labs(title="Aerobic Respiration in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(EnzShort~SampDate, scales="free", space = "free")

ggsave(aero.hm1f,filename = "figures/MGM_Figs/FxnDiv/Aerobic_Respiration/Aerobic_Respiration_KOFxns_MGMs_SampleID_by_Function_SampDate_Aerobic_Respiration_Enzyme_best_heatmap2.png", width=20, height=15, dpi=600)

# aero.hm1g<-ggplot(clr.aero.all, aes(PlotID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.75","0.75","0","-0.5"),breaks=c(1.75,0.75,0,-0.5)) + labs(title="Aerobic Respiration in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(EnzShort~SampDate, scales="free", space = "free")
#
# ggsave(aero.hm1g,filename = "figures/MGM_Figs/FxnDiv/Aerobic_Respiration/Aerobic_Respiration_KOFxns_MGMs_SampleID_by_Function_SampDate_Aerobic_Respiration_Enzyme_best_heatmap.png", width=20, height=15, dpi=600)

# by Depth
# aero.hm1b<-ggplot(clr.aero.all, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.75","0.75","0","-0.5"),breaks=c(1.75,0.75,0,-0.5)) + labs(title="Aerobic Respiration in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))
#
# ggsave(aero.hm1b,filename = "figures/MGM_Figs/FxnDiv/Aerobic_Respiration/Aerobic_Respiration_KOFxns_MGMs_Depth_by_Function_heatmap.png", width=18, height=13, dpi=600)

aero.hm1b2<-ggplot(clr.aero.all, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.75","0.75","0","-0.5"),breaks=c(1.75,0.75,0,-0.5)) + labs(title="Aerobic Respiration in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(EnzShort~.,scales="free_y", space = "free")

ggsave(aero.hm1b2,filename = "figures/MGM_Figs/FxnDiv/Aerobic_Respiration/Aerobic_Respiration_KOFxns_MGMs_Depth_by_Function_Aerobic_Respiration_Enzyme_heatmap.png", width=17, height=15, dpi=600)

aero.hm1c2<-ggplot(clr.aero.all, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.75","0.75","0","-0.5"),breaks=c(1.75,0.75,0,-0.5)) + labs(title="Aerobic Respiration in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.x = element_text(size = 11)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate,scales="free_x", space = "free")

ggsave(aero.hm1c2,filename = "figures/MGM_Figs/FxnDiv/Aerobic_Respiration/Aerobic_Respiration_KOFxns_MGMs_Depth_by_Function_SampDate_best_heatmap.png", width=20, height=13, dpi=600)

aero.hm1d2<-ggplot(clr.aero.all, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.75","0.75","0","-0.5"),breaks=c(1.75,0.75,0,-0.5)) + labs(title="Aerobic Respiration in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(EnzShort~.,scales="free_y", space = "free")

ggsave(aero.hm1d2,filename = "figures/MGM_Figs/FxnDiv/Aerobic_Respiration/Aerobic_Respiration_KOFxns_MGMs_Depth_by_Function_Aerobic_Respiration_Enzyme_heatmap2.png", width=17, height=15, dpi=600)

# aero.hm1e2<-ggplot(clr.aero.all, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.75","0.75","0","-0.5"),breaks=c(1.75,0.75,0,-0.5)) + labs(title="Aerobic Respiration in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Enzyme~SampDate, scales="free", space = "free")
#
# ggsave(aero.hm1e2,filename = "figures/MGM_Figs/FxnDiv/Aerobic_Respiration/Aerobic_Respiration_KOFxns_MGMs_Depth_by_Function_SampDate_Aerobic_Respiration_best_heatmap.png", width=20, height=15, dpi=600)

aero.hm1f2<-ggplot(clr.aero.all, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.75","0.75","0","-0.5"),breaks=c(1.75,0.75,0,-0.5)) + labs(title="Aerobic Respiration in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(EnzShort~SampDate, scales="free", space = "free")

ggsave(aero.hm1f2,filename = "figures/MGM_Figs/FxnDiv/Aerobic_Respiration/Aerobic_Respiration_KOFxns_MGMs_Depth_by_Function_SampDate_Aerobic_Respiration_Enzyme_best_heatmap2.png", width=20, height=15, dpi=600)

# aero.hm1e2<-ggplot(clr.aero.all, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.75","0.75","0","-0.5"),breaks=c(1.75,0.75,0,-0.5)) + labs(title="Aerobic Respiration in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text = element_text(size = 11),strip.text.y=element_text(face="bold")) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(EnzShort~SampDate, scales="free", space = "free")
#
# ggsave(aero.hm1e2,filename = "figures/MGM_Figs/FxnDiv/Aerobic_Respiration/Aerobic_Respiration_KOFxns_MGMs_Depth_by_Function_SampDate_Aerobic_Respiration_Enzyme_best_heatmap.png", width=20, height=15, dpi=600)

#### Pull out Aerobic Respiration Metabolic Fxns from Binary Data ####
aero.ko.bi<-mgm_fxn.binary[,which(colnames(mgm_fxn.binary) %in% aero.fxn$KO_ID)] # merge CLR data w/ N fxns found in contigs from KOFamScan
aero.ko.bi$SampleID<-rownames(aero.ko.bi)
aero.ko.bi.melt<-melt(aero.ko.bi, by="SampleID")
colnames(aero.ko.bi.melt)[which(names(aero.ko.bi.melt) == "variable")] <- "KO_ID"
colnames(aero.ko.bi.melt)[which(names(aero.ko.bi.melt) == "value")] <- "PresAb"
head(aero.ko.bi.melt) #sanity check

clr.aero.ko.bi<-merge(aero.ko.bi.melt,aero.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
head(clr.aero.ko.bi)
colnames(clr.aero.ko.bi)[which(names(clr.aero.ko.bi) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
clr.cov.sum.aero.ko.bi<-as.data.frame(dcast(clr.aero.ko.bi, SampleID~KO_Function.KEGG, value.var="PresAb", fun.aggregate=sum)) ###just dcast, nothing is being added here!
rownames(clr.cov.sum.aero.ko.bi)<-clr.cov.sum.aero.ko.bi$SampleID
clr.cov.sum.aero.ko.bi[1:4,]

# sanity check
clr.cov.sum.aero.ko.bi$`accB, bccP; acetyl-CoA aerooxylase biotin aerooxyl carrier protein`[1:4]
head(clr.cov.sum.aero.ko.bi)

#### Aerobic Respiration Binary Heat Maps ####
# prep for ggplot2 heatmap
clr.aero.ko.bi[1:4,]
clr.aero.all.bi<-merge(clr.aero.ko.bi,meta_scaled,by="SampleID")

head(clr.aero.all.bi)
clr.aero.all.bi$PlotID = factor(clr.aero.all.bi$PlotID, levels=unique(clr.aero.all.bi$PlotID[order(clr.aero.all.bi$SampDate,clr.aero.all.bi$Depth_m)]), ordered=TRUE)
clr.aero.all.bi$SampDate<-gsub("\\."," ",clr.aero.all.bi$SampDate)
clr.aero.all.bi$SampDate<-factor(clr.aero.all.bi$SampDate, levels=c("August 2021","December 2021","April 2022"))

clr.aero.all.bi$EnzShort<-clr.aero.all.bi$Enzyme
clr.aero.all.bi$EnzShort[(clr.aero.all.bi$EnzShort) == "Cytochrome c oxidase"] <- "Cox"
clr.aero.all.bi$EnzShort[(clr.aero.all.bi$EnzShort) == "F-type ATPase"] <- "F-ATPase"
clr.aero.all.bi$EnzShort[(clr.aero.all.bi$EnzShort) == "NADH:quinone oxidoreductase"] <- "NDH-2"
clr.aero.all.bi$EnzShort[(clr.aero.all.bi$EnzShort) == "Fumarate reductase"] <- "FRD"
clr.aero.all.bi$EnzShort[(clr.aero.all.bi$EnzShort) == "Succinate dehydrogenase"] <- "SDH"

clr.aero.all.bi$Enzyme<-factor(clr.aero.all.bi$Enzyme,levels=c("Cytochrome c oxidase","F-type ATPase","NADH:quinone oxidoreductase","Fumarate reductase","Succinate dehydrogenase"))
clr.aero.all.bi$EnzShort<-factor(clr.aero.all.bi$EnzShort,levels=c("Cox","F-ATPase","NDH-2","FRD","SDH"))

clr.aero.all.bi$KO_Function.KEGG = factor(clr.aero.all.bi$KO_Function.KEGG, levels=unique(clr.aero.all.bi$KO_Function.KEGG[order(clr.aero.all.bi$Enzyme)]), ordered=TRUE)

head(clr.aero.all.bi)

# Figures

# By sample ID
aero.bi.hm1a<-ggplot(clr.aero.all.bi, aes(PlotID, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Aerobic Respiration in Salton Seawater Metagenomes",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(aero.bi.hm1a,filename = "figures/MGM_Figs/FxnDiv/Aerobic_Respiration/PresenceAbsence/Aerobic_Respiration_KOFxns_MGMs_by_Function_Binary_heatmap.png", width=18, height=13, dpi=600)

aero.bi.hm1b<-ggplot(clr.aero.all.bi, aes(PlotID, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Aerobic Respiration in Salton Seawater Metagenomes",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(EnzShort~.,scales="free_y", space = "free")

ggsave(aero.bi.hm1b,filename = "figures/MGM_Figs/FxnDiv/Aerobic_Respiration/PresenceAbsence/Aerobic_Respiration_KOFxns_MGMs_by_Function_Enzyme_Binary_heatmap.png", width=17, height=15, dpi=600)

# aero.bi.hm1c<-ggplot(clr.aero.all.bi, aes(PlotID, KO_Function.KEGG, fill=factor(PresAb))) +
#   geom_tile(colour="black",size=0.25) +
#   scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Aerobic Respiration in Salton Seawater Metagenomes",fill="Presence/Absence") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Enzyme~.,scales="free_y", space = "free")
#
# ggsave(aero.bi.hm1c,filename = "figures/MGM_Figs/FxnDiv/Aerobic_Respiration/PresenceAbsence/Aerobic_Respiration_KOFxns_MGMs_by_Function_Aerobic_Respiration_Binary_heatmap.png", width=17, height=15, dpi=600)

aero.bi.hm1d<-ggplot(clr.aero.all.bi, aes(PlotID, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Aerobic Respiration in Salton Seawater Metagenomes",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate,scales="free_x", space = "free")

ggsave(aero.bi.hm1d,filename = "figures/MGM_Figs/FxnDiv/Aerobic_Respiration/PresenceAbsence/Aerobic_Respiration_KOFxns_MGMs_by_Function_SampDate_Binary_best_heatmap.png", width=20, height=13, dpi=600)

aero.bi.hm1e<-ggplot(clr.aero.all.bi, aes(PlotID, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Aerobic Respiration in Salton Seawater Metagenomes",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(EnzShort~SampDate, scales="free", space = "free")

ggsave(aero.bi.hm1e,filename = "figures/MGM_Figs/FxnDiv/Aerobic_Respiration/PresenceAbsence/Aerobic_Respiration_KOFxns_MGMs_by_Function_SampDate_Enzyme_Binary_best_heatmap.png", width=20, height=15, dpi=600)
#
# aero.bi.hm1f<-ggplot(clr.aero.all.bi, aes(PlotID, KO_Function.KEGG, fill=factor(PresAb))) +
#   geom_tile(colour="black",size=0.25) +
#   scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Aerobic Respiration in Salton Seawater Metagenomes",fill="Presence/Absence") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Enzyme~SampDate, scales="free", space = "free")
#
# ggsave(aero.bi.hm1f,filename = "figures/MGM_Figs/FxnDiv/Aerobic_Respiration/PresenceAbsence/Aerobic_Respiration_KOFxns_MGMs_by_Function_SampDate_Aerobic_Respiration_Binary_best_heatmap.png", width=20, height=15, dpi=600)

# aero.bi.hm1g<-ggplot(clr.aero.all.bi, aes(PlotID, KO_Function.KEGG, fill=factor(PresAb))) +
#   geom_tile(colour="black",size=0.25) +
#   scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Aerobic Respiration in Salton Seawater Metagenomes",fill="Presence/Absence") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(EnzShort~SampDate, scales="free", space = "free")
#
# ggsave(aero.bi.hm1g,filename = "figures/MGM_Figs/FxnDiv/Aerobic_Respiration/PresenceAbsence/Aerobic_Respiration_KOFxns_MGMs_by_Function_SampDate_Photo_Enzyme_Binary_heatmap.png", width=20, height=15, dpi=600)

# By depth

aero.bi.hm1b2<-ggplot(clr.aero.all.bi, aes(Depth_m, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Aerobic Respiration in Salton Seawater Metagenomes",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(EnzShort~.,scales="free_y", space = "free")

ggsave(aero.bi.hm1b2,filename = "figures/MGM_Figs/FxnDiv/Aerobic_Respiration/PresenceAbsence/Aerobic_Respiration_KOFxns_MGMs_Depth_by_Function_Enzyme_Binary_heatmap.png", width=17, height=15, dpi=600)

# aero.bi.hm1c2<-ggplot(clr.aero.all.bi, aes(Depth_m, KO_Function.KEGG, fill=factor(PresAb))) +
#   geom_tile(colour="black",size=0.25) +
#   scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Aerobic Respiration in Salton Seawater Metagenomes",fill="Presence/Absence") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Enzyme~.,scales="free_y", space = "free")
#
# ggsave(aero.bi.hm1c2,filename = "figures/MGM_Figs/FxnDiv/Aerobic_Respiration/PresenceAbsence/Aerobic_Respiration_KOFxns_MGMs_Depth_by_Function_Aerobic_Respiration_Binary_heatmap.png", width=17, height=15, dpi=600)

aero.bi.hm1d2<-ggplot(clr.aero.all.bi, aes(Depth_m, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Aerobic Respiration in Salton Seawater Metagenomes",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate,scales="free_x", space = "free")

ggsave(aero.bi.hm1d2,filename = "figures/MGM_Figs/FxnDiv/Aerobic_Respiration/PresenceAbsence/Aerobic_Respiration_KOFxns_MGMs_Depth_by_Function_SampDate_Binary_best_heatmap.png", width=20, height=13, dpi=600)

aero.bi.hm1e2<-ggplot(clr.aero.all.bi, aes(Depth_m, KO_Function.KEGG, fill=factor(PresAb))) +
  geom_tile(colour="black",size=0.25) +
  scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Aerobic Respiration in Salton Seawater Metagenomes",fill="Presence/Absence") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(EnzShort~SampDate, scales="free", space = "free")

ggsave(aero.bi.hm1e2,filename = "figures/MGM_Figs/FxnDiv/Aerobic_Respiration/PresenceAbsence/Aerobic_Respiration_KOFxns_MGMs_Depth_by_Function_SampDate_Enzyme_Binary_best_heatmap.png", width=20, height=15, dpi=600)

# aero.bi.hm1f2<-ggplot(clr.aero.all.bi, aes(Depth_m, KO_Function.KEGG, fill=factor(PresAb))) +
#   geom_tile(colour="black",size=0.25) +
#   scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Aerobic Respiration in Salton Seawater Metagenomes",fill="Presence/Absence") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Enzyme~SampDate, scales="free", space = "free")
#
# ggsave(aero.bi.hm1f2,filename = "figures/MGM_Figs/FxnDiv/Aerobic_Respiration/PresenceAbsence/Aerobic_Respiration_KOFxns_MGMs_Depth_by_Function_SampDate_Aerobic_Respiration_Binary_best_heatmap.png", width=20, height=15, dpi=600)

# aero.bi.hm1g2<-ggplot(clr.aero.all.bi, aes(Depth_m, KO_Function.KEGG, fill=factor(PresAb))) +
#   geom_tile(colour="black",size=0.25) +
#   scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Aerobic Respiration in Salton Seawater Metagenomes",fill="Presence/Absence") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(EnzShort~SampDate, scales="free", space = "free")
#
# ggsave(aero.bi.hm1g2,filename = "figures/MGM_Figs/FxnDiv/Aerobic_Respiration/PresenceAbsence/Aerobic_Respiration_KOFxns_MGMs_Depth_by_Function_SampDate_Photo_Enzyme_Binary_best_heatmap.png", width=20, height=15, dpi=600)

# aero.bi.hm1e0<-ggplot(clr.aero.all.bi[clr.aero.all.bi$Depth_m==0,], aes(EnzShort, KO_Function.KEGG, fill=factor(PresAb))) +
#   geom_tile(colour="black",size=0.25) +
#   scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Aerobic Respiration in Salton Seawater Metagenomes",fill="Presence/Absence") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(~SampDate)
#
# ggsave(aero.bi.hm1e0,filename = "figures/MGM_Figs/FxnDiv/Aerobic_Respiration/PresenceAbsence/Aerobic_Respiration_KOFxns_Aerobic Respirations_Binary_0m_heatmap.png", width=18, height=18, dpi=600)
#
# aero.bi.hm1e5<-ggplot(clr.aero.all.bi[clr.aero.all.bi$Depth_m==5,], aes(EnzShort, KO_Function.KEGG, fill=factor(PresAb))) +
#   geom_tile(colour="black",size=0.25) +
#   scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Aerobic Respiration in Salton Seawater Metagenomes",fill="Presence/Absence") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)
#
# ggsave(aero.bi.hm1e5,filename = "figures/MGM_Figs/FxnDiv/Aerobic_Respiration/PresenceAbsence/Aerobic_Respiration_KOFxns_Aerobic Respirations_Binary_5m_heatmap.png", width=18, height=18, dpi=600)
#
# aero.bi.hm1e6<-ggplot(clr.aero.all.bi[clr.aero.all.bi$Depth_m==10,], aes(EnzShort, KO_Function.KEGG, fill=factor(PresAb))) +
#   geom_tile(colour="black",size=0.25) +
#   scale_fill_manual(values=binary.cols,labels=c("Present","Absent"),breaks=c(1,0)) + labs(title="Aerobic Respiration in Salton Seawater Metagenomes",fill="Presence/Absence") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)
#
# ggsave(aero.bi.hm1e6,filename = "figures/MGM_Figs/FxnDiv/Aerobic_Respiration/PresenceAbsence/Aerobic_Respiration_KOFxns_Aerobic Respirations_Binary_10m_heatmap.png", width=18, height=18, dpi=600)

#### Pull Out ALL Genes of Interest per Pathway/Cycle from CLR data ####
## heatmaps of traits of interest

mgm.clr.na[1:4,1:4]

# pull out sulfur functions from CLR transformed, summed coverages (summed coverage per KO)
All_GOI.ko<-mgm.clr.na[,which(colnames(mgm.clr.na) %in% All_GOI.fxns$KO_ID)] # merge CLR data w/ All_GOI-related fxns found in contigs from KOFamScan
All_GOI.ko$SampleID<-rownames(All_GOI.ko)
All_GOI.ko.melt<-melt(All_GOI.ko, by="SampleID")
colnames(All_GOI.ko.melt)[which(names(All_GOI.ko.melt) == "variable")] <- "KO_ID"
colnames(All_GOI.ko.melt)[which(names(All_GOI.ko.melt) == "value")] <- "CLR_SumCovPerKO"
head(All_GOI.ko.melt) #sanity check

clr.All_GOI.ko<-merge(All_GOI.ko.melt,all_goi.kegg,by.x=c("KO_ID"),by.y=c("KO_ID"))
head(clr.All_GOI.ko)
colnames(clr.All_GOI.ko)[which(names(clr.All_GOI.ko) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
clr.cov.sum.All_GOI.ko<-as.data.frame(dcast(clr.All_GOI.ko, SampleID~KO_Function.KEGG, value.var="CLR_SumCovPerKO", fun.aggregate=sum)) ###
rownames(clr.cov.sum.All_GOI.ko)<-clr.cov.sum.All_GOI.ko$SampleID
clr.cov.sum.All_GOI.ko[1:4,1:4]

### ALL Genes of Interest per Pathway/Cycle Heat Maps ####
# see max & mean of summed
max(clr.cov.sum.All_GOI.ko[,-1])
mean(as.matrix(clr.cov.sum.All_GOI.ko[,-1]))

# first heat map of sulfur KOs
heatmap(as.matrix(clr.cov.sum.All_GOI.ko[,-1]), scale = "none")

colSums(clr.cov.sum.All_GOI.ko[,-1])

heatmap(as.matrix(clr.cov.sum.All_GOI.ko[,-1]), scale = "none")

# prep for ggplot2 heatmap
clr.All_GOI.ko[1:4,]
clr.All_GOI.all<-merge(clr.All_GOI.ko,meta_scaled,by="SampleID")
head(clr.All_GOI.all)
clr.All_GOI.all$PlotID = factor(clr.All_GOI.all$PlotID, levels=unique(clr.All_GOI.all$PlotID[order(clr.All_GOI.all$SampDate,clr.All_GOI.all$Depth_m)]), ordered=TRUE)
clr.All_GOI.all$SampDate<-gsub("\\."," ",clr.All_GOI.all$SampDate)
clr.All_GOI.all$SampDate<-factor(clr.All_GOI.all$SampDate, levels=c("August 2021","December 2021","April 2022"))
unique(clr.All_GOI.all$Pathway)
unique(clr.All_GOI.all$Cycle)
clr.All_GOI.all$Cycle<-factor(clr.All_GOI.all$Cycle, levels=c("Sulfur Cycle","Carbon Cycle","Nitrogen Cycle","Photoheterotrophy"))

unique(clr.All_GOI.all$Pathway)
clr.All_GOI.all$PathShort<-clr.All_GOI.all$Pathway
# vvv can only do this type of renaming if variables are characters, not factors
clr.All_GOI.all$PathShort[(clr.All_GOI.all$PathShort) == "Assimilatory Sulfate Reduction"] <- "A.SO4 Red"
clr.All_GOI.all$PathShort[(clr.All_GOI.all$PathShort) == "Dissimilatory Sulfate Redox"] <- "D.SO4 Redox"
clr.All_GOI.all$PathShort[(clr.All_GOI.all$PathShort) == "Calvin Cycle"] <- "CBB"
clr.All_GOI.all$PathShort[(clr.All_GOI.all$PathShort) == "Reductive Tricarboxylic Acid Cycle"] <- "rTCA"
clr.All_GOI.all$PathShort[(clr.All_GOI.all$PathShort) == "3-Hydroxypropionate Bi-cycle"] <- "3-H BC"
clr.All_GOI.all$PathShort[(clr.All_GOI.all$PathShort) == "Phosphate acetyltransferase-acetate kinase Pathway"] <- "P.A.A.K."
clr.All_GOI.all$PathShort[(clr.All_GOI.all$PathShort) == "Assimilatory Nitrate Reduction"] <- "A.NO3 Red"
clr.All_GOI.all$PathShort[(clr.All_GOI.all$PathShort) == "Photoheterotrophy"] <- "PhotoHet"
clr.All_GOI.all$PathShort[(clr.All_GOI.all$PathShort) == "Methanogenesis"] <- "MethGen"
clr.All_GOI.all$PathShort[(clr.All_GOI.all$PathShort) == "Denitrification"] <- "DeNit"

unique(clr.All_GOI.all$Pathway)
clr.All_GOI.all$Pathway<-factor(clr.All_GOI.all$Pathway,levels=c("Assimilatory Sulfate Reduction","Dissimilatory Sulfate Redox","SOX",
                                                         "Reductive Tricarboxylic Acid Cycle","3-Hydroxypropionate Bi-cycle","Methanogenesis",
                                                         "Calvin Cycle","Denitrification","Assimilatory Nitrate Reduction","Anammox","Photoheterotrophy"))

unique(clr.All_GOI.all$PathShort)
clr.All_GOI.all$PathShort<-factor(clr.All_GOI.all$PathShort,levels=c("A.SO4 Red","D.SO4 Redox","SOX",
                                                             "rTCA","3-H BC","MethGen","P.A.A.K.","CBB",
                                                             "DeNit","A.NO3 Red","Anammox","PhotoHet"))
clr.All_GOI.all$KO_Function.KEGG = factor(clr.All_GOI.all$KO_Function.KEGG, levels=unique(clr.All_GOI.all$KO_Function.KEGG[order(clr.All_GOI.all$Cycle,clr.All_GOI.all$Pathway)]), ordered=TRUE)

head(clr.All_GOI.all)

All_GOI.hm1a<-ggplot(clr.All_GOI.all, aes(PlotID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.5","-0.5"),breaks=c(1.5,0.5,-0.5)) + labs(title="All Genes of Interest by Pathway & Cycle in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(All_GOI.hm1a,filename = "figures/MGM_Figs/FxnDiv/All_GOI/All_GOI_KOFxns_MGMs_SampID_by_Function_heatmap.png", width=18, height=13, dpi=600)

All_GOI.hm1b<-ggplot(clr.All_GOI.all, aes(PlotID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.5","-0.5"),breaks=c(1.5,0.5,-0.5)) + labs(title="All Genes of Interest by Pathway & Cycle in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~.,scales="free_y", space = "free")

ggsave(All_GOI.hm1b,filename = "figures/MGM_Figs/FxnDiv/All_GOI/All_GOI_KOFxns_MGMs_SampID_by_Function_Pathway_heatmap.png", width=20, height=23, dpi=600)

All_GOI.hm1c<-ggplot(clr.All_GOI.all, aes(PlotID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.5","-0.5"),breaks=c(1.5,0.5,-0.5)) + labs(title="All Genes of Interest by Pathway & Cycle in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Cycle~.,scales="free_y", space = "free")

ggsave(All_GOI.hm1c,filename = "figures/MGM_Figs/FxnDiv/All_GOI/All_GOI_KOFxns_MGMs_SampID_by_Function_Cycle_heatmap.png", width=17, height=15, dpi=600)

All_GOI.hm1d<-ggplot(clr.All_GOI.all, aes(PlotID, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.5","-0.5"),breaks=c(1.5,0.5,-0.5)) + labs(title="All Genes of Interest by Pathway & Cycle in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~SampDate,scales="free_y", space = "free")

ggsave(All_GOI.hm1d,filename = "figures/MGM_Figs/FxnDiv/All_GOI/All_GOI_KOFxns_MGMs_SampID_by_Function_Pathway_SampDate_heatmap.png", width=17, height=15, dpi=600)

All_GOI.hm1a1<-ggplot(clr.All_GOI.all, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.5","-0.5"),breaks=c(1.5,0.5,-0.5)) + labs(title="All Genes of Interest by Pathway & Cycle in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~.,scales="free_y", space = "free")

ggsave(All_GOI.hm1a1,filename = "figures/MGM_Figs/FxnDiv/All_GOI/All_GOI_KOFxns_MGMs_Depth_by_Function_Pathway_heatmap.png", width=18, height=13, dpi=600)

All_GOI.hm1b1<-ggplot(clr.All_GOI.all, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.5","-0.5"),breaks=c(1.5,0.5,-0.5)) + labs(title="All Genes of Interest by Pathway & Cycle in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(Cycle~.,scales="free_y", space = "free")

ggsave(All_GOI.hm1b1,filename = "figures/MGM_Figs/FxnDiv/All_GOI/All_GOI_KOFxns_MGMs_Depth_by_Function_Cycle_heatmap.png", width=20, height=23, dpi=600)

All_GOI.hm1c1<-ggplot(clr.All_GOI.all, aes(Depth_m, KO_Function.KEGG, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.15) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8",labels=c("1.5","0.5","-0.5"),breaks=c(1.5,0.5,-0.5)) + labs(title="All Genes of Interest by Pathway & Cycle in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),strip.text.y = element_text(size = 11,face="bold")) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(PathShort~SampDate,scales="free_y", space = "free")

ggsave(All_GOI.hm1c1,filename = "figures/MGM_Figs/FxnDiv/All_GOI/All_GOI_KOFxns_MGMs_Depth_by_Function_Pathway_SampDate_heatmap.png", width=17, height=15, dpi=600)

#### Pull Out Arsenic Metabolic Fxns from CLR data ####
ars.ko<-mgm.clr.na[,which(colnames(mgm.clr.na) %in% arsen.fxns$KO_ID)]
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
clr.ars.all$PlotID = factor(clr.ars.all$PlotID, levels=unique(clr.ars.all$PlotID[order(clr.ars.all$SampDate,clr.ars.all$Depth_m)]), ordered=TRUE)

# For heatmap color gradient
max(clr.ars.all$CLR_SumCovPerKO, na.rm=TRUE)
max(clr.ars.all$CLR_SumCovPerKO, na.rm=TRUE)/2
min(clr.ars.all$CLR_SumCovPerKO, na.rm=TRUE)

ars.hm1<-ggplot(clr.ars.all, aes(PlotID, KO_Function, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Arsenic Functions in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),
        axis.text = element_text(size=12),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17),
        axis.ticks=element_line(size=0.4),panel.border=element_blank()) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(ars.hm1,filename = "figures/MGM_Figs/FxnDiv/Arsenic/Arsenic_KOFxns_MGMs_heatmap1.png", width=18, height=15, dpi=600)

ars.hm2<-ggplot(clr.ars.all, aes(Depth_m, KO_Function, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Arsenic Functions in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),
        axis.text = element_text(size=12),axis.text.x = element_text(),legend.text = element_text(size=12),plot.title = element_text(size=17),
        axis.ticks=element_line(size=0.4),panel.border=element_blank()) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)

ggsave(ars.hm2,filename = "figures/MGM_Figs/FxnDiv/Arsenic/Arsenic_KOFxns_MGMs_heatmap_better.png", width=18, height=12, dpi=600)

#### Pull Out Selenium Metabolism Fxns from CLR data ####
sel.ko<-mgm.clr.na[,which(colnames(mgm.clr.na) %in% selen.fxns$KO_ID)]
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
clr.sel.all$PlotID = factor(clr.sel.all$PlotID, levels=unique(clr.sel.all$PlotID[order(clr.sel.all$SampDate,clr.sel.all$Depth_m)]), ordered=TRUE)

# For heatmap color gradient
max(clr.sel.all$CLR_SumCovPerKO, na.rm=TRUE)
max(clr.sel.all$CLR_SumCovPerKO, na.rm=TRUE)/2
min(clr.sel.all$CLR_SumCovPerKO, na.rm=TRUE)

sel.hm1<-ggplot(clr.sel.all, aes(PlotID, KO_Function, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Selenium Functions in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),
        axis.text = element_text(size=12),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17),
        axis.ticks=element_line(size=0.4),panel.border=element_blank()) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(sel.hm1,filename = "figures/MGM_Figs/FxnDiv/Selenium/Selenium_KOFxns_MGMs_heatmap1.png", width=18, height=15, dpi=600)

sel.hm2<-ggplot(clr.sel.all, aes(Depth_m, KO_Function, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Selenium Functions in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),
        axis.text = element_text(size=12),axis.text.x = element_text(),legend.text = element_text(size=12),plot.title = element_text(size=17),
        axis.ticks=element_line(size=0.4),panel.border=element_blank()) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)

ggsave(sel.hm2,filename = "figures/MGM_Figs/FxnDiv/Selenium/Selenium_KOFxns_MGMs_heatmap_better.png", width=18, height=12, dpi=600)

#### Pull Out Osmoprotectant/tolerance Fxns from CLR data ####
osmo.ko<-mgm.clr.na[,which(colnames(mgm.clr.na) %in% osmo.fxns$KO_ID)]
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
clr.osmo.all$PlotID = factor(clr.osmo.all$PlotID, levels=unique(clr.osmo.all$PlotID[order(clr.osmo.all$SampDate,clr.osmo.all$Depth_m)]), ordered=TRUE)

# For heatmap color gradient
max(clr.osmo.all$CLR_SumCovPerKO, na.rm=TRUE)
max(clr.osmo.all$CLR_SumCovPerKO, na.rm=TRUE)/2
min(clr.osmo.all$CLR_SumCovPerKO, na.rm=TRUE)

osmo.hm1<-ggplot(clr.osmo.all, aes(PlotID, KO_Function, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Osmoprotectant Functions in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),
        axis.text = element_text(size=12),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17),
        axis.ticks=element_line(size=0.4),panel.border=element_blank()) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(osmo.hm1,filename = "figures/MGM_Figs/FxnDiv/OsmoProctection_Tolerance/OsmoProtectant_KOFxns_MGMs_heatmap1.png", width=18, height=12, dpi=600)

osmo.hm2<-ggplot(clr.osmo.all, aes(Depth_m, KO_Function, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Osmoprotectant & Tolerance Functions in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),
        axis.text = element_text(size=12),axis.text.x = element_text(),legend.text = element_text(size=12),plot.title = element_text(size=17),
        axis.ticks=element_line(size=0.4),panel.border=element_blank()) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)

ggsave(osmo.hm2,filename = "figures/MGM_Figs/FxnDiv/OsmoProctection_Tolerance/OsmoProtectant_KOFxns_MGMs_heatmap_better.png", width=18, height=12, dpi=600)

#### Pull Out Metal Resistance/Tolerance Fxns from CLR data ####
met.ko<-mgm.clr.na[,which(colnames(mgm.clr.na) %in% metal.fxns$KO_ID)]
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
clr.met.all$PlotID = factor(clr.met.all$PlotID, levels=unique(clr.met.all$PlotID[order(clr.met.all$SampDate,clr.met.all$Depth_m)]), ordered=TRUE)

met.hm1<-ggplot(clr.met.all, aes(PlotID, KO_Function, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Metal Resistance & Tolerance Functions in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),
        axis.text = element_text(size=12),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17),
        axis.ticks=element_line(size=0.4),panel.border=element_blank()) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(met.hm1,filename = "figures/MGM_Figs/FxnDiv/MetalResist_Tolerance/Metal_ResistToler_KOFxns_MGMs_heatmap1.png", width=18, height=15, dpi=600)

met.hm2<-ggplot(clr.met.all, aes(Depth_m, KO_Function, fill=CLR_SumCovPerKO)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="#ffaf43", high="#5f03f8") + labs(title="Metal Resistance & Tolerance Functions in Salton Seawater Metagenomes",subtitle="Using CLR-Transformed, Gene Coverage Summed by KO",fill="CLR Coverage Per KO") +
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),
        axis.text = element_text(size=12),axis.text.x = element_text(),legend.text = element_text(size=12),plot.title = element_text(size=17),
        axis.ticks=element_line(size=0.4),panel.border=element_blank()) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+ facet_grid(.~SampDate)

ggsave(met.hm2,filename = "figures/MGM_Figs/FxnDiv/MetalResist_Tolerance/Metal_ResistToler_KOFxns_MGMs_heatmap_better.png", width=18, height=15, dpi=600)

#### Sulfur Metabolism PCoA ####
## PCOA with CLR transformed data first
# calculate our Euclidean distance matrix using CLR data

# pull out sulfur functions from CLR transformed, summed coverages (summed coverage per KO)
sulf.ko<-mgm.clr.na[,which(colnames(mgm.clr.na) %in% sulfur.fxns$KO_ID)] # merge CLR data w/ S fxns found in contigs from KOFamScan
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

ggsave(pcoa.s1,filename = "figures/MGM_Figs/FxnDiv/Sulfur/SSW_SulfurOnly_pcoa_CLR_SummedCoverage_Per_KO_sampdate.png", width=12, height=10, dpi=600)

# sample month shape, depth color
pcoa.s2<-ggplot(sulf.pcoa.clr.meta, aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampleMonth), size=5)+theme_bw()+
  labs(title="PCoA: Sulfur Metabolsim in Salton Seawater",subtitle="Using CLR Transformed, Summed Gene Coverage per KO Function",xlab="PC1", ylab="PC2",color="Depth (m)")+
  theme_classic()+ theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),axis.text = element_text(size=12),axis.text.x = element_text(vjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17))+
  scale_color_continuous(low="blue3",high="red",trans = 'reverse') + scale_shape_discrete(labels=c("August 2021","December 2021","April 2022"),name="Sample Date") +
  xlab("PC1 [67.36%]") + ylab("PC2 [16.73%]")

ggsave(pcoa.s2,filename = "figures/MGM_Figs/FxnDiv/Sulfur/SSW_SulfurOnly_pcoa_CLR_SummedCoverage_Per_KO.traits_depth.png", width=12, height=10, dpi=600)

#### Sulfur Fxns ANOVA ####

clr.sulf.ko.meta<-merge(clr.cov.sum.sulf.ko,meta_scaled,by="SampleID")
clr.sulf.ko.meta[1:4,1:4]
clr.sulf.ko.meta[(nrow(clr.sulf.ko.meta)-4):(nrow(clr.sulf.ko.meta)),(ncol(clr.sulf.ko.meta)-4):(ncol(clr.sulf.ko.meta))] # last 4 rows & cols
rownames(clr.sulf.ko.meta)<-clr.sulf.ko.meta$SampleID

sox.f1<-aov(`soxC; sulfane dehydrogenase subunit SoxC` ~ SampDate, data=clr.sulf.ko.meta)
#pairwise.adonis(bac.div.metadat$Bac_Shannon_Diversity, bac.div.metadat$SampDate, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different

summary(sox.f1)
#Df           Sum Sq Mean Sq    F value   Pr(>F)
#SampDate     2 0.2197 0.10984   6.027 0.0367 *
#Residuals    6 0.1094 0.01822
Tuk1<-TukeyHSD(sox.f1)
Tuk1$SampDate
#                             diff        lwr      upr      p adj
#December.2021-August.2021 -0.31513947 -0.653346242 0.02306729 0.06483634
#April.2022-August.2021     0.03047278 -0.307733990 0.36867955 0.95902845
#April.2022-December.2021   0.34561225  0.007405484 0.68381902 0.04604173

# Levene's test with one independent variable
## Levene's tests whether variances of 2 samples are equal
## we want variances to be the same -- want NON SIGNIFICANCE!
## t test assumes that variances are the same, so Levene's test needs to be non significant
## Fligner's test is a Levene's test for data that are not normally distributed
## more here: https://www.geeksforgeeks.org/fligner-killeen-test-in-r-programming/
fligner.test(`soxC; sulfane dehydrogenase subunit SoxC` ~ SampDate, data = clr.sulf.ko.meta)
# Fligner-Killeen:med chi-squared = 1.1963, df = 7, p-value = 0.991
# Which shows that the data do not deviate significantly from homogeneity.
#compare_means(`soxC; sulfane dehydrogenase subunit SoxC` ~ SampDate, data=clr.sulf.ko.meta, method="anova",p.adjust.method = "bonferroni") # won't take as.factor(Elevation) as input

#### DOM Metabolism PCoA ####
## PCOA with CLR transformed data first
# calculate our Euclidean distance matrix using CLR data


mgm.clr.na[1:4,1:4]

# pull out sulfur functions from CLR transformed, summed coverages (summed coverage per KO)
DOM.ko<-mgm.clr.na[,which(colnames(mgm.clr.na) %in% DOM.fxns$KO_ID)] # merge CLR data w/ DOM-related fxns found in contigs from KOFamScan
DOM.ko$SampleID<-rownames(DOM.ko)
DOM.ko.melt<-melt(DOM.ko, by="SampleID")
colnames(DOM.ko.melt)[which(names(DOM.ko.melt) == "variable")] <- "KO_ID"
colnames(DOM.ko.melt)[which(names(DOM.ko.melt) == "value")] <- "CLR_SumCovPerKO"
head(DOM.ko.melt) #sanity check

clr.DOM.ko<-merge(DOM.ko.melt,all_goi.kegg,by.x=c("KO_ID"),by.y=c("KO_ID"))
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

ggsave(pcoa.s1,filename = "figures/MGM_Figs/FxnDiv/DOM/SSW_DOMOnly_pcoa_CLR_SummedCoverage_Per_KO_sampdate.png", width=12, height=10, dpi=600)

# sample month shape, depth color
pcoa.s2<-ggplot(DOM.pcoa.clr.meta, aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampleMonth), size=5)+theme_bw()+
  labs(title="PCoA: DOM Metabolsim in Salton Seawater",subtitle="Using CLR Transformed, Summed Gene Coverage per KO Function",xlab="PC1", ylab="PC2",color="Depth (m)")+
  theme_classic()+ theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),axis.text = element_text(size=12),axis.text.x = element_text(vjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17))+
  scale_color_continuous(low="blue3",high="red",trans = 'reverse') + scale_shape_discrete(labels=c("August 2021","December 2021","April 2022"),name="Sample Date") +
  xlab("PC1 [63.86%]") + ylab("PC2 [19.14%]")

ggsave(pcoa.s2,filename = "figures/MGM_Figs/FxnDiv/DOM/SSW_DOMOnly_pcoa_CLR_SummedCoverage_Per_KO.traits_depth.png", width=12, height=10, dpi=600)


#### Carbon Fixation PCoA ####
## PCOA with CLR transformed data first
# calculate our Euclidean distance matrix using CLR data


mgm.clr.na[1:4,1:4]

# pull out sulfur functions from CLR transformed, summed coverages (summed coverage per KO)
carb.ko<-mgm.clr.na[,which(colnames(mgm.clr.na) %in% carb.fxns$KO_ID)] # merge CLR data w/ DOM-related fxns found in contigs from KOFamScan
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

ggsave(pcoa.s1,filename = "figures/MGM_Figs/FxnDiv/Carbon/SSW_CarbonFixationOnly_pcoa_CLR_SummedCoverage_Per_KO_sampdate.png", width=12, height=10, dpi=600)

# sample month shape, depth color
pcoa.s2<-ggplot(carb.pcoa.clr.meta, aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampleMonth), size=5)+theme_bw()+
  labs(title="PCoA: DOM Metabolsim in Salton Seawater",subtitle="Using CLR Transformed, Summed Gene Coverage per KO Function",xlab="PC1", ylab="PC2",color="Depth (m)")+
  theme_classic()+ theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),axis.text = element_text(size=12),axis.text.x = element_text(vjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17))+
  scale_color_continuous(low="blue3",high="red",trans = 'reverse') + scale_shape_discrete(labels=c("August 2021","December 2021","April 2022"),name="Sample Date") +
  xlab("PC1 [29.64%]") + ylab("PC2 [22.51%]")

ggsave(pcoa.s2,filename = "figures/MGM_Figs/FxnDiv/Carbon/SSW_CarbonFixationOnly_pcoa_CLR_SummedCoverage_Per_KO.traits_depth.png", width=12, height=10, dpi=600)

#### Arsenic Metabolism PCoA ####
## PCOA with CLR transformed data first
# calculate our Euclidean distance matrix using CLR data
ars.ko<-mgm.clr.na[,which(colnames(mgm.clr.na) %in% arsen.kegg$KO_ID)]
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

ggsave(pcoa.a1,filename = "figures/MGM_Figs/FxnDiv/Arsenic/SSW_ArsenicOnly_pcoa_CLR_SummedCoverage_Per_KO_sampdate.png", width=12, height=10, dpi=600)

# sample month shape, depth color
pcoa.a2<-ggplot(ars.pcoa.clr.meta, aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampleMonth), size=5)+theme_bw()+
  labs(title="PCoA: Arsenic Metabolsim in Salton Seawater",subtitle="Using CLR Transformed, Summed Gene Coverage per KO Function",xlab="PC1", ylab="PC2",color="Depth (m)")+
  theme_classic()+ theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),axis.text = element_text(size=12),axis.text.x = element_text(vjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17))+
  scale_color_continuous(low="blue3",high="red",trans = 'reverse') + scale_shape_discrete(labels=c("August 2021","December 2021","April 2022"),name="Sample Date") +
  xlab("PC1 [53.60%]") + ylab("PC2 [30.40%]")

ggsave(pcoa.a2,filename = "figures/MGM_Figs/FxnDiv/Arsenic/SSW_ArsenicOnly_pcoa_CLR_SummedCoverage_Per_KO.traits_depth.png", width=12, height=10, dpi=600)


#### Metal Resistance/Tolerance PCoA ####
## PCOA with CLR transformed data first
# calculate our Euclidean distance matrix using CLR data
met.ko<-mgm.clr.na[,which(colnames(mgm.clr.na) %in% metal.fxns$KO_ID)]
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

ggsave(pcoa.m1,filename = "figures/MGM_Figs/FxnDiv/MetalResist_Tolerance/SSW_MetalResistTolerance_pcoa_CLR_SummedCoverage_Per_KO_sampdate.png", width=12, height=10, dpi=600)

# sample month shape, depth color
pcoa.m2<-ggplot(met.pcoa.clr.meta, aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampleMonth), size=5)+theme_bw()+
  labs(title="PCoA: Metal Resistance & Tolerance in Salton Seawater",subtitle="Using CLR Transformed, Summed Gene Coverage per KO Function",xlab="PC1", ylab="PC2",color="Depth (m)")+
  theme_classic()+ theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),axis.text = element_text(size=12),axis.text.x = element_text(vjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17))+
  scale_color_continuous(low="blue3",high="red",trans = 'reverse') + scale_shape_discrete(labels=c("August 2021","December 2021","April 2022"),name="Sample Date") +
  xlab("PC1 [61.21%]") + ylab("PC2 [19.15]")

ggsave(pcoa.m2,filename = "figures/MGM_Figs/FxnDiv/MetalResist_Tolerance/SSW_MetalResistTolerance_pcoa_CLR_SummedCoverage_Per_KO.traits_depth.png", width=12, height=10, dpi=600)

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
png('figures/MGM_Figs/FxnDiv/SSW_MGM_pcoa_CLR_SummedCoverage_perKO_betadispersion_sampledate.png',width = 700, height = 600, res=100)
plot(mgm.disper,main = "Centroids and Dispersion based on Aitchison Distance (CLR Data)", col=colorset1$SampDate_Color)
dev.off()

png('figures/MGM_Figs/FxnDiv/SSW_MGM_boxplot_CLR_SummedCoverage_perKO_centroid_distance_sampledate.png',width = 700, height = 600, res=100)
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
png('figures/MGM_Figs/FxnDiv/ssw_mgm_pcoa_CLR_SummedCoverage_per_KO_betadispersion_depth.png',width = 700, height = 600, res=100)
plot(mgm.disper2,main = "Centroids and Dispersion based on Aitchison Distance (CLR Data)", col=colfunc(3))
dev.off()

png('figures/MGM_Figs/FxnDiv/ssw_mgm_boxplot_CLR_SummedCoverage_per_KO_centroid_distance_depth.png',width = 700, height = 600, res=100)
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
png('figures/MGM_Figs/FxnDiv/SSW_MGM_pcoa_CLR_SummedCoverage_perKO_SulfurOnly_betadispersion_sampledate.png',width = 700, height = 600, res=100)
plot(sulf.disper1,main = "Centroids and Dispersion based on Aitchison Distance (Sulfur CLR Data)", col=colorset1$SampDate_Color)
dev.off()

png('figures/MGM_Figs/FxnDiv/SSW_MGM_boxplot_CLR_SummedCoverage_perKO_SulfurOnly_centroid_distance_sampledate.png',width = 700, height = 600, res=100)
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
png('figures/MGM_Figs/FxnDiv/ssw_mgm_pcoa_CLR_SummedCoverage_per_KO_SulfOnly_betadispersion_depth.png',width = 700, height = 600, res=100)
plot(sulf.disper2,main = "Centroids and Dispersion based on Aitchison Distance (Sulfur CLR Data)", col=colfunc(3))
dev.off()

png('figures/MGM_Figs/FxnDiv/ssw_mgm_boxplot_CLR_SummedCoverage_per_KO_SulfOnly_centroid_distance_depth.png',width = 700, height = 600, res=100)
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

#### Sulfur Fxns PERMANOVA ####

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

# Create perm variable
perm <- with(meta_scaled, how(nperm = 1000))

# Run PERMANOVAs
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

#### Sulfur Pathway PERMANOVA ####
## what about genes by S metabolic pathway?

# check rownames to make sure sub-dfs of S genes are in same order as meta_scaled
rownames(sox.ko.cov) %in% rownames(meta_scaled)
rownames(asSO4.ko.cov) %in% rownames(meta_scaled)

# includes CLR transformed coverage (summed up per gene per KO) for genes in specific Sulfur Metabolic pathways

# SOX genes
spath0<-adonis2(sox.ko.cov ~ SampDate,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
spath0
#           Df SumOfSqs     R2     F   Pr(>F)
# SampDate  2  0.73160 0.7062 7.211 0.003996 **
# Residual  6  0.30437 0.2938
# Total     8  1.03596 1.0000
p.adjust(spath0$`Pr(>F)`,method="bonferroni",n=3) # adjusted pval = 0.01198801

spath0a<-adonis2(sox.ko.cov ~ SampDate*Depth.num,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
spath0a

sox.euc.dist.mr <- dist(sox.ko.cov, method = "euclidean")
sox.disper1<-betadisper(sox.euc.dist.mr, meta_scaled$SampDate)
sox.disper1

permutest(sox.disper1, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
#Pairwise comparisons:
#  (Observed p-value below diagonal, permuted p-value above diagonal)
#               August.2021 December.2021 April.2022
#August.2021                     0.48200      0.596
#December.2021     0.43936                    0.495
#April.2022        0.55487       0.47369

anova(sox.disper1) # p = 0.5879 --> accept the Null H, spatial medians are NOT significantly difference across sample dates

TukeyHSD(sox.disper1) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other

#                             diff        lwr       upr     p adj
#December.2021-August.2021 -0.10864260 -0.4286886 0.2114034 0.5803442
#April.2022-August.2021    -0.07943587 -0.3994818 0.2406101 0.7382982
#April.2022-December.2021   0.02920673 -0.2908392 0.3492527 0.9579997

# Visualize dispersions
png('figures/MGM_Figs/FxnDiv/SSW_SOX_PCoA_betadispersion_sampledate.png',width = 700, height = 600, res=100)
plot(sox.disper1,main = "SOX Functions - Centroids and Dispersion (Median-Ratio Data)", label=FALSE,col=colorset1$SampDate_Color)
dev.off()

png('figures/MGM_Figs/FxnDiv/SSW_SOX_boxplot_MR_centroid_distance_sampledate.png',width = 700, height = 600, res=100)
boxplot(sox.disper1,xlab="Sample Collection Date", main = "SOX Functions - Distance to Centroid by Category (Median-Ratio Data)", sub="Euclidean Distance of Median-Ratio Transformed Data", col=colorset1$SampDate_Color)
dev.off()

# Assimilatory Sulfate Reduction
spath1<-adonis2(asSO4.ko.cov ~ SampDate,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
spath1
#           Df SumOfSqs      R2      F Pr(>F)
# SampDate  2  0.40909 0.28398 1.1898 0.3506
# Residual  6  1.03147 0.71602
# Total     8  1.44056 1.00000
p.adjust(spath1$`Pr(>F)`,method="bonferroni",n=3)

spath1a<-adonis2(asSO4.ko.cov ~ SampDate*Depth.num,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
spath1a

# Dissimilatory Sulfate Reduction
spath2<-adonis2(disSO4.ko.cov ~ SampDate,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
spath2
#           Df SumOfSqs      R2      F   Pr(>F)
# SampDate  2   3.0096 0.98128 157.29 0.004995 **
# Residual  6   0.0574 0.01872
# Total     8   3.0670 1.00000
p.adjust(spath2$`Pr(>F)`,method="bonferroni",n=3) # adjusted pval= 0.01498501

spath2a<-adonis2(disSO4.ko.cov ~ SampDate*Depth.num,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
spath2a

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
# using PCoA axes as a variable to see if can be predicted by env vars
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