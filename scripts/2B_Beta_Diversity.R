#### Set WD & Load Libraries ####
getwd() # use setwd("path/to/files") if you are not in the right directory
setwd("/Volumes/HLF_SSD/Aronson_Lab_Data/Salton_Sea/SaltonSeaWater")
suppressPackageStartupMessages({ # load packages quietly
  library(devtools)
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
  library(ggbiplot)
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
  library(decontam)
  library(ggvegan)
  library(microbiome)
  library(pairwiseAdonis)
  library(corrplot)
})

#### Load Global Env to Import Count/ASV Tables ####
load("data/SSeawater_Data_Ready.Rdata") # save global env to Rdata file
#load("data/ssw_clr.euc.dist.Rdata")
#save.image("data/SSW_env.seq_analysis.Rdata")

#save.image("data/Env_Seqs_All/env.seq_analysis.Rdata") # save global env to Rdata file
bac.dat.all[1:6,1:6]
bac.ASV_table[1:4,1:4]
bac.ASV_table[(nrow(bac.ASV_table)-4):(nrow(bac.ASV_table)),(ncol(bac.ASV_table)-4):(ncol(bac.ASV_table))] # last 4 rows & cols
head(meta_scaled)

#### Beta Diversity ####
rownames(bac.ASV_table)
bac.ASV_table[1:4,1:4]

# CLR transformation of ASV table
# df must have rownames are SampleIDs, columns are ASV IDs for vegan functions below
b.clr<-decostand(bac.ASV_table[,-1],method = "clr", pseudocount = 1) #CLR transformation
b.clr[1:4,1:4]

# check rownames of CLR transformed ASV data & metadata
rownames(b.clr) %in% rownames(meta_scaled)
meta_scaled=meta_scaled[rownames(b.clr),] ## reorder metadata to match order of CLR data

# calculate our Euclidean distance matrix using CLR data
b.euc_dist <- dist(b.clr, method = "euclidean")

# creating our hierarcical clustering dendrogram
b.euc_clust <- hclust(b.euc_dist, method="ward.D2")

# let's make it a little nicer...
b.euc_dend <- as.dendrogram(b.euc_clust, hang=0.2)
b.dend_cols <- as.character(meta_scaled$SampDate_Color[order.dendrogram(b.euc_dend)])
labels_colors(b.euc_dend) <- b.dend_cols

## DO NOT RUN THIS LINE, THIS IS YOUR COLOR REFERENCE!!!!
(August.2021="#ef781c",December.2021="#03045e",April.2022="#059c3f")

plot(b.euc_dend, ylab="CLR Euclidean Distance",cex = 0.5) + title(main = "Bacteria/Archaea Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 2)
#legend("topright",legend = c("August 2021","December 2021","April 2022"),cex=.8,col = c( "#36ab57","#32cbff","#ff6f00"),pch = 15, bty = "n")
# Control is dark blue ("#218380"), #Alternaria is light blue ("#73d2de")
dev.off()

# PCOA w/ Euclidean distance matrix (of CLR data)
b.pcoa <- pcoa(b.euc_dist) # pcoa of euclidean distance matrix = PCA of euclidean distance matrix
save.image("data/ssw_clr.euc.dist.Rdata")

# The proportion of variances explained is in its element values$Relative_eig
b.pcoa$values

# extract principal coordinates
b.pcoa.vectors<-data.frame(b.pcoa$vectors)
b.pcoa.vectors$SampleID<-rownames(b.pcoa$vectors)

# merge pcoa coordinates w/ metadata
b.pcoa.meta<-merge(b.pcoa.vectors, meta_scaled, by.x="SampleID", by.y="SampleID")
b.pcoa.meta$SampleMonth
b.pcoa.meta$SampDate

head(b.pcoa.meta)

head(b.pcoa$values) # pull out Relative (Relative_eig) variation % to add to axes labels

# create PCoA ggplot fig
pcoa1<-ggplot(b.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(SampDate)), size=4)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Seawater",subtitle="Using Centered-Log Ratio Data",color="Sample Date")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Type",values=unique(b.pcoa.meta$SampDate_Color[order(b.pcoa.meta$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("PC1 [28.88%]") + ylab("PC2 [19.60%]")

ggsave(pcoa1,filename = "figures/BetaDiversity/SSW_16S_pcoa_CLR_sampdate.png", width=12, height=10, dpi=600)

# sample month shape, depth color
pcoa2<-ggplot(b.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(color=as.numeric(as.character(Depth_m)),shape=SampleMonth), size=5)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Seawater",subtitle="Using Centered-Log Ratio Data",color="Depth (m)")+
  theme_classic()+ theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),axis.text = element_text(size=12),axis.text.x = element_text(vjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17))+
  scale_color_continuous(low="blue3",high="red",trans = 'reverse') + scale_shape_discrete(labels=c("August 2021","December 2021","April 2022"),name="Sample Date") +
  xlab("PC1 [28.88%]") + ylab("PC2 [19.60%]")

ggsave(pcoa2,filename = "figures/BetaDiversity/SSW_16S_pcoa_CLR_depth_sampdate.png", width=12, height=10, dpi=600)

#### Homogeneity of Variance & PERMANOVA tests - Composition by Groups ####
## betadisper to look at homogeneity of group dispersions (aka variance) when considering multiple variables
# multivariate analogue to Levene's test of homogeneity of variances
# program finds spatial median or centroid of the group, & compare distances of group to centroid/spatial median via ANOVA

rownames(metadata) %in% rownames(b.clr) #b.clr was used to make the distance matrix b.euc_dist

# first by compare dispersions by sampling date
b.disper1<-betadisper((vegdist(b.clr,method="euclidean")), meta_scaled$SampDate)
b.disper1

permutest(b.disper1, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
#Pairwise comparisons:
#  (Observed p-value below diagonal, permuted p-value above diagonal)
#                            August.2021 December.2021 April.2022
#August.2021                     0.66200      0.191
#December.2021     0.57192                    0.157
#April.2022        0.19575       0.14302

anova(b.disper1) # p = 0.2673 --> accept the Null H, spatial medians (a measure of dispersion) are NOT significantly difference across sample dates

TukeyHSD(b.disper1) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other
# timepoints are not significantly different from each other considering ALL ASVs

# If PERMANOVA is significant but betadisper() IS NOT, then you can infer that there is only a location effect.
# If both tests are significant, then there is a dispersion effect for sure and there might also be (not always) a location effect.
# Dispersion effect means the actual spread of the data points is influencing the significant differences, not the actual data itself

pnova1<-adonis2(b.clr ~ SampDate,data=meta_scaled,method = "euclidean",by="terms",permutations=1000)
pnova1 # p-value = 0.000999

b.clr.dist = (vegdist(b.clr, "euclidean", na.rm = TRUE)) #distance matrix using Bray's dissimilarity index for trait distribution (traits of interest only)
pair.mod1<-pairwise.adonis(b.clr.dist,meta_scaled$SampDate, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different
pair.mod1
#                          pairs Df SumsOfSqs   F.Model        R2 p.value p.adjusted sig
#1  December.2021 vs April.2022  1  9501.431 12.160599 0.2884352   0.001      0.003   *
#2 December.2021 vs August.2021  1  8569.637  9.114434 0.2929327   0.001      0.003   *
#3    April.2022 vs August.2021  1 10311.629 18.675966 0.4591401   0.001      0.003   *

# Visualize dispersions
png('figures/BetaDiversity/pcoa_betadispersion_sampledate.png',width = 700, height = 600, res=100)
plot(b.disper1,main = "Centroids and Dispersion based on Aitchison Distance", col=colorset1$SampDate_Color)
dev.off()

png('boxplot_centroid_distance_sampledate.png',width = 700, height = 600, res=100)
boxplot(b.disper1,xlab="Sample Collection Date", main = "Distance to Centroid by Category", sub="Based on Aitchison Distance", col=colorset1$SampDate_Color)
dev.off()

# Next compare dispersions by depth
b.disper2<-betadisper((vegdist(b.clr,method="euclidean")), meta_scaled$Depth_m)
b.disper2

permutest(b.disper2, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons

anova(b.disper2) # p = 0.7632 --> accept the Null H, spatial medians are NOT significantly difference across sample dates

TukeyHSD(b.disper2) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other
# no sig results

# If PERMANOVA is significant but betadisper() IS NOT, then you can infer that there is only a location effect.
# If both tests are significant, then there is a dispersion effect for sure and there might also be (not always) a location effect.
# Dispersion effect means the actual spread of the data points is influencing the significant differences, not the actual data itself

pnova2<-adonis2(b.clr ~ Depth_m,data=meta_scaled,method = "euclidean",by="terms",permutations=1000)
pnova2 # p-value = 1

#b.clr.dist = (vegdist(b.clr, "euclidean", na.rm = TRUE)) #distance matrix using Bray's dissimilarity index for trait distribution (traits of interest only)
pair.mod2<-pairwise.adonis(b.clr.dist,meta_scaled$Depth_m, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different
pair.mod2
# none are significantly different

col.depth <- colorRampPalette(c("red", "blue"))
col.depth(9)

# Visualize dispersions
png('figures/BetaDiversity/pcoa_betadispersion_depth.png',width = 700, height = 600, res=100)
plot(b.disper2,main = "Centroids and Dispersion based on Aitchison Distance", col=col.depth(9))
dev.off()

png('figures/BetaDiversity/boxplot_centroid_distance_depth.png',width = 700, height = 600, res=100)
boxplot(b.disper2,xlab="Sample Collection Depth", main = "Distance to Centroid by Category", sub="Based on Aitchison Distance", col=col.depth(9))
dev.off()

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
rownames(b.clr) %in% rownames(meta_scaled)
meta_scaled=meta_scaled[rownames(b.clr),] ## reorder metadata to match order of CLR data
perm <- with(meta_scaled, how(nperm = 1000, blocks = SampDate))

pnova1<-adonis2(b.clr ~ DO_Percent_Local*ORP_mV*Temp_DegC*Dissolved_OrganicMatter_RFU*Depth_m*Sulfate_milliM*Sulfide_microM,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
## none are significant

adonis2(b.clr ~ DO_Percent_Local*ORP_mV*Temp_DegC*Dissolved_OrganicMatter_RFU*Depth_m*Sulfate_milliM*Sulfide_microM,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm)
#         Df SumOfSqs     R2    F Pr(>F)
#Model    23    34412 0.73114 1.8918 0.4825
#Residual 16    12654 0.26886
#Total    39    47066 1.00000

# remove categorical variables
pnova2<-adonis2(b.clr ~ DO_Percent_Local*ORP_mV*Temp_DegC*Dissolved_OrganicMatter_RFU*Sulfate_milliM*Sulfide_microM,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
#                                   Df SumOfSqs      R2       F   Pr(>F)
#DO_Percent_Local                   1     2269 0.03771  2.9545 0.036963 *
#ORP_mV                             1     7763 0.12901 10.1061 0.000999 ***
#Temp_DegC                          1     6991 0.11618  9.1010 0.075924 .
#DO_Percent_Local:ORP_mV            1     1023 0.01701  1.3321 0.739261
#DO_Percent_Local:Temp_DegC         1     8964 0.14897 11.6699 0.442557
#ORP_mV:Temp_DegC                   1     2379 0.03954  3.0972 0.007992 **
#DO_Percent_Local:ORP_mV:Temp_DegC  1      827 0.01375  1.0768 0.819181
#Residual                          39    29956 0.49784
#Total                             46    60172 1.00000
adonis2(b.clr ~ DO_Percent_Local*ORP_mV*Temp_DegC*Dissolved_OrganicMatter_RFU*Sulfate_milliM*Sulfide_microM,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm)
#         Df SumOfSqs      R2      F   Pr(>F)
#Model     7    30216 0.50216 5.6197 0.003996 **
#Residual 39    29956 0.49784
#Total    46    60172 1.00000

pnova3<-adonis2(b.clr ~ DO_Percent_Local*ORP_mV*Temp_DegC*Dissolved_OrganicMatter_RFU,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
#                                   Df SumOfSqs      R2       F   Pr(>F)
#ORP_mV                                                         1     2473 0.05255  3.4075 0.02997 *
#Dissolved_OrganicMatter_RFU                                    1     8035 0.17073 11.0705 0.06394 .
#DO_Percent_Local:Dissolved_OrganicMatter_RFU                   1     2471 0.05251  3.4046 0.09191 .
#DO_Percent_Local:ORP_mV:Temp_DegC                              1     1169 0.02483  1.6101 0.02897 *

adonis2(b.clr ~ DO_Percent_Local*ORP_mV*Temp_DegC*Dissolved_OrganicMatter_RFU,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm)
#         Df SumOfSqs      R2      F   Pr(>F)
#Model     15    29646 0.62988 2.7229 0.1518
#Residual 24    17420 0.37012
#Total    39    47066 1.00000

pnova4<-adonis2(b.clr ~ Sulfate_milliM*Sulfide_microM,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
#                                   Df SumOfSqs      R2       F   Pr(>F)
#Sulfate_milliM                 1     6242 0.13262 6.2631 0.37363
#Sulfide_microM                 1     4193 0.08908 4.2069 0.05495 .
#Sulfate_milliM:Sulfide_microM  1      753 0.01600 0.7554 0.61738
#Residual                      36    35878 0.76230
#Total                         39    47066 1.00000

adonis2(b.clr ~ Sulfate_milliM*Sulfide_microM,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm)
#         Df SumOfSqs      R2      F   Pr(>F)
#Model    3    11187 0.2377 3.7418 0.1678
#Residual 36    35878 0.7623
#Total    39    47066 1.0000

pnova5<-adonis2(b.clr ~ ORP_mV*Temp_DegC*Dissolved_OrganicMatter_RFU,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
#                                   Df SumOfSqs      R2       F   Pr(>F)
#ORP_mV                                        1     4239 0.09006 5.7372 0.004995 **
#Temp_DegC                                     1     6087 0.12933 8.2389 0.106893
#Dissolved_OrganicMatter_RFU                   1     5450 0.11579 7.3765 0.063936 .
#ORP_mV:Temp_DegC                              1     3492 0.07419 4.7263 0.180819
#ORP_mV:Dissolved_OrganicMatter_RFU            1     1226 0.02605 1.6596 0.299700
#Temp_DegC:Dissolved_OrganicMatter_RFU         1     1521 0.03231 2.0584 0.039960 *
#  ORP_mV:Temp_DegC:Dissolved_OrganicMatter_RFU  1     1409 0.02994 1.9075 0.057942 .
#Residual                                     32    23642 0.50232
#Total                                        39    47066 1.00000

adonis2(b.clr ~ ORP_mV*Temp_DegC*Dissolved_OrganicMatter_RFU,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm)
#         Df SumOfSqs      R2      F   Pr(>F)
#Model    7    23424 0.49768 4.5292 0.02597 *
#Residual 32    23642 0.50232
#Total    39    47066 1.00000

pnova6a<-adonis2(b.clr ~ ORP_mV*Dissolved_OrganicMatter_RFU,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
adonis2(b.clr ~ ORP_mV*Dissolved_OrganicMatter_RFU,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm) #significant

pnova6b<-adonis2(b.clr ~ ORP_mV*Temp_DegC,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
adonis2(b.clr ~ ORP_mV*Temp_DegC,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm) # insignificant

## BEST MODEL as of 5/11/23:
adonis2(b.clr ~ ORP_mV*Temp_DegC*Dissolved_OrganicMatter_RFU,data=meta_scaled,method = "euclidean",by=NULL,permutations=perm)

#While PERMANOVA tests differences in group means (analogous to MANOVA),
## a related test called PERMDISP can be used to evaluate homogeneity of group dispersion
#(analogous to Levene's test for equal variances). The vegan function for this test is “betadisper”:
## * need a distance matrix!
b.clr.dist = (vegdist(b.clr, "euclidean", na.rm = TRUE)) #distance matrix using Bray's dissimilarity index for trait distribution (traits of interest only)
bac.disper <- betadisper(b.clr.dist, meta_scaled$SampDate)
bac.disper

## Significant differences in homogeneities can be tested using either parametric or permutational tests,
##and parametric post hoc contrasts can also be investigated:
anova(bac.disper)
permutest(bac.disper)
TukeyHSD(bac.disper)

##one issue with adonis is that it doesn't do multiple comparisons *******
# tells us that something is different, but what is different? Which sample/plot/location?
## our four provinces differ, but do all of them differ,or just one?

##random person on the internet to the rescue!
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
meta_scaled=meta_scaled[rownames(b.clr),] ## reorder metadata to match order of CLR data

pair.mod<-pairwise.adonis(b.clr.dist,meta_scaled$Depth_m, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different
pair.mod
pairwise.adonis2(b.clr.dist~Depth_m, data=meta_scaled,strata='SampDate')

pair.mod1<-pairwise.adonis(b.clr.dist,meta_scaled$SampDate, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different
pair.mod1

### SELF REMINDER FOR R^2
### Coefficient of Determination, denoted R2 or r2 and pronounced "R squared"
### is the proportion of the variance in the dependent variable that is predictable from the independent variable(s)

### Pseudo F stat for PERMANOVA
### pseudo F-ratio: It compares the total sum of squared dissimilarities (or ranked dissimilarities) among objects belonging to different groups to that of objects belonging to the same group.
### Larger F-ratios indicate more pronounced group separation, however, the significance of this ratio is usually of more interest than its magnitude.

