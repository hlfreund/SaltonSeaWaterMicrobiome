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
# NOTE: chem data has already been scaled in meta_scaled; raw chem/env data in metadata
# all samples included (outliers from alpha div data were not removed)

bac.dat.all[1:6,1:6]
bac.ASV_table[1:4,1:4]
head(metadata) # not scaled metadata
head(meta_scaled) # centered & scaled metadata
head(meta_scaled_S)

meta_scaled_S<-na.omit(meta_scaled_S)

#### Using Shapiro-Wilk test for Normality ####

shapiro.test(meta_scaled$DO_Percent_Local) # what is the p-value?
# my p-value was p-value =  3.406e-08
# p > 0.05 states distribution of data are not significantly different from normal distribution
# p < 0.05 means that data is significantly different from a normal distribution
hist(meta_scaled$DO_Percent_Local, col="blue")

# visualize Q-Q plot for species richness
qqnorm(meta_scaled$DO_Percent_Local, pch = 1, frame = FALSE)
qqline(meta_scaled$DO_Percent_Local, col = "steelblue", lwd = 2)

shapiro.test(meta_scaled$ORP_mV) # what is the p-value? p-value = 3.277e-11
hist(meta_scaled$ORP_mV, col="blue")
# visualize Q-Q plot for species richness
qqnorm(meta_scaled$ORP_mV, pch = 1, frame = FALSE)
qqline(meta_scaled$ORP_mV, col = "steelblue", lwd = 2)

shapiro.test(meta_scaled$Temp_DegC) # what is the p-value? p-value = 4.186e-05
hist(meta_scaled$Temp_DegC, col="blue")
# visualize Q-Q plot for species richness
qqnorm(meta_scaled$Temp_DegC, pch = 1, frame = FALSE)
qqline(meta_scaled$Temp_DegC, col = "steelblue", lwd = 2)

shapiro.test(meta_scaled$Dissolved_OrganicMatter_RFU) # what is the p-value? p-value = 1.997e-07
hist(meta_scaled$Dissolved_OrganicMatter_RFU, col="blue")
# visualize Q-Q plot for species richness
qqnorm(meta_scaled$Dissolved_OrganicMatter_RFU, pch = 1, frame = FALSE)
qqline(meta_scaled$Dissolved_OrganicMatter_RFU, col = "steelblue", lwd = 2)

shapiro.test(meta_scaled_S$Sulfate_milliM) # what is the p-value?
# my p-value was p-value =  0.006965
# p > 0.05 states distribution of data are not significantly different from normal distribution
# p < 0.05 means that data is significantly different from a normal distribution
hist(meta_scaled_S$Sulfate_milliM, col="blue")

# visualize Q-Q plot for species richness
qqnorm(meta_scaled_S$Sulfate_milliM, pch = 1, frame = FALSE)
qqline(meta_scaled_S$Sulfate_milliM, col = "steelblue", lwd = 2)

shapiro.test(meta_scaled_S$Sulfide_microM) # what is the p-value?
# my p-value was p-value =  5.934e-12
# p > 0.05 states distribution of data are not significantly different from normal distribution
# p < 0.05 means that data is significantly different from a normal distribution
hist(meta_scaled_S$Sulfide_microM, col="blue")

# visualize Q-Q plot for species richness
qqnorm(meta_scaled_S$Sulfide_microM, pch = 1, frame = FALSE)
qqline(meta_scaled_S$Sulfide_microM, col = "steelblue", lwd = 2)

shapiro.test(meta_scaled$Turbidity_FNU) # what is the p-value?  p-value = 0.0005629
hist(meta_scaled$Turbidity_FNU, col="blue")
# visualize Q-Q plot for species richness
qqnorm(meta_scaled$Turbidity_FNU, pch = 1, frame = FALSE)
qqline(meta_scaled$Turbidity_FNU, col = "steelblue", lwd = 2)

shapiro.test(meta_scaled$Chlorophyll_RFU) # what is the p-value? p-value = 1.044e-11
hist(meta_scaled$Chlorophyll_RFU, col="blue")
# visualize Q-Q plot for species richness
qqnorm(meta_scaled$Chlorophyll_RFU, pch = 1, frame = FALSE)
qqline(meta_scaled$Chlorophyll_RFU, col = "steelblue", lwd = 2)

#### PCA w/ Env Variables ####
env.dat<-metadata[,c(8,10:11,15)]

# NOTE: PCA requires normally distributed data, so we are log transforming env data before scaling
# this is because scaled env data alone is not normally distributed
env.log<-decostand(env.dat,method = "log", pseudocount = 1) # log transformation of env data before scaling
env.log[1:4,1:4]

# check rownames of log transformed environmental data & metadata
rownames(env.log) %in% rownames(meta_scaled)
meta_scaled=meta_scaled[rownames(env.log),] ## reorder metadata to match order of log data

# calculate our Euclidean distance matrix using log data
env.euc_dist <- dist(env.log, method = "euclidean")

# creating our hierarcical clustering dendrogram
env.euc_clust <- hclust(env.euc_dist, method="ward.D2")

# let's make it a little nicer...
env.euc_dend <- as.dendrogram(env.euc_clust, hang=0.2)
env.dend_cols <- as.character(metadata$SampDate_Color[order.dendrogram(env.euc_dend)])
labels_colors(env.euc_dend) <- env.dend_cols

plot(env.euc_dend, ylab="Log Euclidean Distance",cex = 0.5) + title(main = "Environmental Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 3)
#legend("topright",legend = c("June 2021","August 2021","December 2021","April 2022"),cex=.8,col = c( "#26547c","#36ab57","#32cbff","#ff6f00"),pch = 15, bty = "n")
# Control is dark blue ("#218380"), #Alternaria is light blue ("#73d2de")
dev.off()

# let's use our Euclidean distance matrix from before
env.pca <- prcomp(env.log, center=TRUE, scale=TRUE) # pca of euclidean distance matrix = PCA of euclidean distance matrix
env.pca$x # where sites fall on PC axes
env.pca$rotation # variables on PC axes
summary(env.pca)$importance
# The proportion of variances

# extract principal coordinates
env.pca.vectors<-data.frame(env.pca$x)
env.pca.vectors$SampleID<-rownames(env.pca$x)

# merge pca coordinates w/ metadata
env.pca.meta<-merge(env.pca.vectors, meta_scaled, by.x="SampleID", by.y="SampleID")
env.pca.meta$SampleMonth
env.pca.meta$SampDate

head(env.pca.meta)
summary(env.pca)$importance # percentage of variation explained for pca below

# create pca ggplot fig
pca1<-ggplot(env.pca.meta, aes(x=PC1, y=PC2)) +geom_point(aes(color=factor(SampDate)), size=2)+theme_bw()+
  labs(title="PCA:DO%, DOM, & ORP in Salton Seawater",subtitle="Using Log Transformed, Scaled Data",color="Sample Type")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Type",values=unique(env.pca.meta$SampDate_Color[order(env.pca.meta$SampDate)]),labels=c("June.2021"="June 2021","August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Axis 1 [62.91%]") + ylab("Axis 2 [24.10%]")

ggsave(pca1,filename = "figures/EnvVariablesOnly/SSW_env_pca_log_sampdate.png", width=12, height=10, dpi=600)

# sample month shape, depth color
pca2<-ggplot(env.pca.meta, aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampleMonth), size=5)+theme_bw()+
  labs(title="pca: Bacteria/Archaea in Salton Seawater",subtitle="Using Centered-Log Ratio Data",xlab="Axis 1", ylab="Axis 2",color="Depth (m)")+
  theme_classic()+ theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),axis.text = element_text(size=12),axis.text.x = element_text(vjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17))+
  scale_color_continuous(low="blue3",high="red",trans = 'reverse') + scale_shape_discrete(labels=c("June 2021","August 2021","December 2021","April 2022"),name="Sample Date") +
  xlab("Axis 1 [22.72%]") + ylab("Axis 2 [17.97%]")

ggsave(pca2,filename = "figures/EnvVariablesOnly/SSW_env_pca_log_depth.png", width=12, height=10, dpi=600)

### what if we include HS/SO4 data?
metadata_S[,c(8,10:11,15:17)] # what env data are we considering?
env.dat.S<-metadata_S[,c(8,10:11,15:17)]
env.dat.S<-na.omit(env.dat.S)

# NOTE: PCA requires normally distributed data, so we are log transforming env data before scaling
# this is because scaled env data alone is not normally distributed
env.log.S<-decostand(env.dat.S,method = "log", pseudocount = 1) # log transformation of env data before scaling
env.log.S[1:4,1:4]

# check rownames of log transformed environmental data & metadata
rownames(env.log.S) %in% rownames(meta_scaled_S)
meta_scaled_S=meta_scaled_S[rownames(env.log.S),] ## reorder metadata to match order of log data

# calculate our Euclidean distance matrix using log data
env.S.euc_dist <- dist(env.log.S, method = "euclidean")

# creating our hierarcical clustering dendrogram
env.S.euc_clust <- hclust(env.S.euc_dist, method="ward.D2")

# let's make it a little nicer...
env.S.euc_dend <- as.dendrogram(env.S.euc_clust, hang=0.2)
env.S.dend_cols <- as.character(metadata_S$SampDate_Color[order.dendrogram(env.S.euc_dend)])
labels_colors(env.S.euc_dend) <- env.S.dend_cols

plot(env.S.euc_dend, ylab="Log Euclidean Distance",cex = 0.5) + title(main = "Environmental Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 3)
#legend("topright",legend = c("June 2021","August 2021","December 2021","April 2022"),cex=.8,col = c( "#26547c","#36ab57","#32cbff","#ff6f00"),pch = 15, bty = "n")
# Control is dark blue ("#218380"), #Alternaria is light blue ("#73d2de")
dev.off()

# let's use our Euclidean distance matrix from before
env.S.pca <- prcomp(env.log.S, center=TRUE, scale=TRUE) # pca of euclidean distance matrix = PCA of euclidean distance matrix
env.S.pca$x # where sites fall on PC axes
env.S.pca$rotation # variables on PC axes
summary(env.S.pca)$importance
# The proportion of variances

# extract principal coordinates
env.S.pca.vectors<-data.frame(env.S.pca$x)
env.S.pca.vectors$SampleID<-rownames(env.S.pca$x)

# merge pca coordinates w/ metadata
env.S.pca.meta<-merge(env.S.pca.vectors, meta_scaled_S, by.x="SampleID", by.y="SampleID")
env.S.pca.meta$SampleMonth
env.S.pca.meta$SampDate

head(env.S.pca.meta)
summary(env.S.pca)$importance # percentage of variation explained for pca below

# create pca ggplot fig
s.pca1<-ggplot(env.S.pca.meta, aes(x=PC1, y=PC2)) +geom_point(aes(color=factor(SampDate)), size=2)+theme_bw()+
  labs(title="PCA:DO%, DOM, ORP, HS, SO4, & Temp in Salton Seawater",subtitle="Using Log Transformed, Scaled Data",color="Sample Type")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Type",values=unique(env.S.pca.meta$SampDate_Color[order(env.S.pca.meta$SampDate)]),labels=c("June.2021"="June 2021","August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Axis 1 [62.87%]") + ylab("Axis 2 [20.64%]")

ggsave(s.pca1,filename = "figures/EnvVariablesOnly/SSW_env_w.S_pca_log_sampdate.png", width=12, height=10, dpi=600)

# sample month shape, depth color
s.pca2<-ggplot(env.S.pca.meta, aes(x=PC1, y=PC2)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampleMonth), size=5)+theme_bw()+
  labs(title="PCA:DO%, DOM, ORP, HS, SO4, & Temp in Salton Seawater",subtitle="Using Log Transformed, Scaled Data",xlab="PC1", ylab="PC2",color="Depth (m)")+
  theme_classic()+ theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),axis.text = element_text(size=12),axis.text.x = element_text(vjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17))+
  scale_color_continuous(low="blue3",high="red",trans = 'reverse') + scale_shape_discrete(labels=c("August 2021","December 2021","April 2022"),name="Sample Date") +
  xlab("Axis 1 [62.87%]") + ylab("Axis 2 [20.64%]")

ggsave(s.pca2,filename = "figures/EnvVariablesOnly/SSW_env_w.S_pca_log_depth.png", width=12, height=10, dpi=600)


#### Compare Env Samples Variance by Sample Date ####
kruskal.test(DO_Percent_Local ~ SampleMonth, data = meta_scaled)
pairwise.wilcox.test(meta_scaled$DO_Percent_Local, meta_scaled$SampleMonth, p.adjust.method = "bonf") # returns p values
#           June    August  December   ##August x December, and December x April are significantly different
#August   1.00000 -       -
#December 0.68634 0.00057 -
#April    0.68634 0.01408 0.00102

kruskal.test(ORP_mV ~ SampleMonth, data = meta_scaled)
pairwise.wilcox.test(meta_scaled$ORP_mV, meta_scaled$SampleMonth, p.adjust.method = "bonf") # returns p values
#           June    August  December
#August   0.42636 -       -
#December 0.00116 0.00192 -
#April    0.00109 0.00054 8.3e-06

kruskal.test(Temp_DegC ~ SampleMonth, data = meta_scaled)
pairwise.wilcox.test(meta_scaled$Temp_DegC, meta_scaled$SampleMonth, p.adjust.method = "bonf") # returns p values
#             June    August  December
#August   0.00817 -       -
#December 0.00116 0.00057 -
#April    0.00116 0.00057 8.6e-06

kruskal.test(Dissolved_OrganicMatter_RFU ~ SampleMonth, data = meta_scaled)
pairwise.wilcox.test(meta_scaled$Dissolved_OrganicMatter_RFU, meta_scaled$SampleMonth, p.adjust.method = "bonf") # returns p values
#         June    August  December
#August   0.87282 -       -
#December 0.68200 0.00055 -
#April    0.68634 0.00057 8.5e-06

kruskal.test(Turbidity_FNU ~ SampleMonth, data = meta_scaled)
pairwise.wilcox.test(meta_scaled$Turbidity_FNU, meta_scaled$SampleMonth, p.adjust.method = "bonf") # returns p values
#           June    August  December
#August   1.00000 -       -
#December 1.00000 0.00055 -
#April    0.67329 0.00052 8e-06

kruskal.test(Chlorophyll_RFU ~ SampleMonth, data = meta_scaled)
pairwise.wilcox.test(meta_scaled$Chlorophyll_RFU, meta_scaled$SampleMonth, p.adjust.method = "bonf") # returns p values
#           June    August  December
#August   0.02597 -       -
#December 1.00000 0.00057 -
#April    0.00530 0.17466 8.6e-06

#### Do Env Variables Correlate? ####
# check for colinearity among env variables themselves
heatmap(abs(cor(meta_scaled[,c(8,10:14)])),
        # Compute pearson correlation (note they are absolute values)
        col = rev(heat.colors(6)),
        Colv = NA, Rowv = NA)
legend("topleft",
       title = "Absolute Pearson R",
       legend =  round(seq(0,1, length.out = 6),1),
       y.intersp = 0.7, bty = "n",
       fill = rev(heat.colors(6)))
dev.off()

# excluding HS & SO4
cor_mat.env1 <- cor(meta_scaled[,c(8,10:15)], method='pearson')
cor_mat.env1

symnum(cor_mat.env1)

corrplot.mixed(cor_mat.env1, tl.pos='lt', tl.cex=0.6, number.cex=0.5, addCoefasPercent=T)
# env variables with a correlation of <|0.7| is a good threshold for determining if predictors correlate

# including HS & SO4 (and not June 2021 data)

cor_mat.env2 <- cor(meta_scaled_S[,c(8,10:17)], method='pearson')
cor_mat.env2

symnum(cor_mat.env2)

corrplot.mixed(cor_mat.env2, tl.pos='lt', tl.cex=0.6, number.cex=0.5, addCoefasPercent=T)
# env variables with a correlation of <|0.7| is a good threshold for determining if predictors correlate

# DO %
cor.test(meta_scaled$DO_Percent_Local, meta_scaled$ORP_mV, method="pearson") # ***
# r = 0.4019, p = 0.005 --> not a strong correlation & it's significant
cor.test(meta_scaled$DO_Percent_Local, meta_scaled$Temp_DegC, method="pearson")
cor.test(meta_scaled$DO_Percent_Local, meta_scaled$Chlorophyll_RFU, method="pearson") # ****
# r = 0.8008, p = 1.397e-11 --> strong correlation & significant
cor.test(meta_scaled$DO_Percent_Local, meta_scaled$Dissolved_OrganicMatter_RFU, method="pearson") # ****
# r = -0.8907, p < 2.2e-16 --> strong correlation & significant

plot(x=meta_scaled$DO_Percent_Local, y=meta_scaled$Dissolved_OrganicMatter_RFU)
summary(lm(DO_Percent_Local ~ Dissolved_OrganicMatter_RFU, data=meta_scaled) %>%
  adjust_pvalue(method="bonferroni"))

# ORP
cor.test(meta_scaled$ORP_mV, meta_scaled$Temp_DegC, method="pearson") # ***
# r = -0.4499, p-value = 0.0015 --> mediocre correlation & significant
cor.test(meta_scaled$ORP_mV, meta_scaled$Chlorophyll_RFU, method="pearson")
cor.test(meta_scaled$ORP_mV, meta_scaled$Dissolved_OrganicMatter_RFU, method="pearson") # ***
# r = -0.380499, p-value = 0.00833 --> not a strong correlation & it's significant

# Chlorophyll
cor.test(meta_scaled$Chlorophyll_RFU, meta_scaled$Dissolved_OrganicMatter_RFU, method="pearson") # ****
# r = -0.70582, p-value = 3.002e-08 --> strong correlation & significant
cor.test(meta_scaled$Chlorophyll_RFU, meta_scaled$Temp_DegC, method="pearson")

# Dissolved Organic Matter
cor.test(meta_scaled$Dissolved_OrganicMatter_RFU, meta_scaled$Temp_DegC, method="pearson")

# Sulfate (milliM)
cor.test(meta_scaled_S$Sulfate_milliM, meta_scaled_S$ORP_mV, method="pearson")
# r = 0.1102182, p = 0.4984 --> not a strong correlation & it's significant
cor.test(meta_scaled_S$Sulfate_milliM, meta_scaled_S$Temp_DegC, method="pearson") # ***
# r = -0.5614189, p = 0.0001639 --> strong negative correlation, significant
cor.test(meta_scaled_S$Sulfate_milliM, meta_scaled_S$DO_Percent_Local, method="pearson") # ****
# r = 0.6452846, p = 6.944e-06 --> strong correlation & significant
cor.test(meta_scaled_S$Sulfate_milliM, meta_scaled_S$Dissolved_OrganicMatter_RFU, method="pearson")
# r = -0.03455734 , p = 0.8323 --> no corr, not sig
cor.test(meta_scaled_S$Sulfate_milliM, meta_scaled_S$Sulfide_microM, method="pearson")
# r = -0.1571489 , p = 0.3328 --> no corr, not sig

# Sulfide (microM)
cor.test(meta_scaled_S$Sulfide_microM, meta_scaled_S$ORP_mV, method="pearson") # ******
# r = -0.979198, p < 2.2e-16 --> STRONG & significant negative correlation
cor.test(meta_scaled_S$Sulfide_microM, meta_scaled_S$Temp_DegC, method="pearson") # ***
# r = 0.5532523 , p = 0.0002134 --> medium correlation, significant
cor.test(meta_scaled_S$Sulfide_microM, meta_scaled_S$DO_Percent_Local, method="pearson") # ****
# r = -0.6286855, p = 1.398e-05 --> medium-strong negative correlation & significant
cor.test(meta_scaled_S$Sulfide_microM, meta_scaled_S$Dissolved_OrganicMatter_RFU, method="pearson") # ****
# r = 0.621356 , p = 1.88e-05 --> medium to strong correlation, significant

#### Do Env Data Vary Significantly By Group?#####

## sampling date

## sampling depth
#### Does Env Predict Alpha Diversity? ####
## linear regression time!

## First let's make the alpha diversity data frame
## Calculate Shannon Diversity (abundance + richness considered in diversity calculation)
# if you have another package loaded that has a diversity function, you can specify that you want to use vegan's diversity function as shown below
Shan_ent.16s<-vegan::diversity(bac.ASV_table[,-1], index="shannon") # Shannon entropy
Shan_div.16s<- exp(Shan_ent.16s) # Shannon Diversity aka Hill number 1

# create data frame with Shannon entropy and Shannon diversity values
div_16s<-data.frame(Bac_Shannon_Entropy=Shan_ent.16s,Bac_Shannon_Diversity=Shan_div.16s)
class(div_16s)
div_16s$SampleID<-rownames(div_16s)
head(div_16s)

# Calculate species richness (number of species per sample)
specnumber(bac.ASV_table[,-1])

# Create a DF with Species Richness
S_16s<-data.frame(Bac_Species_Richness=specnumber(bac.ASV_table[,-1]), SampleID=rownames(bac.ASV_table)) # finds # of species per sample using RAW count data; if MARGIN = 2 it finds frequencies of species

# merge richness and diversity dataframes together
d.r_16s<-merge(div_16s, S_16s, by.x="SampleID", by.y="SampleID")

# merge w/ metadata
bac.div.metadat <- merge(d.r_16s,meta_scaled, by.x="SampleID", by.y="SampleID")
head(bac.div.metadat)
class(bac.div.metadat) # want data frame

unique(bac.div.metadat$SampleMonth) # see how many elements there are in the Group variable
unique(bac.div.metadat$Depth_m) # see how many elements there are in the Group variable

# drop the outliers
bac.div.metadat2<-subset(bac.div.metadat, bac.div.metadat$Bac_Shannon_Diversity<=200)

# Linear Regression time!
## here the focus is comparing dust complexity to alpha diversity, species richness, & elevation
head(bac.div.metadat2)
s.div.lm.fit1<-lm(Bac_Shannon_Diversity ~ DO_Percent_Local, data=bac.div.metadat2) %>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(s.div.lm.fit1)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)        97.679      6.548  14.916   <2e-16 ***
#  DO_Percent_Local   -1.922      6.647  -0.289    0.774

s.div.lm.fit2<-lm(Bac_Shannon_Diversity ~ ORP_mV, data=bac.div.metadat2) %>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(s.div.lm.fit2)

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)   97.674      3.709  26.333  < 2e-16 ***
#  ORP_mV       -15.535      3.720  -4.176 0.000138 ***
## ^^^ the two lms below show that this model is significant only for June & August 2021, not December & April

not_summer_months<-subset(bac.div.metadat2, SampDate=="December.2021" | SampDate=="April.2022" )

s.div.lm.fit2a<-lm(Bac_Shannon_Diversity ~ ORP_mV, data=not_summer_months) %>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(s.div.lm.fit2a)

summer_months<-subset(bac.div.metadat2, SampDate=="June.2021" | SampDate=="August.2021" )

s.div.lm.fit2b<-lm(Bac_Shannon_Diversity ~ ORP_mV, data=summer_months) %>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(s.div.lm.fit2b)

s.div.lm.fit3<-lm(Bac_Shannon_Diversity ~ Temp_DegC, data=bac.div.metadat2) %>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(s.div.lm.fit3)

s.div.lm.fit5<-lm(Bac_Shannon_Diversity ~ Dissolved_OrganicMatter_RFU, data=bac.div.metadat2) %>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(s.div.lm.fit5)




fit1<-aov(Bac_Shannon_Diversity ~ Depth_m, data=bac.div.metadat2)
pairwise.adonis(bac.div.metadat2$Bac_Shannon_Diversity, bac.div.metadat2$Depth_m, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different

summary(fit1)
#Df           Sum Sq Mean Sq    F value   Pr(>F)
#Elevation2   3 0.7277 0.24258   0.084 0.774
#Residuals   27 0.4444 0.01646
Tuk1<-TukeyHSD(fit1)
Tuk1$Depth_m

#plot(DustComplexity ~ Elevation, data=bac.div.metadat2)
#abline(aov(DustComplexity ~ Elevation, data=bac.div.metadat2))

# fit.0<-aov(DustComplexity ~ as.factor(Elevation), data=bac.div.metadat2)
# summary(fit.0)
# TukeyHSD(fit.0)
# Levene's test with one independent variable
## Levene's tests whether variances of 2 samples are equal
## we want variances to be the same -- want NON SIGNIFICANCE!
## t test assumes that variances are the same, so Levene's test needs to be non significant
fligner.test(Bac_Shannon_Diversity ~ Depth_m, data = bac.div.metadat2)
# Levenes Test for Homogeneity of Variance
#        Df  Chi square value  Pr(>F)
# group  3   1.0952   0.7411
# Which shows that the data do not deviate significantly from homogeneity.
elev<-bac.div.metadat2$Elevation2
compare_means(DustComplexity ~ Elevation2, data=bac.div.metadat2, method="anova",p.adjust.method = "bonferroni") # won't take as.factor(Elevation) as input

p.adj.dc.elev<-compare_means(DustComplexity ~ Elevation2, data=bac.div.metadat2, method="t.test",p.adjust.method = "bonferroni") # won't take as.factor(Elevation) as input# Note https://github.com/kassambara/ggpubr/issues/65

fit.test<-ggplot(bac.div.metadat2, aes(x = as.factor(Elevation), y = DustComplexity, fill=as.factor(Elevation))) +
  geom_boxplot() + theme_classic() + guides(fill = guide_legend(reverse=TRUE)) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  labs(title="Dust Complexity x Elevation",fill="Elevation (ft)")+ylab("Dust Complexity")+xlab("Elevation (ft)")+scale_fill_manual(values=saturation(fair_cols, 0.9))+stat_compare_means(method = "anova",label.y=1.5) +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4)), method="t.test", hide.ns = TRUE,label = "p.signif")

ggsave(fit.test,filename = "figures/EnvVariablesOnly/DustComp_by_Elevation_ALL_sigbars_5.24.21.pdf", width=10, height=8, dpi=600)

fit.testa<-ggplot(bac.div.metadat2, aes(x = as.factor(Elevation), y = DustComplexity, fill=as.factor(Elevation))) +
  geom_boxplot() + theme_classic() + guides(fill = guide_legend(reverse=TRUE)) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  labs(title="Dust Complexity x Elevation",fill="Elevation (ft)")+ylab("Dust Complexity")+xlab("Elevation (ft)")+scale_fill_manual(values=saturation(fair_cols, 0.9))+stat_compare_means(method = "anova",label.y=1.5,mapping=aes(label = format.pval(..p.adj.., digits = 3))) +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4)), method="t.test", hide.ns = TRUE,mapping=aes(label = format.pval(..p.adj.., digits = 3)))

ggsave(fit.testa,filename = "figures/EnvVariablesOnly/DustComp_by_Elevation_ALL_sigbars_5.24.21.pdf", width=10, height=8, dpi=600)

fit.test0<-ggplot(bac.div.metadat2, aes(x = as.factor(Elevation), y = DustComplexity, fill=as.factor(Elevation))) +
  geom_boxplot() + theme_classic() + guides(fill = guide_legend(reverse=TRUE)) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  labs(title="Dust Complexity x Elevation",fill="Elevation (ft)")+ylab("Dust Complexity")+xlab("Elevation (ft)")+scale_fill_manual(values=saturation(fair_cols, 0.9))+stat_compare_means(method = "anova",label.y=1.5) +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4)), method="t.test", hide.ns = TRUE,label = "p.signif")

ggsave(fit.test,filename = "figures/EnvVariablesOnly/DustComp_by_Elevation_ALL_sigbars_5.24.21.pdf", width=10, height=8, dpi=600)

fit.testa<-ggplot(bac.div.metadat2, aes(x = as.factor(Elevation), y = DustComplexity, fill=as.factor(Elevation))) +
  geom_boxplot() + theme_classic() + guides(fill = guide_legend(reverse=TRUE)) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  labs(title="Dust Complexity x Elevation",fill="Elevation (ft)")+ylab("Dust Complexity")+xlab("Elevation (ft)")+scale_fill_manual(values=saturation(fair_cols, 0.9))

ggsave(fit.testa,filename = "figures/EnvVariablesOnly/DustComp_by_Elevation_ALL_no.sigbars_5.24.21.pdf", width=10, height=8, dpi=600)

fit.testb<-ggplot(bac.div.metadat2, aes(x = as.factor(Elevation), y = DustComplexity, fill=as.factor(Elevation))) +
  geom_boxplot() + theme_classic() + guides(fill = guide_legend(reverse=TRUE)) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  labs(title="Dust Complexity x Elevation",fill="Elevation (ft)")+ylab("Dust Complexity")+xlab("Elevation (ft)")+scale_fill_grey(start=0.8, end=0.3)

ggsave(fit.testb,filename = "figures/EnvVariablesOnly/DustComp_by_Elevation_ALL_gray_5.24.21.pdf", width=10, height=8, dpi=600)

fit.testb.0<-ggplot(bac.div.metadat2, aes(x = as.factor(Elevation), y = DustComplexity, fill=as.factor(Elevation))) +
  geom_boxplot() + theme_classic() + guides(fill = guide_legend(reverse=TRUE)) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  labs(title="Dust Complexity x Elevation",fill="Elevation (ft)")+ylab("Dust Complexity")+xlab("Elevation (ft)")+scale_fill_grey(start=0.8, end=0.3)+stat_compare_means(method = "anova",label.y=1.5) +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4)), method="t.test", hide.ns = TRUE,label = "p.signif")

ggsave(fit.testb.0,filename = "figures/EnvVariablesOnly/DustComp_by_Elevation_ALL_gray_sigbars_5.24.21.pdf", width=10, height=8, dpi=600)

### Fungi comparisons first
# Dust Comp x ITS1 Shannon diversity
hist(bac.div.metadat2$ITS1_Shannon_Diversity) # NOT normally distributed
hist(bac.div.metadat2$DustComplexity) # somewhat normally distributed
chisq.test(bac.div.metadat2$ITS1_Shannon_Diversity, bac.div.metadat2$DustComplexity)

its1.fit1<-lm(DustComplexity ~ ITS1_Shannon_Diversity, data=bac.div.metadat2) %>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
summary(its1.fit1)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)             0.614082   0.050858  12.074  7.8e-13 ***
#  ITS1_Shannon_Diversity -0.001942   0.001416  -1.372    0.181

its1.fit1.p<-p.adjust(coef(summary(its1.fit1))[8], method="bonferroni") # pvalue

plot(DustComplexity ~ ITS1_Shannon_Diversity, data=bac.div.metadat2)
abline(its1.fit1)


#leveneTest(bac.div.metadat2$DustComplexity,
#            bac.div.metadat2$ITS1_Shannon_Diversity,
#            location = c("median"),
#            trim.alpha = 0.25)
# Levenes Test for Homogeneity of Variance
#        Df  F value  Pr(>F)
# group  3   2.3415   0.0818
# Which shows that the data do not deviate significantly from homogeneity.

fig.its1.fit1<-ggplot(its1_div_meta, aes(x = ITS1_Shannon_Diversity, y = DustComplexity)) +
  geom_point(aes(color=Elev.num), size=3) + theme_classic() + saturation(scale_colour_gradientn(colours=fair_cols,limits=c(400,2700),breaks = c(500,1250,2000,2600),labels=c("400","1100","2000","2700")), 0.9) +
  stat_smooth(method = "lm", col = "black", se=FALSE, size=1)+ labs(color="Elevation (ft)")+ylab("Dust Complexity")+xlab("ITS1 Shannon Diversity")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  stat_cor(label.y = 1, label.x=75) +
  stat_regline_equation(aes(label=paste(..adj.rr.label..)),label.y = 1.20,label.x=75)
## use summary(its1.fit1) to double check that stat_cor gives same p value as linear regression!

ggsave(fig.its1.fit1,filename = "figures/EnvVariablesOnly/DustComp_by_ITS1_ShanDiv_ALL_1.4.22.pdf", width=10, height=8, dpi=600)

fig.its1.fit1<-ggplot(its1_div_meta, aes(x = ITS1_Shannon_Diversity, y = DustComplexity)) +
  geom_point(aes(color=Elev.num),size=3) + theme_classic() + saturation(scale_colour_gradientn(colours=fair_cols,limits=c(400,2700),breaks = c(500,1250,2000,2600),labels=c("400","1100","2000","2700")), 0.9) +
  stat_smooth(method = "lm", col = "black", se=FALSE, size=1)+ labs(title="Dust Complexity x ITS1 Shannon Diversity", color="Elevation (ft)")+ylab("Dust Complexity")+xlab("ITS1 Shannon Diversity")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  stat_cor(label.y = 1, label.x=75) +
  stat_regline_equation(aes(label=paste(..adj.rr.label..)),label.y = 1.05,label.x=75)
## use summary(its1.fit1) to double check that stat_cor gives same p value as linear regression!

#fig.its1.fit2<-ggplot(bac.div.metadat2, aes(x = ITS1_Shannon_Diversity, y = DustComplexity)) +
#  geom_point(aes(color=Elevation), size=3) + theme_classic() + saturation(scale_colour_gradientn(colours=fair_cols,limits=c(400,2700),breaks = c(500,1250,2000,2600),labels=c("400","1100","2000","2700")), 0.9) +
#  stat_smooth(method = "glm", col = "black", se=FALSE, size=1)+ labs(title="Dust Complexity x ITS1 Shannon Diversity", color="Elevation (ft)")+ylab("Dust Complexity")+xlab("ITS1 Shannon Diversity")+
#  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
#  stat_cor(label.y = 1, label.x=75) +
#  stat_regline_equation(label.y = 1.05,label.x=75)

#ggsave(fig.its1.fit2,filename = "figures/EnvVariablesOnly/DustComp_by_ITS1_Shan_Div_ALL_5.19.21.pdf", width=10, height=8, dpi=600)


# DustComp x ITS1 Species Richness
hist(bac.div.metadat2$ITS1_Species_Richness)
its1.sr.fit1<-lm(DustComplexity ~ ITS1_Species_Richness, data=bac.div.metadat2) %>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
summary(its1.sr.fit1)
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)
# (Intercept)            0.7761611  0.0757511  10.246 3.79e-11 ***
#  ITS1_Species_Richness -0.0004551  0.0001475  -3.084  0.00445 **
coef(summary(its1.sr.fit1)) # pvalue

its1.sr.fit2<-lm(DustComplexity ~ ITS1_Species_Richness, data=bac.div.metadat2) %>%
  adjust_pvalue(method="bonferroni")
summary(its1.sr.fit2)
coef(summary(its1.sr.fit2))
p.adjust(coef(summary(its1.sr.fit2))[,4], method="bonferroni") # pvalue

its1.sr.fit2<-glm(DustComplexity ~ ITS1_Species_Richness, data=bac.div.metadat2, family=poisson)
its1.sr.fit3<-glm.nb(DustComplexity ~ ITS1_Species_Richness, data=bac.div.metadat2)

summary(its1.sr.fit2)
dispersiontest(its1.sr.fit2)
# null hypothesis is that equidispersion exists; alternative hypothesis is overdispersion
# if overdispersion, use negative binomial not Poisson
## Poisson distribution implies that the mean and variance are equal --> little dispersion
# negative binomial means # of observations is not fixed, whereas binomial means observations are a fixed #

# z = -16.609, p-value = 1 (cannot reject null)
# alternative hypothesis: true dispersion is greater than 1
# sample estimates:
#   dispersion
# 0.0495281 -- equidispersion exists

plot(DustComplexity ~ ITS1_Species_Richness, data=bac.div.metadat2)
abline(its1.fit2)

fig.its1.sr.fit1<-ggplot(its1_div_meta, aes(x = ITS1_Species_Richness, y = DustComplexity)) +
  geom_point(aes(color=fair_cols), size=3) + theme_classic() + saturation(scale_colour_gradientn(colours=fair_cols,limits=c(400,2700),breaks = c(500,1250,2000,2600),labels=c("400","1100","2000","2700")), 0.9) +
  stat_smooth(method = "lm", col = "black", se=FALSE, size=1)+ labs(title="Dust Complexity x ITS1 Species Richness", color="Elevation (ft)")+ylab("Dust Complexity")+xlab("ITS1 Species Richness")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  stat_cor(label.y = 1, label.x=700) +
  stat_regline_equation(aes(label=paste(..adj.rr.label..)),label.y = 1.05,label.x=700)

## use summary(its1.sr.fit1) to double check that stat_cor gives same p value as linear regression!

ggsave(fig.its1.sr.fit1,filename = "figures/EnvVariablesOnly/DustComp_by_ITS1_Spec_Richness_ALL_1.4.22.pdf", width=10, height=8, dpi=600)


