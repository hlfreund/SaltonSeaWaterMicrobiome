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
load("data/SSeawater_AlphaBetaDiv_Data.Rdata")
#load("data/ssw_clr.euc.dist_2.21.23.Rdata")

#save.image("data/Env_Seqs_All/env.seq_analysis.Rdata") # save global env to Rdata file
bac.dat.all[1:6,1:6]
bac.ASV_table[1:4,1:4]
bac.ASV_table[(nrow(bac.ASV_table)-4):(nrow(bac.ASV_table)),(ncol(bac.ASV_table)-4):(ncol(bac.ASV_table))] # last 4 rows & cols
head(meta_scaled)

#### Rarefaction & Species Accumulation Curves ####
# bacteria/archaea

# Species Accumulation Curve
sc2<-specaccum(bac.ASV_table[,-1],"random")
plot(sc2, ci.type="poly", col="darkgreen", lwd=2, ci.lty=0, ci.col="lightgreen")
boxplot(sc2, col="yellow", add=TRUE, pch=20)

# Prep for Rarefaction Curve
rowSums(bac.ASV_table[,-1]) # total # ASVs per sample, excluding SampleID from calculation
sort(colSums(bac.ASV_table[,-1]))

# Create Rarefaction curve
png('figures/SSW_16S_rarecurve.png')
rarecurve(as.matrix(bac.ASV_table[,-1]),col=metadata$SampDate_Color, step=1000, label=F,ylab="ASVs")
# to show sampel labels per curve, change label=T
dev.off()

# ASVs per Sample
total_asvs<-data.frame(ASV_Total=rowSums(bac.ASV_table[,-1]),metadata)

ggplot(data=total_asvs, aes(x=SampleID, y=ASV_Total,fill=Sample_Type)) +
  geom_bar(stat="identity",colour="black")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

## average ASV per sample month & depth
aggregate(bac.ASV_all$Count, list(bac.ASV_all$SampleMonth), FUN=mean)
aggregate(bac.ASV_all$Count, list(bac.ASV_all$Depth_m), FUN=mean)

#### Alpha Diversity & Species Richness ####

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

## Shannon Diversity by Sample Month & Depth
bac.a.div<-ggplot(bac.div.metadat2, aes(x=SampDate, y=Bac_Shannon_Diversity,fill=SampDate)) +geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values=unique(bac.div.metadat$SampDate_Color[order(bac.div.metadat$SampDate)]), name ="Sample Date")+theme_classic()+
  labs(title = "Bacterial Shannon Diversity by Sample Date", x="Sample Date", y="Shannon Diversity", fill="Sample Month")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1,size=10),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4)), method="t.test", hide.ns = TRUE,label = "p.format")

ggsave(bac.a.div,filename = "figures/SSW_Bacterial_alpha_diversity_samplemonth_2.1.23.png", width=13, height=10, dpi=600)

bac.div.metadat2$Depth_m=as.numeric(levels(bac.div.metadat2$Depth_m))[bac.div.metadat2$Depth_m]
# note: cannot turn numbers that are factors in R into numeric values...
## have to convert factor levels into numeric, then use the numeric "levels" to pull out numbers from Depth_m column in df to make sure the Depth_m columns is now numeric, not a factor

bac.a.div2<-ggplot(bac.div.metadat2, aes(x=Depth_m, y=Bac_Shannon_Diversity)) +geom_boxplot(aes(fill=bac.div.metadat2$Depth_m),color="black")+
  labs(title = "Bacterial Shannon Diversity by Sampling Depth", x="Depth (m)", y="Shannon Diversity", fill="Depth (m)")+
  scale_fill_gradient(low="red",high="blue",guide = guide_colourbar(reverse = TRUE)) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1,,size=10),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))

ggsave(bac.a.div2,filename = "figures/SSW_Bacterial_alpha_diversity_depth_v1_2.1.23.png", width=13, height=10, dpi=600)

bac.a.div3<-ggplot(bac.div.metadat2, aes(x=as.factor(Depth_m), y=Bac_Shannon_Diversity)) +geom_boxplot(aes(fill=bac.div.metadat2$Depth_m),color="black")+
  labs(title = "Bacterial Shannon Diversity by Sampling Depth", x="Depth (m)", y="Shannon Diversity", fill="Depth (m)")+
  scale_fill_gradient(low="red",high="blue",guide = guide_colourbar(reverse = TRUE)) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1,,size=10),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  coord_flip() + scale_x_discrete(limits=rev)

ggsave(bac.a.div3,filename = "figures/SSW_Bacterial_alpha_diversity_depth_v2_2.1.23.png", width=13, height=10, dpi=600)

## Species Richness by Sample Type
bac.a.sr<-ggplot(bac.div.metadat2, aes(x=SampDate, y=Bac_Species_Richness,fill=SampDate)) +geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values=unique(bac.div.metadat$SampDate_Color[order(bac.div.metadat$SampDate)]), name ="Sample Date")+theme_classic()+
  labs(title = "Bacterial Species Richness by Sample Date", x="Sample Date", y="Species Richness")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1,size=10),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4)), method="t.test", hide.ns = TRUE,label = "p.format")
ggsave(bac.a.sr,filename = "figures/SSW_Bacterial_species_richness_samplemonth_2.1.23.png", width=13, height=10, dpi=600)

bac.a.sr2<-ggplot(bac.div.metadat2, aes(x=as.factor(Depth_m), y=Bac_Species_Richness,fill=bac.div.metadat2$Depth_m)) +geom_boxplot(aes(fill=as.numeric(bac.div.metadat2$Depth_m)),color="black")+
  labs(title = "Bacterial Species Richness by Sampling Depth", x="Depth (m)", y="Species Richness", fill="Depth (m)")+
  scale_fill_gradient(low="red",high="blue",guide = guide_colourbar(reverse = TRUE)) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1,,size=10),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))

ggsave(bac.a.sr2,filename = "figures/SSW_Bacterial_species_richness_depth_v1_2.1.23.png", width=13, height=10, dpi=600)

bac.a.sr3<-ggplot(bac.div.metadat2, aes(x=as.factor(Depth_m), y=Bac_Species_Richness,fill=bac.div.metadat2$Depth_m)) +geom_boxplot(aes(fill=as.numeric(bac.div.metadat2$Depth_m)),color="black")+
  labs(title = "Bacterial Species Richness by Sampling Depth", x="Depth (m)", y="Species Richness", fill="Depth (m)")+
  scale_fill_gradient(low="red",high="blue",guide = guide_colourbar(reverse = TRUE)) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1,,size=10),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  coord_flip() + scale_x_discrete(limits=rev)

ggsave(bac.a.sr3,filename = "figures/SSW_Bacterial_species_richness_depth_v2_2.1.23.png", width=13, height=10, dpi=600)

#### Using Shapiro-Wilk test for Normality ####
shapiro.test(bac.div.metadat2$Bac_Shannon_Diversity) # what is the p-value?
# my p-value was p-value = 2.875e-13
# p > 0.05 states distribution of data are not significantly different from normal distribution
# p < 0.05 means that data is significantly different from a normal distribution
hist(bac.div.metadat2$Bac_Shannon_Diversity, col="blue")

# visualize Q-Q plot for species richness
qqnorm(bac.div.metadat2$Bac_Shannon_Diversity, pch = 1, frame = FALSE)
qqline(bac.div.metadat2$Bac_Shannon_Diversity, col = "steelblue", lwd = 2)

shapiro.test(bac.div.metadat2$Bac_Species_Richness) # what is the p-value?
# my p-value was p-value =  0.01702
# p > 0.05 states distribution of data are not significantly different from normal distribution
# p < 0.05 means that data is significantly different from a normal distribution
hist(bac.div.metadat$Bac_Species_Richness, col="blue")

# visualize Q-Q plot for species richness
qqnorm(bac.div.metadat$Bac_Species_Richness, pch = 1, frame = FALSE)
qqline(bac.div.metadat$Bac_Species_Richness, col = "steelblue", lwd = 2)

shapiro.test(bac.div.metadat2$DO_Percent_Local) # what is the p-value?
hist(bac.div.metadat$DO_Percent_Local, col="blue")

shapiro.test(bac.div.metadat2$ORP_mV) # what is the p-value?
hist(bac.div.metadat$ORP_mV, col="blue")

shapiro.test(bac.div.metadat2$Temp_DegC) # what is the p-value?
hist(bac.div.metadat$Temp_DegC, col="blue")

shapiro.test(bac.div.metadat2$Dissolved_OrganicMatter_RFU) # what is the p-value?
hist(bac.div.metadat$Dissolved_OrganicMatter_RFU, col="blue")

shapiro.test(bac.div.metadat2$Turbidity_FNU) # what is the p-value?
hist(bac.div.metadat$Turbidity_FNU, col="blue")

shapiro.test(bac.div.metadat2$Chlorophyll_RFU) # what is the p-value?
hist(bac.div.metadat$Chlorophyll_RFU, col="blue")

#### Linear Regression Comparisons ####
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

ggsave(fit.test,filename = "figures/DustComp_by_Elevation_ALL_sigbars_5.24.21.pdf", width=10, height=8, dpi=600)

fit.testa<-ggplot(bac.div.metadat2, aes(x = as.factor(Elevation), y = DustComplexity, fill=as.factor(Elevation))) +
  geom_boxplot() + theme_classic() + guides(fill = guide_legend(reverse=TRUE)) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  labs(title="Dust Complexity x Elevation",fill="Elevation (ft)")+ylab("Dust Complexity")+xlab("Elevation (ft)")+scale_fill_manual(values=saturation(fair_cols, 0.9))+stat_compare_means(method = "anova",label.y=1.5,mapping=aes(label = format.pval(..p.adj.., digits = 3))) +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4)), method="t.test", hide.ns = TRUE,mapping=aes(label = format.pval(..p.adj.., digits = 3)))

ggsave(fit.testa,filename = "figures/DustComp_by_Elevation_ALL_sigbars_5.24.21.pdf", width=10, height=8, dpi=600)

fit.test0<-ggplot(bac.div.metadat2, aes(x = as.factor(Elevation), y = DustComplexity, fill=as.factor(Elevation))) +
  geom_boxplot() + theme_classic() + guides(fill = guide_legend(reverse=TRUE)) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  labs(title="Dust Complexity x Elevation",fill="Elevation (ft)")+ylab("Dust Complexity")+xlab("Elevation (ft)")+scale_fill_manual(values=saturation(fair_cols, 0.9))+stat_compare_means(method = "anova",label.y=1.5) +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4)), method="t.test", hide.ns = TRUE,label = "p.signif")

ggsave(fit.test,filename = "figures/DustComp_by_Elevation_ALL_sigbars_5.24.21.pdf", width=10, height=8, dpi=600)

fit.testa<-ggplot(bac.div.metadat2, aes(x = as.factor(Elevation), y = DustComplexity, fill=as.factor(Elevation))) +
  geom_boxplot() + theme_classic() + guides(fill = guide_legend(reverse=TRUE)) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  labs(title="Dust Complexity x Elevation",fill="Elevation (ft)")+ylab("Dust Complexity")+xlab("Elevation (ft)")+scale_fill_manual(values=saturation(fair_cols, 0.9))

ggsave(fit.testa,filename = "figures/DustComp_by_Elevation_ALL_no.sigbars_5.24.21.pdf", width=10, height=8, dpi=600)

fit.testb<-ggplot(bac.div.metadat2, aes(x = as.factor(Elevation), y = DustComplexity, fill=as.factor(Elevation))) +
  geom_boxplot() + theme_classic() + guides(fill = guide_legend(reverse=TRUE)) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  labs(title="Dust Complexity x Elevation",fill="Elevation (ft)")+ylab("Dust Complexity")+xlab("Elevation (ft)")+scale_fill_grey(start=0.8, end=0.3)

ggsave(fit.testb,filename = "figures/DustComp_by_Elevation_ALL_gray_5.24.21.pdf", width=10, height=8, dpi=600)

fit.testb.0<-ggplot(bac.div.metadat2, aes(x = as.factor(Elevation), y = DustComplexity, fill=as.factor(Elevation))) +
  geom_boxplot() + theme_classic() + guides(fill = guide_legend(reverse=TRUE)) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  labs(title="Dust Complexity x Elevation",fill="Elevation (ft)")+ylab("Dust Complexity")+xlab("Elevation (ft)")+scale_fill_grey(start=0.8, end=0.3)+stat_compare_means(method = "anova",label.y=1.5) +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4)), method="t.test", hide.ns = TRUE,label = "p.signif")

ggsave(fit.testb.0,filename = "figures/DustComp_by_Elevation_ALL_gray_sigbars_5.24.21.pdf", width=10, height=8, dpi=600)

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

ggsave(fig.its1.fit1,filename = "figures/DustComp_by_ITS1_ShanDiv_ALL_1.4.22.pdf", width=10, height=8, dpi=600)

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

#ggsave(fig.its1.fit2,filename = "figures/DustComp_by_ITS1_Shan_Div_ALL_5.19.21.pdf", width=10, height=8, dpi=600)


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

ggsave(fig.its1.sr.fit1,filename = "figures/DustComp_by_ITS1_Spec_Richness_ALL_1.4.22.pdf", width=10, height=8, dpi=600)


#### Richness, Diversity vs Env Variables ####

## Shannon Diversity & Environmental Variables
# note: R (correlation coefficient) vs R^2 (coefficient of determination): https://towardsdatascience.com/r%C2%B2-or-r%C2%B2-when-to-use-what-4968eee68ed3

ggplot(bac.div.metadat2, aes(x = DO_Percent_Local, y = Bac_Shannon_Diversity)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampDate), size=3) + theme_classic() +
  stat_smooth(method = "lm", col = "black", se=FALSE, size=1)+ labs(title="Dissolved Oxygen x 16S Shannon Diversity", color="Depth (m)")+ylab("Shannon Diversity")+xlab("Dissolved Oxygen (%)")+
  scale_colour_gradient(low="red",high="blue",guide = guide_colourbar(reverse = TRUE)) +
  scale_shape_discrete(labels=c("June 2021","August 2021","December 2021","April 2022"),name="Sample Date") +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  stat_cor(label.y = 3, label.x=1) +
  stat_regline_equation(aes(label=paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),label.y = 3.1,label.x=1)

ggplot(bac.div.metadat2, aes(x = ORP_mV, y = Bac_Shannon_Diversity)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampDate), size=3) + theme_classic() +
  stat_smooth(method = "lm", col = "black", se=FALSE, size=1)+ labs(title="Oxidation-Reduction Potential x 16S Shannon Diversity", color="Depth (m)")+ylab("Shannon Diversity")+xlab("Redox Potential (mV)")+
  scale_colour_gradient(low="red",high="blue",guide = guide_colourbar(reverse = TRUE)) +
  scale_shape_discrete(labels=c("June 2021","August 2021","December 2021","April 2022"),name="Sample Date") +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  stat_cor(label.y = 3, label.x=1) +
  stat_regline_equation(aes(label=paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),label.y = 3.2,label.x=1)

ggplot(bac.div.metadat2, aes(x = Temp_DegC, y = Bac_Shannon_Diversity)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampDate), size=3) + theme_classic() + labs(title="Temperature x 16S Shannon Diversity", color="Depth (m)")+ylab("16S Shannon Diversity")+xlab("Temperature (C)")+
  scale_colour_gradient(low="red",high="blue",guide = guide_colourbar(reverse = TRUE)) +
  scale_shape_discrete(labels=c("June 2021","August 2021","December 2021","April 2022"),name="Sample Date") +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))

## Species Richness & Environmental Variables
# note: R (correlation coefficient) vs R^2 (coefficient of determination): https://towardsdatascience.com/r%C2%B2-or-r%C2%B2-when-to-use-what-4968eee68ed3

ggplot(bac.div.metadat2, aes(x = DO_Percent_Local, y = Bac_Species_Richness)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampDate), size=3) + theme_classic() +
  stat_smooth(method = "lm", col = "black", se=FALSE, size=1)+ labs(title="Dissolved Oxygen x 16S Species Richness", color="Depth (m)")+ylab("Species Richness")+xlab("Dissolved Oxygen (%)")+
  scale_colour_gradient(low="red",high="blue",guide = guide_colourbar(reverse = TRUE)) +
  scale_shape_discrete(labels=c("June 2021","August 2021","December 2021","April 2022"),name="Sample Date") +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  stat_cor(label.y = 3, label.x=1) +
  stat_regline_equation(aes(label=paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),label.y = 3.1,label.x=1)

ggplot(bac.div.metadat2, aes(x = ORP_mV, y = Bac_Species_Richness)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampDate), size=3) + theme_classic() +
  stat_smooth(method = "lm", col = "black", se=FALSE, size=1)+ labs(title="Oxidation-Reduction Potential x 16S Species Richness", color="Depth (m)")+ylab("Species Richness")+xlab("Redox Potential (mV)")+
  scale_colour_gradient(low="red",high="blue",guide = guide_colourbar(reverse = TRUE)) +
  scale_shape_discrete(labels=c("June 2021","August 2021","December 2021","April 2022"),name="Sample Date") +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  stat_cor(label.y = 3, label.x=1) +
  stat_regline_equation(aes(label=paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),label.y = 3.2,label.x=1)

ggplot(bac.div.metadat2, aes(x = Temp_DegC, y = Bac_Species_Richness)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampDate), size=3) + theme_classic() + labs(title="Temperature x 16S Species Richness", color="Depth (m)")+ylab("16S Species Richness")+xlab("Temperature (C)")+
  scale_colour_gradient(low="red",high="blue",guide = guide_colourbar(reverse = TRUE)) +
  scale_shape_discrete(labels=c("June 2021","August 2021","December 2021","April 2022"),name="Sample Date") +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))

#### Beta Diversity ####
rownames(bac.ASV_table)

# df must have rownames are SampleIDs, columns are ASV IDs for vegan functions below\
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
b.dend_cols <- as.character(metadata$SampDate_Color[order.dendrogram(b.euc_dend)])
labels_colors(b.euc_dend) <- b.dend_cols

plot(b.euc_dend, ylab="CLR Euclidean Distance",cex = 0.5) + title(main = "Bacteria/Archaea Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 3)
#legend("topright",legend = c("June 2021","August 2021","December 2021","April 2022"),cex=.8,col = c( "#26547c","#36ab57","#32cbff","#ff6f00"),pch = 15, bty = "n")
# Control is dark blue ("#218380"), #Alternaria is light blue ("#73d2de")
dev.off()

# let's use our Euclidean distance matrix from before
b.pcoa <- pcoa(b.euc_dist) # pcoa of euclidean distance matrix = PCA of euclidean distance matrix
save.image("data/ssw_clr.euc.dist_2.21.23.Rdata")

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
  labs(title="PCoA: Bacteria/Archaea in Salton Seawater",subtitle="Using Centered-Log Ratio Data",xlab="Axis 1 [41.14%]", ylab="Axis 2 [9.04%]",color="Sample Type")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Type",values=unique(b.pcoa.meta$SampDate_Color[order(b.pcoa.meta$SampDate)]),labels=c("June.2021"="June 2021","August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Axis 1 [22.72%]") + ylab("Axis 2 [17.97%]")

ggsave(pcoa1,filename = "figures/SSW_16S_pcoa_CLR_sampdate.png", width=12, height=10, dpi=600)

# sample month shape, depth color
pcoa2<-ggplot(b.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampleMonth), size=5)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Seawater",subtitle="Using Centered-Log Ratio Data",xlab="Axis 1", ylab="Axis 2",color="Depth (m)")+
  theme_classic()+ theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),axis.text = element_text(size=12),axis.text.x = element_text(vjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17))+
  scale_color_continuous(low="blue3",high="red",trans = 'reverse') + scale_shape_discrete(labels=c("June 2021","August 2021","December 2021","April 2022"),name="Sample Date") +
  xlab("Axis 1 [22.72%]") + ylab("Axis 2 [17.97%]")

ggsave(pcoa2,filename = "figures/SSW_16S_pcoa_CLR_depth.png", width=12, height=10, dpi=600)

## betadisper to look at within group variance
rownames(metadata) %in% rownames(b.euc_dist)

# first by sampling date
b.disper1<-betadisper(b.euc_dist, metadata$SampDate)
b.disper1

permutest(b.disper1, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
#Pairwise comparisons:
#  (Observed p-value below diagonal, permuted p-value above diagonal)
#               June.2021  August.2021 December.2021 April.2022
#June.2021                 0.5840000     0.0140000      0.602
#August.2021   0.5614347                 0.0060000      0.418
#December.2021 0.0187594   0.0089359                    0.399
#April.2022    0.5814591   0.3744861     0.3655430

anova(b.disper1) # p = 0.1998 --> accept the Null H, spatial medians are NOT significantly difference across sample dates

TukeyHSD(b.disper1) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other

# Visualize dispersions
png('figures/pcoa_betadispersion_sampledate.png',width = 700, height = 600, res=100)
plot(b.disper,main = "Centroids and Dispersion based on Aitchison Distance", col=colorset1$SampDate_Color)
dev.off()

png('boxplot_centroid_distance_sampledate.png',width = 700, height = 600, res=100)
boxplot(b.disper,xlab="Sample Collection Date", main = "Distance to Centroid by Category", sub="Based on Aitchison Distance", col=colorset1$Sample_Color)
dev.off()

# What about between sampling depths?
b.disper2<-betadisper(b.euc_dist, metadata$Depth_m)
b.disper2

permutest(b.disper2, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons

anova(b.disper2) # p = 0.7251 --> accept the Null H, spatial medians are NOT significantly difference across sample dates

TukeyHSD(b.disper2) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other

colfunc <- colorRampPalette(c("red", "blue"))
colfunc(9)

# Visualize dispersions
png('figures/pcoa_betadispersion_depth.png',width = 700, height = 600, res=100)
plot(b.disper2,main = "Centroids and Dispersion based on Aitchison Distance", col=colfunc(9))
dev.off()

png('figures/boxplot_centroid_distance_depth.png',width = 700, height = 600, res=100)
boxplot(b.disper2,xlab="Sample Collection Depth", main = "Distance to Centroid by Category", sub="Based on Aitchison Distance", col=colfunc(9))
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
meta_scaled=meta_scaled[rownames(b.clr),] ## reorder metadata to match order of CLR data
perm <- with(meta_scaled, how(nperm = 1000, blocks = SampDate))

adonis2(b.clr ~ SampDate,data=meta_scaled,method = "euclidean",by="terms")
adonis2(b.clr ~ Depth_m,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)

adonis2(b.clr ~ DO_Percent_Local*ORP_mV*Temp_DegC*Chlorophyll_RFU*Dissolved_OrganicMatter_RFU*Depth_m,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
pairwise.adonis2(b.clr.dist ~ Temp_DegC,data=meta_scaled,nperm=perm)

# only using significant variables from previous model comparison
adonis2(b.clr ~ ORP_mV*Temp_DegC*Dissolved_OrganicMatter_RFU,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
adonis2(b.clr ~ Dissolved_OrganicMatter_RFU,data=meta_scaled,method = "euclidean",by="terms",permutations=perm) # Dissolved Organic Matter not significant
# maybe Dissolved organic matter correlates with ORP or temp, which is driving this model significance
adonis2(b.clr ~ ORP_mV,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
adonis2(b.clr ~ Temp_DegC,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)

adonis2(b.clr ~ ORP_mV*Temp_DegC,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
anova(rda(b.clr ~ ORP_mV*Temp_DegC,data=meta_scaled,method = "euclidean",permutations=perm),by="terms") ## same results as previous line!

adonis2(b.clr ~ ORP_mV*Temp_DegC*Depth_m,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
anova(rda(b.clr ~ ORP_mV*Temp_DegC*Depth_m,data=meta_scaled,method = "euclidean",permutations=perm),by="terms") ## same results as previous line!

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


