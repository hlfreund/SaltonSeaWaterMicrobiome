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
#load("data/SSeawater_AlphaDiv_Data.Rdata")
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

### NOTE: bac.div.metadat2 has dropped outliers based on Shannon Diversity!

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
## here the focus is comparing dust complexity to alpha diversity, species richness, depth, & sample date
head(bac.div.metadat2)
s.div.lm.fit1<-lm(Bac_Shannon_Diversity ~ DO_Percent_Local, data=bac.div.metadat2) %>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(s.div.lm.fit1)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)        97.679      6.548  14.916   <2e-16 ***
#  DO_Percent_Local   -1.922      6.647  -0.289    0.774

s.div.lm.fit2<-lm(Bac_Shannon_Diversity ~ ORP_mV, data=bac.div.metadat2) %>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity are continuous data, despite not being normally distributed
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
## ^ went with linear regression because Shannon diversity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(s.div.lm.fit2a)

summer_months<-subset(bac.div.metadat2, SampDate=="June.2021" | SampDate=="August.2021" )

s.div.lm.fit2b<-lm(Bac_Shannon_Diversity ~ ORP_mV, data=summer_months) %>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(s.div.lm.fit2b)

s.div.lm.fit3<-lm(Bac_Shannon_Diversity ~ Temp_DegC, data=bac.div.metadat2) %>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(s.div.lm.fit3)

s.div.lm.fit5<-lm(Bac_Shannon_Diversity ~ Dissolved_OrganicMatter_RFU, data=bac.div.metadat2) %>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity are continuous data, despite not being normally distributed
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

plot(Bac_Shannon_Diversity ~ Depth_m, data=bac.div.metadat2)

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
compare_means(Bac_Shannon_Diversity ~ Depth_m, data=bac.div.metadat2, method="anova",p.adjust.method = "bonferroni") # won't take as.factor(Elevation) as input

p.adj.dc.elev<-compare_means(Bac_Shannon_Diversity ~ Depth_m, data=bac.div.metadat2, method="t.test",p.adjust.method = "bonferroni") # won't take as.factor(Elevation) as input# Note https://github.com/kassambara/ggpubr/issues/65

fit.test<-ggplot(bac.div.metadat2, aes(x = Depth_m, y = Bac_Shannon_Diversity, fill=SampDate)) +
  geom_boxplot() + theme_classic() + guides(fill = guide_legend(reverse=TRUE)) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  labs(title="Dust Complexity x Elevation",fill="Elevation (ft)")+ylab("Dust Complexity")+xlab("Elevation (ft)")+scale_fill_manual(values=saturation(fair_cols, 0.9))+stat_compare_means(method = "anova",label.y=1.5) +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4)), method="t.test", hide.ns = TRUE,label = "p.signif")

#ggsave(fit.test,filename = "figures/DustComp_by_Elevation_ALL_sigbars_5.24.21.pdf", width=10, height=8, dpi=600)

## Now for richness...

head(bac.div.metadat2)

s.sr.lm.fit1<-lm(Bac_Species_Richness ~ DO_Percent_Local, data=bac.div.metadat2) %>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Species richness are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(s.sr.lm.fit1)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)        950.80      54.40  17.479   <2e-16 ***
#  DO_Percent_Local   -10.44      55.22  -0.189    0.851

s.sr.lm.fit2<-lm(Bac_Species_Richness ~ ORP_mV, data=bac.div.metadat2) %>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because SSpecies richness are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(s.sr.lm.fit2)

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)   97.674      3.709  26.333  < 2e-16 ***
#  ORP_mV     -93.98      33.63  -2.794  0.00768 **
## ^^^ the two lms below show that this model is significant only for June & August 2021, not December & April

not_summer_months<-subset(bac.div.metadat2, SampDate=="December.2021" | SampDate=="April.2022" )

s.sr.lm.fit2a<-lm(Bac_Species_Richness ~ ORP_mV, data=not_summer_months) %>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon srersity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(s.sr.lm.fit2a)

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)   3409.4      971.9   3.508  0.00149 **
#  ORP_mV       -4694.2     1789.7  -2.623  0.01375 *

summer_months<-subset(bac.div.metadat2, SampDate=="June.2021" | SampDate=="August.2021" )

s.sr.lm.fit2b<-lm(Bac_Species_Richness ~ ORP_mV, data=summer_months) %>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Species richness are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(s.sr.lm.fit2b)

# Coefficients:
#            Estimate Std. Error t value Pr(>|t|)
#(Intercept)  1088.19      70.88  15.353 1.04e-09 ***
#  ORP_mV        -27.43      45.38  -0.605    0.556

s.sr.lm.fit3<-lm(Bac_Species_Richness ~ Temp_DegC, data=bac.div.metadat2) %>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Species richness are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(s.sr.lm.fit3)

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)    639.8      141.3   4.529 4.49e-05 ***
#  Temp_DegC      314.6      142.1   2.214    0.032 *

s.sr.lm.fit5<-lm(Bac_Species_Richness ~ Dissolved_OrganicMatter_RFU, data=bac.div.metadat2) %>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Species richness are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(s.sr.lm.fit5)

fit2<-aov(Bac_Species_Richness ~ Depth_m, data=bac.div.metadat2)
pairwise.adonis(bac.div.metadat2$Bac_Species_Richness, bac.div.metadat2$Depth_m, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different
#adonis2(b.clr ~ DO_Percent_Local*ORP_mV*Temp_DegC*Dissolved_OrganicMatter_RFU*Depth_m,data=meta_scaled,method = "euclidean",by="terms",permutations=perm)
#test<-adonis2(bac.div.metadat2$Bac_Species_Richness ~ Depth_m, data=bac.div.metadat2)

summary(fit2)
#Df           Sum Sq Mean Sq    F value   Pr(>F)
#Elevation2   3 0.7277 0.24258   0.084 0.774
#Residuals   27 0.4444 0.01646
Tuk1<-TukeyHSD(fit2)
Tuk1$Depth_m

plot(Bac_Species_Richness ~ Depth_m, data=bac.div.metadat2)

# fit.0<-aov(DustComplexity ~ as.factor(Elevation), data=bac.div.metadat2)
# summary(fit.0)
# TukeyHSD(fit.0)
# Levene's test with one independent variable
## Levene's tests whether variances of 2 samples are equal
## we want variances to be the same -- want NON SIGNIFICANCE!
## t test assumes that variances are the same, so Levene's test needs to be non significant
fligner.test(Bac_Species_Richness ~ Depth_m, data = bac.div.metadat2)
# Levenes Test for Homogeneity of Variance
#        Df  Chi square value  Pr(>F)
# group  3   1.0952   0.7411
# Which shows that the data do not deviate significantly from homogeneity.
compare_means(Bac_Species_Richness ~ Depth_m, data=bac.div.metadat2, method="anova",p.adjust.method = "bonferroni") # won't take as.factor(Elevation) as input

p.adj.dc.elev<-compare_means(Bac_Species_Richness ~ Depth_m, data=bac.div.metadat2, method="t.test",p.adjust.method = "bonferroni") # won't take as.factor(Elevation) as input# Note https://github.com/kassambara/ggpubr/issues/65

fit.test<-ggplot(bac.div.metadat2, aes(x = Depth_m, y = Bac_Species_Richness, fill=SampDate)) +
  geom_boxplot() + theme_classic() + guides(fill = guide_legend(reverse=TRUE)) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  labs(title="Dust Complexity x Elevation",fill="Elevation (ft)")+ylab("Dust Complexity")+xlab("Elevation (ft)")+scale_fill_manual(values=saturation(fair_cols, 0.9))+stat_compare_means(method = "anova",label.y=1.5) +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4)), method="t.test", hide.ns = TRUE,label = "p.signif")

#ggsave(fit.test,filename = "figures/DustComp_by_Elevation_ALL_sigbars_5.24.21.pdf", width=10, height=8, dpi=600)

#### Richness, Diversity vs Env Variables ####

## Shannon Diversity & Environmental Variables
# note: R (correlation coefficient) vs R^2 (coefficient of determination): https://towardsdatascience.com/r%C2%B2-or-r%C2%B2-when-to-use-what-4968eee68ed3

ggplot(bac.div.metadat2, aes(x = DO_Percent_Local, y = Bac_Shannon_Diversity)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampDate), size=3) + theme_classic() +
  stat_smooth(method = "lm", col = "black", se=FALSE, size=1)+ labs(title="Dissolved Oxygen x 16S Shannon Diversity", color="Depth (m)")+ylab("Shannon Diversity")+xlab("Dissolved Oxygen (%)")+
  scale_colour_gradient(low="red",high="blue",guide = guide_colourbar(reverse = TRUE)) +
  scale_shape_discrete(labels=c("June 2021","August 2021","December 2021","April 2022"),name="Sample Date") +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  stat_cor(label.y = 150, label.x=3) +
  stat_regline_equation(aes(label=paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),label.y = 160,label.x=3)

ggplot(bac.div.metadat2, aes(x = DO_Percent_Local, y = Bac_Shannon_Diversity)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampDate), size=3) + theme_classic() +
  labs(title="Dissolved Oxygen x 16S Shannon Diversity", color="Depth (m)")+ylab("Shannon Diversity")+xlab("Dissolved Oxygen (%)")+
  scale_colour_gradient(low="red",high="blue",guide = guide_colourbar(reverse = TRUE)) +
  scale_shape_discrete(labels=c("June 2021","August 2021","December 2021","April 2022"),name="Sample Date") +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))

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

