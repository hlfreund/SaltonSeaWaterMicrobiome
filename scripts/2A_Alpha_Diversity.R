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
load("data/SSeawater_AlphaDiv_Data.Rdata")
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
png('figures/AlphaDiversity/SSW_16S_rarecurve.png')
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
bac.div.metadat$Depth_m<-factor(bac.div.metadat$Depth_m, levels=c("0","3","4","5","7","9","10","11"))

# drop the outliers
bac.div.metadat2<-subset(bac.div.metadat, bac.div.metadat$Bac_Shannon_Diversity<=200)

# save diversity data
save.image("data/SSeawater_AlphaDiv_Data.Rdata")

#### Visualize Alpha Diversity & Species Richness ####
## Shannon Diversity by Sample Month & Depth
bac.a.div<-ggplot(bac.div.metadat2, aes(x=SampDate, y=Bac_Shannon_Diversity,fill=SampDate)) +geom_boxplot(color="black")+scale_x_discrete(labels=c("August 2021","December 2021","April 2022"))+theme_bw()+
  scale_fill_manual(values=unique(bac.div.metadat$SampDate_Color[order(bac.div.metadat$SampDate)]), name ="Sample Date",labels=c("August 2021","December 2021","April 2022"))+theme_classic()+
  labs(title = "Bacterial Shannon Diversity by Sample Date", x="Sample Date", y="Shannon Diversity", fill="Sample Month")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1,size=10),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3)), method="t.test", hide.ns = TRUE,label = "p.format")

ggsave(bac.a.div,filename = "figures/AlphaDiversity/SSW_Bacterial_alpha_diversity_sampledate.png", width=13, height=10, dpi=600)

bac.div.metadat2$Depth_m=as.numeric(levels(bac.div.metadat2$Depth_m))[bac.div.metadat2$Depth_m]
# ^ note: cannot turn numbers that are factors in R into numeric values...
## have to convert factor levels into numeric, then use the numeric "levels" to pull out numbers from Depth_m column in df to make sure the Depth_m columns is now numeric, not a factor

bac.a.div2<-ggplot(bac.div.metadat2, aes(x=as.factor(Depth_m), y=Bac_Shannon_Diversity)) +geom_boxplot(aes(fill=bac.div.metadat2$Depth_m),color="black")+
  labs(title = "Bacterial Shannon Diversity by Sampling Depth", x="Depth (m)", y="Shannon Diversity", fill="Depth (m)")+
  scale_fill_gradient(low="red",high="blue",guide = guide_colourbar(reverse = TRUE)) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1,,size=10),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  coord_flip() + scale_x_discrete(limits=rev)

ggsave(bac.a.div2,filename = "figures/AlphaDiversity/SSW_Bacterial_alpha_diversity_depth_v1.png", width=13, height=10, dpi=600)

bac.a.div3<-ggplot(bac.div.metadat2, aes(x=as.factor(Depth_m), y=Bac_Shannon_Diversity)) +geom_boxplot(aes(fill=bac.div.metadat2$Depth_m),color="black")+
  labs(title = "Bacterial Shannon Diversity by Sampling Depth", x="Depth (m)", y="Shannon Diversity", fill="Depth (m)")+
  scale_fill_gradient(low="red",high="blue",guide = guide_colourbar(reverse = TRUE)) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1,,size=10),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))

ggsave(bac.a.div3,filename = "figures/AlphaDiversity/SSW_Bacterial_alpha_diversity_depth_v2.png", width=13, height=10, dpi=600)

## Species Richness by Sample Type
bac.a.sr<-ggplot(bac.div.metadat2, aes(x=SampDate, y=Bac_Species_Richness,fill=SampDate)) +geom_boxplot(color="black")+scale_x_discrete(labels=c("August 2021","December 2021","April 2022"))+theme_bw()+
  scale_fill_manual(values=unique(bac.div.metadat$SampDate_Color[order(bac.div.metadat$SampDate)]), name ="Sample Date",labels=c("August 2021","December 2021","April 2022"))+theme_classic()+
  labs(title = "Bacterial Species Richness by Sample Date", x="Sample Date", y="Species Richness")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1,size=10),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3)), method="t.test", hide.ns = TRUE,label = "p.format")

ggsave(bac.a.sr,filename = "figures/AlphaDiversity/SSW_Bacterial_species_richness_samplemonth.png", width=13, height=10, dpi=600)

bac.a.sr2<-ggplot(bac.div.metadat2, aes(x=as.factor(Depth_m), y=Bac_Species_Richness,fill=bac.div.metadat2$Depth_m)) +geom_boxplot(aes(fill=as.numeric(bac.div.metadat2$Depth_m)),color="black")+
  labs(title = "Bacterial Species Richness by Sampling Depth", x="Depth (m)", y="Species Richness", fill="Depth (m)")+
  scale_fill_gradient(low="red",high="blue",guide = guide_colourbar(reverse = TRUE)) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1,,size=10),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  coord_flip() + scale_x_discrete(limits=rev)

ggsave(bac.a.sr2,filename = "figures/AlphaDiversity/SSW_Bacterial_species_richness_depth_v1.png", width=13, height=10, dpi=600)

bac.a.sr3<-ggplot(bac.div.metadat2, aes(x=as.factor(Depth_m), y=Bac_Species_Richness,fill=bac.div.metadat2$Depth_m)) +geom_boxplot(aes(fill=as.numeric(bac.div.metadat2$Depth_m)),color="black")+
  labs(title = "Bacterial Species Richness by Sampling Depth", x="Depth (m)", y="Species Richness", fill="Depth (m)")+
  scale_fill_gradient(low="red",high="blue",guide = guide_colourbar(reverse = TRUE)) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1,,size=10),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))

ggsave(bac.a.sr3,filename = "figures/AlphaDiversity/SSW_Bacterial_species_richness_depth_v2.png", width=13, height=10, dpi=600)

#### Using Shapiro-Wilk test for Normality ####
shapiro.test(bac.div.metadat$Bac_Shannon_Diversity) # what is the p-value?
# p-value = 1.281e-12 including outliers
shapiro.test(bac.div.metadat2$Bac_Shannon_Diversity) # what is the p-value? No outliers **
# p-value = 0.006545; excluding outliars
# p > 0.05 states distribution of data are not significantly different from normal distribution
# p < 0.05 means that data is significantly different from a normal distribution
hist(bac.div.metadat2$Bac_Shannon_Diversity, col="blue")

# visualize Q-Q plot for species richness
qqnorm(bac.div.metadat2$Bac_Shannon_Diversity, pch = 1, frame = FALSE)
qqline(bac.div.metadat2$Bac_Shannon_Diversity, col = "steelblue", lwd = 2)

shapiro.test(bac.div.metadat$Bac_Species_Richness) # what is the p-value?
# p-value = 0.02873 w/ outliars
shapiro.test(bac.div.metadat2$Bac_Species_Richness) # what is the p-value? * No outliars
# p-value =  0.009811; no outliars
# p > 0.05 states distribution of data are not significantly different from normal distribution
# p < 0.05 means that data is significantly different from a normal distribution
hist(bac.div.metadat$Bac_Species_Richness, col="blue")

# visualize Q-Q plot for species richness
qqnorm(bac.div.metadat$Bac_Species_Richness, pch = 1, frame = FALSE)
qqline(bac.div.metadat$Bac_Species_Richness, col = "steelblue", lwd = 2)

### NOTE: bac.div.metadat2 has dropped outliers based on Shannon Diversity!

shapiro.test(bac.div.metadat$DO_Percent_Local) # p-value = 0.0007935
hist(bac.div.metadat$DO_Percent_Local, col="blue")

shapiro.test(bac.div.metadat2$ORP_mV) # p-value = 5.255e-12
hist(bac.div.metadat$ORP_mV, col="blue")

shapiro.test(bac.div.metadat2$Temp_DegC) # p-value = 5.39e-06
hist(bac.div.metadat$Temp_DegC, col="blue")

shapiro.test(bac.div.metadat$Dissolved_OrganicMatter_RFU) #  p-value = 0.0003007
hist(bac.div.metadat$Dissolved_OrganicMatter_RFU, col="blue")

shapiro.test(bac.div.metadat2$Sulfate_milliM) # p-value = 0.01146
hist(bac.div.metadat$Sulfate_milliM, col="blue")

shapiro.test(bac.div.metadat2$Sulfide_microM) # p-value = 9.566e-12
hist(bac.div.metadat$Sulfide_microM, col="blue")

#### Linear Regression Comparisons - Shannon Diversity ####
## here the focus is comparing dust complexity to alpha diversity, species richness, & elevation
head(bac.div.metadat2)
s.div.lm.fit1<-lm(Bac_Shannon_Diversity ~ DO_Percent_Local, data=bac.div.metadat2) %>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(s.div.lm.fit1)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)        89.279      3.587  24.893   <2e-16 ***
#DO_Percent_Local    7.088      3.683   1.925    0.062 .

s.div.lm.fit2<-lm(Bac_Shannon_Diversity ~ ORP_mV, data=bac.div.metadat2) %>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(s.div.lm.fit2)

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)   89.028      3.755  23.712   <2e-16 ***
#ORP_mV         1.142      3.757   0.304    0.763

s.div.lm.fit3<-lm(Bac_Shannon_Diversity ~ Temp_DegC, data=bac.div.metadat2) %>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(s.div.lm.fit3)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)   89.202      3.485  25.596   <2e-16 ***
#  Temp_DegC     -8.659      3.515  -2.464   0.0185 *

s.div.lm.fit5<-lm(Bac_Shannon_Diversity ~ Dissolved_OrganicMatter_RFU, data=bac.div.metadat2) %>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(s.div.lm.fit5)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)                   89.024      3.717  23.949   <2e-16 ***
#  Dissolved_OrganicMatter_RFU   -3.405      3.717  -0.916    0.366

s.div.lm.fit6<-lm(Bac_Shannon_Diversity ~ Sulfate_milliM, data=bac.div.metadat2) %>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(s.div.lm.fit6)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)      89.075      3.745   23.78   <2e-16 ***
#  Sulfate_milliM    2.087      3.795    0.55    0.586

s.div.lm.fit7<-lm(Bac_Shannon_Diversity ~ Sulfide_microM, data=bac.div.metadat2) %>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(s.div.lm.fit7)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)      89.039      3.743   23.79   <2e-16 ***
#  Sulfide_microM   -2.137      3.748   -0.57    0.572

fit1<-aov(Bac_Shannon_Diversity ~ as.factor(Depth_m), data=bac.div.metadat2)
#pairwise.adonis(bac.div.metadat2$Bac_Shannon_Diversity, bac.div.metadat2$Depth_m, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different

summary(fit1)
#Df           Sum Sq Mean Sq    F value   Pr(>F)
#Depth_m      7   4097   585.3   1.114   0.38
#Residuals   31  16294   525.6
Tuk1<-TukeyHSD(fit1)
Tuk1$Depth_m

#plot(Bac_Shannon_Diversity ~ Depth_m, data=bac.div.metadat2)
#abline(aov(DustComplexity ~ Elevation, data=bac.div.metadat2))

# fit.0<-aov(DustComplexity ~ as.factor(Elevation), data=bac.div.metadat2)
# summary(fit.0)
# TukeyHSD(fit.0)
# Levene's test with one independent variable
## Levene's tests whether variances of 2 samples are equal
## we want variances to be the same -- want NON SIGNIFICANCE!
## t test assumes that variances are the same, so Levene's test needs to be non significant
## Fligner's test is a Levene's test for data that are not normally distributed
## more here: https://www.geeksforgeeks.org/fligner-killeen-test-in-r-programming/
fligner.test(Bac_Shannon_Diversity ~ Depth_m, data = bac.div.metadat2)
# Fligner-Killeen:med chi-squared = 4.091, df = 7, p-value = 0.7692
# Which shows that the data do not deviate significantly from homogeneity.
compare_means(Bac_Shannon_Diversity ~ Depth_m, data=bac.div.metadat2, method="anova",p.adjust.method = "bonferroni") # won't take as.factor(Elevation) as input

#### Linear Regression Comparisons - Species Richness ####
## here the focus is comparing dust complexity to alpha diversity, species richness, & elevation
head(bac.div.metadat2)
s.r.lm.fit1<-lm(Bac_Species_Richness ~ DO_Percent_Local, data=bac.div.metadat2) %>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(s.r.lm.fit1)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)        917.14      39.73   23.08   <2e-16 ***
#DO_Percent_Local   -37.11      40.79   -0.91    0.369

s.r.lm.fit2<-lm(Bac_Species_Richness ~ ORP_mV, data=bac.div.metadat2) %>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(s.r.lm.fit2)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)   917.95      37.45  24.514   <2e-16 ***
#  ORP_mV        -88.09      37.47  -2.351   0.0242 *

s.r.lm.fit3<-lm(Bac_Species_Richness ~ Temp_DegC, data=bac.div.metadat2) %>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(s.r.lm.fit3)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)   917.11      38.69   23.71   <2e-16 ***
#Temp_DegC      65.95      39.02    1.69   0.0994 .

s.r.lm.fit5<-lm(Bac_Species_Richness ~ Dissolved_OrganicMatter_RFU, data=bac.div.metadat2) %>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(s.r.lm.fit5)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)                   918.41      35.34  25.989  < 2e-16 ***
#Dissolved_OrganicMatter_RFU   115.85      35.34   3.278  0.00228 **

s.r.lm.fit6<-lm(Bac_Species_Richness ~ Sulfate_milliM, data=bac.div.metadat2) %>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(s.r.lm.fit6)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)      918.89      40.08  22.928   <2e-16 ***
#Sulfate_milliM    15.65      40.61   0.385    0.702

s.r.lm.fit7<-lm(Bac_Species_Richness ~ Sulfide_microM, data=bac.div.metadat2) %>%
  adjust_pvalue(method="bonferroni")
## ^ went with linear regression because Shannon diversity and dust complexity are continuous data, despite not being normally distributed
# model form is response ~ terms (y ~ x) where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
summary(s.r.lm.fit7)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)      917.85      38.09  24.094   <2e-16 ***
#Sulfide_microM    77.15      38.15   2.023   0.0504 .

fit2<-aov(Bac_Species_Richness ~ as.factor(Depth_m), data=bac.div.metadat2)
#pairwise.adonis(bac.div.metadat2$Bac_Species_Richness, bac.div.metadat2$Depth_m, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different

summary(fit2)
#Df           Sum Sq Mean Sq    F value   Pr(>F)
#as.factor(Depth_m)  7  243251   34750   0.517  0.814
#Residuals          31 2082211   67168
Tuk2<-TukeyHSD(fit2)
Tuk2$Depth_m

#plot(DustComplexity ~ Elevation, data=bac.div.metadat2)
#abline(aov(DustComplexity ~ Elevation, data=bac.div.metadat2))

# Levene's test with one independent variable
## Levene's tests whether variances of 2 samples are equal
## we want variances to be the same -- want NON SIGNIFICANCE!
## t test assumes that variances are the same, so Levene's test needs to be non significant
## Fligner's test is a Levene's test for data that are not normally distributed
## more here: https://www.geeksforgeeks.org/fligner-killeen-test-in-r-programming/
fligner.test(Bac_Species_Richness ~ Depth_m, data = bac.div.metadat2)
# Fligner-Killeen:med chi-squared = 3.2235, df = 7, p-value = 0.8636
# Which shows that the data do not deviate significantly from homogeneity.
compare_means(Bac_Species_Richness ~ Depth_m, data=bac.div.metadat2, method="anova",p.adjust.method = "bonferroni")

#### Save Everything ####
save.image("data/SSeawater_AlphaDiv_Data.Rdata")
