#### Set WD & Load Libraries ####
getwd() # use setwd("path/to/files") if you are not in the right directory
setwd("/Volumes/HLF_SSD/Aronson_Lab_Data/Salton_Sea/SaltonSeaWater")
suppressPackageStartupMessages({ # load packages quietly
  library(devtools)
  library(phyloseq)
  library(ggplot2)
  library(vegan)
  library(lme4)
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

# NOTE: chem data has already been Raw in meta_scaled; raw chem/env data in metadata
# all samples included (outliers from alpha div data were not removed)

bac.dat.all[1:6,1:6]
bac.ASV_table[1:4,1:4]
head(metadata) # not Raw metadata
head(meta_scaled) # centered & Raw metadata

# create column for Depth that is a numeric version of this variable, rather than a factor
metadata$Depth.num<-as.numeric(as.character(metadata$Depth_m)) # for env PCA w/ log -transformed variables
meta_scaled$Depth.num<-as.numeric(as.character(meta_scaled$Depth_m))

meta_salinity<-metadata # save metadata with salinity data to just graph salinity data later

# drop salinity from metadata & meta_scaled --> excluding this env variable
metadata<-subset(metadata, select=-c(Salinity_ppt))
head(metadata)

meta_scaled<-subset(meta_scaled, select=-c(Salinity_ppt))
head(meta_scaled)

#### Separate All Data by Timepoints ####
# create metadata df that will contain scaled chemical data
head(metadata)
head(meta_scaled)

site_list<-unique(meta_scaled$SampDate) #define an array of string values
# go through metadata & create a list of data frames
## when metadata$Variable == element in site_list (aka x in this case), subset metadata by said element into elements of a list

# here the function(x) is using site_list aka x to subset metadata, when $Variable column == site_list
# Run the function so it's stored in Global Env
site_subsets<-lapply(site_list, function(x) {subset(meta_scaled, SampDate==x)})

site_subsets # sanity check1 (should see all elements in list)
site_subsets[[1]] # sanity check2 (see 1st element in list)
#rename the list elements

# name each element in list
names(site_subsets)<-site_list # * only do this if the order of names in site_list match order of the elements in site_subsets!
site_subsets$April.2022 # sanity check3 - should be able to pull dataframes by names rather than index now

# example of subsetting
site_subsets[[2]][1:3]
site_subsets$August.2021[1:3] # should produce same ouptut as line above

site_subsets[[2]][1:2,1:2] # another example

# ^ subsetting to [[second dataframe]], [[row #, column #]]
site_subsets[[2]][[1,2]] # [[second dataframe]], [[row 1, column 2]]

# set up the function and run this to store it in our Global environment
df_specific.subset<-function(var_vec,var_subsets){
  # var_vec = vector of variable elements from specific categorical variable;
  ## e.g. vector of names from Site categorical variable (metadata sites)
  # var_subsets = list of dataframes subsetted by column$element from original dataframe;
  ## e.g. list of dataframes (each df = element of list) subsetted from metadata using vector of metadata$Site names
  for(i in seq_along(var_vec)){
    # print(var_vec[i]) -- var_vec[i] = each element in var_vec
    # print(var_subsets[[i]]) -- var_subsets[[i]] = each sub
    df<-paste(var_vec[i])
    #print(df)
    assign(df, var_subsets[[i]], envir = .GlobalEnv)
    print(paste("Dataframe", var_vec[i] ,"done"))

  }

}

# run the function
df_specific.subset(site_list, site_subsets) # used scaled metadata quantitative values

head(August.2021) # sanity check
August.2021[1:5,] # double check that our new Variable (here SampDate) data frames still have scaled chemical data
rownames(August.2021)

#### Using Shapiro-Wilk test for Normality ####

shapiro.test(meta_scaled$DO_Percent_Local) # what is the p-value?
# my p-value was p-value =  0.02586
# p > 0.05 states distribution of data are not significantly different from normal distribution
# p < 0.05 means that data is significantly different from a normal distribution
hist(meta_scaled$DO_Percent_Local, col="blue")

# visualize Q-Q plot for species richness
qqnorm(meta_scaled$DO_Percent_Local, pch = 1, frame = FALSE)
qqline(meta_scaled$DO_Percent_Local, col = "steelblue", lwd = 2)

shapiro.test(meta_scaled$ORP_mV) # what is the p-value? p-value = 1.731e-08
hist(meta_scaled$ORP_mV, col="blue")
# visualize Q-Q plot for species richness
qqnorm(meta_scaled$ORP_mV, pch = 1, frame = FALSE)
qqline(meta_scaled$ORP_mV, col = "steelblue", lwd = 2)

shapiro.test(meta_scaled$Temp_DegC) # what is the p-value? p-value = 0.0002829
hist(meta_scaled$Temp_DegC, col="blue")
# visualize Q-Q plot for species richness
qqnorm(meta_scaled$Temp_DegC, pch = 1, frame = FALSE)
qqline(meta_scaled$Temp_DegC, col = "steelblue", lwd = 2)

shapiro.test(meta_scaled$Dissolved_OrganicMatter_RFU) # what is the p-value? p-value = 0.05411
hist(meta_scaled$Dissolved_OrganicMatter_RFU, col="blue")
# visualize Q-Q plot for species richness
qqnorm(meta_scaled$Dissolved_OrganicMatter_RFU, pch = 1, frame = FALSE)
qqline(meta_scaled$Dissolved_OrganicMatter_RFU, col = "steelblue", lwd = 2)

shapiro.test(meta_scaled$Sulfate_milliM) # what is the p-value? p-value =  0.1912
# p > 0.05 states distribution of data are not significantly different from normal distribution
# p < 0.05 means that data is significantly different from a normal distribution
hist(meta_scaled$Sulfate_milliM, col="blue")

# visualize Q-Q plot for species richness
qqnorm(meta_scaled$Sulfate_milliM, pch = 1, frame = FALSE)
qqline(meta_scaled$Sulfate_milliM, col = "steelblue", lwd = 2)

shapiro.test(meta_scaled$Sulfide_microM) # what is the p-value? p-value =  3.813e-08
# p > 0.05 states distribution of data are not significantly different from normal distribution
# p < 0.05 means that data is significantly different from a normal distribution
hist(meta_scaled$Sulfide_microM, col="blue")

# visualize Q-Q plot for species richness
qqnorm(meta_scaled$Sulfide_microM, pch = 1, frame = FALSE)
qqline(meta_scaled$Sulfide_microM, col = "steelblue", lwd = 2)

shapiro.test(meta_scaled$Turbidity_FNU) # what is the p-value?  p-value = 0.002374
hist(meta_scaled$Turbidity_FNU, col="blue")
# visualize Q-Q plot for species richness
qqnorm(meta_scaled$Turbidity_FNU, pch = 1, frame = FALSE)
qqline(meta_scaled$Turbidity_FNU, col = "steelblue", lwd = 2)

shapiro.test(meta_scaled$Chlorophyll_RFU) # what is the p-value? p-value = 0.1947
hist(meta_scaled$Chlorophyll_RFU, col="blue")
# visualize Q-Q plot for species richness
qqnorm(meta_scaled$Chlorophyll_RFU, pch = 1, frame = FALSE)
qqline(meta_scaled$Chlorophyll_RFU, col = "steelblue", lwd = 2)

#### PCA w/ Env Variables ####
head(metadata)
env.dat<-metadata[,c(8,10:11,14:16)]
head(env.dat)

# NOTE: PCA requires normally distributed data, so we are log transforming env data before scaling
# this is because Raw env data alone is not normally distributed
env.log<-decostand(env.dat,method = "log", pseudocount = 1) # log transformation of env data before scaling
env.log[1:4,1:4]

# check rownames of log transformed environmental data & metadata
rownames(env.log) %in% rownames(meta_scaled)
#meta_scaled=meta_scaled[rownames(env.log),] ## reorder metadata to match order of log data

# calculate our Euclidean distance matrix using log data
env.euc_dist <- dist(env.log, method = "euclidean")

# creating our hierarcical clustering dendrogram
env.euc_clust <- hclust(env.euc_dist, method="ward.D2")

# let's make it a little nicer...
env.euc_dend <- as.dendrogram(env.euc_clust, hang=0.2)
env.dend_cols <- as.character(meta_scaled$SampDate_Color[order.dendrogram(env.euc_dend)])
labels_colors(env.euc_dend) <- env.dend_cols

plot(env.euc_dend, ylab="Log Euclidean Distance",cex = 0.5) + title(main = "Environmental Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 3)
#legend("topright",legend = c("August 2021","December 2021","April 2022"),cex=.8,col = c( "#ef781c","#03045e","#32cbff","#059c3f"),pch = 15, bty = "n")
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
pca1<-ggplot(env.pca.meta, aes(x=PC1, y=PC2)) +geom_point(aes(color=factor(SampDate)), size=4)+theme_bw()+
  labs(title="PCA: Environmental Variables in Salton Seawater",subtitle="Using Log-Transformed Data",color="Sample Type")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Type",values=unique(env.pca.meta$SampDate_Color[order(env.pca.meta$SampDate)]),labels=c("June.2021"="June 2021","August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("PC1 [64.74%]") + ylab("PC2 [19.73%]")

ggsave(pca1,filename = "figures/EnvVariablesOnly/SSW_LogEnvOnly_PCA_SampDate.png", width=12, height=10, dpi=600)

# sample month shape, depth color
pca2<-ggplot(env.pca.meta, aes(x=PC1, y=PC2)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampleMonth), size=5)+theme_bw()+
  labs(title="PCA: Environmental Variables in Salton Seawater",subtitle="Using Log-Transformed Data",color="Depth (m)")+
  theme_classic()+ theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),axis.text = element_text(size=12),axis.text.x = element_text(vjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17))+
  scale_color_continuous(low="blue3",high="red",trans = 'reverse') + scale_shape_discrete(labels=c("August 2021","December 2021","April 2022"),name="Sample Date") +
  xlab("PC1 [64.74%]") + ylab("PC2 [19.73%]")

ggsave(pca2,filename = "figures/EnvVariablesOnly/SSW_LogEnvOnly_PCA_Depth_SampDate.png", width=12, height=10, dpi=600)

#### Do Env Variables Vary by Sample Date ####
# Kruskal-Wallis test = nonparametric one-way ANOVA
kruskal.test(DO_Percent_Local ~ SampDate, data = meta_scaled)
#Kruskal-Wallis chi-squared = 16.266, df = 2, p-value = 0.0002937
pairwise.wilcox.test(meta_scaled$DO_Percent_Local, meta_scaled$SampDate, p.adjust.method = "bonf") # returns p values
#               August  December
#December.2021 0.0028      -
#April.2022    0.0299  0.0299

kruskal.test(ORP_mV ~ SampDate, data = meta_scaled)
#Kruskal-Wallis chi-squared = 19.763, df = 2, p-value = 5.112e-05
pairwise.wilcox.test(meta_scaled$ORP_mV, meta_scaled$SampDate, p.adjust.method = "bonf") # returns p values
#               August.2021 December.2021
# December.2021 0.0068      -
# April.2022    0.0027      0.0027

kruskal.test(Temp_DegC ~ SampDate, data = meta_scaled)
#Kruskal-Wallis chi-squared = 20.507, df = 2, p-value = 3.524e-05
pairwise.wilcox.test(meta_scaled$Temp_DegC, meta_scaled$SampDate, p.adjust.method = "bonf") # returns p values
#               August.2021 December.2021
# December.2021 0.0028      -
#   April.2022  0.0028      0.0028

kruskal.test(Dissolved_OrganicMatter_RFU ~ SampDate, data = meta_scaled)
#Kruskal-Wallis chi-squared = 20.516, df = 2, p-value = 3.508e-05
pairwise.wilcox.test(meta_scaled$Dissolved_OrganicMatter_RFU, meta_scaled$SampDate, p.adjust.method = "bonf") # returns p values
#               August.2021 December.2021
# December.2021 0.0027      -
# April.2022    0.0028      0.0027

#kruskal.test(Sulfate_milliM ~ SampDate, data = meta_scaled)
#Kruskal-Wallis chi-squared = 15.605, df = 2, p-value = 0.0004087
#pairwise.wilcox.test(meta_scaled$Sulfate_milliM, meta_scaled$SampDate, p.adjust.method = "bonf") # returns p values
#               August.2021 December.2021
# December.2021 0.00047     -
# April.2022    1.00000     0.00047

# sulfate is normally distributed -- using ANOVA
sulf.fit1<-aov(Sulfate_milliM ~ SampDate, data=meta_scaled)

summary(sulf.fit1)
#Df           Sum Sq Mean Sq    F value   Pr(>F)
#SampDate     2 16.188   8.094   24.95 2.83e-06 ***
#Residuals   21  6.812   0.324

# Tukey test - tells us which groups are significantly different from each other (more here: https://www.r-bloggers.com/2013/06/anova-and-tukeys-test-on-r/)
Tuk.Sulf<-TukeyHSD(sulf.fit1)
Tuk.Sulf$SampDate
#                                 diff        lwr        upr        p adj
#December.2021-August.2021  1.70414819  0.9863543  2.4219421 1.766873e-05
#April.2022-August.2021    -0.07374444 -0.7915383  0.6440494 9.637716e-01
#April.2022-December.2021  -1.77789263 -2.4956865 -1.0600987 9.888236e-06

kruskal.test(Sulfide_microM ~ SampDate, data = meta_scaled)
#Kruskal-Wallis chi-squared = 10.75, df = 2, p-value = 0.004632
pairwise.wilcox.test(meta_scaled$Sulfide_microM, meta_scaled$SampDate, p.adjust.method = "bonf") # returns p values
#               August.2021 December.2021
# December.2021 0.030       -
# April.2022    0.703       0.016

#### Env Variable Corrplots ####
head(meta_scaled)
# check for colinearity among env variables themselves
heatmap(abs(cor(meta_scaled[,c(8,10:11,14:16)])),
        # Compute pearson correlation (note they are absolute values)
        col = rev(heat.colors(6)),
        Colv = NA, Rowv = NA)
legend("topleft",
       title = "Absolute Pearson R",
       legend =  round(seq(0,1, length.out = 6),1),
       y.intersp = 0.7, bty = "n",
       fill = rev(heat.colors(6)))
dev.off()

# Calculate correlations for corr coefficient & p values
cor(meta_scaled[,c(8,10:11,14:16)],method='pearson')
cor.all.mat = cor.mtest(meta_scaled[,c(8,10:11,14:16)],method='pearson', conf.level = 0.95)

# Visualize with a corrplot
cor_mat.env1 <- cor(meta_scaled[,c(8,10:11,14:16)], method='pearson')
cor_mat.env1

symnum(cor_mat.env1)

tiff('figures/EnvVariablesOnly/SSW_ScaledCentered_EnvVarOnly_AllData_CorrPlot.tiff', width = 7, height = 7, units = 'in', res = 300)
corrplot.mixed(cor_mat.env1, p.mat=cor.all.mat$p, tl.pos='lt', tl.cex=0.7, sig.level = 0.05, order="alphabet",insig='blank',number.cex=0.8,
               diag='n',cl.ratio = 0.2, tl.srt = 45,title="All Timepoints",mar=c(0,0,1,0))
# env variables with a correlation of <|0.7| is a good threshold for determining if predictors correlate
dev.off()

tiff('figures/EnvVariablesOnly/SSW_ScaledCentered_EnvVarOnly_AllData_CorrPlot2.tiff', width = 7, height = 7, units = 'in', res = 300)
corrplot(cor_mat.env1, p.mat = cor.all.mat$p, method = 'square', type = 'lower', insig='blank',
         addCoef.col ='white', number.cex = 0.7, order = 'alphabet', diag=FALSE, tl.cex=0.5,COL2(diverging = c("RdYlBu"), n = 200),
         title="All Timepoints",mar=c(0,0,1,0))
dev.off()

## August Corrplot
# Calculate correlations for corr coefficient & p values
cor(August.2021[,c(8,10:11,14:16)],method='pearson')
cor.aug21.mat = cor.mtest(August.2021[,c(8,10:11,14:16)],method='pearson', conf.level = 0.95)

# Visualize with a corrplot
cor_mat.env.aug <- cor(August.2021[,c(8,10:11,14:16)], method='pearson')
cor_mat.env.aug

symnum(cor_mat.env.aug)

tiff('figures/EnvVariablesOnly/SSW_ScaledCentered_EnvVarOnly_August21_CorrPlot.tiff', width = 7, height = 7, units = 'in', res = 300)
corrplot.mixed(cor_mat.env.aug, p.mat=cor.aug21.mat$p, tl.pos='lt', tl.cex=0.7, sig.level = 0.05, order="alphabet",insig='blank',number.cex=0.8,
                        diag='n',cl.ratio = 0.2, tl.srt = 45,title="August 2021 Env Variables",mar=c(0,0,1,0))
# env variables with a correlation of <|0.7| is a good threshold for determining if predictors correlate
dev.off()

tiff('figures/EnvVariablesOnly/SSW_ScaledCentered_EnvVarOnly_August21_CorrPlot2.tiff', width = 7, height = 7, units = 'in', res = 300)
corrplot(cor_mat.env.aug, p.mat = cor.aug21.mat$p, method = 'square', type = 'lower', insig='blank',
         addCoef.col ='white', number.cex = 0.7, order = 'alphabet', diag=FALSE, tl.cex=0.5,COL2(diverging = c("RdYlBu"), n = 200),
         title="August 2021",mar=c(0,0,1,0))
dev.off()

## December Corrplot
cor(December.2021[,c(8,10:11,14:16)],method='pearson')
cor.dec21.mat = cor.mtest(December.2021[,c(8,10:11,14:16)],method='pearson', conf.level = 0.95)

cor_mat.env.dec <- cor(December.2021[,c(8,10:11,14:16)], method='pearson')
cor_mat.env.dec

symnum(cor_mat.env.dec)

tiff('figures/EnvVariablesOnly/SSW_ScaledCentered_EnvVarOnly_December21_CorrPlot.tiff', width = 7, height = 7, units = 'in', res = 300)
corrplot.mixed(cor_mat.env.dec, p.mat=cor.dec21.mat$p, tl.pos='lt', tl.cex=0.7, sig.level = 0.05, order="alphabet",insig='blank',number.cex=0.8,
                        diag='n',cl.ratio = 0.2, tl.srt = 45,title="December 2021 Env Variables",mar=c(0,0,1,0))
# env variables with a correlation of <|0.7| is a good threshold for determining if predictors correlate
dev.off()

tiff('figures/EnvVariablesOnly/SSW_ScaledCentered_EnvVarOnly_December21_CorrPlot2.tiff', width = 7, height = 7, units = 'in', res = 300)
corrplot(cor_mat.env.dec, p.mat = cor.dec21.mat$p, method = 'square', type = 'lower', insig='blank',
         addCoef.col ='white', number.cex = 0.7, order = 'alphabet', diag=FALSE, tl.cex=0.5,COL2(diverging = c("RdYlBu"), n = 200),
         title="December 2021",mar=c(0,0,1,0))
dev.off()

## April Corrplot
cor(April.2022[,c(8,10:11,14:16)],method='pearson')
cor.apr22.mat = cor.mtest(April.2022[,c(8,10:11,14:16)],method='pearson', conf.level = 0.95)

cor_mat.env.apr <- cor(April.2022[,c(8,10:11,14:16)], method='pearson')
cor_mat.env.apr

symnum(cor_mat.env.apr)

tiff('figures/EnvVariablesOnly/SSW_ScaledCentered_EnvVarOnly_April22_CorrPlot.tiff', width = 7, height = 7, units = 'in', res = 300)
corrplot.mixed(cor_mat.env.apr, p.mat=cor.apr22.mat$p, tl.pos='lt', tl.cex=0.7, sig.level = 0.05, order="alphabet",insig='blank',number.cex=0.8,
               diag='n',cl.ratio = 0.2, tl.srt = 45,title="April 2022 Env Variables",mar=c(0,0,1,0))
# env variables with a correlation of <|0.7| is a good threshold for determining if predictors correlate
dev.off()

tiff('figures/EnvVariablesOnly/SSW_ScaledCentered_EnvVarOnly_April22_CorrPlot2.tiff', width = 7, height = 7, units = 'in', res = 300)
corrplot(cor_mat.env.apr, p.mat = cor.apr22.mat$p, method = 'square', type = 'lower', insig='blank',
         addCoef.col ='white', number.cex = 0.7, order = 'alphabet', diag=FALSE, tl.cex=0.5,COL2(diverging = c("RdYlBu"), n = 200),
         title="April 2022",mar=c(0,0,1,0))
dev.off()

#### Do Env Variables Correlate Across Time Points ####
# DO %
cor.test(meta_scaled$DO_Percent_Local, meta_scaled$ORP_mV, method="pearson") # ***
# r = 0.5619567, p = 0.000161 --> medium corr, & it's significant
cor.test(meta_scaled$DO_Percent_Local, meta_scaled$Temp_DegC, method="pearson") # ***
# r = -0.7973317, p-value = 7.386e-10 --> strong negative correlation, significant
cor.test(meta_scaled$DO_Percent_Local, meta_scaled$Dissolved_OrganicMatter_RFU, method="pearson") # ***
# r = -0.6131097, p-value = 2.6e-05 --> medium correlation & significant
cor.test(meta_scaled$DO_Percent_Local, meta_scaled$Depth.num, method="pearson") #
# r = -0.3850655 , p = 0.01414 --> not strong negative correlation,significant

plot(x=meta_scaled$DO_Percent_Local, y=meta_scaled$Dissolved_OrganicMatter_RFU, col=meta_scaled$SampDate_Color)
plot(x=meta_scaled$DO_Percent_Local, y=meta_scaled$ORP_mV, col=meta_scaled$SampDate_Color)
plot(x=meta_scaled$DO_Percent_Local, y=meta_scaled$Temp_DegC, col=meta_scaled$SampDate_Color)
plot(x=meta_scaled$DO_Percent_Local, y=meta_scaled$Depth.num, col=meta_scaled$SampDate_Color)

# ORP
cor.test(meta_scaled$ORP_mV, meta_scaled$Temp_DegC, method="pearson") # ***
# r = -0.513023, p-value = 0.0007117 --> medium negative correlation & significant
cor.test(meta_scaled$ORP_mV, meta_scaled$Dissolved_OrganicMatter_RFU, method="pearson") # ***
# r = -0.380499, p-value = 0.00833 --> lame negative correlation, it's significant
cor.test(meta_scaled$ORP_mV, meta_scaled$Depth.num, method="pearson")
# r = -0.2568054, p-value =  0.1097 --> not sig, not strong corr

plot(x=meta_scaled$ORP_mV, y=meta_scaled$Dissolved_OrganicMatter_RFU, col=meta_scaled$SampDate_Color)
plot(x=meta_scaled$ORP_mV, y=meta_scaled$Temp_DegC, col=meta_scaled$SampDate_Color)
plot(x=meta_scaled$ORP_mV, y=meta_scaled$Depth.num, col=meta_scaled$SampDate_Color)

# Dissolved Organic Matter
cor.test(meta_scaled$Dissolved_OrganicMatter_RFU, meta_scaled$Temp_DegC, method="pearson") # ***
# r = 0.4585448, p-value = 0.002923 # medium corr, significant
cor.test(meta_scaled$Dissolved_OrganicMatter_RFU, meta_scaled$Depth.num, method="pearson")
# r = 0.3099062 , p = 0.05165 # not sig, lame corr

plot(x=meta_scaled$Dissolved_OrganicMatter_RFU, y=meta_scaled$Temp_DegC, col=meta_scaled$SampDate_Color)
plot(x=meta_scaled$Dissolved_OrganicMatter_RFU, y=meta_scaled$Depth.num, col=meta_scaled$SampDate_Color)

# Sulfate (milliM)
cor.test(meta_scaled$Sulfate_milliM, meta_scaled$ORP_mV, method="pearson")
# r = 0.1102182, p = 0.4984 --> not a strong correlation & it's significant
cor.test(meta_scaled$Sulfate_milliM, meta_scaled$Temp_DegC, method="pearson") # ***
# r = -0.5614189, p = 0.0001639 --> strong negative correlation, significant
cor.test(meta_scaled$Sulfate_milliM, meta_scaled$DO_Percent_Local, method="pearson") # ****
# r = 0.6452846, p = 6.944e-06 --> strong correlation & significant
cor.test(meta_scaled$Sulfate_milliM, meta_scaled$Dissolved_OrganicMatter_RFU, method="pearson")
# r = -0.03455734 , p = 0.8323 --> no corr, not sig
cor.test(meta_scaled$Sulfate_milliM, meta_scaled$Sulfide_microM, method="pearson")
# r = -0.1571489 , p = 0.3328 --> no corr, not sig
cor.test(meta_scaled$Sulfate_milliM, meta_scaled$Depth.num, method="pearson")
# r = -0.1800819, p-value = 0.2662

plot(x=meta_scaled$Sulfate_milliM, y=meta_scaled$DO_Percent_Local, col=meta_scaled$SampDate_Color)
plot(x=meta_scaled$Sulfate_milliM, y=meta_scaled$Dissolved_OrganicMatter_RFU, col=meta_scaled$SampDate_Color)
plot(x=meta_scaled$Sulfate_milliM, y=meta_scaled$ORP_mV, col=meta_scaled$SampDate_Color)
plot(x=meta_scaled$Sulfate_milliM, y=meta_scaled$Temp_DegC, col=meta_scaled$SampDate_Color)
plot(x=meta_scaled$Sulfate_milliM, y=meta_scaled$Depth.num, col=meta_scaled$SampDate_Color)

# Sulfide (microM)
cor.test(meta_scaled$Sulfide_microM, meta_scaled$ORP_mV, method="pearson") # ******
# r = -0.979198, p < 2.2e-16 --> STRONG & significant negative correlation
cor.test(meta_scaled$Sulfide_microM, meta_scaled$Temp_DegC, method="pearson") # ***
# r = 0.5532523 , p = 0.0002134 --> medium correlation, significant
cor.test(meta_scaled$Sulfide_microM, meta_scaled$DO_Percent_Local, method="pearson") # ****
# r = -0.6286855, p = 1.398e-05 --> medium-strong negative correlation & significant
cor.test(meta_scaled$Sulfide_microM, meta_scaled$Dissolved_OrganicMatter_RFU, method="pearson") # ****
# r = 0.621356 , p = 1.88e-05 --> medium to strong correlation, significant
cor.test(meta_scaled$Sulfide_microM, meta_scaled$Depth.num, method="pearson")
# r = 0.2837005 , p-value = 0.07606 # not sig, no corr

plot(x=meta_scaled$Sulfide_microM, y=meta_scaled$DO_Percent_Local, col=meta_scaled$SampDate_Color)
plot(x=meta_scaled$Sulfide_microM, y=meta_scaled$Dissolved_OrganicMatter_RFU, col=meta_scaled$SampDate_Color)
plot(x=meta_scaled$Sulfide_microM, y=meta_scaled$ORP_mV, col=meta_scaled$SampDate_Color)
plot(x=meta_scaled$Sulfide_microM, y=meta_scaled$Temp_DegC, col=meta_scaled$SampDate_Color)
plot(x=meta_scaled$Sulfide_microM, y=meta_scaled$Depth.num, col=meta_scaled$SampDate_Color)

# Chlorophyll
#cor.test(meta_scaled$Chlorophyll_RFU, meta_scaled$Dissolved_OrganicMatter_RFU, method="pearson") # ****
# r = -0.70582, p-value = 3.002e-08 --> strong correlation & significant
#cor.test(meta_scaled$Chlorophyll_RFU, meta_scaled$Temp_DegC, method="pearson")

#### Do Env Variables Correlate within August 2021 ####
head(August.2021)
# create numeric depth variable
August.2021$Depth.num<-as.numeric(as.character(August.2021$Depth_m))

# DO %
cor.test(August.2021$DO_Percent_Local, August.2021$ORP_mV, method="pearson") #
# r = 0.4597742, p = 0.2517 --> little corr, not significant
cor.test(August.2021$DO_Percent_Local, August.2021$Temp_DegC, method="pearson") # ***
# r = 0.8213249, p-value = 0.01242 --> strong pos correlation, significant
cor.test(August.2021$DO_Percent_Local, August.2021$Dissolved_OrganicMatter_RFU, method="pearson") # ***
# r = -0.8891263, p-value = 0.00313 --> strong neg correlation & significant
cor.test(August.2021$DO_Percent_Local, August.2021$Depth.num, method="pearson") #
# r = -0.9295648 , p = 0.0008281 --> not strong negative correlation,significant

plot(x=August.2021$DO_Percent_Local, y=August.2021$Dissolved_OrganicMatter_RFU, col=August.2021$Depth_m) %>% text(x=August.2021$DO_Percent_Local, y=August.2021$Dissolved_OrganicMatter_RFU,labels=August.2021$Depth_m,cex=0.65, pos=3,col="black")
plot(x=August.2021$DO_Percent_Local, y=August.2021$ORP_mV, col=August.2021$Depth_m) %>% text(x=August.2021$DO_Percent_Local, y=August.2021$ORP_mV,labels=August.2021$Depth_m,cex=0.65, pos=3,col="black")
plot(x=August.2021$DO_Percent_Local, y=August.2021$Temp_DegC, col=August.2021$Depth_m) %>% text(x=August.2021$DO_Percent_Local, y=August.2021$Temp_DegC,labels=August.2021$Depth_m,cex=0.65, pos=3,col="black")
plot(x=August.2021$DO_Percent_Local, y=August.2021$Depth.num, col=August.2021$Depth_m) %>% text(x=August.2021$DO_Percent_Local, y=August.2021$Depth.num,labels=August.2021$Depth_m,cex=0.65, pos=3,col="black")

# ORP
cor.test(August.2021$ORP_mV, August.2021$Temp_DegC, method="pearson") # ***
# r = -0.06085622, p-value = 0.8862 --> not sig
cor.test(August.2021$ORP_mV, August.2021$Dissolved_OrganicMatter_RFU, method="pearson") # ***
# r = -0.8044803, p-value = 0.01605 --> strong negative, significant corr
cor.test(August.2021$ORP_mV, August.2021$Depth.num, method="pearson")
# r = -0.719755, p-value =  0.04411 --> neg cor, sig

plot(x=August.2021$ORP_mV, y=August.2021$Dissolved_OrganicMatter_RFU, col=August.2021$Depth_m) %>% text(x=August.2021$ORP_mV, y=August.2021$Dissolved_OrganicMatter_RFU,labels=August.2021$Depth_m,cex=0.65, pos=3,col="black")
plot(x=August.2021$ORP_mV, y=August.2021$Temp_DegC, col=August.2021$Depth_m) %>% text(x=August.2021$ORP_mV, y=August.2021$Temp_DegC,labels=August.2021$Depth_m,cex=0.65, pos=3,col="black")
plot(x=August.2021$ORP_mV, y=August.2021$Depth.num, col=August.2021$Depth_m) %>% text(x=August.2021$ORP_mV, y=August.2021$Depth.num,labels=August.2021$Depth_m,cex=0.65, pos=3,col="black")

# Dissolved Organic Matter
cor.test(August.2021$Dissolved_OrganicMatter_RFU, August.2021$Temp_DegC, method="pearson") # ***
# r = -0.5117475 , p-value = 0.1948 # not sig
cor.test(August.2021$Dissolved_OrganicMatter_RFU, August.2021$Depth.num, method="pearson")
# r = 0.9774045 , p = 2.835e-05 # strong sig corr

plot(x=August.2021$Dissolved_OrganicMatter_RFU, y=August.2021$Temp_DegC, col=August.2021$Depth_m) %>% text(x=August.2021$Dissolved_OrganicMatter_RFU, y=August.2021$Temp_DegC,labels=August.2021$Depth_m,cex=0.65, pos=3,col="black")
plot(x=August.2021$Dissolved_OrganicMatter_RFU, y=August.2021$Depth.num, col=August.2021$Depth_m) %>% text(x=August.2021$Dissolved_OrganicMatter_RFU, y=August.2021$Depth.num,labels=August.2021$Depth_m,cex=0.65, pos=3,col="black")

# Sulfate (milliM)
cor.test(August.2021$Sulfate_milliM, August.2021$ORP_mV, method="pearson")
# r = -0.5448982 , p = 0.1625 --> med cor, not sig
cor.test(August.2021$Sulfate_milliM, August.2021$Temp_DegC, method="pearson") # ***
# r = -0.1507523, p = 0.7216 --> strong negative correlation, significant
cor.test(August.2021$Sulfate_milliM, August.2021$DO_Percent_Local, method="pearson") # ****
# r = -0.4719281, p = 0.2377 --> strong correlation & significant
cor.test(August.2021$Sulfate_milliM, August.2021$Dissolved_OrganicMatter_RFU, method="pearson")
# r = -0.618558 , p = 0.1021 --> no corr, not sig
cor.test(August.2021$Sulfate_milliM, August.2021$Sulfide_microM, method="pearson")
# r = 0.5825559 , p = 0.1297 --> no corr, not sig
cor.test(August.2021$Sulfate_milliM, August.2021$Depth.num, method="pearson")
# r = 0.5011928, p-value = 0.2058

plot(x=August.2021$Sulfate_milliM, y=August.2021$DO_Percent_Local, col=August.2021$Depth_m) %>% text(x=August.2021$Sulfate_milliM, y=August.2021$DO_Percent_Local,labels=August.2021$Depth_m,cex=0.65, pos=3,col="black")
plot(x=August.2021$Sulfate_milliM, y=August.2021$Dissolved_OrganicMatter_RFU, col=August.2021$Depth_m) %>% text(x=August.2021$Sulfate_milliM, y=August.2021$Dissolved_OrganicMatter_RFU,labels=August.2021$Depth_m,cex=0.65, pos=3,col="black")
plot(x=August.2021$Sulfate_milliM, y=August.2021$ORP_mV, col=August.2021$Depth_m) %>% text(x=August.2021$Sulfate_milliM, y=August.2021$ORP_mV,labels=August.2021$Depth_m,cex=0.65, pos=3,col="black")
plot(x=August.2021$Sulfate_milliM, y=August.2021$Temp_DegC, col=August.2021$Depth_m) %>% text(x=August.2021$Sulfate_milliM, y=August.2021$Temp_DegC,labels=August.2021$Depth_m,cex=0.65, pos=3,col="black")
plot(x=August.2021$Sulfate_milliM, y=August.2021$Depth.num, col=August.2021$Depth_m) %>% text(x=August.2021$Sulfate_milliM, y=August.2021$Depth.num,labels=August.2021$Depth_m,cex=0.65, pos=3,col="black")

# Sulfide (microM)
cor.test(August.2021$Sulfide_microM, August.2021$ORP_mV, method="pearson") # ******
# r = -0.9708461, p = 6.06e-05 --> STRONG & significant negative correlation
cor.test(August.2021$Sulfide_microM, August.2021$Temp_DegC, method="pearson") #
# r = -0.04161607 , p = 0.9221 --> no correlation, not sig
cor.test(August.2021$Sulfide_microM, August.2021$DO_Percent_Local, method="pearson") # ****
# r = -0.5503672, p = 0.1575 --> medium negative correlation, not sig
cor.test(August.2021$Sulfide_microM, August.2021$Dissolved_OrganicMatter_RFU, method="pearson") # ****
# r = 0.8537728 , p = 0.006985 --> strong correlation, significant
cor.test(August.2021$Sulfide_microM, August.2021$Depth.num, method="pearson")
# r = 0.7801233 , p-value = 0.02239 # sig sig, medstrong corr

plot(x=August.2021$Sulfide_microM, y=August.2021$DO_Percent_Local, col=August.2021$Depth_m) %>% text(x=August.2021$Sulfide_microM, y=August.2021$DO_Percent_Local,labels=August.2021$Depth_m,cex=0.65, pos=3,col="black")
plot(x=August.2021$Sulfide_microM, y=August.2021$Dissolved_OrganicMatter_RFU, col=August.2021$Depth_m) %>% text(x=August.2021$Sulfide_microM, y=August.2021$Dissolved_OrganicMatter_RFU,labels=August.2021$Depth_m,cex=0.65, pos=3,col="black")
plot(x=August.2021$Sulfide_microM, y=August.2021$ORP_mV, col=August.2021$Depth_m) %>% text(x=August.2021$Sulfide_microM, y=August.2021$ORP_mV,labels=August.2021$Depth_m,cex=0.65, pos=3,col="black")
plot(x=August.2021$Sulfide_microM, y=August.2021$Temp_DegC, col=August.2021$Depth_m) %>% text(x=August.2021$Sulfide_microM, y=August.2021$Temp_DegC,labels=August.2021$Depth_m,cex=0.65, pos=3,col="black")
plot(x=August.2021$Sulfide_microM, y=August.2021$Depth.num, col=August.2021$Depth_m) %>% text(x=August.2021$Sulfide_microM, y=August.2021$Depth.num,labels=August.2021$Depth_m,cex=0.65, pos=3,col="black")

# Chlorophyll
#cor.test(August.2021$Chlorophyll_RFU, August.2021$Dissolved_OrganicMatter_RFU, method="pearson") # ****
# r = -0.70582, p-value = 3.002e-08 --> strong correlation & significant
#cor.test(August.2021$Chlorophyll_RFU, August.2021$Temp_DegC, method="pearson")


#### Do Env Variables Correlate within December 2021 ####
head(December.2021)
# create numeric depth variable
December.2021$Depth.num<-as.numeric(as.character(December.2021$Depth_m))

# DO %
cor.test(December.2021$DO_Percent_Local, December.2021$ORP_mV, method="pearson") #
# r = -0.4285268, p = 0.2895 --> little corr, not significant
cor.test(December.2021$DO_Percent_Local, December.2021$Temp_DegC, method="pearson") # ***
# r = 0.6685309, p-value = 0.06991 --> strong pos correlation, significant
cor.test(December.2021$DO_Percent_Local, December.2021$Dissolved_OrganicMatter_RFU, method="pearson") # ***
# r = -0.2855835, p-value = 0.4929 --> strong neg correlation & significant
cor.test(December.2021$DO_Percent_Local, December.2021$Depth.num, method="pearson") #
# r = -0.7703916 , p = 0.02529 --> not strong negative correlation,significant

plot(x=December.2021$DO_Percent_Local, y=December.2021$Dissolved_OrganicMatter_RFU, col=December.2021$Depth_m) %>% text(x=December.2021$DO_Percent_Local, y=December.2021$Dissolved_OrganicMatter_RFU,labels=December.2021$Depth_m,cex=0.65, pos=3,col="black")
plot(x=December.2021$DO_Percent_Local, y=December.2021$ORP_mV, col=December.2021$Depth_m) %>% text(x=December.2021$DO_Percent_Local, y=December.2021$ORP_mV,labels=December.2021$Depth_m,cex=0.65, pos=3,col="black")
plot(x=December.2021$DO_Percent_Local, y=December.2021$Temp_DegC, col=December.2021$Depth_m) %>% text(x=December.2021$DO_Percent_Local, y=December.2021$Temp_DegC,labels=December.2021$Depth_m,cex=0.65, pos=3,col="black")
plot(x=December.2021$DO_Percent_Local, y=December.2021$Depth.num, col=December.2021$Depth_m) %>% text(x=December.2021$DO_Percent_Local, y=December.2021$Depth.num,labels=December.2021$Depth_m,cex=0.65, pos=3,col="black")

# ORP
cor.test(December.2021$ORP_mV, December.2021$Temp_DegC, method="pearson") # ***
# r = -0.8191958, p-value = 0.01285 --> strong negative correlation & significant
cor.test(December.2021$ORP_mV, December.2021$Dissolved_OrganicMatter_RFU, method="pearson") # ***
# r = -0.5165888, p-value = 0.1899 --> mid neg corr, not sig
cor.test(December.2021$ORP_mV, December.2021$Depth.num, method="pearson")
# r = 0.7949265, p-value =  0.01838 --> midstrong corr, significant

plot(x=December.2021$ORP_mV, y=December.2021$Dissolved_OrganicMatter_RFU, col=December.2021$Depth_m) %>% text(x=December.2021$ORP_mV, y=December.2021$Dissolved_OrganicMatter_RFU,labels=December.2021$Depth_m,cex=0.65, pos=3,col="black")
plot(x=December.2021$ORP_mV, y=December.2021$Temp_DegC, col=December.2021$Depth_m) %>% text(x=December.2021$ORP_mV, y=December.2021$Temp_DegC,labels=December.2021$Depth_m,cex=0.65, pos=3,col="black")
plot(x=December.2021$ORP_mV, y=December.2021$Depth.num, col=December.2021$Depth_m) %>% text(x=December.2021$ORP_mV, y=December.2021$Depth.num,labels=December.2021$Depth_m,cex=0.65, pos=3,col="black")

# Dissolved Organic Matter
cor.test(December.2021$Dissolved_OrganicMatter_RFU, December.2021$Temp_DegC, method="pearson") # ***
# r = 0.5102602 , p-value = 0.1964 # med cor,not sig
cor.test(December.2021$Dissolved_OrganicMatter_RFU, December.2021$Depth.num, method="pearson")
# r = -0.3023965 , p = 0.4666 # no cor

plot(x=December.2021$Dissolved_OrganicMatter_RFU, y=December.2021$Temp_DegC, col=December.2021$Depth_m) %>% text(x=December.2021$Dissolved_OrganicMatter_RFU, y=December.2021$Temp_DegC,labels=December.2021$Depth_m,cex=0.65, pos=3,col="black")
plot(x=December.2021$Dissolved_OrganicMatter_RFU, y=December.2021$Depth.num, col=December.2021$Depth_m) %>% text(x=December.2021$Dissolved_OrganicMatter_RFU, y=December.2021$Depth.num,labels=December.2021$Depth_m,cex=0.65, pos=3,col="black")

# Sulfate (milliM)
cor.test(December.2021$Sulfate_milliM, December.2021$ORP_mV, method="pearson")
# r = -0.04002565 , p = 0.925 --> not sig
cor.test(December.2021$Sulfate_milliM, December.2021$Temp_DegC, method="pearson") # ***
# r = -0.08700684, p = 0.8377 not sig
cor.test(December.2021$Sulfate_milliM, December.2021$DO_Percent_Local, method="pearson") # ****
# r = -0.05199845, p = 0.9027 --> not sig
cor.test(December.2021$Sulfate_milliM, December.2021$Dissolved_OrganicMatter_RFU, method="pearson")
# r = -0.07049322 , p = 0.8683 --> not sig
cor.test(December.2021$Sulfate_milliM, December.2021$Sulfide_microM, method="pearson")
# r = 0.2142547 , p = 0.6104 --> not sig
cor.test(December.2021$Sulfate_milliM, December.2021$Depth.num, method="pearson")
# r = 0.2767267, p-value = 0.507 --> not sig

plot(x=December.2021$Sulfate_milliM, y=December.2021$DO_Percent_Local, col=December.2021$Depth_m) %>% text(x=December.2021$Sulfate_milliM, y=December.2021$DO_Percent_Local,labels=December.2021$Depth_m,cex=0.65, pos=3,col="black")
plot(x=December.2021$Sulfate_milliM, y=December.2021$Dissolved_OrganicMatter_RFU, col=December.2021$Depth_m) %>% text(x=December.2021$Sulfate_milliM, y=December.2021$Dissolved_OrganicMatter_RFU,labels=December.2021$Depth_m,cex=0.65, pos=3,col="black")
plot(x=December.2021$Sulfate_milliM, y=December.2021$ORP_mV, col=December.2021$Depth_m) %>% text(x=December.2021$Sulfate_milliM, y=December.2021$ORP_mV,labels=December.2021$Depth_m,cex=0.65, pos=3,col="black")
plot(x=December.2021$Sulfate_milliM, y=December.2021$Temp_DegC, col=December.2021$Depth_m) %>% text(x=December.2021$Sulfate_milliM, y=December.2021$Temp_DegC,labels=December.2021$Depth_m,cex=0.65, pos=3,col="black")
plot(x=December.2021$Sulfate_milliM, y=December.2021$Depth.num, col=December.2021$Depth_m) %>% text(x=December.2021$Sulfate_milliM, y=December.2021$Depth.num,labels=December.2021$Depth_m,cex=0.65, pos=3,col="black")

# Sulfide (microM)
cor.test(December.2021$Sulfide_microM, December.2021$ORP_mV, method="pearson")
# r = 0.5933017, p = 0.121 --> not sig
cor.test(December.2021$Sulfide_microM, December.2021$Temp_DegC, method="pearson")
# r = -0.529703 , p = 0.177 --> no correlation, not sig
cor.test(December.2021$Sulfide_microM, December.2021$DO_Percent_Local, method="pearson")
# r = -0.08996701, p = 0.1575 --> medium negative correlation, not sig
cor.test(December.2021$Sulfide_microM, December.2021$Dissolved_OrganicMatter_RFU, method="pearson")
# r = 0.8537728 , p = 0.8322 --> no cor, not sig
cor.test(December.2021$Sulfide_microM, December.2021$Depth.num, method="pearson")
# r = 0.593074 , p-value = 0.1212 # not sig, medstrong corr

plot(x=December.2021$Sulfide_microM, y=December.2021$DO_Percent_Local, col=December.2021$Depth_m) %>% text(x=December.2021$Sulfide_microM, y=December.2021$DO_Percent_Local,labels=December.2021$Depth_m,cex=0.65, pos=3,col="black")
plot(x=December.2021$Sulfide_microM, y=December.2021$Dissolved_OrganicMatter_RFU, col=December.2021$Depth_m) %>% text(x=December.2021$Sulfide_microM, y=December.2021$Dissolved_OrganicMatter_RFU,labels=December.2021$Depth_m,cex=0.65, pos=3,col="black")
plot(x=December.2021$Sulfide_microM, y=December.2021$ORP_mV, col=December.2021$Depth_m) %>% text(x=December.2021$Sulfide_microM, y=December.2021$ORP_mV,labels=December.2021$Depth_m,cex=0.65, pos=3,col="black")
plot(x=December.2021$Sulfide_microM, y=December.2021$Temp_DegC, col=December.2021$Depth_m) %>% text(x=December.2021$Sulfide_microM, y=December.2021$Temp_DegC,labels=December.2021$Depth_m,cex=0.65, pos=3,col="black")
plot(x=December.2021$Sulfide_microM, y=December.2021$Depth.num, col=December.2021$Depth_m) %>% text(x=December.2021$Sulfide_microM, y=December.2021$Depth.num,labels=December.2021$Depth_m,cex=0.65, pos=3,col="black")

# Chlorophyll
#cor.test(December.2021$Chlorophyll_RFU, December.2021$Dissolved_OrganicMatter_RFU, method="pearson") # ****
# r = -0.70582, p-value = 3.002e-08 --> strong correlation & significant
#cor.test(December.2021$Chlorophyll_RFU, December.2021$Temp_DegC, method="pearson")


#### Do Env Variables Correlate within April 2022 ####
head(April.2022)
# create numeric depth variable
April.2022$Depth.num<-as.numeric(as.character(April.2022$Depth_m))

# DO %
cor.test(April.2022$DO_Percent_Local, April.2022$ORP_mV, method="pearson") #
# r = 0.8447182, p = 0.008304 --> strong, sig
cor.test(April.2022$DO_Percent_Local, April.2022$Temp_DegC, method="pearson") # ***
# r = 0.916565, p-value = 0.001363 --> strong, sig
cor.test(April.2022$DO_Percent_Local, April.2022$Dissolved_OrganicMatter_RFU, method="pearson") # ***
# r = -0.9415401, p-value = 0.0004778 --> strong neg corr, sig
cor.test(April.2022$DO_Percent_Local, April.2022$Depth.num, method="pearson") #
# r = -0.8727097 , p = 0.004676 -->

plot(x=April.2022$DO_Percent_Local, y=April.2022$Dissolved_OrganicMatter_RFU, col=April.2022$Depth_m) %>% text(x=April.2022$DO_Percent_Local, y=April.2022$Dissolved_OrganicMatter_RFU,labels=April.2022$Depth_m,cex=0.65, pos=3,col="black")
plot(x=April.2022$DO_Percent_Local, y=April.2022$ORP_mV, col=April.2022$Depth_m) %>% text(x=April.2022$DO_Percent_Local, y=April.2022$ORP_mV,labels=April.2022$Depth_m,cex=0.65, pos=3,col="black")
plot(x=April.2022$DO_Percent_Local, y=April.2022$Temp_DegC, col=April.2022$Depth_m) %>% text(x=April.2022$DO_Percent_Local, y=April.2022$Temp_DegC,labels=April.2022$Depth_m,cex=0.65, pos=3,col="black")
plot(x=April.2022$DO_Percent_Local, y=April.2022$Depth.num, col=April.2022$Depth_m) %>% text(x=April.2022$DO_Percent_Local, y=April.2022$Depth.num,labels=April.2022$Depth_m,cex=0.65, pos=3,col="black")

# ORP
cor.test(April.2022$ORP_mV, April.2022$Temp_DegC, method="pearson") # ***
# r = 0.9457778, p-value = 0.0003825 --> strong cor, sig
cor.test(April.2022$ORP_mV, April.2022$Dissolved_OrganicMatter_RFU, method="pearson") # ***
# r = -0.9428313, p-value = 0.0004473 --> strong neg corr, sig
cor.test(April.2022$ORP_mV, April.2022$Depth.num, method="pearson")
# r = -0.6067133, p-value =  0.1107 --> med strong neg corr, not sig

plot(x=April.2022$ORP_mV, y=April.2022$Dissolved_OrganicMatter_RFU, col=April.2022$Depth_m) %>% text(x=April.2022$ORP_mV, y=April.2022$Dissolved_OrganicMatter_RFU,labels=April.2022$Depth_m,cex=0.65, pos=3,col="black")
plot(x=April.2022$ORP_mV, y=April.2022$Temp_DegC, col=April.2022$Depth_m) %>% text(x=April.2022$ORP_mV, y=April.2022$Temp_DegC,labels=April.2022$Depth_m,cex=0.65, pos=3,col="black")
plot(x=April.2022$ORP_mV, y=April.2022$Depth.num, col=April.2022$Depth_m) %>% text(x=April.2022$ORP_mV, y=April.2022$Depth.num,labels=April.2022$Depth_m,cex=0.65, pos=3,col="black")

# Dissolved Organic Matter
cor.test(April.2022$Dissolved_OrganicMatter_RFU, April.2022$Temp_DegC, method="pearson") # ***
# r = -0.9964688 , p-value = 1.098e-07 --> strong neg corr, sig but not real (see plot)
cor.test(April.2022$Dissolved_OrganicMatter_RFU, April.2022$Depth.num, method="pearson")
# r = 0.7382372 , p = 0.0365 # --> strong cor, sig

plot(x=April.2022$Dissolved_OrganicMatter_RFU, y=April.2022$Temp_DegC, col=April.2022$Depth_m) %>% text(x=April.2022$Dissolved_OrganicMatter_RFU, y=April.2022$Temp_DegC,labels=April.2022$Depth_m,cex=0.65, pos=3,col="black")
plot(x=April.2022$Dissolved_OrganicMatter_RFU, y=April.2022$Depth.num, col=April.2022$Depth_m) %>% text(x=April.2022$Dissolved_OrganicMatter_RFU, y=April.2022$Depth.num,labels=April.2022$Depth_m,cex=0.65, pos=3,col="black")

# Sulfate (milliM)
cor.test(April.2022$Sulfate_milliM, April.2022$ORP_mV, method="pearson")
# r = 0.4050627 , p = 0.3195 -->
cor.test(April.2022$Sulfate_milliM, April.2022$Temp_DegC, method="pearson") # ***
# r = 0.3890084, p = 0.3409 -->
cor.test(April.2022$Sulfate_milliM, April.2022$DO_Percent_Local, method="pearson") # ****
# r = 0.5755169, p = 0.1355 -->
cor.test(April.2022$Sulfate_milliM, April.2022$Dissolved_OrganicMatter_RFU, method="pearson")
# r = -0.4332544 , p = 0.2836 -->
cor.test(April.2022$Sulfate_milliM, April.2022$Sulfide_microM, method="pearson")
# r = -0.3599731 , p = 0.3811 -->
cor.test(April.2022$Sulfate_milliM, April.2022$Depth.num, method="pearson")
# r = -0.6734952, p-value = 0.0671

plot(x=April.2022$Sulfate_milliM, y=April.2022$DO_Percent_Local, col=April.2022$Depth_m) %>% text(x=April.2022$Sulfate_milliM, y=April.2022$DO_Percent_Local,labels=April.2022$Depth_m,cex=0.65, pos=3,col="black")
plot(x=April.2022$Sulfate_milliM, y=April.2022$Dissolved_OrganicMatter_RFU, col=April.2022$Depth_m) %>% text(x=April.2022$Sulfate_milliM, y=April.2022$Dissolved_OrganicMatter_RFU,labels=April.2022$Depth_m,cex=0.65, pos=3,col="black")
plot(x=April.2022$Sulfate_milliM, y=April.2022$ORP_mV, col=April.2022$Depth_m) %>% text(x=April.2022$Sulfate_milliM, y=April.2022$ORP_mV,labels=April.2022$Depth_m,cex=0.65, pos=3,col="black")
plot(x=April.2022$Sulfate_milliM, y=April.2022$Temp_DegC, col=April.2022$Depth_m) %>% text(x=April.2022$Sulfate_milliM, y=April.2022$Temp_DegC,labels=April.2022$Depth_m,cex=0.65, pos=3,col="black")
plot(x=April.2022$Sulfate_milliM, y=April.2022$Depth.num, col=April.2022$Depth_m) %>% text(x=April.2022$Sulfate_milliM, y=April.2022$Depth.num,labels=April.2022$Depth_m,cex=0.65, pos=3,col="black")

# Sulfide (microM)
cor.test(April.2022$Sulfide_microM, April.2022$ORP_mV, method="pearson") #
# r = -0.04204334, p = 0.9213 -->
cor.test(April.2022$Sulfide_microM, April.2022$Temp_DegC, method="pearson") #
# r = -0.2711869 , p = 0.5159 -->
cor.test(April.2022$Sulfide_microM, April.2022$DO_Percent_Local, method="pearson") #
# r = -0.4560053, p = 0.2561 -->
cor.test(April.2022$Sulfide_microM, April.2022$Dissolved_OrganicMatter_RFU, method="pearson") #
# r = 0.2883266 , p = 0.4886 -->
cor.test(April.2022$Sulfide_microM, April.2022$Depth.num, method="pearson")
# r = 0.4457843 , p-value = 0.2683 #

plot(x=April.2022$Sulfide_microM, y=April.2022$DO_Percent_Local, col=April.2022$Depth_m) %>% text(x=April.2022$Sulfide_microM, y=April.2022$DO_Percent_Local,labels=April.2022$Depth_m,cex=0.65, pos=3,col="black")
plot(x=April.2022$Sulfide_microM, y=April.2022$Dissolved_OrganicMatter_RFU, col=April.2022$Depth_m) %>% text(x=April.2022$Sulfide_microM, y=April.2022$Dissolved_OrganicMatter_RFU,labels=April.2022$Depth_m,cex=0.65, pos=3,col="black")
plot(x=April.2022$Sulfide_microM, y=April.2022$ORP_mV, col=April.2022$Depth_m) %>% text(x=April.2022$Sulfide_microM, y=April.2022$ORP_mV,labels=April.2022$Depth_m,cex=0.65, pos=3,col="black")
plot(x=April.2022$Sulfide_microM, y=April.2022$Temp_DegC, col=April.2022$Depth_m) %>% text(x=April.2022$Sulfide_microM, y=April.2022$Temp_DegC,labels=April.2022$Depth_m,cex=0.65, pos=3,col="black")
plot(x=April.2022$Sulfide_microM, y=April.2022$Depth.num, col=April.2022$Depth_m) %>% text(x=April.2022$Sulfide_microM, y=April.2022$Depth.num,labels=April.2022$Depth_m,cex=0.65, pos=3,col="black")

# Chlorophyll
#cor.test(April.2022$Chlorophyll_RFU, April.2022$Dissolved_OrganicMatter_RFU, method="pearson") # ****
# r = -0.70582, p-value = 3.002e-08 --> strong correlation & significant
#cor.test(April.2022$Chlorophyll_RFU, April.2022$Temp_DegC, method="pearson")

#### Plots of Env Variables ####

# Compare all variables across Depths
dep.dom<-ggplot(metadata, aes(x=Depth_m, y=Dissolved_OrganicMatter_RFU,color=SampDate,group=SampDate)) +   geom_point(size=5) + geom_line(linewidth=1) + theme_bw()+
  labs(title="Dissolved Organic Matter (DOM) by Depth & Sample Date",subtitle="Using Raw DOM RFU Data",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(metadata$SampDate_Color[order(metadata$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Depth (m)") + ylab("DOM (RFU)")+coord_flip()+ scale_x_discrete(limits=rev)

ggsave(dep.dom,filename = "figures/EnvVariablesOnly/SSW_DOM_Depth_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)

dep.orp<-ggplot(metadata, aes(x=Depth_m, y=ORP_mV,color=SampDate,group=SampDate)) +   geom_point(size=5) + geom_line(linewidth=1) + theme_bw()+
  labs(title="Oxidative-Reduction Potential by Depth & Sample Date",subtitle="Using Raw ORP (mV) Data",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(metadata$SampDate_Color[order(metadata$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Depth (m)") + ylab("ORP (mV)")+coord_flip()+ scale_x_discrete(limits=rev)

ggsave(dep.orp,filename = "figures/EnvVariablesOnly/SSW_ORP_Depth_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)

dep.sulf<-ggplot(metadata, aes(x=Depth_m, y=Sulfate_milliM,color=SampDate,group=SampDate)) +   geom_point(size=5) + geom_line(linewidth=1) + theme_bw()+
  labs(title="Sulfate by Depth & Sample Date",subtitle="Using Raw Sulfate (milliM) Data",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(metadata$SampDate_Color[order(metadata$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Depth (m)") + ylab("Sulfate (milliM)")+coord_flip()+ scale_x_discrete(limits=rev)

ggsave(dep.sulf,filename = "figures/EnvVariablesOnly/SSW_Sulfate_Depth_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)

dep.pdo<-ggplot(metadata, aes(x=Depth_m, y=DO_Percent_Local,color=SampDate,group=SampDate)) +   geom_point(size=5) + geom_line(linewidth=1) + theme_bw()+
  labs(title="Dissolved Oxygen by Depth & Sample Date",subtitle="Using Raw DO (%) Data",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(metadata$SampDate_Color[order(metadata$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Depth (m)") + ylab("DO (%)")+coord_flip()+ scale_x_discrete(limits=rev)

ggsave(dep.pdo,filename = "figures/EnvVariablesOnly/SSW_PercentDO_Depth_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)

dep.do<-ggplot(metadata, aes(x=Depth_m, y=Dissolved_Oxygen_mgL,color=SampDate,group=SampDate)) +   geom_point(size=5) + geom_line(linewidth=1) + theme_bw()+
  labs(title="Dissolved Oxygen by Depth & Sample Date",subtitle="Using Raw DO (mg/L) Data",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(metadata$SampDate_Color[order(metadata$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Depth (m)") + ylab("DO (mg/L)")+coord_flip()+ scale_x_discrete(limits=rev)

ggsave(dep.do,filename = "figures/EnvVariablesOnly/SSW_DO_mgL_Depth_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)

dep.hs<-ggplot(metadata, aes(x=Depth_m, y=Sulfide_microM,color=SampDate,group=SampDate)) +   geom_point(size=5) + geom_line(linewidth=1) + theme_bw()+
  labs(title="Sulfide by Depth & Sample Date",subtitle="Using Raw Sulfide (microM) Data",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(metadata$SampDate_Color[order(metadata$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Depth (m)") + ylab("Sulfide (microM)")+coord_flip()+ scale_x_discrete(limits=rev)

ggsave(dep.hs,filename = "figures/EnvVariablesOnly/SSW_Sulfide_Depth_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)

dep.chlr<-ggplot(metadata, aes(x=Depth_m, y=Chlorophyll_RFU,color=SampDate,group=SampDate)) +   geom_point(size=5) + geom_line(linewidth=1) + theme_bw()+
  labs(title="Chlorophyll by Depth & Sample Date",subtitle="Using Raw Chlorophyll (RFU) Data",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(metadata$SampDate_Color[order(metadata$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Depth (m)") + ylab("Chlorophyll (RFU)")+coord_flip()+ scale_x_discrete(limits=rev)

ggsave(dep.chlr,filename = "figures/EnvVariablesOnly/SSW_Chlorophyll_Depth_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)

dep.temp<-ggplot(metadata, aes(x=Depth_m, y=Temp_DegC,color=SampDate,group=SampDate)) +   geom_point(size=5) + geom_line(linewidth=1) + theme_bw()+
  labs(title="Temperature by Depth & Sample Date",subtitle="Using Raw Data",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(metadata$SampDate_Color[order(metadata$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Depth (m)") + ylab("Temperature (C)")+coord_flip()+ scale_x_discrete(limits=rev)

ggsave(dep.temp,filename = "figures/EnvVariablesOnly/SSW_Temp_Depth_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)

dep.turb<-ggplot(metadata, aes(x=Depth_m, y=Turbidity_FNU,color=SampDate,group=SampDate)) +   geom_point(size=5) + geom_line(linewidth=1) + theme_bw()+
  labs(title="Turbidity by Depth & Sample Date",subtitle="Using Raw Data",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(metadata$SampDate_Color[order(metadata$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Depth (m)") + ylab("Turbidity (FNU)")+coord_flip()+ scale_x_discrete(limits=rev)

ggsave(dep.turb,filename = "figures/EnvVariablesOnly/SSW_Turbidity_Depth_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)

dep.chloro<-ggplot(metadata, aes(x=Depth_m, y=Chlorophyll_RFU,color=SampDate,group=SampDate)) +   geom_point(size=5) + geom_line(linewidth=1) + theme_bw()+
  labs(title="Chlorophyll by Depth & Sample Date",subtitle="Using Raw Data",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(metadata$SampDate_Color[order(metadata$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Depth (m)") + ylab("Chlorophyll (RNU)")+coord_flip()+ scale_x_discrete(limits=rev)

ggsave(dep.chloro,filename = "figures/EnvVariablesOnly/SSW_Chlorophyll_Depth_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)

dep.sal<-ggplot(meta_salinity, aes(x=Depth_m, y=Salinity_ppt,color=SampDate,group=SampDate)) +   geom_point(size=5) + geom_line(linewidth=1) + theme_bw()+
  labs(title="Salinity by Depth & Sample Date",subtitle="Using Raw Data",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(metadata$SampDate_Color[order(metadata$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Depth (m)") + ylab("Salinity (ppt)")+coord_flip()+ scale_x_discrete(limits=rev)

ggsave(dep.sal,filename = "figures/EnvVariablesOnly/SSW_Salinity_Depth_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)

# Compare variables to each other
dom.hs<-ggplot(metadata, aes(x=Dissolved_OrganicMatter_RFU, y=Sulfide_microM,color=SampDate,group=SampDate)) +   geom_point(size=5) + geom_line(linewidth=1) + theme_bw()+
  labs(title="Sulfide ~ DOM",subtitle="Using Raw Data",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(metadata$SampDate_Color[order(metadata$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("DOM (RFU)") + ylab("Sulfide (microM)")

ggsave(dom.hs,filename = "figures/EnvVariablesOnly/SSW_Sulfide_DOM_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)

dom.so4<-ggplot(metadata, aes(x=Dissolved_OrganicMatter_RFU, y=Sulfate_milliM,color=SampDate,group=SampDate)) +   geom_point(size=5) + geom_line(linewidth=1) + theme_bw()+
  labs(title="Sulfate ~ DOM",subtitle="Using Raw Data",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(metadata$SampDate_Color[order(metadata$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("DOM (RFU)") + ylab("Sulfate (milliM)")

ggsave(dom.so4,filename = "figures/EnvVariablesOnly/SSW_Sulfate_DOM_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)

do.hs<-ggplot(metadata, aes(x=DO_Percent_Local, y=Sulfide_microM,color=SampDate,group=SampDate)) +   geom_point(size=5) + geom_line(linewidth=1) + theme_bw()+
  labs(title="Sulfide ~ DO%",subtitle="Using Raw Data",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(metadata$SampDate_Color[order(metadata$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("DO (%)") + ylab("Sulfide (microM)")

ggsave(do.hs,filename = "figures/EnvVariablesOnly/SSW_Sulfide_DO.Percent_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)

do.so4<-ggplot(metadata, aes(x=DO_Percent_Local, y=Sulfate_milliM,color=SampDate,group=SampDate)) +   geom_point(size=5) + geom_line(linewidth=1) + theme_bw()+
  labs(title="Sulfate ~ DO%",subtitle="Using Raw Data",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(metadata$SampDate_Color[order(metadata$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("DO (%)") + ylab("Sulfate (milliM)")

ggsave(do.so4,filename = "figures/EnvVariablesOnly/SSW_Sulfate_DO.Percent_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)

orp.hs<-ggplot(metadata, aes(x=ORP_mV, y=Sulfide_microM,color=SampDate,group=SampDate)) +   geom_point(size=5) + geom_line(linewidth=1) + theme_bw()+
  labs(title="Sulfide ~ ORP (mV)",subtitle="Using Raw Data",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(metadata$SampDate_Color[order(metadata$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("ORP (mV") + ylab("Sulfide (microM)")

ggsave(orp.hs,filename = "figures/EnvVariablesOnly/SSW_Sulfide_ORP_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)

orp.dom<-ggplot(metadata, aes(x=ORP_mV, y=Dissolved_OrganicMatter_RFU,color=SampDate,group=SampDate)) +   geom_point(size=5) + geom_line(linewidth=1) + theme_bw()+
  labs(title="DOM (RFU) ~ ORP (mV)",subtitle="Using Raw Data",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(metadata$SampDate_Color[order(metadata$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("DOM (RFU") + ylab("Sulfide (microM)")

ggsave(orp.dom,filename = "figures/EnvVariablesOnly/SSW_DOM_ORP_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)

orp.do<-ggplot(metadata, aes(x=ORP_mV, y=DO_Percent_Local,color=SampDate,group=SampDate)) +   geom_point(size=5) + geom_line(linewidth=1) + theme_bw()+
  labs(title="DO% ~ ORP (mV)",subtitle="Using Raw Data",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(metadata$SampDate_Color[order(metadata$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("ORP (mV") + ylab("DO %")

ggsave(orp.do,filename = "figures/EnvVariablesOnly/SSW_DO.Percent_ORP_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)

temp.do<-ggplot(metadata, aes(x=Temp_DegC, y=DO_Percent_Local,color=SampDate,group=SampDate)) +   geom_point(size=5) + geom_line(linewidth=1) + theme_bw()+
  labs(title="DO% ~ Temperature (C)",subtitle="Using Raw Data",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(metadata$SampDate_Color[order(metadata$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("Temp (C") + ylab("DO %")

ggsave(temp.do,filename = "figures/EnvVariablesOnly/SSW_DO.Percent_Temp_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)

temp.dom<-ggplot(metadata, aes(x=Temp_DegC, y=Dissolved_OrganicMatter_RFU,color=SampDate,group=SampDate)) +   geom_point(size=5) + geom_line(linewidth=1) + theme_bw()+
  labs(title="DOM (RFU) ~ Temperature (C)",subtitle="Using Raw Data",color="Sample Date")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(metadata$SampDate_Color[order(metadata$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("DOM (RFU") + ylab("Sulfide (microM)")

ggsave(temp.dom,filename = "figures/EnvVariablesOnly/SSW_DOM_ORP_bySampleDate_scatterplot.png", width=12, height=10, dpi=600)

#### Plot Correlations by Sample Date to Check Correlations Above ####

# August
aug1<-ggplot(August.2021, aes(x=Temp_DegC, y=DO_Percent_Local,color=Depth_m)) + geom_point(size=5) + theme_bw()+
  labs(title="Temperature & DO% - August 2021",subtitle="Using Centered & Scaled Data",color="Depth (m)")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  xlab("Temp (C)") + ylab("DO (%)")
ggsave(aug1,filename = "figures/EnvVariablesOnly/Within_SampDates/SSW_DO.Percent_Temp_August2021_scatterplot.png", width=12, height=10, dpi=600)

aug2<-ggplot(August.2021, aes(x=Dissolved_OrganicMatter_RFU, y=DO_Percent_Local,color=Depth_m)) + geom_point(size=5) + theme_bw()+
  labs(title="DOM & DO% - August 2021",subtitle="Using Centered & Scaled Data",color="Depth (m)")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  xlab("DOM (RFU)") + ylab("DO (%)")
ggsave(aug2,filename = "figures/EnvVariablesOnly/Within_SampDates/SSW_DO.Percent_DOM_August2021_scatterplot.png", width=12, height=10, dpi=600)

aug3<-ggplot(August.2021, aes(x=Dissolved_OrganicMatter_RFU, y=ORP_mV,color=Depth_m)) + geom_point(size=5) + theme_bw()+
  labs(title="DOM & ORP - August 2021",subtitle="Using Centered & Scaled Data",color="Depth (m)")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  xlab("DOM (RFU)") + ylab("ORP (mV)")
ggsave(aug3,filename = "figures/EnvVariablesOnly/Within_SampDates/SSW_ORP_DOM_August2021_scatterplot.png", width=12, height=10, dpi=600)

aug4<-ggplot(August.2021, aes(x=Dissolved_OrganicMatter_RFU, y=Sulfide_microM,color=Depth_m)) + geom_point(size=5) + theme_bw()+
  labs(title="DOM & Sulfide - August 2021",subtitle="Using Centered & Scaled Data",color="Depth (m)")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  xlab("DOM (RFU)") + ylab("Sulfide (microM)")
ggsave(aug4,filename = "figures/EnvVariablesOnly/Within_SampDates/SSW_H2S_DOM_August2021_scatterplot.png", width=12, height=10, dpi=600)

aug5<-ggplot(August.2021, aes(x=ORP_mV, y=Sulfide_microM,color=Depth_m)) + geom_point(size=5) + theme_bw()+
  labs(title="ORP & Sulfide - August 2021",subtitle="Using Centered & Scaled Data",color="Depth (m)")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  xlab("ORP (mV)") + ylab("Sulfide (microM)")
ggsave(aug5,filename = "figures/EnvVariablesOnly/Within_SampDates/SSW_H2S_ORP_August2021_scatterplot.png", width=12, height=10, dpi=600)

## December
dec1<-ggplot(December.2021, aes(x=Temp_DegC, y=ORP_mV,color=Depth_m)) + geom_point(size=5) + theme_bw()+
  labs(title="Temperature & ORP - December 2021",subtitle="Using Centered & Scaled Data",color="Depth (m)")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  xlab("Temp (C)") + ylab("ORP (mV)")
ggsave(dec1,filename = "figures/EnvVariablesOnly/Within_SampDates/SSW_ORP_Temp_December2021_scatterplot.png", width=12, height=10, dpi=600)

## April
apr1<-ggplot(April.2022, aes(x=Temp_DegC, y=ORP_mV,color=Depth_m)) + geom_point(size=5) + theme_bw()+
  labs(title="Temperature & ORP - April 2022",subtitle="Using Centered & Scaled Data",color="Depth (m)")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  xlab("Temp (C)") + ylab("ORP (mV)")
ggsave(apr1,filename = "figures/EnvVariablesOnly/Within_SampDates/SSW_ORP_Temp_April2022_scatterplot.png", width=12, height=10, dpi=600)

apr2<-ggplot(April.2022, aes(x=Temp_DegC, y=DO_Percent_Local,color=Depth_m)) + geom_point(size=5) + theme_bw()+
  labs(title="Temperature & DO - April 2022",subtitle="Using Centered & Scaled Data",color="Depth (m)")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  xlab("Temp (C)") + ylab("DO%")
ggsave(apr2,filename = "figures/EnvVariablesOnly/Within_SampDates/SSW_DO_Temp_April2022_scatterplot.png", width=12, height=10, dpi=600)

apr3<-ggplot(April.2022, aes(x=Temp_DegC, y=Dissolved_OrganicMatter_RFU,color=Depth_m)) + geom_point(size=5) + theme_bw()+
  labs(title="Temperature & DOM - April 2022",subtitle="Using Centered & Scaled Data",color="Depth (m)")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  xlab("Temp (C)") + ylab("DOM (RFU)")
ggsave(apr3,filename = "figures/EnvVariablesOnly/Within_SampDates/SSW_DOM_Temp_April2022_scatterplot.png", width=12, height=10, dpi=600)

apr4<-ggplot(April.2022, aes(x=DO_Percent_Local, y=Dissolved_OrganicMatter_RFU,color=Depth_m)) + geom_point(size=5) + theme_bw()+
  labs(title="DO & DOM - April 2022",subtitle="Using Centered & Scaled Data",color="Depth (m)")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  xlab("DO%") + ylab("DOM (RFU)")
ggsave(apr4,filename = "figures/EnvVariablesOnly/Within_SampDates/SSW_DOM_DO_April2022_scatterplot.png", width=12, height=10, dpi=600)

apr5<-ggplot(April.2022, aes(x=ORP_mV, y=Dissolved_OrganicMatter_RFU,color=Depth_m)) + geom_point(size=5) + theme_bw()+
  labs(title="DOM & ORP - April 2022",subtitle="Using Centered & Scaled Data",color="Depth (m)")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  xlab("DOM") + ylab("ORP (mV)")
ggsave(apr5,filename = "figures/EnvVariablesOnly/Within_SampDates/SSW_DOM_ORP_April2022_scatterplot.png", width=12, height=10, dpi=600)

apr6<-ggplot(April.2022, aes(x=ORP_mV, y=DO_Percent_Local,color=Depth_m)) + geom_point(size=5) + theme_bw()+
  labs(title="DO & ORP - April 2022",subtitle="Using Centered & Scaled Data",color="Depth (m)")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  xlab("DO%") + ylab("ORP (mV)")
ggsave(apr6,filename = "figures/EnvVariablesOnly/Within_SampDates/SSW_DOM_ORP_April2022_scatterplot.png", width=12, height=10, dpi=600)

