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
  library(ape)
  library(plyr)
  library(dplyr)
  library(viridis)
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
  library(devtools)
  library(decontam)
  library(ggvegan)
  library(microbiome)
})

#### Load Data & See Info About Data ####
#load("data/Metagenomes/Analysis/SSW_mgm_analysis.Rdata") # load Rdata to global env
load("data/Metagenomes/Analysis/SSW_MGM_FxnBetaDiv.Rdata") # contains separated CLR by function & pathways
load("data/SSW_MGM_Functions_EnvDriver.Rdata")

head(meta_scaled)
arsen.fxns[1:4,]
ko.cov.sum_table[1:4,1:4]
head(mgm.clr.ars)
mgm.clr[1:4,1:4]

# ABOUT THE DATA:
# Before transformations (i.e., VST, CLR, etc) were done, the following was performed:
# featureCounts counted reads that mapped to genes in contigs
# Reads mapped to genes were divided by gene length for all genes across all samples
# Gene coverage was then added together for each KO ID, since multiple genes were assigned the same KO ID
# Summed coverage per KO was then transformed via median-ratio, vst, and clr

## For pathway analyses -- after gene coverage was calculated and added together per KO ID, they were added together for each pathway
## summed coverages per KO ID, then per pathway were transformed by CLR

#### Separate Sulfur Fxns by Timepoints ####
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

# matching data with user defined function -- here is the function, must run to store function in Global env
match_dat<-function(compdata, subset_metadata){
  subset_comp_data = pullrow<-(is.element(row.names(compdata), row.names(subset_metadata)))
  ### * comp data and metadata need to have row names - rownames should be Sample IDs
  subset_comp_data=compdata[pullrow,]
  return(subset_comp_data)
}

# double check that our data frames are ready for this function, aka that they both have the same rownames
## row #s do not have to be the same, but their row names should be in the same format and be able to match up
rownames(clr.cov.sum.sulf.ko)
rownames(August.2021)

# run the function
S.KO.clr_AUG21<-match_dat(clr.cov.sum.sulf.ko,August.2021)
S.KO.clr_DEC21<-match_dat(clr.cov.sum.sulf.ko,December.2021)
S.KO.clr_APR22<-match_dat(clr.cov.sum.sulf.ko,April.2022)

# did the function work the way we wanted it to?

S.KO.clr_AUG21[1:3,1:3]
rownames(August.2021) %in% rownames(S.KO.clr_AUG21) # hopefully all of the rownames match, aka will get output of TRUE

#### Check Sulfur Fxn Data Relationship w/ Env Variables (w/ DCA) ####
## remember, CCA assumes that our species have a unimodal relationship with our variables.
### unimodal = one maximum, think upsidedown bellcurve or something
## RDA assumes a linear relationship
## check the assumption

# ALL data
# add pseudocount so row sums are > 0

S.clr.pseudo<-clr.cov.sum.sulf.ko[,!names(clr.cov.sum.sulf.ko) %in% c("SampleID")]+2
S.ko.dca = decorana(S.clr.pseudo)

#plot(S.ko.dca) # may take too long to load, do not run unless you have to
S.ko.dca #DCA1 axis length = 0.159106; use RDA
## The length of first DCA axis:
## > 4 indicates heterogeneous dataset on which unimodal methods should be used (CCA),
##  < 3 indicates homogeneous dataset for which linear methods are suitable (RDA)
## between 3 and 4 both linear and unimodal methods are OK.

# BY MONTH

S.KO.clr_A21.pseudo<-S.KO.clr_AUG21[,!names(S.KO.clr_AUG21) %in% c("SampleID")]+1
S.ko.A21.dca = decorana(S.KO.clr_A21.pseudo)
S.ko.A21.dca #DCA1 axis length = 0.201019; use RDA

S.KO.clr_D21.pseudo<-S.KO.clr_DEC21[,!names(S.KO.clr_DEC21) %in% c("SampleID")]+2
S.ko.D21.dca = decorana(S.KO.clr_D21.pseudo)
S.ko.D21.dca #DCA1 axis length = 0.093706; use RDA

S.KO.clr_A22.pseudo<-S.KO.clr_APR22[,!names(S.KO.clr_APR22) %in% c("SampleID")]+1
S.ko.A22.dca = decorana(S.KO.clr_A22.pseudo)
S.ko.A22.dca #DCA1 axis length = 0.172273; use RDA

#### RDA w/ Sulfur Fxns Data ####

rownames(meta_scaled) %in% rownames(clr.cov.sum.sulf.ko) # check order of DFs
head(meta_scaled)

rda.S.0<-rda(clr.cov.sum.sulf.ko[,!names(clr.cov.sum.sulf.ko) %in% c("SampleID")] ~ DO_Percent_Local+ORP_mV+Temp_DegC+Dissolved_OrganicMatter_RFU+Depth.num+Sulfate_milliM+Sulfide_microM,data=meta_scaled)

# check summary of RDA
rda.S.0
summary(rda.S.0)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.S.0) # 51.74%
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.S.0, permutations = how(nperm=999)) # p = 0.001, significant

## we can also do a permutation test by RDA axis
#anova(rda.S.0, by = "axis", permutations = how(nperm=999)) ### by RDA axis
## or by terms (aka variables)
anova(rda.S.0, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                           Df Variance      F Pr(>F)
#Dissolved_OrganicMatter_RFU  1  0.31700 8.1622  0.024 *

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.S.0)
#DO_Percent_Local                      ORP_mV                   Temp_DegC Dissolved_OrganicMatter_RFU
# 106.181775                 5715.479253                   60.669845                   52.575771
# Depth.num              Sulfate_milliM              Sulfide_microM
# 6.406329                    4.481245                 5347.192460

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/
head(meta_scaled)
## we can use model selection instead of picking variables we think are important (by p values)
rda.S.a = ordistep(rda(clr.cov.sum.sulf.ko[,!names(clr.cov.sum.sulf.ko) %in% c("SampleID")] ~ 1, data = meta_scaled[,c(8,10:11,15:17,20)]),
                     scope=formula(rda.S.0),
                     direction = "forward",
                     permutations = how(nperm=999))
rda.S.a$anova # see significance of individual terms in model

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.S.a2 = ordiR2step(rda(clr.cov.sum.sulf.ko[,!names(clr.cov.sum.sulf.ko) %in% c("SampleID")] ~ 1, data = meta_scaled[,c(8,10:11,15:17,20)]),
                        scope=formula(rda.S.0),
                        permutations = how(nperm=999))
# clr.cov.sum.sulf.ko ~ Temp_DegC + Dissolved_OrganicMatter_RFU + DO_Percent_Local + Depth.num  = best model
rda.S.a2$anova # see significance of individual terms in model
# R2.adj Df     AIC      F Pr(>F)
# + Dissolved_OrganicMatter_RFU 0.20518  1 -3.4237 3.0652  0.048 *
# + DO_Percent_Local            0.41759  1 -5.6095 3.5529  0.042 *

# check best fit model based on above results
anova(rda.S.a, permutations = how(nperm=999))

# Let's double check by removing the variables with high VIF
rda.S1<-rda(clr.cov.sum.sulf.ko[,!names(clr.cov.sum.sulf.ko) %in% c("SampleID")] ~ DO_Percent_Local+Dissolved_OrganicMatter_RFU+ORP_mV+Sulfide_microM+Sulfate_milliM+Depth.num,data=meta_scaled)
summary(rda.S1)
RsquareAdj(rda.S1) # how much variation is explained by our model? 13.27%
anova(rda.S1, by = "terms", permutations = how(nperm=999)) ### by variables
#Dissolved_OrganicMatter_RFU  1  0.37653 6.7363  0.017 *

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.S1)
#DO_Percent_Local Dissolved_OrganicMatter_RFU                      ORP_mV              Sulfide_microM
# 11.340015                   11.969035                  927.593253                  912.903624
# Sulfate_milliM                   Depth.num
# 4.461883                    3.211721
head(meta_scaled)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.S.b1 = ordistep(rda(clr.cov.sum.sulf.ko[,!names(clr.cov.sum.sulf.ko) %in% c("SampleID")] ~ 1, data = meta_scaled[,c(8,10,15:17,20)]),
                      scope=formula(rda.S1),
                      direction = "forward",
                      permutations = how(nperm=999))
rda.S.b1$anova # see significance of individual terms in model

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.S.b2 = ordiR2step(rda(clr.cov.sum.sulf.ko[,!names(clr.cov.sum.sulf.ko) %in% c("SampleID")] ~ 1, data = meta_scaled[,c(8,10,15:17,20)]),
                        scope=formula(rda.S1),
                        permutations = how(nperm=999))
rda.S.b2$anova # see significance of individual terms in model

# check best fit model based on above results
anova(rda.S.b1, permutations = how(nperm=999))

# compare model fits to each other
anova(rda.S.0, rda.S.b1)

rda.S2<-rda(clr.cov.sum.sulf.ko[,!names(clr.cov.sum.sulf.ko) %in% c("SampleID")] ~ Dissolved_OrganicMatter_RFU + DO_Percent_Local,data=meta_scaled)
summary(rda.S2)
RsquareAdj(rda.S2) # how much variation is explained by our model? 8.14%
anova(rda.S2, by = "terms", permutations = how(nperm=999)) ### by variables
# nothing significant

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.S2)
#Dissolved_OrganicMatter_RFU            DO_Percent_Local
#1.852695                    1.852695

head(meta_scaled)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.S.c1 = ordistep(rda(clr.cov.sum.sulf.ko[,!names(clr.cov.sum.sulf.ko) %in% c("SampleID")] ~ 1, data = meta_scaled[,c(8,15)]),
                      scope=formula(rda.S2),
                      direction = "forward",
                      permutations = how(nperm=999))
# clr.cov.sum.sulf.ko ~ Temp_DegC + Dissolved_OrganicMatter_RFU + DO_Percent_Local +      Depth.num + Sulfate_milliM  = best model
rda.S.c1$anova # see significance of individual terms in model
#                               Df    AIC      F Pr(>F)

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.S.c2 = ordiR2step(rda(clr.cov.sum.sulf.ko[,!names(clr.cov.sum.sulf.ko) %in% c("SampleID")] ~ 1, data = meta_scaled[,c(8,15)]),
                        scope=formula(rda.S2),
                        permutations = how(nperm=999))
# clr.cov.sum.sulf.ko ~ Dissolved_OrganicMatter_RFU + Temp_DegC  + DO_%Local = best model
rda.S.c2$anova # see significance of individual terms in model
# R2.adj Df     AIC      F Pr(>F)
# + Dissolved_OrganicMatter_RFU 0.20518  1 -3.4237 3.0652  0.040 *
#   + DO_Percent_Local            0.41759  1 -5.6095 3.5529  0.034 *
#   <All variables>               0.41759

#### Final Sulfur Fxn RDAs ####
# RDA by sampling timepoint
head(meta_scaled)
head(clr.cov.sum.sulf.ko)
rownames(clr.cov.sum.sulf.ko) %in% rownames(meta_scaled) # sanity check 1

# all data
rda.S2$call # best model for all data

rda.S<-rda(clr.cov.sum.sulf.ko[,!names(clr.cov.sum.sulf.ko) %in% c("SampleID")]  ~ Dissolved_OrganicMatter_RFU + DO_Percent_Local,data=meta_scaled)
rda.S
summary(rda.S)
RsquareAdj(rda.S) # how much variation is explained by our model? 18.8 variation
anova(rda.S, permutations = how(nperm=999)) # p-value = 0.017
anova(rda.S, by = "terms", permutations = how(nperm=999))
#                               Df Variance      F Pr(>F)
#Dissolved_OrganicMatter_RFU  1  0.21592 4.1831  0.018 *
# DO_Percent_Local             1  0.18338 3.5529  0.077 .
# Residual                     6  0.30970

#### Plot RDA - ALL Sulfur data ####
#plot(rda.S.aug2021) # depending on how many species you have, this step may take a while
plot(rda.S, scaling = 1)
## scaling = 1 -> emphasizes relationships among sites
plot(rda.S, scaling = 2)
## scaling = 2 -> emphasizes relationships among species

# check summary of RDA
summary(rda.S)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.S) # 41.76%
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
##anova(rda.S, permutations = how(nperm=999)) # p = 0.001, significant

png('figures/MGM_Figs/SSW_AllData_autoplot_rda_example.png',width = 700, height = 600, res=100)
autoplot(rda.S, arrows = TRUE,data = rda.S ,layers=c("biplot","sites"),label = FALSE, label.size = 3, shape = FALSE, loadings = TRUE, loadings.colour = 'blue', loadings.label = TRUE, loadings.label.size = 3, scale= 0)+theme_classic()
dev.off()
## FOR AUTOPLOT -> must load packagve ggvegan first

rda.sum.S<-summary(rda.S)
rda.sum.S$sites[,1:2]
rda.sum.S$cont #cumulative proportion of variance per axis
# RDA1 = 53.11, RDA2 = 3.21

# create data frame w/ RDA axes for sites
# first check rownames of RDA & metadata, then make df
rownames(rda.sum.S$sites) %in% rownames(meta_scaled)
rda.axes.S<-data.frame(RDA1=rda.sum.S$sites[,1], RDA2=rda.sum.S$sites[,2], SampleID=rownames(rda.sum.S$sites), Depth_m=meta_scaled$Depth_m, SampDate=meta_scaled$SampDate)

# create data frame w/ RDA axes for variables
arrows.S<-data.frame(RDA1=rda.sum.S$biplot[,1], RDA2=rda.sum.S$biplot[,2], Label=rownames(rda.sum.S$biplot))
#arrows.S$Label[(arrows.S$Label) == "ORP_mV"] <- "ORP (mV)"
arrows.S$Label[(arrows.S$Label) == "Dissolved_OrganicMatter_RFU"] <- "DOM (RFU)"
arrows.S$Label[(arrows.S$Label) == "DO_Percent_Local"] <- "DO %"
#arrows.S$Label[(arrows.S$Label) == "Temp_DegC"] <- "Temp (C)"

rda.sum.S$cont #cumulative proportion of variance per axis
# RDA1 = 53.11, RDA2 = 3.21

rda.S.plot1<-ggplot(rda.axes.S, aes(x = RDA1, y = RDA2)) + geom_point(size=2) +
  geom_segment(data = arrows.S,mapping = aes(x = 0, y = 0, xend = RDA1, yend = RDA2),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.5,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.S,aes(label = Label, x = RDA1, y = RDA2, fontface="bold"))+
  coord_fixed() + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1))

rda.S.plot2<-ggplot(rda.axes.S, aes(x = RDA1, y = RDA2)) + geom_point(aes(color=as.numeric(as.character(Depth_m)),shape=SampDate),size=4) +
  geom_segment(data = arrows.S,mapping = aes(x = 0, y = 0, xend = RDA1*1.5, yend = RDA2*1.5),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.8,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.S,aes(label = Label, x = RDA1*1.7, y = RDA2*1.7, fontface="bold"), size=4)+
  coord_fixed(ratio = 1, xlim = c(-2,2), ylim = c(-2,2)) + theme_classic() + scale_color_continuous(low="blue3",high="red",trans = 'reverse') +
  scale_shape_discrete(labels=c("August 2021","December 2021","April 2022"),name="Sample Date") +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Sulfur Metabolism in Metagenomes from Salton Seawater",subtitle="Using Centered-Log Ratio Coverage Data",color="Depth (m)") +
  xlab("RDA1 [53.11%]") + ylab("RDA2 [3.21%]")

ggsave(rda.S.plot2,filename = "figures/MGM_Figs/SSW_MGM_S_Fxns_RDA_AllData.png", width=10, height=10, dpi=600)


rda.S.plot3<-ggplot(rda.axes.S, aes(x = RDA1, y = RDA2)) + geom_point(aes(color=as.numeric(as.character(Depth_m)),shape=SampDate),size=5) +
  geom_segment(data = arrows.S,mapping = aes(x = 0, y = 0, xend = RDA1*6, yend = RDA2*6),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 1,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.S,aes(label = Label, x = RDA1*8, y = RDA2*8, fontface="bold"), size=5)+
  coord_fixed(ratio = 1, xlim = c(-10,10), ylim = c(-10,10)) + theme_classic() + scale_color_continuous(low="blue3",high="red",trans = 'reverse') +
  scale_shape_discrete(labels=c("August 2021","December 2021","April 2022"),name="Sample Date") +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea in Salton Seawater",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [53.11%]") + ylab("RDA2 [3.21%]")

ggsave(rda.S.plot3,filename = "figures/MGM_Figs/SSW_MGM_S_Fxns_RDA_AllData_bigger.png", width=15, height=15, dpi=600)

#### Separate ALL Fxns by Timepoints ####
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

# matching data with user defined function -- here is the function, must run to store function in Global env
match_dat<-function(compdata, subset_metadata){
  subset_comp_data = pullrow<-(is.element(row.names(compdata), row.names(subset_metadata)))
  ### * comp data and metadata need to have row names - rownames should be Sample IDs
  subset_comp_data=compdata[pullrow,]
  return(subset_comp_data)
}

# double check that our data frames are ready for this function, aka that they both have the same rownames
## row #s do not have to be the same, but their row names should be in the same format and be able to match up
rownames(mgm.clr)
rownames(August.2021)

# run the function
mgm.clr_AUG21<-match_dat(mgm.clr,August.2021)
mgm.clr_DEC21<-match_dat(mgm.clr,December.2021)
mgm.clr_APR22<-match_dat(mgm.clr,April.2022)

# did the function work the way we wanted it to?

mgm.clr_AUG21[1:3,1:3]
rownames(August.2021) %in% rownames(mgm.clr_AUG21) # hopefully all of the rownames match, aka will get output of TRUE

#### Check ALL Fxn Data Relationship w/ Env Variables (w/ DCA) ####
## remember, CCA assumes that our species have a unimodal relationship with our variables.
### unimodal = one maximum, think upsidedown bellcurve or something
## RDA assumes a linear relationship
## check the assumption

# ALL data
# add pseudocount so row sums are > 0

mgm.clr.pseudo<-mgm.clr[,!names(mgm.clr) %in% c("SampleID")]+2
mgm.clr.dca = decorana(mgm.clr.pseudo)

#plot(mgm.clr.dca) # may take too long to load, do not run unless you have to
mgm.clr.dca #DCA1 axis length = 0.159106; use RDA
## The length of first DCA axis:
## > 4 indicates heterogeneous dataset on which unimodal methods should be used (CCA),
##  < 3 indicates homogeneous dataset for which linear methods are suitable (RDA)
## between 3 and 4 both linear and unimodal methods are OK.

# BY MONTH

mgm.clr_A21.pseudo<-mgm.clr_AUG21[,!names(mgm.clr_AUG21) %in% c("SampleID")]+1
mgm.clr.A21.dca = decorana(mgm.clr_A21.pseudo)
mgm.clr.A21.dca #DCA1 axis length = 0.201019; use RDA

mgm.clr_D21.pseudo<-mgm.clr_DEC21[,!names(mgm.clr_DEC21) %in% c("SampleID")]+2
mgm.clr.D21.dca = decorana(mgm.clr_D21.pseudo)
mgm.clr.D21.dca #DCA1 axis length = 0.093706; use RDA

mgm.clr_A22.pseudo<-mgm.clr_APR22[,!names(mgm.clr_APR22) %in% c("SampleID")]+1
mgm.clr.A22.dca = decorana(mgm.clr_A22.pseudo)
mgm.clr.A22.dca #DCA1 axis length = 0.172273; use RDA

#### RDA w/ ALL Fxns Data ####

rownames(meta_scaled) %in% rownames(mgm.clr) # check order of DFs
head(meta_scaled)

rda.all.0<-rda(mgm.clr[,!names(mgm.clr) %in% c("SampleID")] ~ DO_Percent_Local+ORP_mV+Temp_DegC+Dissolved_OrganicMatter_RFU+Depth.num+Sulfate_milliM+Sulfide_microM,data=meta_scaled)

# check summary of RDA
rda.all.0
summary(rda.all.0)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.all.0) # 24.93%
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.all.0, permutations = how(nperm=999))

## we can also do a permutation test by RDA axis
#anova(rda.all.0, by = "axis", permutations = how(nperm=999)) ### by RDA axis
## or by terms (aka variables)
anova(rda.all.0, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                           Df Variance      F Pr(>F)
#Dissolved_OrganicMatter_RFU  1   49.234 2.3840  0.036 *

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.all.0)
# DO_Percent_Local                      ORP_mV                   Temp_DegC Dissolved_OrganicMatter_RFU
# 106.181775                 5715.479253                   60.669845                   52.575771
# Depth.num              Sulfate_milliM              Sulfide_microM
# 6.406329                    4.481245                 5347.192460

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/
head(meta_scaled)
## we can use model selection instead of picking variables we think are important (by p values)
rda.all.a = ordistep(rda(mgm.clr[,!names(mgm.clr) %in% c("SampleID")] ~ 1, data = meta_scaled[,c(8,10:11,15:17,20)]),
                     scope=formula(rda.all.0),
                     direction = "forward",
                     permutations = how(nperm=999))
rda.all.a$anova # see significance of individual terms in model
# Df    AIC      F Pr(>F)
# + Dissolved_OrganicMatter_RFU  1 48.954 2.2736  0.004 **
#   + Temp_DegC                    1 48.144 2.1990  0.013 *
#   ---

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.all.a2 = ordiR2step(rda(mgm.clr[,!names(mgm.clr) %in% c("SampleID")] ~ 1, data = meta_scaled[,c(8,10:11,15:17,20)]),
                        scope=formula(rda.all.0),
                        permutations = how(nperm=999))
# Step: mgm.clr ~ Temp_DegC  = best model
rda.all.a2$anova # see significance of individual terms in model
# R2.adj Df     AIC      F Pr(>F)
#+ Temp_DegC     0.15932  1 48.722 2.5161  0.017 *

# check best fit model based on above results
anova(rda.all.a, permutations = how(nperm=999))

# Let's double check by removing the variables with high VIF
rda.all1<-rda(mgm.clr[,!names(mgm.clr) %in% c("SampleID")] ~ Dissolved_OrganicMatter_RFU+Temp_DegC,data=meta_scaled)
summary(rda.all1)
RsquareAdj(rda.all1) # how much variation is explained by our model? 13.27%
anova(rda.all1, by = "terms", permutations = how(nperm=999)) ### by variables
# Dissolved_OrganicMatter_RFU  1   53.956 2.6631  0.005 **
#   Temp_DegC                    1   44.553 2.1990  0.034 *
## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.all1)
#Dissolved_OrganicMatter_RFU                   Temp_DegC              Sulfate_milliM
#1.264684                    2.570619                    2.186851

head(meta_scaled)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.all.b1 = ordistep(rda(mgm.clr[,!names(mgm.clr) %in% c("SampleID")] ~ 1, data = meta_scaled[,c(11,15:16)]),
                      scope=formula(rda.all1),
                      direction = "forward",
                      permutations = how(nperm=999))
rda.all.b1$anova # see significance of individual terms in model
# Df    AIC      F Pr(>F)
# + Dissolved_OrganicMatter_RFU  1 48.954 2.2736  0.004 **
#   + Temp_DegC                    1 48.144 2.1990  0.013 *

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.all.b2 = ordiR2step(rda(mgm.clr[,!names(mgm.clr) %in% c("SampleID")] ~ 1, data = meta_scaled[,c(11,15:16)]),
                        scope=formula(rda.all1),
                        permutations = how(nperm=999))
rda.all.b2$anova # see significance of individual terms in model

# check best fit model based on above results
anova(rda.all.b1, permutations = how(nperm=999))

# compare model fits to each other
anova(rda.all.0, rda.all.b1)

#### Final Fxn RDAs ####
# RDA by sampling timepoint
head(meta_scaled)
head(mgm.clr)
rownames(mgm.clr) %in% rownames(meta_scaled) # sanity check 1

# all data
rda.all1$call # best model for all data

rda.all<-rda(mgm.clr[,!names(mgm.clr) %in% c("SampleID")]  ~ Dissolved_OrganicMatter_RFU + Temp_DegC ,data=meta_scaled)
rda.all
summary(rda.all)
RsquareAdj(rda.all) # how much variation is explained by our model? 18.8 variation
anova(rda.all, permutations = how(nperm=999)) # p-value = 0.003
anova(rda.all, by = "terms", permutations = how(nperm=999))
#                               Df Variance      F Pr(>F)
#Dissolved_OrganicMatter_RFU  1   53.956 2.6631  0.004 **
#Temp_DegC                    1   44.553 2.1990  0.036 *
#Residual                     6  121.565

#### Plot RDA - ALL Fxn data ####
#plot(rda.aug2021) # depending on how many species you have, this step may take a while
plot(rda.all, scaling = 1)
## scaling = 1 -> emphasizes relationships among sites
plot(rda.all, scaling = 2)
## scaling = 2 -> emphasizes relationships among species

# check summary of RDA
summary(rda.all)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.all) # 41.76%
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
##anova(rda.all, permutations = how(nperm=999)) # p = 0.001, significant

png('figures/MGM_Figs/SSW_AllData_autoplot_rda_example.png',width = 700, height = 600, res=100)
autoplot(rda.all, arrows = TRUE,data = rda.all ,layers=c("biplot","sites"),label = FALSE, label.size = 3, shape = FALSE, loadings = TRUE, loadings.colour = 'blue', loadings.label = TRUE, loadings.label.size = 3, scale= 0)+theme_classic()
dev.off()
## FOR AUTOPLOT -> must load packagve ggvegan first

rda.sum.all<-summary(rda.all)
rda.sum.all$sites[,1:2]
rda.sum.all$cont #cumulative proportion of variance per axis
# RDA1 = 29.72, RDA2 = 15.04

# create data frame w/ RDA axes for sites
# first check rownames of RDA & metadata, then make df
rownames(rda.sum.all$sites) %in% rownames(meta_scaled)
rda.axes.all<-data.frame(RDA1=rda.sum.all$sites[,1], RDA2=rda.sum.all$sites[,2], SampleID=rownames(rda.sum.all$sites), Depth_m=meta_scaled$Depth_m, SampDate=meta_scaled$SampDate)

# create data frame w/ RDA axes for variables
arrows.all<-data.frame(RDA1=rda.sum.all$biplot[,1], RDA2=rda.sum.all$biplot[,2], Label=rownames(rda.sum.all$biplot))
#arrows.all$Label[(arrows.all$Label) == "ORP_mV"] <- "ORP (mV)"
arrows.all$Label[(arrows.all$Label) == "Dissolved_OrganicMatter_RFU"] <- "DOM (RFU)"
#arrows.all$Label[(arrows.all$Label) == "DO_Percent_Local"] <- "DO %"
arrows.all$Label[(arrows.all$Label) == "Temp_DegC"] <- "Temp (C)"

rda.sum.all$cont #cumulative proportion of variance per axis
# RDA1 = 29.72, RDA2 = 15.04

rda.all.plot1<-ggplot(rda.axes.all, aes(x = RDA1, y = RDA2)) + geom_point(size=2) +
  geom_segment(data = arrows.all,mapping = aes(x = 0, y = 0, xend = RDA1, yend = RDA2),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.5,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.all,aes(label = Label, x = RDA1, y = RDA2, fontface="bold"))+
  coord_fixed() + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1))

rda.all.plot2<-ggplot(rda.axes.all, aes(x = RDA1, y = RDA2)) + geom_point(aes(color=as.numeric(as.character(Depth_m)),shape=SampDate),size=4) +
  geom_segment(data = arrows.all,mapping = aes(x = 0, y = 0, xend = RDA1*4, yend = RDA2*4),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.8,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.all,aes(label = Label, x = RDA1*4.5, y = RDA2*4.5, fontface="bold"), size=4)+
  coord_fixed(ratio = 1, xlim = c(-5,5), ylim = c(-5,5)) + theme_classic() + scale_color_continuous(low="blue3",high="red",trans = 'reverse') +
  scale_shape_discrete(labels=c("August 2021","December 2021","April 2022"),name="Sample Date") +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Metagenome Functions in Salton Seawater",subtitle="Using CLR-Transformed Coverage Data",color="Depth (m)") +
  xlab("RDA1 [29.72%]") + ylab("RDA2 [15.04%]")

ggsave(rda.all.plot2,filename = "figures/MGM_Figs/SSW_MGM_ALLFxns_RDA.png", width=10, height=10, dpi=600)


rda.all.plot3<-ggplot(rda.axes.all, aes(x = RDA1, y = RDA2)) + geom_point(aes(color=as.numeric(as.character(Depth_m)),shape=SampDate),size=5) +
  geom_segment(data = arrows.all,mapping = aes(x = 0, y = 0, xend = RDA1*6, yend = RDA2*6),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 1,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.all,aes(label = Label, x = RDA1*8, y = RDA2*8, fontface="bold"), size=5)+
  coord_fixed(ratio = 1, xlim = c(-10,10), ylim = c(-10,10)) + theme_classic() + scale_color_continuous(low="blue3",high="red",trans = 'reverse') +
  scale_shape_discrete(labels=c("August 2021","December 2021","April 2022"),name="Sample Date") +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea in Salton Seawater",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [30.80%]") + ylab("RDA2 [23.65%]")

ggsave(rda.all.plot3,filename = "figures/MGM_Figs/SSW_MGM_ALLFxns_RDA_bigger.png", width=15, height=15, dpi=600)





#### Separate Carbon Fxns by Timepoints ####
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

# matching data with user defined function -- here is the function, must run to store function in Global env
match_dat<-function(compdata, subset_metadata){
  subset_comp_data = pullrow<-(is.element(row.names(compdata), row.names(subset_metadata)))
  ### * comp data and metadata need to have row names - rownames should be Sample IDs
  subset_comp_data=compdata[pullrow,]
  return(subset_comp_data)
}

# double check that our data frames are ready for this function, aka that they both have the same rownames
## row #s do not have to be the same, but their row names should be in the same format and be able to match up
rownames(clr.cov.sum.carb.ko)
rownames(August.2021)

# run the function
C.KO.clr_AUG21<-match_dat(clr.cov.sum.carb.ko,August.2021)
C.KO.clr_DEC21<-match_dat(clr.cov.sum.carb.ko,December.2021)
C.KO.clr_APR22<-match_dat(clr.cov.sum.carb.ko,April.2022)

# did the function work the way we wanted it to?

C.KO.clr_AUG21[1:3,1:3]
rownames(August.2021) %in% rownames(C.KO.clr_AUG21) # hopefully all of the rownames match, aka will get output of TRUE

#### Check Carbon Fxn Data Relationship w/ Env Variables (w/ DCA) ####
## remember, CCA assumes that our species have a unimodal relationship with our variables.
### unimodal = one maximum, think upsidedown bellcurve or something
## RDA assumes a linear relationship
## check the assumption

# ALL data
# add pseudocount so row sums are > 0

C.clr.pseudo<-clr.cov.sum.carb.ko[,!names(clr.cov.sum.carb.ko) %in% c("SampleID")]+2
C.ko.dca = decorana(C.clr.pseudo)

#plot(C.ko.dca) # may take too long to load, do not run unless you have to
C.ko.dca #DCA1 axis length = 0.154570; use RDA
## The length of first DCA axis:
## > 4 indicates heterogeneous dataset on which unimodal methods should be used (CCA),
##  < 3 indicates homogeneous dataset for which linear methods are suitable (RDA)
## between 3 and 4 both linear and unimodal methods are OK.

# BY MONTH

C.KO.clr_A21.pseudo<-C.KO.clr_AUG21[,!names(C.KO.clr_AUG21) %in% c("SampleID")]+1
C.ko.A21.dca = decorana(C.KO.clr_A21.pseudo)
C.ko.A21.dca #DCA1 axis length = 0.201019; use RDA

C.KO.clr_D21.pseudo<-C.KO.clr_DEC21[,!names(C.KO.clr_DEC21) %in% c("SampleID")]+2
C.ko.D21.dca = decorana(C.KO.clr_D21.pseudo)
C.ko.D21.dca #DCA1 axis length = 0.185162; use RDA

C.KO.clr_A22.pseudo<-C.KO.clr_APR22[,!names(C.KO.clr_APR22) %in% c("SampleID")]+1
C.ko.A22.dca = decorana(C.KO.clr_A22.pseudo)
C.ko.A22.dca #DCA1 axis length = 0.177835; use RDA

#### RDA w/ Carbon Fxns Data ####

rownames(meta_scaled) %in% rownames(clr.cov.sum.carb.ko) # check order of DFs
head(meta_scaled)

rda.C.0<-rda(clr.cov.sum.carb.ko[,!names(clr.cov.sum.carb.ko) %in% c("SampleID")] ~ DO_Percent_Local+ORP_mV+Temp_DegC+Dissolved_OrganicMatter_RFU+Depth.num+Sulfate_milliM+Sulfide_microM,data=meta_scaled)

# check summary of RDA
rda.C.0
summary(rda.C.0)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.C.0) # 16.85%
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.C.0, permutations = how(nperm=999))

## we can also do a permutation test by RDA axis
#anova(rda.C.0, by = "axis", permutations = how(nperm=999)) ### by RDA axis
## or by terms (aka variables)
anova(rda.C.0, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
# nothing significant

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.C.0)
#DO_Percent_Local                      ORP_mV                   Temp_DegC Dissolved_OrganicMatter_RFU
# 106.181775                 5715.479253                   60.669845                   52.575771
# Depth.num              Sulfate_milliM              Sulfide_microM
# 6.406329                    4.481245                 5347.192460

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/
head(meta_scaled)
## we can use model selection instead of picking variables we think are important (by p values)
rda.C.a = ordistep(rda(clr.cov.sum.carb.ko[,!names(clr.cov.sum.carb.ko) %in% c("SampleID")] ~ 1, data = meta_scaled[,c(8,10:11,15:17,20)]),
                   scope=formula(rda.C.0),
                   direction = "forward",
                   permutations = how(nperm=999))
rda.C.a$anova # see significance of individual terms in model
# Df    AIC      F Pr(>F)
# + Temp_DegC                    1 12.245 1.6795  0.063 .

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.C.a2 = ordiR2step(rda(clr.cov.sum.carb.ko[,!names(clr.cov.sum.carb.ko) %in% c("SampleID")] ~ 1, data = meta_scaled[,c(8,10:11,15:17,20)]),
                      scope=formula(rda.C.0),
                      permutations = how(nperm=999))
# clr.cov.sum.carb.ko ~ Temp_DegC + Dissolved_OrganicMatter_RFU + DO_Percent_Local + Depth.num  = best model
rda.C.a2$anova # see significance of individual terms in model
#             Df    AIC      F Pr(>F)
#+ Temp_DegC  1 12.245 1.6795  0.059 .

# check best fit model based on above results
anova(rda.C.a, permutations = how(nperm=999))

# Let's double check by removing the variables with high VIF
rda.C1<-rda(clr.cov.sum.carb.ko[,!names(clr.cov.sum.carb.ko) %in% c("SampleID")] ~ DO_Percent_Local+Dissolved_OrganicMatter_RFU+Sulfate_milliM+Temp_DegC+ORP_mV,data=meta_scaled)
summary(rda.C1)
RsquareAdj(rda.C1) # how much variation is explained by our model? 17.53%
anova(rda.C1, by = "terms", permutations = how(nperm=999)) ### by variables
# nothing

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.C1)
# DO_Percent_Local Dissolved_OrganicMatter_RFU              Sulfate_milliM                   Temp_DegC
# 3.646772                    2.438925                    2.959137                    3.065149
# ORP_mV
# 2.317548
head(meta_scaled)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.C.b1 = ordistep(rda(clr.cov.sum.carb.ko[,!names(clr.cov.sum.carb.ko) %in% c("SampleID")] ~ 1, data = meta_scaled[,c(8,10,11,15:16)]),
                    scope=formula(rda.C1),
                    direction = "forward",
                    permutations = how(nperm=999))
rda.C.b1$anova # see significance of individual terms in model
# Temp near sig

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.C.b2 = ordiR2step(rda(clr.cov.sum.carb.ko[,!names(clr.cov.sum.carb.ko) %in% c("SampleID")] ~ 1, data = meta_scaled[,c(8,10,11,15:16)]),
                      scope=formula(rda.C1),
                      permutations = how(nperm=999))
rda.C.b2$anova # see significance of individual terms in model
# Temp near sig

# check best fit model based on above results
anova(rda.C.b1, permutations = how(nperm=999))

# compare model fits to each other
anova(rda.C.0, rda.C.b1)

rda.C2<-rda(clr.cov.sum.carb.ko[,!names(clr.cov.sum.carb.ko) %in% c("SampleID")] ~ DO_Percent_Local+Dissolved_OrganicMatter_RFU+Sulfate_milliM+Temp_DegC,data=meta_scaled)
summary(rda.C2)
RsquareAdj(rda.C2) # how much variation is explained by our model? 4.72%
anova(rda.C2, by = "terms", permutations = how(nperm=999)) ### by variables
# nothing significant

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.C2)
#DO_Percent_Local Dissolved_OrganicMatter_RFU              Sulfate_milliM                   Temp_DegC
#3.357882                    2.238385                    2.712690                    2.641885

head(meta_scaled)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.C.c1 = ordistep(rda(clr.cov.sum.carb.ko[,!names(clr.cov.sum.carb.ko) %in% c("SampleID")] ~ 1, data = meta_scaled[,c(8,11,15:16)]),
                    scope=formula(rda.C2),
                    direction = "forward",
                    permutations = how(nperm=999))
# Temp near sig
rda.C.c1$anova # see significance of individual terms in model
#                               Df    AIC      F Pr(>F)

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.C.c2 = ordiR2step(rda(clr.cov.sum.carb.ko[,!names(clr.cov.sum.carb.ko) %in% c("SampleID")] ~ 1, data = meta_scaled[,c(8,11,15:16)]),
                      scope=formula(rda.C2),
                      permutations = how(nperm=999))
# clr.cov.sum.carb.ko ~ Dissolved_OrganicMatter_RFU + Temp_DegC  + DO_%Local = best model
rda.C.c2$anova # see significance of individual terms in model
# Temp near sig

rda.C3<-rda(clr.cov.sum.carb.ko[,!names(clr.cov.sum.carb.ko) %in% c("SampleID")] ~ Temp_DegC,data=meta_scaled)
summary(rda.C3)
RsquareAdj(rda.C3) # how much variation is explained by our model? 4.44%
anova(rda.C3, by = "terms", permutations = how(nperm=999)) ### by variables
# nothing significant

head(meta_scaled)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.C.d1 = ordistep(rda(clr.cov.sum.carb.ko[,!names(clr.cov.sum.carb.ko) %in% c("SampleID")] ~ 1, data = meta_scaled[,c(11)]),
                    scope=formula(rda.C3),
                    direction = "forward",
                    permutations = how(nperm=999))
# Temp near sig
rda.C.d1$anova # see significance of individual terms in model
#                               Df    AIC      F Pr(>F)

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.C.d2 = ordiR2step(rda(clr.cov.sum.carb.ko[,!names(clr.cov.sum.carb.ko) %in% c("SampleID")] ~ 1, data = meta_scaled[,c(11)]),
                      scope=formula(rda.C3),
                      permutations = how(nperm=999))
# clr.cov.sum.carb.ko ~ Dissolved_OrganicMatter_RFU + Temp_DegC  + DO_%Local = best model
rda.C.d2$anova # see significance of individual terms in model
# Temp near sig


#### Final Carbon Fxn RDAs ####
# RDA by sampling timepoint
head(meta_scaled)
head(clr.cov.sum.carb.ko)
rownames(clr.cov.sum.carb.ko) %in% rownames(meta_scaled) # sanity check 1

# all data
rda.C3$call # best model for all data, not significant

rda.C<-rda(clr.cov.sum.carb.ko[,!names(clr.cov.sum.carb.ko) %in% c("SampleID")]  ~ Temp_DegC,data=meta_scaled)
rda.C
summary(rda.C)
RsquareAdj(rda.C) # how much variation is explained by our model? 18.8 variation
anova(rda.C, permutations = how(nperm=999)) # p-value = 0.0
anova(rda.C, by = "terms", permutations = how(nperm=999))

#### Plot RDA - ALL Carbon data ####
#plot(rda.C.aug2021) # depending on how many species you have, this step may take a while
plot(rda.C, scaling = 1)
## scaling = 1 -> emphasizes relationships among sites
plot(rda.C, scaling = 2)
## scaling = 2 -> emphasizes relationships among species

# check summary of RDA
summary(rda.C)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.C) # 41.76%
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
##anova(rda.C, permutations = how(nperm=999)) # p = 0.001, significant

png('figures/MGM_Figs/SSW_AllData_autoplot_rda_example.png',width = 700, height = 600, res=100)
autoplot(rda.C, arrows = TRUE,data = rda.C ,layers=c("biplot","sites"),label = FALSE, label.size = 3, shape = FALSE, loadings = TRUE, loadings.colour = 'blue', loadings.label = TRUE, loadings.label.size = 3, scale= 0)+theme_classic()
dev.off()
## FOR AUTOPLOT -> must load packagve ggvegan first

rda.sum.C<-summary(rda.C)
rda.sum.C$sites[,1:2]
rda.sum.C$cont #cumulative proportion of variance per axis
# RDA1 = 53.11, RDA2 = 3.21

# create data frame w/ RDA axes for sites
# first check rownames of RDA & metadata, then make df
rownames(rda.sum.C$sites) %in% rownames(meta_scaled)
rda.axes.C<-data.frame(RDA1=rda.sum.C$sites[,1], RDA2=rda.sum.C$sites[,2], SampleID=rownames(rda.sum.C$sites), Depth_m=meta_scaled$Depth_m, SampDate=meta_scaled$SampDate)

# create data frame w/ RDA axes for variables
arrows.C<-data.frame(RDA1=rda.sum.C$biplot[,1], RDA2=rda.sum.C$biplot[,2], Label=rownames(rda.sum.C$biplot))
#arrows.C$Label[(arrows.C$Label) == "ORP_mV"] <- "ORP (mV)"
arrows.C$Label[(arrows.C$Label) == "Dissolved_OrganicMatter_RFU"] <- "DOM (RFU)"
arrows.C$Label[(arrows.C$Label) == "DO_Percent_Local"] <- "DO %"
#arrows.C$Label[(arrows.C$Label) == "Temp_DegC"] <- "Temp (C)"

rda.sum.C$cont #cumulative proportion of variance per axis
# RDA1 = 53.11, RDA2 = 3.21

rda.C.plot1<-ggplot(rda.axes.C, aes(x = RDA1, y = RDA2)) + geom_point(size=2) +
  geom_segment(data = arrows.C,mapping = aes(x = 0, y = 0, xend = RDA1, yend = RDA2),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.5,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.C,aes(label = Label, x = RDA1, y = RDA2, fontface="bold"))+
  coord_fixed() + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1))

rda.C.plot2<-ggplot(rda.axes.C, aes(x = RDA1, y = RDA2)) + geom_point(aes(color=as.numeric(as.character(Depth_m)),shape=SampDate),size=4) +
  geom_segment(data = arrows.C,mapping = aes(x = 0, y = 0, xend = RDA1*1.5, yend = RDA2*1.5),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.8,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.C,aes(label = Label, x = RDA1*1.7, y = RDA2*1.7, fontface="bold"), size=4)+
  coord_fixed(ratio = 1, xlim = c(-2,2), ylim = c(-2,2)) + theme_classic() + scale_color_continuous(low="blue3",high="red",trans = 'reverse') +
  scale_shape_discrete(labels=c("August 2021","December 2021","April 2022"),name="Sample Date") +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Carbon Fixation in Metagenomes from Salton Seawater",subtitle="Using Centered-Log Ratio Coverage Data",color="Depth (m)") +
  xlab("RDA1 [53.11%]") + ylab("RDA2 [3.21%]")

ggsave(rda.C.plot2,filename = "figures/MGM_Figs/SSW_MGM_C_Fxns_RDA_AllData.png", width=10, height=10, dpi=600)


rda.C.plot3<-ggplot(rda.axes.C, aes(x = RDA1, y = RDA2)) + geom_point(aes(color=as.numeric(as.character(Depth_m)),shape=SampDate),size=5) +
  geom_segment(data = arrows.C,mapping = aes(x = 0, y = 0, xend = RDA1*6, yend = RDA2*6),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 1,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.C,aes(label = Label, x = RDA1*8, y = RDA2*8, fontface="bold"), size=5)+
  coord_fixed(ratio = 1, xlim = c(-10,10), ylim = c(-10,10)) + theme_classic() + scale_color_continuous(low="blue3",high="red",trans = 'reverse') +
  scale_shape_discrete(labels=c("August 2021","December 2021","April 2022"),name="Sample Date") +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea in Salton Seawater",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [53.11%]") + ylab("RDA2 [3.21%]")

ggsave(rda.C.plot3,filename = "figures/MGM_Figs/SSW_MGM_C_Fxns_RDA_AllData_bigger.png", width=15, height=15, dpi=600)


#### Save Progress ####

save.image("data/SSW_MGM_Functions_EnvDriver.Rdata")
