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

#### Load Global Env to Import Count/ASV Tables ####
load("data/SSeawater_Data_Ready.Rdata") # save global env to Rdata file
load("data/SSW_Amplicon_EnvDriver.Rdata")

bac.dat.all[1:4,1:4]
bac.ASV_table[1:4,1:4]
bac.ASV_table[(nrow(bac.ASV_table)-4):(nrow(bac.ASV_table)),(ncol(bac.ASV_table)-4):(ncol(bac.ASV_table))] # last 4 rows & cols

head(metadata)
head(meta_scaled)

# drop salinity from meta_scaled --> excluding this env variable
meta_scaled<-subset(meta_scaled, select=-c(Salinity_ppt))
head(meta_scaled)

# create column for Depth that is a numeric version of this variable, rather than a factor
meta_scaled$Depth.num<-as.numeric(as.character(meta_scaled$Depth_m))

#### Create Centered Log-Ratio Table from ASV table ####
bac.ASV_table[1:4,1:4]
b.clr<-decostand(bac.ASV_table[,-1],method = "clr", pseudocount = 1) #CLR transformation
b.clr[1:4,1:4]
# df must have rownames are SampleIDs, columns are ASV IDs for vegan functions below\

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

# matching data with user defined function -- here is the function, must run to store function in Global env
match_dat<-function(compdata, subset_metadata){
  subset_comp_data = pullrow<-(is.element(row.names(compdata), row.names(subset_metadata)))
  ### * comp data and metadata need to have row names - rownames should be Sample IDs
  subset_comp_data=compdata[pullrow,]
  return(subset_comp_data)
}

# double check that our data frames are ready for this function, aka that they both have the same rownames
## row #s do not have to be the same, but their row names should be in the same format and be able to match up
rownames(b.clr)
rownames(August.2021)

# run the function
b.clr_AUG21<-match_dat(b.clr,August.2021)
b.clr_DEC21<-match_dat(b.clr,December.2021)
b.clr_APR22<-match_dat(b.clr,April.2022)

# did the function work the way we wanted it to?

b.clr_AUG21[1:4,1:4]
rownames(August.2021) %in% rownames(b.clr_AUG21) # hopefully all of the rownames match, aka will get output of TRUE

#### Check Count Data Relationship w/ Env Variables (w/ DCA) ####
## remember, CCA assumes that our species have a unimodal relationship with our variables.
### unimodal = one maximum, think upsidedown bellcurve or something
## RDA assumes a linear relationship
## check the assumption

# ALL data
# add pseudocount so row sums are > 0
b.clr.pseudo<-b.clr+1
b.dca = decorana(b.clr.pseudo)

#plot(b.dca) # may take too long to load, do not run unless you have to
b.dca #DCA1 axis length = 0.30981; use RDA
## The length of first DCA axis:
## > 4 indicates heterogeneous dataset on which unimodal methods should be used (CCA),
##  < 3 indicates homogeneous dataset for which linear methods are suitable (RDA)
## between 3 and 4 both linear and unimodal methods are OK.

# BY MONTH

b.clr_A21.pseudo<-b.clr_AUG21+1
b.A21.dca = decorana(b.clr_A21.pseudo)
b.A21.dca #DCA1 axis length = 0.211195; use RDA

b.clr_D21.pseudo<-b.clr_DEC21+1
b.D21.dca = decorana(b.clr_D21.pseudo)
b.D21.dca #DCA1 axis length = 0.49828; use RDA

b.clr_A22.pseudo<-b.clr_APR22+1
b.A22.dca = decorana(b.clr_A22.pseudo)
b.A22.dca #DCA1 axis length = 0.215366; use RDA

#### RDA w/ All Data ####

rownames(meta_scaled) %in% rownames(b.clr) # check order of DFs
head(meta_scaled)

rda.all.0<-rda(b.clr ~ DO_Percent_Local+ORP_mV+Temp_DegC+Dissolved_OrganicMatter_RFU+Depth.num+Sulfate_milliM+Sulfide_microM,data=meta_scaled)

# check summary of RDA
rda.all.0
summary(rda.all.0)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.all.0) # 17.61%
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.all.0, permutations = how(nperm=999)) # p = 0.001, significant

## we can also do a permutation test by RDA axis
#anova(rda.all.0, by = "axis", permutations = how(nperm=999)) ### by RDA axis
## or by terms (aka variables)
anova(rda.all.0, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                           Df Variance      F Pr(>F)
# DO_Percent_Local             1    50.47 1.9797  0.060 .
# ORP_mV                       1    44.20 1.7337  0.104
# Temp_DegC                    1    61.63 2.4175  0.031 *
# Dissolved_OrganicMatter_RFU  1   161.55 6.3367  0.001 ***
# Depth.num                    1    28.79 1.1292  0.320
# Sulfate_milliM               1    31.38 1.2310  0.252
# Sulfide_microM               1    12.95 0.5080  0.838
# Residual                    32   815.83

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.all.0)
# DO_Percent_Local              ORP_mV                   Temp_DegC Dissolved_OrganicMatter_RFU
# 10.924631                   29.734615                    6.323630                    2.756809
# Depth.num              Sulfate_milliM              Sulfide_microM
# 2.601581                    2.973335                   33.024737

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/
head(meta_scaled)
## we can use model selection instead of picking variables we think are important (by p values)
rda.all.a = ordistep(rda(b.clr ~ 1, data = meta_scaled[,c(8,10:11,14:16,18)]),
                     scope=formula(rda.all.0),
                     direction = "forward",
                     permutations = how(nperm=999))
# b.clr ~ Dissolved_OrganicMatter_RFU + Sulfate_milliM + Temp_DegC = best model
rda.all.a$anova # see significance of individual terms in model
#                               Df    AIC      F Pr(>F)
#+ Dissolved_OrganicMatter_RFU  1 282.18 4.6688  0.001 ***
#+ Sulfate_milliM               1 279.88 4.2035  0.001 ***
#+ Temp_DegC                    1 279.12 2.5697  0.018 *

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.all.a2 = ordiR2step(rda(b.clr ~ 1, data = meta_scaled[,c(8,10:11,14:16,18)]),
                        scope=formula(rda.all.0),
                        permutations = how(nperm=999))
# b.clr ~ b.clr ~ Dissolved_OrganicMatter_RFU + Sulfate_milliM = best model
rda.all.a2$anova # see significance of individual terms in model
#                               R2.adj Df    AIC      F Pr(>F)
# + Dissolved_OrganicMatter_RFU 0.085983  1 282.18 4.6688  0.001 ***
# + Sulfate_milliM              0.157047  1 279.88 4.2035  0.001 ***
# <All variables>               0.176098

# check best fit model based on above results
anova(rda.all.a, permutations = how(nperm=999)) # p =  0.001, significant

# Let's double check by removing the variables with high VIF
rda.all1<-rda(b.clr ~ Dissolved_OrganicMatter_RFU + Sulfate_milliM + Temp_DegC,data=meta_scaled)
summary(rda.all1)
RsquareAdj(rda.all1) # how much variation is explained by our model? 19.14%
anova(rda.all1, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.all1)
# Dissolved_OrganicMatter_RFU  Sulfate_milliM             Temp_DegC
# 1.394314                     1.607954                    2.033632
head(meta_scaled)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.all.b1 = ordistep(rda(b.clr ~ 1, data = meta_scaled[,c(11,14:15)]),
                      scope=formula(rda.all1),
                      direction = "forward",
                      permutations = how(nperm=999))
# b.clr ~ Dissolved_OrganicMatter_RFU + Sulfate_milliM + Temp_DegC  = best model
rda.all.b1$anova # see significance of individual terms in model
#                               Df    AIC      F Pr(>F)
#+ Dissolved_OrganicMatter_RFU  1 282.18 4.6688  0.001 ***
#+ Sulfate_milliM               1 279.88 4.2035  0.001 ***
#+ Temp_DegC                    1 279.12 2.5697  0.022 *

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.all.b2 = ordiR2step(rda(b.clr ~ 1, data = meta_scaled[,c(11,14:15)]),
                        scope=formula(rda.all1),
                        permutations = how(nperm=999))
# b.clr ~ Dissolved_OrganicMatter_RFU + Sulfate_milliM + Temp_DegC  = best model
rda.all.b2$anova # see significance of individual terms in model
#                               R2.adj Df    AIC      F Pr(>F)
#+ Dissolved_OrganicMatter_RFU 0.085983  1 282.18 4.6688  0.001 ***
#+ Sulfate_milliM              0.157047  1 279.88 4.2035  0.001 ***
#+ Temp_DegC                   0.191353  1 279.12 2.5697  0.027 *

# check best fit model based on above results
anova(rda.all.b1, permutations = how(nperm=999)) # p =  0.001, significant

# compare model fits to each other
anova(rda.all.0, rda.all.b1)

#### RDA - August 2021 ####

rownames(August.2021) %in% rownames(b.clr_AUG21) # check order of DFs
head(August.2021)

rda.aug2021.0<-rda(b.clr_AUG21 ~ DO_Percent_Local+ORP_mV+Temp_DegC+Dissolved_OrganicMatter_RFU+Sulfate_milliM+Sulfide_microM+Depth.num,data=August.2021)

# check summary of RDA
rda.aug2021.0
summary(rda.aug2021.0)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.aug2021.0) # NA -- model fits too well
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.aug2021.0, permutations = how(nperm=999))

## we can also do a permutation test by RDA axis
#anova(rda.aug2021.0, by = "axis", permutations = how(nperm=999)) ### by RDA axis
## or by terms (aka variables)
anova(rda.aug2021.0, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.aug2021.0)
# DO_Percent_Local               ORP_mV                   Temp_DegC Dissolved_OrganicMatter_RFU
# 174.14258                    67.49450                    75.24933                   484.30517
# Sulfate_milliM              Sulfide_microM                   Depth.num
# 13.31191                    25.82085                   340.20799

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/
head(August.2021)
## we can use model selection instead of picking variables we think are important (by p values)
rda.aug2021.a = ordistep(rda(b.clr_AUG21 ~ 1, data = August.2021[,c(8,10:11,14:16,18)]),
                         scope=formula(rda.aug2021.0),
                         direction = "forward",
                         permutations = how(nperm=999))
# b.clr_AUG21 ~ Dissolved_OrganicMatter_RFU = best model
rda.aug2021.a$anova # see significance of individual terms in model

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.aug2021.a2 = ordiR2step(rda(b.clr_AUG21 ~ 1, data = August.2021[,c(8,10:11,14:16,18)]),
                            scope=formula(rda.aug2021.0),
                            permutations = how(nperm=999))
# too many terms

# check best fit model based on above results
anova(rda.aug2021.a, permutations = how(nperm=999)) # p =  0.001, significant
#anova(rda.aug2021.a2, permutations = how(nperm=999)) # p =  0.001, significant

# Let's double check by removing the variables with high VIF, & picking significant variables from ordistep
rda.aug2021.1<-rda(b.clr_AUG21 ~ ORP_mV+Temp_DegC+Dissolved_OrganicMatter_RFU+Sulfate_milliM+Sulfide_microM,data=August.2021)
summary(rda.aug2021.1)
RsquareAdj(rda.aug2021.1) # how much variation is explained by our model? 9.56%
anova(rda.aug2021.1, by = "terms", permutations = how(nperm=999)) ### by variables
#                             Df Variance      F Pr(>F)
# ORP_mV                       1  179.323 2.0690  0.024 *
# Temp_DegC                    1   86.495 0.9980  0.528
# Dissolved_OrganicMatter_RFU  1   57.938 0.6685  0.883
# Sulfate_milliM               1   92.071 1.0623  0.362
# Sulfide_microM               1   81.690 0.9425  0.586
# Residual                     2  173.344

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.aug2021.1)
# ORP_mV                   Temp_DegC Dissolved_OrganicMatter_RFU              Sulfate_milliM
# 32.07755                    11.83552                    39.59135                     1.95635
# Sulfide_microM
# 24.79100

head(August.2021)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.aug2021.b1 = ordistep(rda(b.clr_AUG21 ~ 1, data = August.2021[,c(10:11,14:16)]),
                          scope=formula(rda.aug2021.1),
                          direction = "forward",
                          permutations = how(nperm=999))
# b.clr_AUG21 ~ Dissolved_OrganicMatter_RFU = best model
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.aug2021.b2 = ordiR2step(rda(b.clr_AUG21 ~ 1, data = August.2021[,c(10:11,14:16)]),
                            scope=formula(rda.aug2021.1),
                            permutations = how(nperm=999))
# nothing significant; ORP, Sulfide, & DOM have highest R2

# check best fit model based on above results
anova(rda.aug2021.b1, permutations = how(nperm=999)) # p =  0.001, significant

anova(rda.aug2021.0, rda.aug2021.1) # p =  0.003, significant

# choosing sig variables from ordistep & variables with highest variation
rda.aug2021.2<-rda(b.clr_AUG21 ~ Dissolved_OrganicMatter_RFU+ORP_mV+Sulfide_microM,data=August.2021)
summary(rda.aug2021.2)
RsquareAdj(rda.aug2021.2) # how much variation is explained by our model? 13.26%
anova(rda.aug2021.2, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                             Df Variance      F Pr(>F)
#Dissolved_OrganicMatter_RFU  1   157.93 1.8998  0.005 **
#ORP_mV                       1   102.42 1.2320  0.148
#Sulfide_microM               1    77.98 0.9380  0.531
#Residual                     4   332.53
anova(rda.aug2021.2, by=NULL,permutations = how(nperm=999)) # p =  0.018, significant

vif.cca(rda.aug2021.2)
# Dissolved_OrganicMatter_RFU     ORP_mV              Sulfide_microM
#3.835695                        18.095850                   23.552496

head(August.2021)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.aug2021.c1 = ordistep(rda(b.clr_AUG21 ~ 1, data = August.2021[,c(10,14,16)]),
                          scope=formula(rda.aug2021.2),
                          direction = "forward",
                          permutations = how(nperm=999))
# b.clr_AUG21 ~ Dissolved_OrganicMatter_RFU = best model
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.aug2021.c2 = ordiR2step(rda(b.clr_AUG21 ~ 1, data = August.2021[,c(10,14,16)]),
                            scope=formula(rda.aug2021.2),
                            permutations = how(nperm=999))
# no sig variables

# check best fit model based on above results
anova(rda.aug2021.c1, permutations = how(nperm=999)) # p =  0.001, significant

anova(rda.aug2021.0, rda.aug2021.2) # p =  0.001, significant

# trying a model with DOM, sig variable from above RDA, and variables with lowest VIF
rda.aug2021.3<-rda(b.clr_AUG21 ~ Dissolved_OrganicMatter_RFU+Temp_DegC+Sulfate_milliM,data=August.2021)
summary(rda.aug2021.3)
RsquareAdj(rda.aug2021.3) # how much variation is explained by our model? 13.94%
anova(rda.aug2021.3, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                             Df Variance      F Pr(>F)
#Dissolved_OrganicMatter_RFU  1   157.93 1.9150  0.001 ***
#Temp_DegC                    1    97.43 1.1813  0.188
#Sulfate_milliM               1    85.60 1.0380  0.339
#Residual                     4   329.89
anova(rda.aug2021.3, by=NULL,permutations = how(nperm=999)) # p =  0.005, significant

vif.cca(rda.aug2021.3)
#Dissolved_OrganicMatter_RFU    Temp_DegC              Sulfate_milliM
#2.282208                    1.441769                    1.723704

head(August.2021)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.aug2021.d1 = ordistep(rda(b.clr_AUG21 ~ 1, data = August.2021[,c(11,14:15)]),
                          scope=formula(rda.aug2021.3),
                          direction = "forward",
                          permutations = how(nperm=999))
# b.clr_AUG21 ~ Dissolved_OrganicMatter_RFU = best model
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.aug2021.d2 = ordiR2step(rda(b.clr_AUG21 ~ 1, data = August.2021[,c(11,14:15)]),
                            scope=formula(rda.aug2021.3),
                            permutations = how(nperm=999))
# b.clr_AUG21 ~ Dissolved_OrganicMatter_RFU = best model

anova(rda(b.clr_AUG21 ~ Dissolved_OrganicMatter_RFU,data=August.2021)) # p =  0.001, significant

#### RDA - December 2021 ####

rownames(December.2021) %in% rownames(b.clr_DEC21) # check order of DFs
head(December.2021)

rda.dec2021.0<-rda(b.clr_DEC21 ~ DO_Percent_Local+ORP_mV+Temp_DegC+Dissolved_OrganicMatter_RFU+Sulfate_milliM+Sulfide_microM+Depth.num,data=December.2021)

# check summary of RDA
rda.dec2021.0
summary(rda.dec2021.0)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.dec2021.0) # -0.02232391 -- bad fit
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.dec2021.0, permutations = how(nperm=999)) # not significant

## we can also do a permutation test by RDA axis
#anova(rda.dec2021.0, by = "axis", permutations = how(nperm=999)) ### by RDA axis
## or by terms (aka variables)
anova(rda.dec2021.0, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
# nothing significant

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.dec2021.0)
# DO_Percent_Local              ORP_mV                   Temp_DegC Dissolved_OrganicMatter_RFU
# 686.93160                    38.06763                   202.93243                   212.77436
# Sulfate_milliM              Sulfide_microM                   Depth.num
# 17.68156                    37.94605                   543.52692

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/
head(December.2021)
## we can use model selection instead of picking variables we think are important (by p values)
rda.dec2021.a = ordistep(rda(b.clr_DEC21 ~ 1, data = December.2021[,c(8,10:11,14:16,18)]),
                         scope=formula(rda.dec2021.0),
                         direction = "forward",
                         permutations = how(nperm=999))
# Df    AIC      F Pr(>F)
# + Temp_DegC                    1 112.80 1.5858  0.059 .
# + Dissolved_OrganicMatter_RFU  1 112.74 1.6448  0.068 .
# + DO_Percent_Local             1 111.02 3.4122  0.072 .
# + Depth.num                    1 112.61 1.7702  0.080 .
rda.dec2021.a$anova # see significance of individual terms in model

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.dec2021.a2 = ordiR2step(rda(b.clr_DEC21 ~ 1, data = December.2021[,c(8,10:11,14:16,18)]),
                            scope=formula(rda.dec2021.0),
                            permutations = how(nperm=999))
# DO_Percent_Local, DOM, Temp - not significant

# check best fit model based on above results
anova(rda.dec2021.a, permutations = how(nperm=999)) # not significant
anova(rda.dec2021.a2, permutations = how(nperm=999)) # not significant

# Let's look at near sig variables & ones with low VIF...
rda.dec2021.1<-rda(b.clr_DEC21 ~ DO_Percent_Local+Temp_DegC+Dissolved_OrganicMatter_RFU+Sulfate_milliM+Depth.num,data=December.2021)
summary(rda.dec2021.1)
RsquareAdj(rda.dec2021.1) # how much variation is explained by our model? 3.81%
anova(rda.dec2021.1, by = "terms", permutations = how(nperm=999)) ### by variables
#                            Df Variance      F Pr(>F)
#DO_Percent_Local            1   208.89 2.9693  0.098 .

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.dec2021.1)
# DO_Percent_Local            Temp_DegC      Dissolved_OrganicMatter_RFU              Sulfate_milliM
# 46.367877                   86.958369                   37.179997                    1.900942
# Depth.num
# 20.779679

head(December.2021)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.dec2021.b1 = ordistep(rda(b.clr_DEC21 ~ 1, data = December.2021[,c(8,11,14:15,18)]),
                          scope=formula(rda.dec2021.1),
                          direction = "forward",
                          permutations = how(nperm=999))
#                              Df    AIC      F Pr(>F)
# + Dissolved_OrganicMatter_RFU  1 112.74 1.6448  0.054 .
# + Temp_DegC                    1 112.80 1.5858  0.081 .
# + Depth.num                    1 112.61 1.7702  0.083 .
# + DO_Percent_Local             1 111.02 3.4122  0.094 .
# + Sulfate_milliM               1 114.05 0.4082  0.959

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.dec2021.b2 = ordiR2step(rda(b.clr_DEC21 ~ 1, data = December.2021[,c(8,11,14:15,18)]),
                            scope=formula(rda.dec2021.1),
                            permutations = how(nperm=999))
#   * not significant          R2.adjusted
#+ DO_Percent_Local             0.138536240


# check best fit model based on above results
anova(rda.dec2021.b1, permutations = how(nperm=999)) # p =  0.001, significant

anova(rda.dec2021.0, rda.dec2021.1) # no significant difference

rda.dec2021.2<-rda(b.clr_DEC21 ~ Dissolved_OrganicMatter_RFU+DO_Percent_Local+Temp_DegC+Depth.num,data=December.2021)
summary(rda.dec2021.2)
RsquareAdj(rda.dec2021.2) # how much variation is explained by our model? 8.88%
anova(rda.dec2021.2, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                           Df Variance      F Pr(>F)
#Dissolved_OrganicMatter_RFU  1   112.07 1.7306  0.058 .
#DO_Percent_Local             1   180.94 2.7942  0.085 .

vif.cca(rda.dec2021.2)
#Dissolved_OrganicMatter_RFU   DO_Percent_Local           Temp_DegC                   Depth.num
#34.54139                    45.36586                    68.71426                    11.15461

head(December.2021)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.dec2021.c1 = ordistep(rda(b.clr_DEC21 ~ 1, data = December.2021[,c(8,11,14,18)]),
                          scope=formula(rda.dec2021.2),
                          direction = "forward",
                          permutations = how(nperm=999))
#                              Df    AIC      F Pr(>F)
# + Dissolved_OrganicMatter_RFU  1 112.74 1.6448  0.057 .
# + DO_Percent_Local             1 111.02 3.4122  0.076 .
# + Temp_DegC                    1 112.80 1.5858  0.081 .

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.dec2021.c2 = ordiR2step(rda(b.clr_DEC21 ~ 1, data = December.2021[,c(8,11,14,18)]),
                            scope=formula(rda.dec2021.2),
                            permutations = how(nperm=999))
#           not significant   R2.adjusted
#+ DO_Percent_Local             0.13853624

# check best fit model based on above result
anova(rda.dec2021.0, rda.dec2021.2) # no significant difference

rda.dec2021.3<-rda(b.clr_DEC21 ~ Dissolved_OrganicMatter_RFU+DO_Percent_Local+Temp_DegC,data=December.2021)
summary(rda.dec2021.3)
RsquareAdj(rda.dec2021.3) # how much variation is explained by our model? 16.33%
anova(rda.dec2021.3, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                           Df Variance      F Pr(>F)
#Dissolved_OrganicMatter_RFU  1   112.07 1.8849  0.034 *
#DO_Percent_Local             1   180.94 3.0433  0.050 *
#Residual                    13   772.92

vif.cca(rda.dec2021.3)
#Dissolved_OrganicMatter_RFU    DO_Percent_Local
#1.0888                      1.0888

head(December.2021)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.dec2021.d1 = ordistep(rda(b.clr_DEC21 ~ 1, data = December.2021[,c(8,11,14)]),
                          scope=formula(rda.dec2021.3),
                          direction = "forward",
                          permutations = how(nperm=999))
#                              Df    AIC      F Pr(>F)
#+ Dissolved_OrganicMatter_RFU  1 112.74 1.6448  0.062 .
#+ DO_Percent_Local             1 111.02 3.4122  0.086 .

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.dec2021.d2 = ordiR2step(rda(b.clr_DEC21 ~ 1, data = December.2021[,c(8,11,14)]),
                            scope=formula(rda.dec2021.3),
                            permutations = how(nperm=999))
#           near significant   R2.adjusted
#+ DO_Percent_Local  1 111.02 3.4122  0.076 .

# check best fit model based on above result
anova(rda.dec2021.0, rda.dec2021.3) # no significant difference
anova(rda.dec2021.2, rda.dec2021.3) # no significant difference

rda.dec2021.4<-rda(b.clr_DEC21 ~ Dissolved_OrganicMatter_RFU+DO_Percent_Local,data=December.2021)
summary(rda.dec2021.4)
RsquareAdj(rda.dec2021.4) # how much variation is explained by our model? 16.33%
anova(rda.dec2021.4, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                           Df Variance      F Pr(>F)
#Dissolved_OrganicMatter_RFU  1   112.07 1.8849  0.034 *
#DO_Percent_Local             1   180.94 3.0433  0.050 *
#Residual                    13   772.92

vif.cca(rda.dec2021.4)
#Dissolved_OrganicMatter_RFU    DO_Percent_Local
#1.0888                      1.0888

head(December.2021)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.dec2021.e1 = ordistep(rda(b.clr_DEC21 ~ 1, data = December.2021[,c(8,14)]),
                          scope=formula(rda.dec2021.4),
                          direction = "forward",
                          permutations = how(nperm=999))
#                              Df    AIC      F Pr(>F)
#+ Dissolved_OrganicMatter_RFU  1 112.74 1.6448  0.062 .
#+ DO_Percent_Local             1 111.02 3.4122  0.086 .

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.dec2021.e2 = ordiR2step(rda(b.clr_DEC21 ~ 1, data = December.2021[,c(8,14)]),
                            scope=formula(rda.dec2021.4),
                            permutations = how(nperm=999))
#           near significant   R2.adjusted
#+ DO_Percent_Local  1 111.02 3.4122  0.076 .

#### RDA - April 2022 ####

rownames(April.2022) %in% rownames(b.clr_APR22) # check order of DFs
head(April.2022)

rda.apr2022.0<-rda(b.clr_APR22 ~ DO_Percent_Local+ORP_mV+Temp_DegC+Dissolved_OrganicMatter_RFU+Sulfate_milliM+Sulfide_microM+Depth.num,data=April.2022)

# check summary of RDA
rda.apr2022.0
summary(rda.apr2022.0)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.apr2022.0) # 0.009479891 -- lame fit
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.apr2022.0, permutations = how(nperm=999)) # not significant

## we can also do a permutation test by RDA axis
#anova(rda.apr2022.0, by = "axis", permutations = how(nperm=999)) ### by RDA axis
## or by terms (aka variables)
anova(rda.apr2022.0, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                             Df Variance      F Pr(>F)
# Depth.num                    1   59.922 1.8268  0.039 *

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.apr2022.0)
# DO_Percent_Local            ORP_mV                   Temp_DegC Dissolved_OrganicMatter_RFU
# 74.211673                   65.409468                33173.246389                43723.743178
# Sulfate_milliM              Sulfide_microM           Depth.num
# 3.409798                    7.944354                  523.180072

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/
head(April.2022)
## we can use model selection instead of picking variables we think are important (by p values)
rda.apr2022.a = ordistep(rda(b.clr_APR22 ~ 1, data = April.2022[,c(8,10:11,14:16,18)]),
                         scope=formula(rda.apr2022.0),
                         direction = "forward",
                         permutations = how(nperm=999))
# b.clr_APR22 ~ no significant variables
rda.apr2022.a$anova # see significance of individual terms in model

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.apr2022.a2 = ordiR2step(rda(b.clr_APR22 ~ 1, data = April.2022[,c(8,10:11,14:16,18)]),
                            scope=formula(rda.apr2022.0),
                            permutations = how(nperm=999))
# Depth_m - not significant

# check best fit model based on above results
anova(rda.apr2022.a, permutations = how(nperm=999)) # p =  0.001, significant

# Let's double check by removing variables with high VIF
rda.apr2022.1<-rda(b.clr_APR22 ~ ORP_mV+Sulfate_milliM+Sulfide_microM+Depth.num,data=April.2022)
summary(rda.apr2022.1)
RsquareAdj(rda.apr2022.1) # how much variation is explained by our model? -0.0815
anova(rda.apr2022.1, by = "terms", permutations = how(nperm=999)) ### by variables
#  nothing significant

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.apr2022.1)
# ORP_mV Sulfate_milliM Sulfide_microM      Depth.num
# 1.765386       1.846279       1.403601       2.914192

head(April.2022)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.apr2022.b1 = ordistep(rda(b.clr_APR22 ~ 1, data = April.2022[,c(10,15:16,18)]),
                          scope=formula(rda.apr2022.1),
                          direction = "forward",
                          permutations = how(nperm=999))
# b.clr_APR22 ~ 1 -- nothing significant
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.apr2022.b2 = ordiR2step(rda(b.clr_APR22 ~ 1, data = April.2022[,c(10,15:16,18)]),
                            scope=formula(rda.apr2022.1),
                            permutations = how(nperm=999))
# b.clr_APR22 - nothing significant

# check best fit model based on above results
anova(rda.apr2022.b1, permutations = how(nperm=999)) # p =  0.001, significant

anova(rda.apr2022.0, rda.apr2022.1) # no significant difference

rda.apr2022.2<-rda(b.clr_APR22 ~ ORP_mV+Sulfate_milliM+Sulfide_microM,data=April.2022)
summary(rda.apr2022.2)
RsquareAdj(rda.apr2022.2) # how much variation is explained by our model? -0.04356
anova(rda.apr2022.2, by = "terms", permutations = how(nperm=999)) ### by variables
# ^ nothing significant

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.apr2022.2)
#ORP_mV Sulfate_milliM Sulfide_microM
#1.214250       1.392552       1.166129

head(April.2022)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.apr2022.c1 = ordistep(rda(b.clr_APR22 ~ 1, data = April.2022[,c(10,15:16)]),
                          scope=formula(rda.apr2022.2),
                          direction = "forward",
                          permutations = how(nperm=999))
# b.clr_APR22 ~ 1 = nothing significant but close
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.apr2022.c2 = ordiR2step(rda(b.clr_APR22 ~ 1, data = April.2022[,c(10,15:16)]),
                            scope=formula(rda.apr2022.2),
                            permutations = how(nperm=999))

# check best fit model based on above results
anova(rda.apr2022.c1, permutations = how(nperm=999)) # not significant

anova(rda.apr2022.0, rda.apr2022.2) # no significant difference

# trying variables that had lowest, though still bad, p values
rda.apr2022.3<-rda(b.clr_APR22 ~ ORP_mV+DO_Percent_Local+Dissolved_OrganicMatter_RFU+Temp_DegC+Depth.num,data=April.2022)
summary(rda.apr2022.3)
RsquareAdj(rda.apr2022.3) # how much variation is explained by our model? 0.004298883
anova(rda.apr2022.3, by = "terms", permutations = how(nperm=999)) ### by variables
# ^ nothing significant
anova(rda.apr2022.3, by = NULL, permutations = how(nperm=999)) ### not sig overall

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.apr2022.3)
# ORP_mV            DO_Percent_Local Dissolved_OrganicMatter_RFU                   Temp_DegC
# 36.79253                    24.15112                 13733.50274                 10411.94887
# Depth.num
# 174.50822

head(April.2022)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.apr2022.e1 = ordistep(rda(b.clr_APR22 ~ 1, data = April.2022[,c(8,10,11,14,18)]),
                          scope=formula(rda.apr2022.3),
                          direction = "forward",
                          permutations = how(nperm=999))
# b.clr_APR22 ~ 1 = nothing significant
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.apr2022.e2 = ordiR2step(rda(b.clr_APR22 ~ 1, data = April.2022[,c(8,10,11,14,18)]),
                            scope=formula(rda.apr2022.3),
                            permutations = how(nperm=999))
# nothing significant

rda.apr2022.4<-rda(b.clr_APR22 ~ ORP_mV+DO_Percent_Local,data=April.2022)
summary(rda.apr2022.4)
RsquareAdj(rda.apr2022.4) # how much variation is explained by our model? -0.00949
anova(rda.apr2022.4, by = "terms", permutations = how(nperm=999)) ### by variables
# ^ nothing significant
anova(rda.apr2022.4, by = NULL, permutations = how(nperm=999)) ### not sig overall

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.apr2022.4)
#ORP_mV DO_Percent_Local
#3.490997         3.490997

head(April.2022)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.apr2022.f1 = ordistep(rda(b.clr_APR22 ~ 1, data = April.2022[,c(8,10)]),
                          scope=formula(rda.apr2022.4),
                          direction = "forward",
                          permutations = how(nperm=999))
# b.clr_APR22 ~ 1 = nothing significant
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.apr2022.f2 = ordiR2step(rda(b.clr_APR22 ~ 1, data = April.2022[,c(8,10)]),
                            scope=formula(rda.apr2022.4),
                            permutations = how(nperm=999))
# nothing significant

rda.apr2022.5<-rda(b.clr_APR22 ~ Depth.num,data=April.2022)
summary(rda.apr2022.5)
RsquareAdj(rda.apr2022.5) # how much variation is explained by our model? 0.004826971
anova(rda.apr2022.5, by = "terms", permutations = how(nperm=999)) ### by variables
# ^ nothing significant
anova(rda.apr2022.5, by = NULL, permutations = how(nperm=999)) ### not sig overall

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.apr2022.5)

# no sig variables for april 2022

#### Final RDAs ####
# RDA by sampling timepoint
head(meta_scaled)
head(b.clr)
rownames(b.clr) %in% rownames(meta_scaled) # sanity check 1

# all data
rda.all1$call # best model for all data

rda.all<-rda(b.clr ~ Dissolved_OrganicMatter_RFU+Sulfate_milliM+Temp_DegC,data=meta_scaled)
rda.all
summary(rda.all)
RsquareAdj(rda.all) # how much variation is explained by our model? 18.8 variation
anova(rda.all, permutations = how(nperm=999)) # p-value = 0.001
anova(rda.all, by = "terms", permutations = how(nperm=999))
#                               Df Variance      F Pr(>F)
#ORP_mV                       1    85.86 3.7186  0.002 **
#Dissolved_OrganicMatter_RFU  1    88.11 3.8160  0.001 ***
#DO_Percent_Local             1   141.30 6.1199  0.001 ***
#Residual                    43   992.83

# August 2021
rda.aug2021.3$call # best model

rda.aug2021<-rda(b.clr_AUG21 ~ Dissolved_OrganicMatter_RFU+Temp_DegC+Sulfate_milliM,data=August.2021)
summary(rda.aug2021)
RsquareAdj(rda.aug2021) # how much variation is explained by our model? 14.31%
anova(rda.aug2021, permutations = how(nperm=999)) # p-value = 0.009 **
anova(rda.aug2021, by = "terms", permutations = how(nperm=999))
#                           Df Variance      F Pr(>F)
#Dissolved_OrganicMatter_RFU  1   157.93 1.9213  0.002 **
#ORP_mV                       1   102.54 1.2475  0.133
#Residual                     5   411.02

# December 2021
rda.dec2021.4$call # best model from above

rda.dec2021<-rda(b.clr_DEC21 ~ Dissolved_OrganicMatter_RFU+DO_Percent_Local,data=December.2021)
summary(rda.dec2021)
RsquareAdj(rda.dec2021) # how much variation is explained by our model? 16.34%
anova(rda.dec2021, permutations = how(nperm=999)) # p-value = 0.026
anova(rda.dec2021, by = "terms", permutations = how(nperm=999))
#                             Df Variance      F Pr(>F)
#Dissolved_OrganicMatter_RFU  1   112.13 1.8841  0.025 *
#DO_Percent_Local             1   181.23 3.0453  0.062 .
#Residual                    13   773.64

# April 2022 [no significant env variables]
rda.apr2022.5$call  #best model though no significant variables -- had the higest R^2

rda.apr2022<-rda(b.clr_APR22 ~ Depth.num,data=April.2022)
summary(rda.apr2022)
RsquareAdj(rda.apr2022) # how much variation is explained by our model? -0.959%
anova(rda.apr2022, permutations = how(nperm=999)) # p-value = 0.383
anova(rda.apr2022, by = "terms", permutations = how(nperm=999))

#### Plot RDA - ALL data ####
#plot(rda.jun2021) # depending on how many species you have, this step may take a while
plot(rda.all, scaling = 1)
## scaling = 1 -> emphasizes relationships among sites
plot(rda.all, scaling = 2)
## scaling = 2 -> emphasizes relationships among species

# check summary of RDA
summary(rda.all)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.all) # 19.14%
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.all, permutations = how(nperm=999)) # p = 0.001, significant

png('autoplot_rda_example.png',width = 700, height = 600, res=100)
autoplot(rda.all, arrows = TRUE,data = rda.all ,layers=c("biplot","sites"),label = FALSE, label.size = 3, shape = FALSE, loadings = TRUE, loadings.colour = 'blue', loadings.label = TRUE, loadings.label.size = 3, scale= 0)+theme_classic()
dev.off()
## FOR AUTOPLOT -> must load packagve ggvegan first

rda.sum.all<-summary(rda.all)
rda.sum.all$sites[,1:2]
rda.sum.all$cont #cumulative proportion of variance per axis

# create data frame w/ RDA axes for sites
# first check rownames of RDA & metadata, then make df
rownames(rda.sum.all$sites) %in% rownames(meta_scaled)
rda.axes.all<-data.frame(RDA1=rda.sum.all$sites[,1], RDA2=rda.sum.all$sites[,2], SampleID=rownames(rda.sum.all$sites), Depth_m=meta_scaled$Depth_m, SampDate=meta_scaled$SampDate)

# create data frame w/ RDA axes for variables
arrows.all<-data.frame(RDA1=rda.sum.all$biplot[,1], RDA2=rda.sum.all$biplot[,2], Label=rownames(rda.sum.all$biplot))
#arrows.all$Label[(arrows.all$Label) == "ORP_mV"] <- "ORP (mV)"
arrows.all$Label[(arrows.all$Label) == "Dissolved_OrganicMatter_RFU"] <- "DOM (RFU)"
arrows.all$Label[(arrows.all$Label) == "Sulfate_milliM"] <- "Sulfate (milliM)"
arrows.all$Label[(arrows.all$Label) == "Temp_DegC"] <- "Temp (C)"

rda.sum.all$cont #cumulative proportion of variance per axis
# RDA1=15.2%, RDA2=5.55%

rda.plot1<-ggplot(rda.axes.all, aes(x = RDA1, y = RDA2)) + geom_point(size=2) +
  geom_segment(data = arrows.all,mapping = aes(x = 0, y = 0, xend = RDA1, yend = RDA2),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.5,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.all,aes(label = Label, x = RDA1, y = RDA2, fontface="bold"))+
  coord_fixed() + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1))

rda.plot2<-ggplot(rda.axes.all, aes(x = RDA1, y = RDA2)) + geom_point(aes(color=as.numeric(as.character(Depth_m)),shape=SampDate),size=4) +
  geom_segment(data = arrows.all,mapping = aes(x = 0, y = 0, xend = RDA1*12, yend = RDA2*12),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.8,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.all,aes(label = Label, x = RDA1*13.3, y = RDA2*13.3, fontface="bold"), size=4)+
  coord_fixed(ratio = 1, xlim = c(-12,12), ylim = c(-12,12)) + theme_classic() + scale_color_continuous(low="blue3",high="red",trans = 'reverse') +
  scale_shape_discrete(labels=c("August 2021","December 2021","April 2022"),name="Sample Date") +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea in Salton Seawater",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [15.20%]") + ylab("RDA2 [5.55%]")

ggsave(rda.plot2,filename = "figures/EnvDrivers/SSW_16S_RDA_AllData.png", width=15, height=15, dpi=600)


rda.plot3<-ggplot(rda.axes.all, aes(x = RDA1, y = RDA2)) + geom_point(aes(color=as.numeric(as.character(Depth_m)),shape=SampDate),size=5) +
  geom_segment(data = arrows.all,mapping = aes(x = 0, y = 0, xend = RDA1*12, yend = RDA2*12),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 1,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.all,aes(label = Label, x = RDA1*13.3, y = RDA2*13.3, fontface="bold"), size=5)+
  coord_fixed(ratio = 1, xlim = c(-12,12), ylim = c(-12,12)) + theme_classic() + scale_color_continuous(low="blue3",high="red",trans = 'reverse') +
  scale_shape_discrete(labels=c("August 2021","December 2021","April 2022"),name="Sample Date") +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea in Salton Seawater",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [15.20%]") + ylab("RDA2 [5.55%]")

ggsave(rda.plot3,filename = "figures/EnvDrivers/SSW_16S_RDA_AllData_bigger.png", width=15, height=15, dpi=600)

#### Plot RDA - Aug 2021 ####
#plot(rda.aug2021) # depending on how many species you have, this step may take a while
plot(rda.aug2021, scaling = 1)
## scaling = 1 -> emphasizes relationships among sites
plot(rda.aug2021, scaling = 2)
## scaling = 2 -> emphasizes relationships among species

# check summary of RDA
summary(rda.aug2021)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.aug2021) # 13.94%
## ^^ use this b/c chance correlations can inflate R^2

png('autoplot_rda_example.png',width = 700, height = 600, res=100)
autoplot(rda.aug2021, arrows = TRUE,data = rda.aug2021 ,layers=c("biplot","sites"),label = FALSE, label.size = 3, shape = FALSE, loadings = TRUE, loadings.colour = 'blue', loadings.label = TRUE, loadings.label.size = 3, scale= 0)+theme_classic()
dev.off()
## FOR AUTOPLOT -> must load packagve ggvegan first

rda.sum.a21<-summary(rda.aug2021)
rda.sum.a21$sites[,1:2]
rda.sum.a21$cont #cumulative proportion of variance per axis
# RDA1=26.31%, RDA2=12.53%

# create data frame w/ RDA axes for sites
rda.axes.a21<-data.frame(RDA1=rda.sum.a21$sites[,1], RDA2=rda.sum.a21$sites[,2], SampleID=rownames(rda.sum.a21$sites), Depth_m=August.2021$Depth_m)

# create data frame w/ RDA axes for variables
arrows.a21<-data.frame(RDA1=rda.sum.a21$biplot[,1], RDA2=rda.sum.a21$biplot[,2], Label=rownames(rda.sum.a21$biplot))
#arrows.a21$Label[(arrows.a21$Label) == "ORP_mV"] <- "ORP (mV)"
arrows.a21$Label[(arrows.a21$Label) == "Dissolved_OrganicMatter_RFU"] <- "DOM (RFU)"
arrows.a21$Label[(arrows.a21$Label) == "Sulfate_milliM"] <- "Sulfate (milliM)"
arrows.a21$Label[(arrows.a21$Label) == "Temp_DegC"] <- "Temp (C)"

rda.plot5<-ggplot(rda.axes.a21, aes(x = RDA1, y = RDA2)) + geom_point(size=2) +
  geom_segment(data = arrows.a21,mapping = aes(x = 0, y = 0, xend = RDA1, yend = RDA2),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.5,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.a21,aes(label = Label, x = RDA1, y = RDA2, fontface="bold"))+
  coord_fixed() + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1))

rda.plot6<-ggplot(rda.axes.a21, aes(x = RDA1, y = RDA2)) + geom_point(aes(color=as.numeric(as.character(Depth_m))),size=4) +
  geom_segment(data = arrows.a21,mapping = aes(x = 0, y = 0, xend = RDA1*8, yend = RDA2*8),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.8,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.a21,aes(label = Label, x = RDA1*9, y = RDA2*9, fontface="bold"), size=4)+
  coord_fixed(ratio = 1, xlim = c(-10,10), ylim = c(-10,10)) + theme_classic() + scale_color_continuous(low="blue3",high="red",trans = 'reverse') +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea in Salton Seawater, August 2021",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [26.31%]") + ylab("RDA2 [12.53%]")

ggsave(rda.plot6,filename = "figures/EnvDrivers/SSW_16S_RDA_Aug2021.png", width=16, height=12, dpi=600)

#### Plot RDA - Dec 2021 ####
#plot(rda.dec2021) # depending on how many species you have, this step may take a while
plot(rda.dec2021, scaling = 1)
## scaling = 1 -> emphasizes relationships among sites
plot(rda.dec2021, scaling = 2)
## scaling = 2 -> emphasizes relationships among species

# check summary of RDA
summary(rda.dec2021)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.dec2021) # 16.33%
## ^^ use this b/c chance correlations can inflate R^2

png('autoplot_rda_example.png',width = 700, height = 600, res=100)
autoplot(rda.dec2021, arrows = TRUE,data = rda.dec2021 ,layers=c("biplot","sites"),label = FALSE, label.size = 3, shape = FALSE, loadings = TRUE, loadings.colour = 'blue', loadings.label = TRUE, loadings.label.size = 3, scale= 0)+theme_classic()
dev.off()
## FOR AUTOPLOT -> must load packagve ggvegan first

rda.sum.d21<-summary(rda.dec2021)
rda.sum.d21$sites[,1:2]
rda.sum.d21$cont # cumulative proportion of variation per axis
# RDA1 = 20.35, RDA2 = 7.14

# create data frame w/ RDA axes for sites
rda.axes.d21<-data.frame(RDA1=rda.sum.d21$sites[,1], RDA2=rda.sum.d21$sites[,2], SampleID=rownames(rda.sum.d21$sites), Depth_m=December.2021$Depth_m)

# create data frame w/ RDA axes for variables
arrows.d21<-data.frame(RDA1=rda.sum.d21$biplot[,1], RDA2=rda.sum.d21$biplot[,2], Label=rownames(rda.sum.d21$biplot))
#arrows.d21$Label[(arrows.d21$Label) == "ORP_mV"] <- "ORP (mV)"
arrows.d21$Label[(arrows.d21$Label) == "Dissolved_OrganicMatter_RFU"] <- "DOM (RFU)"
arrows.d21$Label[(arrows.d21$Label) == "DO_Percent_Local"] <- "DO%"
#arrows.d21$Label[(arrows.d21$Label) == "Temp_DegC"] <- "Temp (C)"

rda.plot7<-ggplot(rda.axes.d21, aes(x = RDA1, y = RDA2)) + geom_point(size=2) +
  geom_segment(data = arrows.d21,mapping = aes(x = 0, y = 0, xend = RDA1, yend = RDA2),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.5,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.d21,aes(label = Label, x = RDA1, y = RDA2, fontface="bold"))+
  coord_fixed() + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1))

rda.plot8<-ggplot(rda.axes.d21, aes(x = RDA1, y = RDA2)) + geom_point(aes(color=as.numeric(as.character(Depth_m))),size=3) +
  geom_segment(data = arrows.d21,mapping = aes(x = 0, y = 0, xend = RDA1*9, yend = RDA2*9),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.8,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.d21,aes(label = Label, x = RDA1*10, y = RDA2*10, fontface="bold"), size=4)+
  coord_fixed(ratio = 1, xlim = c(-10,10), ylim = c(-10,10)) + theme_classic() + scale_color_continuous(low="blue3",high="red",trans = 'reverse') +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea in Salton Seawater",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [20.35%]") + ylab("RDA2 [7.14%]")

ggsave(rda.plot8,filename = "figures/EnvDrivers/SSW_16S_RDA_Dec2021.png", width=15, height=12, dpi=600)


#### Plot RDA - Apr 2022 ####
#plot(rda.dec2021) # depending on how many species you have, this step may take a while
plot(rda.apr2022, scaling = 1)
## scaling = 1 -> emphasizes relationships among sites
plot(rda.apr2022, scaling = 2)
## scaling = 2 -> emphasizes relationships among species

# check summary of RDA
summary(rda.apr2022)

png('autoplot_rda_example.png',width = 700, height = 600, res=100)
autoplot(rda.apr2022, arrows = TRUE,data = rda.apr2022 ,layers=c("biplot","sites"),label = FALSE, label.size = 3, shape = FALSE, loadings = TRUE, loadings.colour = 'blue', loadings.label = TRUE, loadings.label.size = 3, scale= 0)+theme_classic()
dev.off()
## FOR AUTOPLOT -> must load packagve ggvegan first

rda.sum.a22<-summary(rda.apr2022)
rda.sum.a22$sites[,1:2]
rda.sum.a22$cont

# create data frame w/ RDA axes for sites
rda.axes.a22<-data.frame(RDA1=rda.sum.a22$sites[,1], RDA2=rda.sum.a22$sites[,2], SampleID=rownames(rda.sum.a22$sites), Depth_m=April.2022$Depth_m)

# create data frame w/ RDA axes for variables
arrows.a22<-data.frame(RDA1=rda.sum.a22$biplot[,1], RDA2=rda.sum.a22$biplot[,2], Label=rownames(rda.sum.a22$biplot))
arrows.a22$Label[(arrows.a22$Label) == "ORP_mV"] <- "ORP (mV)"
arrows.a22$Label[(arrows.a22$Label) == "Dissolved_OrganicMatter_RFU"] <- "DOM (RFU)"
arrows.a22$Label[(arrows.a22$Label) == "DO_Percent_Local"] <- "DO%"
#arrows.a22$Label[(arrows.a22$Label) == "Temp_DegC"] <- "Temp (C)"

rda.plot9<-ggplot(rda.axes.a22, aes(x = RDA1, y = RDA2)) + geom_point(size=2) +
  geom_segment(data = arrows.a22,mapping = aes(x = 0, y = 0, xend = RDA1, yend = RDA2),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.5,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.a22,aes(label = Label, x = RDA1, y = RDA2, fontface="bold"))+
  coord_fixed() + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1))

rda.plot10<-ggplot(rda.axes.d21, aes(x = RDA1, y = PC1)) + geom_point(aes(color=as.numeric(as.character(Depth_m))),size=3) +
  geom_segment(data = arrows.a22,mapping = aes(x = 0, y = 0, xend = RDA1*17, yend = RDA2*17),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.8,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.a22,aes(label = Label, x = RDA1*19, y = RDA2*19, fontface="bold"), size=4)+
  coord_fixed(ratio = 1, xlim = c(-10,10), ylim = c(-10,10)) + theme_classic() + scale_color_continuous(low="blue3",high="red",trans = 'reverse') +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea in Salton Seawater",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [7.56%]") + ylab("RDA2 [6.11%]")

ggsave(rda.plot10,filename = "figures/EnvDrivers/SSW_16S_RDA_apr2022.png", width=15, height=12, dpi=600)

