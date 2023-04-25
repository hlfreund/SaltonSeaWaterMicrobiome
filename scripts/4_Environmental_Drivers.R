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
#load("data/SSW_env.seq_analysis.Rdata") # save global env to Rdata file
bac.dat.all[1:4,1:4]
bac.ASV_table[1:4,1:4]
bac.ASV_table[(nrow(bac.ASV_table)-4):(nrow(bac.ASV_table)),(ncol(bac.ASV_table)-4):(ncol(bac.ASV_table))] # last 4 rows & cols

head(metadata)
rownames(metadata)<-metadata$SampleID

# drop salinity from meta_scaled --> excluding this env variable
meta_scaled<-subset(meta_scaled, select=-c(Salinity_ppt))
head(meta_scaled)

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
b.clr_JUN21<-match_dat(b.clr,June.2021)
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
b.dca #DCA1 axis length = 0.266855; use RDA
## The length of first DCA axis:
## > 4 indicates heterogeneous dataset on which unimodal methods should be used (CCA),
##  < 3 indicates homogeneous dataset for which linear methods are suitable (RDA)
## between 3 and 4 both linear and unimodal methods are OK.

# BY MONTH
b.clr_J21.pseudo<-b.clr_JUN21+1
b.J21.dca = decorana(b.clr_J21.pseudo)
b.J21.dca #DCA1 axis length = 0.263829; use RDA

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

rda.all.0<-rda(b.clr ~ DO_Percent_Local+ORP_mV+Temp_DegC+Dissolved_OrganicMatter_RFU+Depth_m,data=meta_scaled)

# check summary of RDA
rda.all.0
summary(rda.all.0)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.all.0) # 27.03%
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.all.0, permutations = how(nperm=999)) # p = 0.001, significant

## we can also do a permutation test by RDA axis
anova(rda.all.0, by = "axis", permutations = how(nperm=999)) ### by RDA axis
## or by terms (aka variables)
anova(rda.all.0, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.all.0)
# DO_Percent_Local                      ORP_mV                   Temp_DegC Dissolved_OrganicMatter_RFU                    Depth_m2
# 52.040697                    4.454924                    9.498487                    6.248031                   39.501228
# Depth_m3                    Depth_m4                    Depth_m5                    Depth_m7                    Depth_m9
# 3.677240                    4.531623                    7.294175                    5.294976                    5.691828
# Depth_m10                   Depth_m11
# 7.659433                    5.202147

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/
head(meta_scaled)
## we can use model selection instead of picking variables we think are important (by p values)
rda.all.a = ordistep(rda(b.clr ~ 1, data = meta_scaled[,c(6,8,10:11,14)]),
                     scope=formula(rda.all.0),
                     direction = "forward",
                     permutations = how(nperm=999))
# b.clr ~ ORP_mV + Dissolved_OrganicMatter_RFU + DO_Percent_Local +  Depth_m = best model
rda.all.a$anova # see significance of individual terms in model

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.all.a2 = ordiR2step(rda(b.clr ~ 1, data = meta_scaled[,c(6,8,10:11,14)]),
                        scope=formula(rda.all.0),
                        permutations = how(nperm=999))
# b.clr ~ ORP_mV + Dissolved_OrganicMatter_RFU + DO_Percent_Local +  Depth_m  = best model
rda.all.a2$anova # see significance of individual terms in model

# check best fit model based on above results
anova(rda.all.a, permutations = how(nperm=999)) # p =  0.001, significant

# Let's double check by removing the variables with high VIF
rda.all1<-rda(b.clr ~ ORP_mV + Dissolved_OrganicMatter_RFU + Temp_DegC + DO_Percent_Local,data=meta_scaled)
summary(rda.all1)
RsquareAdj(rda.all1) # how much variation is explained by our model? 20.72%
anova(rda.all1, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.all1)
# ORP_mV Dissolved_OrganicMatter_RFU                   Temp_DegC            DO_Percent_Local
# 1.459876                    5.111708                    1.321858                    5.191590
head(meta_scaled)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.all.b1 = ordistep(rda(b.clr ~ 1, data = meta_scaled[,c(8,10:11,14)]),
                      scope=formula(rda.all1),
                      direction = "forward",
                      permutations = how(nperm=999))
# b.clr ~ ORP_mV + Dissolved_OrganicMatter_RFU + DO_Percent_Local +      Temp_DegC  = best model
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.all.b2 = ordiR2step(rda(b.clr ~ 1, data = meta_scaled[,c(8,10:11,14)]),
                        scope=formula(rda.all1),
                        permutations = how(nperm=999))
# b.clr ~ ORP_mV + Dissolved_OrganicMatter_RFU + DO_Percent_Local = best model

# check best fit model based on above results
anova(rda.all.b1, permutations = how(nperm=999)) # p =  0.001, significant

# compare model fits to each other
anova(rda.all.0, rda.all.b1)

rda.all2<-rda(b.clr ~ ORP_mV + Dissolved_OrganicMatter_RFU + DO_Percent_Local,data=meta_scaled)
summary(rda.all2)
RsquareAdj(rda.all2) # how much variation is explained by our model? 11.09%
anova(rda.all2, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.all2)
vif.cca(rda(b.clr ~ ORP_mV + Dissolved_OrganicMatter_RFU + Temp_DegC,data=meta_scaled))

#### RDA - June 2021 ####

rownames(June.2021) %in% rownames(b.clr_JUN21) # check order of DFs
head(June.2021)

rda.jun2021.0<-rda(b.clr_JUN21 ~ DO_Percent_Local+ORP_mV+Temp_DegC+Dissolved_OrganicMatter_RFU+Depth_m,data=June.2021)

# check summary of RDA
rda.jun2021.0
summary(rda.jun2021.0)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.jun2021.0) # 24.88%
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.jun2021.0, permutations = how(nperm=999)) # p = 0.022, significant

## we can also do a permutation test by RDA axis
anova(rda.jun2021.0, by = "axis", permutations = how(nperm=999)) ### by RDA axis
## or by terms (aka variables)
anova(rda.jun2021.0, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.jun2021.0)
# DO_Percent_Local                      ORP_mV                   Temp_DegC Dissolved_OrganicMatter_RFU                    Depth_m5
# 2.900072                    2.900072                          NA                          NA                          NA
# Depth_m10
# NA

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/
head(meta_scaled)
## we can use model selection instead of picking variables we think are important (by p values)
rda.jun2021.a = ordistep(rda(b.clr_JUN21 ~ 1, data = June.2021[,c(6,8,10:11,14)]),
                     scope=formula(rda.jun2021.0),
                     direction = "forward",
                     permutations = how(nperm=999))
# b.clr ~ DO_Percent_Local = best model
rda.jun2021.a$anova # see significance of individual terms in model

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.jun2021.a2 = ordiR2step(rda(b.clr_JUN21 ~ 1, data = June.2021[,c(6,8,10:11,14)]),
                        scope=formula(rda.jun2021.0),
                        permutations = how(nperm=999))
# b.clr ~ Depth_m = best model
rda.jun2021.a2$anova # see significance of individual terms in model

# check best fit model based on above results
anova(rda.jun2021.a, permutations = how(nperm=999)) # p =  0.001, significant
anova(rda.jun2021.a2, permutations = how(nperm=999)) # p =  0.001, significant

# Let's double check by removing the variables with high VIF
rda.jun2021.1<-rda(b.clr_JUN21 ~ DO_Percent_Local+ORP_mV+Temp_DegC+Dissolved_OrganicMatter_RFU,data=June.2021)
summary(rda.jun2021.1)
RsquareAdj(rda.jun2021.1) # how much variation is explained by our model? 24.88%
anova(rda.jun2021.1, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.jun2021.1)
# DO_Percent_Local     ORP_mV
# 2.900072      2.900072

head(June.2021)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.jun2021.b1 = ordistep(rda(b.clr_JUN21 ~ 1, data = June.2021[,c(8,10:11,14)]),
                      scope=formula(rda.jun2021.1),
                      direction = "forward",
                      permutations = how(nperm=999))
# b.clr_JUN21 ~ Temp_DegC = best model
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.jun2021.b2 = ordiR2step(rda(b.clr_JUN21 ~ 1, data = June.2021[,c(8,10:11,14)]),
                        scope=formula(rda.jun2021.1),
                        permutations = how(nperm=999))
# b.clr_JUN21 ~ Temp_DegC = best model

# check best fit model based on above results
anova(rda.jun2021.b1, permutations = how(nperm=999)) # p =  0.001, significant

anova(rda.jun2021.0, rda.jun2021.1) # p =  0.001, significant

rda.jun2021.2<-rda(b.clr_JUN21 ~ DO_Percent_Local+ORP_mV+Temp_DegC,data=June.2021)
summary(rda.jun2021.2)
RsquareAdj(rda.jun2021.2) # how much variation is explained by our model? 24.88%
anova(rda.jun2021.2, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.jun2021.2)
# DO_Percent_Local           ORP_mV        Temp_DegC
# 2.900072         2.900072               NA

head(June.2021)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.jun2021.c1 = ordistep(rda(b.clr_JUN21 ~ 1, data = June.2021[,c(8,10:11)]),
                          scope=formula(rda.jun2021.2),
                          direction = "forward",
                          permutations = how(nperm=999))
# b.clr_JUN21 ~ ORP_mV = best model
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.jun2021.c2 = ordiR2step(rda(b.clr_JUN21 ~ 1, data = June.2021[,c(8,10:11)]),
                            scope=formula(rda.jun2021.2),
                            permutations = how(nperm=999))
# b.clr_JUN21 ~ Temp_DegC = best model

# check best fit model based on above results
anova(rda.jun2021.c1, permutations = how(nperm=999)) # p =  0.001, significant

anova(rda.jun2021.0, rda.jun2021.2) # p =  0.001, significant

rda.jun2021.3<-rda(b.clr_JUN21 ~ ORP_mV+Temp_DegC,data=June.2021)
summary(rda.jun2021.3)
RsquareAdj(rda.jun2021.3) # how much variation is explained by our model? 24.88%
anova(rda.jun2021.3, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.jun2021.3)
# ORP_mV Temp_DegC
# 3.829534  3.829534

head(June.2021)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.jun2021.d1 = ordistep(rda(b.clr_JUN21 ~ 1, data = June.2021[,c(10:11)]),
                          scope=formula(rda.jun2021.3),
                          direction = "forward",
                          permutations = how(nperm=999))
# b.clr_JUN21 ~ Temp_DegC, sometimes ORP_mV = best model
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.jun2021.d2 = ordiR2step(rda(b.clr_JUN21 ~ 1, data = June.2021[,c(10:11)]),
                            scope=formula(rda.jun2021.3),
                            permutations = how(nperm=999))
# b.clr_JUN21 ~ Temp_DegC = best model

# check best fit model based on above results
anova(rda.jun2021.c1, permutations = how(nperm=999)) # p =  0.009, significant

anova(rda.jun2021.0, rda.jun2021.3) # p =  0.001, significant

#### RDA - August 2021 ####

rownames(August.2021) %in% rownames(b.clr_AUG21) # check order of DFs
head(August.2021)

rda.aug2021.0<-rda(b.clr_AUG21 ~ DO_Percent_Local+ORP_mV+Temp_DegC+Dissolved_OrganicMatter_RFU+Depth_m,data=August.2021)

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
anova(rda.aug2021.0, by = "axis", permutations = how(nperm=999)) ### by RDA axis
## or by terms (aka variables)
anova(rda.aug2021.0, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.aug2021.0)

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/
head(August.2021)
## we can use model selection instead of picking variables we think are important (by p values)
rda.aug2021.a = ordistep(rda(b.clr_AUG21 ~ 1, data = August.2021[,c(6,8,10:11,14)]),
                         scope=formula(rda.aug2021.0),
                         direction = "forward",
                         permutations = how(nperm=999))
# b.clr_AUG21 ~ Dissolved_OrganicMatter_RFU = best model
rda.aug2021.a$anova # see significance of individual terms in model

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.aug2021.a2 = ordiR2step(rda(b.clr_AUG21 ~ 1, data = August.2021[,c(6,8,10:11,14)]),
                            scope=formula(rda.aug2021.0),
                            permutations = how(nperm=999))
# too many terms

# check best fit model based on above results
anova(rda.aug2021.a, permutations = how(nperm=999)) # p =  0.001, significant
anova(rda.aug2021.a2, permutations = how(nperm=999)) # p =  0.001, significant

# Let's double check by removing the variables with high VIF
rda.aug2021.1<-rda(b.clr_AUG21 ~ DO_Percent_Local+ORP_mV+Temp_DegC+Dissolved_OrganicMatter_RFU,data=August.2021)
summary(rda.aug2021.1)
RsquareAdj(rda.aug2021.1) # how much variation is explained by our model? 10.36%
anova(rda.aug2021.1, by = "terms", permutations = how(nperm=999)) ### by variables
#                             Df Variance      F Pr(>F)
#DO_Percent_Local             1  115.890 1.3478  0.099 .
#ORP_mV                       1  146.883 1.7082  0.032 *
#Temp_DegC                    1   77.559 0.9020  0.602
#Dissolved_OrganicMatter_RFU  1   73.208 0.8514  0.773
#Residual                     3  257.957

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.aug2021.1)
# DO_Percent_Local        ORP_mV                        Temp_DegC                Dissolved_OrganicMatter_RFU
# 52.92964                    29.11192                    12.28227                    91.79633

head(August.2021)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.aug2021.b1 = ordistep(rda(b.clr_AUG21 ~ 1, data = August.2021[,c(8,10:11,14)]),
                          scope=formula(rda.aug2021.1),
                          direction = "forward",
                          permutations = how(nperm=999))
# b.clr_AUG21 ~ Dissolved_OrganicMatter_RFU = best model
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.aug2021.b2 = ordiR2step(rda(b.clr_AUG21 ~ 1, data = August.2021[,c(8,10:11,14)]),
                            scope=formula(rda.aug2021.1),
                            permutations = how(nperm=999))
# b.clr_AUG21 ~ Dissolved_OrganicMatter_RFU + ORP_mV = best model

# check best fit model based on above results
anova(rda.aug2021.b1, permutations = how(nperm=999)) # p =  0.001, significant

anova(rda.aug2021.0, rda.aug2021.1) # p =  0.001, significant

rda.aug2021.2<-rda(b.clr_AUG21 ~ Dissolved_OrganicMatter_RFU+ORP_mV,data=August.2021)
summary(rda.aug2021.2)
RsquareAdj(rda.aug2021.2) # how much variation is explained by our model? 14.31%
anova(rda.aug2021.2, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.aug2021.2)
# Dissolved_OrganicMatter_RFU  ORP_mV
# 2.834376                    2.834376

head(August.2021)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.aug2021.c1 = ordistep(rda(b.clr_AUG21 ~ 1, data = August.2021[,c(10,14)]),
                          scope=formula(rda.aug2021.2),
                          direction = "forward",
                          permutations = how(nperm=999))
# b.clr_AUG21 ~ Dissolved_OrganicMatter_RFU = best model
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.aug2021.c2 = ordiR2step(rda(b.clr_AUG21 ~ 1, data = August.2021[,c(10,14)]),
                            scope=formula(rda.aug2021.2),
                            permutations = how(nperm=999))
# b.clr_AUG21 ~ ORP_mV = best model

# check best fit model based on above results
anova(rda.aug2021.c1, permutations = how(nperm=999)) # p =  0.001, significant

anova(rda.aug2021.0, rda.aug2021.2) # p =  0.001, significant

#### RDA - December 2021 ####

rownames(December.2021) %in% rownames(b.clr_DEC21) # check order of DFs
head(December.2021)

rda.dec2021.0<-rda(b.clr_DEC21 ~ DO_Percent_Local+ORP_mV+Temp_DegC+Dissolved_OrganicMatter_RFU+Depth_m,data=December.2021)

# check summary of RDA
rda.dec2021.0
summary(rda.dec2021.0)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.dec2021.0) # -0.02242383 -- bad fit
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.dec2021.0, permutations = how(nperm=999)) # not significant

## we can also do a permutation test by RDA axis
anova(rda.dec2021.0, by = "axis", permutations = how(nperm=999)) ### by RDA axis
## or by terms (aka variables)
anova(rda.dec2021.0, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.dec2021.0)

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/
head(December.2021)
## we can use model selection instead of picking variables we think are important (by p values)
rda.dec2021.a = ordistep(rda(b.clr_DEC21 ~ 1, data = December.2021[,c(6,8,10:11,14)]),
                         scope=formula(rda.dec2021.0),
                         direction = "forward",
                         permutations = how(nperm=999))
# b.clr_DEC21 ~ Dissolved_OrganicMatter_RFU + Temp_DegC + DO_Percent_Local = best model but no significant variables
rda.dec2021.a$anova # see significance of individual terms in model

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.dec2021.a2 = ordiR2step(rda(b.clr_DEC21 ~ 1, data = December.2021[,c(6,8,10:11,14)]),
                            scope=formula(rda.dec2021.0),
                            permutations = how(nperm=999))
# DO_Percent_Local, DOM, Temp - not significant

# check best fit model based on above results
anova(rda.dec2021.a, permutations = how(nperm=999)) # p =  0.001, significant
anova(rda.dec2021.a2, permutations = how(nperm=999)) # p =  0.001, significant

# Let's double check by removing the variables with high VIF
rda.dec2021.1<-rda(b.clr_DEC21 ~ Dissolved_OrganicMatter_RFU + Temp_DegC + DO_Percent_Local,data=December.2021)
summary(rda.dec2021.1)
RsquareAdj(rda.dec2021.1) # how much variation is explained by our model? 12.77%
anova(rda.dec2021.1, by = "terms", permutations = how(nperm=999)) ### by variables
#                            Df Variance      F Pr(>F)
#Dissolved_OrganicMatter_RFU  1   112.13 1.8070  0.044 *
#Temp_DegC                    1   178.48 2.8763  0.064 .
#DO_Percent_Local             1    31.77 0.5120  0.901
#Residual                    12   744.62

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.dec2021.1)
# Dissolved_OrganicMatter_RFU  Temp_DegC            DO_Percent_Local
# 33.92260                    56.33309                    45.36584

head(December.2021)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.dec2021.b1 = ordistep(rda(b.clr_DEC21 ~ 1, data = December.2021[,c(8,11,14)]),
                          scope=formula(rda.dec2021.1),
                          direction = "forward",
                          permutations = how(nperm=999))
# b.clr_DEC21 ~ 1 -- nothing significant but all near it
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.dec2021.b2 = ordiR2step(rda(b.clr_DEC21 ~ 1, data = December.2021[,c(8,11,14)]),
                            scope=formula(rda.dec2021.1),
                            permutations = how(nperm=999))
# b.clr_DEC21 ~ DO% = not significant but most variation

# check best fit model based on above results
anova(rda.dec2021.b1, permutations = how(nperm=999)) # p =  0.001, significant

anova(rda.dec2021.0, rda.dec2021.1) # no significant difference

rda.dec2021.2<-rda(b.clr_DEC21 ~ Dissolved_OrganicMatter_RFU+DO_Percent_Local,data=December.2021)
summary(rda.dec2021.2)
RsquareAdj(rda.dec2021.2) # how much variation is explained by our model? 16.34%
anova(rda.dec2021.2, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.dec2021.2)
# Dissolved_OrganicMatter_RFU            DO_Percent_Local
# 1.0888                      1.0888

head(December.2021)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.dec2021.c1 = ordistep(rda(b.clr_DEC21 ~ 1, data = December.2021[,c(8,14)]),
                          scope=formula(rda.dec2021.2),
                          direction = "forward",
                          permutations = how(nperm=999))
# b.clr_DEC21 ~ 1 = nothing significant but close
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.dec2021.c2 = ordiR2step(rda(b.clr_DEC21 ~ 1, data = December.2021[,c(8,14)]),
                            scope=formula(rda.dec2021.2),
                            permutations = how(nperm=999))
# b.clr_DEC21 ~ DO_Percent_Local = best model

# check best fit model based on above results
anova(rda.dec2021.c1, permutations = how(nperm=999)) # p =  0.001, significant

anova(rda.dec2021.0, rda.dec2021.2) # no significant difference

#### RDA - April 2022 ####

rownames(April.2022) %in% rownames(b.clr_APR22) # check order of DFs
head(April.2022)

rda.apr2022.0<-rda(b.clr_APR22 ~ DO_Percent_Local+ORP_mV+Temp_DegC+Dissolved_OrganicMatter_RFU+Depth_m,data=April.2022)

# check summary of RDA
rda.apr2022.0
summary(rda.apr2022.0)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.apr2022.0) # 0.009249572 -- lame fit
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.apr2022.0, permutations = how(nperm=999)) # not significant

## we can also do a permutation test by RDA axis
anova(rda.apr2022.0, by = "axis", permutations = how(nperm=999)) ### by RDA axis
## or by terms (aka variables)
anova(rda.apr2022.0, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.apr2022.0)
#DO_Percent_Local                      ORP_mV                   Temp_DegC Dissolved_OrganicMatter_RFU                    Depth_m3
#8386.86617                    42.61986                  5587.71718                 24533.66496                   737.36967
#Depth_m4                    Depth_m5                    Depth_m7                    Depth_m9                   Depth_m10
#78.16160                          NA                          NA                          NA                     1.75000
#Depth_m11
#NA

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/
head(April.2022)
## we can use model selection instead of picking variables we think are important (by p values)
rda.apr2022.a = ordistep(rda(b.clr_APR22 ~ 1, data = April.2022[,c(6,8,10:11,14)]),
                         scope=formula(rda.apr2022.0),
                         direction = "forward",
                         permutations = how(nperm=999))
# b.clr_APR22 ~ no significant variables
rda.apr2022.a$anova # see significance of individual terms in model

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.apr2022.a2 = ordiR2step(rda(b.clr_APR22 ~ 1, data = April.2022[,c(6,8,10:11,14)]),
                            scope=formula(rda.apr2022.0),
                            permutations = how(nperm=999))
# Depth_m - not significant

# check best fit model based on above results
anova(rda.apr2022.a, permutations = how(nperm=999)) # p =  0.001, significant

# Let's double check by removing the variables with high VIF
rda.apr2022.1<-rda(b.clr_APR22 ~ ORP_mV + Dissolved_OrganicMatter_RFU + Temp_DegC + DO_Percent_Local,data=April.2022)
summary(rda.apr2022.1)
RsquareAdj(rda.apr2022.1) # how much variation is explained by our model? -0.02986758%
anova(rda.apr2022.1, by = "terms", permutations = how(nperm=999)) ### by variables
#  nothing significant

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.apr2022.1)
# ORP_mV Dissolved_OrganicMatter_RFU                   Temp_DegC            DO_Percent_Local
# 10.69370                   528.50030                   347.09750                    23.99276

head(April.2022)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.apr2022.b1 = ordistep(rda(b.clr_APR22 ~ 1, data = April.2022[,c(8,10:11,14)]),
                          scope=formula(rda.apr2022.1),
                          direction = "forward",
                          permutations = how(nperm=999))
# b.clr_APR22 ~ 1 -- nothing significant but all near it
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.apr2022.b2 = ordiR2step(rda(b.clr_APR22 ~ 1, data = April.2022[,c(8,10:11,14)]),
                            scope=formula(rda.apr2022.1),
                            permutations = how(nperm=999))
# b.clr_APR22 - nothing significant

# check best fit model based on above results
anova(rda.apr2022.b1, permutations = how(nperm=999)) # p =  0.001, significant

anova(rda.apr2022.0, rda.apr2022.1) # no significant difference

rda.apr2022.2<-rda(b.clr_APR22 ~ ORP_mV+DO_Percent_Local,data=April.2022)
summary(rda.apr2022.2)
RsquareAdj(rda.apr2022.2) # how much variation is explained by our model? -0.00959%
anova(rda.apr2022.2, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.apr2022.2)
# ORP_mV DO_Percent_Local
# 3.490997         3.490997

head(April.2022)
## we can use model selection instead of picking variables we think are important -- based on p values
rda.apr2022.c1 = ordistep(rda(b.clr_APR22 ~ 1, data = April.2022[,c(8,10)]),
                          scope=formula(rda.apr2022.2),
                          direction = "forward",
                          permutations = how(nperm=999))
# b.clr_APR22 ~ 1 = nothing significant but close
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.apr2022.c2 = ordiR2step(rda(b.clr_APR22 ~ 1, data = April.2022[,c(8,10)]),
                            scope=formula(rda.apr2022.2),
                            permutations = how(nperm=999))

# check best fit model based on above results
anova(rda.apr2022.c1, permutations = how(nperm=999)) # p =  0.001, significant

anova(rda.apr2022.0, rda.apr2022.2) # no significant difference

#### Final RDAs ####
# RDA by sampling timepoint
head(meta_scaled)
head(b.clr)
rownames(b.clr) %in% rownames(meta_scaled) # sanity check 1
rownames(b.clr_JUN21) %in% rownames(June.2021) # sanity check 2

# all data
rda.all<-rda(b.clr ~ ORP_mV + Dissolved_OrganicMatter_RFU + DO_Percent_Local,data=meta_scaled)
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

# June 2021
rda.jun2021<-rda(b.clr_JUN21 ~ ORP_mV+Temp_DegC,data=June.2021)
rda.jun2021
summary(rda.jun2021)
RsquareAdj(rda.jun2021) # how much variation is explained by our model? 24.88%
anova(rda.jun2021, permutations = how(nperm=999)) # p-value = 0.02
anova(rda.jun2021, by = "terms", permutations = how(nperm=999))
#          Df Variance      F Pr(>F)
#ORP_mV     1   229.94 2.4363  0.008 **
#Temp_DegC  1   146.42 1.5514  0.118
#Residual   4   377.53

# August 2021
rda.aug2021<-rda(b.clr_AUG21 ~ Dissolved_OrganicMatter_RFU+ORP_mV,data=August.2021)
summary(rda.aug2021)
RsquareAdj(rda.aug2021) # how much variation is explained by our model? 14.31%
anova(rda.aug2021, permutations = how(nperm=999)) # p-value = 0.009 **
anova(rda.aug2021, by = "terms", permutations = how(nperm=999))
#                           Df Variance      F Pr(>F)
#Dissolved_OrganicMatter_RFU  1   157.93 1.9213  0.002 **
#ORP_mV                       1   102.54 1.2475  0.133
#Residual                     5   411.02

# December 2021
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
rda.apr2022<-rda(b.clr_APR22 ~ ORP_mV+DO_Percent_Local,data=April.2022)
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
RsquareAdj(rda.all) # 18.81%
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
arrows.all$Label[(arrows.all$Label) == "ORP_mV"] <- "ORP (mV)"
arrows.all$Label[(arrows.all$Label) == "Dissolved_OrganicMatter_RFU"] <- "DOM (RFU)"
arrows.all$Label[(arrows.all$Label) == "DO_Percent_Local"] <- "DO%"
#arrows.all$Label[(arrows.all$Label) == "Temp_DegC"] <- "Temp (C)"

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
  geom_label(data = arrows.all,aes(label = Label, x = RDA1*13, y = RDA2*13, fontface="bold"), size=4)+
  coord_fixed(ratio = 1, xlim = c(-10,10), ylim = c(-10,10)) + theme_classic() + scale_color_continuous(low="blue3",high="red",trans = 'reverse') +
  scale_shape_discrete(labels=c("June 2021","August 2021","December 2021","April 2022"),name="Sample Date") +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea in Salton Seawater",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [13.04%]") + ylab("RDA2 [9.88%]")

ggsave(rda.plot2,filename = "figures/EnvDrivers/SSW_16S_RDA_AllData.png", width=15, height=15, dpi=600)


rda.plot3<-ggplot(rda.axes.all, aes(x = RDA1, y = RDA2)) + geom_point(aes(color=as.numeric(as.character(Depth_m)),shape=SampDate),size=5) +
  geom_segment(data = arrows.all,mapping = aes(x = 0, y = 0, xend = RDA1*12, yend = RDA2*12),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 1,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.all,aes(label = Label, x = RDA1*13, y = RDA2*13, fontface="bold"), size=5)+
  coord_fixed(ratio = 1, xlim = c(-10,10), ylim = c(-10,10)) + theme_classic() + scale_color_continuous(low="blue3",high="red",trans = 'reverse') +
  scale_shape_discrete(labels=c("June 2021","August 2021","December 2021","April 2022"),name="Sample Date") +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea in Salton Seawater",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [13.04%]") + ylab("RDA2 [9.88%]")

ggsave(rda.plot3,filename = "figures/EnvDrivers/SSW_16S_RDA_AllData_bigger.png", width=15, height=15, dpi=600)

#### Plot RDA - June 2021 ####
#plot(rda.jun2021) # depending on how many species you have, this step may take a while
plot(rda.jun2021, scaling = 1)
## scaling = 1 -> emphasizes relationships among sites
plot(rda.jun2021, scaling = 2)
## scaling = 2 -> emphasizes relationships among species

# check summary of RDA
summary(rda.jun2021)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.jun2021) # 24.88%
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.jun2021, permutations = how(nperm=999)) # p = 0.022, significant

png('autoplot_rda_example.png',width = 700, height = 600, res=100)
autoplot(rda.jun2021, arrows = TRUE,data = rda.jun2021 ,layers=c("biplot","sites"),label = FALSE, label.size = 3, shape = FALSE, loadings = TRUE, loadings.colour = 'blue', loadings.label = TRUE, loadings.label.size = 3, scale= 0)+theme_classic()
dev.off()
## FOR AUTOPLOT -> must load packagve ggvegan first

rda.sum.j21<-summary(rda.jun2021)
rda.sum.j21$sites[,1:2]
rda.sum.j21$cont #cumulative proportion of variance per axis

# create data frame w/ RDA axes for sites
rda.axes.j21<-data.frame(RDA1=rda.sum.j21$sites[,1], RDA2=rda.sum.j21$sites[,2], SampleID=rownames(rda.sum.j21$sites), Depth_m=June.2021$Depth_m)

# create data frame w/ RDA axes for variables
arrows.j21<-data.frame(RDA1=rda.sum.j21$biplot[,1], RDA2=rda.sum.j21$biplot[,2], Label=rownames(rda.sum.j21$biplot))
arrows.j21$Label[(arrows.j21$Label) == "ORP_mV"] <- "ORP (mV)"
#arrows.j21$Label[(arrows.j21$Label) == "Dissolved_OrganicMatter_RFU"] <- "DOM (RFU)"
#arrows.j21$Label[(arrows.j21$Label) == "DO_Percent_Local"] <- "DO%"
arrows.j21$Label[(arrows.j21$Label) == "Temp_DegC"] <- "Temp (C)"

rda.plot3<-ggplot(rda.axes.j21, aes(x = RDA1, y = RDA2)) + geom_point(size=2) +
  geom_segment(data = arrows.j21,mapping = aes(x = 0, y = 0, xend = RDA1, yend = RDA2),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.5,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.j21,aes(label = Label, x = RDA1, y = RDA2, fontface="bold"))+
  coord_fixed() + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1))

rda.plot4<-ggplot(rda.axes.j21, aes(x = RDA1, y = RDA2)) + geom_point(aes(color=as.numeric(as.character(Depth_m))),size=4) +
  geom_segment(data = arrows.j21,mapping = aes(x = 0, y = 0, xend = RDA1*10, yend = RDA2*10),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.8,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.j21,aes(label = Label, x = RDA1*11.25, y = RDA2*11.25, fontface="bold"), size=4)+
  coord_fixed(ratio = 1, xlim = c(-8,12), ylim = c(-10,10)) + theme_classic() + scale_color_continuous(low="blue3",high="red",trans = 'reverse') +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea in Salton Seawater, June 2021",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [33.56%]") + ylab("RDA2 [16.36%]")

ggsave(rda.plot4,filename = "figures/EnvDrivers/SSW_16S_RDA_Jun2021.png", width=15, height=12, dpi=600)

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
RsquareAdj(rda.aug2021) # 14.31%
## ^^ use this b/c chance correlations can inflate R^2

png('autoplot_rda_example.png',width = 700, height = 600, res=100)
autoplot(rda.aug2021, arrows = TRUE,data = rda.aug2021 ,layers=c("biplot","sites"),label = FALSE, label.size = 3, shape = FALSE, loadings = TRUE, loadings.colour = 'blue', loadings.label = TRUE, loadings.label.size = 3, scale= 0)+theme_classic()
dev.off()
## FOR AUTOPLOT -> must load packagve ggvegan first

rda.sum.a21<-summary(rda.aug2021)
rda.sum.a21$sites[,1:2]
rda.sum.a21$cont #cumulative proportion of variance per axis
# RDA1=27%, RDA2=11.79%

# create data frame w/ RDA axes for sites
rda.axes.a21<-data.frame(RDA1=rda.sum.a21$sites[,1], RDA2=rda.sum.a21$sites[,2], SampleID=rownames(rda.sum.a21$sites), Depth_m=August.2021$Depth_m)

# create data frame w/ RDA axes for variables
arrows.a21<-data.frame(RDA1=rda.sum.a21$biplot[,1], RDA2=rda.sum.a21$biplot[,2], Label=rownames(rda.sum.a21$biplot))
arrows.a21$Label[(arrows.a21$Label) == "ORP_mV"] <- "ORP (mV)"
arrows.a21$Label[(arrows.a21$Label) == "Dissolved_OrganicMatter_RFU"] <- "DOM (RFU)"
#arrows.a21$Label[(arrows.a21$Label) == "DO_Percent_Local"] <- "DO%"
#arrows.a21$Label[(arrows.a21$Label) == "Temp_DegC"] <- "Temp (C)"

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
  xlab("RDA1 [27.00%]") + ylab("RDA2 [11.79%]")

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
RsquareAdj(rda.dec2021) # 16.39%
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

