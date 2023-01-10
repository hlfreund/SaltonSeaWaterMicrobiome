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
#save.image("data/Env_Seqs_All/env.seq_analysis.Rdata") # save global env to Rdata file
bac.dat.all[1:4,1:4]
bac.ASV_table[1:4,1:4]
bac.ASV_table[(nrow(bac.ASV_table)-4):(nrow(bac.ASV_table)),(ncol(bac.ASV_table)-4):(ncol(bac.ASV_table))] # last 4 rows & cols

head(metadata)
rownames(metadata)<-metadata$SampleID

#### Create Centered Log-Ratio Table from ASV table ####
# turning bac.ASV_counts_no.contam into phyloseq object called ASV
dim(bac.ASV_table)
ASV<-otu_table(as.matrix(t(bac.ASV_table[,-1])), taxa_are_rows = TRUE)
## need to have ASVs as rows, sample IDs as columns
head(ASV)
class(ASV) # phyloseq otu_table object

# CLR transformation on phyloseq object
asv_clr<-microbiome::transform(ASV, "clr") # transform function from microbiome package
head(asv_clr)

#testCLR <- apply(pseudo.ASV[,-1], 2, function(x) log(x) - log(rowMeans(pseudo.ASV[,-1])))

testCLR2<-decostand(bac.ASV_table[,-1],method = "clr", pseudocount = 1)
testCLR2[1:4,1:4]

# create CLR Sample x Species matrix for input into dist()
# transform so that rownames are SampleIDs, columns are ASV IDs
b.clr<-as.matrix(t(asv_clr))
rownames(b.clr)
b.clr[1:4,1:4]

#### Check Count Data Relationship w/ Env Variables ####
## remember, CCA assumes that our species have a unimodal relationship with our variables.
### unimodal = one maximum, think upsidedown bellcurve or something
## RDA assumes a linear relationship
## check the assumption
chem_data<-subset(metadata, select=c(DO_Percent_Local, ORP_mV, Salinity_ppt))
head(chem_data)
pairs(c(b.clr[1:2,1:2],chem_data))

b.clr[1:4,1:4]

# add pseudocount so row sums are > 0
b.clr.pseudo<-b.clr+1
b.dca = decorana(b.clr.pseudo)

#plot(b.dca) # may take too long to load, do not run unless you have to
summary(b.dca) #DCA1 axis length = 0.79313; use RDA
## The length of first DCA axis:
## > 4 indicates heterogeneous dataset on which unimodal methods should be used (CCA),
##  < 3 indicates homogeneous dataset for which linear methods are suitable (RDA)
## between 3 and 4 both linear and unimodal methods are OK.

#### Prepare Data for RDA/CCA ####
# create metadata df that will contain scaled chemical data
head(metadata)
meta_scaled<-metadata
meta_scaled$SampDate<-factor(meta_scaled$SampDate, levels=c("June.2021","August.2021","December.2021","April.2022"))

head(meta_scaled)
meta_scaled[,8:11]<-as.data.frame(scale(meta_scaled[,8:11], center=TRUE, scale=TRUE)) #not centering before scaling
head(meta_scaled)
rownames(meta_scaled)<-meta_scaled$SampleID
meta_scaled=meta_scaled[rownames(b.clr),]
head(meta_scaled) # another sanity check

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
b.dca #DCA1 axis length = 0.79313; use RDA
## The length of first DCA axis:
## > 4 indicates heterogeneous dataset on which unimodal methods should be used (CCA),
##  < 3 indicates homogeneous dataset for which linear methods are suitable (RDA)
## between 3 and 4 both linear and unimodal methods are OK.

# BY MONTH
b.clr_J21.pseudo<-b.clr_JUN21+1
b.J21.dca = decorana(b.clr_J21.pseudo)
b.J21.dca #DCA1 axis length = 0.38687; use RDA

b.clr_A21.pseudo<-b.clr_AUG21+1
b.A21.dca = decorana(b.clr_A21.pseudo)
b.A21.dca #DCA1 axis length = 0.37763; use RDA

b.clr_D21.pseudo<-b.clr_DEC21+1
b.D21.dca = decorana(b.clr_D21.pseudo)
b.D21.dca #DCA1 axis length = 0.74789; use RDA

b.clr_A22.pseudo<-b.clr_APR22+1
b.A22.dca = decorana(b.clr_A22.pseudo)
b.A22.dca #DCA1 axis length = 0.367263; use RDA

#### Run the RDAs with Env variables of interest ####
b.clr_JUN21[1:4,1:4]
June.2021[1:4,1:4]
# double check that rownames in metadata & composition data match
rownames(August.2021) %in% rownames(b.clr_AUG21) # TRUE = match

# RDA by sampling timepoint
rda.jun2021<-rda(b.clr_JUN21 ~ DO_Percent_Local+ORP_mV+Temp_DegC+Salinity_ppt,data=June.2021)
summary(rda.jun2021)
RsquareAdj(rda.jun2021) # how much variation is explained by our model? 14.28%

rda.aug2021<-rda(b.clr_AUG21 ~ DO_Percent_Local+ORP_mV+Temp_DegC+Salinity_ppt,data=August.2021)
summary(rda.aug2021)
RsquareAdj(rda.aug2021) # how much variation is explained by our model? 5.77%

rda.dec2021<-rda(b.clr_DEC21 ~ DO_Percent_Local+ORP_mV+Temp_DegC+Salinity_ppt,data=December.2021)
summary(rda.dec2021)
RsquareAdj(rda.dec2021) # how much variation is explained by our model? 8.99%

rda.apr2022<-rda(b.clr_APR22 ~ DO_Percent_Local+ORP_mV+Temp_DegC+Salinity_ppt,data=April.2022)
summary(rda.apr2022)
RsquareAdj(rda.apr2022) # how much variation is explained by our model? - 2.57%, meaning model is insufficient to explain variation


rownames(b.clr) %in% rownames(meta_scaled) # TRUE = match
rda.all<-rda(b.clr ~ DO_Percent_Local+ORP_mV+Temp_DegC+Salinity_ppt,data=meta_scaled)
summary(rda.all)
RsquareAdj(rda.all) # how much variation is explained by our model? 22.98%

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
RsquareAdj(rda.all) # 22.98%
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.all, permutations = how(nperm=999)) # p = 0.001, significant

## we can also do a permutation test by RDA axis
anova(rda.all, by = "axis", permutations = how(nperm=999)) ### by RDA axis
## or by terms (aka variables)
anova(rda.all, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant

# check order of data frames
rownames(b.clr) %in% rownames(meta_scaled)

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.all)
## Temp_DegC & Salinity have VIFs > 14 ... very high --> maybe these two variables are correlated?
## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/

## we can use model selection instead of picking variables we think are important (by p values)
rda.all.a = ordistep(rda(b.clr ~ 1, data = meta_scaled[,8:11]),
                         scope=formula(rda.all),
                         direction = "forward",
                         permutations = how(nperm=999))
# b.clr ~ Temp_DegC + ORP_mV + Salinity_ppt + DO_Percent_Local= best model
rda.all.a$anova # see significance of individual terms in model

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.all.a2 = ordiR2step(rda(b.clr ~ 1, data = meta_scaled[,8:11]),
                        scope=formula(rda.all),
                        permutations = how(nperm=999))
# b.clr ~ Temp_DegC + ORP_mV + Salinity_ppt + DO_Percent_Local = best model
rda.all.a2$anova # see significance of individual terms in model

# check best fit model based on above results
anova(rda.all.a, permutations = how(nperm=999)) # p =  0.001, significant

# Let's double check by removing the variables with high VIF
rda.all1<-rda(b.clr ~ DO_Percent_Local+ORP_mV,data=meta_scaled)
summary(rda.all1)
RsquareAdj(rda.all1) # how much variation is explained by our model? 9.56%
anova(rda.all1, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.all1)

## we can use model selection instead of picking variables we think are important -- based on p values
rda.all.b1 = ordistep(rda(b.clr ~ 1, data = meta_scaled[,8:11]),
                     scope=formula(rda.all1),
                     direction = "forward",
                     permutations = how(nperm=999))
# b.clr ~ ORP_mV + DO_Percent_Local= best model
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.all.b2 = ordiR2step(rda(b.clr ~ 1, data = meta_scaled[,8:11]),
                     scope=formula(rda.all1),
                     permutations = how(nperm=999))
# b.clr ~ ORP_mV = best model

# check best fit model based on above results
anova(rda.all.b1, permutations = how(nperm=999)) # p =  0.001, significant

step.res <- ordiR2step(rda.all, rda.all1, perm.max = 999)
step.res$anova

png('autoplot_rda_example.png',width = 700, height = 600, res=100)
autoplot(rda.all.a, arrows = TRUE,data = rda.all ,layers=c("biplot","sites"),label = FALSE, label.size = 3, shape = FALSE, loadings = TRUE, loadings.colour = 'blue', loadings.label = TRUE, loadings.label.size = 3, scale= 0)+theme_classic()
dev.off()
## FOR AUTOPLOT -> must load packagve ggvegan first

rda.sum.all<-summary(rda.all.a)
rda.sum.all$sites[,1:2]
rda.sum.all$cont #cumulative proportion of variance per axis

# create data frame w/ RDA axes for sites
# first check rownames of RDA & metadata, then make df
rownames(rda.sum.all$sites) %in% rownames(meta_scaled)
rda.axes.all<-data.frame(RDA1=rda.sum.all$sites[,1], RDA2=rda.sum.all$sites[,2], SampleID=rownames(rda.sum.all$sites), Depth_m=meta_scaled$Depth_m, SampDate=meta_scaled$SampDate)

# create data frame w/ RDA axes for variables
arrows.all<-data.frame(RDA1=rda.sum.all$biplot[,1], RDA2=rda.sum.all$biplot[,2], Label=rownames(rda.sum.all$biplot))

rda.plot1<-ggplot(rda.axes.all, aes(x = RDA1, y = RDA2)) + geom_point(size=2) +
  geom_segment(data = arrows.all,mapping = aes(x = 0, y = 0, xend = RDA1, yend = RDA2),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.5,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "#7400b8") +
  geom_label(data = arrows.all,aes(label = Label, x = RDA1, y = RDA2, fontface="bold"))+
  coord_fixed() + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1))

rda.plot2<-ggplot(rda.axes.all, aes(x = RDA1, y = RDA2)) + geom_point(aes(color=as.numeric(Depth_m),shape=SampDate),size=3) +
  geom_segment(data = arrows.all,mapping = aes(x = 0, y = 0, xend = RDA1*17, yend = RDA2*17),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.8,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "#7400b8") +
  geom_label(data = arrows.all,aes(label = Label, x = RDA1*19, y = RDA2*19, fontface="bold"), size=4)+
  coord_fixed() + theme_classic() + scale_color_continuous(low="blue3",high="red",trans = 'reverse') +
  scale_shape_discrete(labels=c("June 2021","August 2021","December 2021","April 2022"),name="Sample Date") +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea in Salton Seawater",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [13.26%]") + ylab("RDA2 [7.34%]")

ggsave(rda.plot2,filename = "figures/SSW_16S_RDA_AllData.png", width=15, height=15, dpi=600)

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
RsquareAdj(rda.jun2021) # 14.28%
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.jun2021, permutations = how(nperm=999)) # p = 0.017, significant

## we can also do a permutation test by RDA axis
anova(rda.jun2021, by = "axis", permutations = how(nperm=999)) ### by RDA axis
## or by terms (aka variables)
anova(rda.jun2021, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant

# check order of data frames
rownames(b.clr_JUN21) %in% rownames(June.2021)
## we can use model selection instead of picking variables we think are important
rda.jun2021.a = ordistep(rda(b.clr_JUN21 ~ 1, data = June.2021[,8:11]),
                     scope=formula(rda.jun2021),
                     direction = "forward",
                     permutations = how(nperm=999))
# b.clr_JUN21 ~ Temp_DegC = best model

# create best fit model based on above results
rda.jun2021.b<-rda(b.clr_JUN21 ~ Temp_DegC,data=June.2021)
anova(rda.jun2021.b, permutations = how(nperm=999)) # p =  0.011, significant

png('autoplot_rda_example.png',width = 700, height = 600, res=100)
autoplot(rda.jun2021.b, arrows = TRUE,data = rda.jun2021 ,layers=c("biplot","sites"),label = FALSE, label.size = 3, shape = FALSE, loadings = TRUE, loadings.colour = 'blue', loadings.label = TRUE, loadings.label.size = 3, scale= 0)+theme_classic()
dev.off()
## FOR AUTOPLOT -> must load packagve ggvegan first

rda.sum.j21<-summary(rda.jun2021.b)
rda.sum.j21$sites[,1:2]
rda.sum.j21$cont #cumulative proportion of variance per axis

# create data frame w/ RDA axes for sites
rda.axes.j21<-data.frame(RDA1=rda.sum.j21$sites[,1], PC1=rda.sum.j21$sites[,2], SampleID=rownames(rda.sum.j21$sites), Depth_m=June.2021$Depth_m)

# create data frame w/ RDA axes for variables
arrows.j21<-data.frame(RDA1=rda.sum.j21$biplot[,1], PC1=rda.sum.j21$biplot[,2], Label=rownames(rda.sum.j21$biplot))

rda.plot3<-ggplot(rda.axes.j21, aes(x = RDA1, y = PC1)) + geom_point(size=2) +
  geom_segment(data = arrows.j21,mapping = aes(x = 0, y = 0, xend = RDA1, yend = PC1),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.5,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "#7400b8") +
  geom_label(data = arrows.j21,aes(label = Label, x = RDA1, y = PC1, fontface="bold"))+
  coord_fixed() + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1))

rda.plot4<-ggplot(rda.axes.j21, aes(x = RDA1, y = PC1)) + geom_point(aes(color=as.numeric(Depth_m)),size=4) +
  geom_segment(data = arrows.j21,mapping = aes(x = 0, y = 0, xend = RDA1*12, yend = PC1*12),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.8,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "#7400b8") +
  geom_label(data = arrows.j21,aes(label = Label, x = RDA1*13, y = PC1*13, fontface="bold"), size=4)+
  coord_fixed() + theme_classic() + scale_color_continuous(low="blue3",high="red",trans = 'reverse') +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea in Salton Seawater, June 2021",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [24.92%]") + ylab("PC1 [21.16%]")

ggsave(rda.plot4,filename = "figures/SSW_16S_RDA_Jun2021.png", width=15, height=12, dpi=600)

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
RsquareAdj(rda.aug2021) # 5.77%
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.aug2021, permutations = how(nperm=999)) # p = 0.277, not significant

## we can also do a permutation test by RDA axis
anova(rda.aug2021, by = "axis", permutations = how(nperm=999)) ### by RDA axis
## or by terms (aka variables)
anova(rda.aug2021, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant

# check order of data frames
rownames(b.clr_AUG21) %in% rownames(August.2021)
## we can use model selection instead of picking variables we think are important
rda.aug2021.a = ordistep(rda(b.clr_AUG21 ~ 1, data = August.2021[,8:11]),
                         scope=formula(rda.aug2021),
                         direction = "forward",
                         permutations = how(nperm=999))
# b.clr_AUG21 ~ Salinity_ppt = best model

# create best fit model based on above results
rda.aug2021.b<-rda(b.clr_AUG21 ~ Salinity_ppt,data=August.2021)
anova(rda.aug2021.b, permutations = how(nperm=999)) # p = 0.01, significant

png('autoplot_rda_example.png',width = 700, height = 600, res=100)
autoplot(rda.aug2021.b, arrows = TRUE,data = rda.aug2021 ,layers=c("biplot","sites"),label = FALSE, label.size = 3, shape = FALSE, loadings = TRUE, loadings.colour = 'blue', loadings.label = TRUE, loadings.label.size = 3, scale= 0)+theme_classic()
dev.off()
## FOR AUTOPLOT -> must load packagve ggvegan first

rda.sum.a21<-summary(rda.aug2021.b)
rda.sum.a21$sites[,1:2]
rda.sum.a21$cont #cumulative proportion of variance per axis
# RDA1=20.88%, PC1=19.26%

# create data frame w/ RDA axes for sites
rda.axes.a21<-data.frame(RDA1=rda.sum.a21$sites[,1], PC1=rda.sum.a21$sites[,2], SampleID=rownames(rda.sum.a21$sites), Depth_m=August.2021$Depth_m)

# create data frame w/ RDA axes for variables
arrows.a21<-data.frame(RDA1=rda.sum.a21$biplot[,1], PC1=rda.sum.a21$biplot[,2], Label=rownames(rda.sum.a21$biplot))

rda.plot5<-ggplot(rda.axes.a21, aes(x = RDA1, y = PC1)) + geom_point(size=2) +
  geom_segment(data = arrows.a21,mapping = aes(x = 0, y = 0, xend = RDA1, yend = PC1),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.5,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "#7400b8") +
  geom_label(data = arrows.a21,aes(label = Label, x = RDA1, y = PC1, fontface="bold"))+
  coord_fixed() + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1))

rda.plot6<-ggplot(rda.axes.a21, aes(x = RDA1, y = PC1)) + geom_point(aes(color=as.numeric(Depth_m)),size=4) +
  geom_segment(data = arrows.a21,mapping = aes(x = 0, y = 0, xend = RDA1*8, yend = PC1*8),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.8,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "#7400b8") +
  geom_label(data = arrows.a21,aes(label = Label, x = RDA1*9, y = PC1*9, fontface="bold"), size=6)+
  coord_fixed() + theme_classic() + scale_color_continuous(low="blue3",high="red",trans = 'reverse') +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea in Salton Seawater, August 2021",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [20.88%]") + ylab("PC1 [19.26%]")

ggsave(rda.plot6,filename = "figures/SSW_16S_RDA_aug2021.png", width=16, height=12, dpi=600)

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
RsquareAdj(rda.dec2021) # 8.99%
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.dec2021, permutations = how(nperm=999)) # p = 0.109, not significant

## we can also do a permutation test by RDA axis
anova(rda.dec2021, by = "axis", permutations = how(nperm=999)) ### by RDA axis
## or by terms (aka variables)
anova(rda.dec2021, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant

# check order of data frames
rownames(b.clr_DEC21) %in% rownames(December.2021)
## we can use model selection instead of picking variables we think are important
rda.dec2021.a = ordistep(rda(b.clr_DEC21 ~ 1, data = December.2021[,8:11]),
                         scope=formula(rda.dec2021),
                         direction = "forward",
                         permutations = how(nperm=999))

# create best fit model based on above results
rda.dec2021.b<-rda(b.clr_DEC21 ~ Salinity_ppt,data=December.2021)
anova(rda.dec2021.b, permutations = how(nperm=999)) # p =  0.046,  significant

png('autoplot_rda_example.png',width = 700, height = 600, res=100)
autoplot(rda.dec2021, arrows = TRUE,data = rda.dec2021 ,layers=c("biplot","sites"),label = FALSE, label.size = 3, shape = FALSE, loadings = TRUE, loadings.colour = 'blue', loadings.label = TRUE, loadings.label.size = 3, scale= 0)+theme_classic()
dev.off()
## FOR AUTOPLOT -> must load packagve ggvegan first

rda.sum.d21<-summary(rda.dec2021.b)
rda.sum.d21$sites[,1:2]
rda.sum.d21$cont # cumulative proportion of variation per axis
# RDA1 = 11.18

# create data frame w/ RDA axes for sites
rda.axes.d21<-data.frame(RDA1=rda.sum.d21$sites[,1], PC1=rda.sum.d21$sites[,2], SampleID=rownames(rda.sum.d21$sites), Depth_m=December.2021$Depth_m)

# create data frame w/ RDA axes for variables
arrows.d21<-data.frame(RDA1=rda.sum.d21$biplot[,1], PC1=rda.sum.d21$biplot[,2], Label=rownames(rda.sum.d21$biplot))

rda.plot7<-ggplot(rda.axes.d21, aes(x = RDA1, y = PC1)) + geom_point(size=2) +
  geom_segment(data = arrows.d21,mapping = aes(x = 0, y = 0, xend = RDA1, yend = PC1),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.5,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "#7400b8") +
  geom_label(data = arrows.d21,aes(label = Label, x = RDA1, y = PC1, fontface="bold"))+
  coord_fixed() + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1))

rda.plot8<-ggplot(rda.axes.d21, aes(x = RDA1, y = PC1)) + geom_point(aes(color=as.numeric(Depth_m)),size=3) +
  geom_segment(data = arrows.d21,mapping = aes(x = 0, y = 0, xend = RDA1*17, yend = PC1*17),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.8,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "#7400b8") +
  geom_label(data = arrows.d21,aes(label = Label, x = RDA1*19, y = PC1*19, fontface="bold"), size=4)+
  coord_fixed() + theme_classic() + scale_color_continuous(low="blue3",high="red",trans = 'reverse') +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea in Salton Seawater",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [20.88%]") + ylab("PC1 [34.81%]")

ggsave(rda.plot8,filename = "figures/SSW_16S_RDA_dec2021.png", width=15, height=12, dpi=600)


#### Plot RDA - Apr 2022 ####
#plot(rda.dec2021) # depending on how many species you have, this step may take a while
plot(rda.apr2022, scaling = 1)
## scaling = 1 -> emphasizes relationships among sites
plot(rda.apr2022, scaling = 2)
## scaling = 2 -> emphasizes relationships among species

# check summary of RDA
summary(rda.apr2022)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.apr2022) # 5.77%
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.apr2022, permutations = how(nperm=999)) # p = 0.277, nit significant

## we can also do a permutation test by RDA axis
anova(rda.apr2022, by = "axis", permutations = how(nperm=999)) ### by RDA axis
## or by terms (aka variables)
anova(rda.apr2022, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant

# check order of data frames
rownames(b.clr_APR22) %in% rownames(April.2022)
## we can use model selection instead of picking variables we think are important
rda.apr2022.a = ordistep(rda(b.clr_APR22 ~ 1, data = April.2022[,8:11]),
                         scope=formula(rda.apr2022),
                         direction = "forward",
                         permutations = how(nperm=999))
# no significant variables

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

rda.plot9<-ggplot(rda.axes.a22, aes(x = RDA1, y = RDA2)) + geom_point(size=2) +
  geom_segment(data = arrows.a22,mapping = aes(x = 0, y = 0, xend = RDA1, yend = RDA2),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.5,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "#7400b8") +
  geom_label(data = arrows.a22,aes(label = Label, x = RDA1, y = RDA2, fontface="bold"))+
  coord_fixed() + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1))

rda.plot10<-ggplot(rda.axes.d21, aes(x = RDA1, y = PC1)) + geom_point(aes(color=as.numeric(Depth_m)),size=3) +
  geom_segment(data = arrows.a22,mapping = aes(x = 0, y = 0, xend = RDA1*17, yend = RDA2*17),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.8,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "#7400b8") +
  geom_label(data = arrows.a22,aes(label = Label, x = RDA1*19, y = RDA2*19, fontface="bold"), size=4)+
  coord_fixed() + theme_classic() + scale_color_continuous(low="blue3",high="red",trans = 'reverse') +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea in Salton Seawater",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [7.56%]") + ylab("RDA2 [6.11%]")

ggsave(rda.plot10,filename = "figures/SSW_16S_RDA_apr2022.png", width=15, height=12, dpi=600)

