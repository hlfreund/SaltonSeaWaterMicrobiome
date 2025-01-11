#### Set WD & Load Libraries ####
getwd() # use setwd("path/to/files") if you are not in the right directory
#setwd("/Volumes/HLF_SSD/Aronson_Lab_Data/Salton_Sea/SaltonSeaWater")
suppressPackageStartupMessages({ # load packages quietly
  library(phyloseq)
  library(ggplot2)
  library(vegan)
  library(ggpubr)
  #library(scales)
  library(grid)
  library(data.table)
  library(ape)
  #library(apeglm)
  library(plyr)
  library(dplyr)
  library(viridis)
  library(readxl)
  library(metagenomeSeq)
  #library(heatmaply)
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
  library(fitdistrplus)
  library(logspline)
  library(shades)
  library(ALDEx2)
  library(rstatix)
  library(devtools)
  library(decontam)
  library(pairwiseAdonis)
})

#### Load Data & See Info About Data ####
load("data/Metagenomes/Analysis/SSW_mgm_analysis.Rdata") # load Rdata to global env
#load("data/Metagenomes/Analysis/SSW_MGM_FxnBetaDiv.Rdata")

head(mgm_meta)
arsen.fxns[1:4,]
ko.cov.sum_table[1:4,1:4] # contains the sum of coverages per gene per KO -- featureCounts was normalized by gene length across samples first to get coverage, then summed up per KO ID
head(sulf.path.clr.ars)

# ABOUT THE DATA:
# Before transformations (i.e., VST, CLR, etc) were done, the following was performed:
# featureCounts counted reads that mapped to genes in contigs
# Reads mapped to genes were divided by gene length for all genes across all samples
# Gene coverage was then added together for each KO ID, since multiple genes were assigned the same KO ID
# Summed coverage per KO was then transformed via median-ratio, vst, and clr

## For pathway analyses -- after gene coverage was calculated and added together per KO ID, they were added together for each pathway
## summed coverages per KO ID, then per pathway were transformed by CLR

# NOTE about CLR transformation:
## uses a pseudocount of 1 to replace 0s, which is why not all 0s are treated equally
## need to look into robustCLR, which uses CLR transformation without 0s. Need more info on this methodology...

#### Sum Gene Coverage by KO ID in Contigs Before Transformations ####
mgm_fxns.cov[1:4,]
mgm_fxns.cov_noNA<-as.data.frame(mgm_fxns.cov[!is.na(mgm_fxns.cov$KO_ID),]) # drop genes with KOs given NA as assignment (aka no KO ID assigned at all)

ko.cov.sum_table<-as.data.frame(dcast(mgm_fxns.cov_noNA, SampleID~KO_ID, value.var="CovPerGene", fun.aggregate=sum)) ###
ko.cov.sum_table[1:4,1:4]
rownames(ko.cov.sum_table)<-ko.cov.sum_table$SampleID
ko.cov.sum_table[1:4,1:4]

# check rownames of summed transformed feature coverage data & metadata
rownames(ko.cov.sum_table) %in% rownames(mgm_meta)

### Pull out traits of interest ####
# create unique list of KO ID and functions
# check for duplicates to make sure each KO_ID has a unique function assignment
head(mgm_fxns.cov)
NA %in% mgm_fxns.cov$CovPerGene # just to ensure there are no NAs for genes in this df

ko_fxns1<-unique(mgm_fxns.cov[,1:3]) # first subset out data based on unique KO functions

n_occur <- data.frame(table(ko_fxns1$KO_ID)) # see how many duplicates there are of KO IDs, compare duplicates
n_occur[n_occur$Freq > 1,] # what traits appear more than once?

ko_ID<-unique(data.frame(KO_ID=ko_fxns1$KO_ID)) # get a list of unique KO IDs in data

ko_fxns<-as.data.frame(ko_fxns1[!is.na(ko_fxns1$KO_ID),])  # use unique KO ID list to subset out KO function data
head(ko_fxns)

## pull out functions of interest
#sulfur.fxns<-as.data.frame(ko_fxns[grep("sulf|thio", ko_fxns$KO_Function), ]) # pull out sulfur functions
sulfur.fxns<-ko_fxns[which(ko_fxns$KO_ID %in% sulf.kegg$KO_ID),]
nitro.fxns<-ko_fxns[which(ko_fxns$KO_ID %in% nitro.kegg$KO_ID),]
carb.fxns<-ko_fxns[which(ko_fxns$KO_ID %in% carb.kegg$KO_ID),]
All_GOI.fxns<-ko_fxns[which(ko_fxns$KO_ID %in% all_goi.kegg$KO_ID),]
osmo.fxns<-ko_fxns[which(ko_fxns$KO_ID %in% osmo.kegg$KO_ID),]
selen.fxns<-ko_fxns[which(ko_fxns$KO_ID %in% selen.kegg$KO_ID),]
arsen.fxns<-ko_fxns[which(ko_fxns$KO_ID %in% arsen.kegg$KO_ID),]
HS.fxns<-ko_fxns[which(ko_fxns$KO_ID %in% heatshock.kegg$KO_ID),]
metal.fxns<-ko_fxns[which(ko_fxns$KO_ID %in% metal.re.kegg$KO_ID),]
photo.fxn<-ko_fxns[which(ko_fxns$KO_ID %in% photo.kegg$KO_ID),]
aero.fxn<-ko_fxns[which(ko_fxns$KO_ID %in% aero.kegg$KO_ID),]

#### Pull Out Sulfur Metabolic Fxns from CLR data - with NAs ####
## heatmaps of traits of interest

ko.cov.sum_table[1:4,1:4]

# pull out sulfur functions from CLR transformed, summed coverages (summed gene coverage per KO)
sulf.ko.covsums<-ko.cov.sum_table[,which(colnames(ko.cov.sum_table) %in% sulfur.fxns$KO_ID)] # merge CLR data w/ S fxns found in contigs from KOFamScan
sulf.ko.covsums$SampleID<-rownames(sulf.ko.covsums)
sulf.ko.covsums.melt<-melt(sulf.ko.covsums, by="SampleID")
colnames(sulf.ko.covsums.melt)[which(names(sulf.ko.covsums.melt) == "variable")] <- "KO_ID"
colnames(sulf.ko.covsums.melt)[which(names(sulf.ko.covsums.melt) == "value")] <- "SumCovPerKO"
head(sulf.ko.covsums.melt) #sanity check

sulf.kofxn.covsums<-merge(sulf.ko.covsums.melt,sulf.kegg,by.x=c("KO_ID"),by.y=c("KO_ID")) # merge data w/ KO assignments from KEGG db
head(sulf.kofxn.covsums)
colnames(sulf.kofxn.covsums)[which(names(sulf.kofxn.covsums) == "KO_Function")] <- "KO_Function.KEGG" # so we know they are KO assignments from KEGG db website
sulf.path.ko.covsums<-as.data.frame(dcast(sulf.kofxn.covsums, SampleID~Pathway, value.var="SumCovPerKO", fun.aggregate=sum)) ###just dcast, nothing is being added here!
rownames(sulf.path.ko.covsums)<-sulf.path.ko.covsums$SampleID
sulf.path.ko.covsums[1:4,]

# sanity check
sulf.path.ko.covsums$`Assimilatory Sulfate Reduction`
dim(sulf.path.ko.covsums)

#### Centered Log Ratio Transformation - by Pathway ####
# for Reviewer 1 comments; wants statistical test to compare pathway differences across metagenomes
## summed up relative gene coverage by KO, then by pathway, then CLR transforming the pathways
sulf.path.ko.covsums[1:4,]

# df must have rownames are SampleIDs, columns are ASV IDs for vegan functions below\
sulf.path.clr<-decostand(sulf.path.ko.covsums[,-1],method = "clr",pseudocount=1) #CLR transformation
sulf.path.clr[1:4,]
# NOTE: CLR transformation does not treat all 0s equally, it has to do with the pseudocount that's added before transformation
# The method can operate only with positive data; a common way to deal with zeroes is to add pseudocount, either by adding it manually to the input data, or by using the argument pseudocount as in decostand(x, method = "clr", pseudocount = 1).
# Adding pseudocount will inevitably introduce some bias; see the rclr method for one available solution

sulf.path.clr$SampleID<-rownames(sulf.path.clr)

#### Compare Variance by All S Pathways ####
sulf.path.clr[1:4,1:4] # sample IDs are rows, genes are columns

# check rownames of CLR & VST transformed feature count data & metadata
rownames(sulf.path.clr) %in% rownames(mgm_meta)

## PCOA with CLR transformed data first
# calculate our Euclidean distance matrix using CLR data
s.path.euc.clr_dist <- dist(sulf.path.clr, method = "euclidean")

# creating our hierarcical clustering dendrogram
s.path.euc.clr_clust <- hclust(s.path.euc.clr_dist, method="ward.D2")

# let's make it a little nicer...
s.path.euc.clr_dend <- as.dendrogram(s.path.euc.clr_clust, hang=0.2)
s.path.dend_cols <- as.character(mgm_meta$SampDate_Color[order.dendrogram(s.path.euc.clr_dend)])
labels_colors(s.path.euc.clr_dend) <- s.path.dend_cols

par(mar=c(1,1,1,1))
plot(s.path.euc.clr_dend, ylab="CLR Euclidean Distance",cex = 0.5) + title(main = "Bacteria/Archaea Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 3)
legend("topright",legend = c("August 2021","December 2021","April 2022"),cex=.8,col = c("#ef781c","#03045e","#059c3f"),pch = 15, bty = "n")
# Control is dark blue ("#218380"), #Alternaria is light blue ("#73d2de")
dev.off()

# let's use our Euclidean distance matrix from before
s.path.pcoa.clr <- pcoa(s.path.euc.clr_dist) # pcoa of euclidean distance matrix = PCA of euclidean distance matrix
##save.image("data/ssw_clr.euc.dist1_3.7.23.Rdata")

# The proportion of variances explained is in its element values$Relative_eig
s.path.pcoa.clr$values

# extract principal coordinates
s.path.pcoa.clr.vectors<-data.frame(s.path.pcoa.clr$vectors)
s.path.pcoa.clr.vectors$SampleID<-rownames(s.path.pcoa.clr$vectors)

# merge pcoa coordinates w/ metadata
s.path.pcoa.clr.meta<-merge(s.path.pcoa.clr.vectors, mgm_meta, by.x="SampleID", by.y="SampleID")
s.path.pcoa.clr.meta$SampleMonth
s.path.pcoa.clr.meta$SampDate

head(s.path.pcoa.clr.meta)

s.path.pcoa.clr$values # pull out Relative (Relative_eig) variation % to add to axes labels

# create PCoA ggplot fig
ggplot(s.path.pcoa.clr.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(SampDate)), size=4)+theme_bw()+
  labs(title="PCoA: Sulfur Cycling in Salton Seawater",subtitle="Using CLR Transformed, Summed Gene Coverage per KO, per Pathway",
       color="Sample Date")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Sample Date",values=unique(s.path.pcoa.clr.meta$SampDate_Color[order(s.path.pcoa.clr.meta$SampDate)]),labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022")) +
  xlab("PC1 [90.02%]") + ylab("PC2 [7.28%]")

ggsave(pcoa5,filename = "figures/Revised/MGM_Figs/FxnDiv/PCoAs/CenterLogRatioTransformation/SSW_MGM_pcoa_CLR_SummedCoverage_Per_KO_sampdate.png", width=12, height=10, dpi=600,create.dir=TRUE)

# sample month shape, depth color
ggplot(s.path.pcoa.clr.meta, aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(color=as.numeric(Depth_m),shape=SampleMonth), size=5)+theme_bw()+
  labs(title="PCoA: Sulfur Cycling in Salton Seawater",subtitle="Using CLR Transformed, Summed Gene Coverage per KO, per Pathway",xlab="PC1", ylab="PC2",color="Depth (m)")+
  theme_classic()+ theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),axis.text = element_text(size=12),axis.text.x = element_text(vjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17))+
  scale_color_continuous(low="blue3",high="red",trans = 'reverse') + scale_shape_discrete(labels=c("August 2021","December 2021","April 2022"),name="Sample Date") +
  xlab("PC1 [90.02%]") + ylab("PC2 [7.28%]")


#### Using Shapiro-Wilk test for Normality - CLR-Transformed MGM (S) Pathway Data ####
shapiro.test(sulf.path.clr$SOX) # what is the p-value?
# p-value = 0.09133
# p > 0.05 states distribution of data are not significantly different from normal distribution
# p < 0.05 means that data is significantly different from a normal distribution
par(mar = c(1,1,1,1)) # Set the margin on all sides to 1
hist(sulf.path.clr$SOX, col="blue") # with outliars

# visualize Q-Q plot for alpha div
# The Q-Q plot, or quantile-quantile plot, is a graphical tool to help us assess if a set of data plausibly came from some theoretical distribution such as a normal or exponential.
# For example, if we run a statistical analysis that assumes our residuals are normally distributed, we can use a normal Q-Q plot to check that assumption
# more on Q-Q plots here: https://data.library.virginia.edu/understanding-q-q-plots/
# more here too: https://grodri.github.io/glms/notes/c2s9#:~:text=8%20The%20Q%2DQ%20Plot,versus%20quantiles%20of%20a%20distribution.
qqnorm(sulf.path.clr$SOX, pch = 1, frame = FALSE)
qqline(sulf.path.clr$SOX, col = "red", lwd = 2)

shapiro.test(sulf.path.clr$`Dissimilatory Sulfate Redox`) # what is the p-value?
# p-value = 0.002686
# p > 0.05 states distribution of data are not significantly different from normal distribution
# p < 0.05 means that data is significantly different from a normal distribution
par(mar = c(1,1,1,1)) # Set the margin on all sides to 1
hist(sulf.path.clr$`Dissimilatory Sulfate Redox`, col="blue")

# visualize Q-Q plot for species richness
qqnorm(sulf.path.clr$`Dissimilatory Sulfate Redox`, pch = 1, frame = FALSE) # with outliars
qqline(sulf.path.clr$`Dissimilatory Sulfate Redox`, col = "red", lwd = 2)

qqnorm(sulf.path.clr$`Dissimilatory Sulfate Redox`, pch = 1, frame = FALSE) # without outliars
qqline(sulf.path.clr$`Dissimilatory Sulfate Redox`, col = "red", lwd = 2)


#### Merge S Pathway & Metadata ####
s.path.meta<-merge(sulf.path.clr,mgm_meta,by="SampleID")
head(s.path.meta)

# subset CLR-transformed S pathway coverages by date
aug.spath<-s.path.meta[s.path.meta$SampDate=="August.2021",]
dec.spath<-s.path.meta[s.path.meta$SampDate=="December.2021",]
apr.spath<-s.path.meta[s.path.meta$SampDate=="April.2022",]

#### Compare Means of Relative Coverage - by Pathway ####
head(sulf.path.clr)

# First run individual t-tests by time point
# SOX vs Ass.Sulf Reduction
sox.assim<-kruskal.test(s.path.meta$SOX,s.path.meta$`Assimilatory Sulfate Reduction`)
sox.assim

# SOX vs Diss.Sulf Reduction
sox.dissim<-kruskal.test(s.path.meta$SOX,s.path.meta$`Dissimilatory Sulfate Redox`)
sox.dissim

# Combine the p-values
wiltest.Div.pvals<-c(wilcox.test.a21.d21$p.value, wilcox.test.a21.a22$p.value, wilcox.test.d21.a22$p.value)

# Adjust the p-values based on the # of comparisons you did
p.adjust(wiltest.Div.pvals, method="bonferroni",n=3)
# ^ matches findings on the figure, which uses the t_test function from the rstatix package (see geom_pwc())

#### Compare Variance - SOX Pathway ####
# use the following statisitcal tests for variance comparisons

fit1<-kruskal.test(SOX ~ SampDate, data=s.path.meta)
# ANOVA is basically a regression but w/ categorical variables more info here https://www.statology.org/anova-vs-regression/
#pairwise.adonis(sulf.path.clr$SOX, sulf.path.clr$SampDate, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different
fit1

p.adjust(summary(fit1)[[1]][["Pr(>F)"]][1],method="bonferroni")

# Instead of using Tukey test, we can use Dunn's test to see which groups significantly vary if Kruskal-Wallis test is significant
# ANOVA + Tukey for normally distributed data, Kruskal-Wallis + Dunn's test for non-normal data
rstatix::dunn_test(s.path.meta, SOX ~ SampDate, p.adjust.method = "bonferroni", detailed = TRUE)

## The Fligner-Killeen test is a non-parametric test for homogeneity of group variances based on ranks. It is useful when the data are non-normally distributed or when problems related to outliers in the dataset cannot be resolved.
### Fligner's test is a Levene's test for data that are not normally distributed
### It is also one of the many tests for homogeneity of variances which is most robust against departures from normality.
## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
## more here: https://www.geeksforgeeks.org/fligner-killeen-test-in-r-programming/
fligner.test(SOX ~ SampDate, data = s.path.meta)
# Fligner-Killeen:med chi-squared = 1.7395, df = 2, p-value = 0.4191
# Which shows that the data DO NOT deviate significantly from homogeneity.

compare_means(SOX ~ SampDate, data=s.path.meta, method="kruskal.test",p.adjust.method = "bonferroni") #

compare_means(SOX ~ SampDate, data=s.path.meta, method="wilcox.test",p.adjust.method = "bonferroni") #
compare_means(SOX ~ SampDate, data=s.path.meta, method="t.test",p.adjust.method = "bonferroni") #

ggplot(s.path.meta, aes(x=SampDate, y=SOX)) +geom_jitter(aes(color=as.numeric(as.character(Depth_m))), size=4, width=0.15, height=0) +
  scale_colour_gradient2(low="red",mid="pink",high="blue3",midpoint=5,guide = guide_colourbar(reverse = TRUE)) +
  geom_boxplot(fill=NA, outlier.color=NA)+scale_x_discrete(labels=c("August 2021","December 2021","April 2022"))+theme_bw()+theme_classic()+
  labs(title = "", subtitle="", x="Sample Date", y="SOX Relative Coverage (Summed)", color="Depth (m)")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45,hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))

#### Compare Variance - Dissimilatory Sulfate Redox Pathway ####
# shannon diversity is normally distributed; used rarefied counts to calculate ShanDiv
# use the following statisitcal tests for variance comparisons
## ANOVA: are variances significantly different between groups
## Tukey test: which groups' variances are significant different from one another
## Levene's test: is variance homogenous aka equal across samples?

fit2<-kruskal.test(`Dissimilatory Sulfate Redox` ~ SampDate, data=s.path.meta)
# ANOVA is basically a regression but w/ categorical variables more info here https://www.statology.org/anova-vs-regression/
#pairwise.adonis(sulf.path.clr$`Dissimilatory Sulfate Redox`, sulf.path.clr$SampDate, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different
fit2

# Instead of using Tukey test, we can use Dunn's test to see which groups significantly vary if Kruskal-Wallis test is significant
# ANOVA + Tukey for normally distributed data, Kruskal-Wallis + Dunn's test for non-normal data
rstatix::dunn_test(s.path.meta, `Dissimilatory Sulfate Redox` ~ SampDate, p.adjust.method = "bonferroni", detailed = TRUE)

## The Fligner-Killeen test is a non-parametric test for homogeneity of group variances based on ranks. It is useful when the data are non-normally distributed or when problems related to outliers in the dataset cannot be resolved.
### Fligner's test is a Levene's test for data that are not normally distributed
### It is also one of the many tests for homogeneity of variances which is most robust against departures from normality.
## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
## more here: https://www.geeksforgeeks.org/fligner-killeen-test-in-r-programming/
fligner.test(`Dissimilatory Sulfate Redox` ~ SampDate, data = s.path.meta)
# Fligner-Killeen:med chi-squared = 0.1724, df = 2, p-value = 0.9174
# Which shows that the data DO NOT deviate significantly from homogeneity.

compare_means(`Dissimilatory Sulfate Redox` ~ SampDate, data=s.path.meta, method="kruskal.test",p.adjust.method = "bonferroni") #

compare_means(`Dissimilatory Sulfate Redox` ~ SampDate, data=s.path.meta, method="wilcox.test",p.adjust.method = "bonferroni") #
compare_means(`Dissimilatory Sulfate Redox` ~ SampDate, data=s.path.meta, method="t.test",p.adjust.method = "bonferroni") #

ggplot(s.path.meta, aes(x=SampDate, y=`Dissimilatory Sulfate Redox`)) +geom_jitter(aes(color=as.numeric(as.character(Depth_m))), size=4, width=0.15, height=0) +
  scale_colour_gradient2(low="red",mid="pink",high="blue3",midpoint=5,guide = guide_colourbar(reverse = TRUE)) +
  geom_boxplot(fill=NA, outlier.color=NA)+scale_x_discrete(labels=c("August 2021","December 2021","April 2022"))+theme_bw()+theme_classic()+
  labs(title = "", subtitle="", x="Sample Date", y="Dissim SO4 Redox Relative Coverage (Summed)", color="Depth (m)")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=45,hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))


#### Compare Variance - Assimilatory Sulfate Reduction Pathway ####
# shannon diversity is normally distributed; used rarefied counts to calculate ShanDiv
# use the following statisitcal tests for variance comparisons
## ANOVA: are variances significantly different between groups
## Tukey test: which groups' variances are significant different from one another
## Levene's test: is variance homogenous aka equal across samples?

fit3<-kruskal.test(`Assimilatory Sulfate Reduction` ~ SampDate, data=s.path.meta)
# ANOVA is basically a regression but w/ categorical variables more info here https://www.statology.org/anova-vs-regression/
#pairwise.adonis(sulf.path.clr$`Assimilatory Sulfate Reduction`, sulf.path.clr$SampDate, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different
fit3

p.adjust(summary(fit1)[[1]][["Pr(>F)"]][1],method="bonferroni")

# Instead of using Tukey test, we can use Dunn's test to see which groups significantly vary if Kruskal-Wallis test is significant
# ANOVA + Tukey for normally distributed data, Kruskal-Wallis + Dunn's test for non-normal data
rstatix::dunn_test(s.path.meta, `Assimilatory Sulfate Reduction` ~ SampDate, p.adjust.method = "bonferroni", detailed = TRUE)

## The Fligner-Killeen test is a non-parametric test for homogeneity of group variances based on ranks. It is useful when the data are non-normally distributed or when problems related to outliers in the dataset cannot be resolved.
### Fligner's test is a Levene's test for data that are not normally distributed
### It is also one of the many tests for homogeneity of variances which is most robust against departures from normality.
## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
## more here: https://www.geeksforgeeks.org/fligner-killeen-test-in-r-programming/
fligner.test(`Assimilatory Sulfate Reduction` ~ SampDate, data = s.path.meta)
# Fligner-Killeen:med chi-squared = 1.7395, df = 2, p-value = 0.4191
# Which shows that the data DO NOT deviate significantly from homogeneity.

compare_means(`Assimilatory Sulfate Reduction` ~ SampDate, data=s.path.meta, method="kruskal.test",p.adjust.method = "bonferroni") #

compare_means(`Assimilatory Sulfate Reduction` ~ SampDate, data=s.path.meta, method="wilcox.test",p.adjust.method = "bonferroni") #
compare_means(`Assimilatory Sulfate Reduction` ~ SampDate, data=s.path.meta, method="t.test",p.adjust.method = "bonferroni") #


