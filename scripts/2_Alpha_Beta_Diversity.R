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
})

#### Load Global Env to Import Count/ASV Tables ####
load("data/SSeawater_Data_Ready.Rdata") # save global env to Rdata file
#save.image("data/Env_Seqs_All/env.seq_analysis.Rdata") # save global env to Rdata file
bac.dat.all[1:4,1:4]
bac.ASV_table[1:4,1:4]
bac.ASV_table[(nrow(bac.ASV_table)-4):(nrow(bac.ASV_table)),(ncol(bac.ASV_table)-4):(ncol(bac.ASV_table))] # last 4 rows & cols

### Alpha Diversity & Species Richness ####
### Rarefaction Curves
# bacteria/archaea

# Species Accumulation Curve
sc2<-specaccum(bac.ASV_table[,-1],"random")
plot(sc2, ci.type="poly", col="darkgreen", lwd=2, ci.lty=0, ci.col="lightgreen")
boxplot(sc2, col="yellow", add=TRUE, pch=20)

# Prep for Rarefaction Curve
rowSums(bac.ASV_table[,-1]) # total # ASVs per sample, excluding SampleID from calculation
sort(colSums(bac.ASV_table[,-1]))

# Create Rarefaction curve
png('results/Env_Seqs_All/rarecurve.png')
rarecurve(as.matrix(bac.ASV_table[,-1]),col=metadata$Sample_Color, step=1000, label=T,ylab="ASVs")
# to show sampel labels per curve, change label=T
dev.off()

# ASVs per Sample
total_asvs<-data.frame(ASV_Total=rowSums(bac.ASV_table[,-1]),metadata)

ggplot(data=total_asvs, aes(x=SampleID, y=ASV_Total,fill=Sample_Type)) +
  geom_bar(stat="identity",colour="black")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(data=total_asvs, aes(x=Sample_Type, y=ASV_Total)) +
  geom_bar(stat="identity")

# ASVs per Sample Type
asv.sampletype <- as.data.frame(dcast(bac.dat.meta, ASV_ID~Sample_Type, value.var="Count", fun.aggregate=sum)) ###
head(asv.sampletype)
## average ASV per sample type
aggregate(bac.dat.meta$Count, list(bac.dat.meta$Sample_Type), FUN=mean)

# if you have another package loaded that has a diversity function, you can specify that you want to use vegan's diversity function as shown below
Shan_ent.16s<-vegan::diversity(bac.ASV_table[,-1], index="shannon") # Shannon entropy
Shan_div.16s<- exp(Shan_ent.16s) # Shannon Diversity aka Hill number 1

# create data frame with Shannon entropy and Shannon diversity values
div_16s<-data.frame(Bac_Shannon_Entropy=Shan_ent.16s,Bac_Shannon_Diversity=Shan_div.16s)
class(div_16s)
div_16s$SampleID<-rownames(div_16s)
head(div_16s)

# create a data frame with species richness
S_16s<-data.frame(Bac_Species_Richness=specnumber(bac.ASV_table[,-1]), SampleID=rownames(bac.ASV_table)) # finds # of species per sample using RAW count data; if MARGIN = 2 it finds frequencies of species

# merge richness and diversity dataframes together
d.r_16s<-merge(div_16s, S_16s, by.x="SampleID", by.y="SampleID")

# merge w/ metadata
bac.div.metadat <- merge(d.r_16s,metadata, by.x="SampleID", by.y="SampleID")
head(bac.div.metadat)
class(bac.div.metadat) # want data frame

unique(bac.div.metadat$Sample_Type) # see how many elements there are in the Group variable
#bac.div.metadat$Sample_Type <- factor(bac.div.metadat$Sample_Type, levels = c("Seawater", "Soil", "Fecal"))

## Shannon Diversity by Sample Type
bac.a.div<-ggplot(bac.div.metadat, aes(x=Sample_Type, y=Bac_Shannon_Diversity, fill=Sample_Type)) +geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual( values=unique(bac.div.metadat$Sample_Color[order(bac.div.metadat$Sample_Type)]), name ="Sample Type")+theme_classic()+
  labs(title = "Bacterial Shannon Diversity by Sample Type", x="Sample Type", y="Shannon Diversity", fill="Sample Type")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1,,size=10),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))
ggsave(bac.a.div,filename = "figures_copy/Env_Seqs_All_Figs/Bacterial_alpha_diversity_envseqs_4.20.22.png", width=13, height=10, dpi=600)

#colorset1 = melt(c(Dust="#ca6702",Seawater="#168aad",Lung="#c9184a"))

bac.b.div<-ggplot(bac.div.metadat[which(bac.div.metadat$Sample_Type=="Seawater"),], aes(x=SampleID, y=Bac_Shannon_Diversity)) +geom_point(aes(color=Sample_Type),size=4)+theme_bw()+
  labs(title = "Bacterial Shannon Diversity in Seawater", x="Sample ID", y="Shannon Diversity", color="Sample Type")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.position="none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=10))+
  guides(legend="none")+scale_color_manual(values=c("#1f547b"))
ggsave(bac.b.div,filename = "figures_copy/Env_Seqs_All_Figs/Bacterial_alpha_diversity_seawater_4.20.22.png", width=13, height=10, dpi=600)

bac.c.div<-ggplot(bac.div.metadat[which(bac.div.metadat$Sample_Type=="Dust"),], aes(x=SampleID, y=Bac_Shannon_Diversity)) +geom_point(aes(color=factor(Sample_Type)), size=4)+theme_bw()+
  labs(title = "Bacterial Shannon Diversity in Dust", x="Sample ID", y="Shannon Diversity", color="Sample Type")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.position="none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=10))+
  scale_color_manual(values=c("#582707"))+guides(legend="none")
ggsave(bac.c.div,filename = "figures_copy/Env_Seqs_All_Figs/Bacterial_alpha_diversity_dust_4.20.22.png", width=13, height=10, dpi=600)

## Species Richness by Sample Type
bac.a.sr<-ggplot(bac.div.metadat, aes(x=Sample_Type, y=Bac_Species_Richness, fill=Sample_Type)) +geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual( values=unique(bac.div.metadat$Sample_Color[order(bac.div.metadat$Sample_Type)]), name ="Sample Type")+theme_classic()+
  labs(title = "Bacterial Species Richness by Sample Type", x="Sample Type", y="Species Richness", fill="Sample Type")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1,,size=10),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))
ggsave(bac.a.sr,filename = "figures_copy/Env_Seqs_All_Figs/Bacterial_species_richness_envseqs_4.22.22.png", width=13, height=10, dpi=600)

#colorset1 = melt(c(Dust="#ca6702",Seawater="#168aad",Lung="#c9184a"))

bac.b.sr<-ggplot(bac.div.metadat[which(bac.div.metadat$Sample_Type=="Seawater"),], aes(x=SampleID, y=Bac_Species_Richness)) +geom_point(aes(color=Sample_Type),size=4)+theme_bw()+
  labs(title = "Bacterial Species Richness in Seawater", x="Sample ID", y="Species Richness", color="Sample Type")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.position="none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=10))+
  guides(legend="none")+scale_color_manual(values=c("#1f547b"))
ggsave(bac.b.sr,filename = "figures_copy/Env_Seqs_All_Figs/Bacterial_species_richness_seawater_4.22.22.png", width=13, height=10, dpi=600)

bac.c.sr<-ggplot(bac.div.metadat[which(bac.div.metadat$Sample_Type=="Dust"),], aes(x=SampleID, y=Bac_Species_Richness)) +geom_point(aes(color=factor(Sample_Type)), size=4)+theme_bw()+
  labs(title = "Bacterial Species Richness in Dust", x="Sample ID", y="Species Richness", color="Sample Type")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.position="none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=10))+
  scale_color_manual(values=c("#582707"))+guides(legend="none")
ggsave(bac.c.sr,filename = "figures_copy/Env_Seqs_All_Figs/Bacterial_species_richness_dust_4.22.22.png", width=13, height=10, dpi=600)

## Using Shapiro-Wilk test for normality
shapiro.test(bac.div.metadat$Bac_Species_Richness) # what is the p-value?
# my p-value was p-value =  1.579e-09
# p > 0.05 states distribution of data are not significantly different from normal distribution
# p < 0.05 means that data is significantly different from a normal distribution

# visualize Q-Q plot for species richness
qqnorm(bac.div.metadat$Bac_Species_Richness, pch = 1, frame = FALSE)
qqline(bac.div.metadat$Bac_Species_Richness, col = "steelblue", lwd = 2)

#### Beta Diversity ####

# turning bac.ASV_counts_no.contam into phyloseq object called ASV
dim(bac.ASV_table)
ASV<-otu_table(as.matrix(bac.ASV_table[,-1]), taxa_are_rows = TRUE)
head(ASV)
class(ASV) # phyloseq otu_table object

# CLR transformation on phyloseq object
asv_clr<-microbiome::transform(ASV, "clr")
head(asv_clr)

# create CLR Sample x Species matrix for input into dist()
b.clr<-as.matrix(t(asv_clr))
rownames(b.clr)
# calculate our Euclidean distance matrix using CLR data
b.euc_dist <- dist(b.clr, method = "euclidean")

# creating our hierarcical clustering dendrogram
b.euc_clust <- hclust(b.euc_dist, method="ward.D2")

# let's make it a little nicer...
b.euc_dend <- as.dendrogram(b.euc_clust, hang=0.2)
b.dend_cols <- as.character(metadata$SampDate_Color[order.dendrogram(b.euc_dend)])
labels_colors(b.euc_dend) <- b.dend_cols

plot(b.euc_dend, ylab="CLR Euclidean Distance",cex = 0.5) + title(main = "Bacteria/Archaea Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 3)
legend("topright",legend = c("Seawater","Soil","Dust","Playa","Fecal","Lung","Control"),cex=.8,col = c("#1f547b","#c44536","#432818","#d00000","#66615f","#47126b","#b13d1e"),pch = 15, bty = "n")
# Control is dark blue ("#218380"), #Alternaria is light blue ("#73d2de")
dev.off()

# let's use our Euclidean distance matrix from before
b.pcoa <- pcoa(b.euc_dist)

# The proportion of variances explained is in its element values$Relative_eig
b.pcoa$values

# extract principal coordinates
b.pcoa.vectors<-data.frame(b.pcoa$vectors)
b.pcoa.vectors$SampleID<-rownames(b.pcoa$vectors)

# merge pcoa coordinates w/ metadata
b.pcoa.meta<-merge(b.pcoa.vectors, metadata, by.x="SampleID", by.y="SampleID")
head(b.pcoa.meta)

b.pcoa$values# pull out Relative (relative eigen) variation % to add to axes labels

# create PCoA ggplot fig
pcoa1<-ggplot(b.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Sample_Type)), size=4)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea",subtitle="Using Centered-Log Ratio Data",xlab="Axis 1 [41.14%]", ylab="Axis 2 [9.04%]",color="Sample Type")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+scale_color_manual(name ="Sample Type",values=unique(b.pcoa.meta$Sample_Color[order(b.pcoa.meta$Sample_Type)])) +xlab("Axis 1 [49.37%]") + ylab("Axis 2 [14.28%]")

#colorset1 = melt(c(Dust="#ca6702",Seawater="#168aad",Lung="#c9184a"))

ggsave(pcoa1,filename = "figures_copy/Env_Seqs_All_Figs/16S_pcoa_CLR.png", width=12, height=10, dpi=600)

ggplot(b.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Sample_Type)), size=4)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea",subtitle="Using Centered-Log Ratio Data",xlab="Axis 1", ylab="Axis 2",color="Sample Type")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+scale_color_manual(name ="Sample Type",values=unique(b.pcoa.meta$Sample_Color[order(b.pcoa.meta$Sample_Type)])) +
  xlab("Axis 1 [41.14%]") + ylab("Axis 2 [9.04%]")+
  geom_text(aes(label=SampleID),hjust=0, vjust=0)

## bac.ASV_counts_no.contam[,-ncol(bac.ASV_counts_no.contam)] allows us to drop the last column in the data frame, which in this case is a column of ASV IDs

## betadisper to look at within group variance
b.disper<-betadisper(b.euc_dist, metadata$Sample_Type)
b.disper

permutest(b.disper, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons

anova(b.disper) # p = 0.0001394 --> reject the Null H, spatial medians are significantly difference across Categories

TukeyHSD(b.disper) # tells us which Category's dispersion MEANS are significantly different than each other

# Visualize dispersions
png('pcoa_betadispersion.png',width = 700, height = 600, res=100)
plot(b.disper,main = "Centroids and Dispersion based on Aitchison Distance", col=colorset1$Sample_Color)
dev.off()

png('boxplot_centroid_distance.png',width = 700, height = 600, res=100)
boxplot(b.disper,xlab="Sample Type", main = "Distance to Centroid by Category", sub="Based on Aitchison Distance", col=colorset1$Sample_Color)
dev.off()

### Seawater only PCoA
sw.asv.meta<-subset(sw.asv.meta, select=-c(Exposure_Duration, Exposure_Type))
sw.asv.meta2<-subset(sw.asv.meta, SampleYear!="2020") # drop 2020 samples
sw.asv.melt<-melt(sw.asv.meta2, id.vars=c("SampleID","Sample_Type","SampleMonth","SampleYear","Depth_m","SampleSource","Deployment","SampleOrControl","ExtractionMethod","LysisType","Sample_Color"))
head(sw.asv.melt)
colnames(sw.asv.melt)[which(names(sw.asv.melt) == "value")] <- "Counts"
colnames(sw.asv.melt)[which(names(sw.asv.melt) == "variable")] <- "ASV_ID"
head(sw.asv.melt)

sw_counts <- as.data.frame(dcast(sw.asv.melt, ASV_ID~SampleID, value.var="Counts", fun.aggregate=sum)) ###
head(sw_counts) # counts by phyla per sample
#sw_counts<-sw_counts[,-which(colnames(sw_counts)=="SCityWater")] ## drop non Salton Seawater samples
rownames(sw_counts)<-sw_counts$ASV_ID

sw.ASV<-otu_table(as.matrix(sw_counts[,-1]), taxa_are_rows = TRUE) # create phyloseq otu table
head(sw.ASV)
class(sw.ASV) # phyloseq otu_table object

# CLR transformation on phyloseq object
sw.asv_clr<-microbiome::transform(sw.ASV, "clr")
head(sw.asv_clr)

# create CLR Sample x Species matrix for input into dist()
sw.clr<-as.matrix(t(sw.asv_clr))

# calculate our Euclidean distance matrix using CLR data
sw.euc_dist <- dist(sw.clr, method = "euclidean")

# creating our hierarcical clustering dendrogram
sw.euc_clust <- hclust(sw.euc_dist, method="ward.D2")

# let's make it a little nicer...
sw.euc_dend <- as.dendrogram(sw.euc_clust, hang=0.2)
sw.dend_cols <- as.character(sw.asv.meta$Sample_Color[order.dendrogram(sw.euc_dend)])
labels_colors(sw.euc_dend) <- sw.dend_cols

#png(file="figures_copy/Env_Seqs_All_Figs/16S_CLR_cluster_SampleType.png",width = 1000, height = 1100, res=90)
pdf(file="figures_copy/Env_Seqs_All_Figs/16S_CLR_cluster_Seawater.pdf", width=8,height=10)
par(cex=.7)
plot(sw.euc_dend, ylab="CLR Euclidean Distance",cex = 0.5) + title(main = "Bacteria/Archaea Clustering Dendrogram", sub = "Only Includes Seawater Samples",cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 3)
legend("topright",legend = c("Seawater"),cex=.8,col = c("#004ACE"),pch = 15, bty = "n")
# Control is dark blue ("#218380"), #Alternaria is light blue ("#73d2de")
dev.off()
#colorset1 = melt(c(Dust="#ca6702",Seawater="#168aad",Lung="#c9184a"))

# let's use our Euclidean distance matrix from before
sw.pcoa <- pcoa(sw.euc_dist)

# The proportion of variances explained is in its element values$Relative_eig
sw.pcoa$values

# extract principal coordinates
sw.pcoa.vectors<-data.frame(sw.pcoa$vectors)
sw.pcoa.vectors$SampleID<-rownames(sw.pcoa$vectors)

# merge pcoa coordinates w/ metadata
sw.meta<-unique(subset(sw.asv.melt, select=-c(ASV_ID, Counts)))
sw.pcoa.meta<-merge(sw.pcoa.vectors, sw.meta, by.x="SampleID", by.y="SampleID")
head(sw.pcoa.meta)
sw.pcoa.meta$Depth.num<-as.numeric(sw.pcoa.meta$Depth_m) # turn Elevation factor column into a numeric version for continuous color scheme -- for plot!
sw.pcoa.meta$SampleMonth<-factor(sw.pcoa.meta$SampleMonth, levels=c("June","August","December","April"))

sw.pcoa$values # pull out relative variation % to add to axes labels

# create PCoA ggplot fig
pcoa.sw1<-ggplot(sw.pcoa.meta, aes(x=Axis.1, y=Axis.2), color=Depth.num, show.legend = TRUE)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Seawater",subtitle="Using Centered-Log Ratio Data",xlab="Axis 1", ylab="Axis 2")+theme_classic()+ geom_text(aes(label=SampleID),fontface=2,check_overlap=FALSE, size=4)+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1))+guides(shape = guide_legend(override.aes = list(size = 5)))+
  xlab("Axis 1 [18.32%]") + ylab("Axis 2 [14.25%]")+xlim(-125,125)+saturation(scale_colour_gradientn(colours=fair_cols), 0.9)

#colorset1 = melt(c(Dust="#ca6702",Seawater="#168aad",Lung="#c9184a"))

ggsave(pcoa.sw1,filename = "figures_copy/Env_Seqs_All_Figs/16S_seawater.only_pcoa_CLR.png", width=12, height=10, dpi=600)

# sample month shape, depth color
pcoa.sw2<-ggplot(sw.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(color=Depth.num,shape=SampleMonth), size=5)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Seawater",subtitle="Using Centered-Log Ratio Data",xlab="Axis 1", ylab="Axis 2",color="Depth (m)")+
  theme_classic()+ theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),axis.text = element_text(size=12),axis.text.x = element_text(vjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17))+
  scale_color_continuous(low="blue3",high="red",trans = 'reverse') + scale_shape_discrete(labels=c("June 2021","August 2021","December 2021","April 2022"),name="Sample Month") +
  xlab("Axis 1 [18.32%]") + ylab("Axis 2 [14.25%]")
ggsave(pcoa.sw2,filename = "figures_copy/Env_Seqs_All_Figs/16S_seawater.only_pcoa_all_10.12.22.png", width=12, height=10, dpi=600)

ggplot(sw.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(color=Depth.num,shape=SampleMonth), size=5)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Seawater",subtitle="Using Centered-Log Ratio Data",xlab="Axis 1", ylab="Axis 2",color="Depth (m)")+
  theme_classic()+ theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),axis.text = element_text(size=12),axis.text.x = element_text(vjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17))+
  scale_color_continuous(low="blue3",high="red",trans = 'reverse') +
  xlab("Axis 1 [18.32%]") + ylab("Axis 2 [14.25%]")

# depth
pcoa.sw3<-ggplot(sw.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(color=factor(Sample_Type),shape=as.factor(Depth_m)), size=4)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Seawater",subtitle="Using Centered-Log Ratio Data",xlab="Axis 1", ylab="Axis 2",color="Sample Type")+
  theme_classic()+ theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),axis.text = element_text(size=12),axis.text.x = element_text(vjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17))+
  scale_color_manual(values=unique(sw.pcoa.meta$Sample_Color[order(sw.pcoa.meta$Sample_Type)])) + scale_shape_manual(values=c(0,1,2,6,5,9,15,17,16),name="Depth (m)") +
  xlab("Axis 1 [18.59%]") + ylab("Axis 2 [14.73%]")
ggsave(pcoa.sw3,filename = "figures_copy/Env_Seqs_All_Figs/16S_seawater.only_pcoa_depth.png", width=12, height=10, dpi=600)

