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

#load("data/SSW_analysis.Rdata") # load Rdata to global env
#save.image("data/SSW_analysis.Rdata") # save global env to Rdata file

#### Import and Prepare Data for Analyses ####

## NOTES ABOUT DATA:
# eukaryotic counts have been removed; contaminant counts removed with decontam(); all non-Lyons-samples excluded from data set before importing into this script
# zero counts & singletons also removed

## Import ALL env plate bacterial ASV count data
bac.ASV_all<-data.frame(readRDS("data/SaltonSeawater_16S_AllData_Robject.rds", refhook = NULL))
dim(bac.ASV_all) ## has count info, ASV IDs, taxa info, and original metadata
bac.ASV_all[1:4,1:4]

# Create ASV table
bac.ASV_table <- as.data.frame(dcast(bac.ASV_all, SampleID~ASV_ID, value.var="Count", fun.aggregate=sum)) ###
head(bac.ASV_table) # counts by asvs per sample
rownames(bac.ASV_table)<-bac.ASV_table$SampleID
bac.ASV_table[1:4,1:4]

# Create taxa table
bac.ASV_taxa<-as.data.frame(unique(subset(bac.ASV_all, select=c(ASV_ID,Kingdom, Phylum, Class, Order, Family, Genus, Species))))

#### Update Metadata ####
# upload geochem data from Lyons lab

chem_meta<-as.data.frame(read_xlsx("data/SaltonSeawater_Lyons_Aronson_Metadata_All.xlsx", sheet="Variables_of_Interest"))
head(chem_meta)
# create color variable(s) to identify variables by colors
## color for sample type
meta1<-unique(subset(bac.ASV_all, select=c(SampleID,Sample_Type,SampleMonth,SampleYear,Depth_m,SampleSource)))
head(meta1)
dim(meta1)
rownames(meta1)<-meta1$SampleID

metadata<-merge(meta1, chem_meta, by=c("SampleID","Sample_Type","SampleMonth","SampleYear", "Depth_m", "SampleSource"))
head(metadata)

metadata$Depth_m<-factor(metadata$Depth_m, levels=c("0","2","3","4","5","7","9","10","11"))

unique(metadata$SampleMonth)
metadata$SampleMonth<-factor(metadata$SampleMonth, levels=c("June","August","December","April"))

metadata$SampDate<-interaction(metadata$SampleMonth, metadata$SampleYear)
head(metadata)

cold2warm1<-get_palette(paste0("#",c("252A52", "66ADE5", "FFC465","BF1B0B")),k=10)
names(cold2warm1) <- levels(metadata$Depth_m)

fair_cols <- paste0("#",c("252A52", "66ADE5", "FFC465","BF1B0B"))
names(fair_cols) <- letters[1:4]
fair_ramp <- scales::colour_ramp(fair_cols)
fair_sat <- saturation(fair_ramp, 1)

# save.image("data/SSW_analysis.Rdata")

### Merge Metadata & Count Data Together ####
bac.dat.all<-merge(bac.ASV_all, metadata, by=c("SampleID","Sample_Type","SampleMonth","SampleYear", "Depth_m", "SampleSource"))
#rownames(bac.dat.all)<-bac.dat.all$SampleID
bac.dat.all[1:4,1:4]
dim(bac.dat.all)

### Export Global Env for Other Scripts ####
save.image("data/SSeawater_Data_Ready.Rdata")

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

unique(bac.div.metadat$Depth_m) # see how many elements there are in the Group variable
#bac.div.metadat$Sample_Type <- factor(bac.div.metadat$Sample_Type, levels = c("Seawater", "Soil", "Fecal"))

## Shannon Diversity by depth
bac.a.div<-ggplot(bac.div.metadat, aes(x=Depth_m, y=Bac_Shannon_Diversity)) +geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual( values=unique(bac.div.metadat$Sample_Color[order(bac.div.metadat$Sample_Type)]), name ="Sample Type")+theme_classic()+
  labs(title = "Bacterial Shannon Diversity by Sample Type", x="Depth", y="Shannon Diversity", fill="Sample Type")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1,,size=10),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))
ggsave(bac.a.div,filename = "figures_copy/Env_Seqs_All_Figs/Bacterial_alpha_diversity_envseqs_4.20.22.png", width=13, height=10, dpi=600)

## Shannon Diversity by month
bac.b.div<-ggplot(bac.div.metadat, aes(x=SampleMonth, y=Bac_Shannon_Diversity)) +geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual( values=unique(bac.div.metadat$Sample_Color[order(bac.div.metadat$Sample_Type)]), name ="Sample Type")+theme_classic()+
  labs(title = "Bacterial Shannon Diversity by Sample Type", x="Depth", y="Shannon Diversity", fill="Sample Type")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1,,size=10),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))
ggsave(bac.a.div,filename = "figures_copy/Env_Seqs_All_Figs/Bacterial_alpha_diversity_envseqs_4.20.22.png", width=13, height=10, dpi=600)

#colorset1 = melt(c(Dust="#ca6702",Seawater="#168aad",Lung="#c9184a"))

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
## Seawater only PCoA
ssw.meta2<-subset(ssw.meta, SampleYear!="2020") # drop 2020 samples
sw.asv.melt<-melt(ssw.meta2, id.vars=c("SampleID","Sample_Type","SampleMonth","SampleYear","Depth_m","SampleSource"))
head(sw.asv.melt)
colnames(sw.asv.melt)[which(names(sw.asv.melt) == "value")] <- "Counts"
colnames(sw.asv.melt)[which(names(sw.asv.melt) == "variable")] <- "ASV_ID"
head(sw.asv.melt)

sw_counts <- as.data.frame(dcast(sw.asv.melt, ASV_ID~SampleID, value.var="Counts", fun.aggregate=sum)) ###
head(sw_counts) # counts by asvs per sample
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
sw.dend_cols <- as.character("#1f547b")
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

sw.pcoa$values # pull out relative variation % to add to axes labels

# create PCoA ggplot fig
pcoa.sw1<-ggplot(sw.pcoa.meta, aes(x=Axis.1, y=Axis.2), color=Depth.num, show.legend = TRUE)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Seawater",subtitle="Using Centered-Log Ratio Data",xlab="Axis 1", ylab="Axis 2")+theme_classic()+ geom_text(aes(label=SampleID),fontface=2,check_overlap=FALSE, size=4)+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1))+guides(shape = guide_legend(override.aes = list(size = 5)))+
  xlab("Axis 1 [18.43%]") + ylab("Axis 2 [14.25%]")+xlim(-125,125)+saturation(scale_colour_gradientn(colours=fair_cols), 0.9)

#colorset1 = melt(c(Dust="#ca6702",Seawater="#168aad",Lung="#c9184a"))

ggsave(pcoa.sw1,filename = "figures_copy/Env_Seqs_All_Figs/16S_seawater.only_pcoa_CLR.png", width=12, height=10, dpi=600)

# sample month shape, depth color
pcoa.sw2<-ggplot(sw.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(color=Depth.num,shape=SampleMonth), size=5)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Seawater",subtitle="Using Centered-Log Ratio Data",xlab="Axis 1", ylab="Axis 2",color="Depth (m)")+
  theme_classic()+ theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),axis.text = element_text(size=12),axis.text.x = element_text(vjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17))+
  scale_color_continuous(low="blue3",high="red",trans = 'reverse') + scale_shape_discrete(labels=c("June 2021","August 2021","December 2021","April 2022"),name="Sample Month") +
  xlab("Axis 1 [18.59%]") + ylab("Axis 2 [14.73%]")
ggsave(pcoa.sw2,filename = "figures_copy/Env_Seqs_All_Figs/16S_seawater.only_pcoa_all_10.11.22.png", width=12, height=10, dpi=600)

# depth
pcoa.sw3<-ggplot(sw.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(color=factor(Sample_Type),shape=as.factor(Depth_m)), size=4)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Seawater",subtitle="Using Centered-Log Ratio Data",xlab="Axis 1", ylab="Axis 2",color="Sample Type")+
  theme_classic()+ theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),legend.title.align=0.5, legend.title = element_text(size=15),axis.text = element_text(size=12),axis.text.x = element_text(vjust=1),legend.text = element_text(size=12),plot.title = element_text(size=17))+
  scale_color_manual(values=unique(sw.pcoa.meta$Sample_Color[order(sw.pcoa.meta$Sample_Type)])) + scale_shape_manual(values=c(0,1,2,6,5,9,15,17,16),name="Depth (m)") +
  xlab("Axis 1 [18.59%]") + ylab("Axis 2 [14.73%]")
ggsave(pcoa.sw3,filename = "figures_copy/Env_Seqs_All_Figs/16S_seawater.only_pcoa_depth.png", width=12, height=10, dpi=600)

#### Relative Abundance -- compare by SampleID, depth, time point ####

#### Genus Relative Abundance ####
head(all_bac)

# by genus + depth
bac.g.typ <- as.data.frame(dcast(all_bac,Depth~Genus, value.var="Counts", fun.aggregate=sum)) ###
head(bac.g.typ) # counts by genus + sample type
rownames(bac.g.typ)<-bac.g.typ$Sample_Type

b.RA_g.typ<-data.frame(decostand(bac.g.typ[,-1], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(b.RA_g.typ) # sanity check
b.RA_g.typ$Sample_Type<-rownames(b.RA_g.typ)
head(b.RA_g.typ)

#melt down relativized data to merge with metadata
b.g.typ_m<-melt(b.RA_g.typ, by="Sample_Type")

head(b.g.typ_m)
colnames(b.g.typ_m)[which(names(b.g.typ_m) == "variable")] <- "Genus"
colnames(b.g.typ_m)[which(names(b.g.typ_m) == "value")] <- "Counts"
head(b.g.typ_m) ## relative abundance based on sum of counts by genus!
b.g.typ_m$Genus<-gsub("X.","", b.g.typ_m$Genus)
b.g.typ_m$Genus<-gsub("(*)\\.\\.(*)","\\1.\\2", b.g.typ_m$Genus)
b.g.typ_m$Genus<-gsub("_"," ", b.g.typ_m$Genus)
head(b.g.typ_m)

# merge metadata and RA data
typ_meta<-unique(data.frame("Sample_Type"=metadata$Sample_Type, "Sample_Color"=metadata$Sample_Color))
g_typ_meta<-merge(typ_meta,b.g.typ_m, by="Sample_Type")
#g_typ_meta<-subset(g_typ_meta, Sample_Type=="Dust" | Sample_Type=="Soil" | Sample_Type=="Seawater")
g_typ_meta<-subset(g_typ_meta, Genus!="Unknown")
g_typ_meta<-subset(g_typ_meta, Counts!=0)

#typ_meta<-unique(data.frame("Sample_Type"=metadata$Sample_Type, "Sample_Color"=metadata$Sample_Color))
#gen_typ_meta<-merge(typ_meta,b.gen.typ_m, by="Sample_Type")

ts1<-ggplot(g_typ_meta, aes(Genus, Counts)) +
  geom_jitter(aes(color=factor(Sample_Type)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Sample Type", values=unique(g_typ_meta$Sample_Color[order(g_typ_meta$Sample_Type)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=0.5,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genus", y="Relative Abundance", title="Microbial Genuses & Sample Type")

ggsave(ts1,filename = "figures_copy/Env_Seqs_All_Figs/16S_Genus.RA_typ.mat.png", width=15, height=10, dpi=600)

typ.1p<-subset(g_typ_meta, Counts>=0.01)

ts2<-ggplot(typ.1p, aes(Genus, Counts)) +
  geom_jitter(aes(color=factor(Sample_Type)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Sample Type", values=unique(typ.1p$Sample_Color[order(typ.1p$Sample_Type)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13, vjust=-1),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genus", y="Relative Abundance", title="Microbial Genuses & Sample Type", subtitle="Includes Taxa with a Relative Abundance of at least 1%")

ggsave(ts2,filename = "figures_copy/Env_Seqs_All_Figs/16S_Genus.RA_1percent_typ.mat.png", width=12, height=10, dpi=600)

env_typ_meta<-subset(g_typ_meta, Sample_Type=="Dust" | Sample_Type=="Soil" | Sample_Type=="Seawater")
typ.1p<-subset(g_typ_meta, Counts>=0.01)

ts3<-ggplot(env_typ_meta, aes(Genus, Counts)) +
  geom_jitter(aes(color=factor(Sample_Type)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Sample Type", values=unique(env_typ_meta$Sample_Color[order(env_typ_meta$Sample_Type)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=0.5,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genus", y="Relative Abundance", title="Microbial Genuses & Sample Type")

ggsave(ts3,filename = "figures_copy/Env_Seqs_All_Figs/16S_Genus.RA_env.typ.mat.png", width=15, height=10, dpi=600)

ts3<-ggplot(env_typ_meta[env_typ_meta$Counts>=0.01,], aes(Genus, Counts)) +
  geom_jitter(aes(color=factor(Sample_Type)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Sample Type", values=unique(env_typ_meta$Sample_Color[order(env_typ_meta$Sample_Type)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=0.5,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genus", y="Relative Abundance", title="Microbial Genuses & Sample Type")

ggsave(ts3,filename = "figures_copy/Env_Seqs_All_Figs/16S_Genus.RA_env.typ.mat.png", width=15, height=10, dpi=600)

ggplot(env_typ_meta[env_typ_meta$Counts>=0.01,], aes(x=Sample_Type, y=Counts, fill=env_typ_meta$Genus[env_typ_meta$Sample_Type=="Dust" | Sample_Type=="Soil" | Sample_Type=="Seawater"]))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genera Relative Abundance", x="SampleID", y="Relative Abundance", fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+guides(fill=guide_legend(ncol=2))+scale_y_continuous(expand = c(0,0),limits = c(0,1))

#### Shared Genus Relative Abundance ####
head(all_bac)

head(b.g.typ_m) # made in section above, ~ line 1277

# merge metadata and RA data
typ_meta<-unique(data.frame("Sample_Type"=metadata$Sample_Type, "Sample_Color"=metadata$Sample_Color))
g_typ_meta<-merge(typ_meta,b.g.typ_m, by="Sample_Type")
g_typ_meta<-subset(g_typ_meta, Sample_Type=="Dust" | Sample_Type=="Soil" | Sample_Type=="Seawater")
g_typ_meta<-subset(g_typ_meta, Genus!="Unknown")
g_typ_meta<-subset(g_typ_meta, Counts!=0)

# finding shared genera...

n_occur <- data.frame(table(g_typ_meta$Genus)) # find frequencies of genera to see which are shared between sample types
n_occur[n_occur$Freq > 1,] # shows us which genera have a greater frequency than 2
g_shared_typ<-g_typ_meta[g_typ_meta$Genus %in% n_occur$Var1[n_occur$Freq > 1],]
#write.csv(g_shared_typ,"16S_Genera_SampleType_Shared.csv",row.names=FALSE)

g_not.shared_typ<-subset(g_typ_meta, !(g_typ_meta$Genus %in% g_shared_typ$Genus)) # subset based off of what is NOT in one dataframe from another data frame

ggplot(g_shared_typ[g_shared_typ$Counts>=0.001,], aes(x=Sample_Type, y=Counts, fill=Genus))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genera Relative Abundance", x="SampleID", y="Relative Abundance", fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+guides(fill=guide_legend(ncol=2))+scale_y_continuous(expand = c(0,0),limits = c(0,1))

sh.t1<-ggplot(g_shared_typ, aes(Genus, Counts)) +
  geom_jitter(aes(color=factor(Sample_Type)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Sample Type", values=unique(g_shared_typ$Sample_Color[order(g_shared_typ$Sample_Type)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13, vjust=-1),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genus", y="Relative Abundance", title="Shared Microbial Genera Across Sample Types")

ggsave(sh.t1,filename = "figures_copy/Env_Seqs_All_Figs/16S_shared_Genera_taxasum_type.png", width=23, height=10, dpi=600)

sh.t1a<-ggplot(g_shared_typ[g_shared_typ$Counts>=0.01,], aes(Genus, Counts)) +
  geom_jitter(aes(color=factor(Sample_Type)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Sample Type", values=unique(g_shared_typ$Sample_Color[order(g_shared_typ$Sample_Type)])) + theme_classic() +
  theme(axis.title.x = element_text(size=13, vjust=-1),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genus", y="Relative Abundance", title="Microbial Genera & Sample Type",subtitle="Only Includes Genera Shared Across All Sample Types")

ggsave(sh.t1a,filename = "figures_copy/Env_Seqs_All_Figs/16S_shared_Genera_taxasum_type_1perc.png", width=23, height=10, dpi=600)

ggplot(g_shared_typ[g_shared_typ$Counts>=0.01,], aes(Genus, Counts)) +
  geom_jitter(aes(color=factor(Sample_Type)), size=3, width=0.15, height=0) +
  scale_color_manual(name ="Sample Type", values=unique(g_shared_typ$Sample_Color[order(g_shared_typ$Sample_Type)])) + theme_classic() +
  geom_boxplot(fill=NA, outlier.color=NA) +
  theme(axis.title.x = element_text(size=13, vjust=-1),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genus", y="Relative Abundance", title="Microbial Genera & Sample Type")

ggplot(g_shared_typ[g_shared_typ$Counts>=0.0005,], aes(Genus, Counts)) +
  geom_jitter(aes(color=factor(Sample_Type)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Sample Type", values=unique(g_shared_typ$Sample_Color[order(g_shared_typ$Sample_Type)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13, vjust=-1),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genus", y="Relative Abundance", title="Microbial Genera & Sample Type")

sh.t2<-ggplot(g_shared_typ[g_shared_typ$Counts>=0.005,], aes(Genus, Counts)) +
  geom_jitter(aes(color=factor(Sample_Type)), size=3, width=0.15, height=0) +
  scale_color_manual(name ="Sample Type", values=unique(g_shared_typ$Sample_Color[order(g_shared_typ$Sample_Type)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13, vjust=-1),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genus", y="Relative Abundance", title="Microbial Genera & Sample Type",subtitle="Relative Abundance > 0.25%")

ggsave(sh.t2,filename = "figures_copy/Env_Seqs_All_Figs/16S_shared_Genera_taxasum_type_0.5perc.png", width=23, height=10, dpi=600)

sh.t3<-ggplot(g_shared_typ[g_shared_typ$Counts>=0.0025,], aes(Genus, Counts)) +
  geom_jitter(aes(color=factor(Sample_Type)), size=3, width=0.15, height=0) +
  scale_color_manual(name ="Sample Type", values=unique(g_shared_typ$Sample_Color[order(g_shared_typ$Sample_Type)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13, vjust=-1),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genus", y="Relative Abundance", title="Microbial Genera & Sample Type",subtitle="Relative Abundance > 0.25%")

ggsave(sh.t3,filename = "figures_copy/Env_Seqs_All_Figs/16S_shared_Genera_taxasum_type_0.25perc.png", width=23, height=10, dpi=600)

## Comparing Genera in Dust vs Seawater
# X Axis Breaks and Labels
lbls = paste0(as.character(c(seq(0.05, 0, -0.01), seq(0.01, 0.05, 0.01)))) # labels
brks=seq(-0.05,0.05,0.01)
g_shared_typ1<-subset(g_shared_typ, Sample_Type!="Soil")
g_shared_typ1$Counts2 <- ifelse(g_shared_typ1$Sample_Type == "Seawater", -1*g_shared_typ1$Counts, g_shared_typ1$Counts)
g_shared_typ1<-g_shared_typ1[order(-g_shared_typ1$Counts2,g_shared_typ1$Genus),]
g_shared_typ1$GenSamp<-interaction(g_shared_typ1$Genus,g_shared_typ1$Sample_Type)
g_shared_typ1$GenSamp<-factor(g_shared_typ1$GenSamp, levels=g_shared_typ1$GenSamp)
class(g_shared_typ1$GenSamp)
g_shared_typ1$Genus<-factor(g_shared_typ1$Genus, levels=unique(g_shared_typ1$Genus[sort(g_shared_typ1$GenSamp)]))

share1<-ggplot(g_shared_typ1, aes(x = Genus, y = -Counts2, fill = Sample_Type)) +
  geom_bar(data = subset(g_shared_typ1[g_shared_typ1$Counts2<=-0.0005,], Sample_Type == "Seawater"), stat = "identity") +
  geom_bar(data = subset(g_shared_typ1[g_shared_typ1$Counts2>=0.0005,], Sample_Type == "Dust"), stat = "identity") +
  coord_flip()+scale_y_continuous(labels = lbls,breaks=brks)+theme_classic()+scale_fill_manual(name ="Sample Type", values=unique(g_shared_typ1$Sample_Color[rev(order(g_shared_typ1$Sample_Type))]))+ylab("Relative Abundance")+
  theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+labs(title="Microbial Genera by Sample Type", subtitle="Includes Taxa with a Relative Abundance of at least 0.05%")

ggsave(share1,filename = "figures_copy/Env_Seqs_All_Figs/16S_shared_Genera_dust.v.sw_population.pyramid.png", width=12, height=10, dpi=600)

pp2<-ggplot(g_shared_typ1[g_shared_typ1$Counts>=0.0005,], aes(x = reorder(Genus,Counts), fill = Sample_Type,y = ifelse(test = Sample_Type == "Dust",yes = Counts, no = -Counts))) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = abs, limits = max(g_shared_typ1$Counts) * c(-1,1)) +
  coord_flip()+scale_fill_manual(name ="Sample Type", values=unique(g_shared_typ1$Sample_Color[order(g_shared_typ1$Sample_Type)]))+ylab("Relative Abundance")+theme_classic()+
  theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
  labs(title="Microbial Genera by Sample Type",subtitle="Includes Only Shared Taxa with a Relative Abundance of at least 0.05%")+xlab("Genus")

ggsave(pp2,filename = "figures_copy/Env_Seqs_All_Figs/16S_shared_Genera_dust.v.sw_population.pyramid.pretty.05percent.png", width=13, height=10, dpi=600)

pp3<-ggplot(g_shared_typ1[g_shared_typ1$Counts>=0.0010,], aes(x = reorder(Genus,Counts), fill = Sample_Type,y = ifelse(test = Sample_Type == "Dust",yes = Counts, no = -Counts))) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = abs, limits = max(g_shared_typ1$Counts) * c(-1,1)) +
  coord_flip()+scale_fill_manual(name ="Sample Type", values=unique(g_shared_typ1$Sample_Color[order(g_shared_typ1$Sample_Type)]))+ylab("Relative Abundance")+theme_classic()+
  theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
  labs(title="Microbial Genera by Sample Type",subtitle="Includes Only Shared Taxa with a Relative Abundance of at least 0.1%")+xlab("Genus")

ggsave(pp3,filename = "figures_copy/Env_Seqs_All_Figs/16S_shared_Genera_dust.v.sw_population.pyramid.pretty.1percent.png", width=13, height=10, dpi=600)

pp4<-ggplot(g_shared_typ1[g_shared_typ1$Counts>=0.010,], aes(x = reorder(Genus,Counts), fill = Sample_Type,y = ifelse(test = Sample_Type == "Dust",yes = Counts, no = -Counts))) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = abs, limits = max(g_shared_typ1$Counts) * c(-1,1)) +
  coord_flip()+scale_fill_manual(name ="Sample Type", values=unique(g_shared_typ1$Sample_Color[order(g_shared_typ1$Sample_Type)]))+ylab("Relative Abundance")+theme_classic()+
  theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
  labs(title="Microbial Genera by Sample Type",subtitle="Includes Only Shared Taxa with a Relative Abundance of at least 1%")+xlab("Genus")
ggsave(pp4,filename = "figures_copy/Env_Seqs_All_Figs/16S_shared_Genera_dust.v.sw_population.pyramid.pretty.1percent_10.14.22.png", width=13, height=10, dpi=600)

## Lollipop chart

lg1<-ggplot(g_shared_typ1[g_shared_typ1$Counts>=0.0005,], aes(x = reorder(Genus,Counts),
                                                              y = ifelse(test = Sample_Type == "Dust",yes = Counts, no = -Counts),color=g_shared_typ1$Sample_Type[g_shared_typ1$Counts>=0.0005])) +
  geom_point(stat='identity',size=3)  +
  geom_segment(aes(y = 0,
                   x = Genus,
                   yend = ifelse(test = Sample_Type == "Dust",yes = Counts, no = -Counts),
                   xend = Genus),color = "black") +
  coord_flip()+scale_color_manual(name ="Sample Type", values=unique(g_shared_typ1$Sample_Color[order(g_shared_typ1$Sample_Type)]))+ylab("Relative Abundance")+theme_classic()+theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
  labs(title="Microbial Genera by Sample Type")+xlab("Genus")

ggsave(lg1,filename = "figures_copy/Env_Seqs_All_Figs/16S_shared_Genera_dust.v.sw_lollipop_chart.png", width=13, height=10, dpi=600)

ggplot(g_shared_typ1[g_shared_typ1$Counts>=0.0005,], aes(x = reorder(Genus,Counts),
                                                         y = ifelse(test = Sample_Type == "Dust",yes = Counts, no = -Counts),color=g_shared_typ1$Sample_Type[g_shared_typ1$Counts>=0.0005])) +
  geom_bar(stat = "identity") +
  coord_flip()+scale_color_manual(name ="Sample Type", values=unique(g_shared_typ1$Sample_Color[order(g_shared_typ1$Sample_Type)]))+ylab("Relative Abundance")+theme_classic()+theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
  labs(title="Microbial Genera by Sample Type")+xlab("Genus")

## Dust vs Soil
# X Axis Breaks and Labels
lbls = paste0(as.character(c(seq(0.05, 0, -0.01), seq(0.01, 0.05, 0.01)))) # labels
brks=seq(-0.05,0.05,0.01)
g_shared_typ2<-subset(g_shared_typ, Sample_Type!="Seawater")
g_shared_typ2$Counts2 <- ifelse(g_shared_typ2$Sample_Type == "Soil", -1*g_shared_typ2$Counts, g_shared_typ2$Counts)
g_shared_typ2<-g_shared_typ2[order(-g_shared_typ2$Counts2,g_shared_typ2$Genus),]
g_shared_typ2$GenSamp<-interaction(g_shared_typ2$Genus,g_shared_typ2$Sample_Type)
g_shared_typ2$GenSamp<-factor(g_shared_typ2$GenSamp, levels=g_shared_typ2$GenSamp)
class(g_shared_typ2$GenSamp)
g_shared_typ2$Genus<-factor(g_shared_typ2$Genus, levels=unique(g_shared_typ2$Genus[sort(g_shared_typ2$GenSamp)]))

share2<-ggplot(g_shared_typ2, aes(x = Genus, y = -Counts2, fill = Sample_Type)) +
  geom_bar(data = subset(g_shared_typ2[g_shared_typ2$Counts2<=-0.0005,], Sample_Type == "Soil"), stat = "identity") +
  geom_bar(data = subset(g_shared_typ2[g_shared_typ2$Counts2>=0.0005,], Sample_Type == "Dust"), stat = "identity") +
  coord_flip()+scale_y_continuous(labels = lbls,breaks=brks)+theme_classic()+scale_fill_manual(name ="Sample Type", values=unique(g_shared_typ2$Sample_Color[order(g_shared_typ2$Sample_Type)]))+ylab("Relative Abundance")+
  theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+labs(title="Microbial Genera by Sample Type", subtitle="Includes Taxa with a Relative Abundance of at least 0.05%")

ggsave(share2,filename = "figures_copy/Env_Seqs_All_Figs/16S_shared_Genera_dust.v.soil_population.pyramid.png", width=12, height=10, dpi=600)


ggplot(g_shared_typ2, aes(x = Genus, y = Counts2, fill = Sample_Type)) +
  geom_bar(data = subset(g_shared_typ2, Sample_Type == "Dust"), stat = "identity") +
  geom_bar(data = subset(g_shared_typ2, Sample_Type == "Soil"), stat = "identity") +
  coord_flip()+scale_y_continuous(labels = lbls,breaks=brks)+theme_classic()+scale_fill_manual(name ="Sample Type", values=unique(g_shared_typ2$Sample_Color[order(g_shared_typ2$Sample_Type)]))+ylab("Relative Abundance")+
  theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
  labs(title="Microbial Genera by Sample Type", subtitle="Includes Taxa with a Relative Abundance of at least 0.05%")

pp2a<-ggplot(g_shared_typ2[g_shared_typ2$Counts>=0.0005,], aes(x = reorder(Genus,Counts), fill = Sample_Type,y = ifelse(test = Sample_Type == "Dust",yes = Counts, no = -Counts))) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = abs, limits = max(g_shared_typ2$Counts) * c(-1,1)) +
  coord_flip()+scale_fill_manual(name ="Sample Type", values=unique(g_shared_typ2$Sample_Color[order(g_shared_typ2$Sample_Type)]))+ylab("Relative Abundance")+theme_classic()+
  theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
  labs(title="Microbial Genera by Sample Type",subtitle="Includes Only Shared Taxa with a Relative Abundance of at least 0.05%")+xlab("Genus")

ggsave(pp2a,filename = "figures_copy/Env_Seqs_All_Figs/16S_shared_Genera_dust.v.soil_population.pyramid.pretty.05percent.png", width=13, height=10, dpi=600)

pp3a<-ggplot(g_shared_typ2[g_shared_typ2$Counts>=0.010,], aes(x = reorder(Genus,Counts), fill = Sample_Type,y = ifelse(test = Sample_Type == "Dust",yes = Counts, no = -Counts))) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = abs, limits = max(g_shared_typ2$Counts) * c(-1,1)) +
  coord_flip()+scale_fill_manual(name ="Sample Type", values=unique(g_shared_typ2$Sample_Color[order(g_shared_typ2$Sample_Type)]))+ylab("Relative Abundance")+theme_classic()+
  theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
  labs(title="Microbial Genera by Sample Type",subtitle="Includes Only Shared Taxa with a Relative Abundance of at least 0.1%")+xlab("Genus")

ggsave(pp3a,filename = "figures_copy/Env_Seqs_All_Figs/16S_shared_Genera_dust.v.soil_population.pyramid.pretty.1percent.png", width=13, height=10, dpi=600)

pp4a<-ggplot(g_shared_typ2[g_shared_typ2$Counts>=0.010,], aes(x = reorder(Genus,Counts), fill = Sample_Type,y = ifelse(test = Sample_Type == "Dust",yes = Counts, no = -Counts))) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = abs, limits = max(g_shared_typ2$Counts) * c(-1,1)) +
  coord_flip()+scale_fill_manual(name ="Sample Type", values=unique(g_shared_typ2$Sample_Color[order(g_shared_typ2$Sample_Type)]))+ylab("Relative Abundance")+theme_classic()+
  theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
  labs(title="Microbial Genera by Sample Type",subtitle="Includes Only Shared Taxa with a Relative Abundance of at least 1%")+xlab("Genus")

ggsave(pp4a,filename = "figures_copy/Env_Seqs_All_Figs/16S_shared_Genera_dust.v.soil_population.pyramid.pretty_1percent_10.14.22.png", width=13, height=10, dpi=600)

## Lollipop chart

lg1a<-ggplot(g_shared_typ2[g_shared_typ2$Counts>=0.0005,], aes(x = reorder(Genus,Counts),
                                                               y = ifelse(test = Sample_Type == "Dust",yes = Counts, no = -Counts),color=g_shared_typ2$Sample_Type[g_shared_typ2$Counts>=0.0005])) +
  geom_point(stat='identity',size=3)  +
  geom_segment(aes(y = 0,
                   x = Genus,
                   yend = ifelse(test = Sample_Type == "Dust",yes = Counts, no = -Counts),
                   xend = Genus),color = "black") +
  coord_flip()+scale_color_manual(name ="Sample Type", values=unique(g_shared_typ2$Sample_Color[order(g_shared_typ2$Sample_Type)]))+ylab("Relative Abundance")+theme_classic()+theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+labs(title="Microbial Genera by Sample Type")+xlab("Genus")

ggsave(lg1a,filename = "figures_copy/Env_Seqs_All_Figs/16S_shared_Genera_dust.v.soil_lollipop_chart.png", width=13, height=10, dpi=600)

## Seawater vs Soil
# X Axis Breaks and Labels
lbls = paste0(as.character(c(seq(0.05, 0, -0.01), seq(0.01, 0.05, 0.01)))) # labels
brks=seq(-0.05,0.05,0.01)
g_shared_typ3<-subset(g_shared_typ, Sample_Type!="Dust")
g_shared_typ3$Counts2 <- ifelse(g_shared_typ3$Sample_Type == "Soil", -1*g_shared_typ3$Counts, g_shared_typ3$Counts)
g_shared_typ3<-g_shared_typ3[order(-g_shared_typ3$Counts2,g_shared_typ3$Genus),]
g_shared_typ3$GenSamp<-interaction(g_shared_typ3$Genus,g_shared_typ3$Sample_Type)
g_shared_typ3$GenSamp<-factor(g_shared_typ3$GenSamp, levels=g_shared_typ3$GenSamp)
class(g_shared_typ3$GenSamp)
g_shared_typ3$Genus<-factor(g_shared_typ3$Genus, levels=unique(g_shared_typ3$Genus[sort(g_shared_typ3$GenSamp)]))

share3<-ggplot(g_shared_typ3, aes(x = Genus, y = -Counts2, fill = Sample_Type)) +
  geom_bar(data = subset(g_shared_typ3[g_shared_typ3$Counts2<=-0.0005,], Sample_Type == "Soil"), stat = "identity") +
  geom_bar(data = subset(g_shared_typ3[g_shared_typ3$Counts2>=0.0005,], Sample_Type == "Seawater"), stat = "identity") +
  coord_flip()+scale_y_continuous(labels = lbls,breaks=brks)+theme_classic()+scale_fill_manual(name ="Sample Type", values=unique(g_shared_typ3$Sample_Color[order(g_shared_typ3$Sample_Type)]))+ylab("Relative Abundance")+
  theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+labs(title="Microbial Genera by Sample Type", subtitle="Includes Taxa with a Relative Abundance of at least 0.05%")

ggsave(share3,filename = "figures_copy/Env_Seqs_All_Figs/16S_shared_Genera_sw.v.soil_population.pyramid.png", width=12, height=10, dpi=600)


ggplot(g_shared_typ3, aes(x = Genus, y = Counts2, fill = Sample_Type)) +
  geom_bar(data = subset(g_shared_typ3, Sample_Type == "Seawater"), stat = "identity") +
  geom_bar(data = subset(g_shared_typ3, Sample_Type == "Soil"), stat = "identity") +
  coord_flip()+scale_y_continuous(labels = lbls,breaks=brks)+theme_classic()+scale_fill_manual(name ="Sample Type", values=unique(g_shared_typ3$Sample_Color[order(g_shared_typ3$Sample_Type)]))+ylab("Relative Abundance")+
  theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
  labs(title="Microbial Genera by Sample Type", subtitle="Includes Taxa with a Relative Abundance of at least 0.05%")

pp2b<-ggplot(g_shared_typ3[g_shared_typ3$Counts>=0.0005,], aes(x = reorder(Genus,Counts), fill = Sample_Type,y = ifelse(test = Sample_Type == "Seawater",yes = Counts, no = -Counts))) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = abs, limits = max(g_shared_typ3$Counts) * c(-1,1)) +
  coord_flip()+scale_fill_manual(name ="Sample Type", values=unique(g_shared_typ3$Sample_Color[order(g_shared_typ3$Sample_Type)]))+ylab("Relative Abundance")+theme_classic()+
  theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
  labs(title="Microbial Genera by Sample Type",subtitle="Includes Only Shared Taxa with a Relative Abundance of at least 0.05%")+xlab("Genus")

ggsave(pp2b,filename = "figures_copy/Env_Seqs_All_Figs/16S_shared_Genera_sw.v.soil_population.pyramid.pretty.05percent.png", width=13, height=10, dpi=600)

pp3b<-ggplot(g_shared_typ3[g_shared_typ3$Counts>=0.0010,], aes(x = reorder(Genus,Counts), fill = Sample_Type,y = ifelse(test = Sample_Type == "Seawater",yes = Counts, no = -Counts))) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = abs, limits = max(g_shared_typ3$Counts) * c(-1,1)) +
  coord_flip()+scale_fill_manual(name ="Sample Type", values=unique(g_shared_typ3$Sample_Color[order(g_shared_typ3$Sample_Type)]))+ylab("Relative Abundance")+theme_classic()+
  theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
  labs(title="Microbial Genera by Sample Type",subtitle="Includes Only Shared Taxa with a Relative Abundance of at least 0.1%")+xlab("Genus")

ggsave(pp3b,filename = "figures_copy/Env_Seqs_All_Figs/16S_shared_Genera_sw.v.soil_population.pyramid.pretty.1percent.png", width=13, height=10, dpi=600)

pp4b<-ggplot(g_shared_typ3[g_shared_typ3$Counts>=0.005,], aes(x = reorder(Genus,Counts), fill = Sample_Type,y = ifelse(test = Sample_Type == "Seawater",yes = Counts, no = -Counts))) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = abs, limits = max(g_shared_typ3$Counts) * c(-1,1)) +
  coord_flip()+scale_fill_manual(name ="Sample Type", values=unique(g_shared_typ3$Sample_Color[order(g_shared_typ3$Sample_Type)]))+ylab("Relative Abundance")+theme_classic()+
  theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
  labs(title="Microbial Genera by Sample Type",subtitle="Includes Only Shared Taxa with a Relative Abundance of at least 0.5%")+xlab("Genus")

ggsave(pp4b,filename = "figures_copy/Env_Seqs_All_Figs/16S_shared_Genera_sw.v.soil_population.pyramid.pretty.1percent.png", width=13, height=10, dpi=600)

pp5b<-ggplot(g_shared_typ3[g_shared_typ3$Counts>=0.010,], aes(x = reorder(Genus,Counts), fill = Sample_Type,y = ifelse(test = Sample_Type == "Seawater",yes = Counts, no = -Counts))) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = abs, limits = max(g_shared_typ3$Counts) * c(-1,1)) +
  coord_flip()+scale_fill_manual(name ="Sample Type", values=unique(g_shared_typ3$Sample_Color[order(g_shared_typ3$Sample_Type)]))+ylab("Relative Abundance")+theme_classic()+
  theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
  labs(title="Microbial Genera by Sample Type",subtitle="Includes Only Shared Taxa with a Relative Abundance of at least 1%")+xlab("Genus")

ggsave(pp5b,filename = "figures_copy/Env_Seqs_All_Figs/16S_shared_Genera_sw.v.soil_population.pyramid.pretty_1percent.png", width=13, height=10, dpi=600)

## Lollipop chart

lg1b<-ggplot(g_shared_typ3[g_shared_typ3$Counts>=0.0005,], aes(x = reorder(Genus,Counts),
                                                               y = ifelse(test = Sample_Type == "Seawater",yes = Counts, no = -Counts),color=g_shared_typ3$Sample_Type[g_shared_typ3$Counts>=0.0005])) +
  geom_point(stat='identity',size=3)  +
  geom_segment(aes(y = 0,
                   x = Genus,
                   yend = ifelse(test = Sample_Type == "Seawater",yes = Counts, no = -Counts),
                   xend = Genus),color = "black") +
  coord_flip()+scale_color_manual(name ="Sample Type", values=unique(g_shared_typ3$Sample_Color[order(g_shared_typ3$Sample_Type)]))+ylab("Relative Abundance")+theme_classic()+theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+labs(title="Microbial Genera by Sample Type")+xlab("Genus")

ggsave(lg1a,filename = "figures_copy/Env_Seqs_All_Figs/16S_shared_Genera_sw.v.soil_lollipop_chart.png", width=13, height=10, dpi=600)

typ_meta<-unique(data.frame("Sample_Type"=metadata$Sample_Type, "Sample_Color"=metadata$Sample_Color))
g_typ_meta<-merge(typ_meta,b.g.typ_m, by="Sample_Type")
g_typ_meta<-subset(g_typ_meta, Genus!="Unknown")
g_typ_meta<-subset(g_typ_meta, Counts!=0)

ts1<-ggplot(g_typ_meta[g_typ_meta$Counts>=0.0005,], aes(Genus,-Counts)) +
  geom_jitter(aes(color=factor(Sample_Type)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Sample Type", values=c(unique(g_typ_meta$Sample_Color[order(g_typ_meta$Sample_Type)])),c("Dust","Seawater")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=0.5,angle=90),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Bacteria/Archaea & Sample Type")

ggsave(ts1,filename = "figures_copy/Env_Seqs_All_Figs/16S_Gen.RA_typ.mat.png", width=20, height=10, dpi=600)

#ts1a<-ggplot(g_typ_meta, aes(Genus, Counts, label=Genus)) +
#  geom_jitter(aes(color=ifelse(Counts>0.01,factor(Sample_Type),"grey")), size=2, width=0.15, height=0) +
#  scale_color_manual(name ="Sample Type", values=c(unique(g_typ_meta$Sample_Color[order(g_typ_meta$Sample_Type)]),"grey"),c("Dust","Seawater","RA < 1%")) +
#  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
#  labs(x="Microbial Genera", y="Relative Abundance", title="Microbial Genera & Sample Type")

#ggsave(ts1a,filename = "figures_copy/Env_Seqs_All_Figs/16S_Phyla.RA_typ.mat2.png", width=15, height=10, dpi=600)

g.typ.05p<-na.omit(subset(g_typ_meta, Counts>=0.005))

g.t.05<-ggplot(g.typ.05p, aes(Genus, Counts)) +
  geom_jitter(aes(color=factor(Sample_Type)), size=2, width=0.15, height=0) + geom_boxplot(fill=NA, outlier.color=NA) +scale_color_manual(name ="Sample Type", values=c(unique(g.typ.05p$Sample_Color[order(g.typ.05p$Sample_Type)])),c("Seawater","Soil", "Dust")) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Bacteria/Archaea & Sample Type")

ggsave(g.t.05,filename = "figures_copy/Env_Seqs_All_Figs/16S_Gen.0.5percRA_typ.mat.png", width=15, height=10, dpi=600)

g.t.05a<-ggplot(g.typ.05p, aes(Genus, Counts)) +
  geom_jitter(aes(color=ifelse(Counts>0.01,factor(Sample_Type),"grey")), size=2, width=0.15, height=0) + geom_boxplot(fill=NA, outlier.color=NA) +scale_color_manual(name ="Sample Type", values=c(unique(g.typ.05p$Sample_Color[order(g.typ.05p$Sample_Type)]),"grey"),c("Seawater","Soil", "Dust","<1% RA")) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Bacteria/Archaea & Sample Type")

ggsave(g.t.05a,filename = "figures_copy/Env_Seqs_All_Figs/16S_Gen.0.5percRA_typ.mat.v2.png", width=15, height=10, dpi=600)

g.typ.1p<-subset(g_typ_meta, Counts>=0.01)
g.typ.1p<-subset(g.typ.1p, Genus!="Unknown")

g2<-ggplot(g.typ.1p, aes(Genus, Counts)) +
  geom_jitter(aes(color=factor(Sample_Type)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Sample Type", values=unique(g.typ.1p$Sample_Color[order(g.typ.1p$Sample_Type)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=0.5),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Microbial Genera by Sample Type", subtitle="Includes Taxa with a Relative Abundance of at least 1%")+coord_flip()

ggsave(g2,filename = "figures_copy/Env_Seqs_All_Figs/16S_Genera.RA_1percent_typ.mat_v1.png", width=12, height=10, dpi=600)

b.gen_RA.st<-ggplot(g.typ.1p, aes(x=Sample_Type, y=Counts, fill=Genus))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genera by Sample Type", x="Sample Type", y="Relative Abundance", fill="Genus", subtitle="Includes Taxa with a Relative Abundance of at least 1%")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=0.5),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+guides(fill=guide_legend(ncol=1))

ggsave(b.gen_RA.st,filename = "figures_copy/Env_Seqs_All_Figs/bacterial_genera_1percent_RA_by_SampleType.png", width=12, height=10, dpi=600)


head(g_typ_meta)
g_typ_meta.no0<-subset(g_typ_meta, Counts!=0)

tg.h1<-ggplot(g_typ_meta.no0, aes(Sample_Type, Genus, fill= Counts)) +geom_tile()+scale_fill_gradient2(low="blue", mid="white",high="red",midpoint=.025)+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.y=element_text(margin = margin(0,0)),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +labs(x="Sample Type", y="Microbial Genera", title="Microbial Genera & Sample Type",fill="Relative Abundance")+theme_classic()+scale_x_discrete(expand = c(0,0))

ggsave(tp.h1,filename = "figures_copy/Env_Seqs_All_Figs/Heatmap_16S_Phyla.RA_typ.mat.png", width=12, height=10, dpi=600)


head(g_typ_meta)
length(g_typ_meta$Sample_Type)

d.gen<-subset(g_typ_meta, Sample_Type=="Dust")
sw.gen<-subset(g_typ_meta, Sample_Type=="Seawater")

colnames(d.gen)[which(names(d.gen) == "Counts")] <- "D.RA"
d.gen$D.RA<-as.numeric(d.gen$D.RA)
#d.gen<-subset(d.gen, D.RA>0.0001)

colnames(sw.gen)[which(names(sw.gen) == "Counts")] <- "SW.RA"
sw.gen$SW.RA<-as.numeric(sw.gen$SW.RA)
#sw.gen<-subset(sw.gen, SW.RA>0.0001)

head(d.gen)
head(sw.gen)

tgen.comp<-merge(d.gen, sw.gen, by="Genus")
colnames(tgen.comp)[which(names(tgen.comp) == "Sample_Type.x")] <- "D_Samp"
colnames(tgen.comp)[which(names(tgen.comp) == "Sample_Type.y")] <- "SW_Samp"
colnames(tgen.comp)[which(names(tgen.comp) == "Sample_Color.x")] <- "D_color"
colnames(tgen.comp)[which(names(tgen.comp) == "Sample_Color.y")] <- "SW_color"
write.csv(tgen.comp,"16S_Genera_SampleType_Shared_RA.Separated.csv",row.names=FALSE)

tgen.comp$Genus[tgen.comp$D.RA==max(tgen.comp$D.RA)] # most relatively abundant genus in dust
tgen.comp$Genus[tgen.comp$SW.RA==max(tgen.comp$SW.RA)] # most relatively abundant genus in seawater
tgen.comp$Gen2 <- ifelse(tgen.comp$SW.RA>0.01 | tgen.comp$D.RA>0.01, tgen.comp$Genus, "Other")
unique(tgen.comp$Gen2) # all genera names
length(unique(tgen.comp$Gen2)) # how many colors do we need

gencol = melt(c(Bacillus="darkgreen",Blastococcus="orange",Geodermatophilus="green3", Halomonas="darkslategray4", Marinobacter="blue", Other="black", Paracoccus="firebrick1", Roseovarius="deeppink3", Salinicoccus="cornflowerblue", Sphingomonas="purple", Spiroplasma="deepskyblue1", Truepera="darkgoldenrod2"))
head(gencol)
gencol$Gen2<-rownames(gencol)
gencol

tgen.comp<-merge(tgen.comp, gencol, by="Gen2")
head(tgen.comp)
tgen.comp$Gen_Col<-as.character(tgen.comp$Gen_Col)

tgen.comp$SumRA<-tgen.comp$SW.RA+tgen.comp$D.RA
tgen.comp<-tgen.comp[order(-tgen.comp$SumRA),]
tgen.comp$SW.RA <- -1*tgen.comp$SW.RA
tgen.comp$Genus<-factor(tgen.comp$Genus, levels=tgen.comp$Genus)
class(tgen.comp$Genus)

test1<-ggplot(tgen.comp, aes(SW.RA, D.RA, color=Gen2))+geom_point(size=2.5) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+scale_colour_manual(values=unique(tgen.comp$Gen_Col[order(tgen.comp$Gen2)]))+
  labs(x="Seawater Relative Abundance", y="Dust Relative Abundance", title="Microbial Genera & Sample Type",subtitle="Points labeled 'Other' include Genera with a Relative Abundance of < 1%")+guides(shape = guide_legend(override.aes = list(size = 3)),col=guide_legend(title="Genus"))

ggsave(test1,filename = "figures_copy/Env_Seqs_All_Figs/16S_Gen.RA_sample.type_scatterplot.png", width=12, height=10, dpi=600)

test2<-ggplot(tgen.comp, aes(SW.RA, D.RA, color=Genus))  +geom_point(aes(color=ifelse(SW.RA>0.01 | D.RA>0.01,Genus,"black")),size=2) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Seawater Relative Abundance", y="Dust Relative Abundance", title="Microbial Genera & Sample Type")

ggsave(test2,filename = "figures_copy/Env_Seqs_All_Figs/16S_Gen.RA_type_linear1_colortest.png", width=12, height=10, dpi=600)

ggplot(tgen.comp, aes(SW.RA, D.RA, color=ifelse(SW.RA>0.01 | D.RA>0.01,Genus,"black")))  +geom_point() + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Seawater Relative Abundance", y="Dust Relative Abundance", title="Microbial Genera & Sample Type")

tgc1<-ggplot(tgen.comp, aes(SW.RA, D.RA, color=Genus)) +geom_point(aes(color=D.RA>0.01 | SW.RA>0.01))+ theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Seawater Relative Abundance", y="Dust Relative Abundance", title="Microbial Genera & Sample Type", subtitle="Only shows shared taxa w/ RA > 0.0001")

ggsave(tgc1,filename = "figures_copy/Env_Seqs_All_Figs/16S_Gen.RA_type_linear1.png", width=12, height=10, dpi=600)
## Date of note: 11/1/21 vvv
# 16S_Gen.RA_type_linear1 - cutoff is 0.0001
# 16S_Gen.RA_type_linear2 - cutoff is 0.00001
# 16S_Gen.RA_type_linear 3 - cutoff is 0.000001
# ** only 2 genera shared at 0.001 cutoff: Halomonas, Truepera

ggplot(g_shared_typ, aes(x = Genus, y = Counts2, fill = Sample_Type)) +
  geom_bar(data = subset(g_shared_typ, Sample_Type == "Dust"), stat = "identity") +
  geom_bar(data = subset(g_shared_typ, Sample_Type == "Seawater"), stat = "identity") +
  coord_flip()+scale_y_continuous(labels = lbls,breaks=brks)+theme_classic()+scale_fill_manual(name ="Sample Type", values=unique(g_shared_typ$Sample_Color[order(g_shared_typ$Sample_Type)]))+ylab("Relative Abundance")+theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+labs(title="Microbial Genera by Sample Type", subtitle="Includes Taxa with a Relative Abundance of at least 0.05%")



save.image("data/Env_Seqs_All/env.seq_analysis.Rdata")


#### Environmental Drivers - CCA/RDA ###

## Run DCA to determine if you should use: CCA or RDA
## *** ONLY USE DCA To determine if you should use a CCA or an RDA

## The length of first DCA axis:
## > 4 indicates heterogeneous dataset on which unimodal methods should be used (CCA),
##  < 3 indicates homogeneous dataset for which linear methods are suitable (RDA)
## between 3 and 4 both linear and unimodal methods are OK.

### all comp data dca
its1.dca = decorana(its1_RA_table)
plot(its1.dca)
summary (its1.dca) #DCA1 axis length = 3.0669; can use RDA or CCA

bac.dca = decorana(bac_RA_table)
plot(bac.dca)
summary (bac.dca) #DCA1 axis length = 5.0610; use CCA



# Version Information
sessionInfo()
