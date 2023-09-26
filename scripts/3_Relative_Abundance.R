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
  library(ape)
  library(plyr)
  library(dplyr)
  library(viridis)
  library(readxl)
  library(metagenomeSeq)
  library(AER)
  library(DESeq2)
  library(dplyr)
  library(magrittr)
  library(MASS)
  library(plotly)
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
  library(pairwiseAdonis)
  library(ggfortify)
})

#### Load Global Env to Import Count/ASV Tables ####
load("data/SSeawater_Data_Ready.Rdata") # save global env to Rdata file
#load("data/SSeawater_RelAb_Data.Rdata")

bac.dat.all[1:4,]
bac.dat.all$Depth_m
head(meta_scaled)

## DO NOT RUN THIS LINE, THIS IS YOUR COLOR REFERENCE!!!!
#(August.2021="#ef781c",December.2021="#03045e",April.2022="#059c3f")

# Create numeric version of depth variable for later correlations
meta_scaled$Depth.num<-as.numeric(as.character(meta_scaled$Depth_m))
unique(meta_scaled$Depth.num)

#### Phyla Relative Abundance ####

# use dcast to count up ASVs within each Phylum across all of the samples
b.phyla_counts <- as.data.frame(dcast(bac.dat.all, SampleID~Phylum, value.var="Count", fun.aggregate=sum)) ###
head(b.phyla_counts) # counts by phyla per sample
dim(b.phyla_counts)
rownames(b.phyla_counts)<-b.phyla_counts$SampleID

b.phyla_RelAb<-data.frame(decostand(b.phyla_counts[,-1], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(b.phyla_RelAb) # sanity check to make sure the transformation worked!

b.phyla_RelAb$SampleID<-rownames(b.phyla_RelAb)
head(b.phyla_RelAb)
#write.csv(b.phyla_RelAb,"16S_Phyla_Relative_Abundance.csv",row.names=TRUE) # good to save just in case

# melt down relativized data to merge with metadata
b.phyla_m<-melt(b.phyla_RelAb)

head(b.phyla_m)
colnames(b.phyla_m)[which(names(b.phyla_m) == "variable")] <- "Phylum"
colnames(b.phyla_m)[which(names(b.phyla_m) == "value")] <- "Count"
head(b.phyla_m) ## relative abundance based on sum of counts by phyla!

b.phyla_RA_meta<-merge(b.phyla_m,metadata, by="SampleID")
head(b.phyla_RA_meta) ## relative abundance based on sum of counts by phyla!

# rename SampleID so it's easy to interpret, then reorder SampleID based on sample date, then depth
b.phyla_RA_meta$PlotID = gsub("^SSW.","",b.phyla_RA_meta$SampleID)
b.phyla_RA_meta$PlotID = factor(b.phyla_RA_meta$PlotID, levels=unique(b.phyla_RA_meta$PlotID[order(b.phyla_RA_meta$SampDate,b.phyla_RA_meta$Depth_m)]), ordered=TRUE)

# Find most abundant phyla per sample per timepoint
aug.phy<-b.phyla_RA_meta[b.phyla_RA_meta$SampDate=="August.2021",]
aug.phy[which.max(aug.phy$Count),] # 48.86%, Actinobacteriota; 7m in Aug 2021
aug.phy[order(aug.phy$Count,decreasing=TRUE),]

dec.phy<-b.phyla_RA_meta[b.phyla_RA_meta$SampDate=="December.2021",]
dec.phy[which.max(dec.phy$Count),] #
dec.phy[order(dec.phy$Count,decreasing=TRUE),]

apr.phy<-b.phyla_RA_meta[b.phyla_RA_meta$SampDate=="April.2022",]
apr.phy[which.max(apr.phy$Count),] #
apr.phy[order(apr.phy$Count,decreasing=TRUE),]

# Barplot by SampleID
b.phy_RA<-ggplot(b.phyla_RA_meta, aes(x=PlotID, y=Count, fill=Phylum))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Phylum Relative Abundance", x="SampleID", y="Relative Abundance", fill="Phylum")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+guides(fill=guide_legend(ncol=1))+scale_y_continuous(expand = c(0,0),limits = c(0,1))

ggsave(b.phy_RA,filename = "figures/RelativeAbundance/Phylum/16S_Phyla.RA_barplot.png", width=12, height=10, dpi=600)

ggplot(b.phyla_RA_meta[b.phyla_RA_meta$Count>0.1,], aes(x=PlotID, y=Count, fill=Phylum))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Phylum Relative Abundance", x="SampleID", y="Relative Abundance", fill="Phylum")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+guides(fill=guide_legend(ncol=1))+scale_y_continuous(expand = c(0,0),limits = c(0,1))

b.phy.a21_RA<-ggplot(aug.phy, aes(x=PlotID, y=Count, fill=Phylum))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Phylum Relative Abundance - August 2021", x="SampleID", y="Relative Abundance", fill="Phylum")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+guides(fill=guide_legend(ncol=1))+scale_y_continuous(expand = c(0,0),limits = c(0,1))

ggsave(b.phy.a21_RA,filename = "figures/RelativeAbundance/Phylum/16S_August2021_Phyla.RA_barplot.png", width=12, height=10, dpi=600)

b.phy.d21_RA<-ggplot(dec.phy, aes(x=PlotID, y=Count, fill=Phylum))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Phylum Relative Abundance - December 2021", x="SampleID", y="Relative Abundance", fill="Phylum")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+guides(fill=guide_legend(ncol=1))+scale_y_continuous(expand = c(0,0),limits = c(0,1))

ggsave(b.phy.d21_RA,filename = "figures/RelativeAbundance/Phylum/16S_December2021_Phyla.RA_barplot.png", width=12, height=10, dpi=600)

b.phy.a22_RA<-ggplot(apr.phy, aes(x=PlotID, y=Count, fill=Phylum))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Phylum Relative Abundance - April 2022", x="SampleID", y="Relative Abundance", fill="Phylum")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+guides(fill=guide_legend(ncol=1))+scale_y_continuous(expand = c(0,0),limits = c(0,1))

ggsave(b.phy.a22_RA,filename = "figures/RelativeAbundance/Phylum/16S_April2022_Phyla.RA_barplot.png", width=12, height=10, dpi=600)

head(b.phyla_RA_meta)

# Heatmap by SampleID
p.h1<-ggplot(b.phyla_RA_meta, aes(PlotID, Phylum, fill= Count)) +geom_tile()+scale_fill_gradient2(low="blue3",mid="white",high="red",midpoint=0.3)+
  theme_classic()+theme(axis.title.x = element_text(size=13,vjust=-0.5),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Sample ID", y="Microbial Phyla", title="Microbial Phyla & Sample Type",fill="Relative Abundance")+scale_x_discrete(expand = c(0,0))

ggsave(p.h1,filename = "figures/RelativeAbundance/Phylum/16S_Phyla.RA_heatmap.png", width=12, height=10, dpi=600)

bac.dat.all[1:4,1:4]

# by Phylum + depth
bac.phy.dep <- as.data.frame(dcast(bac.dat.all,Depth_m~Phylum, value.var="Count", fun.aggregate=sum)) ###
head(bac.phy.dep) # counts by Phylum + sample depe
rownames(bac.phy.dep)<-bac.phy.dep$Depth_m

b.RA_phy.dep<-data.frame(decostand(bac.phy.dep[,-1], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(b.RA_phy.dep) # sanity check
b.RA_phy.dep$Depth_m<-rownames(b.RA_phy.dep) # Depth_m is now a character, not a factor!
head(b.RA_phy.dep)

#melt down relativized data to merge with metadata
b.phy.dep_m<-melt(b.RA_phy.dep, by="Depth_m")

head(b.phy.dep_m)
colnames(b.phy.dep_m)[which(names(b.phy.dep_m) == "variable")] <- "Phylum"
colnames(b.phy.dep_m)[which(names(b.phy.dep_m) == "value")] <- "Count"
head(b.phy.dep_m) ## relative abundance based on sum of counts by Phylum!

#dep_meta<-unique(data.frame("Depth_m"=metadata$Depth_m, "Sample_Color"=metadata$Sample_Color))
#p_dep_meta<-merge(dep_meta,b.phy.dep_m, by="Depth_m")

# Taxonomic Summary by Depth
tp1<-ggplot(b.phy.dep_m, aes(Phylum, Count)) +
  geom_jitter(aes(color=as.numeric(Depth_m)), size=2, width=0.15, height=0) +
  scale_colour_gradient2(low="red",mid="hotpink",high="blue",midpoint=5.25,guide = guide_colourbar(reverse = TRUE)) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Phyla", y="Relative Abundance", title="Microbial Phyla & Depth",color="Depth (m)")

ggsave(tp1,filename = "figures/RelativeAbundance/Phylum/SSW_16S_Phyla.RA_depth_taxasum.png", width=15, height=10, dpi=600)

tp1a<-ggplot(b.phy.dep_m[b.phy.dep_m$Count>0.05,], aes(Phylum, Count)) +
  geom_jitter(aes(color=as.numeric(Depth_m)), size=2, width=0.15, height=0) +
  scale_colour_gradient2(low="red",mid="hotpink",high="blue",midpoint=5.25,guide = guide_colourbar(reverse = TRUE)) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Phyla", y="Relative Abundance", title="Microbial Phyla & Depth",color="Depth (m)")

ggsave(tp1a,filename = "figures/RelativeAbundance/Phylum/SSW_16S_Phyla.RA_depth_taxasum_5percent.png", width=15, height=10, dpi=600)

# by Phylum + Sampling Date
bac.phy.date <- as.data.frame(dcast(bac.dat.all,SampDate~Phylum, value.var="Count", fun.aggregate=sum)) ###
head(bac.phy.date) # counts by Phylum + sample depe
rownames(bac.phy.date)<-bac.phy.date$SampDate

b.RA_phy.date<-data.frame(decostand(bac.phy.date[,-1], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(b.RA_phy.date) # sanity check
b.RA_phy.date$SampDate<-rownames(b.RA_phy.date)
head(b.RA_phy.date)

#melt down relativized data to merge with metadata
b.phy.date_m<-melt(b.RA_phy.date, by="SampDate")

head(b.phy.date_m)
colnames(b.phy.date_m)[which(names(b.phy.date_m) == "variable")] <- "Phylum"
colnames(b.phy.date_m)[which(names(b.phy.date_m) == "value")] <- "Count"
head(b.phy.date_m) ## relative abundance based on sum of counts by Phylum!

b.phy.date_m$SampDate<-factor(b.phy.date_m$SampDate, levels=c("August.2021","December.2021","April.2022"))

# Find most abundant phyla per time point
aug.phy.d<-b.phy.date_m[b.phy.date_m$SampDate=="August.2021",]
#aug.phy.d[which.max(aug.phy.d$Count),]
aug.phy.d[order(aug.phy.d$Count,decreasing=TRUE),]

dec.phy.d<-b.phy.date_m[b.phy.date_m$SampDate=="December.2021",]
#dec.phy.d[which.max(dec.phy.d$Count),] #
dec.phy.d[order(dec.phy.d$Count,decreasing=TRUE),]

apr.phy.d<-b.phy.date_m[b.phy.date_m$SampDate=="April.2022",]
#apr.phy.d[which.max(apr.phy.d$Count),] #
apr.phy.d[order(apr.phy.d$Count,decreasing=TRUE),]

# Barplot by Sample Date

psd0<-ggplot(b.phy.date_m, aes(x=SampDate, y=Count, fill=Phylum))+geom_bar(stat="identity",colour="black")+theme_classic()+
  labs(title = "Relative Abundance of Microbial Phyla", x="SampleID", y="Relative Abundance", fill="Phylum")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1))+ scale_x_discrete(labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022"))

ggsave(psd0,filename = "figures/RelativeAbundance/Phylum/16S_Phyla.RA_barplot_sampdate.png", width=12, height=10, dpi=600)

# Taxonomic Summary by Sample Date

#b.phy.date_m2<-merge(b.phy.date_m, metadata, by="SampDate")

colorset1 # remember which date goes with each color

psd1<-ggplot(b.phy.date_m, aes(Phylum, Count)) +
  geom_jitter(aes(color=factor(SampDate)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Sample Date", values=c("#ef781c","#03045e","#059c3f"), labels=c("August 2021","December 2021","April 2022")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=0.5,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Phyla", y="Relative Abundance", title="Microbial Phyla & Sample Date")

ggsave(psd1,filename = "figures/RelativeAbundance/Phylum/SSW_16S_Phyla.RA_date_taxasum.png", width=15, height=10, dpi=600)

psd1a<-ggplot(b.phy.date_m[b.phy.date_m$Count>0.05,], aes(Phylum, Count)) +
  geom_jitter(aes(color=factor(SampDate)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Sample Date", values=c("#ef781c","#03045e","#059c3f"), labels=c("August 2021","December 2021","April 2022")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=0.5,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Phyla", y="Relative Abundance", title="Microbial Phyla & Sample Date")

ggsave(psd1a,filename = "figures/RelativeAbundance/Phylum/SSW_16S_Phyla.RA_date_taxasum_5percent.png", width=15, height=10, dpi=600)

#### Class Relative Abundance ####

# use dcast to count up ASVs within each Class across all of the samples
b.class_counts <- as.data.frame(dcast(bac.dat.all, SampleID~Class, value.var="Count", fun.aggregate=sum)) ###
head(b.class_counts) # counts by class per sample
dim(b.class_counts)
rownames(b.class_counts)<-b.class_counts$SampleID
b.class_counts<-subset(b.class_counts, select=-c(Unknown))

b.class_RelAb<-data.frame(decostand(b.class_counts[,-1], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(b.class_RelAb) # sanity check to make sure the transformation worked!

b.class_RelAb$SampleID<-rownames(b.class_RelAb)
head(b.class_RelAb)
#write.csv(b.class_RelAb,"16S_class_Relative_Abundance.csv",row.names=TRUE) # good to save just in case

# melt down relativized data to merge with metadata
b.class_m<-melt(b.class_RelAb)

head(b.class_m)
colnames(b.class_m)[which(names(b.class_m) == "variable")] <- "Class"
colnames(b.class_m)[which(names(b.class_m) == "value")] <- "Count"
head(b.class_m) ## relative abundance based on sum of counts by class!

b.class_RA_meta<-merge(b.class_m,metadata, by="SampleID")
head(b.class_RA_meta) ## relative abundance based on sum of counts by class!

# rename SampleID so it's easy to interpret, then reorder SampleID based on sample date, then depth
b.class_RA_meta$PlotID = gsub("^SSW.","",b.class_RA_meta$SampleID)
b.class_RA_meta$PlotID = factor(b.class_RA_meta$PlotID, levels=unique(b.class_RA_meta$PlotID[order(b.class_RA_meta$SampDate,b.class_RA_meta$Depth_m)]), ordered=TRUE)

# Find most abundant class per sample per time point
aug.cls<-b.class_RA_meta[b.class_RA_meta$SampDate=="August.2021",]
aug.cls[which.max(aug.cls$Count),] # 47.08%, Actinobacteria; 7m in Aug 2021
aug.cls[order(aug.cls$Count,decreasing=TRUE),]

dec.cls<-b.class_RA_meta[b.class_RA_meta$SampDate=="December.2021",]
dec.cls[which.max(dec.cls$Count),] # 60.11%, Actinobacteria; m in Dec 2021
dec.cls[order(dec.cls$Count,decreasing=TRUE),]

apr.cls<-b.class_RA_meta[b.class_RA_meta$SampDate=="April.2022",]
apr.cls[which.max(apr.cls$Count),] #
apr.cls[order(apr.cls$Count,decreasing=TRUE),]

# Barplot by SampleID

b.cls_RA<-ggplot(b.class_RA_meta, aes(x=PlotID, y=Count, fill=Class))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Relative Abundance of Microbial Classes", x="SampleID", y="Relative Abundance", fill="Class")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+guides(fill=guide_legend(ncol=1))+scale_y_continuous(expand = c(0,0),limits = c(0,1))

ggsave(b.cls_RA,filename = "figures/RelativeAbundance/Class/16S_Class.RA_barplot.png", width=12, height=10, dpi=600)

ggplot(b.class_RA_meta[b.class_RA_meta$Count>0.1,], aes(x=PlotID, y=Count, fill=Class))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Relative Abundance of Microbial Classes", x="SampleID", y="Relative Abundance", fill="Class")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+guides(fill=guide_legend(ncol=1))+scale_y_continuous(expand = c(0,0),limits = c(0,1))

b.cls.a21_RA<-ggplot(aug.cls, aes(x=PlotID, y=Count, fill=Class))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Classes Relative Abundance - August 2021", x="SampleID", y="Relative Abundance", fill="Class")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+guides(fill=guide_legend(ncol=1))+scale_y_continuous(expand = c(0,0),limits = c(0,1))

ggsave(b.cls.a21_RA,filename = "figures/RelativeAbundance/Class/16S_August2021_Classes.RA_barplot.png", width=12, height=10, dpi=600)

b.cls.d21_RA<-ggplot(dec.cls, aes(x=PlotID, y=Count, fill=Class))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Classes Relative Abundance - December 2021", x="SampleID", y="Relative Abundance", fill="Class")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+guides(fill=guide_legend(ncol=1))+scale_y_continuous(expand = c(0,0),limits = c(0,1))

ggsave(b.cls.d21_RA,filename = "figures/RelativeAbundance/Class/16S_December2021_Classes.RA_barplot.png", width=12, height=10, dpi=600)

b.cls.a22_RA<-ggplot(apr.cls, aes(x=PlotID, y=Count, fill=Class))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Classes Relative Abundance - April 2022", x="SampleID", y="Relative Abundance", fill="Class")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+guides(fill=guide_legend(ncol=1))+scale_y_continuous(expand = c(0,0),limits = c(0,1))

ggsave(b.cls.a22_RA,filename = "figures/RelativeAbundance/Class/16S_April2022_Classes.RA_barplot.png", width=12, height=10, dpi=600)

head(b.class_RA_meta)

# Heatmap by PlotID
c.h1<-ggplot(b.class_RA_meta, aes(PlotID, Class, fill= Count)) +geom_tile()+scale_fill_gradient2(low="blue3",mid="white",high="red",midpoint=0.3)+
  theme_classic()+theme(axis.title.x = element_text(size=13,vjust=-0.5),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Sample ID", y="Microbial Class", title="Microbial Class & Sample Type",fill="Relative Abundance")+scale_x_discrete(expand = c(0,0))

ggsave(c.h1,filename = "figures/RelativeAbundance/Class/16S_class.RA_heatmap.png", width=12, height=10, dpi=600)

bac.dat.all[1:4,1:4]

# by Class + depth
bac.cls.dep <- as.data.frame(dcast(bac.dat.all,Depth_m~Class, value.var="Count", fun.aggregate=sum)) ###
head(bac.cls.dep) # counts by Class + sample depe
rownames(bac.cls.dep)<-bac.cls.dep$Depth_m
bac.cls.dep<-subset(bac.cls.dep, select=-c(Unknown))

b.RA_cls.dep<-data.frame(decostand(bac.cls.dep[,-1], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(b.RA_cls.dep) # sanity check
b.RA_cls.dep$Depth_m<-rownames(b.RA_cls.dep) # Depth_m is now a character, not a factor!
head(b.RA_cls.dep)

#melt down relativized data to merge with metadata
b.cls.dep_m<-melt(b.RA_cls.dep, by="Depth_m")

head(b.cls.dep_m)
colnames(b.cls.dep_m)[which(names(b.cls.dep_m) == "variable")] <- "Class"
colnames(b.cls.dep_m)[which(names(b.cls.dep_m) == "value")] <- "Count"
head(b.cls.dep_m) ## relative abundance based on sum of counts by Class!

b.cls.dep_m$Depth_m<-factor(b.cls.dep_m$Depth_m)
#dep_meta<-unique(data.frame("Depth_m"=metadata$Depth_m, "Sample_Color"=metadata$Sample_Color))
#p_dep_meta<-merge(dep_meta,b.cls.dep_m, by="Depth_m")

# Barplot by Depth

cd1<-ggplot(b.cls.dep_m, aes(x=Depth_m, y=Count, fill=Class))+geom_bar(stat="identity",colour="black")+theme_classic()+
labs(title = "Relative Abundance of Microbial Classes", x="Depth (m)", y="Relative Abundance", fill="Class")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=3))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1))
#+ scale_x_discrete(labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022"))

ggsave(cd1,filename = "figures/RelativeAbundance/Class/SSW_16S_class.RA_barplot_depth.png", width=12, height=10, dpi=600)

cd1a<-ggplot(b.cls.dep_m[b.cls.dep_m$Count>0.05,], aes(x=Depth_m, y=Count, fill=Class))+geom_bar(stat="identity",colour="black")+theme_classic()+
  labs(title = "Relative Abundance of Microbial Classes", x="SampleID", y="Relative Abundance", fill="Class",subtitle="Only Relative Abundance > 5%")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1))
#+ scale_x_discrete(labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022"))

ggsave(cd1a,filename = "figures/RelativeAbundance/Class/SSW_16S_class.RA_barplot_depth_5percent.png", width=12, height=10, dpi=600)

# Taxonomic Summary by Depth

cd2<-ggplot(b.cls.dep_m, aes(Class, Count)) +
  geom_jitter(aes(color=as.numeric(Depth_m)), size=2, width=0.15, height=0) +
  scale_colour_gradient2(high="blue3",low="red",midpoint=5.5,guide = guide_colourbar(reverse = TRUE)) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Class", y="Relative Abundance", title="Microbial Class & Depth",color="Depth (m)")

ggsave(cd2,filename = "figures/RelativeAbundance/Class/SSW_16S_class.RA_depth_taxasum.png", width=15, height=10, dpi=600)

cd2a<-ggplot(b.cls.dep_m[b.cls.dep_m$Count>0.05,], aes(Class, Count)) +
  geom_jitter(aes(color=as.numeric(Depth_m)), size=2, width=0.15, height=0) +
  scale_colour_gradient2(low="red",mid="hotpink",high="blue",midpoint=5.25,guide = guide_colourbar(reverse = TRUE)) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Class", y="Relative Abundance", title="Microbial Class & Depth",color="Depth (m)")

ggsave(cd2a,filename = "figures/RelativeAbundance/Class/SSW_16S_class.RA_depth_5percent_taxasum.png", width=15, height=10, dpi=600)

# by Class + Sampling Date
bac.cls.date <- as.data.frame(dcast(bac.dat.all,SampDate~Class, value.var="Count", fun.aggregate=sum)) ###
head(bac.cls.date) # counts by Class + sample depe
rownames(bac.cls.date)<-bac.cls.date$SampDate
bac.cls.date<-subset(bac.cls.date, select=-c(Unknown))

b.RA_cls.date<-data.frame(decostand(bac.cls.date[,-1], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(b.RA_cls.date) # sanity check
b.RA_cls.date$SampDate<-rownames(b.RA_cls.date)
head(b.RA_cls.date)

#melt down relativized data to merge with metadata
b.cls.date_m<-melt(b.RA_cls.date, by="SampDate")

head(b.cls.date_m)
colnames(b.cls.date_m)[which(names(b.cls.date_m) == "variable")] <- "Class"
colnames(b.cls.date_m)[which(names(b.cls.date_m) == "value")] <- "Count"
head(b.cls.date_m) ## relative abundance based on sum of counts by Class!

b.cls.date_m$SampDate<-factor(b.cls.date_m$SampDate, levels=c("August.2021","December.2021","April.2022"))

# Find most abundant class per time point
aug.cls.d<-b.cls.date_m[b.cls.date_m$SampDate=="August.2021",]
#aug.cls.d[which.max(aug.cls.d$Count),]
aug.cls.d[order(aug.cls.d$Count,decreasing=TRUE),]

dec.cls.d<-b.cls.date_m[b.cls.date_m$SampDate=="December.2021",]
#dec.cls.d[which.max(dec.cls.d$Count),] #
dec.cls.d[order(dec.cls.d$Count,decreasing=TRUE),]

apr.cls.d<-b.cls.date_m[b.cls.date_m$SampDate=="April.2022",]
#apr.cls.d[which.max(apr.cls.d$Count),] #
apr.cls.d[order(apr.cls.d$Count,decreasing=TRUE),]

# Barplot by Sample Date

csd1<-ggplot(b.cls.date_m, aes(x=SampDate, y=Count, fill=Class))+geom_bar(stat="identity",colour="black")+theme_classic()+
  labs(title = "Relative Abundance of Microbial Classes", x="SampleID", y="Relative Abundance", fill="Class")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=3))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1))+ scale_x_discrete(labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022"))

ggsave(csd1,filename = "figures/RelativeAbundance/Class/SSW_16S_class.RA_barplot_sampdate.png", width=12, height=10, dpi=600)

csd1<-ggplot(b.cls.date_m[b.cls.date_m$Count>0.05,], aes(x=SampDate, y=Count, fill=Class))+geom_bar(stat="identity",colour="black")+theme_classic()+
  labs(title = "Relative Abundance of Microbial Classes", x="SampleID", y="Relative Abundance", fill="Class",subtitle="Only Relative Abundance > 5%")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1))+ scale_x_discrete(labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022"))

ggsave(csd1,filename = "figures/RelativeAbundance/Class/SSW_16S_class.RA_barplot_sampdate_5percent.png", width=12, height=10, dpi=600)

# Taxonomic Summary by Sample Date

#b.cls.date_m2<-merge(b.cls.date_m, metadata, by="SampDate")

colorset1 # remember which date goes with each color

csd2<-ggplot(b.cls.date_m, aes(Class, Count)) +
  geom_jitter(aes(color=factor(SampDate)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Sample Date", values=c("#ef781c","#03045e","#059c3f"), labels=c("August 2021","December 2021","April 2022")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=0.5,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Class", y="Relative Abundance", title="Microbial Class & Sample Date")

ggsave(csd2,filename = "figures/RelativeAbundance/Class/SSW_16S_class_RA_date_taxasum.png", width=15, height=10, dpi=600)

csd3<-ggplot(unique(b.cls.date_m[b.cls.date_m$Count>0.1,]), aes(Class, Count)) +
  geom_jitter(aes(color=factor(SampDate)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Sample Date", values=c("#ef781c","#03045e","#059c3f"), labels=c("August 2021","December 2021","April 2022")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=0.5,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Class", y="Relative Abundance", title="Microbial Class & Sample Date")

ggsave(csd3,filename = "figures/RelativeAbundance/Class/SSW_16S_class.RA_date_taxasum_5percent.png", width=15, height=10, dpi=600)

#### Order Relative Abundance ####

# use dcast to count up ASVs within each Order across all of the samples
b.ord_counts <- as.data.frame(dcast(bac.dat.all, SampleID~Order, value.var="Count", fun.aggregate=sum)) ###
head(b.ord_counts) # counts by class per sample
dim(b.ord_counts)
rownames(b.ord_counts)<-b.ord_counts$SampleID
b.ord_counts<-subset(b.ord_counts, select=-c(Unknown))

b.ord_RelAb<-data.frame(decostand(b.ord_counts[,-1], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(b.ord_RelAb) # sanity check to make sure the transformation worked!

b.ord_RelAb$SampleID<-rownames(b.ord_RelAb)
head(b.ord_RelAb)
#write.csv(b.ord_RelAb,"16S_order_Relative_Abundance.csv",row.names=TRUE) # good to save just in case

# melt down relativized data to merge with metadata
b.ord_m<-melt(b.ord_RelAb)

head(b.ord_m)
colnames(b.ord_m)[which(names(b.ord_m) == "variable")] <- "Order"
colnames(b.ord_m)[which(names(b.ord_m) == "value")] <- "Count"
head(b.ord_m) ## relative abundance based on sum of counts by class!

b.ord_RA_meta<-merge(b.ord_m,metadata, by="SampleID")
head(b.ord_RA_meta) ## relative abundance based on sum of counts by class!

# rename SampleID so it's easy to interpret, then reorder SampleID based on sample date, then depth
b.ord_RA_meta$PlotID = gsub("^SSW.","",b.ord_RA_meta$SampleID)
b.ord_RA_meta$PlotID = factor(b.ord_RA_meta$PlotID, levels=unique(b.ord_RA_meta$PlotID[order(b.ord_RA_meta$SampDate,b.ord_RA_meta$Depth_m)]), ordered=TRUE)

# Find most abundant order per sample per time point
aug.ord<-b.ord_RA_meta[b.ord_RA_meta$SampDate=="August.2021",]
aug.ord[which.max(aug.ord$Count),] # 32.81%, Micrococcales; 9m in Aug 2021
aug.ord[order(aug.ord$Count,decreasing=TRUE),]

dec.ord<-b.ord_RA_meta[b.ord_RA_meta$SampDate=="December.2021",]
dec.ord[which.max(dec.ord$Count),] # 36.63%, Micrococcales, 9m
dec.ord[order(dec.ord$Count,decreasing=TRUE),]

apr.ord<-b.ord_RA_meta[b.ord_RA_meta$SampDate=="April.2022",]
apr.ord[which.max(apr.ord$Count),] # Nitriliruptorales
apr.ord[order(apr.ord$Count,decreasing=TRUE),]

# Barplot by SampleID

b.ord_RA<-ggplot(b.ord_RA_meta, aes(x=PlotID, y=Count, fill=Order))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Relative Abundance of Microbial Orders", x="SampleID", y="Relative Abundance", fill="Order")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+guides(fill=guide_legend(ncol=1))+scale_y_continuous(expand = c(0,0),limits = c(0,1))

ggsave(b.ord_RA,filename = "figures/RelativeAbundance/Order/16S_Order.RA_barplot.png", width=12, height=10, dpi=600)

ggplot(b.ord_RA_meta[b.ord_RA_meta$Count>0.1,], aes(x=PlotID, y=Count, fill=Order))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Relative Abundance of Microbial Orders", x="SampleID", y="Relative Abundance", fill="Order")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+guides(fill=guide_legend(ncol=1))+scale_y_continuous(expand = c(0,0),limits = c(0,1))

b.ord.a21_RA<-ggplot(aug.ord, aes(x=PlotID, y=Count, fill=Order))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Orders Relative Abundance - August 2021", x="SampleID", y="Relative Abundance", fill="Order")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+guides(fill=guide_legend(ncol=1))+scale_y_continuous(expand = c(0,0),limits = c(0,1))

ggsave(b.ord.a21_RA,filename = "figures/RelativeAbundance/Order/16S_August2021_Orders.RA_barplot.png", width=12, height=10, dpi=600)

b.ord.d21_RA<-ggplot(dec.ord, aes(x=PlotID, y=Count, fill=Order))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Orders Relative Abundance - December 2021", x="SampleID", y="Relative Abundance", fill="Order")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+guides(fill=guide_legend(ncol=1))+scale_y_continuous(expand = c(0,0),limits = c(0,1))

ggsave(b.ord.d21_RA,filename = "figures/RelativeAbundance/Order/16S_December2021_Orders.RA_barplot.png", width=12, height=10, dpi=600)

b.ord.a22_RA<-ggplot(apr.ord, aes(x=PlotID, y=Count, fill=Order))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Orders Relative Abundance - April 2022", x="SampleID", y="Relative Abundance", fill="Order")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+guides(fill=guide_legend(ncol=1))+scale_y_continuous(expand = c(0,0),limits = c(0,1))

ggsave(b.ord.a22_RA,filename = "figures/RelativeAbundance/Order/16S_April2022_Orders.RA_barplot.png", width=12, height=10, dpi=600)

head(b.ord_RA_meta)

# Heatmap by SampleID
o.h1<-ggplot(b.ord_RA_meta, aes(PlotID, Order, fill= Count)) +geom_tile()+scale_fill_gradient2(low="blue3",mid="white",high="red",midpoint=0.3)+
  theme_classic()+theme(axis.title.x = element_text(size=13,vjust=-0.5),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Sample ID", y="Microbial Order", title="Microbial Order & Sample Type",fill="Relative Abundance")+scale_x_discrete(expand = c(0,0))

ggsave(o.h1,filename = "figures/RelativeAbundance/Order/16S_order.RA_heatmap.png", width=12, height=10, dpi=600)

bac.dat.all[1:4,1:4]

# by Order + depth
bac.ord.dep <- as.data.frame(dcast(bac.dat.all,Depth_m~Order, value.var="Count", fun.aggregate=sum)) ###
head(bac.ord.dep) # counts by Order + sample depe
rownames(bac.ord.dep)<-bac.ord.dep$Depth_m
bac.ord.dep<-subset(bac.ord.dep, select=-c(Unknown))

b.RA_ord.dep<-data.frame(decostand(bac.ord.dep[,-1], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(b.RA_ord.dep) # sanity check
b.RA_ord.dep$Depth_m<-rownames(b.RA_ord.dep) # Depth_m is now a character, not a factor!
head(b.RA_ord.dep)

#melt down relativized data to merge with metadata
b.ord.dep_m<-melt(b.RA_ord.dep, by="Depth_m")

head(b.ord.dep_m)
colnames(b.ord.dep_m)[which(names(b.ord.dep_m) == "variable")] <- "Order"
colnames(b.ord.dep_m)[which(names(b.ord.dep_m) == "value")] <- "Count"
head(b.ord.dep_m) ## relative abundance based on sum of counts by Order!

b.ord.dep_m$Depth_m<-factor(b.ord.dep_m$Depth_m)
#dep_meta<-unique(data.frame("Depth_m"=metadata$Depth_m, "Sample_Color"=metadata$Sample_Color))
#p_dep_meta<-merge(dep_meta,b.ord.dep_m, by="Depth_m")

# Barplot by Depth

od1<-ggplot(b.ord.dep_m, aes(x=Depth_m, y=Count, fill=Order))+geom_bar(stat="identity",colour="black")+theme_classic()+
  labs(title = "Relative Abundance of Microbial Orders", x="Depth (m)", y="Relative Abundance", fill="Order")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=3))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1))
#+ scale_x_discrete(labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022"))

ggsave(od1,filename = "figures/RelativeAbundance/Order/SSW_16S_order.RA_barplot_depth.png", width=12, height=10, dpi=600)

od1a<-ggplot(b.ord.dep_m[b.ord.dep_m$Count>0.05,], aes(x=Depth_m, y=Count, fill=Order))+geom_bar(stat="identity",colour="black")+theme_classic()+
  labs(title = "Relative Abundance of Microbial Orders", x="SampleID", y="Relative Abundance", fill="Order",subtitle="Only Relative Abundance > 5%")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1))
#+ scale_x_discrete(labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022"))

ggsave(od1a,filename = "figures/RelativeAbundance/Order/SSW_16S_order.RA_barplot_depth_5percent.png", width=12, height=10, dpi=600)

# Taxonomic Summary by Depth

od2<-ggplot(b.ord.dep_m, aes(Order, Count)) +
  geom_jitter(aes(color=as.numeric(Depth_m)), size=2, width=0.15, height=0) +
  scale_colour_gradient2(high="blue3",low="red",midpoint=5.5,guide = guide_colourbar(reverse = TRUE)) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Order", y="Relative Abundance", title="Microbial Order & Depth",color="Depth (m)")

ggsave(od2,filename = "figures/RelativeAbundance/Order/SSW_16S_order.RA_depth_taxasum.png", width=15, height=10, dpi=600)

od2a<-ggplot(b.ord.dep_m[b.ord.dep_m$Count>0.05,], aes(Order, Count)) +
  geom_jitter(aes(color=as.numeric(Depth_m)), size=2, width=0.15, height=0) +
  scale_colour_gradient2(low="red",mid="hotpink",high="blue",midpoint=5.25,guide = guide_colourbar(reverse = TRUE)) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Order", y="Relative Abundance", title="Microbial Order & Depth",color="Depth (m)")

ggsave(od2a,filename = "figures/RelativeAbundance/Order/SSW_16S_order.RA_depth_5percent_taxasum.png", width=15, height=10, dpi=600)

# by Order + Sampling Date
bac.ord.date <- as.data.frame(dcast(bac.dat.all,SampDate~Order, value.var="Count", fun.aggregate=sum)) ###
head(bac.ord.date) # counts by Order + sample depe
rownames(bac.ord.date)<-bac.ord.date$SampDate
bac.ord.date<-subset(bac.ord.date, select=-c(Unknown))

b.RA_ord.date<-data.frame(decostand(bac.ord.date[,-1], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(b.RA_ord.date) # sanity check
b.RA_ord.date$SampDate<-rownames(b.RA_ord.date)
head(b.RA_ord.date)

#melt down relativized data to merge with metadata
b.ord.date_m<-melt(b.RA_ord.date, by="SampDate")

head(b.ord.date_m)
colnames(b.ord.date_m)[which(names(b.ord.date_m) == "variable")] <- "Order"
colnames(b.ord.date_m)[which(names(b.ord.date_m) == "value")] <- "Count"
head(b.ord.date_m) ## relative abundance based on sum of counts by Order!

b.ord.date_m$SampDate<-factor(b.ord.date_m$SampDate, levels=c("August.2021","December.2021","April.2022"))

# Find most abundant order per time point
aug.ord.d<-b.ord.date_m[b.ord.date_m$SampDate=="August.2021",]
#aug.ord.d[which.max(aug.ord.d$Count),]
head(aug.ord.d[order(aug.ord.d$Count,decreasing=TRUE),])

dec.ord.d<-b.ord.date_m[b.ord.date_m$SampDate=="December.2021",]
#dec.ord.d[which.max(dec.ord.d$Count),] #
head(dec.ord.d[order(dec.ord.d$Count,decreasing=TRUE),])

apr.ord.d<-b.ord.date_m[b.ord.date_m$SampDate=="April.2022",]
#apr.ord.d[which.max(apr.ord.d$Count),] #
head(apr.ord.d[order(apr.ord.d$Count,decreasing=TRUE),])

# Barplot by Sample Date

osd1<-ggplot(b.ord.date_m, aes(x=SampDate, y=Count, fill=Order))+geom_bar(stat="identity",colour="black")+theme_classic()+
  labs(title = "Relative Abundance of Microbial Orders", x="SampleID", y="Relative Abundance", fill="Order")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=3))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1))+ scale_x_discrete(labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022"))

ggsave(osd1,filename = "figures/RelativeAbundance/Order/SSW_16S_order.RA_barplot_sampdate.png", width=12, height=10, dpi=600)

osd1<-ggplot(b.ord.date_m[b.ord.date_m$Count>0.05,], aes(x=SampDate, y=Count, fill=Order))+geom_bar(stat="identity",colour="black")+theme_classic()+
  labs(title = "Relative Abundance of Microbial Orders", x="SampleID", y="Relative Abundance", fill="Order",subtitle="Only Relative Abundance > 5%")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1))+ scale_x_discrete(labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022"))

ggsave(osd1,filename = "figures/RelativeAbundance/Order/SSW_16S_order.RA_barplot_sampdate_5percent.png", width=12, height=10, dpi=600)

# Taxonomic Summary by Sample Date

#b.ord.date_m2<-merge(b.ord.date_m, metadata, by="SampDate")

colorset1 # remember which date goes with each color

osd2<-ggplot(b.ord.date_m, aes(Order, Count)) +
  geom_jitter(aes(color=factor(SampDate)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Sample Date", values=c("#ef781c","#03045e","#059c3f"), labels=c("August 2021","December 2021","April 2022")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=0.5,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Order", y="Relative Abundance", title="Microbial Order & Sample Date")

ggsave(osd2,filename = "figures/RelativeAbundance/Order/SSW_16S_order_RA_date_taxasum.png", width=15, height=10, dpi=600)

osd3<-ggplot(unique(b.ord.date_m[b.ord.date_m$Count>0.1,]), aes(Order, Count)) +
  geom_jitter(aes(color=factor(SampDate)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Sample Date", values=c("#ef781c","#03045e","#059c3f"), labels=c("August 2021","December 2021","April 2022")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=0.5,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Order", y="Relative Abundance", title="Microbial Order & Sample Date")

ggsave(osd3,filename = "figures/RelativeAbundance/Order/SSW_16S_order.RA_date_taxasum_5percent.png", width=15, height=10, dpi=600)


#### Family Relative Abundance ####

bac.dat.all.f1<-subset(bac.dat.all, bac.dat.all$Family!="Unknown") # drop unknown genera so they don't skew analyses
bac.dat.all.f<-subset(bac.dat.all.f1, bac.dat.all.f1$Family!="Unknown Family") # drop unknown genera so they don't skew analyses
"Unknown Family" %in% bac.dat.all.f

# use dcast to count up ASVs within each Family across all of the samples
b.fam_counts <- as.data.frame(dcast(bac.dat.all.f, SampleID~Family, value.var="Count", fun.aggregate=sum)) ###
head(b.fam_counts) # counts by fam per sample
dim(b.fam_counts)
rownames(b.fam_counts)<-b.fam_counts$SampleID
colnames(b.fam_counts)<-gsub(" ", ".",colnames(b.fam_counts))
#b.fam_counts<-subset(b.fam_counts, select=-c(Unknown, Unknown.Family))

b.fam_RelAb<-data.frame(decostand(b.fam_counts[,-1], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(b.fam_RelAb) # sanity check to make sure the transformation worked!

b.fam_RelAb$SampleID<-rownames(b.fam_RelAb)
head(b.fam_RelAb)
#write.csv(b.fam_RelAb,"16S_fam_Relative_Abundance.csv",row.names=TRUE) # good to save just in case

# melt down relativized data to merge with metadata
b.fam_m<-melt(b.fam_RelAb)

head(b.fam_m)
colnames(b.fam_m)[which(names(b.fam_m) == "variable")] <- "Family"
colnames(b.fam_m)[which(names(b.fam_m) == "value")] <- "Count"
head(b.fam_m) ## relative abundance based on sum of counts by fam!

b.fam_RA_meta<-merge(b.fam_m,metadata, by="SampleID")
head(b.fam_RA_meta) ## relative abundance based on sum of counts by fam!

# rename SampleID so it's easy to interpret, then reorder SampleID based on sample date, then depth
b.fam_RA_meta$PlotID = gsub("^SSW.","",b.fam_RA_meta$SampleID)
b.fam_RA_meta$PlotID = factor(b.fam_RA_meta$PlotID, levels=unique(b.fam_RA_meta$PlotID[order(b.fam_RA_meta$SampDate,b.fam_RA_meta$Depth_m)]), ordered=TRUE)

# Find most abundant family per sample per time point
aug.fam<-b.fam_RA_meta[b.fam_RA_meta$SampDate=="August.2021",]
aug.fam[which.max(aug.fam$Count),] # 32.81%, Micrococcales; 9m in Aug 2021
aug.fam[order(aug.fam$Count,decreasing=TRUE),]

dec.fam<-b.fam_RA_meta[b.fam_RA_meta$SampDate=="December.2021",]
dec.fam[which.max(dec.fam$Count),] # 36.63%, Micrococcales, 9m
dec.fam[order(dec.fam$Count,decreasing=TRUE),]

apr.fam<-b.fam_RA_meta[b.fam_RA_meta$SampDate=="April.2022",]
apr.fam[which.max(apr.fam$Count),] # Nitriliruptorales
apr.fam[order(apr.fam$Count,decreasing=TRUE),]

# Barplot by SampleID

b.fam_RA<-ggplot(b.fam_RA_meta, aes(x=PlotID, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Relative Abundance of Microbial Families", x="SampleID", y="Relative Abundance", fill="Family")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+guides(fill=guide_legend(ncol=1))+scale_y_continuous(expand = c(0,0),limits = c(0,1))

ggsave(b.fam_RA,filename = "figures/RelativeAbundance/Family/16S_fam.RA_barplot.png", width=12, height=10, dpi=600)

ggplot(b.fam_RA_meta[b.fam_RA_meta$Count>0.1,], aes(x=PlotID, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Relative Abundance of Microbial Families", x="SampleID", y="Relative Abundance", fill="Family")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+guides(fill=guide_legend(ncol=1))+scale_y_continuous(expand = c(0,0),limits = c(0,1))

b.fam.a21_RA<-ggplot(aug.fam, aes(x=PlotID, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Families Relative Abundance - August 2021", x="SampleID", y="Relative Abundance", fill="Family")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+guides(fill=guide_legend(ncol=1))+scale_y_continuous(expand = c(0,0),limits = c(0,1))

ggsave(b.fam.a21_RA,filename = "figures/RelativeAbundance/Family/16S_August2021_Families.RA_barplot.png", width=12, height=10, dpi=600)

b.fam.d21_RA<-ggplot(dec.fam, aes(x=PlotID, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Families Relative Abundance - December 2021", x="SampleID", y="Relative Abundance", fill="Family")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+guides(fill=guide_legend(ncol=1))+scale_y_continuous(expand = c(0,0),limits = c(0,1))

ggsave(b.fam.d21_RA,filename = "figures/RelativeAbundance/Family/16S_December2021_Families.RA_barplot.png", width=12, height=10, dpi=600)

b.fam.a22_RA<-ggplot(apr.fam, aes(x=PlotID, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Families Relative Abundance - April 2022", x="SampleID", y="Relative Abundance", fill="Family")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+guides(fill=guide_legend(ncol=1))+scale_y_continuous(expand = c(0,0),limits = c(0,1))

ggsave(b.fam.a22_RA,filename = "figures/RelativeAbundance/Family/16S_April2022_Families.RA_barplot.png", width=12, height=10, dpi=600)

head(b.fam_RA_meta)

# Heatmap by PlotID

f.h1<-ggplot(b.fam_RA_meta, aes(PlotID, Family, fill= Count)) +geom_tile()+scale_fill_gradient2(low="blue3",mid="white",high="red",midpoint=0.3)+
  theme_classic()+theme(axis.title.x = element_text(size=13,vjust=-0.5),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Sample ID", y="Microbial Family", title="Microbial Families & Sample Type",fill="Relative Abundance")+scale_x_discrete(expand = c(0,0))

ggsave(f.h1,filename = "figures/RelativeAbundance/Family/16S_fam.RA_heatmap.png", width=12, height=10, dpi=600)

bac.dat.all.f[1:4,1:4]

# # by Family + depth
# bac.fam.dep <- as.data.frame(dcast(bac.dat.all,Depth_m~Family, value.var="Count", fun.aggregate=sum)) ###
# head(bac.fam.dep) # counts by Family + sample depe
# rownames(bac.fam.dep)<-bac.fam.dep$Depth_m
# colnames(bac.fam.dep)<-gsub(" ", ".",colnames(bac.fam.dep))
# bac.fam.dep<-subset(bac.fam.dep, select=-c(Unknown, Unknown.Family))
#
# b.RA_fam.dep<-data.frame(decostand(bac.fam.dep[,-1], method="total", MARGIN=1, na.rm=TRUE))
# # relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
# rowSums(b.RA_fam.dep) # sanity check
# b.RA_fam.dep$Depth_m<-rownames(b.RA_fam.dep) # Depth_m is now a character, not a factor!
# head(b.RA_fam.dep)
#
# #melt down relativized data to merge with metadata
# b.fam.dep_m<-melt(b.RA_fam.dep, by="Depth_m")
#
# head(b.fam.dep_m)
# colnames(b.fam.dep_m)[which(names(b.fam.dep_m) == "variable")] <- "Family"
# colnames(b.fam.dep_m)[which(names(b.fam.dep_m) == "value")] <- "Count"
# head(b.fam.dep_m) ## relative abundance based on sum of counts by Family!
#
# #dep_meta<-unique(data.frame("Depth_m"=metadata$Depth_m, "Sample_Color"=metadata$Sample_Color))
# #p_dep_meta<-merge(dep_meta,b.fam.dep_m, by="Depth_m")
#
# # Barplot by Depth
#
# fd1<-ggplot(b.fam.dep_m, aes(x=Depth_m, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+theme_classic()+
#   labs(title = "Relative Abundance of Microbial Classes", x="Depth (m)", y="Relative Abundance", fill="Class")+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
#   guides(fill=guide_legend(ncol=3))+
#   scale_y_continuous(expand = c(0,0),limits = c(0,1))
# #+ scale_x_discrete(labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022"))
#
# ggsave(fd1,filename = "figures/RelativeAbundance/Family/SSW_16S_fam.RA_barplot_depth.png", width=12, height=10, dpi=600)
#
# fd1a<-ggplot(b.fam.dep_m[b.fam.dep_m$Count>0.05,], aes(x=Depth_m, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+theme_classic()+
#   labs(title = "Relative Abundance of Microbial Classes", x="SampleID", y="Relative Abundance", fill="Class",subtitle="Only Relative Abundance > 5%")+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
#   scale_y_continuous(expand = c(0,0),limits = c(0,1))
# #+ scale_x_discrete(labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022"))
#
# ggsave(fd1a,filename = "figures/RelativeAbundance/Family/SSW_16S_fam.RA_barplot_depth_5percent.png", width=12, height=10, dpi=600)
#
# # Taxonomic Summary by Depth
#
# fd2<-ggplot(b.fam.dep_m, aes(Family, Count)) +
#   geom_jitter(aes(color=as.numeric(Depth_m)), size=2, width=0.15, height=0) +
#   scale_colour_gradient2(high="blue3",low="red",midpoint=5.5,guide = guide_colourbar(reverse = TRUE)) +
#   geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
#   labs(x="Microbial Family", y="Relative Abundance", title="Microbial Families & Depth",color="Depth (m)")
#
# ggsave(fd2,filename = "figures/RelativeAbundance/Family/SSW_16S_fam.RA_depth_taxasum.png", width=15, height=10, dpi=600)
#
# fd2a<-ggplot(b.fam.dep_m[b.fam.dep_m$Count>0.05,], aes(Family, Count)) +
#   geom_jitter(aes(color=as.numeric(Depth_m)), size=2, width=0.15, height=0) +
#   scale_colour_gradient2(low="red",mid="hotpink",high="blue",midpoint=5.25,guide = guide_colourbar(reverse = TRUE)) +
#   geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
#   labs(x="Microbial Family", y="Relative Abundance", title="Microbial Families & Depth",color="Depth (m)")
#
# ggsave(fd2a,filename = "figures/RelativeAbundance/Family/SSW_16S_fam.RA_depth_5percent_taxasum.png", width=15, height=10, dpi=600)

# by Family + Sampling Date
bac.fam.date <- as.data.frame(dcast(bac.dat.all.f,SampDate~Family, value.var="Count", fun.aggregate=sum)) ###
head(bac.fam.date) # counts by Family + sample depe
rownames(bac.fam.date)<-bac.fam.date$SampDate
colnames(bac.fam.date)<-gsub(" ", ".",colnames(bac.fam.date))
#bac.fam.date<-subset(bac.fam.date, select=-c(Unknown, Unknown.Family))

b.RA_fam.date<-data.frame(decostand(bac.fam.date[,-1], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(b.RA_fam.date) # sanity check
b.RA_fam.date$SampDate<-rownames(b.RA_fam.date)
head(b.RA_fam.date)

#melt down relativized data to merge with metadata
b.fam.date_m<-melt(b.RA_fam.date, by="SampDate")

head(b.fam.date_m)
colnames(b.fam.date_m)[which(names(b.fam.date_m) == "variable")] <- "Family"
colnames(b.fam.date_m)[which(names(b.fam.date_m) == "value")] <- "Count"
head(b.fam.date_m) ## relative abundance based on sum of counts by Family!

b.fam.date_m$SampDate<-factor(b.fam.date_m$SampDate, levels=c("August.2021","December.2021","April.2022"))

# Find most abundant family per time point
aug.fam.d<-b.fam.date_m[b.fam.date_m$SampDate=="August.2021",]
#aug.fam.d[which.max(aug.fam.d$Count),]
head(aug.fam.d[order(aug.fam.d$Count,decreasing=TRUE),])

dec.fam.d<-b.fam.date_m[b.fam.date_m$SampDate=="December.2021",]
#dec.fam.d[which.max(dec.fam.d$Count),] #
head(dec.fam.d[order(dec.fam.d$Count,decreasing=TRUE),])

apr.fam.d<-b.fam.date_m[b.fam.date_m$SampDate=="April.2022",]
#apr.fam.d[which.max(apr.fam.d$Count),] #
head(apr.fam.d[order(apr.fam.d$Count,decreasing=TRUE),])

# Barplot by Sample Date

fsd1<-ggplot(b.fam.date_m, aes(x=SampDate, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+theme_classic()+
  labs(title = "Relative Abundance of Microbial Families", x="SampleID", y="Relative Abundance", fill="Family")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=3))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1))+ scale_x_discrete(labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022"))

ggsave(fsd1,filename = "figures/RelativeAbundance/Family/SSW_16S_fam.RA_barplot_sampdate.png", width=12, height=10, dpi=600)

fsd1a<-ggplot(b.fam.date_m[b.fam.date_m$Count>0.01,], aes(x=SampDate, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+theme_classic()+
  labs(title = "Relative Abundance of Microbial Families", x="SampleID", y="Relative Abundance", fill="Family",subtitle="Only Relative Abundance > 1%")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1))+ scale_x_discrete(labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022"))

ggsave(fsd1a,filename = "figures/RelativeAbundance/Family/SSW_16S_fam.RA_barplot_sampdate_1percent.png", width=12, height=10, dpi=600)

# Taxonomic Summary by Sample Date

#b.fam.date_m2<-merge(b.fam.date_m, metadata, by="SampDate")

colorset1 # remember which date goes with each color

fsd2<-ggplot(b.fam.date_m, aes(Family, Count)) +
  geom_jitter(aes(color=factor(SampDate)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Sample Date", values=c("#ef781c","#03045e","#059c3f"), labels=c("August 2021","December 2021","April 2022")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=0.5,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Family", y="Relative Abundance", title="Microbial Families & Sample Date")

ggsave(fsd2,filename = "figures/RelativeAbundance/Family/SSW_16S_fam.RA_date_taxasum.png", width=15, height=10, dpi=600)

fsd3<-ggplot(b.fam.date_m[b.fam.date_m$Count>0.1,], aes(Family, Count)) +
  geom_jitter(aes(color=factor(SampDate)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Sample Date", values=c("#ef781c","#03045e","#059c3f"), labels=c("August 2021","December 2021","April 2022")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=0.5,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Family", y="Relative Abundance", title="Microbial Families & Sample Date")

ggsave(fsd3,filename = "figures/RelativeAbundance/Family/SSW_16S_fam.RA_date_taxasum_5percent.png", width=15, height=10, dpi=600)

# by Family + Sampling Date + Depth

bac.fam.date.dep <- as.data.frame(dcast(bac.dat.all.f,SampDate+Depth_m~Family, value.var="Count", fun.aggregate=sum)) ###
bac.fam.date.dep[1:5,1:5] # counts by Genus + sample date & depth
rownames(bac.fam.date.dep)<-interaction(bac.fam.date.dep$SampDate,bac.fam.date.dep$Depth_m,sep="_")
bac.fam.date.dep[1:5,1:5]

b.RA_fam.date.dep<-data.frame(decostand(bac.fam.date.dep[,-c(1:2)], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(b.RA_fam.date.dep) # sanity check
b.RA_fam.date.dep$SampDate_Depth<-rownames(b.RA_fam.date.dep)
b.RA_fam.date.dep[1:5,ncol(b.RA_fam.date.dep):5-ncol(b.RA_fam.date.dep)] # first 5 rows, last 5 columns

#melt down relativized data to merge with metadata
b.fam.date.dep_m<-melt(b.RA_fam.date.dep, by="SampDate_Depth")

head(b.fam.date.dep_m)
colnames(b.fam.date.dep_m)[which(names(b.fam.date.dep_m) == "variable")] <- "Family"
colnames(b.fam.date.dep_m)[which(names(b.fam.date.dep_m) == "value")] <- "Count"
head(b.fam.date.dep_m) ## relative abundance based on sum of counts by Family!

# separate SampDate_Depth above by _, then recreate column in new df
b.fam.date.dep_m2<-as.data.frame(separate_wider_delim(data = b.fam.date.dep_m, col=SampDate_Depth, "_", names = c("SampDate", "Depth_m"))) # Separate SampDate & Depth column for Heatmap later
b.fam.date.dep_m2$SampDate_Depth<-interaction(b.fam.date.dep_m2$SampDate,b.fam.date.dep_m2$Depth_m)

b.fam.date.dep_m2$Depth_m<-factor(b.fam.date.dep_m2$Depth_m, levels=c("0","3","4","5","7","9","10","10.5"))
b.fam.date.dep_m2$SampDate<-factor(b.fam.date.dep_m2$SampDate,levels=c("August.2021","December.2021","April.2022"))
b.fam.date.dep_m2$SampDate_Depth<-gsub("(\\..*?)\\.","\\1-",b.fam.date.dep_m2$SampDate_Depth) # instructions on changing second period below
# \\. - first period, .* is any character after that, ? suppresses greedy matching of .*, () - what we want to keep aka \\1
# \\1 is the pattern we want to keep in (), and - is replacing second character we want to replace with -
# more info here: https://stackoverflow.com/questions/43077846/how-to-replace-second-or-more-occurrences-of-a-dot-from-a-column-name
b.fam.date.dep_m2$SampDate_Depth = factor(b.fam.date.dep_m2$SampDate_Depth, levels=unique(b.fam.date.dep_m2$SampDate_Depth[order(b.fam.date.dep_m2$SampDate,b.fam.date.dep_m2$Depth_m)]), ordered=TRUE)

# Heatmaps
max(b.fam.date.dep_m2$Count)
max(b.fam.date.dep_m2$Count)/2 # finding the midpoint

f.sd.d.h1<-ggplot(b.fam.date.dep_m2[b.fam.date.dep_m2$Count>0.01,], aes(SampDate_Depth, Family, fill= Count)) +geom_tile()+scale_fill_gradient2(low="orange",mid="white",high="purple",midpoint=0.18)+
  theme_classic()+theme(axis.title.x = element_text(size=13,vjust=-0.5),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Sampling Date - Depth (m)", y="Microbial Families", title="Microbial Families by Sample Date & Depth",subtitle="Includes taxa with Relative Abundance > 1%",fill="Relative Abundance")+scale_x_discrete(expand = c(0,0))
#+ scale_x_discrete(labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022"))

ggsave(f.sd.d.h1,filename = "figures/RelativeAbundance/Family/16S_Families.RA_heatmap_date_depth_1perc.png", width=20, height=15, dpi=600)

f.sd.d.h2<-ggplot(b.fam.date.dep_m2[b.fam.date.dep_m2$Count>0.05,], aes(SampDate_Depth, Family, fill= Count)) +geom_tile()+scale_fill_gradient2(low="orange",mid="white",high="purple",midpoint=0.18)+
  theme_classic()+theme(axis.title.x = element_text(size=13,vjust=-0.5),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Sampling Date - Depth (m)", y="Microbial Families", title="Microbial Families by Sample Date & Depth",subtitle="Includes taxa with Relative Abundance > 1%",fill="Relative Abundance")+scale_x_discrete(expand = c(0,0))
#+ scale_x_discrete(labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022"))

ggsave(f.sd.d.h2,filename = "figures/RelativeAbundance/Family/16S_Families.RA_heatmap_date_depth_5perc.png", width=20, height=15, dpi=600)

# Taxonomic Summaries
f.sd.d.hm.1<-ggplot(b.fam.date.dep_m2[b.fam.date.dep_m2$Count>0.01,], aes(Family, Count)) +
  geom_jitter(aes(color=as.numeric(as.character(Depth_m)),shape=SampDate), size=3, width=0.15, height=0) +
  scale_colour_gradient2(low="red",mid="hotpink",high="blue",midpoint=5.25,guide = guide_colourbar(reverse = TRUE)) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() + scale_shape_discrete(labels=c("August 2021","December 2021","April 2022"),name="Sample Date") +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Families", y="Relative Abundance", title="Microbial Families by Sample Date & Depth",subtitle="Includes taxa with Relative Abundance > 1%",color="Depth (m)", shape="Sample Date")

ggsave(f.sd.d.hm.1,filename = "figures/RelativeAbundance/Family/SSW_16S_Families.RA_date_depth_taxasum_1perc.png", width=18, height=10, dpi=600)
## ^ this figure includes the relative abundance of each organism by depth & date!!!

f.sd.d.hm.2<-ggplot(b.fam.date.dep_m2[b.fam.date.dep_m2$Count>0.025,], aes(Family, Count)) +
  geom_jitter(aes(color=as.numeric(as.character(Depth_m)),shape=SampDate), size=3, width=0.15, height=0) +
  scale_colour_gradient2(low="red",mid="hotpink",high="blue",midpoint=5.25,guide = guide_colourbar(reverse = TRUE)) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() + scale_shape_discrete(labels=c("August 2021","December 2021","April 2022"),name="Sample Date") +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Families", y="Relative Abundance", title="Microbial Families by Sample Date & Depth",subtitle="Includes taxa with Relative Abundance > 2.5%",color="Depth (m)", shape="Sample Date")

ggsave(f.sd.d.hm.2,filename = "figures/RelativeAbundance/Family/SSW_16S_Families.RA_date_depth_taxasum_2.5perc.png", width=18, height=10, dpi=600)

f.sd.d.hm.3<-ggplot(b.fam.date.dep_m2[b.fam.date.dep_m2$Count>0.05,], aes(Family, Count)) +
  geom_jitter(aes(color=as.numeric(as.character(Depth_m)),shape=SampDate), size=3, width=0.15, height=0) +
  scale_colour_gradient2(low="red",mid="hotpink",high="blue",midpoint=5.25,guide = guide_colourbar(reverse = TRUE)) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() + scale_shape_discrete(labels=c("August 2021","December 2021","April 2022"),name="Sample Date") +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Families", y="Relative Abundance", title="Microbial Families by Sample Date & Depth",subtitle="Includes taxa with Relative Abundance > 5%",color="Depth (m)", shape="Sample Date")

ggsave(f.sd.d.hm.3,filename = "figures/RelativeAbundance/Family/SSW_16S_Families.RA_date_depth_taxasum_5perc.png", width=15, height=10, dpi=600)

#### Genus Relative Abundance ####
# use dcast to count up ASVs within each Genus across all of the samples
head(bac.dat.all)
bac.dat.all.g<-subset(bac.dat.all, bac.dat.all$Genus!="Unknown") # drop unknown genera so they don't skew analyses
"Unknown" %in% bac.dat.all.g$Genus

b.genus_counts <- as.data.frame(dcast(bac.dat.all.g, SampleID~Genus, value.var="Count", fun.aggregate=sum)) ###
head(b.genus_counts) # counts by genus per sample
dim(b.genus_counts)
rownames(b.genus_counts)<-b.genus_counts$SampleID
b.genus_counts[1:4,1:4]

b.genus_RelAb<-data.frame(decostand(b.genus_counts[,-1], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(b.genus_RelAb) # sanity check to make sure the transformation worked!

b.genus_RelAb$SampleID<-rownames(b.genus_RelAb)
head(b.genus_RelAb)
#write.csv(b.genus_RelAb,"16S_Genera_Relative_Abundance.csv",row.names=TRUE) # good to save just in case

# melt down relativized data to merge with metadata
b.genus_m<-melt(b.genus_RelAb)

head(b.genus_m)
colnames(b.genus_m)[which(names(b.genus_m) == "variable")] <- "Genus"
colnames(b.genus_m)[which(names(b.genus_m) == "value")] <- "Count"
head(b.genus_m) ## relative abundance based on sum of counts by genus!
b.genus_m$Genus<-gsub("^X.","",b.genus_m$Genus) # get rid of leading X. in Genus names
b.genus_m$Genus<-gsub("\\.\\."," ",b.genus_m$Genus) # get rid of .. in species name --> . is regex
b.genus_m$Genus<-gsub("_"," ",b.genus_m$Genus) # get rid of .. in species name --> . is regex
head(b.genus_m) ## relative abundance based on sum of counts by genus!

b.genus_RA_meta<-merge(b.genus_m,metadata, by="SampleID")
head(b.genus_RA_meta) ## relative abundance based on sum of counts by genus!
max(b.genus_RA_meta$Count)

# rename SampleID so it's easy to interpret, then reorder SampleID based on sample date, then depth
b.genus_RA_meta$PlotID = gsub("^SSW.","",b.genus_RA_meta$SampleID)
b.genus_RA_meta$PlotID = factor(b.genus_RA_meta$PlotID, levels=unique(b.genus_RA_meta$PlotID[order(b.genus_RA_meta$SampDate,b.genus_RA_meta$Depth_m)]), ordered=TRUE)

# Find most abundant family per sample per time point
aug.gen<-b.genus_RA_meta[b.genus_RA_meta$SampDate=="August.2021",]
aug.gen[which.max(aug.gen$Count),] #
aug.gen[order(aug.gen$Count,decreasing=TRUE),]

dec.gen<-b.genus_RA_meta[b.genus_RA_meta$SampDate=="December.2021",]
dec.gen[which.max(dec.gen$Count),] #
dec.gen[order(dec.gen$Count,decreasing=TRUE),]

apr.gen<-b.genus_RA_meta[b.genus_RA_meta$SampDate=="April.2022",]
apr.gen[which.max(apr.gen$Count),] #
apr.gen[order(apr.gen$Count,decreasing=TRUE),]

# Barplot by PlotID

b.gen_RA0<-ggplot(b.genus_RA_meta[b.genus_RA_meta$Count>0.01,], aes(x=PlotID, y=Count, fill=Genus))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genus Relative Abundance", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 1%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=2))

ggsave(b.gen_RA0,filename = "figures/RelativeAbundance/Genus/SSW_16S_Genera.Spec.RA_barplot_1perc.png", width=12, height=10, dpi=600)

b.gen_RA1<-ggplot(b.genus_RA_meta[b.genus_RA_meta$Count>0.05,], aes(x=PlotID, y=Count, fill=Genus))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genus Relative Abundance", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 5%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))

ggsave(b.gen_RA1,filename = "figures/RelativeAbundance/Genus/SSW_16S_Genera.Spec.RA_barplot_5perc.png", width=12, height=10, dpi=600)

b.gen_RA2<-ggplot(b.genus_RA_meta[b.genus_RA_meta$Count>0.10,], aes(x=PlotID, y=Count, fill=Genus))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genus Relative Abundance", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 10%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))

ggsave(b.gen_RA2,filename = "figures/RelativeAbundance/Genus/SSW_16S_Genera.Spec.RA_barplot_10perc.png", width=12, height=10, dpi=600)

b.gen_RA3<-ggplot(b.genus_RA_meta[b.genus_RA_meta$Count>0.25,], aes(x=PlotID, y=Count, fill=Genus))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genus Relative Abundance", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 25%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))

ggsave(b.gen_RA3,filename = "figures/RelativeAbundance/Genus/SSW_16S_Genera.Spec.RA_barplot_25perc.png", width=12, height=10, dpi=600)

b.gen_RA4<-ggplot(b.genus_RA_meta[b.genus_RA_meta$Count>0.35,], aes(x=PlotID, y=Count, fill=Genus))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genus Relative Abundance", x="SampleID", y="Relative Abundance", subtitle="Includes taxa with Relative Abundance > 35%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))

ggsave(b.gen_RA4,filename = "figures/RelativeAbundance/Genus/SSW_16S_Genera.Spec.RA_barplot_35perc.png", width=12, height=10, dpi=600)

b.gen.a21_RA<-ggplot(aug.gen[aug.gen$Count>0.01,], aes(x=PlotID, y=Count, fill=Genus))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genera Relative Abundance - August 2021", x="SampleID", y="Relative Abundance", subtitle="Includes taxa with Relative Abundance > 1%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))+scale_y_continuous(expand = c(0,0),limits = c(0,1))

ggsave(b.gen.a21_RA,filename = "figures/RelativeAbundance/Genus/16S_August2021_Genera.RA_1perc_barplot.png", width=12, height=10, dpi=600)

b.gen.d21_RA<-ggplot(dec.gen[dec.gen$Count>0.01,], aes(x=PlotID, y=Count, fill=Genus))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genera Relative Abundance - December 2021", x="SampleID", y="Relative Abundance", subtitle="Includes taxa with Relative Abundance > 1%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))+scale_y_continuous(expand = c(0,0),limits = c(0,1))

ggsave(b.gen.d21_RA,filename = "figures/RelativeAbundance/Genus/16S_December2021_Genera.RA_1perc_barplot.png", width=12, height=10, dpi=600)

b.gen.a22_RA<-ggplot(apr.gen[apr.gen$Count>0.01,], aes(x=PlotID, y=Count, fill=Genus))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genera Relative Abundance - April 2022", x="SampleID", y="Relative Abundance", subtitle="Includes taxa with Relative Abundance > 1%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))+scale_y_continuous(expand = c(0,0),limits = c(0,1))

ggsave(b.gen.a22_RA,filename = "figures/RelativeAbundance/Genus/16S_April2022_Genera.RA_1perc_barplot.png", width=12, height=10, dpi=600)

all.g.bymo<-ggarrange(b.gen.a21_RA, b.gen.d21_RA, b.gen.a22_RA,ncol = 3, nrow = 1)
ggsave(all.g.bymo,filename = "figures/RelativeAbundance/Genus/16S_MonthComparison_Genera.RA_1perc_barplot.png", width=25, height=10, dpi=600)

b.gen.a21_RA2<-ggplot(aug.gen[aug.gen$Count>0.05,], aes(x=PlotID, y=Count, fill=Genus))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genera Relative Abundance - August 2021", x="SampleID", y="Relative Abundance", subtitle="Includes taxa with Relative Abundance > 5%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))+scale_y_continuous(expand = c(0,0),limits = c(0,1))

ggsave(b.gen.a21_RA2,filename = "figures/RelativeAbundance/Genus/16S_August2021_Genera.RA_5perc_barplot.png", width=12, height=10, dpi=600)

b.gen.d21_RA2<-ggplot(dec.gen[dec.gen$Count>0.05,], aes(x=PlotID, y=Count, fill=Genus))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genera Relative Abundance - December 2021", x="SampleID", y="Relative Abundance", subtitle="Includes taxa with Relative Abundance > 5%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))+scale_y_continuous(expand = c(0,0),limits = c(0,1))

ggsave(b.gen.d21_RA2,filename = "figures/RelativeAbundance/Genus/16S_December2021_Genera.RA_5perc_barplot.png", width=12, height=10, dpi=600)

b.gen.a22_RA2<-ggplot(apr.gen[apr.gen$Count>0.05,], aes(x=PlotID, y=Count, fill=Genus))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genera Relative Abundance - April 2022", x="SampleID", y="Relative Abundance", subtitle="Includes taxa with Relative Abundance > 5%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))+scale_y_continuous(expand = c(0,0),limits = c(0,1))

ggsave(b.gen.a22_RA2,filename = "figures/RelativeAbundance/Genus/16S_April2022_Genera.RA_5perc_barplot.png", width=12, height=10, dpi=600)

all.g.bymo2<-ggarrange(b.gen.a21_RA2, b.gen.d21_RA2, b.gen.a22_RA2,ncol = 3, nrow = 1)
ggsave(all.g.bymo2,filename = "figures/RelativeAbundance/Genus/16S_MonthComparison_Genera.RA_5perc_barplot.png", width=25, height=10, dpi=600)

# prep for heatmap
max(b.genus_RA_meta$Count)
mean(b.genus_RA_meta$Count)
median(b.genus_RA_meta$Count)
max(b.genus_RA_meta$Count)/2 # what is the mid point of the RA here?

# Heatmap by PlotID

g.h1<-ggplot(b.genus_RA_meta[b.genus_RA_meta$Count>0.01,], aes(PlotID, Genus, fill= Count)) +geom_tile()+scale_fill_gradient2(low="orange",mid="white",high="purple",midpoint=0.3)+
  theme_classic()+theme(axis.title.x = element_text(size=13,vjust=-0.5),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Sample ID", y="Microbial Genera", title="Microbial Genera & Sample Type",subtitle="Includes taxa with Relative Abundance > 1%",fill="Relative Abundance")+scale_x_discrete(expand = c(0,0))

ggsave(g.h1,filename = "figures/RelativeAbundance/Genus/16S_Genera.RA_heatmap_A_1perc.png", width=20, height=15, dpi=600)

g.h2<-ggplot(b.genus_RA_meta[b.genus_RA_meta$Count>0.05,], aes(PlotID, Genus, fill= Count)) +geom_tile()+scale_fill_gradient2(low="orange",mid="white",high="purple",midpoint=0.3)+
  theme_classic()+theme(axis.title.x = element_text(size=13,vjust=-0.5),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Sample ID", y="Microbial Genera", title="Microbial Genera & Sample Type",subtitle="Includes taxa with Relative Abundance > 5%",fill="Relative Abundance")+scale_x_discrete(expand = c(0,0))

ggsave(g.h2,filename = "figures/RelativeAbundance/Genus/16S_Genera.RA_heatmap_B_5perc.png", width=16, height=10, dpi=600)

bac.dat.all[1:4,1:4]
#
# by Genus + Sampling Date
bac.gen.date <- as.data.frame(dcast(bac.dat.all.g,SampDate~Genus, value.var="Count", fun.aggregate=sum)) ###
head(bac.gen.date) # counts by Genus + sample depe
rownames(bac.gen.date)<-bac.gen.date$SampDate

b.RA_gen.date<-data.frame(decostand(bac.gen.date[,-1], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(b.RA_gen.date) # sanity check
b.RA_gen.date$SampDate<-rownames(b.RA_gen.date)
head(b.RA_gen.date)

#melt down relativized data to merge with metadata
b.gen.date_m<-melt(b.RA_gen.date, by="SampDate")

head(b.gen.date_m)
colnames(b.gen.date_m)[which(names(b.gen.date_m) == "variable")] <- "Genus"
colnames(b.gen.date_m)[which(names(b.gen.date_m) == "value")] <- "Count"
head(b.gen.date_m) ## relative abundance based on sum of counts by Genus!
b.gen.date_m$Genus<-gsub("^X.","",b.gen.date_m$Genus) # get rid of leading X. in Genus names
b.gen.date_m$Genus<-gsub("\\.\\."," ",b.gen.date_m$Genus) # get rid of .. in species name --> . is regex
head(b.gen.date_m) ## relative abundance based on sum of counts by genus!

b.gen.date_m$SampDate<-factor(b.gen.date_m$SampDate, levels=c("August.2021","December.2021","April.2022"))

# Find most abundant genus per time point
aug.gen.d<-b.gen.date_m[b.gen.date_m$SampDate=="August.2021",]
#aug.gen.d[which.max(aug.gen.d$Count),]
head(aug.gen.d[order(aug.gen.d$Count,decreasing=TRUE),])

dec.gen.d<-b.gen.date_m[b.gen.date_m$SampDate=="December.2021",]
#dec.gen.d[which.max(dec.gen.d$Count),] #
head(dec.gen.d[order(dec.gen.d$Count,decreasing=TRUE),])

apr.gen.d<-b.gen.date_m[b.gen.date_m$SampDate=="April.2022",]
#apr.gen.d[which.max(apr.gen.d$Count),] #
head(apr.gen.d[order(apr.gen.d$Count,decreasing=TRUE),])

# Barplot by Sample Date

gsd1<-ggplot(b.gen.date_m[b.gen.date_m$Count>0.01,], aes(x=SampDate, y=Count, fill=Genus))+geom_bar(stat="identity",colour="black")+theme_classic()+
  labs(title = "Relative Abundance of Microbial Genera", x="Sample Date", subtitle="Includes taxa with Relative Abundance > 1%",y="Relative Abundance", fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=2))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1))+ scale_x_discrete(labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022"))

ggsave(gsd1,filename = "figures/RelativeAbundance/Genus/SSW_16S_Genera.RA_date_barplot.png", width=15, height=10, dpi=600)

# Taxonomic Summary by Sample Date

#b.gen.date_m2<-merge(b.gen.date_m, metadata, by="SampDate")

colorset1 # remember which date goes with each color

gsd1<-ggplot(b.gen.date_m[b.gen.date_m$Count>0.01,], aes(Genus, Count)) +
  geom_jitter(aes(color=factor(SampDate)), size=3, width=0.15, height=0) +
  scale_color_manual(name ="Sample Date", values=c("#ef781c","#03045e","#059c3f"), labels=c("August 2021","December 2021","April 2022")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=0.5,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Microbial Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 1%")

ggsave(gsd1,filename = "figures/RelativeAbundance/Genus/SSW_16S_Genera.RA_date_taxasum_1perc.png", width=15, height=10, dpi=600)

gsd1a<-ggplot(b.gen.date_m[b.gen.date_m$Count>0.01,], aes(Genus, Count)) +
  geom_jitter(aes(color=factor(SampDate)), size=3, width=0.15, height=0) +
  scale_color_manual(name ="Sample Date", values=c("#ef781c","#03045e","#059c3f"), labels=c("August 2021","December 2021","April 2022")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Microbial Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 1%")+coord_flip()

ggsave(gsd1a,filename = "figures/RelativeAbundance/Genus/SSW_16S_Genera.RA_date_taxasum_1perc_v2.png", width=15, height=10, dpi=600)

gsd1b<-ggplot(b.gen.date_m[b.gen.date_m$Count>0.05,], aes(Genus, Count)) +
  geom_jitter(aes(color=factor(SampDate)), size=3, width=0.15, height=0) +
  scale_color_manual(name ="Sample Date", values=c("#ef781c","#03045e","#059c3f"), labels=c("August 2021","December 2021","April 2022")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=0.5,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Microbial Genera & Sample Date")

ggsave(gsd1b,filename = "figures/RelativeAbundance/Genus/SSW_16S_Genera.RA_date_taxasum_5percent.png", width=15, height=10, dpi=600)

gsd1c<-ggplot(b.gen.date_m[b.gen.date_m$Count>0.05,], aes(Genus, Count)) +
  geom_jitter(aes(color=factor(SampDate)), size=3, width=0.15, height=0) +
  scale_color_manual(name ="Sample Date", values=c("#ef781c","#03045e","#059c3f"), labels=c("August 2021","December 2021","April 2022")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Microbial Genera & Sample Date")+coord_flip()

ggsave(gsd1c,filename = "figures/RelativeAbundance/Genus/SSW_16S_Genera.RA_date_taxasum_5percent_v2.png", width=15, height=10, dpi=600)

# by Genus + Sampling Date + Depth
bac.gen.date.dep <- as.data.frame(dcast(bac.dat.all.g,SampDate+Depth_m~Genus, value.var="Count", fun.aggregate=sum)) ###
bac.gen.date.dep[1:5,1:5] # counts by Genus + sample date & depth
rownames(bac.gen.date.dep)<-interaction(bac.gen.date.dep$SampDate,bac.gen.date.dep$Depth_m,sep="_")
bac.gen.date.dep[1:5,1:5]

b.RA_gen.date.dep<-data.frame(decostand(bac.gen.date.dep[,-c(1:2)], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(b.RA_gen.date.dep) # sanity check
b.RA_gen.date.dep$SampDate_Depth<-rownames(b.RA_gen.date.dep)
b.RA_gen.date.dep[1:5,ncol(b.RA_gen.date.dep):5-ncol(b.RA_gen.date.dep)] # first 5 rows, last 5 columns

#melt down relativized data to merge with metadata
b.gen.date.dep_m<-melt(b.RA_gen.date.dep, by="SampDate_Depth")

head(b.gen.date.dep_m)
colnames(b.gen.date.dep_m)[which(names(b.gen.date.dep_m) == "variable")] <- "Genus"
colnames(b.gen.date.dep_m)[which(names(b.gen.date.dep_m) == "value")] <- "Count"
head(b.gen.date.dep_m) ## relative abundance based on sum of counts by Genus!
b.gen.date.dep_m$Genus<-gsub("^X.","",b.gen.date.dep_m$Genus) # get rid of leading X. in Genus names
b.gen.date.dep_m$Genus<-gsub("\\.\\."," ",b.gen.date.dep_m$Genus) # get rid of .. in species name --> . is regex
head(b.gen.date.dep_m) ## relative abundance based on sum of counts by genus!

# separate SampDate_Depth above by _, then recreate column in new df
b.gen.date.dep_m2<-as.data.frame(separate_wider_delim(data = b.gen.date.dep_m, col=SampDate_Depth, "_", names = c("SampDate", "Depth_m"))) # Separate SampDate & Depth column for Heatmap later
b.gen.date.dep_m2$SampDate_Depth<-interaction(b.gen.date.dep_m2$SampDate,b.gen.date.dep_m2$Depth_m)

b.gen.date.dep_m2$Depth_m<-factor(b.gen.date.dep_m2$Depth_m, levels=c("0","3","4","5","7","9","10","10.5"))
b.gen.date.dep_m2$SampDate<-factor(b.gen.date.dep_m2$SampDate,levels=c("August.2021","December.2021","April.2022"))
b.gen.date.dep_m2$SampDate_Depth<-gsub("(\\..*?)\\.","\\1-",b.gen.date.dep_m2$SampDate_Depth)
# \\. - first period, .* is any character after that, ? suppresses greedy matching of .*, () - what we want to keep aka \\1
# \\1 is the pattern we want to keep in (), and - is replacing second character we want to replace with -
# more info here: https://stackoverflow.com/questions/43077846/how-to-replace-second-or-more-occurrences-of-a-dot-from-a-column-name
b.gen.date.dep_m2$SampDate_Depth = factor(b.gen.date.dep_m2$SampDate_Depth, levels=unique(b.gen.date.dep_m2$SampDate_Depth[order(b.gen.date.dep_m2$SampDate,b.gen.date.dep_m2$Depth_m)]), ordered=TRUE)

# Heatmaps
g.sd.d.h1<-ggplot(b.gen.date.dep_m2[b.gen.date.dep_m2$Count>0.01,], aes(SampDate_Depth, Genus, fill= Count)) +geom_tile()+scale_fill_gradient2(low="orange",mid="white",high="purple",midpoint=0.3)+
  theme_classic()+theme(axis.title.x = element_text(size=13,vjust=-0.5),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Sampling Date - Depth (m)", y="Microbial Genera", title="Microbial Genera by Sample Date & Depth",subtitle="Includes taxa with Relative Abundance > 1%",fill="Relative Abundance")+scale_x_discrete(expand = c(0,0))
#+ scale_x_discrete(labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022"))

ggsave(g.sd.d.h1,filename = "figures/RelativeAbundance/Genus/16S_Genera.RA_heatmap_date_depth_1perc.png", width=20, height=15, dpi=600)

g.sd.d.h2<-ggplot(b.gen.date.dep_m2[b.gen.date.dep_m2$Count>0.05,], aes(SampDate_Depth, Genus, fill= Count)) +geom_tile()+scale_fill_gradient2(low="orange",mid="white",high="purple",midpoint=0.3)+
  theme_classic()+theme(axis.title.x = element_text(size=13,vjust=-0.5),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Sampling Date - Depth (m)", y="Microbial Genera", title="Microbial Genera by Sample Date & Depth",subtitle="Includes taxa with Relative Abundance > 1%",fill="Relative Abundance")+scale_x_discrete(expand = c(0,0))
#+ scale_x_discrete(labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022"))

ggsave(g.sd.d.h2,filename = "figures/RelativeAbundance/Genus/16S_Genera.RA_heatmap_date_depth_5perc.png", width=20, height=15, dpi=600)

# Taxonomic Summaries
g.sd.d.hm.1<-ggplot(b.gen.date.dep_m2[b.gen.date.dep_m2$Count>0.01,], aes(Genus, Count)) +
  geom_jitter(aes(color=as.numeric(as.character(Depth_m)),shape=SampDate), size=3, width=0.15, height=0) +
  scale_colour_gradient2(low="red",mid="hotpink",high="blue",midpoint=5.25,guide = guide_colourbar(reverse = TRUE)) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() + scale_shape_discrete(labels=c("August 2021","December 2021","April 2022"),name="Sample Date") +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Microbial Genera by Sample Date & Depth",subtitle="Includes taxa with Relative Abundance > 1%",color="Depth (m)", shape="Sample Date")

ggsave(g.sd.d.hm.1,filename = "figures/RelativeAbundance/Genus/SSW_16S_Genera.RA_date_depth_taxasum_1perc.png", width=18, height=10, dpi=600)
## ^ this figure includes the relative abundance of each organism by depth & date!!!

g.sd.d.hm.2<-ggplot(b.gen.date.dep_m2[b.gen.date.dep_m2$Count>0.025,], aes(Genus, Count)) +
  geom_jitter(aes(color=as.numeric(as.character(Depth_m)),shape=SampDate), size=3, width=0.15, height=0) +
  scale_colour_gradient2(low="red",mid="hotpink",high="blue",midpoint=5.25,guide = guide_colourbar(reverse = TRUE)) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() + scale_shape_discrete(labels=c("August 2021","December 2021","April 2022"),name="Sample Date") +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Microbial Genera by Sample Date & Depth",subtitle="Includes taxa with Relative Abundance > 2.5%",color="Depth (m)", shape="Sample Date")

ggsave(g.sd.d.hm.2,filename = "figures/RelativeAbundance/Genus/SSW_16S_Genera.RA_date_depth_taxasum_2.5perc.png", width=18, height=10, dpi=600)

g.sd.d.hm.2a<-ggplot(b.gen.date.dep_m2[b.gen.date.dep_m2$Count>0.025,], aes(Genus, Count)) +
  geom_jitter(aes(color=as.numeric(as.character(Depth_m)),shape=SampDate), size=6, width=0.15, height=0) +
  scale_colour_gradient2(low="red",mid="hotpink",high="blue",midpoint=5.25,guide = guide_colourbar(reverse = TRUE)) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() + scale_shape_discrete(labels=c("August 2021","December 2021","April 2022"),name="Sample Date") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),axis.text = element_text(size=16),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01,face="bold"),legend.title.align=0.5, legend.title = element_text(size=16),legend.text = element_text(size=14)) +
  labs(x="Microbial Genera", y="Relative Abundance",color="Depth (m)", shape="Sample Date")

ggsave(g.sd.d.hm.2a,filename = "figures/RelativeAbundance/Genus/SSW_16S_Genera.RA_date_depth_taxasum_2.5perc_poster.png", width=18, height=10, dpi=600)

g.sd.d.hm.3<-ggplot(b.gen.date.dep_m2[b.gen.date.dep_m2$Count>0.05,], aes(Genus, Count)) +
  geom_jitter(aes(color=as.numeric(as.character(Depth_m)),shape=SampDate), size=3, width=0.15, height=0) +
  scale_colour_gradient2(low="red",mid="hotpink",high="blue",midpoint=5.25,guide = guide_colourbar(reverse = TRUE)) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() + scale_shape_discrete(labels=c("August 2021","December 2021","April 2022"),name="Sample Date") +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Microbial Genera by Sample Date & Depth",subtitle="Includes taxa with Relative Abundance > 5%",color="Depth (m)", shape="Sample Date")

ggsave(g.sd.d.hm.3,filename = "figures/RelativeAbundance/Genus/SSW_16S_Genera.RA_date_depth_taxasum_5perc.png", width=15, height=10, dpi=600)


#### Create Most Abundant Genus Color Palette ####
# top 10 bacterial genera overall:
# DS001, Litoricola, Truepera, Synechococcus.CC9902, GKS98.freshwater.group unknown, CL500.3 unknown, Clade.Ia unknown, Candidatus.Aquiluna unknown, Brumimicrobium unknown, Fluviicola unknown

gen.col = melt(c(DS001="#a4133c",Litoricola="#3a86ff",Truepera="#8338ec",Synechococcus.CC9902="#55a630",GKS98.freshwater.group="#ff6d00", CL500.3="violet", Clade.Ia="#4cc9f0",
                 Candidatus.Aquiluna="#ff006e",Brumimicrobium="#006d77",Fluviicola="#ffc300"))

gen.col$Genus<-rownames(gen.col)
colnames(gen.col)[which(names(gen.col) == "value")] <- "Gen_Color"
gen.col

# ^ use this to merge with bacterial relative abundance data later

#### Genus Relative Abundance - Top Taxa Only ####

# Most abundant overall
b.gendecreasing<-b.genus_RA_meta[order(b.genus_RA_meta$Count,decreasing=TRUE),]
b.gendecreasing
# top 10 taxa overall:
# DS001, Litoricola, Truepera, Synechococcus.CC9902, GKS98.freshwater.group unknown, CL500.3 unknown, Clade.Ia unknown, Candidatus.Aquiluna unknown, Brumimicrobium unknown, Fluviicola unknown

head(b.genus_RA_meta)
b.genus_RA_meta2<-merge(b.genus_RA_meta, gen.col, by="Genus")
head(b.genus_RA_meta2)
b.genus_RA_meta2$Gen_Color <- as.character(b.genus_RA_meta2$Gen_Color)

#b.genus_RA_meta2$PlotID<-gsub("(.*\\..*\\..*)\\.(.*\\..*)","\\1-\\2",b.genus_RA_meta2$PlotID) # adjust only labels that end in 10.5m to try and keep period
#b.genus_RA_meta2$PlotID<-gsub("(.*\\..*\\..*)\\.","\\1-",b.genus_RA_meta2$PlotID)
# ^ gsub regex explanation
## () == whatever is in parentheses is what I want to keep after gsub
## .* == any characters; \\. == a period (periods are regex, have to escape with \\)
## \\1- == replacing last period with a -, keeping first pattern described in parentheses in tact with \\1
## turns 4.13.22.0m --> 4.13.22-0m

ten.gen<-ggplot(b.genus_RA_meta2, aes(x=PlotID, y=Count, fill=Genus))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genus Relative Abundance", x="SampleID", y="Relative Abundance", subtitle="Includes Ten Most Relatively Abundant Genera",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))+scale_fill_manual(name ="Genus",values=unique(b.genus_RA_meta2$Gen_Color[order(b.genus_RA_meta2$Genus)]),labels=c(unique(b.genus_RA_meta2$Genus[order(b.genus_RA_meta2$Genus)])))

ggsave(ten.gen,filename = "figures/RelativeAbundance/Genus/SSW_16S_TenMostAbundant_Genera_barplot.png", width=12, height=10, dpi=600)

#### Genus Relative Abundance - August ####
# use dcast to count up ASVs within each Genus across aug of the samples
head(bac.dat.all)
bac.dat.aug<-bac.dat.all[bac.dat.all$SampleMonth=="August",]
bac.dat.aug.g<-subset(bac.dat.aug, bac.dat.aug$Genus!="Unknown") # drop unknown genera so they don't skew analyses
"Unknown" %in% bac.dat.aug.g$Genus

aug.genus_counts <- as.data.frame(dcast(bac.dat.aug.g, SampleID~Genus+Species, value.var="Count", fun.aggregate=sum)) ###
head(aug.genus_counts) # counts by genus per sample
dim(aug.genus_counts)
rownames(aug.genus_counts)<-aug.genus_counts$SampleID
aug.genus_counts[1:4,1:4]

aug.genus_RelAb<-data.frame(decostand(aug.genus_counts[,-1], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(aug.genus_RelAb) # sanity check to make sure the transformation worked!

aug.genus_RelAb$SampleID<-rownames(aug.genus_RelAb)
head(aug.genus_RelAb)
#write.csv(aug.genus_RelAb,"16S_Genera_Relative_Abundance.csv",row.names=TRUE) # good to save just in case

# melt down relativized data to merge with metadata
aug.genus_m<-melt(aug.genus_RelAb)

head(aug.genus_m)
colnames(aug.genus_m)[which(names(aug.genus_m) == "variable")] <- "Genus_species"
colnames(aug.genus_m)[which(names(aug.genus_m) == "value")] <- "Count"
head(aug.genus_m) ## relative abundance based on sum of counts by genus!
aug.genus_m$Genus_species<-gsub("^X.","",aug.genus_m$Genus_species) # get rid of leading X. in Genus_species names
aug.genus_m$Genus_species<-gsub("\\.\\."," ",aug.genus_m$Genus_species) # get rid of .. in species name --> . is regex
aug.genus_m$Genus_species<-gsub("_"," ",aug.genus_m$Genus_species) # get rid of .. in species name --> . is regex
head(aug.genus_m) ## relative abundance based on sum of counts by genus!

aug.genus_RA_meta<-merge(aug.genus_m,metadata, by="SampleID")
head(aug.genus_RA_meta) ## relative abundance based on sum of counts by genus!
max(aug.genus_RA_meta$Count)

# rename SampleID so it's easy to interpret, then reorder SampleID based on sample date, then depth
aug.genus_RA_meta$PlotID = gsub("^SSW.","",aug.genus_RA_meta$SampleID)
aug.genus_RA_meta$PlotID = factor(aug.genus_RA_meta$PlotID, levels=unique(aug.genus_RA_meta$PlotID[order(aug.genus_RA_meta$Depth_m)]), ordered=TRUE)
#
# # Barplot by PlotID
#
b.gen_RA0<-ggplot(aug.genus_RA_meta[aug.genus_RA_meta$Count>0.01,], aes(x=PlotID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genus Relative Abundance", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 1%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=2))

ggsave(b.gen_RA0,filename = "figures/RelativeAbundance/Genus/August2021/SSW_16S_Genera.Spec.RA_barplot_1perc.png", width=12, height=10, dpi=600)

b.gen_RA1<-ggplot(aug.genus_RA_meta[aug.genus_RA_meta$Count>0.05,], aes(x=PlotID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genus Relative Abundance", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 5%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))

ggsave(b.gen_RA1,filename = "figures/RelativeAbundance/Genus/August2021/SSW_16S_Genera.Spec.RA_barplot_5perc.png", width=12, height=10, dpi=600)

# by Genus + depth
aug.gen.dep <- as.data.frame(dcast(bac.dat.aug.g,Depth_m~Genus, value.var="Count", fun.aggregate=sum)) ###
head(aug.gen.dep) # counts by Genus + sample depth
rownames(aug.gen.dep)<-aug.gen.dep$Depth_m

aug.RA_gen.dep<-data.frame(decostand(aug.gen.dep[,-1], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(aug.RA_gen.dep) # sanity check
aug.RA_gen.dep$Depth_m<-rownames(aug.RA_gen.dep) # Depth_m is now a character, not a factor!
head(aug.RA_gen.dep)

#melt down relativized data to merge with metadata
aug.gen.dep_m<-melt(aug.RA_gen.dep, by="Depth_m")

head(aug.gen.dep_m)
colnames(aug.gen.dep_m)[which(names(aug.gen.dep_m) == "variable")] <- "Genus"
colnames(aug.gen.dep_m)[which(names(aug.gen.dep_m) == "value")] <- "Count"
head(aug.gen.dep_m) ## relative abundance based on sum of counts by Genus!
aug.gen.dep_m$Genus<-gsub("^X.","",aug.gen.dep_m$Genus) # get rid of leading X. in Genus names
aug.gen.dep_m$Genus<-gsub("\\.\\."," ",aug.gen.dep_m$Genus) # get rid of .. in species name --> . is regex
head(aug.gen.dep_m) ## relative abundance based on sum of counts by genus!
unique(aug.gen.dep_m$Depth_m)
aug.gen.dep_m$Depth_m<-factor(aug.gen.dep_m$Depth_m, levels=c("0","3","4","5","7","9","10","10.5"))

# Barplot by Depth

gd1<-ggplot(aug.gen.dep_m[aug.gen.dep_m$Count>0.01,], aes(x=Depth_m, y=Count, fill=Genus))+geom_bar(stat="identity",colour="black")+theme_classic()+
  labs(title = "Relative Abundance of Microbial Genera - August 2021", x="Depth (m)", y="Relative Abundance", fill="Class",subtitle="Only Relative Abundance > 1%")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1))+coord_flip() + scale_x_discrete(limits = rev(levels(aug.gen.dep_m$Depth_m)))
#+ scale_x_discrete(labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022"))

ggsave(gd1,filename = "figures/RelativeAbundance/Genus/August2021/SSW_16S_Aug21_Genus.RA_barplot_depth_1percent_B.png", width=12, height=10, dpi=600)

gd1a<-ggplot(aug.gen.dep_m[aug.gen.dep_m$Count>0.05,], aes(x=Depth_m, y=Count, fill=Genus))+geom_bar(stat="identity",colour="black")+theme_classic()+
  labs(title = "Relative Abundance of Microbial Genera - August 2021", x="Depth (m)", y="Relative Abundance", fill="Class",subtitle="Only Relative Abundance > 5%")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1))
#+ scale_x_discrete(labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022"))

ggsave(gd1a,filename = "figures/RelativeAbundance/Genus/August2021/SSW_16S_Aug21_Genus.RA_barplot_depth_5percent_A.png", width=12, height=10, dpi=600)

gd1b<-ggplot(aug.gen.dep_m[aug.gen.dep_m$Count>0.05,], aes(x=Depth_m, y=Count, fill=Genus))+geom_bar(stat="identity",colour="black")+theme_classic()+
  labs(title = "Relative Abundance of Microbial Genera - August 2021", x="Depth (m)", y="Relative Abundance", fill="Class",subtitle="Only Relative Abundance > 5%")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1))+coord_flip() + scale_x_discrete(limits = rev(levels(aug.gen.dep_m$Depth_m)))

ggsave(gd1b,filename = "figures/RelativeAbundance/Genus/August2021/SSW_16S_Aug21_Genus.RA_barplot_depth_5percent_B.png", width=12, height=10, dpi=600)

# Taxonomic Summary by Depth

#dep_meta<-unique(data.frame("Depth_m"=metadata$Depth_m, "Sample_Color"=metadata$Sample_Color))
#p_dep_meta<-merge(dep_meta,aug.gen.dep_m, by="Depth_m")
tg1<-ggplot(aug.gen.dep_m[aug.gen.dep_m$Count>0.01,], aes(Genus, Count)) +
  geom_jitter(aes(color=as.numeric(Depth_m)), size=3, width=0.15, height=0) +
  scale_colour_gradient2(low="red",high="blue3",midpoint=5,guide = guide_colourbar(reverse = TRUE)) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Microbial Genera & Depth - August 2021", subtitle="Includes taxa with Relative Abundance > 1%",color="Depth (m)")

ggsave(tg1,filename = "figures/RelativeAbundance/Genus/August2021/SSW_16S_Aug21_Genera.RA_depth_taxasum_1perc.png", width=15, height=10, dpi=600)

tg1a<-ggplot(aug.gen.dep_m[aug.gen.dep_m$Count>0.01,], aes(Genus, Count)) +
  geom_jitter(aes(color=as.numeric(Depth_m)), size=3, width=0.15, height=0) +
  scale_colour_gradient2(low="red",high="blue3",midpoint=5,guide = guide_colourbar(reverse = TRUE)) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Microbial Genera & Depth - August 2021", subtitle="Includes taxa with Relative Abundance > 1%",color="Depth (m)")+coord_flip()

ggsave(tg1a,filename = "figures/RelativeAbundance/Genus/August2021/SSW_16S_Aug21_Genera.RA_depth_taxasum_1perc_v2.png", width=15, height=10, dpi=600)

tg1b<-ggplot(aug.gen.dep_m[aug.gen.dep_m$Count>0.05,], aes(Genus, Count)) +
  geom_jitter(aes(color=as.numeric(Depth_m)), size=3, width=0.15, height=0) +
  scale_colour_gradient2(low="red",high="blue3",midpoint=5,guide = guide_colourbar(reverse = TRUE)) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Microbial Genera & Depth - August 2021", subtitle="Includes taxa with Relative Abundance > 5%",color="Depth (m)")

ggsave(tg1b,filename = "figures/RelativeAbundance/Genus/August2021/SSW_16S_Aug21_Genera.RA_depth_taxasum_5percent.png", width=15, height=10, dpi=600)

tg1c<-ggplot(aug.gen.dep_m[aug.gen.dep_m$Count>0.05,], aes(Genus, Count)) +
  geom_jitter(aes(color=as.numeric(Depth_m)), size=3, width=0.15, height=0) +
  scale_colour_gradient2(low="red",high="blue3",midpoint=5,guide = guide_colourbar(reverse = TRUE)) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Microbial Genera & Depth - August 2021", subtitle="Includes taxa with Relative Abundance > 5%",color="Depth (m)")+coord_flip()

ggsave(tg1c,filename = "figures/RelativeAbundance/Genus/August2021/SSW_16S_Aug21_Genera.RA_depth_taxasum_5percent_v2.png", width=15, height=10, dpi=600)


#### Genus Relative Abundance - December ####
# use dcast to count up ASVs within each Genus across dec of the samples
head(bac.dat.all)
bac.dat.dec<-bac.dat.all[bac.dat.all$SampleMonth=="December",]
bac.dat.dec.g<-subset(bac.dat.dec, bac.dat.dec$Genus!="Unknown") # drop unknown genera so they don't skew analyses
"Unknown" %in% bac.dat.dec.g$Genus

dec.genus_counts <- as.data.frame(dcast(bac.dat.dec.g, SampleID~Genus+Species, value.var="Count", fun.aggregate=sum)) ###
head(dec.genus_counts) # counts by genus per sample
dim(dec.genus_counts)
rownames(dec.genus_counts)<-dec.genus_counts$SampleID
dec.genus_counts[1:4,1:4]

dec.genus_RelAb<-data.frame(decostand(dec.genus_counts[,-1], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(dec.genus_RelAb) # sanity check to make sure the transformation worked!

dec.genus_RelAb$SampleID<-rownames(dec.genus_RelAb)
head(dec.genus_RelAb)
#write.csv(dec.genus_RelAb,"16S_Genera_Relative_Abundance.csv",row.names=TRUE) # good to save just in case

# melt down relativized data to merge with metadata
dec.genus_m<-melt(dec.genus_RelAb)

head(dec.genus_m)
colnames(dec.genus_m)[which(names(dec.genus_m) == "variable")] <- "Genus_species"
colnames(dec.genus_m)[which(names(dec.genus_m) == "value")] <- "Count"
head(dec.genus_m) ## relative abundance based on sum of counts by genus!
dec.genus_m$Genus_species<-gsub("^X.","",dec.genus_m$Genus_species) # get rid of leading X. in Genus_species names
dec.genus_m$Genus_species<-gsub("\\.\\."," ",dec.genus_m$Genus_species) # get rid of .. in species name --> . is regex
dec.genus_m$Genus_species<-gsub("_"," ",dec.genus_m$Genus_species) # get rid of .. in species name --> . is regex
head(dec.genus_m) ## relative abundance based on sum of counts by genus!

dec.genus_RA_meta<-merge(dec.genus_m,metadata, by="SampleID")
head(dec.genus_RA_meta) ## relative abundance based on sum of counts by genus!
max(dec.genus_RA_meta$Count)

# rename SampleID so it's easy to interpret, then reorder SampleID based on sample date, then depth
dec.genus_RA_meta$PlotID = gsub("^SSW.","",dec.genus_RA_meta$SampleID)
dec.genus_RA_meta$PlotID = factor(dec.genus_RA_meta$PlotID, levels=unique(dec.genus_RA_meta$PlotID[order(dec.genus_RA_meta$SampDate,dec.genus_RA_meta$Depth_m)]), ordered=TRUE)
#
dec.genus_RA_meta[which.max(dec.genus_RA_meta$Count),] #
dec.gen[order(dec.genus_RA_meta$Count,decreasing=TRUE),]
# # Barplot by PlotID
#
# b.gen_RA0<-ggplot(dec.genus_RA_meta[dec.genus_RA_meta$Count>0.01,], aes(x=PlotID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
#   labs(title = "Microbial Genus Relative Abundance", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 1%",fill="Genus")+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
#   guides(fill=guide_legend(ncol=2))
#
# ggsave(b.gen_RA0,filename = "figures/RelativeAbundance/Genus/SSW_16S_Genera.Spec.RA_barplot_1perc.png", width=12, height=10, dpi=600)
#
# b.gen_RA1<-ggplot(dec.genus_RA_meta[dec.genus_RA_meta$Count>0.05,], aes(x=PlotID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
#   labs(title = "Microbial Genus Relative Abundance", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 5%",fill="Genus")+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
#   guides(fill=guide_legend(ncol=1))
#
# ggsave(b.gen_RA1,filename = "figures/RelativeAbundance/Genus/SSW_16S_Genera.Spec.RA_barplot_5perc.png", width=12, height=10, dpi=600)
#
# b.gen_RA2<-ggplot(dec.genus_RA_meta[dec.genus_RA_meta$Count>0.10,], aes(x=PlotID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
#   labs(title = "Microbial Genus Relative Abundance", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 10%",fill="Genus")+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
#   guides(fill=guide_legend(ncol=1))
#
# ggsave(b.gen_RA2,filename = "figures/RelativeAbundance/Genus/SSW_16S_Genera.Spec.RA_barplot_10perc.png", width=12, height=10, dpi=600)
#
# b.gen_RA3<-ggplot(dec.genus_RA_meta[dec.genus_RA_meta$Count>0.25,], aes(x=PlotID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
#   labs(title = "Microbial Genus Relative Abundance", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 25%",fill="Genus")+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
#   guides(fill=guide_legend(ncol=1))
#
# ggsave(b.gen_RA3,filename = "figures/RelativeAbundance/Genus/SSW_16S_Genera.Spec.RA_barplot_25perc.png", width=12, height=10, dpi=600)
#
# b.gen_RA4<-ggplot(dec.genus_RA_meta[dec.genus_RA_meta$Count>0.35,], aes(x=PlotID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
#   labs(title = "Microbial Genus Relative Abundance", x="SampleID", y="Relative Abundance", subtitle="Includes taxa with Relative Abundance > 35%",fill="Genus")+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
#   guides(fill=guide_legend(ncol=1))
#
# ggsave(b.gen_RA4,filename = "figures/RelativeAbundance/Genus/SSW_16S_Genera.Spec.RA_barplot_35perc.png", width=12, height=10, dpi=600)
#
# # prep for heatmap
# max(dec.genus_RA_meta$Count)
# mean(dec.genus_RA_meta$Count)
# median(dec.genus_RA_meta$Count)
# max(dec.genus_RA_meta$Count)/2 # what is the mid point of the RA here?
#
# # Heatmap by PlotID
#
# g.h1<-ggplot(dec.genus_RA_meta[dec.genus_RA_meta$Count>0.01,], aes(PlotID, Genus_species, fill= Count)) +geom_tile()+scale_fill_gradient2(low="orange",mid="white",high="purple",midpoint=0.3)+
#   theme_classic()+theme(axis.title.x = element_text(size=13,vjust=-0.5),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
#   labs(x="Sample ID", y="Microbial Genera", title="Microbial Genera & Sample Type",subtitle="Includes taxa with Relative Abundance > 1%",fill="Relative Abundance")+scale_x_discrete(expand = c(0,0))
#
# ggsave(g.h1,filename = "figures/RelativeAbundance/Genus/16S_Genera.RA_heatmap_A_1perc.png", width=20, height=15, dpi=600)
#
# g.h2<-ggplot(dec.genus_RA_meta[dec.genus_RA_meta$Count>0.05,], aes(PlotID, Genus_species, fill= Count)) +geom_tile()+scale_fill_gradient2(low="orange",mid="white",high="purple",midpoint=0.3)+
#   theme_classic()+theme(axis.title.x = element_text(size=13,vjust=-0.5),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
#   labs(x="Sample ID", y="Microbial Genera", title="Microbial Genera & Sample Type",subtitle="Includes taxa with Relative Abundance > 5%",fill="Relative Abundance")+scale_x_discrete(expand = c(0,0))
#
# ggsave(g.h2,filename = "figures/RelativeAbundance/Genus/16S_Genera.RA_heatmap_B_5perc.png", width=16, height=10, dpi=600)
#
# bac.dat.dec[1:4,1:4]

# by Genus + depth
dec.gen.dep <- as.data.frame(dcast(bac.dat.dec.g,Depth_m~Genus, value.var="Count", fun.aggregate=sum)) ###
head(dec.gen.dep) # counts by Genus + sample depth
rownames(dec.gen.dep)<-dec.gen.dep$Depth_m

dec.RA_gen.dep<-data.frame(decostand(dec.gen.dep[,-1], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(dec.RA_gen.dep) # sanity check
dec.RA_gen.dep$Depth_m<-rownames(dec.RA_gen.dep) # Depth_m is now a character, not a factor!
head(dec.RA_gen.dep)

#melt down relativized data to merge with metadata
dec.gen.dep_m<-melt(dec.RA_gen.dep, by="Depth_m")

head(dec.gen.dep_m)
colnames(dec.gen.dep_m)[which(names(dec.gen.dep_m) == "variable")] <- "Genus"
colnames(dec.gen.dep_m)[which(names(dec.gen.dep_m) == "value")] <- "Count"
head(dec.gen.dep_m) ## relative abundance based on sum of counts by Genus!
dec.gen.dep_m$Genus<-gsub("^X.","",dec.gen.dep_m$Genus) # get rid of leading X. in Genus names
dec.gen.dep_m$Genus<-gsub("\\.\\."," ",dec.gen.dep_m$Genus) # get rid of .. in species name --> . is regex
head(dec.gen.dep_m) ## relative abundance based on sum of counts by genus!
unique(dec.gen.dep_m$Depth_m)
dec.gen.dep_m$Depth_m<-factor(dec.gen.dep_m$Depth_m, levels=c("0","3","4","5","7","9","10","10.5"))

# Barplot by Depth

gd1<-ggplot(dec.gen.dep_m[dec.gen.dep_m$Count>0.01,], aes(x=Depth_m, y=Count, fill=Genus))+geom_bar(stat="identity",colour="black")+theme_classic()+
  labs(title = "Relative Abundance of Microbial Genera - December 2021", x="Depth (m)", y="Relative Abundance", fill="Class",subtitle="Only Relative Abundance > 1%")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1))+coord_flip() + scale_x_discrete(limits = rev(levels(dec.gen.dep_m$Depth_m)))

ggsave(gd1,filename = "figures/RelativeAbundance/Genus/December2021/SSW_16S_Dec21_Genus.RA_barplot_depth_1percent_B.png", width=12, height=10, dpi=600)

gd1a<-ggplot(dec.gen.dep_m[dec.gen.dep_m$Count>0.05,], aes(x=Depth_m, y=Count, fill=Genus))+geom_bar(stat="identity",colour="black")+theme_classic()+
  labs(title = "Relative Abundance of Microbial Genera - December 2021", x="Depth (m)", y="Relative Abundance", fill="Class",subtitle="Only Relative Abundance > 5%")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1))

ggsave(gd1a,filename = "figures/RelativeAbundance/Genus/December2021/SSW_16S_Dec21_Genus.RA_barplot_depth_5percent_A.png", width=12, height=10, dpi=600)

gd1b<-ggplot(dec.gen.dep_m[dec.gen.dep_m$Count>0.05,], aes(x=Depth_m, y=Count, fill=Genus))+geom_bar(stat="identity",colour="black")+theme_classic()+
  labs(title = "Relative Abundance of Microbial Genera - December 2021", x="Depth (m)", y="Relative Abundance", fill="Class",subtitle="Only Relative Abundance > 5%")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1))+coord_flip() + scale_x_discrete(limits = rev(levels(dec.gen.dep_m$Depth_m)))

ggsave(gd1b,filename = "figures/RelativeAbundance/Genus/December2021/SSW_16S_Dec21_Genus.RA_barplot_depth_5percent_B.png", width=12, height=10, dpi=600)

# Taxonomic Summary by Depth

#dep_meta<-unique(data.frame("Depth_m"=metadata$Depth_m, "Sample_Color"=metadata$Sample_Color))
#p_dep_meta<-merge(dep_meta,dec.gen.dep_m, by="Depth_m")
tg1<-ggplot(dec.gen.dep_m[dec.gen.dep_m$Count>0.01,], aes(Genus, Count)) +
  geom_jitter(aes(color=as.numeric(Depth_m)), size=3, width=0.15, height=0) +
  scale_colour_gradient2(low="red",high="blue3",midpoint=5,guide = guide_colourbar(reverse = TRUE)) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Microbial Genera & Depth - December 2021", subtitle="Includes taxa with Relative Abundance > 1%",color="Depth (m)")

ggsave(tg1,filename = "figures/RelativeAbundance/Genus/December2021/SSW_16S_Dec21_Genera.RA_depth_taxasum_1perc.png", width=15, height=10, dpi=600)

tg1a<-ggplot(dec.gen.dep_m[dec.gen.dep_m$Count>0.01,], aes(Genus, Count)) +
  geom_jitter(aes(color=as.numeric(Depth_m)), size=3, width=0.15, height=0) +
  scale_colour_gradient2(low="red",high="blue3",midpoint=5,guide = guide_colourbar(reverse = TRUE)) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Microbial Genera & Depth - December 2021", subtitle="Includes taxa with Relative Abundance > 1%",color="Depth (m)")+coord_flip()

ggsave(tg1a,filename = "figures/RelativeAbundance/Genus/December2021/SSW_16S_Dec21_Genera.RA_depth_taxasum_1perc_v2.png", width=15, height=10, dpi=600)

tg1b<-ggplot(dec.gen.dep_m[dec.gen.dep_m$Count>0.05,], aes(Genus, Count)) +
  geom_jitter(aes(color=as.numeric(Depth_m)), size=3, width=0.15, height=0) +
  scale_colour_gradient2(low="red",high="blue3",midpoint=5,guide = guide_colourbar(reverse = TRUE)) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Microbial Genera & Depth - December 2021", subtitle="Includes taxa with Relative Abundance > 5%",color="Depth (m)")

ggsave(tg1b,filename = "figures/RelativeAbundance/Genus/December2021/SSW_16S_Dec21_Genera.RA_depth_taxasum_5percent.png", width=15, height=10, dpi=600)

tg1c<-ggplot(dec.gen.dep_m[dec.gen.dep_m$Count>0.05,], aes(Genus, Count)) +
  geom_jitter(aes(color=as.numeric(Depth_m)), size=3, width=0.15, height=0) +
  scale_colour_gradient2(low="red",high="blue3",midpoint=5,guide = guide_colourbar(reverse = TRUE)) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Microbial Genera & Depth - December 2021", subtitle="Includes taxa with Relative Abundance > 5%",color="Depth (m)")+coord_flip()

ggsave(tg1c,filename = "figures/RelativeAbundance/Genus/December2021/SSW_16S_Dec21_Genera.RA_depth_taxasum_5percent_v2.png", width=15, height=10, dpi=600)


#### Genus Relative Abundance - April ####
# use dcast to count up ASVs within each Genus across apr of the samples
head(bac.dat.all)
bac.dat.apr<-bac.dat.all[bac.dat.all$SampleMonth=="April",]
bac.dat.apr.g<-subset(bac.dat.apr, bac.dat.apr$Genus!="Unknown") # drop unknown genera so they don't skew analyses
"Unknown" %in% bac.dat.apr.g$Genus

apr.genus_counts <- as.data.frame(dcast(bac.dat.apr.g, SampleID~Genus+Species, value.var="Count", fun.aggregate=sum)) ###
head(apr.genus_counts) # counts by genus per sample
dim(apr.genus_counts)
rownames(apr.genus_counts)<-apr.genus_counts$SampleID
apr.genus_counts[1:4,1:4]

apr.genus_RelAb<-data.frame(decostand(apr.genus_counts[,-1], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(apr.genus_RelAb) # sanity check to make sure the transformation worked!

apr.genus_RelAb$SampleID<-rownames(apr.genus_RelAb)
head(apr.genus_RelAb)
#write.csv(apr.genus_RelAb,"16S_Genera_Relative_Abundance.csv",row.names=TRUE) # good to save just in case

# melt down relativized data to merge with metadata
apr.genus_m<-melt(apr.genus_RelAb)

head(apr.genus_m)
colnames(apr.genus_m)[which(names(apr.genus_m) == "variable")] <- "Genus_species"
colnames(apr.genus_m)[which(names(apr.genus_m) == "value")] <- "Count"
head(apr.genus_m) ## relative abundance based on sum of counts by genus!
apr.genus_m$Genus_species<-gsub("^X.","",apr.genus_m$Genus_species) # get rid of leading X. in Genus_species names
apr.genus_m$Genus_species<-gsub("\\.\\."," ",apr.genus_m$Genus_species) # get rid of .. in species name --> . is regex
apr.genus_m$Genus_species<-gsub("_"," ",apr.genus_m$Genus_species) # get rid of .. in species name --> . is regex
head(apr.genus_m) ## relative abundance based on sum of counts by genus!

apr.genus_RA_meta<-merge(apr.genus_m,metadata, by="SampleID")
head(apr.genus_RA_meta) ## relative abundance based on sum of counts by genus!
max(apr.genus_RA_meta$Count)

# rename SampleID so it's easy to interpret, then reorder SampleID based on sample date, then depth
apr.genus_RA_meta$PlotID = gsub("^SSW.","",apr.genus_RA_meta$SampleID)
apr.genus_RA_meta$PlotID = factor(apr.genus_RA_meta$PlotID, levels=unique(apr.genus_RA_meta$PlotID[order(apr.genus_RA_meta$SampDate,apr.genus_RA_meta$Depth_m)]), ordered=TRUE)
#
# # Barplot by PlotID
#
# b.gen_RA0<-ggplot(apr.genus_RA_meta[apr.genus_RA_meta$Count>0.01,], aes(x=PlotID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
#   labs(title = "Microbial Genus Relative Abundance", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 1%",fill="Genus")+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
#   guides(fill=guide_legend(ncol=2))
#
# ggsave(b.gen_RA0,filename = "figures/RelativeAbundance/Genus/SSW_16S_Genera.Spec.RA_barplot_1perc.png", width=12, height=10, dpi=600)
#
# b.gen_RA1<-ggplot(apr.genus_RA_meta[apr.genus_RA_meta$Count>0.05,], aes(x=PlotID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
#   labs(title = "Microbial Genus Relative Abundance", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 5%",fill="Genus")+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
#   guides(fill=guide_legend(ncol=1))
#
# ggsave(b.gen_RA1,filename = "figures/RelativeAbundance/Genus/SSW_16S_Genera.Spec.RA_barplot_5perc.png", width=12, height=10, dpi=600)
#
# b.gen_RA2<-ggplot(apr.genus_RA_meta[apr.genus_RA_meta$Count>0.10,], aes(x=PlotID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
#   labs(title = "Microbial Genus Relative Abundance", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 10%",fill="Genus")+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
#   guides(fill=guide_legend(ncol=1))
#
# ggsave(b.gen_RA2,filename = "figures/RelativeAbundance/Genus/SSW_16S_Genera.Spec.RA_barplot_10perc.png", width=12, height=10, dpi=600)
#
# b.gen_RA3<-ggplot(apr.genus_RA_meta[apr.genus_RA_meta$Count>0.25,], aes(x=PlotID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
#   labs(title = "Microbial Genus Relative Abundance", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 25%",fill="Genus")+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
#   guides(fill=guide_legend(ncol=1))
#
# ggsave(b.gen_RA3,filename = "figures/RelativeAbundance/Genus/SSW_16S_Genera.Spec.RA_barplot_25perc.png", width=12, height=10, dpi=600)
#
# b.gen_RA4<-ggplot(apr.genus_RA_meta[apr.genus_RA_meta$Count>0.35,], aes(x=PlotID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
#   labs(title = "Microbial Genus Relative Abundance", x="SampleID", y="Relative Abundance", subtitle="Includes taxa with Relative Abundance > 35%",fill="Genus")+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
#   guides(fill=guide_legend(ncol=1))
#
# ggsave(b.gen_RA4,filename = "figures/RelativeAbundance/Genus/SSW_16S_Genera.Spec.RA_barplot_35perc.png", width=12, height=10, dpi=600)
#
# # prep for heatmap
# max(apr.genus_RA_meta$Count)
# mean(apr.genus_RA_meta$Count)
# median(apr.genus_RA_meta$Count)
# max(apr.genus_RA_meta$Count)/2 # what is the mid point of the RA here?
#
# # Heatmap by PlotID
#
# g.h1<-ggplot(apr.genus_RA_meta[apr.genus_RA_meta$Count>0.01,], aes(PlotID, Genus_species, fill= Count)) +geom_tile()+scale_fill_gradient2(low="orange",mid="white",high="purple",midpoint=0.3)+
#   theme_classic()+theme(axis.title.x = element_text(size=13,vjust=-0.5),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
#   labs(x="Sample ID", y="Microbial Genera", title="Microbial Genera & Sample Type",subtitle="Includes taxa with Relative Abundance > 1%",fill="Relative Abundance")+scale_x_discrete(expand = c(0,0))
#
# ggsave(g.h1,filename = "figures/RelativeAbundance/Genus/16S_Genera.RA_heatmap_A_1perc.png", width=20, height=15, dpi=600)
#
# g.h2<-ggplot(apr.genus_RA_meta[apr.genus_RA_meta$Count>0.05,], aes(PlotID, Genus_species, fill= Count)) +geom_tile()+scale_fill_gradient2(low="orange",mid="white",high="purple",midpoint=0.3)+
#   theme_classic()+theme(axis.title.x = element_text(size=13,vjust=-0.5),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
#   labs(x="Sample ID", y="Microbial Genera", title="Microbial Genera & Sample Type",subtitle="Includes taxa with Relative Abundance > 5%",fill="Relative Abundance")+scale_x_discrete(expand = c(0,0))
#
# ggsave(g.h2,filename = "figures/RelativeAbundance/Genus/16S_Genera.RA_heatmap_B_5perc.png", width=16, height=10, dpi=600)
#
# bac.dat.apr[1:4,1:4]

# by Genus + depth
apr.gen.dep <- as.data.frame(dcast(bac.dat.apr.g,Depth_m~Genus, value.var="Count", fun.aggregate=sum)) ###
head(apr.gen.dep) # counts by Genus + sample depth
rownames(apr.gen.dep)<-apr.gen.dep$Depth_m

apr.RA_gen.dep<-data.frame(decostand(apr.gen.dep[,-1], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(apr.RA_gen.dep) # sanity check
apr.RA_gen.dep$Depth_m<-rownames(apr.RA_gen.dep) # Depth_m is now a character, not a factor!
head(apr.RA_gen.dep)

#melt down relativized data to merge with metadata
apr.gen.dep_m<-melt(apr.RA_gen.dep, by="Depth_m")

head(apr.gen.dep_m)
colnames(apr.gen.dep_m)[which(names(apr.gen.dep_m) == "variable")] <- "Genus"
colnames(apr.gen.dep_m)[which(names(apr.gen.dep_m) == "value")] <- "Count"
head(apr.gen.dep_m) ## relative abundance based on sum of counts by Genus!
apr.gen.dep_m$Genus<-gsub("^X.","",apr.gen.dep_m$Genus) # get rid of leading X. in Genus names
apr.gen.dep_m$Genus<-gsub("\\.\\."," ",apr.gen.dep_m$Genus) # get rid of .. in species name --> . is regex
head(apr.gen.dep_m) ## relative abundance based on sum of counts by genus!
unique(apr.gen.dep_m$Depth_m)
apr.gen.dep_m$Depth_m<-factor(apr.gen.dep_m$Depth_m, levels=c("0","3","4","5","7","9","10","10.5"))

# Barplot by Depth

gd1<-ggplot(apr.gen.dep_m[apr.gen.dep_m$Count>0.01,], aes(x=Depth_m, y=Count, fill=Genus))+geom_bar(stat="identity",colour="black")+theme_classic()+
  labs(title = "Relative Abundance of Microbial Genera - April 2022", x="Depth (m)", y="Relative Abundance", fill="Class",subtitle="Only Relative Abundance > 1%")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1))+coord_flip() + scale_x_discrete(limits = rev(levels(apr.gen.dep_m$Depth_m)))

ggsave(gd1,filename = "figures/RelativeAbundance/Genus/April2022/SSW_16S_Apr22_Genus.RA_barplot_depth_1percent_B.png", width=12, height=10, dpi=600)

gd1a<-ggplot(apr.gen.dep_m[apr.gen.dep_m$Count>0.05,], aes(x=Depth_m, y=Count, fill=Genus))+geom_bar(stat="identity",colour="black")+theme_classic()+
  labs(title = "Relative Abundance of Microbial Genera - April 2022", x="Depth (m)", y="Relative Abundance", fill="Class",subtitle="Only Relative Abundance > 5%")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1))

ggsave(gd1a,filename = "figures/RelativeAbundance/Genus/April2022/SSW_16S_Apr22_Genus.RA_barplot_depth_5percent_A.png", width=12, height=10, dpi=600)

gd1b<-ggplot(apr.gen.dep_m[apr.gen.dep_m$Count>0.05,], aes(x=Depth_m, y=Count, fill=Genus))+geom_bar(stat="identity",colour="black")+theme_classic()+
  labs(title = "Relative Abundance of Microbial Genera - April 2022", x="Depth (m)", y="Relative Abundance", fill="Class",subtitle="Only Relative Abundance > 5%")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1))+coord_flip() + scale_x_discrete(limits = rev(levels(apr.gen.dep_m$Depth_m)))

ggsave(gd1b,filename = "figures/RelativeAbundance/Genus/April2022/SSW_16S_Apr22_Genus.RA_barplot_depth_5percent_B.png", width=12, height=10, dpi=600)

# Taxonomic Summary by Depth

#dep_meta<-unique(data.frame("Depth_m"=metadata$Depth_m, "Sample_Color"=metadata$Sample_Color))
#p_dep_meta<-merge(dep_meta,apr.gen.dep_m, by="Depth_m")
tg1<-ggplot(apr.gen.dep_m[apr.gen.dep_m$Count>0.01,], aes(Genus, Count)) +
  geom_jitter(aes(color=as.numeric(Depth_m)), size=3, width=0.15, height=0) +
  scale_colour_gradient2(low="red",high="blue3",midpoint=5,guide = guide_colourbar(reverse = TRUE)) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Microbial Genera & Depth - April 2022", subtitle="Includes taxa with Relative Abundance > 1%",color="Depth (m)")

ggsave(tg1,filename = "figures/RelativeAbundance/Genus/April2022/SSW_16S_Apr22_Genera.RA_depth_taxasum_1perc.png", width=15, height=10, dpi=600)

tg1a<-ggplot(apr.gen.dep_m[apr.gen.dep_m$Count>0.01,], aes(Genus, Count)) +
  geom_jitter(aes(color=as.numeric(Depth_m)), size=3, width=0.15, height=0) +
  scale_colour_gradient2(low="red",high="blue3",midpoint=5,guide = guide_colourbar(reverse = TRUE)) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Microbial Genera & Depth - April 2022", subtitle="Includes taxa with Relative Abundance > 1%",color="Depth (m)")+coord_flip()

ggsave(tg1a,filename = "figures/RelativeAbundance/Genus/April2022/SSW_16S_Apr22_Genera.RA_depth_taxasum_1perc_v2.png", width=15, height=10, dpi=600)

tg1b<-ggplot(apr.gen.dep_m[apr.gen.dep_m$Count>0.05,], aes(Genus, Count)) +
  geom_jitter(aes(color=as.numeric(Depth_m)), size=3, width=0.15, height=0) +
  scale_colour_gradient2(low="red",high="blue3",midpoint=5,guide = guide_colourbar(reverse = TRUE)) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Microbial Genera & Depth - April 2022", subtitle="Includes taxa with Relative Abundance > 5%",color="Depth (m)")

ggsave(tg1b,filename = "figures/RelativeAbundance/Genus/April2022/SSW_16S_Apr22_Genera.RA_depth_taxasum_5percent.png", width=15, height=10, dpi=600)

tg1c<-ggplot(apr.gen.dep_m[apr.gen.dep_m$Count>0.05,], aes(Genus, Count)) +
  geom_jitter(aes(color=as.numeric(Depth_m)), size=3, width=0.15, height=0) +
  scale_colour_gradient2(low="red",high="blue3",midpoint=5,guide = guide_colourbar(reverse = TRUE)) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Microbial Genera & Depth - April 2022", subtitle="Includes taxa with Relative Abundance > 5%",color="Depth (m)")+coord_flip()

ggsave(tg1c,filename = "figures/RelativeAbundance/Genus/April2022/SSW_16S_Apr22_Genera.RA_depth_taxasum_5percent_v2.png", width=15, height=10, dpi=600)



#### Looking at Most Abundant Genus Only - DS001 ####
bac.dat.all.g<-subset(bac.dat.all, bac.dat.all$Genus!="Unknown") # drop unknown genera so they don't skew analyses

# by Genus + Sampling Date + Depth
bac.gen.date.dep <- as.data.frame(dcast(bac.dat.all.g,SampDate+Depth_m~Genus, value.var="Count", fun.aggregate=sum)) ###
bac.gen.date.dep[1:5,1:5] # counts by Genus + sample date & depth
rownames(bac.gen.date.dep)<-interaction(bac.gen.date.dep$SampDate,bac.gen.date.dep$Depth_m,sep="_")
bac.gen.date.dep[1:5,1:5]

b.RA_gen.date.dep<-data.frame(decostand(bac.gen.date.dep[,-c(1:2)], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(b.RA_gen.date.dep) # sanity check
b.RA_gen.date.dep$SampDate_Depth<-rownames(b.RA_gen.date.dep)
b.RA_gen.date.dep[1:5,ncol(b.RA_gen.date.dep):5-ncol(b.RA_gen.date.dep)] # first 5 rows, last 5 columns

ds001.date.dep<-subset(b.RA_gen.date.dep, select=c(DS001, SampDate_Depth))

ds001.date.dep.2<-as.data.frame(separate_wider_delim(data = ds001.date.dep, col=SampDate_Depth, "_", names = c("SampDate", "Depth_m"))) # Separate SampDate & Depth column for Heatmap later
ds001.date.dep.2$SampDate_Depth<-interaction(ds001.date.dep.2$SampDate,ds001.date.dep.2$Depth_m)

ds001.date.dep.2$Depth_m<-factor(ds001.date.dep.2$Depth_m, levels=c("0","3","4","5","7","9","10","10.5"))
ds001.date.dep.2$SampDate<-factor(ds001.date.dep.2$SampDate,levels=c("August.2021","December.2021","April.2022"))
ds001.date.dep.2$SampDate_Depth = factor(ds001.date.dep.2$SampDate_Depth, levels=unique(ds001.date.dep.2$SampDate_Depth[order(ds001.date.dep.2$Depth_m,ds001.date.dep.2$SampDate)]), ordered=TRUE)

ds001.date.dep.2

ds001_meta<-merge(ds001.date.dep.2,meta_scaled, by=c("SampDate","Depth_m"))
rownames(ds001_meta)<-ds001_meta$SampleID
head(ds001_meta)

# is this genus normally distributed?
shapiro.test(ds001_meta$DS001) # what is the p-value?
# p-value = 0.002568
# p > 0.05 states distribution of data are not significantly different from normal distribution
# p < 0.05 means that data is significantly different from a normal distribution
hist(ds001_meta$DS001, col="blue3")

# visualize Q-Q plot for species richness
qqnorm(ds001_meta$DS001, pch = 1, frame = FALSE)
qqline(ds001_meta$DS001, col = "red", lwd = 2)

# do env variables & RelAb of DS001 correlate?
cor.test(ds001_meta$DS001, ds001_meta$ORP_mV, method="pearson")
cor.test(ds001_meta$DS001, ds001_meta$Dissolved_OrganicMatter_RFU, method="pearson") # r = 0.6476844, p-value =  0.0006222
cor.test(ds001_meta$DS001, ds001_meta$DO_Percent_Local, method="pearson")
cor.test(ds001_meta$DS001, ds001_meta$Temp_DegC, method="pearson")
cor.test(ds001_meta$DS001, ds001_meta$Depth.num, method="pearson")
cor.test(ds001_meta$DS001, ds001_meta$Sulfate_milliM, method="pearson") # r = 0.4249563, p-value = 0.03845
cor.test(ds001_meta$DS001, ds001_meta$Sulfide_microM, method="pearson")

# does RelAb of DS001 change with sample date??
# Kruskal-Wallis test = nonparametric one-way ANOVA
kruskal.test(DS001 ~ SampDate, data = ds001_meta)
#Kruskal-Wallis chi-squared = 15.38, df = 2, p-value = 0.0004574
pairwise.wilcox.test(ds001_meta$DS001, ds001_meta$SampDate, p.adjust.method = "bonf") # returns p values
#               August.2021 December.2021
#December.2021 1.00000     -
#April.2022    0.00047     0.00047

# Levene's test with one independent variable
## Levene's tests whether variances of 2 samples are equal
## we want variances to be the same -- want NON SIGNIFICANCE!
## t test assumes that variances are the same, so Levene's test needs to be non significant
fligner.test(DS001 ~ SampDate, data = ds001_meta)
# Levenes Test for Homogeneity of Variance
#  Fligner-Killeen:med chi-squared = 7.9189, df = 2, p-value = 0.01907
# Which shows that the data DO deviate significantly from homogeneity.
compare_means(DS001 ~ SampDate, data=ds001_meta, method="kruskal.test",p.adjust.method = "bonferroni") # not significant
# .y.       p p.adj p.format p.signif method
# <chr> <dbl> <dbl> <chr>    <chr>    <chr>
# DS001 0.000457 0.00046 0.00046  ***      Kruskal-Wallis
plot(DS001 ~ SampDate, data=ds001_meta)

ds001.date.ts1<-ggplot(ds001_meta, aes(SampDate, DS001)) +
  geom_jitter(aes(color=factor(SampDate)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Sample Date", values=c("#ef781c","#03045e","#059c3f"), labels=c("August 2021","December 2021","April 2022")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Sample Date", y="Relative Abundance", title="Genus DS001 & Sample Date")+
  scale_x_discrete(labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022"))

ggsave(ds001.date.ts1,filename = "figures/RelativeAbundance/Genus/DS001/SSW_16S_DS001_RA_bydate_taxasum.png", width=15, height=10, dpi=600)

ds001.date.ts2<-ggplot(ds001_meta, aes(SampDate, DS001)) +
  geom_jitter(aes(color=as.numeric(Depth_m)), size=2, width=0.15, height=0) +
  scale_colour_gradient2(low="red",high="blue3",midpoint=5,guide = guide_colourbar(reverse = TRUE)) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=14),axis.title.y = element_text(size=14),axis.text = element_text(size=12),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=14),legend.text = element_text(size=13),plot.title = element_text(size=16)) +
  labs(x="Sample Date", y="Relative Abundance", title="Genus DS001 & Sample Date", color="Depth (m)")+scale_x_discrete(labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022"))

ggsave(ds001.date.ts2,filename = "figures/RelativeAbundance/Genus/DS001/SSW_16S_DS001_RA_bydate_depth_taxasum.png", width=15, height=10, dpi=600)

# does RelAb of DS001 change with depth??
# Kruskal-Wallis test = nonparametric one-way ANOVA
kruskal.test(DS001 ~ Depth_m, data = ds001_meta)
#Kruskal-Wallis chi-squared = 3.64, df = 7, p-value = 0.8202
pairwise.wilcox.test(ds001_meta$DS001, ds001_meta$Depth_m, p.adjust.method = "bonf") # returns p values

# Levene's test with one independent variable
## Levene's tests whether variances of 2 samples are equal
## we want variances to be the same -- want NON SIGNIFICANCE!
## t test assumes that variances are the same, so Levene's test needs to be non significant
fligner.test(DS001 ~ Depth_m, data = ds001_meta)
# Levenes Test for Homogeneity of Variance
#  Fligner-Killeen:med chi-squared = 1.3502, df = 7, p-value = 0.987
# Which shows that the data DO deviate significantly from homogeneity.
compare_means(DS001 ~ Depth_m, data=ds001_meta, method="kruskal.test",p.adjust.method = "bonferroni") # not significant
# .y.       p p.adj p.format p.signif method
# <chr> <dbl> <dbl> <chr>    <chr>    <chr>
#   1 DS001 0.820  0.82 0.82     ns       Kruskal-Wallis
plot(DS001 ~ SampDate, data=ds001_meta)

# ds001.dep.ts1<-ggplot(ds001_meta, aes(Depth_m, DS001)) +
#   geom_jitter(aes(color=as.numeric(Depth_m)), size=2, width=0.15, height=0) +
#   scale_colour_gradient2(low="red",high="blue3",midpoint=5,guide = guide_colourbar(reverse = TRUE)) +
#   geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
#   labs(x="Depth (m)", y="Relative Abundance", title="Genus DS001 & Depth", color="Depth (m)")
#
# ggsave(ds001.dep.ts1,filename = "figures/RelativeAbundance/Genus/DS001/SSW_16S_DS001_RA_bydepth_taxasum.png", width=15, height=10, dpi=600)

ds001.dep.ts2<-ggplot(ds001_meta, aes(Depth_m, DS001)) +
  geom_jitter(aes(color=factor(SampDate)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Sample Date", values=c("#ef781c","#03045e","#059c3f"), labels=c("August 2021","December 2021","April 2022")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Depth (m)", y="Relative Abundance", title="Genus DS001 by Depth & Sample Date", color="Sample Date")

ggsave(ds001.dep.ts2,filename = "figures/RelativeAbundance/Genus/DS001/SSW_16S_DS001_RA_bydepth_date_taxasum.png", width=15, height=10, dpi=600)

ggplot(metadata, aes(x=Depth_m, y=Sulfate_milliM,color=SampDate,group=SampDate)) + geom_point(size=3) + geom_line()

ds001.dep.ts3<-ggplot(ds001_meta, aes(x=Depth_m,y=DS001,color=SampDate,group=SampDate)) +
  geom_point(size=3) + geom_line() +
  scale_color_manual(name ="Sample Date", values=c("#ef781c","#03045e","#059c3f"), labels=c("August 2021","December 2021","April 2022")) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Depth (m)", y="Relative Abundance", title="Genus DS001 by Depth & Sample Date", color="Sample Date")+coord_flip()+ scale_x_discrete(limits=rev)

ggsave(ds001.dep.ts3,filename = "figures/RelativeAbundance/Genus/DS001/SSW_16S_DS001_RA_bydepth_date_scatterplot.png", width=15, height=10, dpi=600)

# can we use env variables to predict RelAb of DS001?!
## correlation gives us the status of their relationship, but linear regression will tell us if these env variables predict DS001 relative abundance
## more on that here https://www.researchgate.net/post/Two_variables_are_correlated_but_regression_is_not_significant#:~:text=The%20simple%20answer%20is%20yes,opposite%20direction%20(negative%20correlation).

# first lets make a big model and see what we find
step.ds001<-step(glm(formula = DS001 ~ ., data=ds001_meta[,c(3,10,12:13,17:19,21)]))
summary(step.ds001)
#                             Estimate Std. Error t value Pr(>|t|)
# (Intercept)                 0.235790   0.042692   5.523 3.72e-05 ***
#   DO_Percent_Local            0.112263   0.043520   2.580   0.0195 *
#   ORP_mV                      0.052997   0.019824   2.673   0.0160 *
#   Temp_DegC                   0.089822   0.038534   2.331   0.0323 *
#   Dissolved_OrganicMatter_RFU 0.147940   0.025322   5.842 1.96e-05 ***
#   Sulfate_milliM              0.063009   0.022376   2.816   0.0119 *
#   Depth.num                   0.011092   0.006612   1.677   0.1117

# Remove sulfide & rerun
ds001.dom.fit1<-glm(DS001 ~ Dissolved_OrganicMatter_RFU+DO_Percent_Local+ORP_mV+Temp_DegC+Sulfate_milliM+Depth.num, data=ds001_meta) %>%
  adjust_pvalue(method="bonferroni")
summary(ds001.dom.fit1)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)                 0.235790   0.042692   5.523 3.72e-05 ***
# Dissolved_OrganicMatter_RFU 0.147940   0.025322   5.842 1.96e-05 ***
#   DO_Percent_Local            0.112263   0.043520   2.580   0.0195 *
#   ORP_mV                      0.052997   0.019824   2.673   0.0160 *
#   Temp_DegC                   0.089822   0.038534   2.331   0.0323 *
#   Sulfate_milliM              0.063009   0.022376   2.816   0.0119 *
#   Depth.num                   0.011092   0.006612   1.677   0.1117

anova(ds001.dom.fit1, test = "LR")
# Terms added sequentially (first to last)
#                             Df Deviance Resid. Df Resid. Dev  Pr(>Chi)
# NULL                                           23    0.45248
# Dissolved_OrganicMatter_RFU  1 0.189814        22    0.26267 4.466e-11 ***
#   DO_Percent_Local             1 0.128567        21    0.13410 5.903e-08 ***
#   ORP_mV                       1 0.013252        20    0.12085  0.081749 .
# Temp_DegC                    1 0.003926        19    0.11692  0.343415
# Sulfate_milliM               1 0.030262        18    0.08666  0.008529 **
#   Depth.num                    1 0.012307        17    0.07435  0.093450 .
coef(summary(ds001.dom.fit1))
#p.adjust(coef(summary(ds001.dom.fit1))[,4], method="bonferroni") # pvalue

autoplot(ds001.dom.fit1, which = 1:6, label.size = 3,colour='SampDate',data=ds001_meta)

#dispersiontest(ds001.dom.fit1)
# null hypothesis is that equidispersion exists; alternative hypothesis is overdispersion
# if overdispersion, use negative binomial not Poisson
## Poisson distribution implies that the mean and variance are equal --> little dispersion
# negative binomial means # of observations is not fixed, whereas binomial means observations are a fixed #

ds001.dom.fit2<-glm(DS001 ~ Dissolved_OrganicMatter_RFU+DO_Percent_Local+Sulfate_milliM, data=ds001_meta) %>%
  adjust_pvalue(method="bonferroni")
summary(ds001.dom.fit2)

anova(ds001.dom.fit2, test = "LR")

ds001.dom.fit3<-glm(DS001 ~ Dissolved_OrganicMatter_RFU+DO_Percent_Local, data=ds001_meta) %>%
  adjust_pvalue(method="bonferroni")
summary(ds001.dom.fit3)

anova(ds001.dom.fit3, test = "LR")

ds001.dom.fit4<-glm(DS001 ~ Dissolved_OrganicMatter_RFU, data=ds001_meta) %>%
  adjust_pvalue(method="bonferroni")
summary(ds001.dom.fit4)

anova(ds001.dom.fit4, test = "LR")

plot(DS001 ~ Dissolved_OrganicMatter_RFU, data=ds001_meta,col=SampDate_Color)
# June 2021 = green, orange = August 2021, dark blue = December 2021, light blue = April 2022

plot(DS001 ~ DO_Percent_Local, data=ds001_meta,col=SampDate_Color)
plot(DS001 ~ Sulfate_milliM, data=ds001_meta,col=SampDate_Color)

# Plot DS001 RelAb x DOM

fig.ds001.dom.fit1<-ggplot(ds001_meta, aes(x = Dissolved_OrganicMatter_RFU, y = DS001)) +
  geom_point(aes(color=SampDate), size=3) + theme_classic() +
  scale_color_manual(name ="Sample Date", values=c("#ef781c","#03045e","#059c3f"), labels=c("August 2021","December 2021","April 2022")) +
  stat_smooth(method = "glm", col = "black", se=FALSE, linewidth=1)+ labs(title="Relative Abundance of Genus DS001 vs. Dissolved Organic Matter", color="Sampling Date")+ylab("Relative Abundance")+xlab("Dissolved Organic Matter (RFU)")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  stat_cor(label.y = 0.53, label.x=0.8) +
  stat_regline_equation(aes(label=paste(after_stat(adj.rr.label))),label.y = 0.55,label.x=0.8)

## use summary(ds001.so4.fit1) to double check that stat_cor gives same p value as linear regression! it does here :)
ggsave(fig.ds001.dom.fit1,filename = "figures/RelativeAbundance/Genus/DS001/DS001_Genus_RelAb_vs_DOM_scatterplot.pdf", width=10, height=8, dpi=600)

fig.ds001.dom.fit1a<-ggplot(ds001_meta, aes(x = Dissolved_OrganicMatter_RFU, y = DS001)) +
  geom_point(aes(color=SampDate), size=3) + theme_classic() +
  scale_color_manual(name ="Sample Date", values=c("#ef781c","#03045e","#059c3f"), labels=c("August 2021","December 2021","April 2022")) + labs(title="Relative Abundance of Genus DS001 vs. Dissolved Organic Matter", color="Sampling Date")+ylab("Relative Abundance")+xlab("Dissolved Organic Matter (RFU)")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  stat_cor(label.y = 0.53, label.x=0.8) +
  stat_regline_equation(aes(label=paste(after_stat(adj.rr.label))),label.y = 0.55,label.x=0.8)

## use summary(ds001.so4.fit1) to double check that stat_cor gives same p value as linear regression! it does here :)
ggsave(fig.ds001.dom.fit1a,filename = "figures/RelativeAbundance/Genus/DS001/DS001_Genus_RelAb_vs_DOM_scatterplot_noline.pdf", width=10, height=8, dpi=600)

## next Sulfate

ds001.so4.fit1<-glm(DS001 ~ Sulfate_milliM, data=ds001_meta) %>%
  adjust_pvalue(method="bonferroni")
summary(ds001.so4.fit1)
#Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)
#(Intercept)     0.26851    0.01905  14.097  < 2e-16 ***
#  Sulfate_milliM  0.05906    0.01929   3.061  0.00403 **

coef(summary(ds001.so4.fit1))
p.adjust(coef(summary(ds001.so4.fit1))[,4], method="bonferroni") # pvalue

plot(DS001 ~ Sulfate_milliM, data=ds001_meta,col=SampDate_Color)
# orange = August 2021, dark blue = December 2021, light blue = April 2022

# Plot RelAb of DS001 x Sulfate

fig.ds001.so4.fit1<-ggplot(ds001_meta, aes(x = Sulfate_milliM, y = DS001)) +
  geom_point(aes(color=SampDate), size=3) + theme_classic() +
  scale_color_manual(name ="Sample Date", values=c("#ef781c","#03045e","#059c3f"), labels=c("August 2021","December 2021","April 2022")) +
  stat_smooth(method = "glm", col = "black", se=FALSE, size=1)+ labs(title="Relative Abundance of Genus DS001 vs. Sulfate (milliM)", color="Sampling Date")+ylab("Relative Abundance")+xlab("Sulfate (milliM)")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  stat_cor(label.y = 0.5, label.x=0.6) +
  stat_regline_equation(aes(label=paste(after_stat(adj.rr.label))),label.y = 0.53,label.x=0.6)

## use summary(ds001.so4.fit1) to double check that stat_cor gives same p value as linear regression! it does here :)
ggsave(fig.ds001.so4.fit1,filename = "figures/RelativeAbundance/Genus/DS001/DS001_Genus_RelAb_vs_Sulfate_scatterplot.pdf", width=10, height=8, dpi=600)

#### Shared Genus Relative Abundance ####
head(bac.dat.all)

head(b.g.typ_m) # made in section above, ~ line 1277

# merge metadata and RA data
typ_meta<-unique(data.frame("Sample_Type"=metadata$Sample_Type, "Sample_Color"=metadata$Sample_Color))
g_typ_meta<-merge(typ_meta,b.g.typ_m, by="Sample_Type")
g_typ_meta<-subset(g_typ_meta, Sample_Type=="Dust" | Sample_Type=="Soil" | Sample_Type=="Seawater")
g_typ_meta<-subset(g_typ_meta, Genus!="Unknown")
g_typ_meta<-subset(g_typ_meta, Count!=0)

# finding shared genera...

n_occur <- data.frame(table(g_typ_meta$Genus)) # find frequencies of genera to see which are shared between sample types
n_occur[n_occur$Freq > 1,] # shows us which genera have a greater frequency than 2
g_shared_typ<-g_typ_meta[g_typ_meta$Genus %in% n_occur$Var1[n_occur$Freq > 1],]
#write.csv(g_shared_typ,"16S_Genera_SampleType_Shared.csv",row.names=FALSE)

g_not.shared_typ<-subset(g_typ_meta, !(g_typ_meta$Genus %in% g_shared_typ$Genus)) # subset based off of what is NOT in one dataframe from another data frame

ggplot(g_shared_typ[g_shared_typ$Count>=0.001,], aes(x=Sample_Type, y=Count, fill=Genus))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genera Relative Abundance", x="SampleID", y="Relative Abundance", fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+guides(fill=guide_legend(ncol=2))+scale_y_continuous(expand = c(0,0),limits = c(0,1))

sh.t1<-ggplot(g_shared_typ, aes(Genus, Count)) +
  geom_jitter(aes(color=factor(Sample_Type)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Sample Type", values=unique(g_shared_typ$Sample_Color[order(g_shared_typ$Sample_Type)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13, vjust=-1),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genus", y="Relative Abundance", title="Shared Microbial Genera Across Sample Types")

ggsave(sh.t1,filename = "figures/RelativeAbundance/16S_shared_Genera_taxasum_type.png", width=23, height=10, dpi=600)

sh.t1a<-ggplot(g_shared_typ[g_shared_typ$Count>=0.01,], aes(Genus, Count)) +
  geom_jitter(aes(color=factor(Sample_Type)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Sample Type", values=unique(g_shared_typ$Sample_Color[order(g_shared_typ$Sample_Type)])) + theme_classic() +
  theme(axis.title.x = element_text(size=13, vjust=-1),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genus", y="Relative Abundance", title="Microbial Genera & Sample Type",subtitle="Only Includes Genera Shared Across All Sample Types")

ggsave(sh.t1a,filename = "figures/RelativeAbundance/16S_shared_Genera_taxasum_type_1perc.png", width=23, height=10, dpi=600)

ggplot(g_shared_typ[g_shared_typ$Count>=0.01,], aes(Genus, Count)) +
  geom_jitter(aes(color=factor(Sample_Type)), size=3, width=0.15, height=0) +
  scale_color_manual(name ="Sample Type", values=unique(g_shared_typ$Sample_Color[order(g_shared_typ$Sample_Type)])) + theme_classic() +
  geom_boxplot(fill=NA, outlier.color=NA) +
  theme(axis.title.x = element_text(size=13, vjust=-1),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genus", y="Relative Abundance", title="Microbial Genera & Sample Type")

ggplot(g_shared_typ[g_shared_typ$Count>=0.0005,], aes(Genus, Count)) +
  geom_jitter(aes(color=factor(Sample_Type)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Sample Type", values=unique(g_shared_typ$Sample_Color[order(g_shared_typ$Sample_Type)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13, vjust=-1),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genus", y="Relative Abundance", title="Microbial Genera & Sample Type")

sh.t2<-ggplot(g_shared_typ[g_shared_typ$Count>=0.005,], aes(Genus, Count)) +
  geom_jitter(aes(color=factor(Sample_Type)), size=3, width=0.15, height=0) +
  scale_color_manual(name ="Sample Type", values=unique(g_shared_typ$Sample_Color[order(g_shared_typ$Sample_Type)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13, vjust=-1),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genus", y="Relative Abundance", title="Microbial Genera & Sample Type",subtitle="Relative Abundance > 0.25%")

ggsave(sh.t2,filename = "figures/RelativeAbundance/16S_shared_Genera_taxasum_type_0.5perc.png", width=23, height=10, dpi=600)

sh.t3<-ggplot(g_shared_typ[g_shared_typ$Count>=0.0025,], aes(Genus, Count)) +
  geom_jitter(aes(color=factor(Sample_Type)), size=3, width=0.15, height=0) +
  scale_color_manual(name ="Sample Type", values=unique(g_shared_typ$Sample_Color[order(g_shared_typ$Sample_Type)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13, vjust=-1),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genus", y="Relative Abundance", title="Microbial Genera & Sample Type",subtitle="Relative Abundance > 0.25%")

ggsave(sh.t3,filename = "figures/RelativeAbundance/16S_shared_Genera_taxasum_type_0.25perc.png", width=23, height=10, dpi=600)

## Comparing Genera in Dust vs Seawater
# X Axis Breaks and Labels
lbls = paste0(as.character(c(seq(0.05, 0, -0.01), seq(0.01, 0.05, 0.01)))) # labels
brks=seq(-0.05,0.05,0.01)
g_shared_typ1<-subset(g_shared_typ, Sample_Type!="Soil")
g_shared_typ1$Count2 <- ifelse(g_shared_typ1$Sample_Type == "Seawater", -1*g_shared_typ1$Count, g_shared_typ1$Count)
g_shared_typ1<-g_shared_typ1[order(-g_shared_typ1$Count2,g_shared_typ1$Genus),]
g_shared_typ1$GenSamp<-interaction(g_shared_typ1$Genus,g_shared_typ1$Sample_Type)
g_shared_typ1$GenSamp<-factor(g_shared_typ1$GenSamp, levels=g_shared_typ1$GenSamp)
class(g_shared_typ1$GenSamp)
g_shared_typ1$Genus<-factor(g_shared_typ1$Genus, levels=unique(g_shared_typ1$Genus[sort(g_shared_typ1$GenSamp)]))

share1<-ggplot(g_shared_typ1, aes(x = Genus, y = -Count2, fill = Sample_Type)) +
  geom_bar(data = subset(g_shared_typ1[g_shared_typ1$Count2<=-0.0005,], Sample_Type == "Seawater"), stat = "identity") +
  geom_bar(data = subset(g_shared_typ1[g_shared_typ1$Count2>=0.0005,], Sample_Type == "Dust"), stat = "identity") +
  coord_flip()+scale_y_continuous(labels = lbls,breaks=brks)+theme_classic()+scale_fill_manual(name ="Sample Type", values=unique(g_shared_typ1$Sample_Color[rev(order(g_shared_typ1$Sample_Type))]))+ylab("Relative Abundance")+
  theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+labs(title="Microbial Genera by Sample Type", subtitle="Includes Taxa with a Relative Abundance of at least 0.05%")

ggsave(share1,filename = "figures/RelativeAbundance/16S_shared_Genera_dust.v.sw_population.pyramid.png", width=12, height=10, dpi=600)

pp2<-ggplot(g_shared_typ1[g_shared_typ1$Count>=0.0005,], aes(x = reorder(Genus,Count), fill = Sample_Type,y = ifelse(test = Sample_Type == "Dust",yes = Count, no = -Count))) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = abs, limits = max(g_shared_typ1$Count) * c(-1,1)) +
  coord_flip()+scale_fill_manual(name ="Sample Type", values=unique(g_shared_typ1$Sample_Color[order(g_shared_typ1$Sample_Type)]))+ylab("Relative Abundance")+theme_classic()+
  theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
  labs(title="Microbial Genera by Sample Type",subtitle="Includes Only Shared Taxa with a Relative Abundance of at least 0.05%")+xlab("Genus")

ggsave(pp2,filename = "figures/RelativeAbundance/16S_shared_Genera_dust.v.sw_population.pyramid.pretty.05percent.png", width=13, height=10, dpi=600)

pp3<-ggplot(g_shared_typ1[g_shared_typ1$Count>=0.0010,], aes(x = reorder(Genus,Count), fill = Sample_Type,y = ifelse(test = Sample_Type == "Dust",yes = Count, no = -Count))) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = abs, limits = max(g_shared_typ1$Count) * c(-1,1)) +
  coord_flip()+scale_fill_manual(name ="Sample Type", values=unique(g_shared_typ1$Sample_Color[order(g_shared_typ1$Sample_Type)]))+ylab("Relative Abundance")+theme_classic()+
  theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
  labs(title="Microbial Genera by Sample Type",subtitle="Includes Only Shared Taxa with a Relative Abundance of at least 0.1%")+xlab("Genus")

ggsave(pp3,filename = "figures/RelativeAbundance/16S_shared_Genera_dust.v.sw_population.pyramid.pretty.1percent.png", width=13, height=10, dpi=600)

pp4<-ggplot(g_shared_typ1[g_shared_typ1$Count>=0.010,], aes(x = reorder(Genus,Count), fill = Sample_Type,y = ifelse(test = Sample_Type == "Dust",yes = Count, no = -Count))) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = abs, limits = max(g_shared_typ1$Count) * c(-1,1)) +
  coord_flip()+scale_fill_manual(name ="Sample Type", values=unique(g_shared_typ1$Sample_Color[order(g_shared_typ1$Sample_Type)]))+ylab("Relative Abundance")+theme_classic()+
  theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
  labs(title="Microbial Genera by Sample Type",subtitle="Includes Only Shared Taxa with a Relative Abundance of at least 1%")+xlab("Genus")
ggsave(pp4,filename = "figures/RelativeAbundance/16S_shared_Genera_dust.v.sw_population.pyramid.pretty.1percent_10.14.22.png", width=13, height=10, dpi=600)

## Lollipop chart

lg1<-ggplot(g_shared_typ1[g_shared_typ1$Count>=0.0005,], aes(x = reorder(Genus,Count),
                                                              y = ifelse(test = Sample_Type == "Dust",yes = Count, no = -Count),color=g_shared_typ1$Sample_Type[g_shared_typ1$Count>=0.0005])) +
  geom_point(stat='identity',size=3)  +
  geom_segment(aes(y = 0,
                   x = Genus,
                   yend = ifelse(test = Sample_Type == "Dust",yes = Count, no = -Count),
                   xend = Genus),color = "black") +
  coord_flip()+scale_color_manual(name ="Sample Type", values=unique(g_shared_typ1$Sample_Color[order(g_shared_typ1$Sample_Type)]))+ylab("Relative Abundance")+theme_classic()+theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
  labs(title="Microbial Genera by Sample Type")+xlab("Genus")

ggsave(lg1,filename = "figures/RelativeAbundance/16S_shared_Genera_dust.v.sw_lollipop_chart.png", width=13, height=10, dpi=600)

ggplot(g_shared_typ1[g_shared_typ1$Count>=0.0005,], aes(x = reorder(Genus,Count),
                                                         y = ifelse(test = Sample_Type == "Dust",yes = Count, no = -Count),color=g_shared_typ1$Sample_Type[g_shared_typ1$Count>=0.0005])) +
  geom_bar(stat = "identity") +
  coord_flip()+scale_color_manual(name ="Sample Type", values=unique(g_shared_typ1$Sample_Color[order(g_shared_typ1$Sample_Type)]))+ylab("Relative Abundance")+theme_classic()+theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
  labs(title="Microbial Genera by Sample Type")+xlab("Genus")

## Dust vs Soil
# X Axis Breaks and Labels
lbls = paste0(as.character(c(seq(0.05, 0, -0.01), seq(0.01, 0.05, 0.01)))) # labels
brks=seq(-0.05,0.05,0.01)
g_shared_typ2<-subset(g_shared_typ, Sample_Type!="Seawater")
g_shared_typ2$Count2 <- ifelse(g_shared_typ2$Sample_Type == "Soil", -1*g_shared_typ2$Count, g_shared_typ2$Count)
g_shared_typ2<-g_shared_typ2[order(-g_shared_typ2$Count2,g_shared_typ2$Genus),]
g_shared_typ2$GenSamp<-interaction(g_shared_typ2$Genus,g_shared_typ2$Sample_Type)
g_shared_typ2$GenSamp<-factor(g_shared_typ2$GenSamp, levels=g_shared_typ2$GenSamp)
class(g_shared_typ2$GenSamp)
g_shared_typ2$Genus<-factor(g_shared_typ2$Genus, levels=unique(g_shared_typ2$Genus[sort(g_shared_typ2$GenSamp)]))

share2<-ggplot(g_shared_typ2, aes(x = Genus, y = -Count2, fill = Sample_Type)) +
  geom_bar(data = subset(g_shared_typ2[g_shared_typ2$Count2<=-0.0005,], Sample_Type == "Soil"), stat = "identity") +
  geom_bar(data = subset(g_shared_typ2[g_shared_typ2$Count2>=0.0005,], Sample_Type == "Dust"), stat = "identity") +
  coord_flip()+scale_y_continuous(labels = lbls,breaks=brks)+theme_classic()+scale_fill_manual(name ="Sample Type", values=unique(g_shared_typ2$Sample_Color[order(g_shared_typ2$Sample_Type)]))+ylab("Relative Abundance")+
  theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+labs(title="Microbial Genera by Sample Type", subtitle="Includes Taxa with a Relative Abundance of at least 0.05%")

ggsave(share2,filename = "figures/RelativeAbundance/16S_shared_Genera_dust.v.soil_population.pyramid.png", width=12, height=10, dpi=600)


ggplot(g_shared_typ2, aes(x = Genus, y = Count2, fill = Sample_Type)) +
  geom_bar(data = subset(g_shared_typ2, Sample_Type == "Dust"), stat = "identity") +
  geom_bar(data = subset(g_shared_typ2, Sample_Type == "Soil"), stat = "identity") +
  coord_flip()+scale_y_continuous(labels = lbls,breaks=brks)+theme_classic()+scale_fill_manual(name ="Sample Type", values=unique(g_shared_typ2$Sample_Color[order(g_shared_typ2$Sample_Type)]))+ylab("Relative Abundance")+
  theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
  labs(title="Microbial Genera by Sample Type", subtitle="Includes Taxa with a Relative Abundance of at least 0.05%")

pp2a<-ggplot(g_shared_typ2[g_shared_typ2$Count>=0.0005,], aes(x = reorder(Genus,Count), fill = Sample_Type,y = ifelse(test = Sample_Type == "Dust",yes = Count, no = -Count))) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = abs, limits = max(g_shared_typ2$Count) * c(-1,1)) +
  coord_flip()+scale_fill_manual(name ="Sample Type", values=unique(g_shared_typ2$Sample_Color[order(g_shared_typ2$Sample_Type)]))+ylab("Relative Abundance")+theme_classic()+
  theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
  labs(title="Microbial Genera by Sample Type",subtitle="Includes Only Shared Taxa with a Relative Abundance of at least 0.05%")+xlab("Genus")

ggsave(pp2a,filename = "figures/RelativeAbundance/16S_shared_Genera_dust.v.soil_population.pyramid.pretty.05percent.png", width=13, height=10, dpi=600)

pp3a<-ggplot(g_shared_typ2[g_shared_typ2$Count>=0.010,], aes(x = reorder(Genus,Count), fill = Sample_Type,y = ifelse(test = Sample_Type == "Dust",yes = Count, no = -Count))) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = abs, limits = max(g_shared_typ2$Count) * c(-1,1)) +
  coord_flip()+scale_fill_manual(name ="Sample Type", values=unique(g_shared_typ2$Sample_Color[order(g_shared_typ2$Sample_Type)]))+ylab("Relative Abundance")+theme_classic()+
  theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
  labs(title="Microbial Genera by Sample Type",subtitle="Includes Only Shared Taxa with a Relative Abundance of at least 0.1%")+xlab("Genus")

ggsave(pp3a,filename = "figures/RelativeAbundance/16S_shared_Genera_dust.v.soil_population.pyramid.pretty.1percent.png", width=13, height=10, dpi=600)

pp4a<-ggplot(g_shared_typ2[g_shared_typ2$Count>=0.010,], aes(x = reorder(Genus,Count), fill = Sample_Type,y = ifelse(test = Sample_Type == "Dust",yes = Count, no = -Count))) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = abs, limits = max(g_shared_typ2$Count) * c(-1,1)) +
  coord_flip()+scale_fill_manual(name ="Sample Type", values=unique(g_shared_typ2$Sample_Color[order(g_shared_typ2$Sample_Type)]))+ylab("Relative Abundance")+theme_classic()+
  theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
  labs(title="Microbial Genera by Sample Type",subtitle="Includes Only Shared Taxa with a Relative Abundance of at least 1%")+xlab("Genus")

ggsave(pp4a,filename = "figures/RelativeAbundance/16S_shared_Genera_dust.v.soil_population.pyramid.pretty_1percent_10.14.22.png", width=13, height=10, dpi=600)

## Lollipop chart

lg1a<-ggplot(g_shared_typ2[g_shared_typ2$Count>=0.0005,], aes(x = reorder(Genus,Count),
                                                               y = ifelse(test = Sample_Type == "Dust",yes = Count, no = -Count),color=g_shared_typ2$Sample_Type[g_shared_typ2$Count>=0.0005])) +
  geom_point(stat='identity',size=3)  +
  geom_segment(aes(y = 0,
                   x = Genus,
                   yend = ifelse(test = Sample_Type == "Dust",yes = Count, no = -Count),
                   xend = Genus),color = "black") +
  coord_flip()+scale_color_manual(name ="Sample Type", values=unique(g_shared_typ2$Sample_Color[order(g_shared_typ2$Sample_Type)]))+ylab("Relative Abundance")+theme_classic()+theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+labs(title="Microbial Genera by Sample Type")+xlab("Genus")

ggsave(lg1a,filename = "figures/RelativeAbundance/16S_shared_Genera_dust.v.soil_lollipop_chart.png", width=13, height=10, dpi=600)

## Seawater vs Soil
# X Axis Breaks and Labels
lbls = paste0(as.character(c(seq(0.05, 0, -0.01), seq(0.01, 0.05, 0.01)))) # labels
brks=seq(-0.05,0.05,0.01)
g_shared_typ3<-subset(g_shared_typ, Sample_Type!="Dust")
g_shared_typ3$Count2 <- ifelse(g_shared_typ3$Sample_Type == "Soil", -1*g_shared_typ3$Count, g_shared_typ3$Count)
g_shared_typ3<-g_shared_typ3[order(-g_shared_typ3$Count2,g_shared_typ3$Genus),]
g_shared_typ3$GenSamp<-interaction(g_shared_typ3$Genus,g_shared_typ3$Sample_Type)
g_shared_typ3$GenSamp<-factor(g_shared_typ3$GenSamp, levels=g_shared_typ3$GenSamp)
class(g_shared_typ3$GenSamp)
g_shared_typ3$Genus<-factor(g_shared_typ3$Genus, levels=unique(g_shared_typ3$Genus[sort(g_shared_typ3$GenSamp)]))

share3<-ggplot(g_shared_typ3, aes(x = Genus, y = -Count2, fill = Sample_Type)) +
  geom_bar(data = subset(g_shared_typ3[g_shared_typ3$Count2<=-0.0005,], Sample_Type == "Soil"), stat = "identity") +
  geom_bar(data = subset(g_shared_typ3[g_shared_typ3$Count2>=0.0005,], Sample_Type == "Seawater"), stat = "identity") +
  coord_flip()+scale_y_continuous(labels = lbls,breaks=brks)+theme_classic()+scale_fill_manual(name ="Sample Type", values=unique(g_shared_typ3$Sample_Color[order(g_shared_typ3$Sample_Type)]))+ylab("Relative Abundance")+
  theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+labs(title="Microbial Genera by Sample Type", subtitle="Includes Taxa with a Relative Abundance of at least 0.05%")

ggsave(share3,filename = "figures/RelativeAbundance/16S_shared_Genera_sw.v.soil_population.pyramid.png", width=12, height=10, dpi=600)


ggplot(g_shared_typ3, aes(x = Genus, y = Count2, fill = Sample_Type)) +
  geom_bar(data = subset(g_shared_typ3, Sample_Type == "Seawater"), stat = "identity") +
  geom_bar(data = subset(g_shared_typ3, Sample_Type == "Soil"), stat = "identity") +
  coord_flip()+scale_y_continuous(labels = lbls,breaks=brks)+theme_classic()+scale_fill_manual(name ="Sample Type", values=unique(g_shared_typ3$Sample_Color[order(g_shared_typ3$Sample_Type)]))+ylab("Relative Abundance")+
  theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
  labs(title="Microbial Genera by Sample Type", subtitle="Includes Taxa with a Relative Abundance of at least 0.05%")

pp2b<-ggplot(g_shared_typ3[g_shared_typ3$Count>=0.0005,], aes(x = reorder(Genus,Count), fill = Sample_Type,y = ifelse(test = Sample_Type == "Seawater",yes = Count, no = -Count))) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = abs, limits = max(g_shared_typ3$Count) * c(-1,1)) +
  coord_flip()+scale_fill_manual(name ="Sample Type", values=unique(g_shared_typ3$Sample_Color[order(g_shared_typ3$Sample_Type)]))+ylab("Relative Abundance")+theme_classic()+
  theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
  labs(title="Microbial Genera by Sample Type",subtitle="Includes Only Shared Taxa with a Relative Abundance of at least 0.05%")+xlab("Genus")

ggsave(pp2b,filename = "figures/RelativeAbundance/16S_shared_Genera_sw.v.soil_population.pyramid.pretty.05percent.png", width=13, height=10, dpi=600)

pp3b<-ggplot(g_shared_typ3[g_shared_typ3$Count>=0.0010,], aes(x = reorder(Genus,Count), fill = Sample_Type,y = ifelse(test = Sample_Type == "Seawater",yes = Count, no = -Count))) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = abs, limits = max(g_shared_typ3$Count) * c(-1,1)) +
  coord_flip()+scale_fill_manual(name ="Sample Type", values=unique(g_shared_typ3$Sample_Color[order(g_shared_typ3$Sample_Type)]))+ylab("Relative Abundance")+theme_classic()+
  theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
  labs(title="Microbial Genera by Sample Type",subtitle="Includes Only Shared Taxa with a Relative Abundance of at least 0.1%")+xlab("Genus")

ggsave(pp3b,filename = "figures/RelativeAbundance/16S_shared_Genera_sw.v.soil_population.pyramid.pretty.1percent.png", width=13, height=10, dpi=600)

pp4b<-ggplot(g_shared_typ3[g_shared_typ3$Count>=0.005,], aes(x = reorder(Genus,Count), fill = Sample_Type,y = ifelse(test = Sample_Type == "Seawater",yes = Count, no = -Count))) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = abs, limits = max(g_shared_typ3$Count) * c(-1,1)) +
  coord_flip()+scale_fill_manual(name ="Sample Type", values=unique(g_shared_typ3$Sample_Color[order(g_shared_typ3$Sample_Type)]))+ylab("Relative Abundance")+theme_classic()+
  theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
  labs(title="Microbial Genera by Sample Type",subtitle="Includes Only Shared Taxa with a Relative Abundance of at least 0.5%")+xlab("Genus")

ggsave(pp4b,filename = "figures/RelativeAbundance/16S_shared_Genera_sw.v.soil_population.pyramid.pretty.1percent.png", width=13, height=10, dpi=600)

pp5b<-ggplot(g_shared_typ3[g_shared_typ3$Count>=0.010,], aes(x = reorder(Genus,Count), fill = Sample_Type,y = ifelse(test = Sample_Type == "Seawater",yes = Count, no = -Count))) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = abs, limits = max(g_shared_typ3$Count) * c(-1,1)) +
  coord_flip()+scale_fill_manual(name ="Sample Type", values=unique(g_shared_typ3$Sample_Color[order(g_shared_typ3$Sample_Type)]))+ylab("Relative Abundance")+theme_classic()+
  theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
  labs(title="Microbial Genera by Sample Type",subtitle="Includes Only Shared Taxa with a Relative Abundance of at least 1%")+xlab("Genus")

ggsave(pp5b,filename = "figures/RelativeAbundance/16S_shared_Genera_sw.v.soil_population.pyramid.pretty_1percent.png", width=13, height=10, dpi=600)

## Lollipop chart

lg1b<-ggplot(g_shared_typ3[g_shared_typ3$Count>=0.0005,], aes(x = reorder(Genus,Count),
                                                               y = ifelse(test = Sample_Type == "Seawater",yes = Count, no = -Count),color=g_shared_typ3$Sample_Type[g_shared_typ3$Count>=0.0005])) +
  geom_point(stat='identity',size=3)  +
  geom_segment(aes(y = 0,
                   x = Genus,
                   yend = ifelse(test = Sample_Type == "Seawater",yes = Count, no = -Count),
                   xend = Genus),color = "black") +
  coord_flip()+scale_color_manual(name ="Sample Type", values=unique(g_shared_typ3$Sample_Color[order(g_shared_typ3$Sample_Type)]))+ylab("Relative Abundance")+theme_classic()+theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+labs(title="Microbial Genera by Sample Type")+xlab("Genus")

ggsave(lg1a,filename = "figures/RelativeAbundance/16S_shared_Genera_sw.v.soil_lollipop_chart.png", width=13, height=10, dpi=600)

typ_meta<-unique(data.frame("Sample_Type"=metadata$Sample_Type, "Sample_Color"=metadata$Sample_Color))
g_typ_meta<-merge(typ_meta,b.g.typ_m, by="Sample_Type")
g_typ_meta<-subset(g_typ_meta, Genus!="Unknown")
g_typ_meta<-subset(g_typ_meta, Count!=0)

ts1<-ggplot(g_typ_meta[g_typ_meta$Count>=0.0005,], aes(Genus,-Count)) +
  geom_jitter(aes(color=factor(Sample_Type)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Sample Type", values=c(unique(g_typ_meta$Sample_Color[order(g_typ_meta$Sample_Type)])),c("Dust","Seawater")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=0.5,angle=90),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Bacteria/Archaea & Sample Type")

ggsave(ts1,filename = "figures/RelativeAbundance/16S_Gen.RA.png", width=20, height=10, dpi=600)

#ts1a<-ggplot(g_typ_meta, aes(Genus, Count, label=Genus)) +
#  geom_jitter(aes(color=ifelse(Count>0.01,factor(Sample_Type),"grey")), size=2, width=0.15, height=0) +
#  scale_color_manual(name ="Sample Type", values=c(unique(g_typ_meta$Sample_Color[order(g_typ_meta$Sample_Type)]),"grey"),c("Dust","Seawater","RA < 1%")) +
#  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
#  labs(x="Microbial Genera", y="Relative Abundance", title="Microbial Genera & Sample Type")

#ggsave(ts1a,filename = "figures/RelativeAbundance/16S_Phyla.RA2.png", width=15, height=10, dpi=600)

g.typ.05p<-na.omit(subset(g_typ_meta, Count>=0.005))

g.t.05<-ggplot(g.typ.05p, aes(Genus, Count)) +
  geom_jitter(aes(color=factor(Sample_Type)), size=2, width=0.15, height=0) + geom_boxplot(fill=NA, outlier.color=NA) +scale_color_manual(name ="Sample Type", values=c(unique(g.typ.05p$Sample_Color[order(g.typ.05p$Sample_Type)])),c("Seawater","Soil", "Dust")) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Bacteria/Archaea & Sample Type")

ggsave(g.t.05,filename = "figures/RelativeAbundance/16S_Gen.0.5percRA.png", width=15, height=10, dpi=600)

g.t.05a<-ggplot(g.typ.05p, aes(Genus, Count)) +
  geom_jitter(aes(color=ifelse(Count>0.01,factor(Sample_Type),"grey")), size=2, width=0.15, height=0) + geom_boxplot(fill=NA, outlier.color=NA) +scale_color_manual(name ="Sample Type", values=c(unique(g.typ.05p$Sample_Color[order(g.typ.05p$Sample_Type)]),"grey"),c("Seawater","Soil", "Dust","<1% RA")) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Bacteria/Archaea & Sample Type")

ggsave(g.t.05a,filename = "figures/RelativeAbundance/16S_Gen.0.5percRA.v2.png", width=15, height=10, dpi=600)

g.typ.1p<-subset(g_typ_meta, Count>=0.01)
g.typ.1p<-subset(g.typ.1p, Genus!="Unknown")

g2<-ggplot(g.typ.1p, aes(Genus, Count)) +
  geom_jitter(aes(color=factor(Sample_Type)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Sample Type", values=unique(g.typ.1p$Sample_Color[order(g.typ.1p$Sample_Type)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=0.5),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Microbial Genera by Sample Type", subtitle="Includes Taxa with a Relative Abundance of at least 1%")+coord_flip()

ggsave(g2,filename = "figures/RelativeAbundance/16S_Genera.RA_1percent_v1.png", width=12, height=10, dpi=600)

b.gen_RA.st<-ggplot(g.typ.1p, aes(x=Sample_Type, y=Count, fill=Genus))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genera by Sample Type", x="Sample Type", y="Relative Abundance", fill="Genus", subtitle="Includes Taxa with a Relative Abundance of at least 1%")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=0.5),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+guides(fill=guide_legend(ncol=1))

ggsave(b.gen_RA.st,filename = "figures/RelativeAbundance/bacterial_genera_1percent_RA_by_SampleType.png", width=12, height=10, dpi=600)


head(g_typ_meta)
g_typ_meta.no0<-subset(g_typ_meta, Count!=0)

tg.h1<-ggplot(g_typ_meta.no0, aes(Sample_Type, Genus, fill= Count)) +geom_tile()+scale_fill_gradient2(low="blue3", mid="white",high="red",midpoint=.025)+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.y=element_text(margin = margin(0,0)),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Sample Type", y="Microbial Genera", title="Microbial Genera & Sample Type",fill="Relative Abundance")+theme_classic()+scale_x_discrete(expand = c(0,0))

ggsave(tp.h1,filename = "figures/RelativeAbundance/16S_Phyla.RA_heatmap.png", width=12, height=10, dpi=600)


head(g_typ_meta)
length(g_typ_meta$Sample_Type)

d.gen<-subset(g_typ_meta, Sample_Type=="Dust")
sw.gen<-subset(g_typ_meta, Sample_Type=="Seawater")

colnames(d.gen)[which(names(d.gen) == "Count")] <- "D.RA"
d.gen$D.RA<-as.numeric(d.gen$D.RA)
#d.gen<-subset(d.gen, D.RA>0.0001)

colnames(sw.gen)[which(names(sw.gen) == "Count")] <- "SW.RA"
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

gencol = melt(c(Bacillus="darkgreen",Blastococcus="orange",Geodermatophilus="green3", Halomonas="darkslategray4", Marinobacter="blue3", Other="black", Paracoccus="firebrick1", Roseovarius="deeppink3", Salinicoccus="cornflowerblue", Sphingomonas="purple", Spiroplasma="deepskyblue1", Truepera="darkgoldenrod2"))
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

ggsave(test1,filename = "figures/RelativeAbundance/16S_Gen.RA_sample.type_scatterplot.png", width=12, height=10, dpi=600)

test2<-ggplot(tgen.comp, aes(SW.RA, D.RA, color=Genus))  +geom_point(aes(color=ifelse(SW.RA>0.01 | D.RA>0.01,Genus,"black")),size=2) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Seawater Relative Abundance", y="Dust Relative Abundance", title="Microbial Genera & Sample Type")

ggsave(test2,filename = "figures/RelativeAbundance/16S_Gen.RA_type_linear1_colortest.png", width=12, height=10, dpi=600)

ggplot(tgen.comp, aes(SW.RA, D.RA, color=ifelse(SW.RA>0.01 | D.RA>0.01,Genus,"black")))  +geom_point() + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Seawater Relative Abundance", y="Dust Relative Abundance", title="Microbial Genera & Sample Type")

tgc1<-ggplot(tgen.comp, aes(SW.RA, D.RA, color=Genus)) +geom_point(aes(color=D.RA>0.01 | SW.RA>0.01))+ theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Seawater Relative Abundance", y="Dust Relative Abundance", title="Microbial Genera & Sample Type", subtitle="Only shows shared taxa w/ RA > 0.0001")

ggsave(tgc1,filename = "figures/RelativeAbundance/16S_Gen.RA_type_linear1.png", width=12, height=10, dpi=600)
## Date of note: 11/1/21 vvv
# 16S_Gen.RA_type_linear1 - cutoff is 0.0001
# 16S_Gen.RA_type_linear2 - cutoff is 0.00001
# 16S_Gen.RA_type_linear 3 - cutoff is 0.000001
# ** only 2 genera shared at 0.001 cutoff: Halomonas, Truepera

ggplot(g_shared_typ, aes(x = Genus, y = Count2, fill = Sample_Type)) +
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

#### Save Everything ####
save.image("data/SSeawater_RelAb_Data.Rdata")

# Version Information
sessionInfo()