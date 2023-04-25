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
})

#### Load Global Env to Import Count/ASV Tables ####
load("data/SSeawater_Data_Ready.Rdata") # save global env to Rdata file
#save.image("data/Env_Seqs_All/env.seq_analysis.Rdata") # save global env to Rdata file
bac.dat.all[1:4,1:4]
bac.dat.all$Depth_m
head(meta_scaled)
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

# Barplot by SampleID
b.phy_RA<-ggplot(b.phyla_RA_meta, aes(x=SampleID, y=Count, fill=Phylum))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Phylum Relative Abundance", x="SampleID", y="Relative Abundance", fill="Phylum")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+guides(fill=guide_legend(ncol=2))+scale_y_continuous(expand = c(0,0),limits = c(0,1))

ggsave(b.phy_RA,filename = "figures/RelativeAbundance/16S_Phyla.RA_barplot.png", width=12, height=10, dpi=600)

head(b.phyla_RA_meta)

# Heatmap by SampleID
p.h1<-ggplot(b.phyla_RA_meta, aes(SampleID, Phylum, fill= Count)) +geom_tile()+scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=0.3)+
  theme_classic()+theme(axis.title.x = element_text(size=13,vjust=-0.5),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Sample ID", y="Microbial Phyla", title="Microbial Phyla & Sample Type",fill="Relative Abundance")+scale_x_discrete(expand = c(0,0))

ggsave(p.h1,filename = "figures/RelativeAbundance/16S_Phyla.RA_heatmap.png", width=12, height=10, dpi=600)

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
  scale_colour_gradient2(low="darkred",high="blue",midpoint=5.5,guide = guide_colourbar(reverse = TRUE)) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Phyla", y="Relative Abundance", title="Microbial Phyla & Depth",color="Depth (m)")

ggsave(tp1,filename = "figures/RelativeAbundance/SSW_16S_Phyla.RA_depth_taxasum.png", width=15, height=10, dpi=600)

tp1a<-ggplot(b.phy.dep_m[b.phy.dep_m$Count>0.05,], aes(Phylum, Count)) +
  geom_jitter(aes(color=as.numeric(Depth_m)), size=2, width=0.15, height=0) +
  scale_colour_gradient2(low="darkred",high="blue",midpoint=5.5,guide = guide_colourbar(reverse = TRUE)) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Phyla", y="Relative Abundance", title="Microbial Phyla & Depth",color="Depth (m)")

ggsave(tp1a,filename = "figures/RelativeAbundance/SSW_16S_Phyla.RA_depth_taxasum_5percent.png", width=15, height=10, dpi=600)

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

b.phy.date_m$SampDate<-factor(b.phy.date_m$SampDate, levels=c("June.2021","August.2021","December.2021","April.2022"))

# Barplot by Sample Date

psd0<-ggplot(b.phy.date_m, aes(x=SampDate, y=Count, fill=Phylum))+geom_bar(stat="identity",colour="black")+theme_classic()+
  labs(title = "Relative Abundance of Microbial Phyla", x="SampleID", y="Relative Abundance", fill="Phylum")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=2))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1))+ scale_x_discrete(labels=c("June.2021"="June 2021","August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022"))

ggsave(psd0,filename = "figures/RelativeAbundance/16S_Phyla.RA_barplot_sampdate.png", width=12, height=10, dpi=600)

# Taxonomic Summary by Sample Date

#b.phy.date_m2<-merge(b.phy.date_m, metadata, by="SampDate")

colorset1 # remember which date goes with each color

psd1<-ggplot(b.phy.date_m, aes(Phylum, Count)) +
  geom_jitter(aes(color=factor(SampDate)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Sample Date", values=c("#36ab57","#ff6f00","#26547c","#32cbff"), labels=c("June 2021","August 2021","December 2021","April 2022")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=0.5,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Phyla", y="Relative Abundance", title="Microbial Phyla & Sample Date")

ggsave(psd1,filename = "figures/RelativeAbundance/SSW_16S_Phyla.RA_date_taxasum.png", width=15, height=10, dpi=600)

psd1a<-ggplot(b.phy.date_m[b.phy.date_m$Count>0.05,], aes(Phylum, Count)) +
  geom_jitter(aes(color=factor(SampDate)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Sample Date", values=c("#36ab57","#ff6f00","#26547c","#32cbff"), labels=c("June 2021","August 2021","December 2021","April 2022")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=0.5,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Phyla", y="Relative Abundance", title="Microbial Phyla & Sample Date")

ggsave(psd1a,filename = "figures/RelativeAbundance/SSW_16S_Phyla.RA_date_taxasum_5percent.png", width=15, height=10, dpi=600)

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

# Barplot by SampleID

b.cls_RA<-ggplot(b.class_RA_meta, aes(x=SampleID, y=Count, fill=Class))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Relative Abundance of Microbial Classes", x="SampleID", y="Relative Abundance", fill="Class")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+guides(fill=guide_legend(ncol=2))+scale_y_continuous(expand = c(0,0),limits = c(0,1))

ggsave(b.cls_RA,filename = "figures/RelativeAbundance/16S_Class.RA_barplot.png", width=12, height=10, dpi=600)

head(b.class_RA_meta)

# Heatmap by SampleID
c.h1<-ggplot(b.class_RA_meta, aes(SampleID, Class, fill= Count)) +geom_tile()+scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=0.3)+
  theme_classic()+theme(axis.title.x = element_text(size=13,vjust=-0.5),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Sample ID", y="Microbial Class", title="Microbial Class & Sample Type",fill="Relative Abundance")+scale_x_discrete(expand = c(0,0))

ggsave(c.h1,filename = "figures/RelativeAbundance/16S_class.RA_heatmap.png", width=12, height=10, dpi=600)

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
#+ scale_x_discrete(labels=c("June.2021"="June 2021","August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022"))

ggsave(cd1,filename = "figures/RelativeAbundance/SSW_16S_class.RA_barplot_depth.png", width=12, height=10, dpi=600)

cd1a<-ggplot(b.cls.dep_m[b.cls.dep_m$Count>0.05,], aes(x=Depth_m, y=Count, fill=Class))+geom_bar(stat="identity",colour="black")+theme_classic()+
  labs(title = "Relative Abundance of Microbial Classes", x="SampleID", y="Relative Abundance", fill="Class",subtitle="Only Relative Abundance > 5%")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1))
#+ scale_x_discrete(labels=c("June.2021"="June 2021","August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022"))

ggsave(cd1a,filename = "figures/RelativeAbundance/SSW_16S_class.RA_barplot_depth_5percent.png", width=12, height=10, dpi=600)

# Taxonomic Summary by Depth

cd2<-ggplot(b.cls.dep_m, aes(Class, Count)) +
  geom_jitter(aes(color=as.numeric(Depth_m)), size=2, width=0.15, height=0) +
  scale_colour_gradient2(high="blue",low="darkred",midpoint=5.5,guide = guide_colourbar(reverse = TRUE)) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Class", y="Relative Abundance", title="Microbial Class & Depth",color="Depth (m)")

ggsave(cd2,filename = "figures/RelativeAbundance/SSW_16S_class.RA_depth_taxasum.png", width=15, height=10, dpi=600)

cd2a<-ggplot(b.cls.dep_m[b.cls.dep_m$Count>0.05,], aes(Class, Count)) +
  geom_jitter(aes(color=as.numeric(Depth_m)), size=2, width=0.15, height=0) +
  scale_colour_gradient2(low="darkred",high="blue",midpoint=5.5,guide = guide_colourbar(reverse = TRUE)) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Class", y="Relative Abundance", title="Microbial Class & Depth",color="Depth (m)")

ggsave(cd2a,filename = "figures/RelativeAbundance/SSW_16S_class.RA_depth_5percent_taxasum.png", width=15, height=10, dpi=600)

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

b.cls.date_m$SampDate<-factor(b.cls.date_m$SampDate, levels=c("June.2021","August.2021","December.2021","April.2022"))

# Barplot by Sample Date

csd1<-ggplot(b.cls.date_m, aes(x=SampDate, y=Count, fill=Class))+geom_bar(stat="identity",colour="black")+theme_classic()+
  labs(title = "Relative Abundance of Microbial Classes", x="SampleID", y="Relative Abundance", fill="Class")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=3))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1))+ scale_x_discrete(labels=c("June.2021"="June 2021","August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022"))

ggsave(csd1,filename = "figures/RelativeAbundance/SSW_16S_class.RA_barplot_sampdate.png", width=12, height=10, dpi=600)

csd1<-ggplot(b.cls.date_m[b.cls.date_m$Count>0.05,], aes(x=SampDate, y=Count, fill=Class))+geom_bar(stat="identity",colour="black")+theme_classic()+
  labs(title = "Relative Abundance of Microbial Classes", x="SampleID", y="Relative Abundance", fill="Class",subtitle="Only Relative Abundance > 5%")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1))+ scale_x_discrete(labels=c("June.2021"="June 2021","August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022"))

ggsave(csd1,filename = "figures/RelativeAbundance/SSW_16S_class.RA_barplot_sampdate_5percent.png", width=12, height=10, dpi=600)

# Taxonomic Summary by Sample Date

#b.cls.date_m2<-merge(b.cls.date_m, metadata, by="SampDate")

colorset1 # remember which date goes with each color

csd2<-ggplot(b.cls.date_m, aes(Class, Count)) +
  geom_jitter(aes(color=factor(SampDate)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Sample Date", values=c("#36ab57","#ff6f00","#26547c","#32cbff"), labels=c("June 2021","August 2021","December 2021","April 2022")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=0.5,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Class", y="Relative Abundance", title="Microbial Class & Sample Date")

ggsave(csd2,filename = "figures/RelativeAbundance/SSW_16S_class_RA_date_taxasum.png", width=15, height=10, dpi=600)

csd3<-ggplot(unique(b.cls.date_m[b.cls.date_m$Count>0.1,]), aes(Class, Count)) +
  geom_jitter(aes(color=factor(SampDate)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Sample Date", values=c("#36ab57","#ff6f00","#26547c","#32cbff"), labels=c("June 2021","August 2021","December 2021","April 2022")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=0.5,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Class", y="Relative Abundance", title="Microbial Class & Sample Date")

ggsave(csd3,filename = "figures/RelativeAbundance/SSW_16S_class.RA_date_taxasum_5percent.png", width=15, height=10, dpi=600)

#### Family Relative Abundance ####

# use dcast to count up ASVs within each Family across all of the samples
b.fam_counts <- as.data.frame(dcast(bac.dat.all, SampleID~Family, value.var="Count", fun.aggregate=sum)) ###
head(b.fam_counts) # counts by fam per sample
dim(b.fam_counts)
rownames(b.fam_counts)<-b.fam_counts$SampleID
colnames(b.fam_counts)<-gsub(" ", ".",colnames(b.fam_counts))
b.fam_counts<-subset(b.fam_counts, select=-c(Unknown, Unknown.Family))

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

# Barplot by SampleID

b.fam_RA<-ggplot(b.fam_RA_meta, aes(x=SampleID, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Relative Abundance of Microbial Families", x="SampleID", y="Relative Abundance", fill="Family")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+guides(fill=guide_legend(ncol=2))+scale_y_continuous(expand = c(0,0),limits = c(0,1))

ggsave(c.h1,filename = "figures/RelativeAbundance/16S_fam.RA_barplot.png", width=12, height=10, dpi=600)

head(b.fam_RA_meta)

# Heatmap by SampleID

p.h1<-ggplot(b.fam_RA_meta, aes(SampleID, Family, fill= Count)) +geom_tile()+scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=0.3)+
  theme_classic()+theme(axis.title.x = element_text(size=13,vjust=-0.5),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Sample ID", y="Microbial Family", title="Microbial Families & Sample Type",fill="Relative Abundance")+scale_x_discrete(expand = c(0,0))

ggsave(c.h1,filename = "figures/RelativeAbundance/16S_fam.RA_heatmap.png", width=12, height=10, dpi=600)

bac.dat.all[1:4,1:4]

# by Family + depth
bac.fam.dep <- as.data.frame(dcast(bac.dat.all,Depth_m~Family, value.var="Count", fun.aggregate=sum)) ###
head(bac.fam.dep) # counts by Family + sample depe
rownames(bac.fam.dep)<-bac.fam.dep$Depth_m
colnames(bac.fam.dep)<-gsub(" ", ".",colnames(bac.fam.dep))
bac.fam.dep<-subset(bac.fam.dep, select=-c(Unknown, Unknown.Family))

b.RA_fam.dep<-data.frame(decostand(bac.fam.dep[,-1], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(b.RA_fam.dep) # sanity check
b.RA_fam.dep$Depth_m<-rownames(b.RA_fam.dep) # Depth_m is now a character, not a factor!
head(b.RA_fam.dep)

#melt down relativized data to merge with metadata
b.fam.dep_m<-melt(b.RA_fam.dep, by="Depth_m")

head(b.fam.dep_m)
colnames(b.fam.dep_m)[which(names(b.fam.dep_m) == "variable")] <- "Family"
colnames(b.fam.dep_m)[which(names(b.fam.dep_m) == "value")] <- "Count"
head(b.fam.dep_m) ## relative abundance based on sum of counts by Family!

#dep_meta<-unique(data.frame("Depth_m"=metadata$Depth_m, "Sample_Color"=metadata$Sample_Color))
#p_dep_meta<-merge(dep_meta,b.fam.dep_m, by="Depth_m")

# Barplot by Depth

fd1<-ggplot(b.fam.dep_m, aes(x=Depth_m, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+theme_classic()+
  labs(title = "Relative Abundance of Microbial Classes", x="Depth (m)", y="Relative Abundance", fill="Class")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=3))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1))
#+ scale_x_discrete(labels=c("June.2021"="June 2021","August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022"))

ggsave(fd1,filename = "figures/RelativeAbundance/SSW_16S_fam.RA_barplot_depth.png", width=12, height=10, dpi=600)

fd1a<-ggplot(b.fam.dep_m[b.fam.dep_m$Count>0.05,], aes(x=Depth_m, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+theme_classic()+
  labs(title = "Relative Abundance of Microbial Classes", x="SampleID", y="Relative Abundance", fill="Class",subtitle="Only Relative Abundance > 5%")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1))
#+ scale_x_discrete(labels=c("June.2021"="June 2021","August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022"))

ggsave(fd1a,filename = "figures/RelativeAbundance/SSW_16S_fam.RA_barplot_depth_5percent.png", width=12, height=10, dpi=600)

# Taxonomic Summary by Depth

fd2<-ggplot(b.fam.dep_m, aes(Family, Count)) +
  geom_jitter(aes(color=as.numeric(Depth_m)), size=2, width=0.15, height=0) +
  scale_colour_gradient2(high="blue",low="darkred",midpoint=5.5,guide = guide_colourbar(reverse = TRUE)) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Family", y="Relative Abundance", title="Microbial Families & Depth",color="Depth (m)")

ggsave(fd2,filename = "figures/RelativeAbundance/SSW_16S_fam.RA_depth_taxasum.png", width=15, height=10, dpi=600)

fd2a<-ggplot(b.fam.dep_m[b.fam.dep_m$Count>0.05,], aes(Family, Count)) +
  geom_jitter(aes(color=as.numeric(Depth_m)), size=2, width=0.15, height=0) +
  scale_colour_gradient2(low="darkred",high="blue",midpoint=5.5,guide = guide_colourbar(reverse = TRUE)) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Family", y="Relative Abundance", title="Microbial Families & Depth",color="Depth (m)")

ggsave(fd2a,filename = "figures/RelativeAbundance/SSW_16S_fam.RA_depth_5percent_taxasum.png", width=15, height=10, dpi=600)

# by Family + Sampling Date
bac.fam.date <- as.data.frame(dcast(bac.dat.all,SampDate~Family, value.var="Count", fun.aggregate=sum)) ###
head(bac.fam.date) # counts by Family + sample depe
rownames(bac.fam.date)<-bac.fam.date$SampDate
colnames(bac.fam.date)<-gsub(" ", ".",colnames(bac.fam.date))
bac.fam.date<-subset(bac.fam.date, select=-c(Unknown, Unknown.Family))

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

b.fam.date_m$SampDate<-factor(b.fam.date_m$SampDate, levels=c("June.2021","August.2021","December.2021","April.2022"))

# Barplot by Sample Date

fsd1<-ggplot(b.fam.date_m, aes(x=SampDate, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+theme_classic()+
  labs(title = "Relative Abundance of Microbial Families", x="SampleID", y="Relative Abundance", fill="Family")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=3))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1))+ scale_x_discrete(labels=c("June.2021"="June 2021","August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022"))

ggsave(fsd1,filename = "figures/RelativeAbundance/SSW_16S_fam.RA_barplot_sampdate.png", width=12, height=10, dpi=600)

fsd1a<-ggplot(b.fam.date_m[b.fam.date_m$Count>0.01,], aes(x=SampDate, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+theme_classic()+
  labs(title = "Relative Abundance of Microbial Families", x="SampleID", y="Relative Abundance", fill="Family",subtitle="Only Relative Abundance > 1%")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=2))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1))+ scale_x_discrete(labels=c("June.2021"="June 2021","August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022"))

ggsave(fsd1a,filename = "figures/RelativeAbundance/SSW_16S_fam.RA_barplot_sampdate_1percent.png", width=12, height=10, dpi=600)

# Taxonomic Summary by Sample Date

#b.fam.date_m2<-merge(b.fam.date_m, metadata, by="SampDate")

colorset1 # remember which date goes with each color

fsd2<-ggplot(b.fam.date_m, aes(Family, Count)) +
  geom_jitter(aes(color=factor(SampDate)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Sample Date", values=c("#36ab57","#ff6f00","#26547c","#32cbff"), labels=c("June 2021","August 2021","December 2021","April 2022")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=0.5,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Family", y="Relative Abundance", title="Microbial Families & Sample Date")

ggsave(fsd2,filename = "figures/RelativeAbundance/SSW_16S_fam.RA_date_taxasum.png", width=15, height=10, dpi=600)

fsd3<-ggplot(b.fam.date_m[b.fam.date_m$Count>0.1,], aes(Family, Count)) +
  geom_jitter(aes(color=factor(SampDate)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Sample Date", values=c("#36ab57","#ff6f00","#26547c","#32cbff"), labels=c("June 2021","August 2021","December 2021","April 2022")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=0.5,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Family", y="Relative Abundance", title="Microbial Families & Sample Date")

ggsave(fsd3,filename = "figures/RelativeAbundance/SSW_16S_fam.RA_date_taxasum_5percent.png", width=15, height=10, dpi=600)

#### Genus Relative Abundance ####
# use dcast to count up ASVs within each Genus across all of the samples
bac.dat.all.g<-subset(bac.dat.all, bac.dat.all$Genus!="Unknown")
"Unknown" %in% bac.dat.all.g$Genus

b.genus_counts <- as.data.frame(dcast(bac.dat.all.g, SampleID~Genus+Species, value.var="Count", fun.aggregate=sum)) ###
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
colnames(b.genus_m)[which(names(b.genus_m) == "variable")] <- "Genus_species"
colnames(b.genus_m)[which(names(b.genus_m) == "value")] <- "Count"
head(b.genus_m) ## relative abundance based on sum of counts by genus!
b.genus_m$Genus_species<-gsub("^X.","",b.genus_m$Genus_species) # get rid of leading X. in Genus_species names
b.genus_m$Genus_species<-gsub("\\.\\."," ",b.genus_m$Genus_species) # get rid of .. in species name --> . is regex
head(b.genus_m) ## relative abundance based on sum of counts by genus!

b.genus_RA_meta<-merge(b.genus_m,metadata, by="SampleID")
head(b.genus_RA_meta) ## relative abundance based on sum of counts by genus!
max(b.genus_RA_meta$Count)

# Barplot by SampleID

b.gen_RA1<-ggplot(b.genus_RA_meta[b.genus_RA_meta$Count>0.05,], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genus Relative Abundance", x="SampleID", y="Relative Abundance", subtitle="Includes taxa with Relative Abundance > 5%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))

ggsave(b.gen_RA1,filename = "figures/RelativeAbundance/SSW_16S_Genera.Spec.RA_barplot_5perc.png", width=12, height=10, dpi=600)

b.gen_RA2<-ggplot(b.genus_RA_meta[b.genus_RA_meta$Count>0.10,], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genus Relative Abundance", x="SampleID", y="Relative Abundance", subtitle="Includes taxa with Relative Abundance > 10%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))

ggsave(b.gen_RA2,filename = "figures/RelativeAbundance/SSW_16S_Genera.Spec.RA_barplot_10perc.png", width=12, height=10, dpi=600)

b.gen_RA3<-ggplot(b.genus_RA_meta[b.genus_RA_meta$Count>0.25,], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genus Relative Abundance", x="SampleID", y="Relative Abundance", subtitle="Includes taxa with Relative Abundance > 25%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))

ggsave(b.gen_RA3,filename = "figures/RelativeAbundance/SSW_16S_Genera.Spec.RA_barplot_25perc.png", width=12, height=10, dpi=600)

b.gen_RA4<-ggplot(b.genus_RA_meta[b.genus_RA_meta$Count>0.35,], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Genus Relative Abundance", x="SampleID", y="Relative Abundance", subtitle="Includes taxa with Relative Abundance > 35%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))

ggsave(b.gen_RA4,filename = "figures/RelativeAbundance/SSW_16S_Genera.Spec.RA_barplot_35perc.png", width=12, height=10, dpi=600)

# prep for heatmap
max(b.genus_RA_meta$Count)
mean(b.genus_RA_meta$Count)
median(b.genus_RA_meta$Count)
max(b.genus_RA_meta$Count)/2 # what is the mid point of the RA here?

# Heatmap by SampleID

g.h1<-ggplot(b.genus_RA_meta[b.genus_RA_meta$Count>0.01,], aes(SampleID, Genus_species, fill= Count)) +geom_tile()+scale_fill_gradient2(low="lightblue",mid="white",high="orange",midpoint=0.3)+
  theme_classic()+theme(axis.title.x = element_text(size=13,vjust=-0.5),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Sample ID", y="Microbial Genera", title="Microbial Genera & Sample Type",subtitle="Includes taxa with Relative Abundance > 1%",fill="Relative Abundance")+scale_x_discrete(expand = c(0,0))

ggsave(g.h1,filename = "figures/RelativeAbundance/16S_Genera.RA_heatmap_A_1perc.png", width=20, height=15, dpi=600)

g.h2<-ggplot(b.genus_RA_meta[b.genus_RA_meta$Count>0.05,], aes(SampleID, Genus_species, fill= Count)) +geom_tile()+scale_fill_gradient2(low="lightblue",mid="white",high="orange",midpoint=0.3)+
  theme_classic()+theme(axis.title.x = element_text(size=13,vjust=-0.5),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Sample ID", y="Microbial Genera", title="Microbial Genera & Sample Type",subtitle="Includes taxa with Relative Abundance > 5%",fill="Relative Abundance")+scale_x_discrete(expand = c(0,0))

ggsave(g.h2,filename = "figures/RelativeAbundance/16S_Genera.RA_heatmap_B_5perc.png", width=16, height=10, dpi=600)

bac.dat.all[1:4,1:4]

# by Genus + depth
bac.gen.dep <- as.data.frame(dcast(bac.dat.all.g,Depth_m~Genus, value.var="Count", fun.aggregate=sum)) ###
head(bac.gen.dep) # counts by Genus + sample depth
rownames(bac.gen.dep)<-bac.gen.dep$Depth_m

b.RA_gen.dep<-data.frame(decostand(bac.gen.dep[,-1], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(b.RA_gen.dep) # sanity check
b.RA_gen.dep$Depth_m<-rownames(b.RA_gen.dep) # Depth_m is now a character, not a factor!
head(b.RA_gen.dep)

#melt down relativized data to merge with metadata
b.gen.dep_m<-melt(b.RA_gen.dep, by="Depth_m")

head(b.gen.dep_m)
colnames(b.gen.dep_m)[which(names(b.gen.dep_m) == "variable")] <- "Genus"
colnames(b.gen.dep_m)[which(names(b.gen.dep_m) == "value")] <- "Count"
head(b.gen.dep_m) ## relative abundance based on sum of counts by Genus!
b.gen.dep_m$Genus<-gsub("^X.","",b.gen.dep_m$Genus) # get rid of leading X. in Genus names
b.gen.dep_m$Genus<-gsub("\\.\\."," ",b.gen.dep_m$Genus) # get rid of .. in species name --> . is regex
head(b.gen.dep_m) ## relative abundance based on sum of counts by genus!

# Barplot by Depth

gd1<-ggplot(b.gen.dep_m, aes(x=Depth_m, y=Count, fill=Genus))+geom_bar(stat="identity",colour="black")+theme_classic()+
  labs(title = "Relative Abundance of Microbial Classes", x="Depth (m)", y="Relative Abundance", fill="Class")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=3))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1))
#+ scale_x_discrete(labels=c("June.2021"="June 2021","August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022"))

ggsave(gd1,filename = "figures/RelativeAbundance/SSW_16S_Genus.RA_barplot_depth.png", width=12, height=10, dpi=600)

gd1a<-ggplot(b.gen.dep_m[b.gen.dep_m$Count>0.05,], aes(x=Depth_m, y=Count, fill=Genus))+geom_bar(stat="identity",colour="black")+theme_classic()+
  labs(title = "Relative Abundance of Microbial Classes", x="SampleID", y="Relative Abundance", fill="Class",subtitle="Only Relative Abundance > 5%")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1))
#+ scale_x_discrete(labels=c("June.2021"="June 2021","August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022"))

ggsave(gd1a,filename = "figures/RelativeAbundance/SSW_16S_Genus.RA_barplot_depth_5percent.png", width=12, height=10, dpi=600)

# Taxonomic Summary by Depth

#dep_meta<-unique(data.frame("Depth_m"=metadata$Depth_m, "Sample_Color"=metadata$Sample_Color))
#p_dep_meta<-merge(dep_meta,b.gen.dep_m, by="Depth_m")
tg1<-ggplot(b.gen.dep_m[b.gen.dep_m$Count>0.01,], aes(Genus, Count)) +
  geom_jitter(aes(color=as.numeric(Depth_m)), size=2, width=0.15, height=0) +
  scale_colour_gradient2(low="darkred",high="blue",midpoint=5.5,guide = guide_colourbar(reverse = TRUE)) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Microbial Genera & Depth", subtitle="Includes taxa with Relative Abundance > 1%",color="Depth (m)")

ggsave(tg1,filename = "figures/RelativeAbundance/SSW_16S_Genera.RA_depth_taxasum_1perc.png", width=15, height=10, dpi=600)

tg1a<-ggplot(b.gen.dep_m[b.gen.dep_m$Count>0.05,], aes(Genus, Count)) +
  geom_jitter(aes(color=as.numeric(Depth_m)), size=2, width=0.15, height=0) +
  scale_colour_gradient2(low="darkred",high="blue",midpoint=5.5,guide = guide_colourbar(reverse = TRUE)) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Microbial Genera & Depth", subtitle="Includes taxa with Relative Abundance > 5%",color="Depth (m)")

ggsave(tg1a,filename = "figures/RelativeAbundance/SSW_16S_Genera.RA_depth_taxasum_5percent.png", width=15, height=10, dpi=600)

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
b.gen.date_m$Genus<-gsub("^X.","",b.gen.date_m$Genus) # get rid of leading X. in Genus_species names
b.gen.date_m$Genus<-gsub("\\.\\."," ",b.gen.date_m$Genus) # get rid of .. in species name --> . is regex
head(b.genus_m) ## relative abundance based on sum of counts by genus!

b.gen.date_m$SampDate<-factor(b.gen.date_m$SampDate, levels=c("June.2021","August.2021","December.2021","April.2022"))

# Barplot by Sample Date

gsd1<-ggplot(b.gen.date_m[b.gen.date_m$Count>0.01,], aes(x=SampDate, y=Count, fill=Genus))+geom_bar(stat="identity",colour="black")+theme_classic()+
  labs(title = "Relative Abundance of Microbial Genera", x="SampleID", subtitle="Includes taxa with Relative Abundance > 1%",y="Relative Abundance", fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=2))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1))+ scale_x_discrete(labels=c("June.2021"="June 2021","August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022"))

ggsave(gsd1,filename = "figures/RelativeAbundance/SSW_16S_Genera.RA_date_barplot.png", width=15, height=10, dpi=600)

# Taxonomic Summary by Sample Date

#b.gen.date_m2<-merge(b.gen.date_m, metadata, by="SampDate")

colorset1 # remember which date goes with each color

gsd1<-ggplot(b.gen.date_m[b.gen.date_m$Count>0.01,], aes(Genus, Count)) +
  geom_jitter(aes(color=factor(SampDate)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Sample Date", values=c("#36ab57","#ff6f00","#26547c","#32cbff"), labels=c("June 2021","August 2021","December 2021","April 2022")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=0.5,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Microbial Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 1%")

ggsave(gsd1,filename = "figures/RelativeAbundance/SSW_16S_Genera.RA_date_taxasum_1perc.png", width=15, height=10, dpi=600)

gsd1a<-ggplot(b.gen.date_m[b.gen.date_m$Count>0.05,], aes(Genus, Count)) +
  geom_jitter(aes(color=factor(SampDate)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Sample Date", values=c("#36ab57","#ff6f00","#26547c","#32cbff"), labels=c("June 2021","August 2021","December 2021","April 2022")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=0.5,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Microbial Genera & Sample Date")

ggsave(gsd1a,filename = "figures/RelativeAbundance/SSW_16S_Genera.RA_date_taxasum_5percent.png", width=15, height=10, dpi=600)

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
b.gen.date.dep_m$Genus<-gsub("^X.","",b.gen.date.dep_m$Genus) # get rid of leading X. in Genus_species names
b.gen.date.dep_m$Genus<-gsub("\\.\\."," ",b.gen.date.dep_m$Genus) # get rid of .. in species name --> . is regex
head(b.gen.date.dep_m) ## relative abundance based on sum of counts by genus!

b.gen.date.dep_m2<-as.data.frame(separate_wider_delim(data = b.gen.date.dep_m, col=SampDate_Depth, "_", names = c("SampDate", "Depth_m"))) # Separate SampDate & Depth column for Heatmap later
b.gen.date.dep_m2$SampDate_Depth<-interaction(b.gen.date.dep_m2$SampDate,b.gen.date.dep_m2$Depth_m)

b.gen.date.dep_m2$Depth_m<-factor(b.gen.date.dep_m2$Depth_m, levels=c("0","2","3","4","5","7","9","10","11"))
b.gen.date.dep_m2$SampDate<-factor(b.gen.date.dep_m2$SampDate,levels=c("June.2021","August.2021","December.2021","April.2022"))
b.gen.date.dep_m2$SampDate_Depth = factor(b.gen.date.dep_m2$SampDate_Depth, levels=unique(b.gen.date.dep_m2$SampDate_Depth[order(b.gen.date.dep_m2$Depth_m,b.gen.date.dep_m2$SampDate)]), ordered=TRUE)

g.sd.d.h1<-ggplot(b.gen.date.dep_m2[b.gen.date.dep_m2$Count>0.01,], aes(SampDate_Depth, Genus, fill= Count)) +geom_tile()+scale_fill_gradient2(low="lightblue",mid="white",high="orange",midpoint=0.3)+
  theme_classic()+theme(axis.title.x = element_text(size=13,vjust=-0.5),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Sample ID", y="Microbial Genera", title="Microbial Genera & Sample Type",subtitle="Includes taxa with Relative Abundance > 1%",fill="Relative Abundance")+scale_x_discrete(expand = c(0,0))

ggsave(g.sd.d.h1,filename = "figures/RelativeAbundance/16S_Genera.RA_heatmap_date_depth_1perc.png", width=20, height=15, dpi=600)

g.sd.d.h2<-ggplot(b.gen.date.dep_m2[b.gen.date.dep_m2$Count>0.05,], aes(SampDate_Depth, Genus, fill= Count)) +geom_tile()+scale_fill_gradient2(low="lightblue",mid="white",high="orange",midpoint=0.3)+
  theme_classic()+theme(axis.title.x = element_text(size=13,vjust=-0.5),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Sample ID", y="Microbial Genera", title="Microbial Genera & Sample Type",subtitle="Includes taxa with Relative Abundance > 1%",fill="Relative Abundance")+scale_x_discrete(expand = c(0,0))

ggsave(g.sd.d.h2,filename = "figures/RelativeAbundance/16S_Genera.RA_heatmap_date_depth_5perc.png", width=20, height=15, dpi=600)

g.sd.d.hm.1<-ggplot(b.gen.date.dep_m2[b.gen.date.dep_m2$Count>0.05,], aes(Genus, Count)) +
  geom_jitter(aes(color=as.numeric(as.character(Depth_m)),shape=SampDate), size=3, width=0.15, height=0) +
  scale_colour_gradient2(low="darkred",high="blue",midpoint=5.5,guide = guide_colourbar(reverse = TRUE)) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() + scale_shape_discrete(labels=c("June 2021","August 2021","December 2021","April 2022"),name="Sample Date") +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Microbial Genera at Depth by Sample Date",color="Depth (m)", shape="Sample Date")

ggsave(g.sd.d.hm.1,filename = "figures/RelativeAbundance/SSW_16S_Genera.RA_date_depth_taxasum.png", width=15, height=10, dpi=600)
## ^ this figure includes the relative abundance of each organism by depth & date!!!

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

tg.h1<-ggplot(g_typ_meta.no0, aes(Sample_Type, Genus, fill= Count)) +geom_tile()+scale_fill_gradient2(low="blue", mid="white",high="red",midpoint=.025)+
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



# Version Information
sessionInfo()