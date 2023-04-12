library(phyloseq)
library(ggplot2)
library (vegan)
library(ggpubr)
library(scales)
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
library("dendextend")
library("tidyr")
library("viridis")
library("reshape")
library(reshape2)
library(wesanderson)
library(nationalparkcolors)
library(shades)

getwd()


#### Import metadata ####
metadata<-as.data.frame(read.csv("data/Mouse_Seqs/mapzym.csv"), header=TRUE)
head(metadata)
metadata<-na.omit(metadata)

metadata$Description<-gsub("\xca", "", metadata$Description) ## gsub is global sub (does not just remove first instance of pattern, but multiple)
rownames(metadata)<-metadata$SampleID
head(metadata)

metadata$lungTissue<-sub("yes", "lung", metadata$lungTissue)
metadata$lungTissue<-sub("no", "BALF", metadata$lungTissue)
head(metadata)

# create colors for later figures
colorset = melt(c(Alt="#540b0e",Con="#4c956c",Silica="#98c1d9"))
colorset$Exposed<-rownames(colorset)
colnames(colorset)[which(names(colorset) == "value")] <- "color"
colorset

metadata<-merge(metadata, colorset, by="Exposed")
head(metadata)
metadata$color <- as.character(metadata$color)
rownames(metadata)<-metadata$SampleID

#### Import and reformat the data ####
# bacteria first
bac.ASV_counts<-data.frame(readRDS("data/Mouse_Seqs/Results_16S/16S_ASVs_Counts_dada2_6.20.2021_Robject.rds", refhook = NULL))
dim(bac.ASV_counts)
head(bac.ASV_counts)
rownames(bac.ASV_counts)<-bac.ASV_counts$ASV_ID

# fungi next
its2.ASV_counts<-data.frame(readRDS("data/Mouse_Seqs/Results_ITS2/ITS2_ASVs_Counts_dada2_6.20.2021_Robject.rds", refhook = NULL))
dim(its2.ASV_counts)
head(its2.ASV_counts)
rownames(its2.ASV_counts)<-its2.ASV_counts$ASV_ID

### ^^^ *** Note to self - before you run stats, make sure your table is in a SAMPLE x ASV/SPECIES orientation!!

#### Import ASV to taxa sheet for ASV identification ####

# bacteria first
bac.ASV_tax<-data.frame(readRDS("data/Mouse_Seqs/Results_16S/ASVs_Taxonomy_dada2_6.16.2021_Robject.rds", refhook = NULL))
head(bac.ASV_tax)
bac.ASV_tax$Genus<-gsub("\\[(.*)\\]", "\\1", bac.ASV_tax$Genus) ## drop brackets around Eubacterium while not losing the string inside brackets

bac.ASV_tax[is.na(bac.ASV_tax)]<- "Unknown"
bac.ASV_tax$Species<-gsub("Unknown", "unknown", bac.ASV_tax$Species) ## drop brackets around Eubacterium while not losing the string inside brackets

head(bac.ASV_tax)
class(bac.ASV_tax)
bac.ASV_tax$ASV_ID<-rownames(bac.ASV_tax)
head(bac.ASV_tax)

# fungi next
its2.ASV_tax<-data.frame(readRDS("data/Mouse_Seqs/Results_ITS2/ITS2_ASVs_Taxonomy_dada2_6.16.2021_Robject.rds", refhook = NULL))
head(its2.ASV_tax)
its2.ASV_tax$ASV_ID<-rownames(its2.ASV_tax)
its2.ASV_tax<-as.data.frame(lapply(its2.ASV_tax, function(x) gsub("[a-z]__", "", x)))

its2.ASV_tax[is.na(its2.ASV_tax)]<- "Unknown"
its2.ASV_tax$Species<-gsub("Unknown", "unknown", its2.ASV_tax$Species) ## drop brackets around Eubacterium while not losing the string inside brackets

head(its2.ASV_tax)
rownames(its2.ASV_tax)<-its2.ASV_tax$ASV_ID
head(its2.ASV_tax)

#### Data Formatting and Transformation ####

# bacteria first
bac.ASV_dat<-merge(bac.ASV_counts,bac.ASV_tax, by="ASV_ID")
head(bac.ASV_dat)

bac.ASV_dat<-subset(bac.ASV_dat, Kingdom!="Unknown") ## keep only bacteria and archaean -- drop Unknowns
bac.ASV_dat<-subset(bac.ASV_dat, Phylum!="Unknown") ## keep only bacteria and archaean -- drop Unknowns=
head(bac.ASV_dat)
bac.ASV_dat<-subset(bac.ASV_dat, Class!="Chloroplast") ## keep only bacteria -- exclude Chloroplast sequences
bac.ASV_dat<-subset(bac.ASV_dat, Order!="Chloroplast") ## keep only bacteria -- exclude Chloroplast sequences
bac.ASV_dat<-subset(bac.ASV_dat, Family!="Mitochondria") ## keep only bacteria -- exclude Chloroplast sequences

'Chloroplast' %in% bac.ASV_dat # check if Chloroplast counts are still in df, should be false because they've been removed
'Mitochondria' %in% bac.ASV_dat # check if Chloroplast counts are still in df, should be false because they've been removed

head(bac.ASV_dat)
rownames(bac.ASV_dat)<-bac.ASV_dat$ASV_ID
head(bac.ASV_dat)

"Undetermined" %in% bac.ASV_dat

## Create Sample x Species (ASV) table from counts
b.dat.m<-melt(bac.ASV_dat)
head(b.dat.m)
colnames(b.dat.m)[which(names(b.dat.m) == "variable")] <- "SampleID"
colnames(b.dat.m)[which(names(b.dat.m) == "value")] <- "Count"

bac.ASV_table<-as.data.frame(dcast(b.dat.m, SampleID~ASV_ID, value.var="Count", fun.aggregate=sum)) ###
head(bac.ASV_table)
rownames(bac.ASV_table)<-bac.ASV_table$SampleID
bac.ASV_table<-subset(bac.ASV_table, select=-c(SampleID))
head(bac.ASV_table)

# fungi next
its2.ASV_dat<-merge(its2.ASV_counts,its2.ASV_tax, by="ASV_ID")
head(its2.ASV_dat)

its2.ASV_dat<-subset(its2.ASV_dat, Kingdom!="Unknown") ## keep only its2teria and archaean -- drop Unknowns
its2.ASV_dat<-subset(its2.ASV_dat, Phylum!="Unknown") ## keep only its2teria and archaean -- drop Unknowns=
head(its2.ASV_dat)
its2.ASV_dat<-subset(its2.ASV_dat, Class!="Chloroplast") ## keep only its2teria -- exclude Chloroplast sequences
its2.ASV_dat<-subset(its2.ASV_dat, Order!="Chloroplast") ## keep only its2teria -- exclude Chloroplast sequences
its2.ASV_dat<-subset(its2.ASV_dat, Family!="Mitochondria") ## keep only its2teria -- exclude Chloroplast sequences

'Chloroplast' %in% its2.ASV_dat # check if Chloroplast counts are still in df, should be false because they've been removed
'Mitochondria' %in% its2.ASV_dat # check if Chloroplast counts are still in df, should be false because they've been removed

head(its2.ASV_dat)
rownames(its2.ASV_dat)<-its2.ASV_dat$ASV_ID
head(its2.ASV_dat)

"Undetermined" %in% its2.ASV_dat

## Create Sample x Species (ASV) table from counts
f.dat.m<-melt(its2.ASV_dat)
head(f.dat.m)
colnames(f.dat.m)[which(names(f.dat.m) == "variable")] <- "SampleID"
colnames(f.dat.m)[which(names(f.dat.m) == "value")] <- "Count"

its2.ASV_table<-as.data.frame(dcast(f.dat.m, SampleID~ASV_ID, value.var="Count", fun.aggregate=sum)) ###
head(its2.ASV_table)
rownames(its2.ASV_table)<-its2.ASV_table$SampleID
its2.ASV_table<-subset(its2.ASV_table, select=-c(SampleID))
head(its2.ASV_table)

## Reorder metadata to have same rows as ASV tables
metadata=metadata[rownames(its2.ASV_table),]
# ** ^ this indexing method will only work if the two dfs have the same # of rows AND the same row names!
bac.ASV_table=bac.ASV_table[rownames(its2.ASV_table),]

#### Checking out some color pallettes before we visualize ####
wes1<-wes_palette("Chevalier1")
wes2<-wes_palette("Moonrise3")
wes3<-wes_palette("IsleofDogs1")
wes4<-wes_palette("GrandBudapest1")
wes5<-wes_palette("GrandBudapest2")
#scale_fill_manual(values = wes_palette("IsleofDogs1"))

SM_pal <- park_palette("SmokyMountains") # create a palette and specify # of colors youw ant
Arc_pal <- park_palette("Arches") # create a palette and specify # of colors youw ant
CL_pal <- park_palette("CraterLake") # create a palette and specify # of colors youw ant
Sag_pal <- park_palette("Saguaro") # create a palette and specify # of colors youw ant
Aca_pal <- park_palette("Acadia") # create a palette and specify # of colors youw ant
DV_pal <- park_palette("DeathValley") # create a palette and specify # of colors youw ant
CI_pal <- park_palette("ChannelIslands") # create a palette and specify # of colors youw ant
Bad_pal <- park_palette("Badlands") # create a palette and specify # of colors youw ant
MR_pal <- park_palette("MtRainier") # create a palette and specify # of colors youw ant
HI_pal <- park_palette("Hawaii") # create a palette and specify # of colors youw ant

fair_cols <- c("#38170B","#BF1B0B", "#FFC465", "#66ADE5", "#252A52")
names(fair_cols) <- letters[1:5]
fair_ramp <- scales::colour_ramp(fair_cols)
fair_sat <- saturation(fair_ramp, 1)

# custom color palette
head(metadata)

warm2cold1<-get_palette(paste0("#", c("720026", "D14D60", "0077B6","03045E")),k=4)
names(warm2cold1) <- levels(metadata$Elevation)

warm2cold2<-get_palette(paste0("#", c("720026", "D14D60", "0077B6","03045E")),k=4)
names(warm2cold2) <- levels(metadata$Site)

warm2cold2<-get_palette(paste0("#", c("720026", "D14D60", "0077B6","03045E")),k=4)
names(warm2cold2) <- levels(metadata$Site)

#### Bray-Curtis Distance Matrices w/ relativized data ####
# *** using relativized data for Bray-Curtis distance matrices
## ITS1 first

class(its2.ASV_table) # raw counts
dim(its2.ASV_table)
row.names(its2.ASV_table)
head(its2.ASV_table)

## Relativize ASV Table
its2_RA_table<-data.frame(decostand(its2.ASV_table[,-1], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(its2_RA_table)
its2_RA_table<-its2_RA_table[which(rowSums(its2_RA_table)>0),]
rowSums(its2_RA_table)

## Create Bray-Curtis dissimilarity distance matrix
its2_RA_bray<-vegdist(its2_RA_table,method="bray")     ####### distance matrix with RELATIVIZED data!!!
# **** in vegan: ROWS need to be SITES/samples; COLUMNS are SPECIES (OTUs, ASVs) -- vegdist!
class(its2_RA_bray) # needs to be class dist for pcoa()

## 16S next
class(bac.ASV_table) # raw counts
dim(bac.ASV_table)
row.names(bac.ASV_table)

## Relativize ASV table
bac_RA_table<-data.frame(decostand(bac.ASV_table, method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(bac_RA_table)

## Create Bray-Curtis dissimilarity distance matrix
bac_RA_bray<-vegdist(bac_RA_table,method="bray")     ####### distance matrix with RELATIVIZED data!!!
# **** in vegan: ROWS need to be SITES/samples; COLUMNS are SPECIES (OTUs, ASVs) -- vegdist!
class(bac_RA_bray) # needs to be class dist for pcoa()

#### PCoA w/ relativized, BrayCurtis matrix ####

## fungi first
its2_pcoa1 = pcoa(its2_RA_bray) # pcoa of Bray-Curtis dissimilarity matrix made w/ relativized abundance (site x species) data
# The proportion of variances explained is in its element values$Relative_eig
#cmdscale(bac.bray, k=(nrow(bac.ASV_table)-1), eig=TRUE)

str(its2_pcoa1)
biplot(its2_pcoa1)
head(metadata)
#biplot(its2_pcoa1, meta_qdat[,2:9], main="ITS1 PCoA + Scaled Chemical Data", xlab = "Axis 1", ylab = "Axis 2")
#biplot(its2_pcoa1, meta_qdat[,2:9],display = c("sites", "species"),type = c("text","points"), main="ITS1 PCoA + Scaled Chemical Data", xlab = "Axis 1", ylab = "Axis 2")

## bacteria + archaea next

bac_pcoa1 = pcoa(bac_RA_bray) # pcoa of Bray-Curtis dissimilarity matrix made w/ relativized abundance (site x species) data
# The proportion of variances explained is in its element values$Relative_eig
#cmdscale(bac.bray, k=(nrow(bac.ASV_table)-1), eig=TRUE)

biplot(bac_pcoa1)
#biplot(bac_pcoa1, meta_qdat[,2:9], main="16S PCoA + Scaled Chemical Data", xlab = "Axis 1", ylab = "Axis 2")
#biplot(bac_pcoa1, meta_qdat[,2:9],display = c("sites", "species"),type = c("text","points"))

#### Function for custom biplot ####
PCbiplot <- function(PC, x="PC1", y="PC2") {
  # PC being a prcomp object
  data <- data.frame(obsnames=row.names(PC$x), PC$x)
  plot <- ggplot(data, aes_string(x=x, y=y)) + geom_text(alpha=.4, size=4, aes(label=obsnames))
  plot <- plot + geom_hline(aes(0), size=.2) + geom_vline(aes(0), size=.2)
  datapc <- data.frame(varnames=rownames(PC$rotation), PC$rotation)
  mult <- min(
    (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
    (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
  )
  datapc <- transform(datapc,
                      v1 = .7 * mult * (get(x)),
                      v2 = .7 * mult * (get(y))
  )
  plot <- plot + coord_equal() + geom_text(data=datapc, aes(x=v1, y=v2, label=varnames), size = 5, vjust=1, color="red")
  plot <- plot + geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="red")
  plot
}

fit <- prcomp(USArrests, scale=T)
PCbiplot(fit)
#### ITS2 PCoA Visualization (relativized Bray-Curtis counts) ####

### fungi first

#tiff(file="figures/its2_pcoa_biplot.tiff",width=10, height=10, units="in", res=600)
biplot(its2_pcoa1)
dev.off()

# setting up the dataframe for ggplot2 visualization

its2_pcoa1.vectors<-data.frame(its2_pcoa1$vectors)
its2_pcoa1.vectors$SampleID<-rownames(its2_pcoa1.vectors)

its2_pcoa1_meta<-merge(its2_pcoa1.vectors, metadata, by.x="SampleID", by.y="SampleID")

its2_pcoa1_meta$runDuration <- factor(its2_pcoa1_meta$runDuration, levels = c("week","month"))

head(its2_pcoa1_meta)

# plot time


f.f1<-ggplot(its2_pcoa1_meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Exposed)), size=4)+theme_bw()+
  labs(title="PCoA: Fungi",subtitle="Using Relativized Bray-Curtis Dissimilarity",xlab="Axis 1", ylab="Axis 2",color="Exposure Material")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Exposure Material", labels=c("Alt"="Alternaria", "Con"="Control", "Silica"="Silica"),
                     values=unique(its2_pcoa1_meta$color[order(its2_pcoa1_meta$Exposed)])) +
  xlab("Axis 1") + ylab("Axis 2")
ggsave(f.f1,filename = "figures/its2_pcoa_6.20.21.pdf", width=10, height=6, dpi=600)

f.f1<-ggplot(its2_pcoa1_meta, aes(x=Axis.1, y=Axis.2, label=SampleID)) +geom_text(aes(color=factor(Exposed)), size=4)+theme_bw()+
  labs(title="PCoA: Fungi",subtitle="Using Relativized Bray-Curtis Dissimilarity",xlab="Axis 1", ylab="Axis 2",color="Exposure Material")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Exposure Material", labels=c("Alt"="Alternaria", "Con"="Control", "Silica"="Silica"),
                     values=unique(its2_pcoa1_meta$color[order(its2_pcoa1_meta$Exposed)])) +
  xlab("Axis 1") + ylab("Axis 2")
ggsave(f.f1,filename = "figures/its2_pcoa.labeled_6.20.21.pdf", width=10, height=6, dpi=600)

its2.pcoa.fig.1<-ggplot(its2_pcoa1_meta, aes(x=Axis.1, y=Axis.2, col=Exposed)) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(name = 'Exposure Material', labels = c('Alternaria', 'Control', 'Silica'),values =  saturation(SM_pal, 1))+
  labs(xlab="Axis 1", ylab="Axis 2",color="Exposure Material")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))
#coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))
ggsave(its2.pcoa.fig.1,filename = "figures/fungal_pcoa1_expos.mat_6.19.21.pdf", width=8, height=6, dpi=600)

its2.pcoa.fig.1a<-ggplot(its2_pcoa1_meta, aes(x=Axis.1, y=Axis.2, col=Exposed)) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(name = 'Exposure Material', labels = c('Alternaria', 'Control', 'Silica'),values =  saturation(SM_pal, 1))+
  labs(xlab="Axis 1", ylab="Axis 2",color="Exposure Material")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+stat_ellipse()
#coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))
ggsave(its2.pcoa.fig.1a,filename = "figures/fungal_pcoa1a_expos.mat_ellipses_6.19.21.pdf", width=8, height=6, dpi=600)

its2.pcoa.fig.2<-ggplot(its2_pcoa1_meta, aes(x=Axis.1, y=Axis.2, col=Exposed, shape=runDuration)) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(name = 'Exposure Material', labels = c('Alternaria', 'Control', 'Silica'),values =  saturation(SM_pal, 1))+
  labs(xlab="Axis 1", ylab="Axis 2",color="Exposure Material")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5))) + scale_shape_manual(name='Run Time', labels=c('Week', 'Month'), values=c(16,17))
#coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))
ggsave(its2.pcoa.fig.2,filename = "figures/fungal_pcoa2_exp.mat_runtime_6.19.21.pdf", width=8, height=6, dpi=600)

its2.pcoa.fig.3<-ggplot(its2_pcoa1_meta, aes(x=Axis.1, y=Axis.2, col=exposed_y.n)) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(name = 'Exposed Condition', labels = c('Y'='Yes','N'='No'),values =  saturation(Arc_pal, 1))+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Exposed Condition")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))

ggsave(its2.pcoa.fig.3,filename = "figures/fungal_pcoa3_expose.yes.no_6.19.21.pdf", width=8, height=6, dpi=600)

its2.pcoa.fig.4<-ggplot(its2_pcoa1_meta, aes(x=Axis.1, y=Axis.2, col=sampleType)) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(name = 'Sample Type', labels = c('DNA'='DNA','supernatant'='Supernatant', 'pellet'='Pellet', 'lungTissue'='Lung Tissue'),values =  saturation(wes5, 1))+
  labs(xlab="Axis 1", ylab="Axis 2",color="Sample Type")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))

ggsave(its2.pcoa.fig.4,filename = "figures/fungal_pcoa4_sample.type_6.19.21.pdf", width=8, height=6, dpi=600)

its2.pcoa.fig.4a<-ggplot(its2_pcoa1_meta, aes(x=Axis.1, y=Axis.2, col=sampleType, shape=runDuration)) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(name = 'Sample Type', labels = c('DNA'='DNA','supernatant'='Supernatant', 'pellet'='Pellet', 'lungTissue'='Lung Tissue'),values =  saturation(wes5, 1))+
  labs(xlab="Axis 1", ylab="Axis 2",color="Sample Type", shape="Run Time")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+ scale_shape_manual(name='Run Time', labels=c('Week', 'Month'), values=c(16,17))

ggsave(its2.pcoa.fig.4a,filename = "figures/fungal_pcoa4a_sample.type_runtime_6.19.21.pdf", width=8, height=6, dpi=600)

its2.pcoa.fig.4b<-ggplot(its2_pcoa1_meta, aes(x=Axis.1, y=Axis.2, col=sampleType, shape=animalBALFd)) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(name = 'Sample Type', labels = c('DNA'='DNA','supernatant'='Supernatant', 'pellet'='Pellet', 'lungTissue'='Lung Tissue'),values =  saturation(wes5, 1))+
  labs(xlab="Axis 1", ylab="Axis 2",color="Sample Type", shape="Run Time")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+ scale_shape_manual(name = 'DNA Extraction Material', labels = c('BALFdLung'='BALF-lung','na'='BALF', 'rawLung' = 'Lung'), values=c(15,16,17))

ggsave(its2.pcoa.fig.4b,filename = "figures/fungal_pcoa4b_sample.type_dna.ext.mat_6.19.21.pdf", width=8, height=6, dpi=600)

its2.pcoa.fig.5<-ggplot(its2_pcoa1_meta, aes(x=Axis.1, y=Axis.2, col=lungTissue, shape=exposed_y.n)) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(name = 'Lung Sample Type', labels = c('lung'='Lung','BALF'='BALF'),values =  saturation(Arc_pal, 1))+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Lung Sample Type", shape="Exposed Condition")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+ scale_shape_manual(name='Exposed Condition', labels=c('Y'='Yes', 'N'='No'), values=c(16,17))

ggsave(its2.pcoa.fig.5,filename = "figures/fungal_pcoa5_lung.tiss_expose.yes.no_6.19.21.pdf", width=8, height=6, dpi=600)

its2.pcoa.fig.6<-ggplot(its2_pcoa1_meta, aes(x=Axis.1, y=Axis.2, col=Exposed, shape=exposed_y.n)) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(name = 'Exposure Material', labels = c('Alternaria', 'Control', 'Silica'),values =  saturation(SM_pal, 1))+
  labs(xlab="Axis 1", ylab="Axis 2",color="Exposure Material", shape="Exposed Condition")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5))) + scale_shape_manual(name='Exposed Condition', labels=c('Y'='Yes', 'N'='No'), values=c(16,17))
#coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))
ggsave(its2.pcoa.fig.6,filename = "figures/fungal_pcoa6_exp.mat_exposed.yes.no_6.19.21.pdf", width=8, height=6, dpi=600)

its2.pcoa.fig.6a<-ggplot(its2_pcoa1_meta, aes(x=Axis.1, y=Axis.2, col=Exposed, shape=lungTissue)) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(name = 'Exposure Material', labels = c('Alternaria', 'Control', 'Silica'),values =  saturation(SM_pal, 1))+
  labs(xlab="Axis 1", ylab="Axis 2",color="Exposure Material", shape="Exposed Condition")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5))) + scale_shape_manual(name = 'Lung Sample Type', labels = c('lung'='Lung','BALF'='BALF'), values=c(16,17))
#coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))
ggsave(its2.pcoa.fig.6a,filename = "figures/fungal_pcoa6a_exp.mat_lung.tiss.type_6.19.21.pdf", width=8, height=6, dpi=600)

its2.pcoa.fig.7<-ggplot(its2_pcoa1_meta, aes(x=Axis.1, y=Axis.2, col=Exposed, shape=animalBALFd)) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(name = 'Exposure Material', labels = c('Alternaria', 'Control', 'Silica'),values =  saturation(SM_pal, 1))+
  labs(xlab="Axis 1", ylab="Axis 2",color="Exposure Material", shape="DNA Extraction Material")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5))) + scale_shape_manual(name = 'DNA Extraction Material', labels = c('BALFdLung'='BALF-lung','na'='BALF', 'rawLung' = 'Lung'), values=c(15,16,17))
#coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))
ggsave(its2.pcoa.fig.7,filename = "figures/fungal_pcoa7_exp.mat_dna.ext.mat_6.19.21.pdf", width=8, height=6, dpi=600)

its2.pcoa.fig.8<-ggplot(its2_pcoa1_meta, aes(x=Axis.1, y=Axis.2, col=animalBALFd, shape=runDuration)) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_shape_manual(name = 'Run Time', labels = c('Week', 'Month'), values=c(16,17))+
  labs(xlab="Axis 1", ylab="Axis 2",color="Exposure Material", shape="DNA Extraction Material")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5))) + scale_colour_manual(name = 'DNA Extraction Material', labels = c('BALFdLung'='BALF-lung','na'='BALF', 'rawLung' = 'Lung'),values =  saturation(wes2, 1))
#coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))
ggsave(its2.pcoa.fig.8,filename = "figures/fungal_pcoa8_dna.ext.mat_runtime_6.19.21.pdf", width=8, height=6, dpi=600)

#### 16S PCoA Visualization (relativized Bray-Curtis counts) ####


#tiff(file="figures/bac_pcoa_biplot.tiff",width=10, height=10, units="in", res=600)
biplot(bac_pcoa1)
dev.off()

# setting up the dataframe for ggplot2 visualization

bac_pcoa1.vectors<-data.frame(bac_pcoa1$vectors)
bac_pcoa1.vectors$SampleID<-rownames(bac_pcoa1.vectors)

bac_pcoa1_meta<-merge(bac_pcoa1.vectors, metadata, by.x="SampleID", by.y="SampleID")
bac_pcoa1_meta$runDuration <- factor(bac_pcoa1_meta$runDuration, levels = c("week","month"))

head(bac_pcoa1_meta)

# plot time

b.f1<-ggplot(bac_pcoa1_meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Exposed)), size=4)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea",subtitle="Using Relativized Bray-Curtis Dissimilarity",xlab="Axis 1", ylab="Axis 2",color="Exposure Material")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Exposure Material", labels=c("Alt"="Alternaria", "Con"="Control", "Silica"="Silica"),
                     values=unique(bac_pcoa1_meta$color[order(bac_pcoa1_meta$Exposed)])) +
  xlab("Axis 1") + ylab("Axis 2")
ggsave(b.f1,filename = "figures/bac_pcoa_6.20.21.pdf", width=10, height=6, dpi=600)

b.f1<-ggplot(bac_pcoa1_meta, aes(x=Axis.1, y=Axis.2, label=SampleID)) +geom_text(aes(color=factor(Exposed)), size=4)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea",subtitle="Using Relativized Bray-Curtis Dissimilarity",xlab="Axis 1", ylab="Axis 2",color="Exposure Material")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Exposure Material", labels=c("Alt"="Alternaria", "Con"="Control", "Silica"="Silica"),
                     values=unique(bac_pcoa1_meta$color[order(bac_pcoa1_meta$Exposed)])) +
  xlab("Axis 1") + ylab("Axis 2")
ggsave(b.f1,filename = "figures/bac_pcoa.labeled_6.20.21.pdf", width=10, height=6, dpi=600)

bac.pcoa.fig.1<-ggplot(bac_pcoa1_meta, aes(x=Axis.1, y=Axis.2, col=Exposed)) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(name = 'Exposure Material', labels = c('Alternaria', 'Control', 'Silica'),values =  saturation(SM_pal, 1))+
  labs(xlab="Axis 1", ylab="Axis 2",color="Exposure Material")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))
#coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))
ggsave(bac.pcoa.fig.1,filename = "figures/microb_pcoa1_expos.mat_6.19.21.pdf", width=8, height=6, dpi=600)

bac.pcoa.fig.1a<-ggplot(bac_pcoa1_meta, aes(x=Axis.1, y=Axis.2, col=Exposed)) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(name = 'Exposure Material', labels = c('Alternaria', 'Control', 'Silica'),values =  saturation(SM_pal, 1))+
  labs(xlab="Axis 1", ylab="Axis 2",color="Exposure Material")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+stat_ellipse()
#coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))
ggsave(bac.pcoa.fig.1a,filename = "figures/microb_pcoa1a_expos.mat_ellipses_6.19.21.pdf", width=8, height=6, dpi=600)

bac.pcoa.fig.2<-ggplot(bac_pcoa1_meta, aes(x=Axis.1, y=Axis.2, col=Exposed, shape=runDuration)) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(name = 'Exposure Material', labels = c('Alternaria', 'Control', 'Silica'),values =  saturation(SM_pal, 1))+
  labs(xlab="Axis 1", ylab="Axis 2",color="Exposure Material")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5))) + scale_shape_manual(name='Run Time', labels=c('Week', 'Month'), values=c(16,17))
#coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))
ggsave(bac.pcoa.fig.2,filename = "figures/microb_pcoa2_exp.mat_runDuration_6.19.21.pdf", width=8, height=6, dpi=600)

bac.pcoa.fig.3<-ggplot(bac_pcoa1_meta, aes(x=Axis.1, y=Axis.2, col=exposed_y.n)) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(name = 'Exposed Condition', labels = c('Y'='Yes','N'='No'),values =  saturation(Arc_pal, 1))+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Exposed Condition")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))

ggsave(bac.pcoa.fig.3,filename = "figures/microb_pcoa3_expose.yes.no_6.19.21.pdf", width=8, height=6, dpi=600)

bac.pcoa.fig.4<-ggplot(bac_pcoa1_meta, aes(x=Axis.1, y=Axis.2, col=sampleType)) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(name = 'Sample Type', labels = c('DNA'='DNA','supernatant'='Supernatant', 'pellet'='Pellet', 'lungTissue'='Lung Tissue'),values =  saturation(wes5, 1))+
  labs(xlab="Axis 1", ylab="Axis 2",color="Sample Type")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))

ggsave(bac.pcoa.fig.4,filename = "figures/microb_pcoa4_sample.type_6.19.21.pdf", width=8, height=6, dpi=600)

bac.pcoa.fig.4a<-ggplot(bac_pcoa1_meta, aes(x=Axis.1, y=Axis.2, col=sampleType, shape=runDuration)) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(name = 'Sample Type', labels = c('DNA'='DNA','supernatant'='Supernatant', 'pellet'='Pellet', 'lungTissue'='Lung Tissue'),values =  saturation(wes5, 1))+
  labs(xlab="Axis 1", ylab="Axis 2",color="Sample Type", shape="Run Time")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+ scale_shape_manual(name='Run Time', labels=c('Week', 'Month'), values=c(16,17))

ggsave(bac.pcoa.fig.4a,filename = "figures/microb_pcoa4a_sample.type_runDuration_6.19.21.pdf", width=8, height=6, dpi=600)

bac.pcoa.fig.4b<-ggplot(bac_pcoa1_meta, aes(x=Axis.1, y=Axis.2, col=sampleType, shape=animalBALFd)) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(name = 'Sample Type', labels = c('DNA'='DNA','supernatant'='Supernatant', 'pellet'='Pellet', 'lungTissue'='Lung Tissue'),values =  saturation(wes5, 1))+
  labs(xlab="Axis 1", ylab="Axis 2",color="Sample Type", shape="Run Time")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+ scale_shape_manual(name = 'DNA Extraction Material', labels = c('BALFdLung'='BALF-lung','na'='BALF', 'rawLung' = 'Lung'), values=c(15,16,17))

ggsave(bac.pcoa.fig.4b,filename = "figures/microb_pcoa4b_sample.type_dna.ext.mat_6.19.21.pdf", width=8, height=6, dpi=600)

bac.pcoa.fig.5<-ggplot(bac_pcoa1_meta, aes(x=Axis.1, y=Axis.2, col=lungTissue, shape=exposed_y.n)) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(name = 'Lung Sample Type', labels = c('lung'='Lung','BALF'='BALF'),values =  saturation(Arc_pal, 1))+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Lung Sample Type", shape="Exposed Condition")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+ scale_shape_manual(name='Exposed Condition', labels=c('Y'='Yes', 'N'='No'), values=c(16,17))

ggsave(bac.pcoa.fig.5,filename = "figures/microb_pcoa5_lung.tiss_expose.yes.no_6.19.21.pdf", width=8, height=6, dpi=600)

bac.pcoa.fig.6<-ggplot(bac_pcoa1_meta, aes(x=Axis.1, y=Axis.2, col=Exposed, shape=exposed_y.n)) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(name = 'Exposure Material', labels = c('Alternaria', 'Control', 'Silica'),values =  saturation(SM_pal, 1))+
  labs(xlab="Axis 1", ylab="Axis 2",color="Exposure Material", shape="Exposed Condition")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5))) + scale_shape_manual(name='Exposed Condition', labels=c('Y'='Yes', 'N'='No'), values=c(16,17))
#coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))
ggsave(bac.pcoa.fig.6,filename = "figures/microb_pcoa6_exp.mat_exposed.yes.no_6.19.21.pdf", width=8, height=6, dpi=600)

bac.pcoa.fig.6a<-ggplot(bac_pcoa1_meta, aes(x=Axis.1, y=Axis.2, col=Exposed, shape=lungTissue)) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(name = 'Exposure Material', labels = c('Alternaria', 'Control', 'Silica'),values =  saturation(SM_pal, 1))+
  labs(xlab="Axis 1", ylab="Axis 2",color="Exposure Material", shape="Exposed Condition")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5))) + scale_shape_manual(name = 'Lung Sample Type', labels = c('lung'='Lung','BALF'='BALF'), values=c(16,17))
#coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))
ggsave(bac.pcoa.fig.6a,filename = "figures/microb_pcoa6a_exp.mat_lung.tiss.type_6.19.21.pdf", width=8, height=6, dpi=600)

bac.pcoa.fig.7<-ggplot(bac_pcoa1_meta, aes(x=Axis.1, y=Axis.2, col=Exposed, shape=animalBALFd)) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(name = 'Exposure Material', labels = c('Alternaria', 'Control', 'Silica'),values =  saturation(SM_pal, 1))+
  labs(xlab="Axis 1", ylab="Axis 2",color="Exposure Material", shape="DNA Extraction Material")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5))) + scale_shape_manual(name = 'DNA Extraction Material', labels = c('BALFdLung'='BALF-lung','na'='BALF', 'rawLung' = 'Lung'), values=c(15,16,17))
#coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))
ggsave(bac.pcoa.fig.7,filename = "figures/microb_pcoa7_exp.mat_dna.ext.mat_6.19.21.pdf", width=8, height=6, dpi=600)

bac.pcoa.fig.8<-ggplot(bac_pcoa1_meta, aes(x=Axis.1, y=Axis.2, col=animalBALFd, shape=runDuration)) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_shape_manual(name = 'Run Time', labels = c('Week', 'Month'), values=c(16,17))+
  labs(xlab="Axis 1", ylab="Axis 2",color="Exposure Material", shape="DNA Extraction Material")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5))) + scale_colour_manual(name = 'DNA Extraction Material', labels = c('BALFdLung'='BALF-lung','na'='BALF', 'rawLung' = 'Lung'),values =  saturation(wes2, 1))
#coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))
ggsave(bac.pcoa.fig.8,filename = "figures/microb_pcoa8_dna.ext.mat_runDuration_6.19.21.pdf", width=8, height=6, dpi=600)


#### Prep for PERMANOVAs ####
head(metadata)
head(its2.RA_table)
head(bac.RA_table)

metadata<-subset(metadata, select=c(Exposed, sampleType, exposed_y.n, lungTissue, balfSup, runDuration, animalBALFd))

## reorder dataframes of metadata and tables by row so all dfs have same order
# * must have same rownames and same # of rows to do this!
metadata=metadata[rownames(its2.RA_table),]
bac.RA_table=bac.RA_table[rownames(its2.RA_table),]

samp.type<-metadata$sampleType
exp.mat<-metadata$Exposed
lung.tissue<-metadata$lungTissue
exp.yes.no<-metadata$exposed_y.n
run.time<-metadata$runDuration
balf.type<-metadata$animalBALFd
# BALFDlung - lung from animal that endured BALF; BALF - combining BALF pellet + supernatant together; raw lung - lung from animal that did NOT endure BALF
samp<-metadata$balfSup

#### ITS2 PERMANOVAs ####

## The currently preferred analysis for evaluating differences among groups is PERMANOVA.
## This analysis partitions sums of squares using dissimilarities,
##  evaluating differences in the centroids of groups in multivariate space.
##  The vegan functions “adonis” and “adonis2” are used to compute PERMANOVA in R.

help(adonis)

# Fungi first

##  The following models are equivalent:
#adonis(its2.RA_table ~ prov,  permutations = 999, method = "bray")
adonis2(its2.RA_table ~ exp.mat, permutations = 999, method = "bray")
adonis2(its2.RA_table ~ lung.tissue, permutations = 999, method = "bray")
adonis2(its2.RA_table ~ exp.yes.no, permutations = 999, method = "bray")
adonis2(its2.RA_table ~ run.time, permutations = 999, method = "bray")
adonis2(its2.RA_table ~ balf.type, permutations = 999, method = "bray")
adonis2(its2.RA_table ~ samp, permutations = 999, method = "bray")
adonis2(its2.RA_table ~ samp.type, permutations = 999, method = "bray")

## The above models specify dataframes for analysis, but we can alternatively specify a dissimilarity matrix:
#adonis2(compDS ~ soil, permutations = 999)

#Other advantages of using PERMANOVA are that we can test for interactions between predictor variables,
## and we can use both categorical and continuous predictor variables.
adonis2(its2.RA_table ~ Exposed*lungTissue, data = metadata, permutations = 999, method = "bray", by='terms') # looks for interactions between predictor variables
#These models provide the significance for statistical interactions and the main effects.
## An advantage of adonis2 is that we can also test for overall model fit, using the “by” command:
adonis2(its2.RA_table ~ Exposed*lungTissue, data = metadata, permutations = 999, method = "bray", by = NULL)

adonis2(its2.RA_table ~ lungTissue*sampleType, data = metadata, permutations = 999, method = "bray", by='terms') # looks for interactions between predictor variables
adonis2(its2.RA_table ~ runDuration*lungTissue, data = metadata, permutations = 999, method = "bray", by='terms') # looks for interactions between predictor variables
adonis2(its2.RA_table ~ Exposed*balfSup, data = metadata, permutations = 999, method = "bray", by='terms') # looks for interactions between predictor variables

adonis2(its2.RA_table ~ Exposed*lungTissue*sampleType*balfSup*runDuration*animalBALFd, data = metadata, permutations = 999, method = "bray", by='terms') # looks for interactions between predictor variables

#While PERMANOVA tests differences in group means (analogous to MANOVA),
## a related test called PERMDISP can be used to evaluate homogeneity of group dispersion
#(analogous to Levene's test for equal variances). The vegan function for this test is “betadisper”:
its2.disper <- betadisper(its2.RA_bray, metadata$Exposed)
its2.disper

## Significant differences in homogeneities can be tested using either parametric or permutational tests,
##and parametric post hoc contrasts can also be investigated:
anova(its2.disper)
permutest(its2.disper)
TukeyHSD(its2.disper)

##one issue with adonis is that it doesn't do multiple comparisons *******
# tells us that something is different, but what is different? Which sample/plot/location?
## our four provinces differ, but do all of them differ,or just one?

##random person on the internet to the rescue!
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

pair.mod<-pairwise.adonis(its2.RA_table,exp.mat, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different
pair.mod

pair.mod1<-pairwise.adonis(its2.RA_table,lung.tissue, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different
pair.mod1

pair.mod2<-pairwise.adonis(its2.RA_table,balf.type, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different
pair.mod2

pair.mod3<-pairwise.adonis(its2.RA_table,metadata$balfSup, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different
pair.mod3

### SELF REMINDER FOR R^2
### Coefficient of Determination, denoted R2 or r2 and pronounced "R squared"
### is the proportion of the variance in the dependent variable that is predictable from the independent variable(s)

### Pseudo F stat for PERMANOVA
### pseudo F-ratio: It compares the total sum of squared dissimilarities (or ranked dissimilarities) among objects belonging to different groups to that of objects belonging to the same group.
### Larger F-ratios indicate more pronounced group separation, however, the significance of this ratio is usually of more interest than its magnitude.

#### 16S PERMANOVAs ####

## The currently preferred analysis for evaluating differences among groups is PERMANOVA.
## This analysis partitions sums of squares using dissimilarities,
##  evaluating differences in the centroids of groups in multivariate space.
##  The vegan functions “adonis” and “adonis2” are used to compute PERMANOVA in R.

help(adonis)

##  The following models are equivalent:
#adonis(bac.RA_table ~ prov,  permutations = 999, method = "bray")
adonis2(bac.RA_table ~ exp.mat, permutations = 999, method = "bray")
adonis2(bac.RA_table ~ lung.tissue, permutations = 999, method = "bray")
adonis2(bac.RA_table ~ exp.yes.no, permutations = 999, method = "bray")
adonis2(bac.RA_table ~ run.time, permutations = 999, method = "bray")
adonis2(bac.RA_table ~ balf.type, permutations = 999, method = "bray")
adonis2(bac.RA_table ~ samp, permutations = 999, method = "bray")
adonis2(bac.RA_table ~ samp.type, permutations = 999, method = "bray")

## The above models specify dataframes for analysis, but we can alternatively specify a dissimilarity matrix:
#adonis2(compDS ~ soil, permutations = 999)

#Other advantages of using PERMANOVA are that we can test for interactions between predictor variables,
## and we can use both categorical and continuous predictor variables.
adonis2(bac.RA_table ~ Exposed*lungTissue, data = metadata, permutations = 999, method = "bray", by='terms') # looks for interactions between predictor variables
#These models provide the significance for statistical interactions and the main effects.
## An advantage of adonis2 is that we can also test for overall model fit, setting by=NULL

adonis2(bac.RA_table ~ lungTissue*sampleType, data = metadata, permutations = 999, method = "bray", by='terms') # looks for interactions between predictor variables
adonis2(bac.RA_table ~ lungTissue*sampleType, data = metadata, permutations = 999, method = "bray", by=NULL) # looks for interactions between predictor variables

adonis2(bac.RA_table ~ runDuration*lungTissue, data = metadata, permutations = 999, method = "bray", by='terms') # looks for interactions between predictor variables
adonis2(bac.RA_table ~ runDuration*sampleType, data = metadata, permutations = 999, method = "bray", by='terms') # looks for interactions between predictor variables

adonis2(bac.RA_table ~ Exposed*balfSup, data = metadata, permutations = 999, method = "bray", by='terms') # looks for interactions between predictor variables

adonis2(bac.RA_table ~ Exposed*lungTissue*sampleType*balfSup*runDuration*animalBALFd, data = metadata, permutations = 999, method = "bray", by='terms') # looks for interactions between predictor variables

#While PERMANOVA tests differences in group means (analogous to MANOVA),
## a related test called PERMDISP can be used to evaluate homogeneity of group dispersion
#(analogous to Levene's test for equal variances). The vegan function for this test is “betadisper”:
bac.disper <- betadisper(bac.RA_bray, metadata$Exposed)
bac.disper

## Significant differences in homogeneities can be tested using either parametric or permutational tests,
##and parametric post hoc contrasts can also be investigated:
anova(bac.disper)
permutest(bac.disper)
TukeyHSD(bac.disper)

##one issue with adonis is that it doesn't do multiple comparisons *******
# tells us that something is different, but what is different? Which sample/plot/location?
## our four provinces differ, but do all of them differ,or just one?

##random person on the internet to the rescue!
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

pair.mod<-pairwise.adonis(bac.RA_table,exp.mat, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different
pair.mod

pair.mod1<-pairwise.adonis(bac.RA_table,lung.tissue, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different
pair.mod1

pair.mod2<-pairwise.adonis(bac.RA_table,balf.type, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different
pair.mod2

pair.mod3<-pairwise.adonis(bac.RA_table,metadata$balfSup, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different
pair.mod3

pair.mod4<-pairwise.adonis(bac.RA_table,metadata$sampleType, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different
pair.mod4

### SELF REMINDER FOR R^2
### Coefficient of Determination, denoted R2 or r2 and pronounced "R squared"
### is the proportion of the variance in the dependent variable that is predictable from the independent variable(s)

### Pseudo F stat for PERMANOVA
### pseudo F-ratio: It compares the total sum of squared dissimilarities (or ranked dissimilarities) among objects belonging to different groups to that of objects belonging to the same group.
### Larger F-ratios indicate more pronounced group separation, however, the significance of this ratio is usually of more interest than its magnitude.


#### Jaccard Distance Matrices ####
## using binary tables for these Jaccard

## fungi first
class(its2_binary_table) # raw counts
dim(its2_binary_table)

its2_jac<-vegdist(its2_binary_table,method="jaccard")     ####### distance matrix with RELATIVIZED data!!!
# **** in vegan: ROWS need to be SITES/samples; COLUMNS are SPECIES (OTUs, ASVs) -- vegdist!
class(its2_jac) # needs to be class dist for pcoa()

## bacterial/archaeal next

class(bac_binary_table) # raw counts
dim(bac_binary_table)

bac_jac<-vegdist(bac_binary_table,method="jaccard")     ####### distance matrix with RELATIVIZED data!!!
# **** in vegan: ROWS need to be SITES/samples; COLUMNS are SPECIES (OTUs, ASVs) -- vegdist!
class(bac_jac) # needs to be class dist for pcoa()



#### PCoA -- w/ Jaccard matrix (from binary counts) ####

## fungi first

its2_pcoa2 = pcoa(its2_jac) # pcoa of Bray-Curtis dissimilarity matrix made w/ relativized abundance (site x species) data

biplot(its2_pcoa2, scale(meta_quant))

its2_pcoa2.vectors<-data.frame(its2_pcoa2$vectors)
its2_pcoa2.vectors$SampleID<-rownames(its2_pcoa2.vectors)

its2_pcoa2_meta<-merge(its2_pcoa2.vectors, metadata, by="SampleID")
head(its2_pcoa2_meta)

its2_pcoa2_df=data.frame(SampleID = its2_pcoa2_meta$SampleID, Axis1 = its2_pcoa2_meta$Axis.1, Axis2 = its2_pcoa2_meta$Axis.2, Site = its2_pcoa2_meta$Site, Elevation = its2_pcoa2_meta$Elevation, Month=its2_pcoa2_meta$Month, Year= its2_pcoa2_meta$Year) #create dataframe with PCoA axes and some metadata
head(its2_pcoa2_df)

its2_pcoa2_df$Month2 <- factor(its2_pcoa2_df$Month, levels = c("July","August","October"))
its2_pcoa2_df$Elevation2 <- factor(its2_pcoa2_df$Elevation, levels = c("400","1100","2000","2700"))

head(its2_pcoa2_df)

its2.pcoa2.fig.1<-ggplot(its2_pcoa2_df, aes(x=Axis1, y=Axis2, col=Month2, shape=factor(Year))) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(values =  saturation(SM_pal, 1))+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Month", shape="Year")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))

ggsave(its2.pcoa2.fig.1,filename = "figures/fungal_pcoa1_Jaccard_Sierra_3.29.21.pdf", width=8, height=6, dpi=600)

its2.pcoa2.fig.2<-ggplot(its2_pcoa2_df, aes(x=Axis1, y=Axis2, col=Month2, shape=factor(Year), size=Elevation2)) +geom_point(alpha=0.5)+theme_bw()+scale_colour_manual(values =  saturation(SM_pal, 1))+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Month", shape="Year",size="Elevation")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(color = guide_legend(override.aes = list(size = 5)),shape = guide_legend(override.aes = list(size = 5)))

ggsave(its2.pcoa2.fig.2,filename = "figures/fungal_pcoa2_Jaccard_Sierra_3.29.21.pdf", width=8, height=6, dpi=600)

## bacteria + archaea next

bac_pcoa2 = pcoa(bac_jac) # pcoa of Bray-Curtis dissimilarity matrix made w/ relativized abundance (site x species) data

biplot(bac_pcoa2, scale(meta_quant))

bac_pcoa2.vectors<-data.frame(bac_pcoa2$vectors)
bac_pcoa2.vectors$SampleID<-rownames(bac_pcoa2.vectors)

bac_pcoa2_meta<-merge(bac_pcoa2.vectors, metadata, by="SampleID")
head(bac_pcoa2_meta)

bac_pcoa2_df=data.frame(SampleID = bac_pcoa2_meta$SampleID, Axis1 = bac_pcoa2_meta$Axis.1, Axis2 = bac_pcoa2_meta$Axis.2, Site = bac_pcoa2_meta$Site, Elevation = bac_pcoa2_meta$Elevation, Month=bac_pcoa2_meta$Month, Year= bac_pcoa2_meta$Year) #create dataframe with PCoA axes and some metadata
head(bac_pcoa2_df)

bac_pcoa2_df$Month2 <- factor(bac_pcoa2_df$Month, levels = c("July","August","October"))
bac_pcoa2_df$Elevation2 <- factor(bac_pcoa2_df$Elevation, levels = c("400","1100","2000","2700"))

head(bac_pcoa2_df)

bac.pcoa2.fig.1<-ggplot(bac_pcoa2_df, aes(x=Axis1, y=Axis2, col=Month2, shape=factor(Year))) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(values =  saturation(SM_pal, 1))+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Month", shape="Year")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))

ggsave(bac.pcoa2.fig.1,filename = "figures/microbial_pcoa1_Jaccard_Sierra_3.29.21.pdf", width=8, height=6, dpi=600)

bac.pcoa2.fig.2<-ggplot(bac_pcoa2_df, aes(x=Axis1, y=Axis2, col=Month2, shape=factor(Year), size=Elevation2)) +geom_point(alpha=0.5)+theme_bw()+scale_colour_manual(values =  saturation(SM_pal, 1))+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Month", shape="Year",size="Elevation")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(color = guide_legend(override.aes = list(size = 5)),shape = guide_legend(override.aes = list(size = 5)))

ggsave(bac.pcoa2.fig.2,filename = "figures/microbial_pcoa2_Jaccard_Sierra_3.29.21.pdf", width=8, height=6, dpi=600)





#### PCoA ####

# PCoA produces a set of uncorrelated (orthogonal) axes to summarise the variability in the data set.
# Each axis has an eigenvalue whose magnitude indicates the amount of variation captured in that axis
# The proportion of a given eigenvalue to the sum of all eigenvalues reveals the relative ‚importance‘ of each axis.
# A successful PCoA will generate a few (2-3) axes with relatively large eigenvalues, capturing above 50% of the variation in the input data, with all other axes having small eigenvalues
# Each object has a ’score‘ along each axis. The object scores provide the object coordinates in the ordination plot
# ^^^ source: https://archetypalecology.wordpress.com/2018/02/19/principal-coordinates-analysis-pcoa-in-r/

## Let’s compute a matrix of Bray-Curtis similarities among sites, and subject this matrix to PCoA.
## If the metric is non-euclidean (as in our case), then the PCoA may produce several negative eigenvalues in
## addition to the positive ones. In most applications, this does not affect the representation of the first
## several axes. You will still receive a warning message, though, when this occurs. You also will receive a
## warning about species scores not being available; there is a way to project weighted averages of species
## abundances on a PCoA plot using the function wascores, we will do that too.

PCOA = pcoa(BrayC_ASV) # pcoa of Bray-Curtis dissimilarity matrix made w/ relativized abundance (site x species) data

## Check to see if negative eigenvalues affect the interpretation of the first several axes
PCOA$values

biplot(PCOA)

## project species on the PCoA ordination
biplot(PCOA, metadata$Type)

## project soils on the PCoA ordination
biplot(PCOA, soils)

## so just CA and MG? or should we scale our variables
biplot(PCOA, scale(soils))

## these look a little better, but you sill wouldn't want to publish this
## I personally use a different program, but if you want to do it in R I
## suggest you continue to explore ggplot
PCOA$values # Eigenvalues
PCOA$vectors # Eigenvectors
PCOA$correction # No correction used

## unlike NMDS, we can use PCoA sores for other analyses and
str(PCOA)
PCoAscores=as.data.frame(PCOA$vectors[,1:2])

## plot to check
plot(PCoAscores$Axis.2 ~ PCoAscores$Axis.1 )
## these can be used as a predictor variable in other analyses. NMDS scores cannot!!!

PCOA.vectors<-data.frame(PCOA$vectors)

pcoa.df = data.frame(Axis1 = PCOA.vectors$Axis.1, Axis2 = PCOA.vectors$Axis.2, SampleID = sample_names$SampleID) #create dataframe with PCoA axes and some metadata
pcoa.df
pcoa.meta<-merge(pcoa.df, metadata, by="SampleID")
pcoa.meta$Type2 <- factor(pcoa.meta$Type, levels = c("air","soil","water","control"))

pcoa.1<-ggplot(pcoa.meta, aes(x=Axis1, y=Axis2, col=Type2)) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(values =saturation(ewf, 1), name ="Sample Type", labels=c("air"="Aeolian", "soil"="Soil", "water"="Water", "control"="Control"))+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(title="PCoA Based on Relativized, Bray-Curtis Dissimilarity Index",xlab="Axis 1", ylab="Axis 2",color="Sample Type")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11), plot.title=element_text(hjust=0.5, size=15))+
  guides(shape = guide_legend(override.aes = list(size = 5)))

ggsave(pcoa.1,filename = "bacteria_pcoa_SaltonSea_wArchaea_6.16.21.pdf", width=10, height=10, dpi=600)

pcoa.2<-ggplot(pcoa.meta, aes(x=Axis1, y=Axis2, col=Month2, shape=factor(Year), size=Elevation2)) +geom_point(alpha=0.5)+theme_bw()+scale_colour_manual(values =  saturation(SM_pal, 1))+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Month", shape="Year",size="Elevation")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(color = guide_legend(override.aes = list(size = 5)),shape = guide_legend(override.aes = list(size = 5)))

ggsave(pcoa.2,filename = "bacteria_pcoa.2_Sierra_wArchaea_9.13.2020.pdf", width=8, height=6, dpi=600)


# pcoa biplot with log transformed environmental variables!
metadat_env<-data.frame(subset(metadata, select=-c(NdEp, Year, X87Sr86Sr, Elevation, SampleID, Month, Site, DateCode)))
metadat_env=metadat_env[rownames(metadata),] ## reorder metadata.dnv to match metadata df, which matches same order as original OTU table
head(metadat_env)
metadat_env.n<-sapply(metadat_env, as.numeric) # turn variables into numeric variables
rownames(metadat_env.n)<-rownames(metadat_env)
metadat_env.scale<-scale(metadat_env)
#metadat_env.log<-data.frame(decostand(metadat_env.n, method="log", MARGIN=1, na.rm=TRUE))  # log transformation of environmental data
#metadat_env.sqrt<-sqrt(metadat_env.n) # square root transformation of environmental data
#metadat_env.norm<-data.frame(decostand(metadat_env.n, method="normalize", MARGIN=1, na.rm=TRUE))  # standardized environmental data

biplot.pcoa(otu.braysq.pcoa, scale(metadat_env), plot.axes = c(1, 2)) # PCOA w/ sq rt transformed Bray-curtis dissimilarity matrix (fron OTU table) + log transformed environmental data

biplot.pcoa(otu.braysq.pcoa, metadat_env.scale, plot.axes = c(1, 2)) # PCOA w/ sq rt transformed Bray-curtis dissimilarity matrix (fron OTU table) + log transformed environmental data

biplot.pcoa(otu.braysq.pcoa, metadat_env.n, plot.axes = c(1, 2)) # PCOA w/ sq rt transformed Bray-curtis dissimilarity matrix (fron OTU table) + raw environmental data

biplot.pcoa(otu.braysq.pcoa, metadat_env.sqrt, plot.axes = c(1, 2)) # PCOA w/ sq rt transformed Bray-curtis dissimilarity matrix (fron OTU table) + square-root transformed environmental data

biplot.pcoa(otu.braysq.pcoa, metadat_env.norm, plot.axes = c(1, 2)) # PCOA w/ sq rt transformed Bray-curtis dissimilarity matrix (fron OTU table) + normalized environmental data

biplot.pcoa(otu.braysq.pcoa, metadat_env.log, plot.axes = c(1, 2)) # PCOA w/ sq rt transformed Bray-curtis dissimilarity matrix (fron OTU table) + log transformed environmental data

metadat_env.2<-data.frame(subset(metadata, select=c(NdEp, X87Sr86Sr))) ### ******* PLOT COMES OUT WEIRD 9.10
metadat_env.2=metadat_env.2[rownames(metadata),] ## reorder metadata.dnv to match metadata df, which matches same order as original OTU table

metadat_env.n.2<-sapply(metadat_env.2, as.numeric) # turn variables into numeric variables
rownames(metadat_env.n.2)<-rownames(metadat_env.2)
metadat_env.2.log<-data.frame(decostand(metadat_env.n.2, method="log", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows, which are SAMPLES in this case)

biplot.pcoa(otu.braysq.pcoa, metadat_env.n.2, plot.axes = c(1, 2)) # PCOA w/ sq rt transformed Bray-curtis dissimilarity matrix (fron OTU table) + log transformed environmental data

# if you decide to plot with biplot, rerun the pcoa and skip the line above ^^ biplot needs class dist!

## PCOA with Relative Abundances of OTUs (RElAb transformed data) ##

otu.bray.RA.pcoa<-pcoa(bray.otu.RelAb) # PCoA command!
class(otu.bray.RA.pcoa)

otu.bray.RA.pcoa$values # Eigenvalues
otu.bray.RA.pcoa$vectors # Eigenvectors
otu.bray.RA.pcoa$correction # No correction used

otu.bray.RA.pcoa.vectors = data.frame(otu.bray.RA.pcoa$vectors) #turn PCoA vectors into a dataframe so that you can access the Axes

head(metadata)

pcoa.meta_RA = data.frame(Axis1 = otu.bray.RA.pcoa.vectors$Axis.1, Axis2 = otu.bray.RA.pcoa.vectors$Axis.2, Site = metadata$Site,Elevation = metadata$Elevation,Month=metadata$Month, Year= metadata$Year) #create dataframe with PCoA axes and some metadata
pcoa.meta_RA
pcoa.meta_RA$Month2 <- factor(pcoa.meta_RA$Month, levels = c("July","August","October"))
pcoa.meta_RA$Elevation2 <- factor(pcoa.meta_RA$Elevation, levels = c("400","1100","2000","2700"))


pcoa.RA.1<-ggplot(pcoa.meta_RA, aes(x=Axis1, y=Axis2, col=Month2, shape=factor(Year))) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(values =  saturation(SM_pal, 1))+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Month", shape="Year")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))
ggsave(pcoa.RA.1,filename = "bacteria_pcoa_Sierra_wArchaea_RelativeAbundance_9.13.2020.pdf", width=8, height=6, dpi=600)

pcoa.RA.2<-ggplot(pcoa.meta_RA, aes(x=Axis1, y=Axis2, col=Month2, shape=factor(Year), size=Elevation2)) +geom_point(alpha=0.5)+theme_bw()+scale_colour_manual(values =  saturation(SM_pal, 1))+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Month", shape="Year",size="Elevation")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(color = guide_legend(override.aes = list(size = 5)),shape = guide_legend(override.aes = list(size = 5)))
ggsave(pcoa.RA.2,filename = "bacteria_pcoa_Sierra_wArchaea_RelativeAbundance_2_9.13.2020.pdf", width=8, height=6, dpi=600)


