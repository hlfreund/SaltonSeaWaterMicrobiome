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

# load Rdata to global env
#save.image("data/Env_Seqs_All/env.seq_analysis.Rdata") # save global env to Rdata file

bac.ASV_table[1:4,1:4]
bac.ASV_table[(nrow(bac.ASV_table)-4):(nrow(bac.ASV_table)),(ncol(bac.ASV_table)-4):(ncol(bac.ASV_table))] # last 4 rows & cols

#### 16S: Relative Abundance by Taxa Level ####

### * decostand(df, method="total) is the function (with argument total) used to get relative abundance of OTU table

## ASV_ID...
bac.ASV_table[1:4,1:4]
b.ASV_RelAb<-data.frame(decostand(b.ASV_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(b.ASV_RelAb) # sanity check to make sure the transformation worked!
b.ASV_RelAb$SampleID<-rownames(b.ASV_RelAb)
head(b.ASV_RelAb)
b.ASV_m<-melt(b.ASV_RelAb)

head(b.ASV_m)
colnames(b.ASV_m)[which(names(b.ASV_m) == "variable")] <- "ASV"
colnames(b.ASV_m)[which(names(b.ASV_m) == "value")] <- "Count"
head(b.ASV_m) ## relative abundance based on sum of counts by phyla!
#b.ASV_m<-subset(b.ASV_m, SampleID!="Undetermined_S0") ## keep only bacteria -- exclude Chloroplast sequences

b.ASV_RA_meta<-merge(b.ASV_m,metadata, by="SampleID")
head(b.ASV_RA_meta) ## relative abundance based on sum of counts by phyla!

## phylum ....
b.phyla_counts <- as.data.frame(dcast(b.dat.m, SampleID~Phylum, value.var="Count", fun.aggregate=sum)) ###
head(b.phyla_counts) # counts by phyla per sample
rownames(b.phyla_counts)<-b.phyla_counts$SampleID
b.phyla_counts$SampleID<-NULL
b.phyla_RelAb<-data.frame(decostand(b.phyla_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(b.phyla_RelAb) # sanity check to make sure the transformation worked!
b.phyla_RelAb$SampleID<-rownames(b.phyla_RelAb)
head(b.phyla_RelAb)
b.phyla_m<-melt(b.phyla_RelAb)

head(b.phyla_m)
colnames(b.phyla_m)[which(names(b.phyla_m) == "variable")] <- "Phylum"
colnames(b.phyla_m)[which(names(b.phyla_m) == "value")] <- "Count"
head(b.phyla_m) ## relative abundance based on sum of counts by phyla!
b.phyla_m<-subset(b.phyla_m, SampleID!="Undetermined_S0") ## keep only bacteria -- exclude Chloroplast sequences

b.phyla_RA_meta<-merge(b.phyla_m,metadata, by="SampleID")
head(b.phyla_RA_meta) ## relative abundance based on sum of counts by phyla!

## class...

b.class_counts <- as.data.frame(dcast(b.dat.m, SampleID~Class, value.var="Count", fun.aggregate=sum)) ###
head(b.class_counts) # counts by class per sample
rownames(b.class_counts)<-b.class_counts$SampleID
b.class_counts$SampleID<-NULL
b.class_RelAb<-data.frame(decostand(b.class_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(b.class_RelAb) # sanity check to make sure the transformation worked!
b.class_RelAb$SampleID<-rownames(b.class_RelAb)
head(b.class_RelAb)
b.class_m<-melt(b.class_RelAb)

head(b.class_m)
colnames(b.class_m)[which(names(b.class_m) == "variable")] <- "Class"
colnames(b.class_m)[which(names(b.class_m) == "value")] <- "Count"
head(b.class_m) ## relative abundance based on sum of counts by class!
b.class_m<-subset(b.class_m, SampleID!="Undetermined_S0") ## keep only bacteria -- exclude Chloroplast sequences

b.class_RA_meta<-merge(b.class_m,metadata, by="SampleID")
head(b.class_RA_meta) ## relative abundance based on sum of counts by class!

## order ...

b.order_counts <- as.data.frame(dcast(b.dat.m, SampleID~Order, value.var="Count", fun.aggregate=sum)) ###
head(b.order_counts) # counts by order per sample
rownames(b.order_counts)<-b.order_counts$SampleID
b.order_counts$SampleID<-NULL
b.order_RelAb<-data.frame(decostand(b.order_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(b.order_RelAb) # sanity check to make sure the transformation worked!
b.order_RelAb$SampleID<-rownames(b.order_RelAb)
head(b.order_RelAb)
b.order_m<-melt(b.order_RelAb)

head(b.order_m)
colnames(b.order_m)[which(names(b.order_m) == "variable")] <- "Order"
colnames(b.order_m)[which(names(b.order_m) == "value")] <- "Count"
head(b.order_m) ## relative abundance based on sum of counts by order!
b.order_m<-subset(b.order_m, SampleID!="Undetermined_S0") ## keep only bacteria -- exclude Chloroplast sequences

b.order_RA_meta<-merge(b.order_m,metadata, by="SampleID")
head(b.order_RA_meta) ## relative abundance based on sum of counts by order!

## family ...

b.family_counts <- as.data.frame(dcast(b.dat.m, SampleID~Family, value.var="Count", fun.aggregate=sum)) ###
head(b.family_counts) # counts by family per sample
rownames(b.family_counts)<-b.family_counts$SampleID
b.family_counts$SampleID<-NULL
b.family_RelAb<-data.frame(decostand(b.family_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(b.family_RelAb) # sanity check to make sure the transformation worked!
b.family_RelAb$SampleID<-rownames(b.family_RelAb)
head(b.family_RelAb)
b.family_m<-melt(b.family_RelAb)

head(b.family_m)
colnames(b.family_m)[which(names(b.family_m) == "variable")] <- "Family"
colnames(b.family_m)[which(names(b.family_m) == "value")] <- "Count"
head(b.family_m) ## relative abundance based on sum of counts by family!
b.family_m<-subset(b.family_m, SampleID!="Undetermined_S0") ## keep only bacteria -- exclude Chloroplast sequences

b.family_RA_meta<-merge(b.family_m,metadata, by="SampleID")
head(b.family_RA_meta) ## relative abundance based on sum of counts by family!

## genus ...

b.genus_counts <- as.data.frame(dcast(b.dat.m, SampleID~Genus, value.var="Count", fun.aggregate=sum)) ###
head(b.genus_counts) # counts by genus per sample
rownames(b.genus_counts)<-b.genus_counts$SampleID
b.genus_counts$SampleID<-NULL
b.genus_RelAb<-data.frame(decostand(b.genus_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(b.genus_RelAb) # sanity check to make sure the transformation worked!
b.genus_RelAb$SampleID<-rownames(b.genus_RelAb)
head(b.genus_RelAb)
b.genus_m<-melt(b.genus_RelAb)

head(b.genus_m)
colnames(b.genus_m)[which(names(b.genus_m) == "variable")] <- "Genus"
colnames(b.genus_m)[which(names(b.genus_m) == "value")] <- "Count"
head(b.genus_m) ## relative abundance based on sum of counts by genus!
b.genus_m<-subset(b.genus_m, SampleID!="Undetermined_S0") ## keep only bacteria -- exclude Chloroplast sequences

b.genus_RA_meta<-merge(b.genus_m,metadata, by="SampleID")
head(b.genus_RA_meta) ## relative abundance based on sum of counts by genus!

## species ...

b.species_counts <- as.data.frame(dcast(b.dat.m, SampleID~Genus+Species, value.var="Count", fun.aggregate=sum)) ###
head(b.species_counts) # counts by species per sample
rownames(b.species_counts)<-b.species_counts$SampleID
b.species_counts$SampleID<-NULL
b.species_RelAb<-data.frame(decostand(b.species_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(b.species_RelAb) # sanity check to make sure the transformation worked!
b.species_RelAb$SampleID<-rownames(b.species_RelAb)
head(b.species_RelAb)
b.species_m<-melt(b.species_RelAb)

head(b.species_m)
colnames(b.species_m)[which(names(b.species_m) == "variable")] <- "Genus_Species"
colnames(b.species_m)[which(names(b.species_m) == "value")] <- "Count"
b.species_m$Genus_Species<-gsub("_", " ", b.species_m$Genus_Species) ## gsub is global sub (does not just remove first instance of pattern, but multiple)
head(b.species_m) ## relative abundance based on sum of counts by species!
b.species_m<-subset(b.species_m, SampleID!="Undetermined_S0") ## keep only bacteria -- exclude Chloroplast sequences

b.species_RA_meta<-merge(b.species_m,metadata, by="SampleID")
head(b.species_RA_meta) ## relative abundance based on sum of counts by species!

#### ITS2: Relative Abundance by Taxa Level ####
head(f.dat.m)
### * below we use the dcast() function to "cast" the data into a wide format based on given elements (column names), taking sum of "Count"
### * decostand(df, method="total) is the function (with argument total) used to get relative abundance of OTU table

## phylum ....
f.phyla_counts <- as.data.frame(dcast(f.dat.m, SampleID~Phylum, value.var="Count", fun.aggregate=sum)) ###
head(f.phyla_counts) # counts by phyla per sample
rownames(f.phyla_counts)<-f.phyla_counts$SampleID
f.phyla_counts$SampleID<-NULL
f.phyla_RelAb<-data.frame(decostand(f.phyla_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(f.phyla_RelAb) # sanity check to make sure the transformation worked!
f.phyla_RelAb$SampleID<-rownames(f.phyla_RelAb)
head(f.phyla_RelAb)
f.phyla_m<-melt(f.phyla_RelAb)

head(f.phyla_m)
colnames(f.phyla_m)[which(names(f.phyla_m) == "variable")] <- "Phylum"
colnames(f.phyla_m)[which(names(f.phyla_m) == "value")] <- "Count"
head(f.phyla_m) ## relative abundance based on sum of counts by phyla!

f.phyla_RA_meta<-merge(f.phyla_m,metadata, by="SampleID")
head(f.phyla_RA_meta) ## relative abundance based on sum of counts by phyla!

## class...

f.class_counts <- as.data.frame(dcast(f.dat.m, SampleID~Class, value.var="Count", fun.aggregate=sum)) ###
head(f.class_counts) # counts by class per sample
rownames(f.class_counts)<-f.class_counts$SampleID
f.class_counts$SampleID<-NULL
f.class_RelAb<-data.frame(decostand(f.class_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(f.class_RelAb) # sanity check to make sure the transformation worked!
f.class_RelAb$SampleID<-rownames(f.class_RelAb)
head(f.class_RelAb)
f.class_m<-melt(f.class_RelAb)

head(f.class_m)
colnames(f.class_m)[which(names(f.class_m) == "variable")] <- "Class"
colnames(f.class_m)[which(names(f.class_m) == "value")] <- "Count"
head(f.class_m) ## relative abundance based on sum of counts by class!

f.class_RA_meta<-merge(f.class_m,metadata, by="SampleID")
head(f.class_RA_meta) ## relative abundance based on sum of counts by class!

## order ...

f.order_counts <- as.data.frame(dcast(f.dat.m, SampleID~Order, value.var="Count", fun.aggregate=sum)) ###
head(f.order_counts) # counts by order per sample
rownames(f.order_counts)<-f.order_counts$SampleID
f.order_counts$SampleID<-NULL
f.order_RelAb<-data.frame(decostand(f.order_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(f.order_RelAb) # sanity check to make sure the transformation worked!
f.order_RelAb$SampleID<-rownames(f.order_RelAb)
head(f.order_RelAb)
f.order_m<-melt(f.order_RelAb)

head(f.order_m)
colnames(f.order_m)[which(names(f.order_m) == "variable")] <- "Order"
colnames(f.order_m)[which(names(f.order_m) == "value")] <- "Count"
head(f.order_m) ## relative abundance based on sum of counts by order!

f.order_RA_meta<-merge(f.order_m,metadata, by="SampleID")
head(f.order_RA_meta) ## relative abundance based on sum of counts by order!

## family ...

f.family_counts <- as.data.frame(dcast(f.dat.m, SampleID~Family, value.var="Count", fun.aggregate=sum)) ###
head(f.family_counts) # counts by family per sample
rownames(f.family_counts)<-f.family_counts$SampleID
f.family_counts$SampleID<-NULL
f.family_RelAb<-data.frame(decostand(f.family_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(f.family_RelAb) # sanity check to make sure the transformation worked!
f.family_RelAb$SampleID<-rownames(f.family_RelAb)
head(f.family_RelAb)
f.family_m<-melt(f.family_RelAb)

head(f.family_m)
colnames(f.family_m)[which(names(f.family_m) == "variable")] <- "Family"
colnames(f.family_m)[which(names(f.family_m) == "value")] <- "Count"
head(f.family_m) ## relative abundance based on sum of counts by family!

f.family_RA_meta<-merge(f.family_m,metadata, by="SampleID")
head(f.family_RA_meta) ## relative abundance based on sum of counts by family!

## genus ...

f.genus_counts <- as.data.frame(dcast(f.dat.m, SampleID~Genus, value.var="Count", fun.aggregate=sum)) ###
head(f.genus_counts) # counts by genus per sample
rownames(f.genus_counts)<-f.genus_counts$SampleID
f.genus_counts$SampleID<-NULL
f.genus_RelAb<-data.frame(decostand(f.genus_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(f.genus_RelAb) # sanity check to make sure the transformation worked!
f.genus_RelAb$SampleID<-rownames(f.genus_RelAb)
head(f.genus_RelAb)
f.genus_m<-melt(f.genus_RelAb)

head(f.genus_m)
colnames(f.genus_m)[which(names(f.genus_m) == "variable")] <- "Genus"
colnames(f.genus_m)[which(names(f.genus_m) == "value")] <- "Count"
head(f.genus_m) ## relative abundance based on sum of counts by genus!

f.genus_RA_meta<-merge(f.genus_m,metadata, by="SampleID")
head(f.genus_RA_meta) ## relative abundance based on sum of counts by genus!

## species ...

f.species_counts <- as.data.frame(dcast(f.dat.m, SampleID~Genus+Species, value.var="Count", fun.aggregate=sum)) ###
head(f.species_counts) # counts by species per sample
rownames(f.species_counts)<-f.species_counts$SampleID
f.species_counts$SampleID<-NULL
f.species_RelAb<-data.frame(decostand(f.species_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(f.species_RelAb) # sanity check to make sure the transformation worked!
f.species_RelAb$SampleID<-rownames(f.species_RelAb)
head(f.species_RelAb)
f.species_m<-melt(f.species_RelAb)

head(f.species_m)
colnames(f.species_m)[which(names(f.species_m) == "variable")] <- "Genus_Species"
colnames(f.species_m)[which(names(f.species_m) == "value")] <- "Count"
f.species_m$Genus_Species<-gsub("_", " ", f.species_m$Genus_Species) ## gsub is global sub (does not just remove first instance of pattern, but multiple)

head(f.species_m) ## relative abundance based on sum of counts by species!

f.species_RA_meta<-merge(f.species_m,metadata, by="SampleID")
head(f.species_RA_meta) ## relative abundance based on sum of counts by species!

#### Checking out some color pallettes before we visualize ####
wes1<-wes_palette("Chevalier1")
wes2<-wes_palette("Moonrise3")
wes3<-wes_palette("IsleofDogs1")
wes4<-wes_palette("GrandBudapest1")
wes5<-wes_palette("GrandBudapest2")

SM_pal <- park_palette("SmokyMountains") # create a palette and specify # of colors youw ant
Arc_pal <- park_palette("Arches") # create a palette and specify # of colors youw ant
CL_pal <- park_palette("CraterLake") # create a palette and specify # of colors youw ant
Sag_pal <- park_palette("Saguaro") # create a palette and specify # of colors youw ant
Aca_pal <- park_palette("Acadia") # create a palette and specify # of colors youw ant
DV_pal <- park_palette("DeathValley") # create a palette and specify # of colors youw ant
DV_pal2 <- park_palette("DeathValley", 1) # create a palette and specify # of colors youw ant
CI_pal <- park_palette("ChannelIslands") # create a palette and specify # of colors youw ant
Bad_pal <- park_palette("Badlands") # create a palette and specify # of colors youw ant
MR_pal2 <- park_palette("MtRainier", 1) # create a palette and specify # of colors youw ant
MR_pal <- park_palette("MtRainier") # create a palette and specify # of colors youw ant

HI_pal <- park_palette("Hawaii") # create a palette and specify # of colors youw ant

bright_pal<-get_palette(palette = paste0("#", c("FFBE0B","D03325","8338EC","3A86FF")), k = 4)

ewf<-get_palette(palette = paste0("#", c("E9AE3D","694344","97A7C0","132138")), k = 4)

colorset = c(air="#E9AE3D",soil="#694344",water="#97A7C0",control="#132138")
# colorset1 = c(Actinobacteriota="#FCAF58",Asgardarchaeota="#3A86FF",Bacteroidota="#58A4B0",Chloroflexi="#43AA8B",
#               Desulfobacterota="#7400B8",Firmicutes="#90BE6D",Patescibacteria="#990066",
#               Planctomycetota="#70A8FF",Proteobacteria="#9A031E",Spirochaetota="#CC3399",Unknown="#74DABC", Verrucomicrobiota="#EEB625")

#fungal_pal<-get_palette(palette = paste0("#", c("222529","c72421","da7049","eedc2b","c5debf","b6f2fa","513a56","fafccf")), k = 8)
#fungal_pal2<-get_palette(palette = paste0("#", c("07060f","dd1f19","eda964","bba20d","44be45","74dabc","5e2479","5f3523")), k = 8)

#scale_fill_manual(values = wes_palette("IsleofDogs1"))


#### Visualize 16S Relative Abundances ####

# phyla
b.phyla_RA_2<-subset(b.phyla_m, c(Count)>(1/100)) ## DROP BACTERIAL PHYLA that are less than 1% abundant!!!!!!1
b.phyla_RA_3<-subset(b.phyla_m, c(Count)>(10/100)) ## DROP BACTERIAL PHYLA that are less than 10% abundant!!!!!!1

b.phyla_f1<-ggplot(b.phyla_m, aes(x=SampleID, y=Count, fill=Phylum)) + geom_bar(stat="identity", colour="black")+theme_bw()+theme_classic()+scale_x_discrete()+labs(title="Microbial Phyla",subtitle="16S (Bacteria & Archaea) Ampicon Sequencing",x="Microbial Phyla", y="Relative Abundance")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5,
        legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11),
        plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),plot.caption=element_text(hjust=0.5))
ggsave(b.phyla_f1,filename = "figures/16S_phyla_RelativeAbundance_wArchaea_6.16.21.pdf", width=15, height=10, dpi=600)

b.phyla_f2<-ggplot(b.phyla_RA_2, aes(x=SampleID, y=Count, fill=Phylum)) + geom_bar(stat="identity",colour="black")+theme_bw()+theme_classic()+scale_x_discrete()+labs(title="Microbial Phyla",subtitle="16S (Bacteria & Archaea) Ampicon Sequencing",x="Microbial Phyla", y="Relative Abundance",caption="Includes Taxa > 1% Relative Abundance")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5,
        legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11),
        plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),plot.caption=element_text(hjust=0.5))
ggsave(b.phyla_f2,filename = "figures/16S_phyla_1percent_RelativeAbundance_wArchaea_6.16.21.pdf", width=15, height=10, dpi=600)

b.phyla_f3<-ggplot(b.phyla_RA_3, aes(x=SampleID, y=Count, fill=Phylum)) + geom_bar(stat="identity",colour="black")+theme_bw()+theme_classic()+scale_x_discrete()+labs(title="Microbial Phyla",subtitle="16S (Bacteria & Archaea) Ampicon Sequencing",x="Microbial Phyla", y="Relative Abundance",caption="Includes Taxa > 10% Relative Abundance")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5,
        legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11),
        plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),plot.caption=element_text(hjust=0.5))
ggsave(b.phyla_f2,filename = "figures/16S_phyla_10percent_RelativeAbundance_wArchaea_6.16.21.pdf", width=15, height=10, dpi=600)

# class
b.class_RA_melt2<-subset(b.class_m, c(Count)>(1/100)) ## DROP BACTERIAL CLASS that are less than 1% abundant!!!!!!1\
b.class_RA_melt3<-subset(b.class_m, c(Count)>(10/100)) ## DROP BACTERIAL CLASS that are less than 10% abundant!!!!!!1

b.class_f1<-ggplot(b.class_m, aes(x=SampleID, y=Count, fill=Class)) + geom_bar(stat="identity",colour="black")+theme_bw()+theme_classic()+scale_x_discrete()+labs(title="Microbial Classes",subtitle="16S (Bacteria & Archaea) Ampicon Sequencing",x="Microbial Class", y="Relative Abundance")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5,
        legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11),
        plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),plot.caption=element_text(hjust=0.5))
ggsave(b.class_f1,filename = "figures/16S_class_RelativeAbundance_wArchaea_6.16.21.pdf", width=20, height=10, dpi=600)

b.class_f2<-ggplot(b.class_RA_melt2, aes(x=SampleID, y=Count, fill=Class)) + geom_bar(stat="identity",colour="black")+theme_bw()+theme_classic()+scale_x_discrete()+labs(title="Microbial Classes",subtitle="16S (Bacteria & Archaea) Ampicon Sequencing",x="Microbial Class", y="Relative Abundance",caption="Includes Taxa > 1% Relative Abundance")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5,
        legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11),
        plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),plot.caption=element_text(hjust=0.5))
ggsave(b.class_f2,filename = "figures/16S_class_1percent_RelativeAbundance_wArchaea_6.16.21.pdf", width=20, height=10, dpi=600)

b.class_f3<-ggplot(b.class_RA_melt3, aes(x=SampleID, y=Count, fill=Class)) + geom_bar(stat="identity",colour="black")+theme_bw()+theme_classic()+scale_x_discrete()+labs(title="Microbial Classes",subtitle="16S (Bacteria & Archaea) Ampicon Sequencing",x="Microbial Class", y="Relative Abundance",caption="Includes Taxa > 10% Relative Abundance")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5,
        legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11),
        plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),plot.caption=element_text(hjust=0.5))
ggsave(b.class_f2,filename = "figures/16S_class_10percent_RelativeAbundance_wArchaea_6.16.21.pdf", width=20, height=10, dpi=600)

#+scale_y_continuous(limits=c(0,0.75))

# order
b.order_f1<-ggplot(b.order_m, aes(x=SampleID, y=Count, fill=Order)) + geom_bar(stat="identity",colour="black")+theme_bw()+theme_classic()+scale_x_discrete()+labs(title="Microbial Orders",subtitle="16S (Bacteria & Archaea) Ampicon Sequencing",x="Microbial Order", y="Relative Abundance")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5,
        legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11),
        plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),plot.caption=element_text(hjust=0.5))
ggsave(b.order_f1,filename = "figures/16S_order_RelativeAbundance_wArchaea_6.16.21.pdf", width=25, height=10, dpi=600)

b.order_RA_melt2<-subset(b.order_m, c(Count)>(1/100)) ## DROP BACTERIAL ORDER that are less than 1% abundant!!!!!!1\
b.order_RA_melt3<-subset(b.order_m, c(Count)>(10/100)) ## DROP BACTERIAL ORDER that are less than 10% abundant!!!!!!1

b.order_f2<-ggplot(b.order_RA_melt2, aes(x=SampleID, y=Count, fill=Order)) + geom_bar(stat="identity",colour="black")+theme_bw()+theme_classic()+scale_x_discrete()+labs(title="Microbial Orders",subtitle="16S (Bacteria & Archaea) Ampicon Sequencing",x="Microbial Order", y="Relative Abundance",caption="Includes Taxa > 1% Relative Abundance")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5,
        legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11),
        plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),plot.caption=element_text(hjust=0.5))
ggsave(b.order_f2,filename = "figures/16S_order_1percent_RelativeAbundance_wArchaea_6.16.21.pdf", width=20, height=10, dpi=600)

b.order_f3<-ggplot(b.order_RA_melt3, aes(x=SampleID, y=Count, fill=Order)) + geom_bar(stat="identity",colour="black")+theme_bw()+theme_classic()+scale_x_discrete()+labs(title="Microbial Orders",subtitle="16S (Bacteria & Archaea) Ampicon Sequencing",x="Microbial Order", y="Relative Abundance",caption="Includes Taxa > 10% Relative Abundance")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5,
        legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11),
        plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),plot.caption=element_text(hjust=0.5))
ggsave(b.order_f2,filename = "figures/16S_order_10percent_RelativeAbundance_wArchaea_6.16.21.pdf", width=20, height=10, dpi=600)

# family
b.family_f1<-ggplot(b.family_m, aes(x=SampleID, y=Count, fill=Family)) + geom_bar(stat="identity",colour="black")+theme_bw()+theme_classic()+scale_x_discrete()+labs(title="Microbial Families",subtitle="16S (Bacteria & Archaea) Ampicon Sequencing",x="Microbial Family", y="Relative Abundance")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5,
        legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11),
        plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),plot.caption=element_text(hjust=0.5))
ggsave(b.family_f1,filename = "figures/16S_family_RelativeAbundance_wArchaea_6.16.21.pdf", width=25, height=10, dpi=600)

b.family_RA_melt2<-subset(b.family_m, c(Count)>(1/100)) ## DROP BACTERIAL FAMILY that are less than 1% abundant!!!!!!1\
b.family_RA_melt3<-subset(b.family_m, c(Count)>(10/100)) ## DROP BACTERIAL FAMILY that are less than 10% abundant!!!!!!1

b.family_f2<-ggplot(b.family_RA_melt2, aes(x=SampleID, y=Count, fill=Family)) + geom_bar(stat="identity",colour="black")+theme_bw()+theme_classic()+scale_x_discrete()+labs(title="Microbial Families",subtitle="16S (Bacteria & Archaea) Ampicon Sequencing",x="Microbial Family", y="Relative Abundance",caption="Includes Taxa > 1% Relative Abundance")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5,
        legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11),
        plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),plot.caption=element_text(hjust=0.5))
ggsave(b.family_f2,filename = "figures/16S_family_1percent_RelativeAbundance_wArchaea_6.16.21.pdf", width=20, height=10, dpi=600)

b.family_f3<-ggplot(b.family_RA_melt3, aes(x=SampleID, y=Count, fill=Family)) + geom_bar(stat="identity",colour="black")+theme_bw()+theme_classic()+scale_x_discrete()+labs(title="Microbial Families",subtitle="16S (Bacteria & Archaea) Ampicon Sequencing",x="Microbial Family", y="Relative Abundance",caption="Includes Taxa > 10% Relative Abundance")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5,
        legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11),
        plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),plot.caption=element_text(hjust=0.5))
ggsave(b.family_f2,filename = "figures/16S_family_10percent_RelativeAbundance_wArchaea_6.16.21.pdf", width=20, height=10, dpi=600)

# genus
b.genus_f1<-ggplot(b.genus_m, aes(x=SampleID, y=Count, fill=Genus)) + geom_bar(stat="identity",colour="black")+theme_bw()+theme_classic()+scale_x_discrete()+labs(title="Microbial Genera",subtitle="16S (Bacteria & Archaea) Ampicon Sequencing",x="Microbial Genus", y="Relative Abundance")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5,
        legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11),
        plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),plot.caption=element_text(hjust=0.5))
ggsave(b.genus_f1,filename = "figures/16S_genus_RelativeAbundance_wArchaea_6.16.21.pdf", width=35, height=10, dpi=600)

b.genus_RA_melt2<-subset(b.genus_m, c(Count)>(1/100)) ## DROP BACTERIAL GENUS that are less than 1% abundant!!!!!!1\
b.genus_RA_melt3<-subset(b.genus_m, c(Count)>(10/100)) ## DROP BACTERIAL GENUS that are less than 10% abundant!!!!!!1

b.genus_f2<-ggplot(b.genus_RA_melt2, aes(x=SampleID, y=Count, fill=Genus)) + geom_bar(stat="identity",colour="black")+theme_bw()+theme_classic()+scale_x_discrete()+labs(title="Microbial Genera",subtitle="16S (Bacteria & Archaea) Ampicon Sequencing",x="Microbial Genus", y="Relative Abundance",caption="Includes Taxa > 1% Relative Abundance")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5,
        legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11),
        plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),plot.caption=element_text(hjust=0.5))
ggsave(b.genus_f2,filename = "figures/16S_genus_1percent_RelativeAbundance_wArchaea_6.16.21.pdf", width=30, height=10, dpi=600)

b.genus_f3<-ggplot(b.genus_RA_melt3, aes(x=SampleID, y=Count, fill=Genus)) + geom_bar(stat="identity",colour="black")+theme_bw()+theme_classic()+scale_x_discrete()+labs(title="Microbial Genera",subtitle="16S (Bacteria & Archaea) Ampicon Sequencing",x="Microbial Genus", y="Relative Abundance",caption="Includes Taxa > 10% Relative Abundance")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5,
        legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11),
        plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),plot.caption=element_text(hjust=0.5))
ggsave(b.genus_f3,filename = "figures/16S_genus_10percent_RelativeAbundance_wArchaea_6.16.21.pdf", width=20, height=10, dpi=600)

# species
# b.species_f1<-ggplot(b.species_m, aes(x=SampleID, y=Count, fill=Genus_Species)) + geom_bar(stat="identity",colour="black")+theme_bw()+theme_classic()+scale_x_discrete()+labs(title="Microbial Species",subtitle="16S (Bacteria & Archaea) Ampicon Sequencing",x="Microbial Species", y="Relative Abundance")+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5,
#         legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11),
#         plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),plot.caption=element_text(hjust=0.5))
# ggsave(b.species_f1,filename = "figures/16S_species_RelativeAbundance_wArchaea_6.16.21.pdf", width=45, height=10, dpi=600)

b.species_RA_melt2<-subset(b.species_m, c(Count)>(1/100)) ## DROP BACTERIAL SPECIES that are less than 1% abundant!!!!!!1\
b.species_RA_melt3<-subset(b.species_m, c(Count)>(10/100)) ## DROP BACTERIAL SPECIES that are less than 10% abundant!!!!!!1

b.species_f2<-ggplot(b.species_RA_melt2, aes(x=SampleID, y=Count, fill=Genus_Species)) + geom_bar(stat="identity",colour="black")+theme_bw()+theme_classic()+scale_x_discrete()+labs(title="Microbial Species",subtitle="16S (Bacteria & Archaea) Ampicon Sequencing",x="Microbial Species", y="Relative Abundance",caption="Includes Taxa > 1% Relative Abundance")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5,
        legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11),
        plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),plot.caption=element_text(hjust=0.5))
ggsave(b.species_f2,filename = "figures/16S_species_1percent_RelativeAbundance_wArchaea_6.16.21.pdf", width=45, height=10, dpi=600)

b.species_f3<-ggplot(b.species_RA_melt3, aes(x=SampleID, y=Count, fill=Genus_Species)) + geom_bar(stat="identity",colour="black")+theme_bw()+theme_classic()+scale_x_discrete()+labs(title="Microbial Species",subtitle="16S (Bacteria & Archaea) Ampicon Sequencing",x="Microbial Species", y="Relative Abundance",caption="Includes Taxa > 10% Relative Abundance")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5,
        legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11),
        plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),plot.caption=element_text(hjust=0.5))
ggsave(b.species_f3,filename = "figures/16S_species_10percent_RelativeAbundance_wArchaea_6.16.21.pdf", width=20, height=10, dpi=600)

### ^^^^ plotting with relative abundance data (decostand "total" transformation where all counts are divided by total hits across taxa per sample)

#### Visualize ITS2 Relative Abundances ####

# phyla
f.phyla_RA_2<-subset(f.phyla_m, c(Count)>(1/100)) ## DROP FUNGAL PHYLA that are less than 1% abundant!!!!!!1
f.phyla_RA_3<-subset(f.phyla_m, c(Count)>(10/100)) ## DROP FUNGAL PHYLA that are less than 10% abundant!!!!!!1

f.phyla_f1<-ggplot(f.phyla_m, aes(x=SampleID, y=Count, fill=Phylum)) + geom_bar(stat="identity", colour="black")+theme_bw()+theme_classic()+scale_x_discrete()+labs(title="Fungal Phyla",subtitle="ITS2 Ampicon Sequencing",x="Fungal Phyla", y="Relative Abundance")+
 theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5,
    legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11),
    plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),plot.caption=element_text(hjust=0.5))
ggsave(f.phyla_f1,filename = "figures/ITS2_phyla_RelativeAbundance_6.16.21.pdf", width=15, height=10, dpi=600)

f.phyla_f2<-ggplot(f.phyla_RA_2, aes(x=SampleID, y=Count, fill=Phylum)) + geom_bar(stat="identity",colour="black")+theme_bw()+theme_classic()+scale_x_discrete()+labs(title="Fungal Phyla",subtitle="ITS2 Ampicon Sequencing",x="Fungal Phyla", y="Relative Abundance",caption="Includes Taxa > 1% Relative Abundance")+
 theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5,
    legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11),
    plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),plot.caption=element_text(hjust=0.5))
ggsave(f.phyla_f2,filename = "figures/ITS2_phyla_1percent_RelativeAbundance_6.16.21.pdf", width=15, height=10, dpi=600)

f.phyla_f3<-ggplot(f.phyla_RA_3, aes(x=SampleID, y=Count, fill=Phylum)) + geom_bar(stat="identity",colour="black")+theme_bw()+theme_classic()+scale_x_discrete()+labs(title="Fungal Phyla",subtitle="ITS2 Ampicon Sequencing",x="Fungal Phyla", y="Relative Abundance",caption="Includes Taxa > 10% Relative Abundance")+
 theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5,
    legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11),
    plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),plot.caption=element_text(hjust=0.5))
ggsave(f.phyla_f2,filename = "figures/ITS2_phyla_10percent_RelativeAbundance_6.16.21.pdf", width=15, height=10, dpi=600)

# class
f.class_RA_melt2<-subset(f.class_m, c(Count)>(1/100)) ## DROP FUNGAL CLASS that are less than 1% abundant!!!!!!1\
f.class_RA_melt3<-subset(f.class_m, c(Count)>(10/100)) ## DROP FUNGAL CLASS that are less than 10% abundant!!!!!!1

f.class_f1<-ggplot(f.class_m, aes(x=SampleID, y=Count, fill=Class)) + geom_bar(stat="identity",colour="black")+theme_bw()+theme_classic()+scale_x_discrete()+labs(title="Fungal Classes",subtitle="ITS2 Ampicon Sequencing",x="Fungal Class", y="Relative Abundance")+
 theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5,
    legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11),
    plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),plot.caption=element_text(hjust=0.5))
ggsave(f.class_f1,filename = "figures/ITS2_class_RelativeAbundance_6.16.21.pdf", width=20, height=10, dpi=600)

f.class_f2<-ggplot(f.class_RA_melt2, aes(x=SampleID, y=Count, fill=Class)) + geom_bar(stat="identity",colour="black")+theme_bw()+theme_classic()+scale_x_discrete()+labs(title="Fungal Classes",subtitle="ITS2 Ampicon Sequencing",x="Fungal Class", y="Relative Abundance",caption="Includes Taxa > 1% Relative Abundance")+
 theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5,
    legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11),
    plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),plot.caption=element_text(hjust=0.5))
ggsave(f.class_f2,filename = "figures/ITS2_class_1percent_RelativeAbundance_6.16.21.pdf", width=20, height=10, dpi=600)

f.class_f3<-ggplot(f.class_RA_melt3, aes(x=SampleID, y=Count, fill=Class)) + geom_bar(stat="identity",colour="black")+theme_bw()+theme_classic()+scale_x_discrete()+labs(title="Fungal Classes",subtitle="ITS2 Ampicon Sequencing",x="Fungal Class", y="Relative Abundance",caption="Includes Taxa > 10% Relative Abundance")+
 theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5,
    legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11),
    plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),plot.caption=element_text(hjust=0.5))
ggsave(f.class_f2,filename = "figures/ITS2_class_10percent_RelativeAbundance_6.16.21.pdf", width=20, height=10, dpi=600)

#+scale_y_continuous(limits=c(0,0.75))

# order
f.order_f1<-ggplot(f.order_m, aes(x=SampleID, y=Count, fill=Order)) + geom_bar(stat="identity",colour="black")+theme_bw()+theme_classic()+scale_x_discrete()+labs(title="Fungal Orders",subtitle="ITS2 Ampicon Sequencing",x="Fungal Order", y="Relative Abundance")+
 theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5,
    legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11),
    plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),plot.caption=element_text(hjust=0.5))
ggsave(f.order_f1,filename = "figures/ITS2_order_RelativeAbundance_6.16.21.pdf", width=15, height=10, dpi=600)

f.order_RA_melt2<-subset(f.order_m, c(Count)>(1/100)) ## DROP FUNGAL ORDER that are less than 1% abundant!!!!!!1\
f.order_RA_melt3<-subset(f.order_m, c(Count)>(10/100)) ## DROP FUNGAL ORDER that are less than 10% abundant!!!!!!1

f.order_f2<-ggplot(f.order_RA_melt2, aes(x=SampleID, y=Count, fill=Order)) + geom_bar(stat="identity",colour="black")+theme_bw()+theme_classic()+scale_x_discrete()+labs(title="Fungal Orders",subtitle="ITS2 Ampicon Sequencing",x="Fungal Order", y="Relative Abundance",caption="Includes Taxa > 1% Relative Abundance")+
 theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5,
    legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11),
    plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),plot.caption=element_text(hjust=0.5))
ggsave(f.order_f2,filename = "figures/ITS2_order_1percent_RelativeAbundance_6.16.21.pdf", width=15, height=10, dpi=600)

f.order_f3<-ggplot(f.order_RA_melt3, aes(x=SampleID, y=Count, fill=Order)) + geom_bar(stat="identity",colour="black")+theme_bw()+theme_classic()+scale_x_discrete()+labs(title="Fungal Orders",subtitle="ITS2 Ampicon Sequencing",x="Fungal Order", y="Relative Abundance",caption="Includes Taxa > 10% Relative Abundance")+
 theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5,
    legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11),
    plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),plot.caption=element_text(hjust=0.5))
ggsave(f.order_f2,filename = "figures/ITS2_order_10percent_RelativeAbundance_6.16.21.pdf", width=15, height=10, dpi=600)

# family
f.family_f1<-ggplot(f.family_m, aes(x=SampleID, y=Count, fill=Family)) + geom_bar(stat="identity",colour="black")+theme_bw()+theme_classic()+scale_x_discrete()+labs(title="Fungal Families",subtitle="ITS2 Ampicon Sequencing",x="Fungal Family", y="Relative Abundance")+
 theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5,
    legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11),
    plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),plot.caption=element_text(hjust=0.5))
ggsave(f.family_f1,filename = "figures/ITS2_family_RelativeAbundance_6.16.21.pdf", width=15, height=10, dpi=600)

f.family_RA_melt2<-subset(f.family_m, c(Count)>(1/100)) ## DROP FUNGAL FAMILY that are less than 1% abundant!!!!!!1\
f.family_RA_melt3<-subset(f.family_m, c(Count)>(10/100)) ## DROP FUNGAL FAMILY that are less than 10% abundant!!!!!!1

f.family_f2<-ggplot(f.family_RA_melt2, aes(x=SampleID, y=Count, fill=Family)) + geom_bar(stat="identity",colour="black")+theme_bw()+theme_classic()+scale_x_discrete()+labs(title="Fungal Families",subtitle="ITS2 Ampicon Sequencing",x="Fungal Family", y="Relative Abundance",caption="Includes Taxa > 1% Relative Abundance")+
 theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5,
    legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11),
    plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),plot.caption=element_text(hjust=0.5))
ggsave(f.family_f2,filename = "figures/ITS2_family_1percent_RelativeAbundance_6.16.21.pdf", width=15, height=10, dpi=600)

f.family_f3<-ggplot(f.family_RA_melt3, aes(x=SampleID, y=Count, fill=Family)) + geom_bar(stat="identity",colour="black")+theme_bw()+theme_classic()+scale_x_discrete()+labs(title="Fungal Families",subtitle="ITS2 Ampicon Sequencing",x="Fungal Family", y="Relative Abundance",caption="Includes Taxa > 10% Relative Abundance")+
 theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5,
    legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11),
    plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),plot.caption=element_text(hjust=0.5))
ggsave(f.family_f2,filename = "figures/ITS2_family_10percent_RelativeAbundance_6.16.21.pdf", width=15, height=10, dpi=600)

# genus
f.genus_f1<-ggplot(f.genus_m, aes(x=SampleID, y=Count, fill=Genus)) + geom_bar(stat="identity",colour="black")+theme_bw()+theme_classic()+scale_x_discrete()+labs(title="Fungal Genera",subtitle="ITS2 Ampicon Sequencing",x="Fungal Genus", y="Relative Abundance")+
 theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5,
    legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11),
    plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),plot.caption=element_text(hjust=0.5))
ggsave(f.genus_f1,filename = "figures/ITS2_genus_RelativeAbundance_6.16.21.pdf", width=15, height=10, dpi=600)

f.genus_RA_melt2<-subset(f.genus_m, c(Count)>(1/100)) ## DROP FUNGAL GENUS that are less than 1% abundant!!!!!!1\
f.genus_RA_melt3<-subset(f.genus_m, c(Count)>(10/100)) ## DROP FUNGAL GENUS that are less than 10% abundant!!!!!!1

f.genus_f2<-ggplot(f.genus_RA_melt2, aes(x=SampleID, y=Count, fill=Genus)) + geom_bar(stat="identity",colour="black")+theme_bw()+theme_classic()+scale_x_discrete()+labs(title="Fungal Genera",subtitle="ITS2 Ampicon Sequencing",x="Fungal Genus", y="Relative Abundance",caption="Includes Taxa > 1% Relative Abundance")+
 theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5,
    legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11),
    plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),plot.caption=element_text(hjust=0.5))
ggsave(f.genus_f2,filename = "figures/ITS2_genus_1percent_RelativeAbundance_6.16.21.pdf", width=15, height=10, dpi=600)

f.genus_f3<-ggplot(f.genus_RA_melt3, aes(x=SampleID, y=Count, fill=Genus)) + geom_bar(stat="identity",colour="black")+theme_bw()+theme_classic()+scale_x_discrete()+labs(title="Fungal Genera",subtitle="ITS2 Ampicon Sequencing",x="Fungal Genus", y="Relative Abundance",caption="Includes Taxa > 10% Relative Abundance")+
 theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5,
    legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11),
    plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),plot.caption=element_text(hjust=0.5))
ggsave(f.genus_f3,filename = "figures/ITS2_genus_10percent_RelativeAbundance_6.16.21.pdf", width=15, height=10, dpi=600)

# species
# f.species_f1<-ggplot(f.species_m, aes(x=SampleID, y=Count, fill=Genus_Species)) + geom_bar(stat="identity",colour="black")+theme_bw()+theme_classic()+scale_x_discrete()+labs(title="Fungal Species",subtitle="ITS2 Ampicon Sequencing",x="Fungal Species", y="Relative Abundance")+
#  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5,
#     legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11),
#     plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),plot.caption=element_text(hjust=0.5))
# ggsave(f.species_f1,filename = "figures/ITS2_species_RelativeAbundance_6.16.21.pdf", width=45, height=10, dpi=600)

f.species_RA_melt2<-subset(f.species_m, c(Count)>(1/100)) ## DROP FUNGAL SPECIES that are less than 1% abundant!!!!!!1\
f.species_RA_melt3<-subset(f.species_m, c(Count)>(10/100)) ## DROP FUNGAL SPECIES that are less than 10% abundant!!!!!!1

f.species_f2<-ggplot(f.species_RA_melt2, aes(x=SampleID, y=Count, fill=Genus_Species)) + geom_bar(stat="identity",colour="black")+theme_bw()+theme_classic()+scale_x_discrete()+labs(title="Fungal Species",subtitle="ITS2 Ampicon Sequencing",x="Fungal Species", y="Relative Abundance",caption="Includes Taxa > 1% Relative Abundance")+
 theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5,
    legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11),
    plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),plot.caption=element_text(hjust=0.5))
ggsave(f.species_f2,filename = "figures/ITS2_species_1percent_RelativeAbundance_6.16.21.pdf", width=45, height=10, dpi=600)

f.species_f3<-ggplot(f.species_RA_melt3, aes(x=SampleID, y=Count, fill=Genus_Species)) + geom_bar(stat="identity",colour="black")+theme_bw()+theme_classic()+scale_x_discrete()+labs(title="Fungal Species",subtitle="ITS2 Ampicon Sequencing",x="Fungal Species", y="Relative Abundance",caption="Includes Taxa > 10% Relative Abundance")+
 theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5,
    legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11),
    plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),plot.caption=element_text(hjust=0.5))
ggsave(f.species_f3,filename = "figures/ITS2_species_10percent_RelativeAbundance_6.16.21.pdf", width=20, height=10, dpi=600)

### ^^^^ plotting with relative abundance data (decostand "total" transformation where all counts are divided by total hits across taxa per sample)


#### 16S Relative Abundance by Exposure Group ####
head(b.dat.m)
bac.asv_meta<-merge(b.dat.m, metadata, by="SampleID")

## phylum ....

head(bac.asv_meta)

# by phylum + exposure
bac_phy_exp<- as.data.frame(dcast(bac.asv_meta,Exposed~Phylum, value.var="Count", fun.aggregate=sum)) ###
head(bac_phy_exp)
rownames(bac_phy_exp)<-bac_phy_exp$Exposed
head(bac_phy_exp)

b.phy.RA_exp<-data.frame(decostand(bac_phy_exp[,-c(1)], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(b.phy.RA_exp)
head(b.phy.RA_exp)

b.phy.RA_exp$Exposed<-rownames(b.phy.RA_exp)
b.phy.RA_exp.m<-melt(b.phy.RA_exp)
head(b.phy.RA_exp.m)
colnames(b.phy.RA_exp.m)[which(names(b.phy.RA_exp.m) == "variable")] <- "Phylum"
colnames(b.phy.RA_exp.m)[which(names(b.phy.RA_exp.m) == "value")] <- "Count"
head(b.phy.RA_exp.m) ## relative abundance based on sum of counts by phyla in exposure group!

# by class + exposure
bac_cls_exp <- as.data.frame(dcast(bac.asv_meta,Exposed~Class, value.var="Count", fun.aggregate=sum)) ###
head(bac_cls_exp) # counts by class + elevation
rownames(bac_cls_exp)<-bac_cls_exp$Exposed

b.cls.RA_exp<-data.frame(decostand(bac_cls_exp[,-c(1)], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(b.cls.RA_exp)
head(b.cls.RA_exp)

b.cls.RA_exp$Exposed<-rownames(b.cls.RA_exp)
b.cls.RA_exp.m<-melt(b.cls.RA_exp)
head(b.cls.RA_exp.m)
colnames(b.cls.RA_exp.m)[which(names(b.cls.RA_exp.m) == "variable")] <- "Class"
colnames(b.cls.RA_exp.m)[which(names(b.cls.RA_exp.m) == "value")] <- "Count"
head(b.cls.RA_exp.m) ## relative abundance based on sum of counts by class in exposure group!

# by order + exposure
bac_ord_exp <- as.data.frame(dcast(bac.asv_meta,Exposed~Order, value.var="Count", fun.aggregate=sum)) ###
head(bac_ord_exp) # counts by Order + elevation
rownames(bac_ord_exp)<-bac_ord_exp$Exposed

b.ord.RA_exp<-data.frame(decostand(bac_ord_exp[,-c(1)], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(b.ord.RA_exp)
head(b.ord.RA_exp)

b.ord.RA_exp$Exposed<-rownames(b.ord.RA_exp)
b.ord.RA_exp.m<-melt(b.ord.RA_exp)
head(b.ord.RA_exp.m)
colnames(b.ord.RA_exp.m)[which(names(b.ord.RA_exp.m) == "variable")] <- "Order"
colnames(b.ord.RA_exp.m)[which(names(b.ord.RA_exp.m) == "value")] <- "Count"
head(b.ord.RA_exp.m) ## relative abundance based on sum of counts by order by exposed group!

# by Family + exposure
bac_fam_exp <- as.data.frame(dcast(bac.asv_meta,Exposed~Family, value.var="Count", fun.aggregate=sum)) ###
head(bac_fam_exp) # counts by Family + elevation
rownames(bac_fam_exp)<-bac_fam_exp$Exposed

b.fam.RA_exp<-data.frame(decostand(bac_fam_exp[,-c(1)], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(b.fam.RA_exp)
head(b.fam.RA_exp)

b.fam.RA_exp$Exposed<-rownames(b.fam.RA_exp)
b.fam.RA_exp.m<-melt(b.fam.RA_exp)
head(b.fam.RA_exp.m)
colnames(b.fam.RA_exp.m)[which(names(b.fam.RA_exp.m) == "variable")] <- "Family"
colnames(b.fam.RA_exp.m)[which(names(b.fam.RA_exp.m) == "value")] <- "Count"
head(b.fam.RA_exp.m) ## relative abundance based on sum of counts by family by exposed group!

# by Genus + exposure
bac_gen_exp <- as.data.frame(dcast(bac.asv_meta,Exposed~Genus, value.var="Count", fun.aggregate=sum)) ###
head(bac_gen_exp) # counts by genus + elevation
rownames(bac_gen_exp)<-bac_gen_exp$Exposed

b.gen.RA_exp<-data.frame(decostand(bac_gen_exp[,-c(1)], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(b.gen.RA_exp)
head(b.gen.RA_exp)

b.gen.RA_exp$Exposed<-rownames(b.gen.RA_exp)
b.gen.RA_exp.m<-melt(b.gen.RA_exp)
head(b.gen.RA_exp.m)
colnames(b.gen.RA_exp.m)[which(names(b.gen.RA_exp.m) == "variable")] <- "Genus"
colnames(b.gen.RA_exp.m)[which(names(b.gen.RA_exp.m) == "value")] <- "Count"
head(b.gen.RA_exp.m) ## relative abundance based on sum of counts by genus by exposed group!

# by Genus/Species + exposure
bac_spec_exp <- as.data.frame(dcast(bac.asv_meta,Exposed~Genus+Species, value.var="Count", fun.aggregate=sum)) ###
head(bac_spec_exp) # counts by genus_species + elevation
rownames(bac_spec_exp)<-bac_spec_exp$Exposed

b.spec.RA_exp<-data.frame(decostand(bac_spec_exp[,-c(1)], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(b.spec.RA_exp)
head(b.spec.RA_exp)

b.spec.RA_exp$Exposed<-rownames(b.spec.RA_exp)
b.spec.RA_exp.m<-melt(b.spec.RA_exp)
head(b.spec.RA_exp.m)
colnames(b.spec.RA_exp.m)[which(names(b.spec.RA_exp.m) == "variable")] <- "Genus_species"
colnames(b.spec.RA_exp.m)[which(names(b.spec.RA_exp.m) == "value")] <- "Count"
b.spec.RA_exp.m$Genus_species<-gsub("_", " ", b.spec.RA_exp.m$Genus_species) ## gsub is global sub (does not just remove first instance of pattern, but multiple)

head(b.spec.RA_exp.m) ## relative abundance based on sum of counts by specus by exposed group!

#### Visualize 16S Relative Abundance by Exposure Group ####

# phylum
b.phy.RA_exp.m1<-subset(b.phy.RA_exp.m, c(Count)>(0.01)) ## DROP BACTERIAL PHYLA that are less than 1% abundant!!!!!!1
b.phy.RA_exp.m5<-subset(b.phy.RA_exp.m, c(Count)>(0.05)) ## DROP BACTERIAL PHYLA that are less than 5% abundant!!!!!!1
b.phy.RA_exp.m10<-subset(b.phy.RA_exp.m, c(Count)>(0.10)) ## DROP BACTERIAL PHYLA that are less than 10% abundant!!!!!!1

b.exp1<-ggplot(b.phy.RA_exp.m, aes(x=Exposed, y=Count, fill=Phylum))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Microbial Phylum Relative Abundance by Exposure Group", x="Exposure Material", y="Relative Abundance", fill="Phylum")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(b.exp1,filename = "figures/microb_phyla_RA_Expo_6.18.21.pdf", width=8, height=6, dpi=600)

b.exp1a<-ggplot(b.phy.RA_exp.m1, aes(x=Exposed, y=Count, fill=Phylum))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Microbial Phylum Relative Abundance by Exposure Group", subtitle="Taxa > 1% Relative Abundance", x="Exposure Material", y="Relative Abundance", fill="Phylum")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(b.exp1a,filename = "figures/microb_phyla_1percent_RA_Expo_6.18.21.pdf", width=8, height=6, dpi=600)

b.exp1b<-ggplot(b.phy.RA_exp.m5, aes(x=Exposed, y=Count, fill=Phylum))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Microbial Phylum Relative Abundance by Exposure Group", subtitle="Taxa > 5% Relative Abundance", x="Exposure Material", y="Relative Abundance", fill="Phylum")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(b.exp1b,filename = "figures/microb_phyla_5percent_RA_Expo_6.18.21.pdf", width=8, height=6, dpi=600)

b.exp1c<-ggplot(b.phy.RA_exp.m10, aes(x=Exposed, y=Count, fill=Phylum))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Microbial Phylum Relative Abundance by Exposure Group", subtitle="Taxa > 10% Relative Abundance", x="Exposure Material", y="Relative Abundance", fill="Phylum")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(b.exp1c,filename = "figures/microb_phyla_10percent_RA_Expo_6.18.21.pdf", width=8, height=6, dpi=600)

#class
b.cls.RA_exp.m1<-subset(b.cls.RA_exp.m, c(Count)>(0.01)) ## DROP BACTERIAL CLASSES that are less than 1% abundant!!!!!!1
b.cls.RA_exp.m5<-subset(b.cls.RA_exp.m, c(Count)>(0.05)) ## DROP BACTERIAL CLASSES that are less than 5% abundant!!!!!!1
b.cls.RA_exp.m10<-subset(b.cls.RA_exp.m, c(Count)>(0.10)) ## DROP BACTERIAL CLASSES that are less than 10% abundant!!!!!!1

b.exp2<-ggplot(b.cls.RA_exp.m, aes(x=Exposed, y=Count, fill=Class))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Microbial Class Relative Abundance by Exposure Group", x="Exposure Material", y="Relative Abundance", fill="Class")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(b.exp2,filename = "figures/microb_class_RA_Expo_6.18.21.pdf", width=10, height=6, dpi=600)

b.exp2a<-ggplot(b.cls.RA_exp.m1, aes(x=Exposed, y=Count, fill=Class))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Microbial Class Relative Abundance by Exposure Group", subtitle="Taxa > 1% Relative Abundance", x="Exposure Material", y="Relative Abundance", fill="Class")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(b.exp2a,filename = "figures/microb_class_1percent_RA_Expo_6.18.21.pdf", width=8, height=6, dpi=600)

b.exp2b<-ggplot(b.cls.RA_exp.m5, aes(x=Exposed, y=Count, fill=Class))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Microbial Class Relative Abundance by Exposure Group", subtitle="Taxa > 5% Relative Abundance", x="Exposure Material", y="Relative Abundance", fill="Class")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(b.exp2b,filename = "figures/microb_class_5percent_RA_Expo_6.18.21.pdf", width=8, height=6, dpi=600)

b.exp2c<-ggplot(b.cls.RA_exp.m10, aes(x=Exposed, y=Count, fill=Class))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Microbial Class Relative Abundance by Exposure Group", subtitle="Taxa > 10% Relative Abundance", x="Exposure Material", y="Relative Abundance", fill="Class")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(b.exp2c,filename = "figures/microb_class_10percent_RA_Expo_6.18.21.pdf", width=8, height=6, dpi=600)

#order
b.ord.RA_exp.m1<-subset(b.ord.RA_exp.m, c(Count)>(0.01)) ## DROP BACTERIAL ORDERS that are less than 1% abundant!!!!!!1
b.ord.RA_exp.m5<-subset(b.ord.RA_exp.m, c(Count)>(0.05)) ## DROP BACTERIAL ORDERS that are less than 5% abundant!!!!!!1
b.ord.RA_exp.m10<-subset(b.ord.RA_exp.m, c(Count)>(0.10)) ## DROP BACTERIAL ORDERS that are less than 10% abundant!!!!!!1

b.exp3<-ggplot(b.ord.RA_exp.m, aes(x=Exposed, y=Count, fill=Order))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Microbial Order Relative Abundance by Exposure Group", x="Exposure Material", y="Relative Abundance", fill="Order")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(b.exp3,filename = "figures/microb_order_RA_Expo_6.18.21.pdf", width=20, height=6, dpi=600)

b.exp3a<-ggplot(b.ord.RA_exp.m1, aes(x=Exposed, y=Count, fill=Order))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Microbial Order Relative Abundance by Exposure Group", subtitle="Taxa > 1% Relative Abundance", x="Exposure Material", y="Relative Abundance", fill="Order")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(b.exp3a,filename = "figures/microb_order_1percent_RA_Expo_6.18.21.pdf", width=10, height=6, dpi=600)

b.exp3b<-ggplot(b.ord.RA_exp.m5, aes(x=Exposed, y=Count, fill=Order))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Microbial Order Relative Abundance by Exposure Group", subtitle="Taxa > 5% Relative Abundance", x="Exposure Material", y="Relative Abundance", fill="Order")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(b.exp3b,filename = "figures/microb_order_5percent_RA_Expo_6.18.21.pdf", width=8, height=6, dpi=600)

b.exp3c<-ggplot(b.ord.RA_exp.m10, aes(x=Exposed, y=Count, fill=Order))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Microbial Order Relative Abundance by Exposure Group", subtitle="Taxa > 10% Relative Abundance", x="Exposure Material", y="Relative Abundance", fill="Order")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(b.exp3c,filename = "figures/microb_order_10percent_RA_Expo_6.18.21.pdf", width=8, height=6, dpi=600)

# family
b.fam.RA_exp.m1<-subset(b.fam.RA_exp.m, c(Count)>(0.01)) ## DROP BACTERIAL FamilyS that are less than 1% abundant!!!!!!1
b.fam.RA_exp.m5<-subset(b.fam.RA_exp.m, c(Count)>(0.05)) ## DROP BACTERIAL FamilyS that are less than 5% abundant!!!!!!1
b.fam.RA_exp.m10<-subset(b.fam.RA_exp.m, c(Count)>(0.10)) ## DROP BACTERIAL FamilyS that are less than 10% abundant!!!!!!1

b.exp4<-ggplot(b.fam.RA_exp.m, aes(x=Exposed, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Microbial Family Relative Abundance by Exposure Group", x="Exposure Material", y="Relative Abundance", fill="Family")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(b.exp4,filename = "figures/microb_Family_RA_Expo_6.18.21.pdf", width=20, height=6, dpi=600)

b.exp4a<-ggplot(b.fam.RA_exp.m1, aes(x=Exposed, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Microbial Family Relative Abundance by Exposure Group", subtitle="Taxa > 1% Relative Abundance", x="Exposure Material", y="Relative Abundance", fill="Family")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(b.exp4a,filename = "figures/microb_Family_1percent_RA_Expo_6.18.21.pdf", width=10, height=6, dpi=600)

b.exp4b<-ggplot(b.fam.RA_exp.m5, aes(x=Exposed, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Microbial Family Relative Abundance by Exposure Group", subtitle="Taxa > 5% Relative Abundance", x="Exposure Material", y="Relative Abundance", fill="Family")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(b.exp4b,filename = "figures/microb_Family_5percent_RA_Expo_6.18.21.pdf", width=8, height=6, dpi=600)

b.exp4c<-ggplot(b.fam.RA_exp.m10, aes(x=Exposed, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Microbial Family Relative Abundance by Exposure Group", subtitle="Taxa > 10% Relative Abundance", x="Exposure Material", y="Relative Abundance", fill="Family")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(b.exp4c,filename = "figures/microb_Family_10percent_RA_Expo_6.18.21.pdf", width=8, height=6, dpi=600)

# Genus
b.gen.RA_exp.m1<-subset(b.gen.RA_exp.m, c(Count)>(0.01)) ## DROP BACTERIAL GENERA that are less than 1% abundant!!!!!!1
b.gen.RA_exp.m5<-subset(b.gen.RA_exp.m, c(Count)>(0.05)) ## DROP BACTERIAL GENERA that are less than 5% abundant!!!!!!1
b.gen.RA_exp.m10<-subset(b.gen.RA_exp.m, c(Count)>(0.10)) ## DROP BACTERIAL GENERA that are less than 10% abundant!!!!!!1

# b.exp5<-ggplot(b.gen.RA_exp.m, aes(x=Exposed, y=Count, fill=Genus))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
#   labs(title = "Microbial Genera Relative Abundance by Exposure Group", x="Exposure Material", y="Relative Abundance", fill="Genus")+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
# ggsave(b.exp5,filename = "figures/microb_Genus_RA_Expo_6.18.21.pdf", width=50, height=6, dpi=600)

b.exp5a<-ggplot(b.gen.RA_exp.m1, aes(x=Exposed, y=Count, fill=Genus))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Microbial Genera Relative Abundance by Exposure Group", subtitle="Taxa > 1% Relative Abundance", x="Exposure Material", y="Relative Abundance", fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(b.exp5a,filename = "figures/microb_Genus_1percent_RA_Expo_6.18.21.pdf", width=20, height=6, dpi=600)

b.exp5b<-ggplot(b.gen.RA_exp.m5, aes(x=Exposed, y=Count, fill=Genus))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Microbial Genera Relative Abundance by Exposure Group", subtitle="Taxa > 5% Relative Abundance", x="Exposure Material", y="Relative Abundance", fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(b.exp5b,filename = "figures/microb_Genus_5percent_RA_Expo_6.18.21.pdf", width=10, height=6, dpi=600)

b.exp5c<-ggplot(b.gen.RA_exp.m10, aes(x=Exposed, y=Count, fill=Genus))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Microbial Genera Relative Abundance by Exposure Group", subtitle="Taxa > 10% Relative Abundance", x="Exposure Material", y="Relative Abundance", fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(b.exp5c,filename = "figures/microb_Genus_10percent_RA_Expo_6.18.21.pdf", width=8, height=6, dpi=600)

# species
b.spec.RA_exp.m1<-subset(b.spec.RA_exp.m, c(Count)>(0.01)) ## DROP BACTERIAL Species that are less than 1% abundant!!!!!!1
b.spec.RA_exp.m5<-subset(b.spec.RA_exp.m, c(Count)>(0.05)) ## DROP BACTERIAL Species that are less than 5% abundant!!!!!!1
b.spec.RA_exp.m10<-subset(b.spec.RA_exp.m, c(Count)>(0.10)) ## DROP BACTERIAL Species that are less than 10% abundant!!!!!!1

# b.exp5<-ggplot(b.spec.RA_exp.m, aes(x=Exposed, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
#   labs(title = "Microbial Species Relative Abundance by Exposure Group", x="Exposure Material", y="Relative Abundance", fill="Species")+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
# ggsave(b.exp5,filename = "figures/microb_Species_RA_Expo_6.18.21.pdf", width=, height=6, dpi=600)

b.exp5a<-ggplot(b.spec.RA_exp.m1, aes(x=Exposed, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Microbial Species Relative Abundance by Exposure Group", subtitle="Taxa > 1% Relative Abundance", x="Exposure Material", y="Relative Abundance", fill="Species")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(b.exp5a,filename = "figures/microb_Species_1percent_RA_Expo_6.18.21.pdf", width=20, height=6, dpi=600)

b.exp5b<-ggplot(b.spec.RA_exp.m5, aes(x=Exposed, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Microbial Species Relative Abundance by Exposure Group", subtitle="Taxa > 5% Relative Abundance", x="Exposure Material", y="Relative Abundance", fill="Species")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(b.exp5b,filename = "figures/microb_Species_5percent_RA_Expo_6.18.21.pdf", width=8, height=6, dpi=600)

b.exp5c<-ggplot(b.spec.RA_exp.m10, aes(x=Exposed, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Microbial Species Relative Abundance by Exposure Group", subtitle="Taxa > 10% Relative Abundance", caption="* Species' relative abundance in Controls did not surpass 10%",x="Exposure Material", y="Relative Abundance", fill="Species")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(b.exp5c,filename = "figures/microb_Species_10percent_RA_Expo_6.18.21.pdf", width=8, height=6, dpi=600)


#### ITS2 Relative Abundance by Exposure Group ####
head(f.dat.m)
its2.asv_meta<-merge(f.dat.m, metadata, by="SampleID")

## phylum ....
head(its2.asv_meta)

# by phylum + exposure
its2_phy_exp<- as.data.frame(dcast(its2.asv_meta,Exposed~Phylum, value.var="Count", fun.aggregate=sum)) ###
head(its2_phy_exp)
rownames(its2_phy_exp)<-its2_phy_exp$Exposed
head(its2_phy_exp)

f.phy.RA_exp<-data.frame(decostand(its2_phy_exp[,-c(1)], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(f.phy.RA_exp)
head(f.phy.RA_exp)

f.phy.RA_exp$Exposed<-rownames(f.phy.RA_exp)
f.phy.RA_exp.m<-melt(f.phy.RA_exp)
head(f.phy.RA_exp.m)
colnames(f.phy.RA_exp.m)[which(names(f.phy.RA_exp.m) == "variable")] <- "Phylum"
colnames(f.phy.RA_exp.m)[which(names(f.phy.RA_exp.m) == "value")] <- "Count"
head(f.phy.RA_exp.m) ## relative abundance based on sum of counts by phyla in exposure group!

# by class + exposure
its2_cls_exp <- as.data.frame(dcast(its2.asv_meta,Exposed~Class, value.var="Count", fun.aggregate=sum)) ###
head(its2_cls_exp) # counts by class + elevation
rownames(its2_cls_exp)<-its2_cls_exp$Exposed

f.cls.RA_exp<-data.frame(decostand(its2_cls_exp[,-c(1)], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(f.cls.RA_exp)
head(f.cls.RA_exp)

f.cls.RA_exp$Exposed<-rownames(f.cls.RA_exp)
f.cls.RA_exp.m<-melt(f.cls.RA_exp)
head(f.cls.RA_exp.m)
colnames(f.cls.RA_exp.m)[which(names(f.cls.RA_exp.m) == "variable")] <- "Class"
colnames(f.cls.RA_exp.m)[which(names(f.cls.RA_exp.m) == "value")] <- "Count"
head(f.cls.RA_exp.m) ## relative abundance based on sum of counts by class in exposure group!

# by order + exposure
its2_ord_exp <- as.data.frame(dcast(its2.asv_meta,Exposed~Order, value.var="Count", fun.aggregate=sum)) ###
head(its2_ord_exp) # counts by Order + elevation
rownames(its2_ord_exp)<-its2_ord_exp$Exposed

f.ord.RA_exp<-data.frame(decostand(its2_ord_exp[,-c(1)], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(f.ord.RA_exp)
head(f.ord.RA_exp)

f.ord.RA_exp$Exposed<-rownames(f.ord.RA_exp)
f.ord.RA_exp.m<-melt(f.ord.RA_exp)
head(f.ord.RA_exp.m)
colnames(f.ord.RA_exp.m)[which(names(f.ord.RA_exp.m) == "variable")] <- "Order"
colnames(f.ord.RA_exp.m)[which(names(f.ord.RA_exp.m) == "value")] <- "Count"
head(f.ord.RA_exp.m) ## relative abundance based on sum of counts by order by exposed group!

# by Family + exposure
its2_fam_exp <- as.data.frame(dcast(its2.asv_meta,Exposed~Family, value.var="Count", fun.aggregate=sum)) ###
head(its2_fam_exp) # counts by Family + elevation
rownames(its2_fam_exp)<-its2_fam_exp$Exposed

f.fam.RA_exp<-data.frame(decostand(its2_fam_exp[,-c(1)], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(f.fam.RA_exp)
head(f.fam.RA_exp)

f.fam.RA_exp$Exposed<-rownames(f.fam.RA_exp)
f.fam.RA_exp.m<-melt(f.fam.RA_exp)
head(f.fam.RA_exp.m)
colnames(f.fam.RA_exp.m)[which(names(f.fam.RA_exp.m) == "variable")] <- "Family"
colnames(f.fam.RA_exp.m)[which(names(f.fam.RA_exp.m) == "value")] <- "Count"
head(f.fam.RA_exp.m) ## relative abundance based on sum of counts by family by exposed group!

# by Genus + exposure
its2_gen_exp <- as.data.frame(dcast(its2.asv_meta,Exposed~Genus, value.var="Count", fun.aggregate=sum)) ###
head(its2_gen_exp) # counts by genus + elevation
rownames(its2_gen_exp)<-its2_gen_exp$Exposed

f.gen.RA_exp<-data.frame(decostand(its2_gen_exp[,-c(1)], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(f.gen.RA_exp)
head(f.gen.RA_exp)

f.gen.RA_exp$Exposed<-rownames(f.gen.RA_exp)
f.gen.RA_exp.m<-melt(f.gen.RA_exp)
head(f.gen.RA_exp.m)
colnames(f.gen.RA_exp.m)[which(names(f.gen.RA_exp.m) == "variable")] <- "Genus"
colnames(f.gen.RA_exp.m)[which(names(f.gen.RA_exp.m) == "value")] <- "Count"
head(f.gen.RA_exp.m) ## relative abundance based on sum of counts by genus by exposed group!

# by Genus/Species + exposure
its2_spec_exp <- as.data.frame(dcast(its2.asv_meta,Exposed~Genus+Species, value.var="Count", fun.aggregate=sum)) ###
head(its2_spec_exp) # counts by genus_species + elevation
rownames(its2_spec_exp)<-its2_spec_exp$Exposed

f.spec.RA_exp<-data.frame(decostand(its2_spec_exp[,-c(1)], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(f.spec.RA_exp)
head(f.spec.RA_exp)

f.spec.RA_exp$Exposed<-rownames(f.spec.RA_exp)
f.spec.RA_exp.m<-melt(f.spec.RA_exp)
head(f.spec.RA_exp.m)
colnames(f.spec.RA_exp.m)[which(names(f.spec.RA_exp.m) == "variable")] <- "Genus_species"
colnames(f.spec.RA_exp.m)[which(names(f.spec.RA_exp.m) == "value")] <- "Count"
f.spec.RA_exp.m$Genus_species<-gsub("_", " ", f.spec.RA_exp.m$Genus_species) ## gsub is global sub (does not just remove first instance of pattern, but multiple)

head(f.spec.RA_exp.m) ## relative abundance based on sum of counts by specus by exposed group!

#### Visualize ITS2 Relative Abundance by Exposure Group ####

# phylum
f.phy.RA_exp.m1<-subset(f.phy.RA_exp.m, c(Count)>(0.01)) ## DROP FUNGAL PHYLA that are less than 1% abundant!!!!!!1
f.phy.RA_exp.m5<-subset(f.phy.RA_exp.m, c(Count)>(0.05)) ## DROP FUNGAL PHYLA that are less than 5% abundant!!!!!!1
f.phy.RA_exp.m10<-subset(f.phy.RA_exp.m, c(Count)>(0.10)) ## DROP FUNGAL PHYLA that are less than 10% abundant!!!!!!1

f.exp1<-ggplot(f.phy.RA_exp.m, aes(x=Exposed, y=Count, fill=Phylum))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Fungal Phylum Relative Abundance by Exposure Group", x="Exposure Material", y="Relative Abundance", fill="Phylum")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(f.exp1,filename = "figures/fungi_phyla_RA_Expo_6.18.21.pdf", width=8, height=6, dpi=600)

f.exp1a<-ggplot(f.phy.RA_exp.m1, aes(x=Exposed, y=Count, fill=Phylum))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Fungal Phylum Relative Abundance by Exposure Group", subtitle="Taxa > 1% Relative Abundance", x="Exposure Material", y="Relative Abundance", fill="Phylum")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(f.exp1a,filename = "figures/fungi_phyla_1percent_RA_Expo_6.18.21.pdf", width=8, height=6, dpi=600)

f.exp1b<-ggplot(f.phy.RA_exp.m5, aes(x=Exposed, y=Count, fill=Phylum))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Fungal Phylum Relative Abundance by Exposure Group", subtitle="Taxa > 5% Relative Abundance", x="Exposure Material", y="Relative Abundance", fill="Phylum")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(f.exp1b,filename = "figures/fungi_phyla_5percent_RA_Expo_6.18.21.pdf", width=8, height=6, dpi=600)

f.exp1c<-ggplot(f.phy.RA_exp.m10, aes(x=Exposed, y=Count, fill=Phylum))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Fungal Phylum Relative Abundance by Exposure Group", subtitle="Taxa > 10% Relative Abundance", x="Exposure Material", y="Relative Abundance", fill="Phylum")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(f.exp1c,filename = "figures/fungi_phyla_10percent_RA_Expo_6.18.21.pdf", width=8, height=6, dpi=600)

#class
f.cls.RA_exp.m1<-subset(f.cls.RA_exp.m, c(Count)>(0.01)) ## DROP FUNGAL CLASSES that are less than 1% abundant!!!!!!1
f.cls.RA_exp.m5<-subset(f.cls.RA_exp.m, c(Count)>(0.05)) ## DROP FUNGAL CLASSES that are less than 5% abundant!!!!!!1
f.cls.RA_exp.m10<-subset(f.cls.RA_exp.m, c(Count)>(0.10)) ## DROP FUNGAL CLASSES that are less than 10% abundant!!!!!!1

f.exp2<-ggplot(f.cls.RA_exp.m, aes(x=Exposed, y=Count, fill=Class))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Fungal Class Relative Abundance by Exposure Group", x="Exposure Material", y="Relative Abundance", fill="Class")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(f.exp2,filename = "figures/fungi_class_RA_Expo_6.18.21.pdf", width=8, height=6, dpi=600)

f.exp2a<-ggplot(f.cls.RA_exp.m1, aes(x=Exposed, y=Count, fill=Class))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Fungal Class Relative Abundance by Exposure Group", subtitle="Taxa > 1% Relative Abundance", x="Exposure Material", y="Relative Abundance", fill="Class")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(f.exp2a,filename = "figures/fungi_class_1percent_RA_Expo_6.18.21.pdf", width=8, height=6, dpi=600)

f.exp2b<-ggplot(f.cls.RA_exp.m5, aes(x=Exposed, y=Count, fill=Class))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Fungal Class Relative Abundance by Exposure Group", subtitle="Taxa > 5% Relative Abundance", x="Exposure Material", y="Relative Abundance", fill="Class")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(f.exp2b,filename = "figures/fungi_class_5percent_RA_Expo_6.18.21.pdf", width=8, height=6, dpi=600)

f.exp2c<-ggplot(f.cls.RA_exp.m10, aes(x=Exposed, y=Count, fill=Class))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Fungal Class Relative Abundance by Exposure Group", subtitle="Taxa > 10% Relative Abundance", x="Exposure Material", y="Relative Abundance", fill="Class")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(f.exp2c,filename = "figures/fungi_class_10percent_RA_Expo_6.18.21.pdf", width=8, height=6, dpi=600)

#order
f.ord.RA_exp.m1<-subset(f.ord.RA_exp.m, c(Count)>(0.01)) ## DROP FUNGAL ORDERS that are less than 1% abundant!!!!!!1
f.ord.RA_exp.m5<-subset(f.ord.RA_exp.m, c(Count)>(0.05)) ## DROP FUNGAL ORDERS that are less than 5% abundant!!!!!!1
f.ord.RA_exp.m10<-subset(f.ord.RA_exp.m, c(Count)>(0.10)) ## DROP FUNGAL ORDERS that are less than 10% abundant!!!!!!1

f.exp3<-ggplot(f.ord.RA_exp.m, aes(x=Exposed, y=Count, fill=Order))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Fungal Order Relative Abundance by Exposure Group", x="Exposure Material", y="Relative Abundance", fill="Order")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(f.exp3,filename = "figures/fungi_order_RA_Expo_6.18.21.pdf", width=8, height=6, dpi=600)

f.exp3a<-ggplot(f.ord.RA_exp.m1, aes(x=Exposed, y=Count, fill=Order))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Fungal Order Relative Abundance by Exposure Group", subtitle="Taxa > 1% Relative Abundance", x="Exposure Material", y="Relative Abundance", fill="Order")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(f.exp3a,filename = "figures/fungi_order_1percent_RA_Expo_6.18.21.pdf", width=8, height=6, dpi=600)

f.exp3b<-ggplot(f.ord.RA_exp.m5, aes(x=Exposed, y=Count, fill=Order))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Fungal Order Relative Abundance by Exposure Group", subtitle="Taxa > 5% Relative Abundance", x="Exposure Material", y="Relative Abundance", fill="Order")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(f.exp3b,filename = "figures/fungi_order_5percent_RA_Expo_6.18.21.pdf", width=8, height=6, dpi=600)

f.exp3c<-ggplot(f.ord.RA_exp.m10, aes(x=Exposed, y=Count, fill=Order))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Fungal Order Relative Abundance by Exposure Group", subtitle="Taxa > 10% Relative Abundance", x="Exposure Material", y="Relative Abundance", fill="Order")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(f.exp3c,filename = "figures/fungi_order_10percent_RA_Expo_6.18.21.pdf", width=8, height=6, dpi=600)

# family
f.fam.RA_exp.m1<-subset(f.fam.RA_exp.m, c(Count)>(0.01)) ## DROP FUNGAL FamilyS that are less than 1% abundant!!!!!!1
f.fam.RA_exp.m5<-subset(f.fam.RA_exp.m, c(Count)>(0.05)) ## DROP FUNGAL FamilyS that are less than 5% abundant!!!!!!1
f.fam.RA_exp.m10<-subset(f.fam.RA_exp.m, c(Count)>(0.10)) ## DROP FUNGAL FamilyS that are less than 10% abundant!!!!!!1

f.exp4<-ggplot(f.fam.RA_exp.m, aes(x=Exposed, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Fungal Family Relative Abundance by Exposure Group", x="Exposure Material", y="Relative Abundance", fill="Family")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(f.exp4,filename = "figures/fungi_Family_RA_Expo_6.18.21.pdf", width=8, height=6, dpi=600)

f.exp4a<-ggplot(f.fam.RA_exp.m1, aes(x=Exposed, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Fungal Family Relative Abundance by Exposure Group", subtitle="Taxa > 1% Relative Abundance", x="Exposure Material", y="Relative Abundance", fill="Family")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(f.exp4a,filename = "figures/fungi_Family_1percent_RA_Expo_6.18.21.pdf", width=8, height=6, dpi=600)

f.exp4b<-ggplot(f.fam.RA_exp.m5, aes(x=Exposed, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Fungal Family Relative Abundance by Exposure Group", subtitle="Taxa > 5% Relative Abundance", x="Exposure Material", y="Relative Abundance", fill="Family")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(f.exp4b,filename = "figures/fungi_Family_5percent_RA_Expo_6.18.21.pdf", width=8, height=6, dpi=600)

f.exp4c<-ggplot(f.fam.RA_exp.m10, aes(x=Exposed, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Fungal Family Relative Abundance by Exposure Group", subtitle="Taxa > 10% Relative Abundance", x="Exposure Material", y="Relative Abundance", fill="Family")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(f.exp4c,filename = "figures/fungi_Family_10percent_RA_Expo_6.18.21.pdf", width=8, height=6, dpi=600)

# Genus
f.gen.RA_exp.m1<-subset(f.gen.RA_exp.m, c(Count)>(0.01)) ## DROP FUNGAL GENERA that are less than 1% abundant!!!!!!1
f.gen.RA_exp.m5<-subset(f.gen.RA_exp.m, c(Count)>(0.05)) ## DROP FUNGAL GENERA that are less than 5% abundant!!!!!!1
f.gen.RA_exp.m10<-subset(f.gen.RA_exp.m, c(Count)>(0.10)) ## DROP FUNGAL GENERA that are less than 10% abundant!!!!!!1

f.exp5<-ggplot(f.gen.RA_exp.m, aes(x=Exposed, y=Count, fill=Genus))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Fungal Genera Relative Abundance by Exposure Group", x="Exposure Material", y="Relative Abundance", fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(f.exp5,filename = "figures/fungi_Genus_RA_Expo_6.18.21.pdf", width=8, height=6, dpi=600)

f.exp5a<-ggplot(f.gen.RA_exp.m1, aes(x=Exposed, y=Count, fill=Genus))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Fungal Genera Relative Abundance by Exposure Group", subtitle="Taxa > 1% Relative Abundance", x="Exposure Material", y="Relative Abundance", fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(f.exp5a,filename = "figures/fungi_Genus_1percent_RA_Expo_6.18.21.pdf", width=8, height=6, dpi=600)

f.exp5b<-ggplot(f.gen.RA_exp.m5, aes(x=Exposed, y=Count, fill=Genus))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Fungal Genera Relative Abundance by Exposure Group", subtitle="Taxa > 5% Relative Abundance", x="Exposure Material", y="Relative Abundance", fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(f.exp5b,filename = "figures/fungi_Genus_5percent_RA_Expo_6.18.21.pdf", width=8, height=6, dpi=600)

f.exp5c<-ggplot(f.gen.RA_exp.m10, aes(x=Exposed, y=Count, fill=Genus))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Fungal Genera Relative Abundance by Exposure Group", subtitle="Taxa > 10% Relative Abundance", x="Exposure Material", y="Relative Abundance", fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(f.exp5c,filename = "figures/fungi_Genus_10percent_RA_Expo_6.18.21.pdf", width=8, height=6, dpi=600)

# species
f.spec.RA_exp.m1<-subset(f.spec.RA_exp.m, c(Count)>(0.01)) ## DROP FUNGAL Species that are less than 1% abundant!!!!!!1
f.spec.RA_exp.m5<-subset(f.spec.RA_exp.m, c(Count)>(0.05)) ## DROP FUNGAL Species that are less than 5% abundant!!!!!!1
f.spec.RA_exp.m10<-subset(f.spec.RA_exp.m, c(Count)>(0.10)) ## DROP FUNGAL Species that are less than 10% abundant!!!!!!1

f.exp5<-ggplot(f.spec.RA_exp.m, aes(x=Exposed, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Fungal Species Relative Abundance by Exposure Group", x="Exposure Material", y="Relative Abundance", fill="Species")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(f.exp5,filename = "figures/fungi_Species_RA_Expo_6.18.21.pdf", width=8, height=6, dpi=600)

f.exp5a<-ggplot(f.spec.RA_exp.m1, aes(x=Exposed, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Fungal Species Relative Abundance by Exposure Group", subtitle="Taxa > 1% Relative Abundance", x="Exposure Material", y="Relative Abundance", fill="Species")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(f.exp5a,filename = "figures/fungi_Species_1percent_RA_Expo_6.18.21.pdf", width=8, height=6, dpi=600)

f.exp5b<-ggplot(f.spec.RA_exp.m5, aes(x=Exposed, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Fungal Species Relative Abundance by Exposure Group", subtitle="Taxa > 5% Relative Abundance", x="Exposure Material", y="Relative Abundance", fill="Species")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(f.exp5b,filename = "figures/fungi_Species_5percent_RA_Expo_6.18.21.pdf", width=8, height=6, dpi=600)

f.exp5c<-ggplot(f.spec.RA_exp.m10, aes(x=Exposed, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
  labs(title = "Fungal Species Relative Abundance by Exposure Group", subtitle="Taxa > 10% Relative Abundance", x="Exposure Material", y="Relative Abundance", fill="Species")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(f.exp5c,filename = "figures/fungi_Species_10percent_RA_Expo_6.18.21.pdf", width=8, height=6, dpi=600)



#### 16S Relative Abundance by Exposed: Yes or No ####
# head(fbdat.m)
# bac.asv_meta<-merge(b.dat.m, metadata, by="SampleID")

## phylum ....
head(bac.asv_meta)

# by phylum + exposure
bac_phy_e.y.n<- as.data.frame(dcast(bac.asv_meta,exposed_y.n~Phylum, value.var="Count", fun.aggregate=sum)) ###
head(bac_phy_e.y.n)
rownames(bac_phy_e.y.n)<-bac_phy_e.y.n$exposed_y.n
head(bac_phy_e.y.n)

b.phy.RA_e.y.n<-data.frame(decostand(bac_phy_e.y.n[,-c(1)], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(b.phy.RA_e.y.n)
head(b.phy.RA_e.y.n)

b.phy.RA_e.y.n$exposed_y.n<-rownames(b.phy.RA_e.y.n)
b.phy.RA_e.y.n.m<-melt(b.phy.RA_e.y.n)
head(b.phy.RA_e.y.n.m)
colnames(b.phy.RA_e.y.n.m)[which(names(b.phy.RA_e.y.n.m) == "variable")] <- "Phylum"
colnames(b.phy.RA_e.y.n.m)[which(names(b.phy.RA_e.y.n.m) == "value")] <- "Count"
head(b.phy.RA_e.y.n.m) ## relative abundance based on sum of counts by phyla in exposed/not exposed!

# by class + exposure
bac_cls_e.y.n <- as.data.frame(dcast(bac.asv_meta,exposed_y.n~Class, value.var="Count", fun.aggregate=sum)) ###
head(bac_cls_e.y.n) # counts by class + elevation
rownames(bac_cls_e.y.n)<-bac_cls_e.y.n$exposed_y.n

b.cls.RA_e.y.n<-data.frame(decostand(bac_cls_e.y.n[,-c(1)], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(b.cls.RA_e.y.n)
head(b.cls.RA_e.y.n)

b.cls.RA_e.y.n$exposed_y.n<-rownames(b.cls.RA_e.y.n)
b.cls.RA_e.y.n.m<-melt(b.cls.RA_e.y.n)
head(b.cls.RA_e.y.n.m)
colnames(b.cls.RA_e.y.n.m)[which(names(b.cls.RA_e.y.n.m) == "variable")] <- "Class"
colnames(b.cls.RA_e.y.n.m)[which(names(b.cls.RA_e.y.n.m) == "value")] <- "Count"
head(b.cls.RA_e.y.n.m) ## relative abundance based on sum of counts by class in exposed/not exposed!

# by order + exposure
bac_ord_e.y.n <- as.data.frame(dcast(bac.asv_meta,exposed_y.n~Order, value.var="Count", fun.aggregate=sum)) ###
head(bac_ord_e.y.n) # counts by Order + elevation
rownames(bac_ord_e.y.n)<-bac_ord_e.y.n$exposed_y.n

b.ord.RA_e.y.n<-data.frame(decostand(bac_ord_e.y.n[,-c(1)], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(b.ord.RA_e.y.n)
head(b.ord.RA_e.y.n)

b.ord.RA_e.y.n$exposed_y.n<-rownames(b.ord.RA_e.y.n)
b.ord.RA_e.y.n.m<-melt(b.ord.RA_e.y.n)
head(b.ord.RA_e.y.n.m)
colnames(b.ord.RA_e.y.n.m)[which(names(b.ord.RA_e.y.n.m) == "variable")] <- "Order"
colnames(b.ord.RA_e.y.n.m)[which(names(b.ord.RA_e.y.n.m) == "value")] <- "Count"
head(b.ord.RA_e.y.n.m) ## relative abundance based on sum of counts by order by exposed group!

# by Family + exposure
bac_fam_e.y.n <- as.data.frame(dcast(bac.asv_meta,exposed_y.n~Family, value.var="Count", fun.aggregate=sum)) ###
head(bac_fam_e.y.n) # counts by Family + elevation
rownames(bac_fam_e.y.n)<-bac_fam_e.y.n$exposed_y.n

b.fam.RA_e.y.n<-data.frame(decostand(bac_fam_e.y.n[,-c(1)], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(b.fam.RA_e.y.n)
head(b.fam.RA_e.y.n)

b.fam.RA_e.y.n$exposed_y.n<-rownames(b.fam.RA_e.y.n)
b.fam.RA_e.y.n.m<-melt(b.fam.RA_e.y.n)
head(b.fam.RA_e.y.n.m)
colnames(b.fam.RA_e.y.n.m)[which(names(b.fam.RA_e.y.n.m) == "variable")] <- "Family"
colnames(b.fam.RA_e.y.n.m)[which(names(b.fam.RA_e.y.n.m) == "value")] <- "Count"
head(b.fam.RA_e.y.n.m) ## relative abundance based on sum of counts by family by exposed group!

# by Genus + exposure
bac_gen_e.y.n <- as.data.frame(dcast(bac.asv_meta,exposed_y.n~Genus, value.var="Count", fun.aggregate=sum)) ###
head(bac_gen_e.y.n) # counts by genus + elevation
rownames(bac_gen_e.y.n)<-bac_gen_e.y.n$exposed_y.n

b.gen.RA_e.y.n<-data.frame(decostand(bac_gen_e.y.n[,-c(1)], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(b.gen.RA_e.y.n)
head(b.gen.RA_e.y.n)

b.gen.RA_e.y.n$exposed_y.n<-rownames(b.gen.RA_e.y.n)
b.gen.RA_e.y.n.m<-melt(b.gen.RA_e.y.n)
head(b.gen.RA_e.y.n.m)
colnames(b.gen.RA_e.y.n.m)[which(names(b.gen.RA_e.y.n.m) == "variable")] <- "Genus"
colnames(b.gen.RA_e.y.n.m)[which(names(b.gen.RA_e.y.n.m) == "value")] <- "Count"
b.gen.RA_e.y.n.m$Genus_Species<-gsub("_", " ", b.gen.RA_e.y.n.m$Genus_Species) ## gsub is global sub (does not just remove first instance of pattern, but multiple)

head(b.gen.RA_e.y.n.m) ## relative abundance based on sum of counts by genus by exposed group!

# by Genus/Species + exposure
bac_spec_e.y.n <- as.data.frame(dcast(bac.asv_meta,exposed_y.n~Genus+Species, value.var="Count", fun.aggregate=sum)) ###
head(bac_spec_e.y.n) # counts by genus_species + elevation
rownames(bac_spec_e.y.n)<-bac_spec_e.y.n$exposed_y.n

b.spec.RA_e.y.n<-data.frame(decostand(bac_spec_e.y.n[,-c(1)], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(b.spec.RA_e.y.n)
head(b.spec.RA_e.y.n)

b.spec.RA_e.y.n$exposed_y.n<-rownames(b.spec.RA_e.y.n)
b.spec.RA_e.y.n.m<-melt(b.spec.RA_e.y.n)
head(b.spec.RA_e.y.n.m)
colnames(b.spec.RA_e.y.n.m)[which(names(b.spec.RA_e.y.n.m) == "variable")] <- "Genus_species"
colnames(b.spec.RA_e.y.n.m)[which(names(b.spec.RA_e.y.n.m) == "value")] <- "Count"
b.spec.RA_e.y.n.m$Genus_species<-gsub("_", " ", b.spec.RA_e.y.n.m$Genus_species)
head(b.spec.RA_e.y.n.m) ## relative abundance based on sum of counts by specus by exposed group!

#### Visualize 16S Relative Abundance by Exposed/ Not Exposed ####

# phylum
b.phy.RA_e.y.n.m1<-subset(b.phy.RA_e.y.n.m, c(Count)>(0.01)) ## DROP BACTERIAL PHYLA that are less than 1% abundant!!!!!!1
b.phy.RA_e.y.n.m5<-subset(b.phy.RA_e.y.n.m, c(Count)>(0.05)) ## DROP BACTERIAL PHYLA that are less than 5% abundant!!!!!!1
b.phy.RA_e.y.n.m10<-subset(b.phy.RA_e.y.n.m, c(Count)>(0.10)) ## DROP BACTERIAL PHYLA that are less than 10% abundant!!!!!!1

b.e.y.n1<-ggplot(b.phy.RA_e.y.n.m, aes(x=exposed_y.n, y=Count, fill=Phylum))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
  labs(title = "Microbial Phylum Relative Abundance by Exposed Condition", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Phylum")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(b.e.y.n1,filename = "figures/microb_phyla_RA_Exp.YesNo_6.18.21.pdf", width=8, height=6, dpi=600)

b.e.y.n1a<-ggplot(b.phy.RA_e.y.n.m1, aes(x=exposed_y.n, y=Count, fill=Phylum))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
  labs(title = "Microbial Phylum Relative Abundance by Exposed Condition", subtitle="Taxa > 1% Relative Abundance", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Phylum")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(b.e.y.n1a,filename = "figures/microb_phyla_1percent_RA_Exp.YesNo_6.18.21.pdf", width=8, height=6, dpi=600)

b.e.y.n1b<-ggplot(b.phy.RA_e.y.n.m5, aes(x=exposed_y.n, y=Count, fill=Phylum))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
  labs(title = "Microbial Phylum Relative Abundance by Exposed Condition", subtitle="Taxa > 5% Relative Abundance", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Phylum")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(b.e.y.n1b,filename = "figures/microb_phyla_5percent_RA_Exp.YesNo_6.18.21.pdf", width=8, height=6, dpi=600)

b.e.y.n1c<-ggplot(b.phy.RA_e.y.n.m10, aes(x=exposed_y.n, y=Count, fill=Phylum))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
  labs(title = "Microbial Phylum Relative Abundance by Exposed Condition", subtitle="Taxa > 10% Relative Abundance", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Phylum")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(b.e.y.n1c,filename = "figures/microb_phyla_10percent_RA_Exp.YesNo_6.18.21.pdf", width=8, height=6, dpi=600)

#class
b.cls.RA_e.y.n.m1<-subset(b.cls.RA_e.y.n.m, c(Count)>(0.01)) ## DROP BACTERIAL CLASSES that are less than 1% abundant!!!!!!1
b.cls.RA_e.y.n.m5<-subset(b.cls.RA_e.y.n.m, c(Count)>(0.05)) ## DROP BACTERIAL CLASSES that are less than 5% abundant!!!!!!1
b.cls.RA_e.y.n.m10<-subset(b.cls.RA_e.y.n.m, c(Count)>(0.10)) ## DROP BACTERIAL CLASSES that are less than 10% abundant!!!!!!1

b.e.y.n2<-ggplot(b.cls.RA_e.y.n.m, aes(x=exposed_y.n, y=Count, fill=Class))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
  labs(title = "Microbial Class Relative Abundance by Exposed Condition", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Class")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(b.e.y.n2,filename = "figures/microb_class_RA_Exp.YesNo_6.18.21.pdf", width=10, height=6, dpi=600)

b.e.y.n2a<-ggplot(b.cls.RA_e.y.n.m1, aes(x=exposed_y.n, y=Count, fill=Class))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
  labs(title = "Microbial Class Relative Abundance by Exposed Condition", subtitle="Taxa > 1% Relative Abundance", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Class")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(b.e.y.n2a,filename = "figures/microb_class_1percent_RA_Exp.YesNo_6.18.21.pdf", width=8, height=6, dpi=600)

b.e.y.n2b<-ggplot(b.cls.RA_e.y.n.m5, aes(x=exposed_y.n, y=Count, fill=Class))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
  labs(title = "Microbial Class Relative Abundance by Exposed Condition", subtitle="Taxa > 5% Relative Abundance", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Class")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(b.e.y.n2b,filename = "figures/microb_class_5percent_RA_Exp.YesNo_6.18.21.pdf", width=8, height=6, dpi=600)

b.e.y.n2c<-ggplot(b.cls.RA_e.y.n.m10, aes(x=exposed_y.n, y=Count, fill=Class))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
  labs(title = "Microbial Class Relative Abundance by Exposed Condition", subtitle="Taxa > 10% Relative Abundance", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Class")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(b.e.y.n2c,filename = "figures/microb_class_10percent_RA_Exp.YesNo_6.18.21.pdf", width=8, height=6, dpi=600)

#order
b.ord.RA_e.y.n.m1<-subset(b.ord.RA_e.y.n.m, c(Count)>(0.01)) ## DROP BACTERIAL ORDERS that are less than 1% abundant!!!!!!1
b.ord.RA_e.y.n.m5<-subset(b.ord.RA_e.y.n.m, c(Count)>(0.05)) ## DROP BACTERIAL ORDERS that are less than 5% abundant!!!!!!1
b.ord.RA_e.y.n.m10<-subset(b.ord.RA_e.y.n.m, c(Count)>(0.10)) ## DROP BACTERIAL ORDERS that are less than 10% abundant!!!!!!1

b.e.y.n3<-ggplot(b.ord.RA_e.y.n.m, aes(x=exposed_y.n, y=Count, fill=Order))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
  labs(title = "Microbial Order Relative Abundance by Exposed Condition", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Order")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(b.e.y.n3,filename = "figures/microb_order_RA_Exp.YesNo_6.18.21.pdf", width=18, height=6, dpi=600)

b.e.y.n3a<-ggplot(b.ord.RA_e.y.n.m1, aes(x=exposed_y.n, y=Count, fill=Order))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
  labs(title = "Microbial Order Relative Abundance by Exposed Condition", subtitle="Taxa > 1% Relative Abundance", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Order")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(b.e.y.n3a,filename = "figures/microb_order_1percent_RA_Exp.YesNo_6.18.21.pdf", width=10, height=6, dpi=600)

b.e.y.n3b<-ggplot(b.ord.RA_e.y.n.m5, aes(x=exposed_y.n, y=Count, fill=Order))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
  labs(title = "Microbial Order Relative Abundance by Exposed Condition", subtitle="Taxa > 5% Relative Abundance", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Order")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(b.e.y.n3b,filename = "figures/microb_order_5percent_RA_Exp.YesNo_6.18.21.pdf", width=8, height=6, dpi=600)

b.e.y.n3c<-ggplot(b.ord.RA_e.y.n.m10, aes(x=exposed_y.n, y=Count, fill=Order))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
  labs(title = "Microbial Order Relative Abundance by Exposed Condition", subtitle="Taxa > 10% Relative Abundance", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Order")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(b.e.y.n3c,filename = "figures/microb_order_10percent_RA_Exp.YesNo_6.18.21.pdf", width=8, height=6, dpi=600)

# family
b.fam.RA_e.y.n.m1<-subset(b.fam.RA_e.y.n.m, c(Count)>(0.01)) ## DROP BACTERIAL FamilyS that are less than 1% abundant!!!!!!1
b.fam.RA_e.y.n.m5<-subset(b.fam.RA_e.y.n.m, c(Count)>(0.05)) ## DROP BACTERIAL FamilyS that are less than 5% abundant!!!!!!1
b.fam.RA_e.y.n.m10<-subset(b.fam.RA_e.y.n.m, c(Count)>(0.10)) ## DROP BACTERIAL FamilyS that are less than 10% abundant!!!!!!1

b.e.y.n4<-ggplot(b.fam.RA_e.y.n.m, aes(x=exposed_y.n, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
  labs(title = "Microbial Family Relative Abundance by Exposed Condition", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Family")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(b.e.y.n4,filename = "figures/microb_Family_RA_Exp.YesNo_6.18.21.pdf", width=22, height=6, dpi=600)

b.e.y.n4a<-ggplot(b.fam.RA_e.y.n.m1, aes(x=exposed_y.n, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
  labs(title = "Microbial Family Relative Abundance by Exposed Condition", subtitle="Taxa > 1% Relative Abundance", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Family")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(b.e.y.n4a,filename = "figures/microb_Family_1percent_RA_Exp.YesNo_6.18.21.pdf", width=10, height=6, dpi=600)

b.e.y.n4b<-ggplot(b.fam.RA_e.y.n.m5, aes(x=exposed_y.n, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
  labs(title = "Microbial Family Relative Abundance by Exposed Condition", subtitle="Taxa > 5% Relative Abundance", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Family")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(b.e.y.n4b,filename = "figures/microb_Family_5percent_RA_Exp.YesNo_6.18.21.pdf", width=8, height=6, dpi=600)

b.e.y.n4c<-ggplot(b.fam.RA_e.y.n.m10, aes(x=exposed_y.n, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
  labs(title = "Microbial Family Relative Abundance by Exposed Condition", subtitle="Taxa > 10% Relative Abundance", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Family")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(b.e.y.n4c,filename = "figures/microb_Family_10percent_RA_Exp.YesNo_6.18.21.pdf", width=8, height=6, dpi=600)

# Genus
b.gen.RA_e.y.n.m1<-subset(b.gen.RA_e.y.n.m, c(Count)>(0.01)) ## DROP BACTERIAL GENERA that are less than 1% abundant!!!!!!1
b.gen.RA_e.y.n.m5<-subset(b.gen.RA_e.y.n.m, c(Count)>(0.05)) ## DROP BACTERIAL GENERA that are less than 5% abundant!!!!!!1
b.gen.RA_e.y.n.m10<-subset(b.gen.RA_e.y.n.m, c(Count)>(0.10)) ## DROP BACTERIAL GENERA that are less than 10% abundant!!!!!!1

# b.e.y.n5<-ggplot(b.gen.RA_e.y.n.m, aes(x=exposed_y.n, y=Count, fill=Genus))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
#   labs(title = "Microbial Genera Relative Abundance by Exposed Condition", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Genus")+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
# ggsave(b.e.y.n5,filename = "figures/microb_Genus_RA_Exp.YesNo_6.18.21.pdf", width=8, height=6, dpi=600)

b.e.y.n5a<-ggplot(b.gen.RA_e.y.n.m1, aes(x=exposed_y.n, y=Count, fill=Genus))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
  labs(title = "Microbial Genera Relative Abundance by Exposed Condition", subtitle="Taxa > 1% Relative Abundance", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(b.e.y.n5a,filename = "figures/microb_Genus_1percent_RA_Exp.YesNo_6.18.21.pdf", width=15, height=6, dpi=600)

b.e.y.n5b<-ggplot(b.gen.RA_e.y.n.m5, aes(x=exposed_y.n, y=Count, fill=Genus))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
  labs(title = "Microbial Genera Relative Abundance by Exposed Condition", subtitle="Taxa > 5% Relative Abundance", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(b.e.y.n5b,filename = "figures/microb_Genus_5percent_RA_Exp.YesNo_6.18.21.pdf", width=10, height=6, dpi=600)

b.e.y.n5c<-ggplot(b.gen.RA_e.y.n.m10, aes(x=exposed_y.n, y=Count, fill=Genus))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
  labs(title = "Microbial Genera Relative Abundance by Exposed Condition", subtitle="Taxa > 10% Relative Abundance", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(b.e.y.n5c,filename = "figures/microb_Genus_10percent_RA_Exp.YesNo_6.18.21.pdf", width=8, height=6, dpi=600)

# species
b.spec.RA_e.y.n.m1<-subset(b.spec.RA_e.y.n.m, c(Count)>(0.01)) ## DROP BACTERIAL Species that are less than 1% abundant!!!!!!1
b.spec.RA_e.y.n.m5<-subset(b.spec.RA_e.y.n.m, c(Count)>(0.05)) ## DROP BACTERIAL Species that are less than 5% abundant!!!!!!1
b.spec.RA_e.y.n.m10<-subset(b.spec.RA_e.y.n.m, c(Count)>(0.10)) ## DROP BACTERIAL Species that are less than 10% abundant!!!!!!1

# b.e.y.n5<-ggplot(b.spec.RA_e.y.n.m, aes(x=exposed_y.n, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
#   labs(title = "Microbial Species Relative Abundance by Exposed Condition", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Species")+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
# ggsave(b.e.y.n5,filename = "figures/microb_Species_RA_Exp.YesNo_6.18.21.pdf", width=8, height=6, dpi=600)

b.e.y.n5a<-ggplot(b.spec.RA_e.y.n.m1, aes(x=exposed_y.n, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
  labs(title = "Microbial Species Relative Abundance by Exposed Condition", subtitle="Taxa > 1% Relative Abundance", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Species")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(b.e.y.n5a,filename = "figures/microb_Species_1percent_RA_Exp.YesNo_6.18.21.pdf", width=20, height=6, dpi=600)

b.e.y.n5b<-ggplot(b.spec.RA_e.y.n.m5, aes(x=exposed_y.n, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
  labs(title = "Microbial Species Relative Abundance by Exposed Condition", subtitle="Taxa > 5% Relative Abundance", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Species")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(b.e.y.n5b,filename = "figures/microb_Species_5percent_RA_Exp.YesNo_6.18.21.pdf", width=8, height=6, dpi=600)

# b.exp5c<-ggplot(b.spec.RA_e.y.n.m10, aes(x=exposed_y.n, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Exposure Material", labels=c("Alt"="Alternaria", "Silica"="Silica", "Con"="Control"))+theme_classic()+
#   labs(title = "Microbial Species Relative Abundance by Exposure Group", subtitle="Taxa > 10% Relative Abundance", x="Exposure Material", y="Relative Abundance", fill="Species")+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
# ggsave(b.exp5c,filename = "figures/microb_Species_10percent_RA_Expo_6.18.21.pdf", width=8, height=6, dpi=600)



#### ITS2 Relative Abundance by Exposed: Yes or No ####
# head(f.dat.m)
# its2.asv_meta<-merge(f.dat.m, metadata, by="SampleID")

## phylum ....
head(its2.asv_meta)

# by phylum + exposure
its2_phy_e.y.n<- as.data.frame(dcast(its2.asv_meta,exposed_y.n~Phylum, value.var="Count", fun.aggregate=sum)) ###
head(its2_phy_e.y.n)
rownames(its2_phy_e.y.n)<-its2_phy_e.y.n$exposed_y.n
head(its2_phy_e.y.n)

f.phy.RA_e.y.n<-data.frame(decostand(its2_phy_e.y.n[,-c(1)], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(f.phy.RA_e.y.n)
head(f.phy.RA_e.y.n)

f.phy.RA_e.y.n$exposed_y.n<-rownames(f.phy.RA_e.y.n)
f.phy.RA_e.y.n.m<-melt(f.phy.RA_e.y.n)
head(f.phy.RA_e.y.n.m)
colnames(f.phy.RA_e.y.n.m)[which(names(f.phy.RA_e.y.n.m) == "variable")] <- "Phylum"
colnames(f.phy.RA_e.y.n.m)[which(names(f.phy.RA_e.y.n.m) == "value")] <- "Count"
head(f.phy.RA_e.y.n.m) ## relative abundance based on sum of counts by phyla in exposed/not exposed!

# by class + exposure
its2_cls_e.y.n <- as.data.frame(dcast(its2.asv_meta,exposed_y.n~Class, value.var="Count", fun.aggregate=sum)) ###
head(its2_cls_e.y.n) # counts by class + elevation
rownames(its2_cls_e.y.n)<-its2_cls_e.y.n$exposed_y.n

f.cls.RA_e.y.n<-data.frame(decostand(its2_cls_e.y.n[,-c(1)], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(f.cls.RA_e.y.n)
head(f.cls.RA_e.y.n)

f.cls.RA_e.y.n$exposed_y.n<-rownames(f.cls.RA_e.y.n)
f.cls.RA_e.y.n.m<-melt(f.cls.RA_e.y.n)
head(f.cls.RA_e.y.n.m)
colnames(f.cls.RA_e.y.n.m)[which(names(f.cls.RA_e.y.n.m) == "variable")] <- "Class"
colnames(f.cls.RA_e.y.n.m)[which(names(f.cls.RA_e.y.n.m) == "value")] <- "Count"
head(f.cls.RA_e.y.n.m) ## relative abundance based on sum of counts by class in exposed/not exposed!

# by order + exposure
its2_ord_e.y.n <- as.data.frame(dcast(its2.asv_meta,exposed_y.n~Order, value.var="Count", fun.aggregate=sum)) ###
head(its2_ord_e.y.n) # counts by Order + elevation
rownames(its2_ord_e.y.n)<-its2_ord_e.y.n$exposed_y.n

f.ord.RA_e.y.n<-data.frame(decostand(its2_ord_e.y.n[,-c(1)], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(f.ord.RA_e.y.n)
head(f.ord.RA_e.y.n)

f.ord.RA_e.y.n$exposed_y.n<-rownames(f.ord.RA_e.y.n)
f.ord.RA_e.y.n.m<-melt(f.ord.RA_e.y.n)
head(f.ord.RA_e.y.n.m)
colnames(f.ord.RA_e.y.n.m)[which(names(f.ord.RA_e.y.n.m) == "variable")] <- "Order"
colnames(f.ord.RA_e.y.n.m)[which(names(f.ord.RA_e.y.n.m) == "value")] <- "Count"
head(f.ord.RA_e.y.n.m) ## relative abundance based on sum of counts by order by exposed group!

# by Family + exposure
its2_fam_e.y.n <- as.data.frame(dcast(its2.asv_meta,exposed_y.n~Family, value.var="Count", fun.aggregate=sum)) ###
head(its2_fam_e.y.n) # counts by Family + elevation
rownames(its2_fam_e.y.n)<-its2_fam_e.y.n$exposed_y.n

f.fam.RA_e.y.n<-data.frame(decostand(its2_fam_e.y.n[,-c(1)], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(f.fam.RA_e.y.n)
head(f.fam.RA_e.y.n)

f.fam.RA_e.y.n$exposed_y.n<-rownames(f.fam.RA_e.y.n)
f.fam.RA_e.y.n.m<-melt(f.fam.RA_e.y.n)
head(f.fam.RA_e.y.n.m)
colnames(f.fam.RA_e.y.n.m)[which(names(f.fam.RA_e.y.n.m) == "variable")] <- "Family"
colnames(f.fam.RA_e.y.n.m)[which(names(f.fam.RA_e.y.n.m) == "value")] <- "Count"
head(f.fam.RA_e.y.n.m) ## relative abundance based on sum of counts by family by exposed group!

# by Genus + exposure
its2_gen_e.y.n <- as.data.frame(dcast(its2.asv_meta,exposed_y.n~Genus, value.var="Count", fun.aggregate=sum)) ###
head(its2_gen_e.y.n) # counts by genus + elevation
rownames(its2_gen_e.y.n)<-its2_gen_e.y.n$exposed_y.n

f.gen.RA_e.y.n<-data.frame(decostand(its2_gen_e.y.n[,-c(1)], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(f.gen.RA_e.y.n)
head(f.gen.RA_e.y.n)

f.gen.RA_e.y.n$exposed_y.n<-rownames(f.gen.RA_e.y.n)
f.gen.RA_e.y.n.m<-melt(f.gen.RA_e.y.n)
head(f.gen.RA_e.y.n.m)
colnames(f.gen.RA_e.y.n.m)[which(names(f.gen.RA_e.y.n.m) == "variable")] <- "Genus"
colnames(f.gen.RA_e.y.n.m)[which(names(f.gen.RA_e.y.n.m) == "value")] <- "Count"
head(f.gen.RA_e.y.n.m) ## relative abundance based on sum of counts by genus by exposed group!

# by Genus/Species + exposure
its2_spec_e.y.n <- as.data.frame(dcast(its2.asv_meta,exposed_y.n~Genus+Species, value.var="Count", fun.aggregate=sum)) ###
head(its2_spec_e.y.n) # counts by genus_species + elevation
rownames(its2_spec_e.y.n)<-its2_spec_e.y.n$exposed_y.n

f.spec.RA_e.y.n<-data.frame(decostand(its2_spec_e.y.n[,-c(1)], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(f.spec.RA_e.y.n)
head(f.spec.RA_e.y.n)

f.spec.RA_e.y.n$exposed_y.n<-rownames(f.spec.RA_e.y.n)
f.spec.RA_e.y.n.m<-melt(f.spec.RA_e.y.n)
head(f.spec.RA_e.y.n.m)
colnames(f.spec.RA_e.y.n.m)[which(names(f.spec.RA_e.y.n.m) == "variable")] <- "Genus_species"
colnames(f.spec.RA_e.y.n.m)[which(names(f.spec.RA_e.y.n.m) == "value")] <- "Count"
f.spec.RA_e.y.n.m$Genus_species<-gsub("_", " ", f.spec.RA_e.y.n.m$Genus_species) ## gsub is global sub (does not just remove first instance of pattern, but multiple)

head(f.spec.RA_e.y.n.m) ## relative abundance based on sum of counts by specus by exposed group!

#### Visualize ITS2 Relative Abundance by Exposed/ Not Exposed ####

# phylum
f.phy.RA_e.y.n.m1<-subset(f.phy.RA_e.y.n.m, c(Count)>(0.01)) ## DROP FUNGAL PHYLA that are less than 1% abundant!!!!!!1
f.phy.RA_e.y.n.m5<-subset(f.phy.RA_e.y.n.m, c(Count)>(0.05)) ## DROP FUNGAL PHYLA that are less than 5% abundant!!!!!!1
f.phy.RA_e.y.n.m10<-subset(f.phy.RA_e.y.n.m, c(Count)>(0.10)) ## DROP FUNGAL PHYLA that are less than 10% abundant!!!!!!1

f.e.y.n1<-ggplot(f.phy.RA_e.y.n.m, aes(x=exposed_y.n, y=Count, fill=Phylum))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
  labs(title = "Fungal Phylum Relative Abundance by Exposed Condition", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Phylum")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(f.e.y.n1,filename = "figures/fungi_phyla_RA_Exp.YesNo_6.18.21.pdf", width=8, height=6, dpi=600)

f.e.y.n1a<-ggplot(f.phy.RA_e.y.n.m1, aes(x=exposed_y.n, y=Count, fill=Phylum))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
  labs(title = "Fungal Phylum Relative Abundance by Exposed Condition", subtitle="Taxa > 1% Relative Abundance", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Phylum")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(f.e.y.n1a,filename = "figures/fungi_phyla_1percent_RA_Exp.YesNo_6.18.21.pdf", width=8, height=6, dpi=600)

f.e.y.n1b<-ggplot(f.phy.RA_e.y.n.m5, aes(x=exposed_y.n, y=Count, fill=Phylum))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
  labs(title = "Fungal Phylum Relative Abundance by Exposed Condition", subtitle="Taxa > 5% Relative Abundance", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Phylum")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(f.e.y.n1b,filename = "figures/fungi_phyla_5percent_RA_Exp.YesNo_6.18.21.pdf", width=8, height=6, dpi=600)

f.e.y.n1c<-ggplot(f.phy.RA_e.y.n.m10, aes(x=exposed_y.n, y=Count, fill=Phylum))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
  labs(title = "Fungal Phylum Relative Abundance by Exposed Condition", subtitle="Taxa > 10% Relative Abundance", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Phylum")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(f.e.y.n1c,filename = "figures/fungi_phyla_10percent_RA_Exp.YesNo_6.18.21.pdf", width=8, height=6, dpi=600)

#class
f.cls.RA_e.y.n.m1<-subset(f.cls.RA_e.y.n.m, c(Count)>(0.01)) ## DROP FUNGAL CLASSES that are less than 1% abundant!!!!!!1
f.cls.RA_e.y.n.m5<-subset(f.cls.RA_e.y.n.m, c(Count)>(0.05)) ## DROP FUNGAL CLASSES that are less than 5% abundant!!!!!!1
f.cls.RA_e.y.n.m10<-subset(f.cls.RA_e.y.n.m, c(Count)>(0.10)) ## DROP FUNGAL CLASSES that are less than 10% abundant!!!!!!1

f.e.y.n2<-ggplot(f.cls.RA_e.y.n.m, aes(x=exposed_y.n, y=Count, fill=Class))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
  labs(title = "Fungal Class Relative Abundance by Exposed Condition", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Class")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(f.e.y.n2,filename = "figures/fungi_class_RA_Exp.YesNo_6.18.21.pdf", width=8, height=6, dpi=600)

f.e.y.n2a<-ggplot(f.cls.RA_e.y.n.m1, aes(x=exposed_y.n, y=Count, fill=Class))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
  labs(title = "Fungal Class Relative Abundance by Exposed Condition", subtitle="Taxa > 1% Relative Abundance", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Class")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(f.e.y.n2a,filename = "figures/fungi_class_1percent_RA_Exp.YesNo_6.18.21.pdf", width=8, height=6, dpi=600)

f.e.y.n2b<-ggplot(f.cls.RA_e.y.n.m5, aes(x=exposed_y.n, y=Count, fill=Class))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
  labs(title = "Fungal Class Relative Abundance by Exposed Condition", subtitle="Taxa > 5% Relative Abundance", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Class")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(f.e.y.n2b,filename = "figures/fungi_class_5percent_RA_Exp.YesNo_6.18.21.pdf", width=8, height=6, dpi=600)

f.e.y.n2c<-ggplot(f.cls.RA_e.y.n.m10, aes(x=exposed_y.n, y=Count, fill=Class))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
  labs(title = "Fungal Class Relative Abundance by Exposed Condition", subtitle="Taxa > 10% Relative Abundance", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Class")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(f.e.y.n2c,filename = "figures/fungi_class_10percent_RA_Exp.YesNo_6.18.21.pdf", width=8, height=6, dpi=600)

#order
f.ord.RA_e.y.n.m1<-subset(f.ord.RA_e.y.n.m, c(Count)>(0.01)) ## DROP FUNGAL ORDERS that are less than 1% abundant!!!!!!1
f.ord.RA_e.y.n.m5<-subset(f.ord.RA_e.y.n.m, c(Count)>(0.05)) ## DROP FUNGAL ORDERS that are less than 5% abundant!!!!!!1
f.ord.RA_e.y.n.m10<-subset(f.ord.RA_e.y.n.m, c(Count)>(0.10)) ## DROP FUNGAL ORDERS that are less than 10% abundant!!!!!!1

f.e.y.n3<-ggplot(f.ord.RA_e.y.n.m, aes(x=exposed_y.n, y=Count, fill=Order))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
  labs(title = "Fungal Order Relative Abundance by Exposed Condition", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Order")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(f.e.y.n3,filename = "figures/fungi_order_RA_Exp.YesNo_6.18.21.pdf", width=8, height=6, dpi=600)

f.e.y.n3a<-ggplot(f.ord.RA_e.y.n.m1, aes(x=exposed_y.n, y=Count, fill=Order))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
  labs(title = "Fungal Order Relative Abundance by Exposed Condition", subtitle="Taxa > 1% Relative Abundance", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Order")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(f.e.y.n3a,filename = "figures/fungi_order_1percent_RA_Exp.YesNo_6.18.21.pdf", width=8, height=6, dpi=600)

f.e.y.n3b<-ggplot(f.ord.RA_e.y.n.m5, aes(x=exposed_y.n, y=Count, fill=Order))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
  labs(title = "Fungal Order Relative Abundance by Exposed Condition", subtitle="Taxa > 5% Relative Abundance", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Order")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(f.e.y.n3b,filename = "figures/fungi_order_5percent_RA_Exp.YesNo_6.18.21.pdf", width=8, height=6, dpi=600)

f.e.y.n3c<-ggplot(f.ord.RA_e.y.n.m10, aes(x=exposed_y.n, y=Count, fill=Order))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
  labs(title = "Fungal Order Relative Abundance by Exposed Condition", subtitle="Taxa > 10% Relative Abundance", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Order")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(f.e.y.n3c,filename = "figures/fungi_order_10percent_RA_Exp.YesNo_6.18.21.pdf", width=8, height=6, dpi=600)

# family
f.fam.RA_e.y.n.m1<-subset(f.fam.RA_e.y.n.m, c(Count)>(0.01)) ## DROP FUNGAL FamilyS that are less than 1% abundant!!!!!!1
f.fam.RA_e.y.n.m5<-subset(f.fam.RA_e.y.n.m, c(Count)>(0.05)) ## DROP FUNGAL FamilyS that are less than 5% abundant!!!!!!1
f.fam.RA_e.y.n.m10<-subset(f.fam.RA_e.y.n.m, c(Count)>(0.10)) ## DROP FUNGAL FamilyS that are less than 10% abundant!!!!!!1

f.e.y.n4<-ggplot(f.fam.RA_e.y.n.m, aes(x=exposed_y.n, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
  labs(title = "Fungal Family Relative Abundance by Exposed Condition", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Family")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(f.e.y.n4,filename = "figures/fungi_Family_RA_Exp.YesNo_6.18.21.pdf", width=8, height=6, dpi=600)

f.e.y.n4a<-ggplot(f.fam.RA_e.y.n.m1, aes(x=exposed_y.n, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
  labs(title = "Fungal Family Relative Abundance by Exposed Condition", subtitle="Taxa > 1% Relative Abundance", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Family")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(f.e.y.n4a,filename = "figures/fungi_Family_1percent_RA_Exp.YesNo_6.18.21.pdf", width=8, height=6, dpi=600)

f.e.y.n4b<-ggplot(f.fam.RA_e.y.n.m5, aes(x=exposed_y.n, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
  labs(title = "Fungal Family Relative Abundance by Exposed Condition", subtitle="Taxa > 5% Relative Abundance", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Family")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(f.e.y.n4b,filename = "figures/fungi_Family_5percent_RA_Exp.YesNo_6.18.21.pdf", width=8, height=6, dpi=600)

f.e.y.n4c<-ggplot(f.fam.RA_e.y.n.m10, aes(x=exposed_y.n, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
  labs(title = "Fungal Family Relative Abundance by Exposed Condition", subtitle="Taxa > 10% Relative Abundance", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Family")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(f.e.y.n4c,filename = "figures/fungi_Family_10percent_RA_Exp.YesNo_6.18.21.pdf", width=8, height=6, dpi=600)

# Genus
f.gen.RA_e.y.n.m1<-subset(f.gen.RA_e.y.n.m, c(Count)>(0.01)) ## DROP FUNGAL GENERA that are less than 1% abundant!!!!!!1
f.gen.RA_e.y.n.m5<-subset(f.gen.RA_e.y.n.m, c(Count)>(0.05)) ## DROP FUNGAL GENERA that are less than 5% abundant!!!!!!1
f.gen.RA_e.y.n.m10<-subset(f.gen.RA_e.y.n.m, c(Count)>(0.10)) ## DROP FUNGAL GENERA that are less than 10% abundant!!!!!!1

f.e.y.n5<-ggplot(f.gen.RA_e.y.n.m, aes(x=exposed_y.n, y=Count, fill=Genus))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
  labs(title = "Fungal Genera Relative Abundance by Exposed Condition", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(f.e.y.n5,filename = "figures/fungi_Genus_RA_Exp.YesNo_6.18.21.pdf", width=8, height=6, dpi=600)

f.e.y.n5a<-ggplot(f.gen.RA_e.y.n.m1, aes(x=exposed_y.n, y=Count, fill=Genus))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
  labs(title = "Fungal Genera Relative Abundance by Exposed Condition", subtitle="Taxa > 1% Relative Abundance", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(f.e.y.n5a,filename = "figures/fungi_Genus_1percent_RA_Exp.YesNo_6.18.21.pdf", width=8, height=6, dpi=600)

f.e.y.n5b<-ggplot(f.gen.RA_e.y.n.m5, aes(x=exposed_y.n, y=Count, fill=Genus))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
  labs(title = "Fungal Genera Relative Abundance by Exposed Condition", subtitle="Taxa > 5% Relative Abundance", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(f.e.y.n5b,filename = "figures/fungi_Genus_5percent_RA_Exp.YesNo_6.18.21.pdf", width=8, height=6, dpi=600)

f.e.y.n5c<-ggplot(f.gen.RA_e.y.n.m10, aes(x=exposed_y.n, y=Count, fill=Genus))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
  labs(title = "Fungal Genera Relative Abundance by Exposed Condition", subtitle="Taxa > 10% Relative Abundance", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(f.e.y.n5c,filename = "figures/fungi_Genus_10percent_RA_Exp.YesNo_6.18.21.pdf", width=8, height=6, dpi=600)

# species
f.spec.RA_e.y.n.m1<-subset(f.spec.RA_e.y.n.m, c(Count)>(0.01)) ## DROP FUNGAL Species that are less than 1% abundant!!!!!!1
f.spec.RA_e.y.n.m5<-subset(f.spec.RA_e.y.n.m, c(Count)>(0.05)) ## DROP FUNGAL Species that are less than 5% abundant!!!!!!1
f.spec.RA_e.y.n.m10<-subset(f.spec.RA_e.y.n.m, c(Count)>(0.10)) ## DROP FUNGAL Species that are less than 10% abundant!!!!!!1

f.e.y.n5<-ggplot(f.spec.RA_e.y.n.m, aes(x=exposed_y.n, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
  labs(title = "Fungal Species Relative Abundance by Exposed Condition", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Species")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(f.e.y.n5,filename = "figures/fungi_Species_RA_Exp.YesNo_6.18.21.pdf", width=8, height=6, dpi=600)

f.e.y.n5a<-ggplot(f.spec.RA_e.y.n.m1, aes(x=exposed_y.n, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
  labs(title = "Fungal Species Relative Abundance by Exposed Condition", subtitle="Taxa > 1% Relative Abundance", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Species")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(f.e.y.n5a,filename = "figures/fungi_Species_1percent_RA_Exp.YesNo_6.18.21.pdf", width=8, height=6, dpi=600)

f.e.y.n5b<-ggplot(f.spec.RA_e.y.n.m5, aes(x=exposed_y.n, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
  labs(title = "Fungal Species Relative Abundance by Exposed Condition", subtitle="Taxa > 5% Relative Abundance", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Species")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(f.e.y.n5b,filename = "figures/fungi_Species_5percent_RA_Exp.YesNo_6.18.21.pdf", width=8, height=6, dpi=600)

f.e.y.n5c<-ggplot(f.spec.RA_e.y.n.m10, aes(x=exposed_y.n, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete(name ="Not Exposed vs. Exposed", labels=c("N"="No", "Y"="Yes"))+theme_classic()+
  labs(title = "Fungal Species Relative Abundance by Exposure Condition", subtitle="Taxa > 10% Relative Abundance", x="Not Exposed vs. Exposed", y="Relative Abundance", fill="Species")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13), legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=12))
ggsave(f.e.y.n5c,filename = "figures/fungi_Species_10percent_RA_Exp.YesNo_6.18.21.pdf", width=8, height=6, dpi=600)











#### Cumulative sum scaling ####
## CSS for beta diversity!
bac.ASV_counts_1<-newMRexperiment(bac.ASV_counts) # samples - columns, OTUs - rows ; newMRexperiment creates MRexperiment object
ASV_count_css<-as.data.frame(cumNormMat(bac.ASV_counts_1, p = cumNormStatFast(bac.ASV_counts_1), sl = 1000)) # Cumulative sum scaling factors: Calculates each column's quantile and calculates the sum up to and including that quantile.
### cumNormMat - Returns a matrix normalized by scaling counts up to and including the pth quantile.
## cumNormStatFast --> Calculates the percentile for which to sum counts up to and scale by.
head(ASV_count_css)

## create a distance matrix of biological data for hierarchical clustering
dist_asvs = as.matrix((vegdist(ASV_count_css, "bray", na.rm = TRUE))) #distance matrix using Bray's dissimilarity index for trait distribution (traits of interest only)

dist_asvs<-na.omit(dist_asvs) ## getting rid of NAs, different than making Nas = 0 **

plot(hclust(vegdist(ASV_count_css,method="bray")))     #####clustering with all the taxa!
plot(hclust(dist_asvs))
#### Rarefaction ####
# RAREFACTION in vegan: ROWS need to be SITES/samples; COLUMNS are SPECIES (OTUs, ASVs)
# see "Numerical Ecology with R", page 13-14
min<-min(rowSums(ASV_table))
rar_ASV_table<-rrarefy(ASV_table,min) ## THIS IS RETURNING A DF OF 0s, FIGURE OUT WTF IS WRONG HERE (4.9.2020)
head(rar_ASV_table)
#### Alpha Diversity + Species Richness ####
# ## H <- diversity() -- Shannon entropy function
## exp(H) -- Shannon diversity (exp of entropy)
# use rarefied data for alpha div; raw counts for species richness
Shan_ent.16s<-vegan::diversity(ASV_table, index="shannon") # Shannon entropy
Shan_div.16s<- exp(Shan_ent.16s) # Shannon Diversity aka Hill number 1
Shannon.raw.16s<-data.frame(Shan_ent.16s,Shan_div.16s)
class(Shannon.raw.16s)
Shannon.raw.16s$SampleID<-rownames(Shannon.raw.16s)

Shan_ent.rar.16S<-vegan::diversity(rar_ASV_table, index = "shannon") #### Shannon ENTROPY (alpha diversity)
Shan_div.rar.16s<- exp(Shan_ent.rar.16S) # Shannon Diversity aka Hill number 1
Shannon.rar.16s<-data.frame(Shan_ent.rar.16S, Shan_div.rar.16s)
class(Shannon.rar.16s)
Shannon.rar.16s

Shannon.rar.16s$SampleID<-rownames(Shannon.rar.16s)

# save as R objects to you can load them into fungi R script and combine with fungi shannon info!
#saveRDS(Shannon.raw.16s, file = "16S_Sierra_Shannon_from.rawdata_9.13.2020_Robject.rds", ascii = FALSE, version = NULL,
#        compress = TRUE, refhook = NULL)

#saveRDS(Shannon.rar.16s, file = "16S_Sierra_Shannon_from.rarefieddata_9.13.2020_Robject.rds", ascii = FALSE, version = NULL,
#        compress = TRUE, refhook = NULL)

#write.csv(Shannon.raw.16s,"16S_Sierra_Shannon_from.rawdata_9.13.2020.csv",row.names=TRUE)
#write.csv(Shannon.rar.16s,"16S_Sierra_Shannon_from.rarefieddata_9.13.2020.csv",row.names=TRUE)


S<-specnumber(ASV_table) # finds # of species per sample using RAW count data; if MARGIN = 2 it finds frequencies of species
S # tells you number of ASVs
S.freq<-specnumber(ASV_table, MARGIN = 2) # # finds how many times each ASV appeared across samples (frequency)
S.freq


#### Merge Alpha div w/ count/taxa data + metadata ####

alpha.div.metadat <- merge(Shannon.raw.16s,metadata, by="SampleID")
head(alpha.div.metadat)
class(alpha.div.metadat)
#write.csv(alpha.div.metadat,"SierraDust_16S_Metadata_AlphaDiversity_9.13.2020.csv",row.names=FALSE)

alpha.div.metadat$Type2 <- factor(alpha.div.metadat$Type, levels = c("air","soil","water","control")) ## reeorder month column so that it will be listed in chronological order in x axis

#### Alpha Diversity Visualization ####
head(alpha.div.metadat)

# shannon entropy across sample type
alpha.ent.raw.elev<-ggplot(alpha.div.metadat, aes(x=Type2, y=Shan_ent.16s, fill=Type2)) +geom_boxplot(color="black")+scale_x_discrete(name ="Sample Type", labels=c("air"="Aeolian", "soil"="Soil", "water"="Water","control"="Control"))+theme_bw()+scale_fill_manual(values = saturation(ewf, 0.6), name ="Sample Type", labels=c("air"="Aeolian", "soil"="Soil", "water"="Water", "control"="Control"))+theme_classic()+
  labs(title = "Bacterial Shannon Entropy by Elevation", x="Elevation", y="Shannon Entropy", fill="Sample Type")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(alpha.ent.raw.elev,filename = "16S_shannon_entropy_raw_by_elevation_6.16.2021.pdf", width=8, height=6, dpi=600)

# shannon diversity across sample type
alpha.div.raw.elev<-ggplot(alpha.div.metadat, aes(x=Type2, y=Shan_div.16s, fill=Type2)) +geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(ewf, 0.6), name ="Sample Type", labels=c("air"="Aeolian", "soil"="Soil", "water"="Water", "control"="Control"))+theme_classic()+
  labs(title = "Bacterial Shannon Diversity by Elevation", x="Elevation", y="Shannon Diversity", fill="Sample Type")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(alpha.div.raw.elev,filename = "16S_shannon_diversity_raw_by_elevation_6.16.2021.pdf", width=8, height=6, dpi=600)

# shannon diversity across elevation + year (no box fill, color outlines only)
#alpha.div.raw.elev.yr.2<-ggplot(alpha.div.metadat, aes(x=factor(Elevation), y=Shan_div.16s, col=factor(Year)))+geom_boxplot()+scale_x_discrete()+theme_bw()+scale_colour_manual(values = saturation(wes1, 1))+theme_classic()+
#  labs(title = "Bacterial Shannon Diversity by Elevation and Year", x="Elevation", y="Shannon Diversity", col="Year")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
#ggsave(alpha.div.raw.elev.yr.2,filename = "bacteria_shannon_diversity_raw_by_elevation.year.2_9.13.2020.pdf", width=8, height=6, dpi=600)

dev.off()


#### Relativize OTU table, then make Bray-Curtis dissimilarity distance matrix ####
# use this for PCoA -- this bypasses the need for using a Cailliez or Lingoes correction w/ PCoA (need to get rid of negative eigenvalues)
head(ASV_table) # raw counts
# t.otu_counts$SampleID<-NULL
# t.otu_counts<-data.frame(t.otu_counts)
ASV_table<-ASV_table[ which(rowSums(ASV_table) > 0),] # drop 0s before distance matrix is made!!!
rowSums(ASV_table)

ASV_table_RA<-data.frame(decostand(ASV_table, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples

rowSums(ASV_table_RA)
rownames(ASV_table_RA)
sample_names<-data.frame(SampleID=rownames(ASV_table_RA))

BrayC_ASV<-vegdist(ASV_table_RA,method="bray")     ####### Bray-Curtis dissimliarity distance matrix with raw data!!!

#t.otu_counts$SampleID<-rownames(t.otu_counts)

# bray.otu.counts<-vegdist(t.otu_counts,method="bray")     ####### distance matrix with raw data!!!
# # **** in vegan: ROWS need to be SITES/samples; COLUMNS are SPECIES (OTUs, ASVs) -- vegdist!
# sqrt.bray.otu<-sqrt(bray.otu.counts) # square root transformation of Bray-curtis dissimilarity matrix -- eliminates negative eigenvalues, making it "true" Euclidean distances for the ordination
# class(bray.otu.counts) # needs to be class dist for pcoa()

#### Saving OTU counts (Bray-Curtis dissimilarity) to compare to fungal diversity via mantel test ####

#saveRDS(bray.otu.counts, file = "16S_Sierra_braycurtis.dissimilarity_wArchaea_9.13.20_Robject.rds", ascii = FALSE, version = NULL,
#        compress = TRUE, refhook = NULL)

#### Standardizing Environmental Data (metadata) by creating Euclidean distance matrices (to compare with diversity, etc) ####
head(metadata)
rownames(ASV_table)
metadata=metadata[rownames(ASV_table),] ## reorder metadata to have same rows as original OTU table
# 1-18-21: can't standardize these types of env variables yet, need more metadata

meta_chem1<-subset(metadata, select=c(SampleID, Cu, Fe, Mg, Mn, Ni, P, S, Zn)) # subset only part of the metadata we need
rownames(meta_chem1)<-meta_chem1$SampleID
meta_chem1$SampleID<-NULL
head(meta_chem1)

log.chem1<-dist(log(meta_chem1))

std.chem1<-vegdist(scale(meta_chem1), "euclid")

meta_chem2<-subset(metadata, select=c(SampleID, NdEp, X87Sr86Sr)) # subset only part of the metadata we need
rownames(meta_chem2)<-meta_chem2$SampleID
meta_chem2$SampleID<-NULL
head(meta_chem2)
meta_chem2<-sapply(meta_chem2, as.numeric)
head(meta_chem2)
rownames(meta_chem2)<-metadata$SampleID
head(meta_chem2)

log.chem2<-dist(log(meta_chem2))

std.chem2<-vegdist(scale(meta_chem2), "euclid", na.rm=TRUE)


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


