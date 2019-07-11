rm(list=ls()) #to clean R environment  
ctrll# to clean console in R
R.Version()
setwd("D:/Rcuracao") #set working directory

################################################################install packages first time########################################################
#devtools
install.packages("devtools")
library(devtools)
source("https://bioconductor.org/biocLite.R")#use this command if "install.packages" command doesn't work 
biocLite("Biostrings") #Memory efficient string containers, to manipulate large biological sequences or sets of sequences
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")

#Since microbiome does not look for all the pacakges it depends on, you need to install the following packages first

source("https://bioconductor.org/biocLite.R")
bioclite("ggplot2") #to declare the input data frame for a graphic
#to read, write, plot, and manipulate phylogenetic trees
install.packages("picante",repos="http://R-Forge.R-project.org") 
install.packages("plyr")#to  Split, apply and combine Data 
install.packages("vegan") # to run all common ordination methods (CCA, PCA, PCOA., etc)
install.packages("tensorA")
install.packages("bayesm")
install.packages("energy")
install.packages("robustbase")
install.packages("plotly")
install.packages("DT")   #to display R data object into tables on html pages
install.packages("RColorBrewer")  # to choose sensible colour schemes for figures in R
install.packages("reshape2")  # to transform data between wide and long formats.
install.packages("scales")  #to provide methods for automatically determining breaks and labels for axes and legends
install.packages("data.table") #to manipulate data operations such as subset, group, update, join etc
install.packages("DESeq2")  #to use negative binomial for differential OTU analysis
install.packages("dplyr")   #to transform and summarize tabular data with rows and columns
source("https://raw.githubusercontent.com/joey711/phyloseq/master/inst/scripts/installer.R",
       local = TRUE)
install.packages("ggthemes")
install.packages("ggdendro")
install.packages("picante")
install.packages("maptree")
install.packages("rmdformats")
install.packages("dynamicTreeCut")
install.packages("fastcluster")
install.packages("Rcpp")# better to install from right box install button.
source("https://bioconductor.org/biocLite.R")
biocLite("metagenomeSeq")
source("https://bioconductor.org/biocLite.R")
biocLite("GO.db")
source("http://www.bioconductor.org/biocLite.R")
biocLite(c("ade4", "fastcluster", "df2json", "rjson", "gplots", "devtools", "ggplot2","MASS","minet","mixOmics", "plyr","qvalue","reshape2","RPA","svDialogs","vegan","WGCNA"))
library(BiocInstaller)
source("http://www.bioconductor.org/biocLite.R")
useDevel()
biocLite("microbiome")
# Test installation by loading the microbiome package in R
library("microbiome") #to analyse tasks in microbiome profiling studies
#load required packages: 
library(ggplot2) 
library(ape)
library(plyr)
library(vegan)
library(Biostrings)
library(RColorBrewer)
library(reshape2)
library(scales)
library(data.table)
library(DESeq2)
library(microbiome)
library(dplyr)
library(phyloseq)
library(DT)
library(Rcpp)
library(picante)
library(pheatmap)
###################################################Getting Started#########################################################################
theme_set(theme_bw()) # Set plotting theme
#create a phyloseq file 
pseq1 <- read_phyloseq(otu.file = "DG_final.biom", taxonomy.file = NULL, metadata.file = "MappingFile_DG.csv", type = "biom")
#Read three file
library(ape)
treefile_p1 <- read.tree("DG.tre")
#Merge into phyloseq object
ps1 <-merge_phyloseq(pseq1,treefile_p1)
ps1 <- prune_taxa(taxa_sums(ps1) > 0, ps1) # remove row with total value of 0
psHeatMap <- ps1 # phyloseq object to generate heatmap
###
otu_table_ps1 <- as.data.frame(ps1@otu_table)
write.csv(otu_table_ps1, file = "otu_table_ps1.csv", fileEncoding = "UTF-16LE")
#the table is interactive 
datatable(tax_table(ps1))
rank_names(ps1)#check rank of samples
#Taxonomic table
tax_table(ps1)
#with this command, you can investigate the number of reads for each sample:
sort(sample_sums(ps1))

##investigate read of each sample
reads_per_OTU <- taxa_sums(ps1)
print(sum(reads_per_OTU))

#####################################################Phylogenetic Diversity#####################################################################
#Calculate Phylogenetic Distance (PD) of the dataset
my_colors <- c("#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD", "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "steelblue" )
otu_table_ps1 <- as.data.frame(ps1@otu_table)
metadata_table_ps1  <- as.data.frame(ps1@sam_data)
df.pd <- pd(t(otu_table_ps1), treefile_p1,include.root=F) 
datatable(df.pd)
ps1.metadata <- as(sample_data(ps1), "data.frame")
head(ps1.metadata)
#plot Phylogenetic Diversity of all samples
ps1.metadata$Phyogenetic_diversity <- df.pd$PD
colnames(ps1.metadata)
Reads <- sort(sample_sums(ps1)) 
#plot Phylogenetic Diversity of all samples
metadata_table_ps1$Phyogenetic_diversity <- df.pd$PD 
plot.pd1 <- ggplot(metadata_table_ps1, aes(Sample_types, Phyogenetic_diversity)) + geom_boxplot(aes(fill = Depth)) + ggtitle("Phylogenetic Diversity of All Samples") + theme(axis.text.x = element_text(size=14, angle = 90)) 
print(plot.pd1)

##test of significance based on PD ##
##all samples
ps1.adiv <- estimate_richness(ps1, measures = c("Chao1", "Shannon", "Observed", "InvSimpson"))
ps1.metadata <- as(sample_data(ps1), "data.frame")
head(ps1.metadata)
ps1.metadata$Phyogenetic_diversity <- df.pd$PD####add PD value
ps1.metadata$Observed <- ps1.adiv$Observed 
ps1.metadata$Shannon <- ps1.adiv$Shannon
ps1.metadata$InvSimpson <- ps1.adiv$InvSimpson
colnames(ps1.metadata)
write.table(ps1.metadata, file="ps1_metadata.csv",sep=",",row.names=F)# import metadata with PD value

##Some of the taxonomy levels may include "<empty>" entries. If you want to remove them for plotting then run the following code.  
ps1.com <- ps1 # create a new pseq object for generating boxplot
tax.mat <- tax_table(ps1.com)
tax.df <- as.data.frame(tax.mat)
tax.df[tax.df == "__*~"] <- "__"
tax.df[tax.df == "k__<empty>"] <- "k__"
tax.df[tax.df == "p__<empty>"] <- "p__"
tax.df[tax.df == "c__<empty>"] <- "c__"
tax.df[tax.df == "o__<empty>"] <- "o__"
tax.df[tax.df == "f__<empty>"] <- "f__"
tax.df[tax.df == "g__<empty>"] <- "g__"
tax_table(ps1.com) <- tax_table(as.matrix(tax.df))

######Barplot composition based on Phylum level###########
# We need to set Palette
taxic <- as.data.frame(ps1.com@tax_table)  
colourCount = length(unique(taxic$Phylum))  #define number of variable colors based on number of Phylum
getPalette = colorRampPalette(brewer.pal(9, "Set1"))  # change the palette as well as the number of colors will change according to palette.
otu.df <- as.data.frame(otu_table(ps1.com))
head(otu.df)
taxic$OTU <- row.names.data.frame(otu.df)
colnames(taxic)
library(knitr)
head(kable(taxic))  # check the table.
taxmat <- as.matrix(taxic)  # convert it into a matrix.
new.tax <- tax_table(taxmat)  # convert into phyloseq compatible file.
tax_table(ps1.com) <- new.tax  # incroporate into phyloseq Object
guide_italics <- guides(fill = guide_legend(label.theme = element_text(size = 15,face = "italic", colour = "Black", angle = 0)))
ps1.com@phy_tree <- NULL # remove tree 
ps1.com.phyl <- aggregate_taxa(ps1.com, "Phylum") #aggregate to Phylum level
ps1.com.phyl.rel <- microbiome::transform(ps1.com.phyl, "compositional")# change counts to relative abudance
otu_table.phyl.rel <- otu_table(ps1.com.phyl.rel)#generate OTU table 
write.csv(otu_table.phyl.rel, file = "otu_table_phylum_rel.csv", fileEncoding = "UTF-16LE")# check result in excel 
ps1.com.phyl.rel.fr <- filter_taxa(ps1.com.phyl.rel, function(x) mean(x) >= 0.0025, TRUE) #select Phyla that have an average 0.25% in all samples 
otu_table.phyl.rel0.25<- otu_table(ps1.com.phyl.rel.fr)
write.csv(otu_table.phyl.rel0.25, file = "otu_table_phylum0.25%.csv", fileEncoding = "UTF-16LE")
#generate box plot
plot.composition.relAbun.phylum <- plot_composition(ps1.com.phyl.rel.fr) + theme(legend.position = "bottom") + 
  scale_fill_manual(values = getPalette(colourCount)) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + 
  ggtitle("Relative abundance based on Phylum level") + guide_italics 
plot.composition.relAbun.phylum

#generate original box plot phylum level
plot.composition.relAbun.phylum <- plot_composition(ps1.com.phyl.rel) + theme(legend.position = "bottom") + 
  scale_fill_manual(values = getPalette(colourCount)) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + 
  ggtitle("Relative abundance based on Phylum level") + guide_italics 
plot.composition.relAbun.phylum

########Beta diversity######
#Bray Curtis dissimilarity, transformed based on Hellinger (square root of relative abundance, but for this command values were generated from counts data)
ps1.hel <- transform(ps1, "hellinger") 
otu_table_ps1.hel <- as.data.frame(ps1.hel@otu_table)#otu_table
View(otu_table_ps1.hel)#to check how transformation to relative abundance looks like 
write.csv(otu_table_ps1.hel, file = "otu_table_ps1.hel.csv", row.names=TRUE)# save in csv
ordu.bray.hel <- ordinate(ps1.hel, "PCoA", "bray", weighted=F)
bray.hel <- plot_ordination(ps1.hel, ordu.bray.hel, color="Sample_types", shape="Depth", label ="DepthMeter")
bray_hel <- bray.hel + scale_fill_manual(values = c("#CBD588", "#5F7FC7", "orange","#DA5724", "#508578")) + ggtitle("A. Bray-Curtis dissimilarity of all samples") + geom_point(size = 3) + theme_bw()
print(bray_hel)

#Bray Curtis new all samples
ps1.hel <- transform(ps1, "hellinger") 
otu_table_ps1.hel <- as.data.frame(ps1.hel@otu_table)#otu_table
View(otu_table_ps1.hel)#to check how transformation to relative abundance looks like 
write.csv(otu_table_ps1.hel, file = "otu_table_ps1.hel.csv", row.names=TRUE)# save in csv
ordu.bray.hel <- ordinate(ps1.hel, "PCoA", "bray", weighted=F)
bray.hel <- plot_ordination(ps1.hel, ordu.bray.hel, color="Sample_types", shape="Depth")
bray_hel <- bray.hel + scale_fill_manual(values = c("#CBD588", "#5F7FC7", "orange","#DA5724", "#508578")) + geom_point(size = 6) + theme_bw()
print(bray_hel)


#Bray Curtis per samples
#X.muta
psXmuta <-subset_samples(ps1, Sample_types == "X.muta") 
psXmuta = prune_taxa(taxa_sums(psXmuta) > 0, psXmuta)
print(psXmuta)
psXmuta.hel <- transform(psXmuta, "hellinger")
otu_table_psXmuta.hel <- as.data.frame(psXmuta.hel@otu_table)
write.csv(otu_table_psXmuta.hel, file = "otu_table_psXmuta.hel.csv", row.names=TRUE)
ordu.bray.Xmuta.hel = ordinate(psXmuta.hel, "PCoA", "bray", weighted=F)
bray.Xmuta.hel <- plot_ordination(psXmuta.hel, ordu.bray.Xmuta.hel,color="Depth", shape="Depth")
bray_Xmuta.hel <- bray.Xmuta.hel + scale_colour_manual(values = c("royalblue4", "royalblue1", "steelblue1")) + geom_point(size = 6) + theme_bw()
print(bray_Xmuta.hel)

#A.sventres_new
psAsventres <- subset_samples(ps1, Sample_types == "A.sventres") 
psAsventres <- prune_taxa(taxa_sums(psAsventres) > 0, psAsventres)
print(psAsventres)
psAsventres.hel <- transform(psAsventres, "hellinger")
otu_table_psAsventres.hel <- as.data.frame(psAsventres.hel@otu_table)
write.csv(otu_table_psAsventres.hel, file = "otu_table_psAsventres.hel.csv", row.names=TRUE)
ordu.bray.Asventres.hel <- ordinate(psAsventres.hel, "PCoA", "bray", weighted=F)
bray.Asventres.hel <- plot_ordination(psAsventres.hel, ordu.bray.Asventres.hel, color="Depth", shape="Depth")
bray_Asventres.hel <- bray.Asventres.hel + scale_colour_manual(values = c("red3", "red")) + geom_point(size = 6) + theme_bw()
print(bray_Asventres.hel)

#SeaWater
psSeawater <-subset_samples(ps1, Sample_types == "Seawater") 
psSeawater  <- prune_taxa(taxa_sums(psSeawater ) > 0, psSeawater)
print(psSeawater)
psSeawater.hel <- transform(psSeawater, "hellinger")
otu_table_psSeawater.hel <- as.data.frame(psSeawater.hel@otu_table)
write.csv(otu_table_psSeawater.hel, file = "otu_table_psSeawater.hel.csv", row.names=TRUE)
ordu.bray.Seawater.hel <- ordinate(psSeawater.hel, "PCoA", "bray", weighted=F)
bray.Seawater.hel <- plot_ordination(psSeawater.hel, ordu.bray.Seawater.hel, color="Depth", shape="Depth")
bray_Seawater.hel <- bray.Seawater.hel + scale_colour_manual(values = c("darkgreen", "green3", "green")) + geom_point(size = 6) + theme_bw()
print(bray_Seawater.hel)

####################Statistical analysis#####################
####Kruskall Wallis based on PD for all samples 
ps1.metadata$Phyogenetic_diversity <- df.pd$PD
#based on parameter Depth
kruskal.PD_Depth <- kruskal.test(ps1.metadata$Phyogenetic_diversity ~ ps1.metadata$Depth)
print(kruskal.PD_Depth)
#based on parameter Samples
kruskal.PD_Samples <- kruskal.test(ps1.metadata$Phyogenetic_diversity ~ ps1.metadata$Sample_types)
print(kruskal.PD_Samples)
#Kruskall Wallis subset per sample identity
###X. muta
psXmuta <- subset_samples(ps1, Sample_types == "X.muta") #subset X. muta
psXmuta <- prune_taxa(taxa_sums(psXmuta) > 0, psXmuta)
otu_table_psXmuta <- as.data.frame(psXmuta@otu_table) 
metadata_table_Xmuta  <- as.data.frame(psXmuta@sam_data)
df.pdXmuta <- pd(t(otu_table_psXmuta), treefile_p1,include.root=F) 
datatable(df.pdXmuta)
psXmuta.metadata <- as(sample_data(psXmuta), "data.frame")
psXmuta.metadata$Phyogenetic_diversity <- df.pdXmuta$PD
#based on parameter Depth
kruskal.Xmuta.PD_Depth <- kruskal.test(psXmuta.metadata$Phyogenetic_diversity ~ psXmuta.metadata$Depth)
print(kruskal.Xmuta.PD_Depth)
### A. sventres
psAsventres <- subset_samples(ps1, Sample_types == "A.sventres") #subset A.sventres
psAsventres <- prune_taxa(taxa_sums(psAsventres) > 0, psAsventres)
otu_table_psAsventres <- as.data.frame(psAsventres@otu_table) 
metadata_table_psAsventres  <- as.data.frame(psAsventres@sam_data)
df.pdAsventres <- pd(t(otu_table_psAsventres), treefile_p1,include.root=F) 
datatable(df.pdAsventres)
psAsventres.metadata <- as(sample_data(psAsventres), "data.frame")
psAsventres.metadata$Phyogenetic_diversity <- df.pdAsventres$PD
#based on parameter Depth
kruskal.Asventres.PD_Depth <- kruskal.test(psAsventres.metadata$Phyogenetic_diversity ~ psAsventres.metadata$Depth)
print(kruskal.Asventres.PD_Depth)
### SeaWater 
psSeaWater <- subset_samples(ps1, Sample_types == "Seawater") #subset SeaWater
psSeaWater <- prune_taxa(taxa_sums(psSeawater) > 0, psSeawater)
otu_table_psSeaWater <- as.data.frame(psSeaWater@otu_table) 
metadata_table_psSeaWater  <- as.data.frame(psSeaWater@sam_data)
df.pdSeaWater <- pd(t(otu_table_psSeaWater), treefile_p1,include.root=F) 
datatable(df.pdSeaWater)
psSeaWater.metadata <- as(sample_data(psSeaWater), "data.frame")
psSeaWater.metadata$Phyogenetic_diversity <- df.pdSeaWater$PD
#based on parameter SampleDepth

#####pairwise Wkruskal.SeaWater.PD_Depth <- kruskal.test(psSeaWater.metadata$Phyogenetic_diversity ~ psSeaWater.metadata$Depth)
print(kruskal.SeaWater.PD_Depth)ilcox.test####
#all samples
pairwise.wilcox.test(ps1.metadata$Phyogenetic_diversity, ps1.metadata$Samples, p.adj = "BH")
pairwise.wilcox.test(ps1.metadata$Phyogenetic_diversity, ps1.metadata$Depth, p.adj = "BH")
#X. muta
pairwise.wilcox.test(psXmuta.metadata$Phyogenetic_diversity, psXmuta.metadata$Depth, p.adj = "BH")
#A.sventres
pairwise.wilcox.test(psAsventres.metadata$Phyogenetic_diversity, psAsventres.metadata$Depth, p.adj = "BH")
#SeaWater
pairwise.wilcox.test(psSeaWater.metadata$Phyogenetic_diversity, psSeaWater.metadata$Depth, p.adj = "BH")

####Beta diversity statistical test###
###ADONIS and betadisper test based on bray curtis####
ps1.hel <- transform(ps1, "hellinger") 
otu_table_ps1.hel <- as.data.frame(ps1.hel@otu_table)
metadf.hel <- data.frame(sample_data(ps1.hel))
#testing influence of depth 
permanova <- adonis(distance(ps1.hel, method="bray") ~ Depth,
       data = metadf.hel)
dist.hel <- vegdist(t(otu_table_ps1.hel))
ps.disper.hel <- betadisper(dist.hel, metadf.hel$Depth)
permutest(ps.disper.hel)

###subset for each each samples 
##X.muta
psXmuta_hel <- transform(psXmuta, "hellinger") 
otu_table_psXmuta.hel <- as.data.frame(psXmuta_hel@otu_table)
metadf.Xmuta_hel <- data.frame(sample_data(psXmuta_hel))
otu_table_psXmuta_hel <- as.data.frame(psXmuta_hel@otu_table)
adonis(distance(psXmuta_hel, method="bray") ~ Depth,
       data = metadf.Xmuta_hel)
dist.helXmuta <- vegdist(t(otu_table_psXmuta.hel)) 
ps.disper.Xmutahel <- betadisper(dist.helXmuta, metadf.Xmuta_hel$Depth) 
permutest(ps.disper.Xmutahel)
##A.sventres
psAsventres_hel <- transform(psAsventres, "hellinger") 
otu_table_psAsventres.hel <- as.data.frame(psAsventres_hel@otu_table)
metadf.Asventres_hel <- data.frame(sample_data(psAsventres_hel))
otu_table_psAsventres_hel <- as.data.frame(psAsventres_hel@otu_table)
adonis(distance(psAsventres_hel, method="bray") ~ Depth,
       data = metadf.Asventres_hel)
dist.helAsventres <- vegdist(t(otu_table_psAsventres_hel)) 
ps.disper.Asventres <- betadisper(dist.helAsventres, metadf.Asventres_hel$Depth) 
permutest(ps.disper.Asventres)
###SeaWater
psSeaWater_hel <- transform(psSeaWater, "hellinger") 
otu_table_psSeaWater_hel <- as.data.frame(psSeaWater_hel@otu_table)
metadf.SeaWater_hel <- data.frame(sample_data(psSeaWater_hel))
otu_table_psSeaWater_hel <- as.data.frame(psSeaWater_hel@otu_table)
adonis(distance(psSeaWater_hel, method="bray") ~ Depth,
       data = metadf.SeaWater_hel)
dist.helSeaWater <- vegdist(t(otu_table_psSeaWater_hel)) 
ps.disper.SeaWater <- betadisper(dist.helSeaWater, metadf.SeaWater_hel$Depth) 
permutest(ps.disper.SeaWater)

####script from Georg Pairwise Adonis#################
#Use the function pairwise.adonis() with following arguments
#
#x = community table
#
#factors = a column or vector with all factors to be tested pairwise
#
#sim.function = which function to calculate the similarity matrix. eg 'daisy' or 'vegdist' default is 'vegdist'
# NOTE that if you wnat to use daisy, you need to install library 'cluster'
#
#sim.method = similarity method from daisy or vegdist: default is 'bray' 
#
#p.adjust.m = the p.value correction method, one of the methods supported by p.adjust(); default is 'bonferroni'
#
#The function will return a table with the pairwise factors, F-values, R^2, p.value and adjusted p.value
#
# load in runnig R session with: 
# source('pairwise.adonis.txt')
#
# example:
# data(iris)
# pairwise.adonis(iris[,1:4],iris$Species)
#
#[1] "Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1"
#                    pairs    F.Model        R2 p.value p.adjusted sig
#1    setosa vs versicolor  552.51266 0.8493496   0.001      0.003   *
#2     setosa vs virginica 1001.54509 0.9108722   0.001      0.003   *
#3 versicolor vs virginica   91.82959 0.4837475   0.001      0.003   *

# similarity euclidean from vegdist and holm correction
# pairwise.adonis(x=iris[,1:4],factors=iris$Species,sim.function='vegdist',sim.method='euclidian',p.adjust.m='holm')

# similarity manhattan from daisy and bonferroni correction
# pairwise.adonis(x=iris[,1:4],factors=iris$Species,sim.function='daisy',sim.method='manhattan',p.adjust.m='bonferroni')
##############################Pairwise adonis all data together####################################
####prepare/import the data#####
ps1.hel <- transform(ps1, "hellinger") 
otu_table_ps1.hel <- as.data.frame(ps1.hel@otu_table)
#transpose abundance data
abundance <- as.data.frame(t(otu_table_ps1.hel))
# remove 0 OTUs / columns & rows
abundance <- droplevels(abundance[rowSums(abundance) != 0,])
abundance<- droplevels(abundance[,colSums(abundance) != 0])
read.csv(file="MappingFile.hab", header=TRUE, sep="\t", row.names=1) -> hab

##start copy here for function pairwise.adonis()
pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='bonferroni')
{
  library(vegan)
  
  co = combn(unique(as.character(factors)),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  
  for(elem in 1:ncol(co)){
    if(sim.function == 'daisy'){
      library(cluster); x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
    } else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
    
    ad = adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])] );
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'.'
  sig[p.adjusted <= 0.01] <-'*'
  sig[p.adjusted <= 0.001] <-'**'
  sig[p.adjusted <= 0.0001] <-'***'
  
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted,sig)
  print("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1")
  return(pairw.res)
  
} 

attach(hab)
detach(hab)

#pairwise adonis parameter "Depth"
pairwise.adonis(abundance, Depth) # this is the actual analysis "all" is your abundance table, "Species" your parameter, like "depth"

#pairwise adonis parameter "Sample_types"
pairwise.adonis(abundance, Sample_types)

####Pairwise Adonis subset X. muta##################
psXmuta <- subset_samples(ps1, Sample_types == "X.muta") #subset X. muta
psXmuta <- prune_taxa(taxa_sums(psXmuta) > 0, psXmuta)

otu_table_psXmuta<- as.data.frame(psXmuta@otu_table)
#transpose abundance data X. muta
abundance.Xmuta <- as.data.frame(t(otu_table_psXmuta.hel))
# remove 0 OTUs / columns & rows
abundance.Xmuta <- droplevels(abundance.Xmuta[rowSums(abundance.Xmuta) != 0,])
abundance.Xmuta<- droplevels(abundance.Xmuta[,colSums(abundance.Xmuta) != 0])

read.csv(file="MappingXmuta.hab", header=TRUE, sep="\t", row.names=1) -> hab.Xmuta

# re-run pairwise adonis (copy paste the script)
##start copy here for function pairwise.adonis()
pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='bonferroni')
{
  library(vegan)
  
  co = combn(unique(as.character(factors)),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  
  for(elem in 1:ncol(co)){
    if(sim.function == 'daisy'){
      library(cluster); x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
    } else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
    
    ad = adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])] );
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'.'
  sig[p.adjusted <= 0.01] <-'*'
  sig[p.adjusted <= 0.001] <-'**'
  sig[p.adjusted <= 0.0001] <-'***'
  
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted,sig)
  print("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1")
  return(pairw.res)
  
} 

attach(hab.Xmuta)
#pairwise adonis parameter "Depth" in X. muta
pairwise.adonis(abundance.Xmuta, Depth)

####Pairwise Adonis subset A. sventres##################
psAsventres <- subset_samples(ps1, Sample_types == "A.sventres") #subset A. sventres
psAsventres <- prune_taxa(taxa_sums(psAsventres) > 0, psAsventres)
psAsventres.hel <- transform(psAsventres, "hellinger") 
otu_table_psAsventres.hel <- as.data.frame(psAsventres.hel@otu_table)
read.csv(file="MappingAsventres.hab", header=TRUE, sep="\t", row.names=1) -> hab.Asventres
#transpose abundance data A.sventres
abundance.Asventres <- as.data.frame(t(otu_table_psAsventres.hel))
# remove 0 OTUs / columns & rows
abundance.Asventres <- droplevels(abundance.Asventres[rowSums(abundance.Asventres) != 0,])
abundance.Asventres <- droplevels(abundance.Asventres[,colSums(abundance.Asventres) != 0])
##start copy here for function pairwise.adonis()
pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='bonferroni')
{
  library(vegan)
  
  co = combn(unique(as.character(factors)),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  
  for(elem in 1:ncol(co)){
    if(sim.function == 'daisy'){
      library(cluster); x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
    } else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
    
    ad = adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])] );
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'.'
  sig[p.adjusted <= 0.01] <-'*'
  sig[p.adjusted <= 0.001] <-'**'
  sig[p.adjusted <= 0.0001] <-'***'
  
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted,sig)
  print("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1")
  return(pairw.res)
  
} 

attach(hab.Asventres)
detach(hab.Asventres)
#pairwise adonis parameter "Depth" in A.sventres
pairwise.adonis(abundance.Asventres, Depth)

####Pairwise Adonis subset Sea water##################
psSeaWater <- subset_samples(ps1, Sample_types == "Seawater") #subset SeaWater
psSeaWater <- prune_taxa(taxa_sums(psSeaWater) > 0, psSeaWater)
psSeaWater.hel <- transform(psSeaWater, "hellinger")
otu_table_psSeawater.hel <- as.data.frame(psSeaWater.hel@otu_table)
read.csv(file="MappingSW.hab", header=TRUE, sep="\t", row.names=1) -> hab.Seawater
#transpose abundance data Seawater
abundance.SW <- as.data.frame(t(otu_table_psSeawater.hel))
# remove 0 OTUs / columns & rows
abundance.SW <- droplevels(abundance.SW[rowSums(abundance.SW) != 0,])
abundance.SW <- droplevels(abundance.SW[,colSums(abundance.SW) != 0])
##start copy here for function pairwise.adonis()
pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='bonferroni')
{
  library(vegan)
  
  co = combn(unique(as.character(factors)),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  
  for(elem in 1:ncol(co)){
    if(sim.function == 'daisy'){
      library(cluster); x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
    } else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
    
    ad = adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])] );
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'.'
  sig[p.adjusted <= 0.01] <-'*'
  sig[p.adjusted <= 0.001] <-'**'
  sig[p.adjusted <= 0.0001] <-'***'
  
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted,sig)
  print("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1")
  return(pairw.res)
  
} 

attach(hab.Seawater)
detach(hab.Seawater)
#pairwise adonis parameter "Depth" in Sea water
pairwise.adonis(abundance.SW, Depth)


###ADONIS and betadisper test based on Unweighted UniFrac####
ps1.rel <- transform(ps1, "relative.abundance") 
metadf <- data.frame(sample_data(ps1.rel))
set.seed(28567)
uwUF.dist = UniFrac(ps1.rel, weighted=FALSE, normalized=TRUE)
#test based on depth 
UwUF.betadisper <- betadisper(uwUF.dist, metadf$Depth)
permutest(UwUF.betadisper)
#test based on samples 
adonis(uwUF.dist ~ Samples, data=metadf)
UwUF.betadisper <- betadisper(uwUF.dist, metadf$Samples)
permutest(UwUF.betadisper)
###subset for each each samples 
##X.muta
psXmuta_rel <- transform(psXmuta, "relative.abundance") 
metadf.Xmuta_rel <- data.frame(sample_data(psXmuta_rel))
set.seed(28567)
uwUF.Xmuta.dist <- UniFrac(psXmuta_rel, weighted=FALSE, normalized=TRUE)
#test based on depth 
adonis(uwUF.Xmuta.dist ~ Depth, data=metadf.Xmuta_rel)
UwUF.betadisper.Xmuta <- betadisper(uwUF.Xmuta.dist, metadf.Xmuta_rel$Depth)
permutest(UwUF.betadisper.Xmuta )
##A.sventres
psAsventres_rel <- transform(psAsventres, "relative.abundance") 
metadf.Asventres_rel <- data.frame(sample_data(psAsventres_rel))
set.seed(28567)
uwUF.Asventres.dist <- UniFrac(psAsventres_rel, weighted=FALSE, normalized=TRUE)
#test based on depth 
adonis(uwUF.Asventres.dist ~ Depth, data=metadf.Asventres_rel)
UwUF.betadisper.Asventres <- betadisper(uwUF.Asventres.dist, metadf.Asventres_rel$Depth)
permutest(UwUF.betadisper.Asventres)
##SeaWater
psSeaWater_rel <- transform(psSeaWater, "relative.abundance") 
metadf.SeaWater_rel <- data.frame(sample_data(psSeaWater_rel))
set.seed(28567)
uwUF.SeaWater.dist <- UniFrac(psSeaWater_rel, weighted=FALSE, normalized=TRUE)
#test based on depth 
adonis(uwUF.Asventres.dist ~ Depth, data=metadf.Asventres_rel)
UwUF.betadisper.Asventres <- betadisper(uwUF.Asventres.dist, metadf.Asventres_rel$Depth)
permutest(UwUF.betadisper.Asventres)
###ADONIS and betadisper test based on Weighted UniFrac####
ps1.rel <- transform(ps1, "relative.abundance") 
metadf <- data.frame(sample_data(ps1.rel))
set.seed(28567)
wUF.dist = UniFrac(ps1.rel, weighted=TRUE, normalized=TRUE)
#test based on depth 
adonis(wUF.dist ~ Depth, data=metadf)
wUF.betadisper <- betadisper(wUF.dist, metadf$Depth)
permutest(wUF.betadisper)
#test based on samples 
adonis(wUF.dist ~ Samples, data=metadf)
wUF.betadisper <- betadisper(wUF.dist, metadf$Samples)
permutest(wUF.betadisper)
###subset for each each samples 
##X.muta
psXmuta_rel <- transform(psXmuta, "relative.abundance") 
metadf.Xmuta_rel <- data.frame(sample_data(psXmuta_rel))
set.seed(28567)
wUF.Xmuta.dist <- UniFrac(psXmuta_rel, weighted=FALSE, normalized=TRUE)
#test based on depth 
adonis(wUF.Xmuta.dist ~ Depth, data=metadf.Xmuta_rel)
wUF.betadisper.Xmuta <- betadisper(wUF.Xmuta.dist, metadf.Xmuta_rel$Depth)
permutest(wUF.betadisper.Xmuta )
##Asventres
psAsventres_rel <- transform(psAsventres, "relative.abundance") 
metadf.Asventres_rel <- data.frame(sample_data(psAsventres_rel))
set.seed(28567)
wUF.Asventres.dist <- UniFrac(psAsventres_rel, weighted=FALSE, normalized=TRUE)
#test based on depth 
adonis(wUF.Asventres.dist ~ Depth, data=metadf.Asventres_rel)
wUF.betadisper.Asventres <- betadisper(wUF.Asventres.dist, metadf.Asventres_rel$Depth)
permutest(wUF.betadisper.Asventres)
##SeaWater
psSeaWater_rel <- transform(psSeaWater, "relative.abundance") 
metadf.SeaWater_rel <- data.frame(sample_data(psSeaWater_rel))
set.seed(28567)
wUF.SeaWater.dist <- UniFrac(psSeaWater_rel, weighted=FALSE, normalized=TRUE)
#test based on depth 
adonis(wUF.SeaWater.dist ~ Depth, data=metadf.SeaWater_rel)
wUF.betadisper.SeaWater <- betadisper(wUF.SeaWater.dist, metadf.SeaWater_rel$Depth)
permutest(wUF.betadisper.SeaWater)

###HeatMap####
#Sort the OTUs by abundance and pick the top 20
psHeatMap.rel <- transform(psHeatMap, "relative.abundance")
#Heat map based on top 100 otus / cut of average total relative abundance >= 0.25%
top100OTU.names <-  names(sort(taxa_sums(psHeatMap.rel), TRUE)[1:100])
top100OTU <-  prune_taxa(top100OTU.names, psHeatMap.rel) 
heatmap <- otu_table(top100OTU, taxa_are_rows)#generate OTU_table containing information of selected OTU, in this case otu with total average >=0.25%
write.csv(heatmap, file = "HeatMap.csv", row.names=TRUE)# save in csv and edit manually e.g. order otu based on same Phylum level
HeatMap <-read.csv("HeatMap.csv", row.names=1, check.names=F, stringsAsFactors=F, header = TRUE)#import heatmap file after re-arrangenment (edited version)
color = c("white", "gray", "green","skyblue", "royalblue", "blue", "purple", "magenta","violet", "red")
breaks = c(-0.5, 0, 0.005, 0.01, 0.05, 0.09, 0.14, 0.20, 0.40, 0.60, 0.80)# 0 = white, 0-0.5% = gray, 0.5%-1% = green, 1%-5% = skyblue, 5%-9%=royalblue, 9%-14% = blue, 14%-19% =purple, 19%-25% = magenta, 25%-45%= violet, 45%-75% = red
pheatmap(HeatMap, cluster_col = FALSE, cluster_row = FALSE, color=color, breaks=breaks, fontsize_row = 8, main = "HeatMap of OTU based total average relative abundance >= 0.25% among all samples")#edit picture in adobe for a better clarity of color scale


##newheatmap
HeatMap <-read.csv("HeatMap.csv", row.names=1, check.names=F, stringsAsFactors=F, header = TRUE)#import heatmap file after re-arrangenment (edited version)
color = c("white", "snow2", "snow3","snow4", "steelblue2", "royalblue2", "royalblue4","tomato", "orangered3", "red4")
breaks = c(-0.5, 0, 0.005, 0.01, 0.05, 0.09, 0.14, 0.20, 0.40, 0.60, 0.80)# 0 = white, 0-0.5% = gray, 0.5%-1% = green, 1%-5% = skyblue, 5%-9%=royalblue, 9%-14% = blue, 14%-19% =purple, 19%-25% = magenta, 25%-45%= violet, 45%-75% = red
pheatmap(HeatMap, cluster_col = FALSE, cluster_row = FALSE, color=color, breaks=breaks, fontsize_row = 8, main = "HeatMap of OTU based total average relative abundance >= 0.25% among all samples")#edit picture in adobe for a better clarity of color scale



#######Wilcoxon test antimicrobial####
###E. coli 
X.muta_Ecoli <-c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3.72, 0, 3.27,0,3.28)
A.sventres_Ecoli <-c(NA, NA, NA, NA, NA, 0, 0, 0, 0, 3.22, 4.15, 0, 0, 3.96, 0) 
#create data frame
data_Ecoli <- data.frame( 
  group = rep(c("X.muta", "A.sventres"), each = 15),
  inhibition = c(X.muta_Ecoli,  A.sventres_Ecoli)
)
#check data
print(data_Ecoli)
#Wilcoxon test E. coli
res_Ecoli <- wilcox.test(inhibition ~ group, data = data_Ecoli,
                   exact = FALSE)
res_Ecoli

###A. salmonicida 
X.muta_salmonicida <-c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3.47, 0)
A.sventres_salmonicida <-c(NA, NA, NA, NA, NA, 0, 0, 0, 0, 4.47, 0, 0, 0, 6.32, 0) 

#create data frame
data_salmonicida <- data.frame( 
  group = rep(c("X.muta", "A.sventres"), each = 15),
  inhibition = c(X.muta_salmonicida,  A.sventres_salmonicida)
)
#check data
print(data_salmonicida)

##Wilcoxon test Salmonicida
res_Salmonicida <- wilcox.test(inhibition ~ group, data = data_salmonicida,
                         exact = FALSE)
res_Salmonicida

###B. subtilis
X.muta_subtilis <-c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3.14, 0)
A.sventres_subtilis <-c(NA, NA, NA, NA, NA, 0, 0, 0, 0, 5.27, 5.59, 0, 4.01, 5.93, 5) 

#create data frame
data_subtilis <- data.frame( 
  group = rep(c("X.muta", "A.sventres"), each = 15),
  inhibition = c(X.muta_subtilis,  A.sventres_subtilis)
)
#check data
print(data_subtilis)

##Wilcoxon test Salmonicida
res_subtilis<- wilcox.test(inhibition ~ group, data = data_subtilis,
                               exact = FALSE)
res_subtilis

###S. simulans
X.muta_simulans <-c(0, 0, 0, 0, 0, 0, 3.2, 0, 0, 0, 0, 0, 0, 3.1, 0)
A.sventres_simulans <-c(NA, NA, NA, NA, NA, 0, 0, 0, 4.27, 5.86, 6.86, 3.15, 4.68, 5.5, 5.63) 

#create data frame
data_simulans <- data.frame( 
  group = rep(c("X.muta", "A.sventres"), each = 15),
  inhibition = c(X.muta_simulans,  A.sventres_simulans)
)
#check data
print(data_simulans)

##Wilcoxon test Simulans
res_simulans<- wilcox.test(inhibition ~ group, data = data_simulans,
                           exact = FALSE)
res_simulans

##Wilcoxon test S.parasitica
X.muta_parasitica <-c(0, 5.55, 13.43, 6.51, 9.38, 0, 0, 0, 0, 0, 8.14, 8.09, 0, 0, 0)
A.sventres_parasitica <-c(NA, NA, NA, NA, NA, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0) 

#create data frame
data_parasitica <- data.frame( 
  group = rep(c("X.muta", "A.sventres"), each = 15),
  inhibition = c(X.muta_parasitica,  A.sventres_parasitica)
)
#check data
print(data_parasitica)

res_parasitica<- wilcox.test(inhibition ~ group, data = data_parasitica,
                           exact = FALSE)
res_parasitica

######Wilcoxon test antimicrobial depth vs strains ####
###E. coli 
Deep_Ecoli <-c(0, 0, 0, 0, 0, NA, NA, NA, NA, NA)
Middle_Ecoli <-c(0, 0, 0, 0, 0, 0, 0, 0, 0, 3.22) 
Shallow_Ecoli <-c(3.72, 0, 3.27, 0, 3.28, 4.15, 0, 0, 3.96, 0)

#create data frame Deep_Middle
data_depthEcoli1<- data.frame( 
  group = rep(c("Deep", "Middle"), each = 10),
  inhibition = c(Deep_Ecoli,  Middle_Ecoli)
)
#check data 
print(data_depthEcoli1)

#wilcoxon Deep_Middle Ecoli
res_depthEcoli1<- wilcox.test(inhibition ~ group, data = data_depthEcoli1,
                             exact = FALSE)
res_depthEcoli1

#create data frame Deep_Shallow
data_depthEcoli2<- data.frame(group = rep(c("Deep", "Shallow"), each = 10),
  inhibition = c(Deep_Ecoli,  Shallow_Ecoli)) 

#check data 
print(data_depthEcoli2)

#wilcoxon Deep_Shallow Ecoli
res_depthEcoli2<- wilcox.test(inhibition ~ group, data = data_depthEcoli2,
                              exact = FALSE)
res_depthEcoli2

#create data frame Middle_Shallow
data_depthEcoli3<- data.frame(group = rep(c("Middle", "Shallow"), each = 10),
                              inhibition = c(Middle_Ecoli,  Shallow_Ecoli)) 

#check data 
print(data_depthEcoli3)

#wilcoxon Deep_Shallow Ecoli
res_depthEcoli3<- wilcox.test(inhibition ~ group, data = data_depthEcoli3,
                              exact = FALSE)
res_depthEcoli3

###A. salmonicida 
Deep_salmonicida <-c(0, 0, 0, 0, 0, NA, NA, NA, NA, NA)
Middle_salmonicida <-c(0, 0, 0, 0, 0, 0, 0, 0, 0, 4.47) 
Shallow_salmonicida <-c(0, 0, 0, 3.47, 0, 0, 0, 0, 6.32, 0)

#create data frame Deep_Middle
data_depthsalmonicida1<- data.frame( 
  group = rep(c("Deep", "Middle"), each = 10),
  inhibition = c(Deep_salmonicida,  Middle_salmonicida )
)
#check data 
print(data_depthsalmonicida1)

#wilcoxon Deep_Middle Salmonicida
res_depthSalmonicida1<- wilcox.test(inhibition ~ group, data = data_depthsalmonicida1,
                              exact = FALSE)
res_depthSalmonicida1

#create data frame Deep_Shallow
data_depthsalmonicida2<- data.frame( 
  group = rep(c("Deep", "Shallow"), each = 10),
  inhibition = c(Deep_salmonicida,  Shallow_salmonicida )
)
#check data 
print(data_depthsalmonicida2)

#wilcoxon Deep_Shallow Salmonicida
res_depthSalmonicida2<- wilcox.test(inhibition ~ group, data = data_depthsalmonicida2,
                                    exact = FALSE)

####create data frame Middle_Shallow
data_depthsalmonicida3<- data.frame( 
  group = rep(c("Middle", "Shallow"), each = 10),
  inhibition = c(Middle_salmonicida,  Shallow_salmonicida )
)
#check data 
print(data_depthsalmonicida3)

#wilcoxon Middle_Shallow Salmonicida
res_depthSalmonicida3<- wilcox.test(inhibition ~ group, data = data_depthsalmonicida3,
                                    exact = FALSE)
#####depth vs B.subtilis 
Deep_subtilis <-c(0, 0, 0, 0, 0, NA, NA, NA, NA, NA)
Middle_subtilis <-c(0, 0, 0, 0, 0, 0, 0, 0, 0, 5.27) 
Shallow_subtilis <-c(0, 0, 0, 3.14, 0, 5.59, 0, 4.01, 5.93, 5)

#create data frame Deep_Middle
data_depthsubtilis1<- data.frame( 
  group = rep(c("Deep", "Middle"), each = 10),
  inhibition = c(Deep_subtilis,  Middle_subtilis )
)
#check data 
print(data_depthsubtilis1)

#wilcoxon Deep_Middle Subtilis
res_depthsubtilis1<- wilcox.test(inhibition ~ group, data = data_depthsubtilis1,
                                    exact = FALSE)

#create data frame Deep_shallow
data_depthsubtilis2<- data.frame( 
  group = rep(c("Deep", "Shallow"), each = 10),
  inhibition = c(Deep_subtilis,  Shallow_subtilis )
)
#check data 
print(data_depthsubtilis2)

#wilcoxon Deep_Shallow Subtilis
res_depthsubtilis2<- wilcox.test(inhibition ~ group, data = data_depthsubtilis2,
                                 exact = FALSE)

#create data frame Middle_Shallow
data_depthsubtilis3<- data.frame( 
  group = rep(c("Middle", "Shallow"), each = 10),
  inhibition = c(Middle_subtilis,  Shallow_subtilis )
)
#check data 
print(data_depthsubtilis3)

#wilcoxon Middle_Shallow Subtilis
res_depthsubtilis3<- wilcox.test(inhibition ~ group, data = data_depthsubtilis3,
                                 exact = FALSE)

##depth vs S. simulans 
Deep_simulans <-c(0, 0, 0, 0, 0, NA, NA, NA, NA, NA)
Middle_simulans <-c(0, 3.2, 0, 0, 0, 0, 0, 0, 4.07, 5.86) 
Shallow_simulans <-c(0, 0, 0, 3.1, 0, 6.86, 3.15, 4.68, 5.5, 5.63)

#create data frame Deep_Middle
data_depthsimulans1<- data.frame( 
  group = rep(c("Deep", "Middle"), each = 10),
  inhibition = c(Deep_simulans,  Middle_simulans )
)
#check data 
print(data_depthsimulans1)

#wilcoxon Deep_Middle Subtilis
res_depthsimulans1<- wilcox.test(inhibition ~ group, data = data_depthsimulans1,
                                 exact = FALSE)

#create data frame Deep_Shallow
data_depthsimulans2<- data.frame( 
  group = rep(c("Deep", "Shallow"), each = 10),
  inhibition = c(Deep_simulans,  Shallow_simulans )
)
#check data 
print(data_depthsimulans2)

#wilcoxon Deep_Middle Subtilis
res_depthsimulans2<- wilcox.test(inhibition ~ group, data = data_depthsimulans2,
                                 exact = FALSE)

#create data frame Middle_Shallow
data_depthsimulans3<- data.frame( 
  group = rep(c("Middle", "Shallow"), each = 10),
  inhibition = c(Middle_simulans,  Shallow_simulans )
)
#check data 
print(data_depthsimulans3)

#wilcoxon Deep_Middle Subtilis
res_depthsimulans3<- wilcox.test(inhibition ~ group, data = data_depthsimulans3,
                                 exact = FALSE)

##depth vs S. parasitica 
Deep_parasitica <-c(0, 5.55, 13.43, 6.51, 9.38, NA, NA, NA, NA, NA)
Middle_parasitica <-c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0) 
Shallow_parasitica <-c(8.14, 8.09, 0, 0, 0, 0, 0, 0, 0, 0)

#create data frame Deep_Middle
data_depthparasitica1<- data.frame( 
  group = rep(c("Deep", "Middle"), each = 10),
  inhibition = c(Deep_parasitica,  Middle_parasitica)
)
#check data 
print(data_depthparasitica1)

#wilcoxon Deep_Middle Subtilis
res_depthparasitica1<- wilcox.test(inhibition ~ group, data = data_depthparasitica1,
                                 exact = FALSE)
#create data frame Deep_Shallow
data_depthparasitica2<- data.frame( 
  group = rep(c("Deep", "Shallow"), each = 10),
  inhibition = c(Deep_parasitica,  Shallow_parasitica)
)
#check data 
print(data_depthparasitica2)

#wilcoxon Deep_Middle Subtilis
res_depthparasitica2<- wilcox.test(inhibition ~ group, data = data_depthparasitica2,
                                   exact = FALSE)

#create data frame Middle_Shallow
data_depthparasitica3<- data.frame( 
  group = rep(c("Middle", "Shallow"), each = 10),
  inhibition = c(Middle_parasitica,  Shallow_parasitica)
)
#check data 
print(data_depthparasitica3)

#wilcoxon Middle_Shallow Subtilis
res_depthparasitica3<- wilcox.test(inhibition ~ group, data = data_depthparasitica3,
                                   exact = FALSE)

####Wilcoxon X. muta###
###E. coli
DeepMuta_Ecoli <-c(0, 0, 0, 0, 0)
MiddleMuta_Ecoli <-c(0, 0, 0, 0,0) 
ShallowMuta_Ecoli <-c(3.72, 0, 3.27, 0, 3.28)

#create data frame Deep_Middle X. muta E. coli
data_XmutaEcoli<- data.frame( 
  group = rep(c("Deep", "Middle"), each = 5),
  inhibition = c(DeepMuta_Ecoli,  MiddleMuta_Ecoli)
)
#check data 
print(data_XmutaEcoli)

#wilcoxon Xmuta Deep_Shallow E coli
res_XmutaEcoli<- wilcox.test(inhibition ~ group, data = data_XmutaEcoli,
                                 exact = FALSE)

#create data frame Deep_Shallow X. muta E. coli
data_XmutaEcoli<- data.frame( 
  group = rep(c("Deep", "Shallow"), each = 5),
  inhibition = c(DeepMuta_Ecoli,  ShallowMuta_Ecoli)
)
#check data 
print(data_XmutaEcoli)

#wilcoxon E.coli Deep_Shallow
res_XmutaEcoli<- wilcox.test(inhibition ~ group, data = data_XmutaEcoli,
                             exact = FALSE)

#create data frame Middle_Shallow
data_XmutaEcoli<- data.frame( 
  group = rep(c("Middle", "Shallow"), each = 5),
  inhibition = c(MiddleMuta_Ecoli,  ShallowMuta_Ecoli)
)
#check data 
print(data_XmutaEcoli)

#wilcoxon Middle_Shallow Subtilis
res_XmutaEcoli<- wilcox.test(inhibition ~ group, data = data_XmutaEcoli,
                             exact = FALSE)

#####Salmonicida
DeepMuta_Salmonicida <-c(0, 0, 0, 0, 0)
MiddleMuta_Salmonicida <-c(0, 0, 0, 0,0) 
ShallowMuta_Salmonicida <-c(0, 0, 0, 3.47, 0)

#create data frame Deep_Middle X. muta Salmonicida
data_XmutaSalmonicida<- data.frame( 
  group = rep(c("Deep", "Middle"), each = 5),
  inhibition = c(DeepMuta_Salmonicida,  MiddleMuta_Salmonicida)
)
#check data 
print(data_XmutaSalmonicida)

#wilcoxon Xmuta Deep_Shallow Salmonicida
res_XmutaSalmonicida<- wilcox.test(inhibition ~ group, data = data_XmutaSalmonicida,
                             exact = FALSE)

#create data frame Deep_Shallow X. muta Salmonicida
data_XmutaSalmonicida<- data.frame( 
  group = rep(c("Middle", "Shallow"), each = 5),
  inhibition = c(MiddleMuta_Salmonicida,  ShallowMuta_Salmonicida)
)
#check data 
print(data_XmutaSalmonicida)

#wilcoxon Xmuta Deep_Shallow E coli
res_XmutaSalmonicida<- wilcox.test(inhibition ~ group, data = data_XmutaSalmonicida,
                                   exact = FALSE)
#####Subtilis
DeepMuta_Subtilis <-c(0, 0, 0, 0, 0)
MiddleMuta_Subtilis <-c(0, 0, 0, 0,0) 
ShallowMuta_Subtilis <-c(0, 0, 0, 0, 3.14)

#create data frame Deep_Middle X. muta Salmonicida
data_XmutaSubtilis<- data.frame( 
  group = rep(c("Deep", "Middle"), each = 5),
  inhibition = c(DeepMuta_Subtilis,  MiddleMuta_Subtilis)
)
#check data 
print(data_XmutaSubtilis)

#wilcoxon Xmuta Deep_Middle Salmonicida
res_XmutaSubtilis<- wilcox.test(inhibition ~ group, data = data_XmutaSubtilis,
                                   exact = FALSE)

#create data frame Deep_Middle X. muta Salmonicida
data_XmutaSubtilis<- data.frame( 
  group = rep(c("Deep", "Shallow"), each = 5),
  inhibition = c(DeepMuta_Subtilis,  ShallowMuta_Subtilis)
)
#check data 
print(data_XmutaSubtilis)

#wilcoxon Xmuta Deep_Middle Salmonicida
res_XmutaSubtilis<- wilcox.test(inhibition ~ group, data = data_XmutaSubtilis,
                                exact = FALSE)

#create data frame Middle_Shallow X. muta Salmonicida
data_XmutaSubtilis<- data.frame( 
  group = rep(c("Middle", "Shallow"), each = 5),
  inhibition = c(MiddleMuta_Subtilis,  ShallowMuta_Subtilis)
)
#check data 
print(data_XmutaSubtilis)

#wilcoxon Xmuta Deep_Shallow Salmonicida
res_XmutaSubtilis<- wilcox.test(inhibition ~ group, data = data_XmutaSubtilis,
                                exact = FALSE)

#####Simulans
DeepMuta_Simulans <-c(0, 0, 0, 0, 0)
MiddleMuta_Simulans <-c(0, 3.2, 0, 0,0) 
ShallowMuta_Simulans <-c(0, 0, 0, 3.1, 0)

#create data frame Deep_Middle X. muta Simulans
data_XmutaSimulans<- data.frame( 
  group = rep(c("Deep", "Middle"), each = 5),
  inhibition = c(DeepMuta_Simulans,  MiddleMuta_Simulans)
)
#check data 
print(data_XmutaSimulans)

#wilcoxon Xmuta Deep_Middle Simulans
res_XmutaSimulans<- wilcox.test(inhibition ~ group, data = data_XmutaSimulans,
                                exact = FALSE)

#create data frame Deep_Middle X. muta Simulans
data_XmutaSimulans<- data.frame( 
  group = rep(c("Deep", "Shallow"), each = 5),
  inhibition = c(DeepMuta_Simulans,  ShallowMuta_Simulans)
)
#check data 
print(data_XmutaSimulans)

#wilcoxon Xmuta Deep_Middle Simulans
res_XmutaSimulans<- wilcox.test(inhibition ~ group, data = data_XmutaSimulans,
                                exact = FALSE)
#create data frame Middle_Shallow X. muta Simulans
data_XmutaSimulans<- data.frame( 
  group = rep(c("Middle", "Shallow"), each = 5),
  inhibition = c(MiddleMuta_Simulans,  ShallowMuta_Simulans)
)
#check data 
print(data_XmutaSimulans)

#wilcoxon Xmuta Deep_Middle Simulans
res_XmutaSimulans<- wilcox.test(inhibition ~ group, data = data_XmutaSimulans,
                                exact = FALSE)
#####parasitica
DeepMuta_parasitica <-c(0, 5.55, 13.43, 6.51, 9.38)
MiddleMuta_parasitica <-c(0, 0, 0, 0,0) 
ShallowMuta_parasitica <-c(8.14, 8.09, 0, 0, 0)

#create data frame Deep_Middle X. muta parasitica
data_Xmutaparasitica<- data.frame( 
  group = rep(c("Deep", "Middle"), each = 5),
  inhibition = c(DeepMuta_parasitica,  MiddleMuta_parasitica)
)
#check data 
print(data_Xmutaparasitica)

#wilcoxon Xmuta Deep_Middle Simulans
res_Xmutaparasitica<- wilcox.test(inhibition ~ group, data = data_Xmutaparasitica,
                                exact = FALSE)

#create data frame Deep_Shallow X. muta parasitica
data_Xmutaparasitica<- data.frame( 
  group = rep(c("Deep", "Shallow"), each = 5),
  inhibition = c(DeepMuta_parasitica,  ShallowMuta_parasitica)
)
#check data 
print(data_Xmutaparasitica)

#wilcoxon Xmuta Deep_Middle Simulans
res_Xmutaparasitica<- wilcox.test(inhibition ~ group, data = data_Xmutaparasitica,
                                  exact = FALSE)

#create data frame Middle_Shallow X. muta parasitica
data_Xmutaparasitica<- data.frame( 
  group = rep(c("Middle", "Shallow"), each = 5),
  inhibition = c(MiddleMuta_parasitica,  ShallowMuta_parasitica)
)
#check data 
print(data_Xmutaparasitica)

#wilcoxon Xmuta Middle_Shallow Simulans
res_Xmutaparasitica<- wilcox.test(inhibition ~ group, data = data_Xmutaparasitica,
                                  exact = FALSE)

###A. sventres
#####E. coli
MiddleSventres_Ecoli <-c(0, 0, 0, 0,3.22) 
ShallowSventres_Ecoli <-c(4.15, 0, 0, 3.96, 0)

#create data frame Deep_Middle X. muta parasitica
data_SventresEcoli<- data.frame( 
  group = rep(c("Middle", "Shallow"), each = 5),
  inhibition = c(MiddleSventres_Ecoli,  ShallowSventres_Ecoli)
)
#check data 
print(data_SventresEcoli)

#wilcoxon A.sventres Deep_Middle Simulans
res_SventresEcoli<- wilcox.test(inhibition ~ group, data = data_SventresEcoli,
                                  exact = FALSE)

#####A. salmonicida
MiddleSventres_Salmonicida <-c(0, 0, 0, 0,4.47) 
ShallowSventres_Salmonicida <-c(0, 0, 0, 6.32, 0)

#create data frame Deep_Middle A. sventres Salmonicida
data_SventresSalmonicida <- data.frame( 
  group = rep(c("Middle", "Shallow"), each = 5),
  inhibition = c(MiddleSventres_Salmonicida,  ShallowSventres_Salmonicida)
)
#check data 
print(data_SventresSalmonicida)

#wilcoxon Xmuta Deep_Middle Simulans
res_SventresSalmonicida<- wilcox.test(inhibition ~ group, data = data_SventresSalmonicida,
                                exact = FALSE)
#####B. subtilis
MiddleSventres_Subtilis <-c(0, 0, 0, 0,5.27) 
ShallowSventres_Subtilis <-c(5.59, 0, 4.01, 5.93, 5)

#create data frame Deep_Middle A. sventres Subtilis
data_SventresBsubtilis <- data.frame( 
  group = rep(c("Middle", "Shallow"), each = 5),
  inhibition = c(MiddleSventres_Subtilis,  ShallowSventres_Subtilis)
)
#check data 
print(data_SventresBsubtilis)

#wilcoxon Xmuta Deep_Middle Simulans
res_SventresBsubtilis<- wilcox.test(inhibition ~ group, data = data_SventresBsubtilis,
                                      exact = FALSE)
#####S. simulans
MiddleSventres_Simulans <-c(0, 0, 0, 4.07,5.86) 
ShallowSventres_Simulans <-c(6.86, 3.15, 4.68, 5.5, 5)

#create data frame Deep_Middle A. sventres Subtilis
data_SventresSsimulans <- data.frame( 
  group = rep(c("Middle", "Shallow"), each = 5),
  inhibition = c(MiddleSventres_Simulans,  ShallowSventres_Simulans)
)
#check data 
print(data_SventresSsimulans)

#wilcoxon Xmuta Deep_Middle Simulans
res_SventresSsimulans <- wilcox.test(inhibition ~ group, data = data_SventresSsimulans,
                                    exact = FALSE)

######ANOVA sponge extract####
#############X.muta###########
#import file
read.csv(file="Xmuta_Ecoli.txt", header=TRUE, sep="\t", row.names=1) -> my_data
#show the levels
levels(my_data$Depth)
library(dplyr)
group_by(my_data, Depth) %>%
  summarise(
    count = n(),
    mean = mean(Inhibition, na.rm = TRUE),
    sd = sd(Inhibition, na.rm = TRUE)
  )
# Compute the analysis of variance
res.aov <- aov(Inhibition ~ Depth, data = my_data)

# Summary of the analysis
summary(res.aov)
TukeyHSD(res.aov)

##Salmonicida
read.csv(file="XSalmonicida.txt", header=TRUE, sep="\t", row.names=1) -> my_data
#show the levels
levels(my_data$Depth)
library(dplyr)
group_by(my_data, Depth) %>%
  summarise(
    count = n(),
    mean = mean(Inhibition, na.rm = TRUE),
    sd = sd(Inhibition, na.rm = TRUE)
  )
# Compute the analysis of variance
res.aov <- aov(Inhibition ~ Depth, data = my_data)

# Summary of the analysis
summary(res.aov)
TukeyHSD(res.aov)

###Xmuta B.subtilis
read.csv(file="XSubtilis.txt", header=TRUE, sep="\t", row.names=1) -> my_data
#show the levels
levels(my_data$Depth)
library(dplyr)
group_by(my_data, Depth) %>%
  summarise(
    count = n(),
    mean = mean(Inhibition, na.rm = TRUE),
    sd = sd(Inhibition, na.rm = TRUE)
  )
# Compute the analysis of variance
res.aov <- aov(Inhibition ~ Depth, data = my_data)

# Summary of the analysis
summary(res.aov)
TukeyHSD(res.aov)

###Xmuta S.simulans
read.csv(file="XSimulans.txt", header=TRUE, sep="\t", row.names=1) -> my_data
#show the levels
levels(my_data$Depth)
library(dplyr)
group_by(my_data, Depth) %>%
  summarise(
    count = n(),
    mean = mean(Inhibition, na.rm = TRUE),
    sd = sd(Inhibition, na.rm = TRUE)
  )
# Compute the analysis of variance
res.aov <- aov(Inhibition ~ Depth, data = my_data)

# Summary of the analysis
summary(res.aov)
TukeyHSD(res.aov)

###Xmuta parasitica
read.csv(file="XParasitica.txt", header=TRUE, sep="\t", row.names=1) -> my_data
#show the levels
levels(my_data$Depth)
library(dplyr)
group_by(my_data, Depth) %>%
  summarise(
    count = n(),
    mean = mean(Inhibition, na.rm = TRUE),
    sd = sd(Inhibition, na.rm = TRUE)
  )
# Compute the analysis of variance
res.aov <- aov(Inhibition ~ Depth, data = my_data)

# Summary of the analysis
summary(res.aov)
TukeyHSD(res.aov)

####A. sventres####
Asventres_Ecoli
read.csv(file="AEcoli.txt", header=TRUE, sep="\t", row.names=1) -> my_data
#show the levels
levels(my_data$Depth)
library(dplyr)
group_by(my_data, Depth) %>%
  summarise(
    count = n(),
    mean = mean(Inhibition, na.rm = TRUE),
    sd = sd(Inhibition, na.rm = TRUE)
  )
# Compute the analysis of variance
res.aov <- aov(Inhibition ~ Depth, data = my_data)

# Summary of the analysis
summary(res.aov)
TukeyHSD(res.aov)

#Asventres_salmonicida
read.csv(file="ASalmonicida.txt", header=TRUE, sep="\t", row.names=1) -> my_data
#show the levels
levels(my_data$Depth)
library(dplyr)
group_by(my_data, Depth) %>%
  summarise(
    count = n(),
    mean = mean(Inhibition, na.rm = TRUE),
    sd = sd(Inhibition, na.rm = TRUE)
  )
# Compute the analysis of variance
res.aov <- aov(Inhibition ~ Depth, data = my_data)

# Summary of the analysis
summary(res.aov)
TukeyHSD(res.aov)

#Asventres_subtilis
read.csv(file="ASubtilis.txt", header=TRUE, sep="\t", row.names=1) -> my_data
#show the levels
levels(my_data$Depth)
library(dplyr)
group_by(my_data, Depth) %>%
  summarise(
    count = n(),
    mean = mean(Inhibition, na.rm = TRUE),
    sd = sd(Inhibition, na.rm = TRUE)
  )
# Compute the analysis of variance
res.aov <- aov(Inhibition ~ Depth, data = my_data)

# Summary of the analysis
summary(res.aov)
TukeyHSD(res.aov)

##Asventres_simulans
read.csv(file="ASimulans.txt", header=TRUE, sep="\t", row.names=1) -> my_data
#show the levels
levels(my_data$Depth)
library(dplyr)
group_by(my_data, Depth) %>%
  summarise(
    count = n(),
    mean = mean(Inhibition, na.rm = TRUE),
    sd = sd(Inhibition, na.rm = TRUE)
  )
# Compute the analysis of variance
res.aov <- aov(Inhibition ~ Depth, data = my_data)

# Summary of the analysis
summary(res.aov)
TukeyHSD(res.aov)

###C. oleophila
read.csv(file="XOleophila.txt", header=TRUE, sep="\t", row.names=1) -> my_data
#show the levels
levels(my_data$Depth)
library(dplyr)
group_by(my_data, Depth) %>%
  summarise(
    count = n(),
    mean = mean(Inhibition, na.rm = TRUE),
    sd = sd(Inhibition, na.rm = TRUE)
  )
# Compute the analysis of variance
res.aov <- aov(Inhibition ~ Depth, data = my_data)

# Summary of the analysis
summary(res.aov)
TukeyHSD(res.aov)


#######PERMANOVA pairwise seawater
trans <- "hellinger"
ps1.hel <- microbiome::transform(ps1, trans)
otu_table_ps1.hel <- as.data.frame(ps1.hel@otu_table)#otu_table
View(otu_table_ps1.hel)#to check how transformation to relative abundance looks like 
write.csv(otu_table_ps1.hel, file = "otu_table_ps1.hel.csv", row.names=TRUE)# save in csv
ordu.bray.hel <- ordinate(ps1.hel, "PCoA", "bray", weighted=F)
bray.hel <- plot_ordination(ps1.hel, ordu.bray.hel, color="Sample_types", shape="Depth")
bray_hel <- bray.hel + scale_fill_manual(values = c("#CBD588", "#5F7FC7", "orange","#DA5724", "#508578")) + geom_point(size = 6) + theme_bw()
print(bray_hel)
otu_table_ps1.hel <- as.data.frame(ps1.hel@otu_table)
metadf.hel <- data.frame(sample_data(ps1.hel))
#testing influence of depth 
permanova <- adonis(distance(ps1.hel, method="bray") ~ Depth,
                    data = metadf.hel)
dist.hel <- vegdist(t(otu_table_ps1.hel))
ps.disper.hel <- betadisper(dist.hel, metadf.hel$Depth)
permutest(ps.disper.hel)
###seawater
psSeaWater_hel <-subset_samples(ps1.hel, Sample_types == "Seawater") 
psSeaWater_hel <- prune_taxa(taxa_sums(psSeaWater_hel) > 0, psSeaWater_hel)
otu_table_psSeaWater_hel <- as.data.frame(psSeaWater_hel@otu_table)
metadf.SeaWater_hel <- data.frame(sample_data(psSeaWater_hel))
otu_table_psSeaWater_hel <- as.data.frame(psSeaWater_hel@otu_table)
adonis(distance(psSeaWater_hel, method="bray") ~ Depth,
       data = metadf.SeaWater_hel)
dist.helSeaWater <- vegdist(t(otu_table_psSeaWater_hel)) 
ps.disper.SeaWater <- betadisper(dist.helSeaWater, metadf.SeaWater_hel$Depth) 
permutest(ps.disper.SeaWater)

###pairwise Adonis seawater_deep_middle#####
psSeaWater_hel1 <-subset_samples(psSeaWater_hel, Depth != "shallow")
psSeaWater_hel1 <- prune_taxa(taxa_sums(psSeaWater_hel1) > 0, psSeaWater_hel1)
otu_table_psSeaWater_hel1 <- as.data.frame(psSeaWater_hel1@otu_table)
metadf.SeaWater_hel1 <- data.frame(sample_data(psSeaWater_hel1))
otu_table_psSeaWater_hel1 <- as.data.frame(psSeaWater_hel1@otu_table)
adonis(distance(psSeaWater_hel1, method="bray") ~ Depth,
       data = metadf.SeaWater_hel1)
dist.helSeaWater1 <- vegdist(t(otu_table_psSeaWater_hel1)) 
ps.disper.SeaWater1 <- betadisper(dist.helSeaWater1, metadf.SeaWater_hel1$Depth) 
permutest(ps.disper.SeaWater1)

####pairwise Adonis seawater_deep_shallow####
psSeaWater_hel2 <-subset_samples(psSeaWater_hel, Depth != "Middle")
psSeaWater_hel2 <- prune_taxa(taxa_sums(psSeaWater_hel2) > 0, psSeaWater_hel2)
otu_table_psSeaWater_hel2 <- as.data.frame(psSeaWater_hel2@otu_table)
metadf.SeaWater_hel2 <- data.frame(sample_data(psSeaWater_hel2))
otu_table_psSeaWater_hel2 <- as.data.frame(psSeaWater_hel2@otu_table)
adonis(distance(psSeaWater_hel2, method="bray") ~ Depth,
       data = metadf.SeaWater_hel2)
dist.helSeaWater2 <- vegdist(t(otu_table_psSeaWater_hel2)) 
ps.disper.SeaWater2 <- betadisper(dist.helSeaWater2, metadf.SeaWater_hel2$Depth) 
permutest(ps.disper.SeaWater2)

####pairwise Adonis seawater_middle_shallow####
psSeaWater_hel3 <-subset_samples(psSeaWater_hel, Depth != "Deep")
psSeaWater_hel3 <- prune_taxa(taxa_sums(psSeaWater_hel3) > 0, psSeaWater_hel3)
otu_table_psSeaWater_hel3 <- as.data.frame(psSeaWater_hel3@otu_table)
metadf.SeaWater_hel3 <- data.frame(sample_data(psSeaWater_hel3))
otu_table_psSeaWater_hel3 <- as.data.frame(psSeaWater_hel3@otu_table)
adonis(distance(psSeaWater_hel3, method="bray") ~ Depth,
       data = metadf.SeaWater_hel3)
dist.helSeaWater3 <- vegdist(t(otu_table_psSeaWater_hel3)) 
ps.disper.SeaWater3 <- betadisper(dist.helSeaWater3, metadf.SeaWater_hel3$Depth) 
permutest(ps.disper.SeaWater3)
