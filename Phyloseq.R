library(phyloseq)
library("tidyverse") 
library(plyr)
library(vegan)
library("ggplot2")
library("grid")

setwd("/Users/Emma/Documents/UCT\ /Honours\ module/honours\ module/Datasets")

#Import Taxonomy and OTU files: 
TAXfile="TO_TAX.txt"
OTUfile="TO_OTU.txt"

otu_TO= read.table(OTUfile, header=TRUE)
tax_TO= read.delim(TAXfile, header=TRUE)
tax_TO=column_to_rownames(tax_TO, "X")  
tax_TO<- as.matrix(tax_TO)

otu_TO = otu_table(otu_TO, taxa_are_rows = TRUE)
tax_TO= tax_table(tax_TO)

TO<- merge_phyloseq(otu_TO, tax_TO)
TO


###Metadata####
META2 <- read.delim("META2.txt", header=TRUE)
META <- read.delim("ALL_META.txt", header=TRUE)
META_TO <- merge(META2, META, by="PANGAEA")

META_TO=column_to_rownames(META_TO, "PANGAEA")  
MAP <- sample_data(META_TO)

TO1 <-merge_phyloseq(TO, MAP)
TO1

#Bacteria fraction:
TO.bacteria = subset_taxa(TO1, Kingdom=="Bacteria")
TO.bacteria

#Archaea fraction:
TO.archaea = subset_taxa(TO1, Kingdom=="Archaea")
TO.archaea

#Euk fraction:
TO.euk = subset_taxa(TO1, Kingdom=="Eukaryota")
TO.euk

#Coastal, Indian ocean data set
#All together:
TO.trades.bac = subset_samples(TO.bacteria, Marine_pelagic_biomes_Longhurst_2007=="Trades")
TO.trades.bac.io = subset_samples(TO.trades.bac, Ocean_sea_regions=="(IO)")
TO.trades.bac.sao = subset_samples(TO.trades.bac, Ocean_sea_regions=="(SAO)")
#TO.coastal = subset_samples(TO.coastal, Ocean_sea_regions=="(IO)")
#TO.coastal
TO.trades.bac
TO.trades.bac.io
TO.trades.bac.sao
SAO.IO <-merge_phyloseq(TO.trades.bac.io, TO.trades.bac.sao)
SAO.IO

#Antarctic versus SAO
TO.so.bac = subset_samples(TO.bacteria, Ocean_sea_regions=="(SO)")
TO.coastal.bac = subset_samples(TO.bacteria, Marine_pelagic_biomes_Longhurst_2007=="Coastal")
TO.coastal.bac.sao = subset_samples(TO.coastal.bac, Ocean_sea_regions=="(SAO)")
SAO.SO.coastal<-merge_phyloseq(TO.so.bac, TO.coastal.bac.sao)
SAO.SO.coastal = subset_samples(SAO.SO.coastal, Sample_label != "TARA_066_DCM_0.22-3" & Sample_label != "TARA_066_SRF_0.22-3" & Sample_label !="TARA_067_SRF_0.22-3")
SAO.SO.coastal
TO.coastal.bac.sao
sample_data(SAO.SO.coastal)
sample_data(TO.bacteria)

#IO coastal versus IO trades
TO.coastal.bac = subset_samples(TO.bacteria, Marine_pelagic_biomes_Longhurst_2007=="Coastal")
TO.coastal.io.bac = subset_samples(TO.coastal.bac, Ocean_sea_regions=="(IO)")
TO.trades.bac.io
TO.coastal.io.bac 
TO.IO <-merge_phyloseq(TO.trades.bac.io, TO.coastal.io.bac)
TO.IO

#Euks:
TO.euk.coastal = subset_samples(TO.euk, Marine_pelagic_biomes_Longhurst_2007=="Coastal")
TO.euk.coastal = subset_samples(TO.euk.coastal, Ocean_sea_regions=="(IO)")
TO.euk.coastal

#Bacteria:
TO.bac.coastal = subset_samples(TO.bacteria, Marine_pelagic_biomes_Longhurst_2007=="Coastal")
TO.bac.coastal = subset_samples(TO.bac.coastal, Ocean_sea_regions=="(IO)")
TO.bac.coastal

#Archaea:
TO.arch.coastal = subset_samples(TO.archaea, Marine_pelagic_biomes_Longhurst_2007=="Coastal")
TO.arch.coastal = subset_samples(TO.arch.coastal, Ocean_sea_regions=="(IO)")
TO.arch.coastal

#relative abundance:
SAO.SO.coastal.rel = transform_sample_counts(SAO.SO.coastal, function(x) x / sum(x) )
SAO.SO.coastal.rel

#archaea:
#TO.arch.coastal.rel = transform_sample_counts(TO.arch.coastal, function(x) x / sum(x) )
#TO.arch.coastal.rel
#bacteria:
#TO.bac.coastal.rel = transform_sample_counts(TO.bac.coastal, function(x) x / sum(x) )
#TO.bac.coastal.rel
#euks:
#TO.euk.coastal.rel = transform_sample_counts(TO.euk.coastal, function(x) x / sum(x) )
#TO.euk.coastal.rel

# Assuming your metadata has a numeric "Latitude" column
TO_polar <- subset_samples(TO1, Latitude >  66 | Latitude < -60)
# drop any taxa that have zero counts in that subset
TO_polar <- prune_taxa(taxa_sums(TO_polar) > 0, ps_polar)



#Westerlies North Atlantic
#all together:
TO.Westerlies = subset_samples(TO1, Marine_pelagic_biomes_Longhurst_2007=="Westerlies")
TO.Westerlies = subset_samples(TO.Westerlies, Ocean_sea_regions=="(NAO)")
TO.Westerlies

#Euks:
TO.euk.Westerlies = subset_samples(TO.euk, Marine_pelagic_biomes_Longhurst_2007=="Westerlies")
TO.euk.Westerlies = subset_samples(TO.euk.Westerlies, Ocean_sea_regions=="(NAO)")
TO.euk.Westerlies 

#Bacteria:
TO.bac.Westerlies = subset_samples(TO.bacteria, Marine_pelagic_biomes_Longhurst_2007=="Westerlies")
TO.bac.Westerlies = subset_samples(TO.bac.Westerlies, Ocean_sea_regions=="(NAO)")
TO.bac.Westerlies

#Archaea:
TO.arch.Westerlies = subset_samples(TO.archaea, Marine_pelagic_biomes_Longhurst_2007=="Westerlies")
TO.arch.Westerlies = subset_samples(TO.arch.Westerlies, Ocean_sea_regions=="(NAO)")
TO.arch.Westerlies

#Relative abundance:
#archaea:
TO.arch.westerlies.rel = transform_sample_counts(TO.arch.Westerlies, function(x) x / sum(x) )
TO.arch.westerlies.rel
#bacteria:
TO.bac.westerlies.rel = transform_sample_counts(TO.bac.Westerlies, function(x) x / sum(x) )
TO.bac.westerlies.rel
#euk:
TO.euk.westerlies.rel = transform_sample_counts(TO.euk.Westerlies, function(x) x / sum(x) )
TO.euk.westerlies.rel

#Trades Indian Ocean:
#all together:
TO.Trades = subset_samples(TO1, Marine_pelagic_biomes_Longhurst_2007=="Trades")
TO.Trades = subset_samples(TO.Trades, Ocean_sea_regions=="(IO)")
TO.Trades

#Euks:
TO.euk.Trades = subset_samples(TO.euk, Marine_pelagic_biomes_Longhurst_2007=="Trades")
TO.euk.Trades = subset_samples(TO.euk.Trades, Ocean_sea_regions=="(IO)")
TO.euk.Trades

#Bacteria:
TO.bac.Trades = subset_samples(TO.bacteria, Marine_pelagic_biomes_Longhurst_2007=="Trades")
TO.bac.Trades = subset_samples(TO.bac.Trades, Ocean_sea_regions=="(IO)")
TO.bac.Trades

#Archaea:
TO.arch.Trades = subset_samples(TO.archaea, Marine_pelagic_biomes_Longhurst_2007=="Trades")
TO.arch.Trades = subset_samples(TO.arch.Trades, Ocean_sea_regions=="(IO)")
TO.arch.Trades

#relative abundance:
#TO.bac.Trades.rel = transform_sample_counts(TO.bac.Trades, function(x) x / sum(x) )
#TO.coastal.rel = transform_sample_counts(TO.coastal, function(x) x / sum(x) )

#only OTUs with a mean greater than 10^-5 are kept (this will not work for Archaea)
SAO.SO.coastal.rel.mean= filter_taxa(SAO.SO.coastal.rel, function(x) mean(x) > 0, TRUE)
SAO.SO.coastal.rel.mean


#subset by taxa
Coastal_Rhizaria <- subset_taxa(TO.coastal.rel.mean, Phylum=="Rhizaria")
Coastal_Rhizaria

Westerlies_Archaepl <- subset_taxa(TO.westerlies.euk.rel.mean, Phylum=="Archaeplastida")
Westerlies_Archaepl

SO_CO_Flavo <- subset_taxa(SAO.SO.coastal.rel.mean, Class=="Flavobacteriia")
SO_CO_Flavo

# Colour scheme -----------------------------------------------------------
library(RColorBrewer)
library("viridis")
library(pals)

nb.cols <- 39
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)

# Barplots ----------------------------------------------------------------
#Simple bar plot example:
A <- plot_bar(SAO.SO.coastal.rel.mean, "Sample_label", fill="Phylum") 
A + scale_fill_manual(values=unname(polychrome()))

#OR:
title = "plot_bar; Antarctic versus SAO coastal"
plot_bar(SAO.SO.coastal.rel.mean, "Sample_label", fill="Class", title=title)

#barplot with facets:
plot_bar(SAO.SO.coastal.rel.mean, "Sample_label", fill="Class", facet_grid=~Marine_pelagic_biomes_Longhurst_2007) + scale_fill_manual(values=unname(polychrome()))
plot_bar(SAO.SO.coastal.rel.mean, "Sample_label", fill="Class", facet_grid=~Marine_pelagic_biomes_Longhurst_2007) + scale_fill_manual(values= mycolors)
plot_bar(SAO.SO.coastal.rel.mean, "Sample_label", fill="Class", facet_grid=~Marine_pelagic_biomes_Longhurst_2007) + scale_fill_viridis_d()

# remove those spesky black lines:
phylumGlommed = tax_glom(SAO.SO.coastal.rel.mean, "Phylum")
ClassGlommed = tax_glom(SAO.SO.coastal.rel.mean, "Class")
plot_bar(phylumGlommed, "Sample_label", fill="Phylum") + scale_fill_manual(values=unname(polychrome()))
plot_bar(ClassGlommed, "Sample_label", fill="Class", facet_grid=~Marine_pelagic_biomes_Longhurst_2007) + scale_fill_manual(values= mycolors)
# alpha diversity ---------------------------------------------------------
#You must run this on your raw data, NOT on the relative abundance data!

#example: IO coastal versus trades
plot_richness(TO.IO)

#Chao and Shannon only:
plot_richness(SAO.SO.coastal, "Sample_label",measures=c("Chao1", "Shannon"))

#OR by a catagorical variable. so first we need to know what is available:
sample_variables(SAO.SO.coastal)

#then pick one. lets say sampling depth:
plot_richness(TO.IO, x="Sampling_depth", measures=c("Chao1", "Shannon"))

#or latitude:
plot_richness(TO.IO, x="Longitude", measures=c("Chao1", "Shannon"))

#or colour by type:
plot_richness(SAO.SO.coastal, x="Latitude", color="Mean_Temperature_degC", measures=c("Chao1", "Shannon"))

# Ordination --------------------------------------------------------------
#Remember to use your normalised data object here!

#ordination:
GP.ord <- ordinate(SAO.SO.coastal.rel.mean, "PCoA", "bray")
p1 = plot_ordination(SAO.SO.coastal.rel.mean, GP.ord, type="taxa", color="Class", title="taxa")
print(p1)
p1 + facet_wrap(~Phylum, 2)

GP.ord <- ordinate(SAO.SO.coastal.rel.mean, "PCoA", "bray")
p4 = plot_ordination(SAO.SO.coastal.rel.mean, GP.ord, type="split", color="Phylum", label="Sample_label", title="split") 
p4



# Ordination with vector overlay ------------------------------------------
# This ordination shows how much each sample differs from the other composition wise. When we superimpose
#vectors on it this will tell us how these metadata are positively (pointing towards the sample) or negatively 
# (pointing directly away from a sample dot) influencing those samples (dots)

psv = data.frame(sample_data(SAO.SO.coastal.rel.mean)) #change TO.trades to your region of interest

# convert the sample_data() within a phyloseq object to a vegan compatible data object
pssd2veg <- function(physeq) {
  sd <- sample_data(physeq)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a phyloseq object to a vegan compatible data object
psotu2veg <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

ps.sd.veg <-  pssd2veg(SAO.SO.coastal.rel.mean) #remember to change TO.trades to your region of interest here
ps.otu.veg <-  psotu2veg(SAO.SO.coastal.rel.mean) # Same here

vare.dis <- vegdist(ps.otu.veg, method="bray")
library(MASS)
vare.mds <- metaMDS(ps.otu.veg, trace = FALSE)
vare.mds

data.scores = as.data.frame(scores(vare.mds))

ef <- envfit(vare.mds, psv, permu = 999, na.rm=TRUE)
ef
plot(vare.mds, display = "sites")
plot(ef, p.max = 0.001)
ordiplot (vare.mds, display = 'sites', type = 'n')
orditorp (vare.mds, display = 'sites')

#colour sites by a variable:
ORD = plot_ordination(SAO.SO.coastal.rel.mean, GP.ord, type="samples", color="EnvironmentalFeature") 
ORD1 = plot_ordination(SAO.SO.coastal.rel.mean, GP.ord, type="samples", color="Mean_Depth_Max.O2_m") 

plot(ORD1)
##Scatterplots - metadata: --------------------------

#convert your OTU table to a data frame:
OTU_Arch = data.frame(otu_table(Westerlies_Archaepl))
TAX_Arch = data.frame(tax_table(Westerlies_Archaepl))

#convert your Metadata table into a dataframe:
META_Arch = data.frame(sample_data(Westerlies_Archaepl))


ggplot(META_trades, aes(x=Mean_Depth_Nitracline_m, y=miTAG.SILVA.Shannon, color=Marine_pelagic_biomes)) +
  geom_point()

#How to plot OTU to a metadata:
#Calculate a total row to your OTU table:
Arch_sum <- colSums(OTU_Arch)

#Now add this line to your META file:
META_Arch1 <- rbind(t(META_Arch), Arch_sum)

META_Arch1 <- t(META_Arch1)
  
META_Arch1 <- as.data.frame(META_Arch1)
sapply(META_Arch1, mode)
META_Arch1 <- transform(META_Arch1, Arch_sum = as.numeric(Arch_sum))
META_Arch1 <- transform(META_Arch1, FC._autotrophs_cells.mL = as.numeric(FC._autotrophs_cells.mL))

#plot total Arch OTU versus Fluorescence using ggplot:
ggplot(META_Arch1, aes(x=Mean_Salinity_PSU, y=Arch_sum, color=FC._autotrophs_cells.mL)) +
  geom_point()


