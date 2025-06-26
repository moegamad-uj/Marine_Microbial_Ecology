# Load libraries
library(phyloseq)
library(tidyverse)
library(vegan)
library(ggplot2)
library(grid)
library(RColorBrewer)
library(viridis)
library(pals)
library(MASS)

# Set working directory
setwd("C:/Users/Uzair/Documents/R/Marine_Microbial_Ecology")

# Import OTU and Taxonomy tables
otu_TO <- read.table("TO_OTU.txt", header=TRUE)
tax_TO <- read.delim("TO_TAX.txt", header=TRUE)
tax_TO <- column_to_rownames(tax_TO, "X")
tax_TO <- as.matrix(tax_TO)

# Create Phyloseq OTU and Taxonomy objects
otu_TO <- otu_table(otu_TO, taxa_are_rows = TRUE)
tax_TO <- tax_table(tax_TO)

# Merge into single Phyloseq object
TO <- merge_phyloseq(otu_TO, tax_TO)

# Load metadata
META2 <- read.delim("META2.txt", header=TRUE)
META <- read.delim("ALL_META.txt", header=TRUE)
META_TO <- merge(META2, META, by="PANGAEA")
META_TO <- column_to_rownames(META_TO, "PANGAEA")
MAP <- sample_data(META_TO)

# Add metadata to Phyloseq object
TO1 <- merge_phyloseq(TO, MAP)

# Subset bacteria only
TO.bacteria <- subset_taxa(TO1, Kingdom == "Bacteria")

# Antarctic vs South Atlantic Ocean coastal comparison
TO.so.bac <- subset_samples(TO.bacteria, Ocean_sea_regions == "(SO)")
TO.coastal.bac <- subset_samples(TO.bacteria, Marine_pelagic_biomes_Longhurst_2007 == "Coastal")
TO.coastal.bac.sao <- subset_samples(TO.coastal.bac, Ocean_sea_regions == "(SAO)")
SAO.SO.coastal <- merge_phyloseq(TO.so.bac, TO.coastal.bac.sao)
SAO.SO.coastal <- subset_samples(SAO.SO.coastal, 
                                 !(Sample_label %in% c("TARA_066_DCM_0.22-3", "TARA_066_SRF_0.22-3", "TARA_067_SRF_0.22-3")))

# Normalise data to relative abundance
SAO.SO.coastal.rel <- transform_sample_counts(SAO.SO.coastal, function(x) x / sum(x))
SAO.SO.coastal.rel.mean <- filter_taxa(SAO.SO.coastal.rel, function(x) mean(x) > 0, TRUE)

# Alpha diversity (on raw data)
plot_richness(SAO.SO.coastal, x = "Latitude", color = "Mean_Temperature_degC", measures = c("Shannon", "Chao1"))

# Ordination
GP.ord <- ordinate(SAO.SO.coastal.rel.mean, method = "PCoA", distance = "bray")
plot_ordination(SAO.SO.coastal.rel.mean, GP.ord, type = "samples", color = "Mean_Temperature_degC")

# Environmental fitting for ordination
psv <- data.frame(sample_data(SAO.SO.coastal.rel.mean))
psotu2veg <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}
ps.otu.veg <- psotu2veg(SAO.SO.coastal.rel.mean)
vare.mds <- metaMDS(ps.otu.veg, trace = FALSE)
ef <- envfit(vare.mds, psv, permu = 999, na.rm = TRUE)
plot(vare.mds, display = "sites")
plot(ef, p.max = 0.01)

# Barplot of dominant classes
ClassGlommed <- tax_glom(SAO.SO.coastal.rel.mean, "Class")
plot_bar(ClassGlommed, "Sample_label", fill = "Class") + 
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Set1"))(39))