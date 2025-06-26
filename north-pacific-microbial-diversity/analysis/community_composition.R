# Load Required Libraries
library(phyloseq)
library(ggplot2)
library(dplyr)
library(reshape2)

# Set Working Directory
setwd("C:/Users/Uzair/Documents/R/north-pacific-microbial-diversity")

# Load Data
otu_table <- read.table("data/raw/TO_OTU.txt", header=TRUE, row.names=1)
tax_table <- read.delim("data/raw/TO_TAX.txt", header=TRUE, row.names=1)
meta_data <- read.delim("data/raw/META2.txt", header=TRUE)

# Create Phyloseq Object
otu <- otu_table(otu_table, taxa_are_rows=TRUE)
tax <- tax_table(as.matrix(tax_table))
sample_data <- sample_data(meta_data)
physeq <- merge_phyloseq(otu, tax, sample_data)

# Subset for North Pacific Gyre Samples
npg_samples <- subset_samples(physeq, Ocean_sea_regions == "NPO")

# Analyze Community Composition
# Agglomerate at Class Level
npg_class <- tax_glom(npg_samples, "Class")

# Get Top 15 Classes
top_classes <- names(sort(taxa_sums(npg_class), decreasing=TRUE)[1:15])
npg_class_filtered <- prune_taxa(top_classes, npg_class)

# Plot Class Composition
class_plot <- ggplot(data=psmelt(npg_class_filtered), 
                     aes(x=Sample_label, y=Abundance, fill=Class)) +
  geom_bar(stat="identity", position="stack", width=0.8) +
  scale_y_continuous(labels=scales::percent_format()) +
  theme_bw() +
  labs(title="Bacterial Class Composition in North Pacific Gyre Samples",
       x="Sample",
       y="Relative Abundance") +
  theme(axis.text.x=element_text(angle=45, hjust=1))

# Save Plot
ggsave("plots/npo_class_composition.png", plot=class_plot, width=12, height=8, dpi=300)

# Print Summary of Findings
print("Community Composition Analysis Complete:")
print(paste("Total Samples Analyzed:", nsamples(npg_samples)))
print(paste("Top Classes Identified:", paste(top_classes, collapse=", ")))