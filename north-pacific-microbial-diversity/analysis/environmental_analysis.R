# Load necessary libraries
library(phyloseq)
library(tidyverse)
library(vegan)
library(ggplot2)

# Set working directory to the project folder
setwd("C:/Users/Uzair/Documents/R/north-pacific-microbial-diversity")

# Load processed NPO analysis results
npo_data <- readRDS("data/processed/npo_analysis_results.rds")

# Extract environmental metadata
npo_meta <- data.frame(sample_data(npo_data))

# Define environmental variables of interest
env_vars <- c("Mean_Lat", "Mean_Long", "Mean_Depth_m", "Mean_Temperature_degC",
              "Mean_Salinity_PSU", "Mean_Oxygen_umolkg", "Mean_Nitrates_umolL")

# Subset the metadata for the environmental variables
npo_env <- npo_meta[, env_vars]
rownames(npo_env) <- npo_meta$Sample_label

# Print environmental conditions
print("Environmental Conditions for NPO Samples:")
print(npo_env)

# Calculate summary statistics for environmental variables
env_summary <- npo_env %>%
  summarise(across(everything(), list(mean = mean, sd = sd), na.rm = TRUE))

print("Summary Statistics for Environmental Variables:")
print(env_summary)

# Visualize environmental differences between the two sampling points
npo_env_long <- pivot_longer(npo_env, cols = everything(), names_to = "Variable", values_to = "Value")

ggplot(npo_env_long, aes(x = Sample_label, y = Value, fill = Sample_label)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Variable, scales = "free") +
  theme_bw() +
  labs(title = "Environmental Conditions Comparison: NPO Samples",
       x = "Sample",
       y = "Value") +
  scale_fill_brewer(palette = "Set2")

# Save the environmental comparison plot
ggsave("plots/npo_environmental_comparison.png", width = 12, height = 8, dpi = 300)

# Conduct correlation analysis between environmental variables
correlation_matrix <- cor(npo_env, use = "complete.obs")
print("Correlation Matrix of Environmental Variables:")
print(correlation_matrix)

# Visualize the correlation matrix
library(corrplot)
corrplot(correlation_matrix, method = "circle", type = "upper", tl.col = "black", tl.srt = 45)

# Save the correlation plot
ggsave("plots/npo_environmental_correlation.png", width = 10, height = 8, dpi = 300)