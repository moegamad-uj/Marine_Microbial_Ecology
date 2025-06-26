# Load necessary libraries
library(phyloseq)
library(vegan)
library(dplyr)

# Load processed NPO analysis results
npo_results <- readRDS("data/processed/npo_analysis_results.rds")

# Extract relevant data for statistical testing
npo_samples <- subset_samples(npo_results, Ocean_sea_regions == "(NPO)")
npo_meta <- data.frame(sample_data(npo_samples))

# Define a function to perform statistical tests
perform_statistical_tests <- function(data, variable) {
  # Check if there are at least two groups to compare
  if(length(unique(data[[variable]])) < 2) {
    stop("Not enough groups to perform statistical tests.")
  }
  
  # Perform ANOVA if there are more than two groups
  if(length(unique(data[[variable]])) > 2) {
    model <- aov(Mean_Oxygen_umolkg ~ get(variable), data = data)
    summary_model <- summary(model)
    return(list(type = "ANOVA", results = summary_model))
  } else {
    # Perform t-test if there are exactly two groups
    t_test_result <- t.test(Mean_Oxygen_umolkg ~ get(variable), data = data)
    return(list(type = "t-test", results = t_test_result))
  }
}

# Example: Perform statistical tests on oxygen levels by sample label
oxygen_test_results <- perform_statistical_tests(npo_meta, "Sample_label")

# Print results
print("Statistical Test Results for Oxygen Levels:")
print(oxygen_test_results)

# Save results to a file for later reference
saveRDS(oxygen_test_results, "data/processed/oxygen_test_results.rds")

# Additional tests can be added here for other variables (e.g., Nitrates)