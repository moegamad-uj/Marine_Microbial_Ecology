# North Pacific Microbial Diversity Project

## Overview
This project aims to assess microbial diversity within the North Pacific Gyre, focusing on two sampling points that, despite their proximity, exhibit significant differences in microbial community composition. Through a series of analyses, we will explore environmental variables, community composition, and statistical significance of observed differences.

## Objectives
1. **Characterize Microbial Diversity**: Analyze the operational taxonomic units (OTUs) and their relative abundances to understand the microbial community structure at the two sampling points.
2. **Compare Environmental Conditions**: Investigate the environmental variables associated with each sampling point, including temperature, salinity, oxygen levels, and nutrient concentrations.
3. **Statistical Analysis**: Conduct statistical tests to determine the significance of differences in microbial diversity and environmental conditions between the two points.

## Data Sources
- **OTU Table**: `data/raw/TO_OTU.txt` - Contains abundance data for different microbial taxa across samples.
- **Taxonomy Information**: `data/raw/TO_TAX.txt` - Provides classification details for each OTU.
- **Sample Metadata**: `data/raw/META2.txt` and `data/raw/ALL_META.txt` - Include environmental conditions and sample identifiers.

## Analysis Workflow
1. **Data Loading and Preprocessing**: The main analysis is conducted in `analysis/NPO_samples.Rmd`, where data is loaded, processed, and visualized.
2. **Environmental Analysis**: The script `analysis/environmental_analysis.R` focuses on comparing environmental variables between the two sampling points.
3. **Community Composition Analysis**: The script `analysis/community_composition.R` analyzes the composition of microbial communities and visualizes taxonomic distributions.
4. **Statistical Testing**: The script `analysis/statistical_tests.R` performs statistical tests to assess the significance of observed differences.

## Results
- **Plots**: Visualizations of environmental comparisons, alpha diversity metrics, community composition, and geographic locations are stored in the `plots` directory.
- **Processed Results**: The results of the analysis are saved in `data/processed/npo_analysis_results.rds` for further exploration.

## Manuscript and Presentation
- A draft manuscript summarizing the findings is available in `manuscript/draft_manuscript.docx`.
- A presentation, including a speech draft and slides, is prepared in the `presentation` directory.

## References
All literature cited in the project can be found in `references/literature.bib`.

## Project Requirements
For a successful run of the project, please refer to `project_requirements.txt` for a list of dependencies and requirements.

## Conclusion
This project provides insights into the microbial diversity of the North Pacific Gyre, highlighting the importance of environmental factors and geographic context in shaping microbial communities.