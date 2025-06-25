# Bacterial Community Analysis in Trade Wind Biomes

This project aims to compare bacterial communities in trade wind biomes at the Deep Chlorophyll Maximum (DCM) under varying salinities, using salinity as a proxy for pH. The analysis focuses on understanding how salinity influences bacterial diversity and community composition in marine environments.

## Project Structure

- **data/**: Contains all the necessary data files for the analysis.
  - **TO_OTU.txt**: Operational Taxonomic Unit (OTU) table with bacterial abundance across samples.
  - **TO_TAX.txt**: Taxonomy table providing classifications for the OTUs.
  - **META2.txt**: Metadata related to the samples, including environmental and sampling information.
  - **ALL_META.txt**: Additional metadata for broader context or supplementary information.

- **scripts/**: Contains the R Markdown document for analysis.
  - **analysis.Rmd**: R Markdown file with code for data loading, processing, statistical analysis, and visualization.

- **results/**: Documentation for the results generated from the analysis.
  - **README.md**: Explains the findings and how to interpret the results.

- **plots/**: Documentation for the plots generated during the analysis.
  - **README.md**: Details the visualizations and their significance.

## Objectives

1. To analyze the bacterial diversity in trade wind biomes at the DCM.
2. To investigate the relationship between salinity and bacterial community composition.
3. To visualize the findings through various statistical and graphical methods.

## Methodology

1. **Data Loading**: Import OTU, taxonomy, and metadata files.
2. **Data Processing**: Normalize data and filter for relevant samples.
3. **Statistical Analysis**: Perform alpha and beta diversity analyses, including PERMANOVA.
4. **Visualization**: Create plots to illustrate the findings, including ordination plots and bar plots of dominant phyla.

## Instructions

To run the analysis, ensure that all data files are placed in the `data/` directory. Open the `scripts/analysis.Rmd` file in RStudio or an appropriate R environment, and execute the code chunks sequentially. The results and plots will be generated and saved in their respective directories.

## Conclusion

This project provides insights into the dynamics of bacterial communities in marine ecosystems, highlighting the impact of environmental factors such as salinity on microbial diversity.