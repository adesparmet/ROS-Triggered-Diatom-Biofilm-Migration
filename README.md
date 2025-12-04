# Data processing for the article "ROS-Triggered Diatom Biofilm Migration"

This repository contains the R scripts used to analyze the data from the article:  
"Reactive oxygen species drive downward vertical migration in diatom microphytobenthic biofilms as a strategy to cope with oxidative stress", Desparmet et al.

## Repository structure (SCRIPT Folder)

- `ROS_MIGRATION.Rproj`: RStudio project file.
- `1_Data_formatting`: all R scripts used for data processing.
- `2_Figures`: R script for figures.

### Raw data

The raw data used in this project are available on ZENODO:  
https://doi.org/10.5281/zenodo.15835853 ;
https://doi.org/10.5281/zenodo.15847025 (not required for scripts)

Please place the downloaded content of `data` (from zenodo) into the `/SCRIPT/1_Data_formatting/data` folder.

The raw metabarcoding data used in this project are available on ENA:
https://www.ebi.ac.uk/ena/browser/view/ERA33152396

#### How to use the scripts

The typical order of execution for the scripts in the `scripts/` folder is:

1. Run the five `_formatting.R` scripts to import, clean and transform raw data.
2. Run `Fig_article.R`: for statistical analyses and article figure generation.

All paths are relative to the project root.

##### Dependencies

This project uses the following main R packages:

-ade4
-ape
-BiocManager
-cowplot
-devtools
-dplyr
-factoextra
-FactoMineR
-ggplot2
-ggpubr
-ggrepel
-gt
-here
-hms
-iNEXT
-jsonlite
-lattice
-lubridate
-metagMisc
-microbiome
-patchwork
-permute
-phyloseq
-plotly
-readr
-readxl
-reshape2
-rnaturalearth
-sf
-stringr
-thejesus
-tidyr
-utils
-vegan
-webshot2

> Install them with `install.packages()`.

##  License

All code in this repository is licensed under the MIT License.  
All data are shared under the Creative Commons Attribution 4.0 International (CC BY 4.0) License.  

# Professional Network

LinkedIn: https://www.linkedin.com/in/alexandre-desparmet-bb5341171/
ResearchGate: https://www.researchgate.net/profile/Alexandre-Desparmet
