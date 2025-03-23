# Chapter 4: GLMs and trait analyses  

This chapter contains the R script and data used for exploring trends in community structure and ecosystem function using Generaized Linear Models (GLMs) on output data from [Chapter 2](https://github.com/ShrimpFather7/Koster_Deep-Learning_Ecology/tree/main/chapter2 "Chapter 2 – Annotations to GBIF data"). Community structure trends are explored through abundance trends of all investigated taxa, ecosystem function trends are explored by assigning trait categories (relating to size, lifestyle, and temperature preference) to taxa (in part using output from [Chapter 3](https://github.com/ShrimpFather7/Koster_Deep-Learning_Ecology/tree/main/chapter3 "Chapter 3 – Depth exploration and temp assignments")) and investigating abundance trends of these categories. Finally the impact of traits on taxon-specific abundance change coefficients will be investigated with Wilcoxon rank sum tests and a linear regression.

## **Citation**
Please use the following citation for use of any data or code provided in this repository:
- Zenodo blablabla Nilsson 2025

Example dataset is taken from [Koster Historical Biodiversity Assessment](https://doi.org/10.15468/rzhmef "GBIF – Koster Historical Biodiversity Assessment"). For any use of this data, please use the following citation:
- Nilsson CL, Anton V, Burman E, Germishuys J, Obst M. (2025). Koster historical biodiversity assessment. Version 1.7. Wildlife.ai. Occurrence dataset. https://ipt.gbif.org.nz/resource?r=koster_historical_assessment&v=1.7 https://doi.org/10.15468/rzhmef.

This script was originally used for data processing in:
- Nilsson et al. 2025

## **Contents**  
- `GLMs_and_trait_analyses.R` – R-script that investigates abundance trends of investigated taxa and trait categories, as well as the impact of traits on abundance trends on taxa.

- `koster_historical_biodiversity_assessment_occurrence.txt` – Tab-separated .txt file with Darwin Core formatted abundance data sorted by taxon and depth, utilized in data analyses by Nilsson et al. (2025).

- `chp_4_traits_df.csv` – A .csv-file with all investigated taxa and their assigned traits as assigned by Nilsson et al. (2025)

- `stats_output.txt` – Example output from Nilsson et al. (2025), copied from the console to a .txt file.

## **Dependencies**
- dplyr

## **How to Use**  
1. Load the script into RStudio.
2. Follow instructions to:
    - Investigate abundance trends of taxa from data formatted as [Koster Historical Biodiversity Assessment](https://doi.org/10.15468/rzhmef "GBIF – Koster Historical Biodiversity Assessment").
    - Assign trait categories to taxa and investigate abundance trends of trait categories.
    - Investigate the impact of trait categories on taxon-specific abundance trends with Wilcoxon rank sum tests and a linear regression.
3. Copy console output from GLMs and all other statistical tests to .txt files for later use/analysis.

## **Output**
As the statistical tests' output include many parameters that can be of interest, users are advised to copy console output of statistical tests to a .txt file for later use. `stats_output.txt` is provided as an example from results of Nilsson et al. (2025).

## **Contact**
Questions, bugs, and suggestions can be submitted via email to christian.nilsson7@outlook.com.

## **License**
- All scripts and data are available under the [MIT license](https://mit-license.org/), please cite any use or derivatives of this work.
