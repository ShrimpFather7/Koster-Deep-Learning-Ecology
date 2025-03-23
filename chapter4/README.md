# Chapter 4: GLMs and trait analyses  

This chapter contains the R script and data used for exploring abundance trends of taxa investigated by Nilsson et al. (2025) using output data from [Chapter 2](https://github.com/ShrimpFather7/Koster_Deep-Learning_Ecology/tree/main/chapter2 "Chapter 2 – Annotations to GBIF data"). Afterwards, occurrences can be intersected with temperature values from [Bio-ORACLE](https://doi.org/10.1111/geb.13813 "Bio-ORACLE v3.0 original publication") and filtered to estimate median temperature and temperature ranges of investigated taxa.

## **Citation**
Please use the following citation for use of any data or code provided in this repository:
- Zenodo blablabla Nilsson 2025

Example dataset is taken from [Koster Historical Biodiversity Assessment](https://doi.org/10.15468/rzhmef "GBIF – Koster Historical Biodiversity Assessment"). For any use of this data, please use the following citation:
- Nilsson CL, Anton V, Burman E, Germishuys J, Obst M. (2025). Koster historical biodiversity assessment. Version 1.7. Wildlife.ai. Occurrence dataset. https://ipt.gbif.org.nz/resource?r=koster_historical_assessment&v=1.7 https://doi.org/10.15468/rzhmef.

This script was originally used for data processing in:
- Nilsson et al. 2025

## **Contents**  
- `Depth_comparison_and_temp_assignments.R` – R-script that extracts depth ranges from investigated taxa in the original dataset and a GBIF occurrence download, and assigns investigated taxa with a median temperature based on Bio-ORACLE data and GBIF occurrence data.

- `koster_historical_biodiversity_assessment_occurrence.txt` – Tab-separated .txt file with Darwin Core formatted abundance data sorted by taxon and depth, utilized in data analyses by Nilsson et al. (2025).

- `chp_3_GBIF_Swedata.csv` – The same .csv-file with [GBIF occurrences](https://doi.org/10.15468/dl.rcne77 "GBIF Occurrence download (Sweden)") utilized by Nilsson et al. (2025) for comparison of their modeled results with depth distributions of the same species they investigated throughout Sweden. [See note 1](#note-1)

- `chp_3_GBIF_Globaldata.csv` – The same .csv-file with [GBIF occurrences](https://doi.org/10.15468/dl.azec6t "GBIF Occurrence download (Global)") utilized by Nilsson et al. (2025) for intersection of temperatures to worldwide occurrences between the years 2000-2019 of the same species they investigated. [See note 1](#note-1)

- `chp_3_publication_temp_occurrences.csv` – The filtered .csv-file utilized by Nilsson et al. (2025) with GBIF occurrences intersected with Bio-ORACLE temperature values. [See note 2](#note-2)

## **Dependencies**
- dplyr
- terra
- raster
- biooracler
- geosphere

## **How to Use**  
1. Load the script into RStudio.
2. Follow instructions to:
    - Extract depth ranges of taxa from data formatted as [Koster Historical Biodiversity Assessment](https://doi.org/10.15468/rzhmef "GBIF – Koster Historical Biodiversity Assessment") and from a GBIF occurrence download.
    - Assign temperature values to GBIF occurrences and filter to estimate taxa's temperature ranges.
    - Do both of the above :)

## **Output**
- `my_depth_ranges.csv` – Depth ranges of investigated taxa from data formatted like [Koster Historical Biodiversity Assessment](https://doi.org/10.15468/rzhmef "GBIF – Koster Historical Biodiversity Assessment"), i.e. from chapter 2 output `test_gbif.csv`.

- `gbif_depth_ranges.csv` – Depth ranges of investigated taxa from a GBIF occurrence download.
- `intersected_temps.csv` – GBIF occurrences intersected with temperature values from Bio-ORACLE.
- `filtered_temps.csv` – Filtered occurrences intersected with temperature.
- `central_ranges.csv` – Central 68% of temperature values and median temperatures assigned to all investigated taxa.

## **Contact**
Questions, bugs, and suggestions can be submitted via email to christian.nilsson7@outlook.com.

## **License**
- All scripts and data are available under the [MIT license](https://mit-license.org/), please cite any use or derivatives of this work.

## **Notes**
### <a id="note-1"></a>Note 1
Some taxa investigated by Nilsson et al. (2025) may have consisted of several species. Therefore, species of these taxa that have been suggested by the Swedish biodiversity portal [Artdatabanken](https://artfakta.se/ "Artfakta – the facts part about species") to reside at the study site of Nilsson et al. (2025) have been included in the occurrence downloads.

### <a id="note-2"></a>Note 2
As one filtering step is random, output may not be identical to `chp_3_publication_temp_occurrences.csv`. Therefore, this file has been included for reproducibility purposes.

