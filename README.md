# MY PAPER GitHub repository

## Koster-Deep-Learning-Ecology

This is the official GitHub repository for the scientific publication [MY PAPER](link.com "MY PAPER"), detailing the data, processing workflows, and statistical analyses utilized to generate results shown in this article. While the primary objective of this repository is to enhance the reproducibility of the study, code here can be highly useful in analyses of data from the Swedish platform for subsea image analysis ([SUBSIM](https://subsim.se/)) and may also be used for unrelated purposes as well.

## Description

MY PAPER trained a [YOLOv8 object-detection model](https://doi.org/10.5281/zenodo.13589902) to annotate prominent benthic invertebrate fauna from ROV footage spanning across 26 years (1997-2023) along a depth gradient in Kosterhavets national park (situated in the northernmost part of Sweden's west coast). 

Analysis of the [data](https://doi.org/10.15468/rzhmef) produced by the object-detection model highlighted changes in community structure along the depth gradient as well as significant temporal trends in community structure and ecosystem function. Temporal trends in community structure were primarily driven by size and temperature preference of the investigated taxa, with large, cold-loving taxa declining. However, the area's designation as a national park appears to benefit the remaining community.

## Repository Structure
The repository is divided into four chapters, following the order of methodology utilized by Nilsson et al. (2025). Each chapter contains a script to carry out the chapter's tasks, example data, and instructions for use (in a README file and the script).

```
/ (root)
|-- /chapter1 # Depth Extraction 
|-- /chapter2 # Annotations to GBIF data 
|-- /chapter3 # Depth exploration and temperature assignments 
|-- /chapter4 # GLMs and trait analyses
|-- README.md # Main repository README 
|-- LICENSE # License file
```

## How to Use

The primary idea of the repository is to follow scripts 1-4 using example data to follow data processing/analysis steps of Nilsson et al. (2025), but users may use any chapter irregardless of order to perform their own analyses if deemed relevant. Chapters include:

1. Using OCR software to extract depth values from a folder of image frames taken from ROV footage (with overlay visible).
   
2. Tying YOLOv8 annotations to depth values from chapter 1, gather abundance data for each investigated taxon and formatting data for GBIF publication.
   
3. Extracting depth ranges of investigated taxa from data generated in chapter 2 and a GBIF occurrence download for comparison. Finding temperature preferences of investigated taxa using GBIF occurrences and seafloor temperature data packages.
   
4. Investigating trends in community structure and ecosystem functioning using data generated in chapter 2 and temperature preferences from chapter 3. Investigating the impact of traits on abundance trends of individual taxa.

For further details on what is done in the different chapters, check chapter READMEs and scripts.

## Dependencies

The original workflow was executed using the programming languages Python (Jupyter Notebooks) and R (RStudio).

Python packages required:
- os
- csv
- PIL
- numpy
- pandas
- easyocr
- matplotlib.pyplot

R packages required:
- zoo
- dplyr
- tidyr
- ggplot2
- terra
- raster
- biooracler
- geosphere

For more information on chapter-specific dependencies, check chapter READMEs.

## Citation

Please use the following citation for use of any data or code provided in this repository:
- Zenodo blablabla Nilsson 2025

Associated article citation:
- Nilsson CL, Faurby S, Burman E, Germishuys J, Obst M. (2025). Title. Journal blabla

Associated dataset citation:
- Nilsson CL, Anton V, Burman E, Germishuys J, Obst M. (2025). Koster historical biodiversity assessment. Version 1.7. Wildlife.ai. Occurrence dataset. https://ipt.gbif.org.nz/resource?r=koster_historical_assessment&v=1.7 https://doi.org/10.15468/rzhmef.

Associated YOLOv8 model citation:
- Nilsson, CL, Germishuys, J, Burman, E, Anton, V, White, J, Obst, M (2024). Koster historical invertebrate model - SUBSIM 17tx. Zenodo. https://doi.org/10.5281/zenodo.13589902

  

## **Contact**
Questions, bugs, and suggestions can be submitted via email to christian.nilsson7@outlook.com.


## **License**
All scripts and data are available under the [MIT license](https://mit-license.org/), please cite any use or derivatives of this work.


