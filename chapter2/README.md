# Chapter 2: Annotations to GBIF data  

This chapter contains the R script and data used for assigning depth values to annotations made by running a YOLOv8 model on footage through [SUBSIM Notebook "Evaluate machine learning models"](https://github.com/ocean-data-factory-sweden/kso/blob/dev/notebooks/classify/Upload_subjects_to_Zooniverse.ipynb "Evaluate models (GitHub)"). Afterward it collects average and max count of annotations per class per depth and reformats data according to [Darwin Core terminology](https://dwc.tdwg.org/list/ "Darwin Core List of Terms") for publication to [The Global Biodiversity Information Facility](https://www.gbif.org/ "GBIF.org").

## **Citation**
Please use the following citation for use of any data or code provided in this repository:
- Zenodo blablabla Nilsson 2025

Originally used for data processing in:
- Nilsson et al. 2025

## **Contents**  
- `Annotations_to_GBIF.R` – R-script that assigns depth values to annotations, collects annotation stats per depth, and formats data for GBIF publication.
- `chp_2_annotations_transect_1.zip` .zip file containing example raw annotations data from YOLOv8 application to footage using SUBSIM.
- `chp_2_depth_df.csv` – A .csv-file including the depth values, corresponding frame numbers, and name of the files from which each frame was taken. [See note 1](#note-1)


## **Dependencies**
- zoo
- dplyr
- tidyr
- ggplot2

## **How to Use**  
1. Load the script into R.
2. Follow instructions until the desired level of data is achieved, i.e:
    - Expanded raw data containing the number of annotations (including absences) of each class per frame).
    - Summarized data of max/mean number of annotations of each class per depth unit, using Darwin Core terminology formatted for GBIF publication.
4. When applying Darwin Core terminology, make sure to check the [term list](https://dwc.tdwg.org/list/ "Darwin Core List of Terms") and use the correct format for values (examples provided in code).

## **Output**
- `rawdata.csv` – Expanded raw data containing the number of annotations (including absences) of each class per frame).
- `test_gbif.csv` – Summarized data of max/mean number of annotations of each class per depth unit, using Darwin Core terminology formatted. This data will be formatted in the same way as [Koster Historical Biodiversity Assessment](https://doi.org/10.15468/rzhmef "GBIF – Koster Historical Biodiversity Assessment") data which is used in chapters 3-4 (but in tab-separated .txt format).

## **Contact**
Questions, bugs, and suggestions can be submitted via email to christian.nilsson7@outlook.com.

## **License**
- All scripts and data are available under the [MIT license](https://mit-license.org/), please cite any use or derivatives of this work.

## **Notes**
### <a id="note-1"></a>Note 1
`depth_df` (generated in [Chapter 1](https://github.com/ShrimpFather7/Koster_Deep-Learning_Ecology/tree/main/chapter1 "Chapter 1 – Depth Extraction")) and `annotations` must stem from the same video files with the same file names.
