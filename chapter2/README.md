# Chapter 2: Annotations to GBIF data  

This chapter contains the R script and data used for assigning depth values to annotations made by running a YOLOv8 model on footage through [SUBSIM Notebook "Evaluate machine learning models"](https://github.com/ocean-data-factory-sweden/kso/blob/dev/notebooks/classify/Upload_subjects_to_Zooniverse.ipynb "Evaluate models (GitHub)"). Afterward it collects average and max count of annotations per class per depth and reformats data according to [Darwin Core terminology](https://dwc.tdwg.org/list/ "Darwin Core List of Terms") for publication to [The Global Biodiversity Information Facility](https://www.gbif.org/ "GBIF.org").

## **Citation**
Please use the following citation for use of any data or code provided in this repository:
- Zenodo blablabla Nilsson 2025

Originally used for data processing in:
- Nilsson et al. 2025

## **Contents**  
- `Annotations_to_GBIF.R` – Script that assigns depth values to annotations, collects annotation stats per depth, and formats data for GBIF publication.
- `chp_2_annotations_transect_1.zip` .zip file containing example raw annotations data from YOLOv8 application to footage using SUBSIM.
- `chp_2_depth_df.csv` – A .csv-file including the depth values, corresponding frame numbers, and name of the files from which each frame was taken. [See note 1](#note-1)


## **Dependencies**
- zoo
- dplyr
- tidyr
- ggplot2

## **How to Use**  
1. Ensure all images have their file name and frame number in the file name (preferably separated by '_').
2. Load the script in into Jupyter Notebooks.
3. Follow instructions in the script to edit and then either crop and read or just read your images. [See note 2](#note-2)
4. Remove misread values with attached boolean terminology.

## **Output**
- `chp_2_depth_df.csv` – A .csv-file including the depth read, frame number, and name of the file which the frame was taken from.

## **License**
- All scripts and data are available under the [MIT license](https://mit-license.org/), please cite any use or derivatives of this work.

## **Contact**
Questions, bugs, and suggestions can be submitted via email to christian.nilsson7@outlook.com.

## **Notes**
### <a id="note-1"></a>Note 1
`depth_df` (generated in [Chapter 1](https://github.com/ShrimpFather7/Koster_Deep-Learning_Ecology/tree/main/chapter1 "Chapter 1 – Depth Extraction")) and `annotations` must stem from the same video files with the same file names.

### <a id="note-2"></a>Note 2
Cropping will overwrite all images in the folder, make sure you only have images you want to crop in the folder, and that you have a backup. 
