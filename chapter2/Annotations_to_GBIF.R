

# Hello! Welcome to this script that allows you to assign depths to annotations from the
# annotations.csv output file from SUBSIM (www.subsim.se) Jupyter notebooks "Evaluate models". 
# Adapted to depth dataframes as created in Chapter 1 - "Extract_Depths". Afterward the 
# script collects statistics on mean/max count per depth of each analyzed taxon, and finally 
# reformat to a Darwin-Core compatible .csv file fit for GBIF publication and OBIS 
# integration. Originally used by Nilsson et al. (2025) to make the Koster Historical
# Biodiversity Assessment https://doi.org/10.15468/rzhmef.


# Load the necessary packages.

library(dplyr)
library(zoo)
library(ggplot2)
library(tidyr)


# Load your annotations data

data_path <- "Path/to/chp_2_annotations_transect_1.csv"
data <- read.csv(data_path)

# Start with checking the filename column
data$filename
# Looks messy right?

# Clean the "filename" column to include only the filename and without path and filetype

# Try different character lengths to only include the video filename without path, frame 
# number, or file extension. i.e. from e.g. "Path/to/video/file_1337.txt" to "file".
first_char <- 34
last_char <- 41

# Check that you got it right. i.e. output from previous example would only give "file".
substr(data$filename, first_char, last_char) # If working with example data, you'll see
                                             # a series of "01747001" if correctly done.

# Apply filename cleanup
data$filename <- substr(data$filename, first_char, last_char)


# OPTIONAL - Filter out annotations by a confidence threshold

# Specify the threshold for confidence scores of model annotations
threshold <- 0.6

# Filter the data based on the threshold
data <- data %>% filter(conf > threshold)


# Create a dataframe including all annotated frames and the number of annotations for each 
# taxon for each frame
create_distribution <- function(df){
  df <- df %>%
    group_by(filename, class_id, frame_no) %>%
    summarize(occurrences = n()) %>%
    ungroup() %>%
    group_by(filename, class_id) %>%
    arrange(filename, class_id, frame_no)
  
  return(df)
}

# Apply it to your data
distribution <- create_distribution(data)


# Load your depth df.

depth_path <- "Path/to/chp_2_depth_df.csv"
depth <- read.csv(depth_path)%>%
  mutate(filename = sub("^'", "", filename)) # Remove the leading apostrophe from "filename"


# OPTIONAL (but usually necessary) - Run this section if your annotations.csv consists of
# multiple files where the frame number resets at the start of each new movie.

# Create a function that ensures that if you have multiple movies from a single dive
# where the frame number count resets that the frame numbers become consecutive.
# (provided your filenames are arranged in chronological order).

# Function that makes frame number consecutive in both your distribution df and
# your depth dataframes.

make_frames_consecutive <- function(df,depthframe) {
  
  filenames <- unique(df$filename)
  
  for(file in filenames) {
    
    print(paste0("The max frame number of file: ", file, " was ", max(df$frame_no[which(df$filename == file)])))
    
    df$frame_no[which(df$filename == file)] <- df$frame_no[which(df$filename == file)] + 
      min(depthframe$frame_no[which(depthframe$filename == file)]) # As the filename in your
                                                              # depth df should be the same 
                                                              # as in your distribution
                                                              # frames can be made consecutive
                                                              # by taking the minimum frame
                                                              # number from the same file in
                                                              # the depthframe (which is
                                                              # already consecutive)!
    
    print(paste0("It is now ", max(df$frame_no[which(df$filename == file)])))
    
  }
  
  return(df)
}

# Run it! (If you get -Inf instead of a revised frame number, something is wrong)
consecutive_distribution <- make_frames_consecutive(distribution,depth)


# Remove the first frames if there's any unsuitable habitat in the footage.

start_threshold <- 16905 # 0 if you want to keep all frames from start

    # Activate if you want to get the 500th frame of the last vid rather than the 35428th 
    # frame of the transect
end_threshold <- 500 + min(depth$frame_no[which(depth$filename==unique(depth$filename)[length(unique(depth$filename))])])

# Activate if you don't want to filter
#end_threshold <- max(consecutive_distribution$frame_no)

corrected_distribution <- consecutive_distribution %>% 
  filter(frame_no >= start_threshold) %>% 
  filter(frame_no <= end_threshold)


# Make a function that includes all occurrences, including absences of each taxon at each
# frame number.

expansivize <- function(df){
  
  all_combinations <- expand.grid(class_id = unique(df$class_id), frame_no = seq(start_threshold,end_threshold))
  all_combinations <- all_combinations[order(all_combinations$class_id,all_combinations$frame_no),]
  all_combinations$occurrences <- 0
  all_combinations$filename <- NA
  
  merged_df <- merge(all_combinations, df, by = c("frame_no", "class_id"), all.x = TRUE)
  merged_df <- merged_df %>%
    mutate(occurrences = ifelse(is.na(occurrences.y), occurrences.x, occurrences.y)) %>%
    dplyr::select(-occurrences.y, -occurrences.x) %>%
    mutate(filename = ifelse(is.na(filename.y), filename.x, filename.y)) %>%
    dplyr::select(-`filename.y`, -`filename.x`)
  
  return(merged_df)
  
}

# Run it!
expanded_distribution <- expansivize(corrected_distribution)


# Create a function that matches depth values to the original dataframe by frame number.
# The function uses linear interpolation for any framenumbers without depth values.

match_depths <- function(original_distribution,depthframe){
  
  # Join depthframe to original_distribution and drop filename from depthframe.
  original_distribution <- original_distribution %>%
    left_join(select(depthframe, -filename), by = "frame_no")
  
  # Interpolate missing depth values.
  original_distribution$depth_m <- na.approx(original_distribution$depth_m, rule = 2)
  
  # Round depth values
  original_distribution$depth_m <- round(original_distribution$depth_m)
  
  return(original_distribution)
  
}

# Run it!!
distribution_and_depth <- match_depths(expanded_distribution,depth)

# Check output! Filename is NA for all non-annotated taxa per frame
distribution_and_depth


# OPTIONAL - If you need to filter out any snippets of the video

file <- "01747002"
snippet_start <- 4900 #+ min(depth$Frame_no[which(depth$filename=='file')]) # Unhide the '+'
snippet_end <- 6230 #+ min(depth$Frame_no[which(depth$filename=='file')])  # part only if
                                                                          # you don't know
                                                                      # the consecutive frame
# Apply filters                                                       # number
distribution_and_depth <- distribution_and_depth %>% 
  filter(!between(frame_no, snippet_start, snippet_end))


# CHECKPOINT - Well done! You can now export this rawdata for e.g. combination with data
# from other transects for multi-transect analyses
write.csv(distribution_and_depth,
          "Path/to/where/you/want/the/rawdata.csv", row.names = FALSE)



# Now, let's prepare the data for GBIF/OBIS publication

# First let's summarize each depth for each taxon with statistics of mean/max count per
# depth.
summarize_data <- function(df){
  
  summarized_data <- df %>%
    group_by(class_id, depth_m) %>%
    summarise(meancount = mean(occurrences),
              maxcount = max(occurrences))
  
  return(summarized_data)
  
}

# Run it
summarized_distribution <- summarize_data(distribution_and_depth)


# OPTIONAL - Make histograms for distributions of your taxa

# Reshape data to long format for faceting
df_plot <- summarized_distribution %>%
  pivot_longer(cols = c(meancount, maxcount), 
               names_to = "count_type", 
               values_to = "count_value")%>%
  mutate(depth_m = factor(depth_m, levels = rev(unique(depth_m))))

# Plot histograms of count distribution across depth. Very simple visualization.
ggplot(df_plot, aes(x = factor(depth_m), y = count_value, fill = count_type)) +
  geom_col(position = "dodge", color = "black", width = 0.8) +  # Wider bars
  facet_wrap(~ class_id, scales = "free_y") +
  labs(x = "Depth (m)", y = "Count", fill = "Count Type", 
       title = "Histograms of Mean and Max Count Across Depth") +
  coord_flip() +  # Flipping the x-axis
  theme_minimal()


# Alright, back to business. Let's make that GBIF dataframe, renaming columns to Darwin-
# core terms.

# List all your taxa, must be in same order as class_id.
taxa <- c("Protanthea simplex", "Sabella pavonina", "Geodia barretti", 
          "Mycale lingua", "Polycarpa pomaria", "Serpulidae", 
          "Molgula", "Ascidia", "Ciona intestinalis",
          "Porania pulvillus","Acesta excavata", "Actiniidae", 
          "Munida","Phakellia ventilabrum", "Lithodes maja",
          "Caryophyllia smithii", "Alcyonium digitatum")

# List all taxonomic ranks, must be in same order as class_id.
ranks <- c("Species", "Species", "Species", "Species",
           "Species", "Family", "Genus", "Genus",
           "Species", "Species", "Species", "Family",
           "Genus", "Species", "Species", "Species",
           "Species")

# List all taxa again using underscores for spaces, must be in same order as class_id. 
# This is for the occurrence ID column.
ids <- c("Protanthea_simplex", "Sabella_pavonina", "Geodia_baretti",
         "Mycale_lingua", "Polycarpa_pomaria", "Serpulidae",
         "Molgula", "Ascidia", "Ciona_intestinalis",
         "Porania_pulvillus", "Acesta_excavata", "Actiniidae",
         "Munida", "Phakellia_ventilabrum", "Lithodes_maja",
         "Caryophyllia_smithii", "Alcyonium_digitatum")

# List all id numbers from the World Register of Marine Species (WoRMS), 
# must be in same order as class_id.
worms_ids <- c("urn:lsid:marinespecies.org:taxname:100913", "urn:lsid:marinespecies.org:taxname:130967",
               "urn:lsid:marinespecies.org:taxname:134023", "urn:lsid:marinespecies.org:taxname:133289",
               "urn:lsid:marinespecies.org:taxname:103909", "urn:lsid:marinespecies.org:taxname:988",
               "urn:lsid:marinespecies.org:taxname:103509", "urn:lsid:marinespecies.org:taxname:103483",
               "urn:lsid:marinespecies.org:taxname:103732", "urn:lsid:marinespecies.org:taxname:125166",
               "urn:lsid:marinespecies.org:taxname:140232", "urn:lsid:marinespecies.org:taxname:100653",
               "urn:lsid:marinespecies.org:taxname:106835", "urn:lsid:marinespecies.org:taxname:132511",
               "urn:lsid:marinespecies.org:taxname:107205", "urn:lsid:marinespecies.org:taxname:135144",
               "urn:lsid:marinespecies.org:taxname:125333")

# Run all the list functions before running this. This will create a GBIF/Darwin-core compatible
# dataframe that includes columns: occurrenceID, basisofRecord, scientificName, taxonRank,
# scientificNameID, verbatimDepth, individualCount, organismQuantity, organismQuantityType,
# eventDate, decimalLatitude, decimalLongitude, geodeticDatum, coordinateUncertaintyInMeters,
# datasetName, institutionCode, and occurrenceStatus.

# For details of these terms, check https://dwc.tdwg.org/list/.
create_gbif_df <- function(og_df, date, id_date, time, lat, lon, coordinatetype, coordinateuncertainty, datasetname, institutioncode, occurrencestatus){
  df <- data.frame(
    occurrenceID = og_df$class_id,
    basisOfRecord = rep("MachineObservation",length(og_df$class_id)),
    scientificName = og_df$class_id,
    taxonRank = og_df$class_id,
    scientificNameID = og_df$class_id,
    verbatimDepth = og_df$depth_m,
    individualCount = og_df$maxcount,
    organismQuantity = og_df$meancount,
    organismQuantityType = rep("Average count per frame",length(og_df$class_id)),
    eventDate = rep(date,length(og_df$class_id)),
    decimalLatitude = rep(lat,length(og_df$class_id)),
    decimalLongitude = rep(lon,length(og_df$class_id)),
    geodeticDatum = rep(coordinatetype,length(og_df$class_id)),
    coordinateUncertaintyInMeters = rep(coordinateuncertainty,length(og_df$class_id)),
    datasetName = rep(datasetname, length(og_df$class_id)),
    institutionCode = rep(institutioncode, length(og_df$class_id))
  )
  
  df$occurrenceStatus <- ifelse(df$individualCount==0, "absent", "present")
  
  for(i in 1:17){
    df$scientificName[which(df$scientificName==(i-1))] <- taxa[[i]]
  }
  
  for(i in 1:17){
    df$taxonRank[which(df$taxonRank==(i-1))] <- ranks[[i]]
  }
  
  for(i in 1:17){
    df$scientificNameID[which(df$scientificNameID==(i-1))] <- worms_ids[[i]]
  }
  
  for(i in 1:17){
    df$occurrenceID[which(df$occurrenceID==(i-1))] <- paste0("SUBSIM_Koster_", df$decimalLatitude[1], "_", df$decimalLongitude[1], "_", ids[i], "_")
  }
  
  for(i in 1:length(df$verbatimDepth)){
    df$occurrenceID[i] <- paste0(df$occurrenceID[i], df$verbatimDepth[i], "_", id_date, "_", time)
  }
  
  return(df)
}

# Parameter info: 

date = "2015-08-12" #date formatted as YYYY-MM-DD 

id_date = "120815" #date formatted as DDMMYY

time = "1100" #time formatted as HHMM (is completely arbitrary and is primarily used to separate
# observations from the same date/depth but different transects)

lat = 58.9 #decimal latitude

lon = 11.1 #decimal longitude

coordinatetype = "WGS84" #type of coordinate measurement

coordinateuncertainty = 3000 #coordinate uncertainty in meters

datasetname = "test_dataset" #name of your dataset

institutioncode = "CLN" #code of your institution

# RUN IT!!!
gbif_data <- create_gbif_df(summarized_distribution, date = date, id_date = id_date,
                            time = time, lat = lat, lon = lon, coordinatetype = coordinatetype, 
                            coordinateuncertainty = coordinateuncertainty, datasetname = datasetname,
                            institutioncode = institutioncode)

# Check it out! Looks alright?
gbif_data

# Beautifully done! Export your data and enjoy your soon to be published GBIF dataset
write.csv(gbif_data,"test_gbif.csv", row.names = FALSE, quote = FALSE)
