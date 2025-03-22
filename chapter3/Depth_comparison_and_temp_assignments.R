
# Hello! Welcome to this script that allows you to extract descriptive statistics on depth 
# distribution from the study's main dataset "Koster Historical Biodiversity Assessment" 
# (https://doi.org/10.15468/rzhmef) for later comparison with GBIF download of the same 
# taxa* throughout Sweden (https://doi.org/10.15468/dl.rcne77). Later on you can utilize a 
# GBIF download of the same taxa worldwide (https://doi.org/10.15468/dl.azec6t) to assign 
# different occurrences with temperature values from Bio-ORACLE (Assis et al., 2024) 
# https://doi.org/10.1111/geb.13813. From this, median temperature can be assigned to each 
# taxon for later traits analysis.


# *Species utilized to represent multispecies taxa were those likely to occur at the site
# based on occurrences and depth ranges from the Swedish biodiversity portal, Artdatabanken 
# (https://artfakta.se/).

# Load necessary packages.

library(dplyr)
library(biooracler)
library(raster)
library(terra)
library(geosphere)

# Load your main data

my_data_path <- "Path/to/koster_historical_biodiversity_assessment_occurrence.txt"
my_data <- read.delim(my_data_path, sep = "\t")


# Specify the preferred percentiles of distribution here. Default is 16th and 84th 
# percentile representing one degree of standard deviation from a Gaussian distribution.

percentiles <- c(0.16, 0.84)


# This script calculates the mean abundance (organismQuantity) across all years at each
# depth for each taxon, then it calculates the cumulative percentage of annotations at 
# each depth for each taxon, from which it extracts the depths where the previously
# specified percentiles and medians of annotations lie.

extract_my_depth_info <- function(df){
  
  df_with_meanQ <- my_data %>%
    group_by(scientificName, verbatimDepth) %>%
    summarize(meanquantity = mean(organismQuantity))
  
  depth_ranges <- df_with_meanQ %>%
    group_by(scientificName) %>%
    mutate(
      cumulative_percent = cumsum(meanquantity) / sum(meanquantity)  # Calculate cumulative distribution of meanquantity
    ) %>%
    summarize(
      lower_depth = verbatimDepth[which.min(abs(cumulative_percent - percentiles[1]))],
      median_depth = verbatimDepth[which.min(abs(cumulative_percent - 0.5))],
      upper_depth = verbatimDepth[which.min(abs(cumulative_percent - percentiles[2]))])
  
  return(depth_ranges)
  
}

# Run it!
my_depth_ranges <- extract_my_depth_info(my_data)

# Check it!
my_depth_ranges

# Write it!
write.csv(my_depth_ranges, "Path/to/where/I/want/it.csv", row.names = FALSE, quote = FALSE)


# Load the GBIF data.

swebif_data_path <- "Path/to/chp_3_GBIF_Swedata.csv"
swebif_data <- read.csv(swebif_data_path, sep = "\t")


# OPTIONAL here you can filter out datasets by their datasetKey, as the gbif_data includes
# data from Nilsson et al. (2024), the same as used for "my_data", this code will filter out
# that data.

swebif_data <- swebif_data[swebif_data$datasetKey!="51d0bd32-e215-45ea-a04d-47a474336125",]

# OPTIONAL If you're working with multiple taxonomic levels, you can add a "taxon" column
# here where you can map names of the different taxonomic levels with their respective names.
# e.g. if the species is Protanthea simplex its taxon will be Protanthea simplex, if the
# family is Serpulidae its taxon will be Serpulidae, etc. With data from Nilsson et al.
# scientificName is sufficient. But with gbif data you will have to specify taxonomic
# rank.

swebif_data <- swebif_data %>%
  mutate(taxon = case_when(
    species == "Protanthea simplex" ~ "Protanthea simplex",
    species == "Sabella pavonina" ~ "Sabella pavonina",
    species == "Geodia barretti" ~ "Geodia barretti",
    species == "Mycale lingua" ~ "Mycale lingua",
    species == "Polycarpa pomaria" ~ "Polycarpa pomaria",
    family == "Serpulidae" ~ "Serpulidae",
    genus == "Molgula" ~ "Molgula",
    genus == "Ascidia" ~ "Ascidia",
    species == "Ciona intestinalis" ~ "Ciona intestinalis",
    species == "Porania pulvillus" ~ "Porania pulvillus",
    species == "Acesta excavata" ~ "Acesta excavata",
    family == "Actiniidae" ~ "Actiniidae",
    genus == "Munida" ~ "Munida",
    species == "Phakellia ventilabrum" ~ "Phakellia ventilabrum",
    species == "Lithodes maja" ~ "Lithodes maja",
    species == "Caryophyllia smithii" ~ "Caryophyllia smithii",
    species == "Alcyonium digitatum" ~ "Alcyonium digitatum",
    TRUE ~ NA_character_  # Default if no condition is met
  ))


# Some GBIF datapoints lack individualCount info even though they have depth info. As
# each occurrence should have a minimum of one individual, we can assign NA values
# in individualCount with 1.

swebif_data$individualCount[which(is.na(swebif_data$individualCount))] <- 1

# This script will identify a taxon by name, count the number of occurrences per depth
# and list them from smallest to largest, extracting both previously specified percentiles 
# and medians from the distribution.

extract_gbif_depth_info <- function(df){
  
  individuals_per_depth <- df %>%
    group_by(taxon, depth) %>%
    summarise(individuals = sum(individualCount, na.rm = TRUE))
  
  depth_info <- individuals_per_depth %>%
    filter(!is.na(depth)) %>%
    group_by(taxon) %>%
    summarise(
      depth_repeated = list(rep(depth, individuals)),
      occurrences = sum(individuals)
    ) %>%
    rowwise() %>%
    mutate(
      lower_depth = round(quantile(unlist(depth_repeated), probs = percentiles[1], na.rm = TRUE)),
      median_depth = round(quantile(unlist(depth_repeated), probs = 0.5, na.rm = TRUE)),
      upper_depth = round(quantile(unlist(depth_repeated), probs = percentiles[2], na.rm = TRUE))
    ) %>%
    select(taxon, lower_depth, median_depth, upper_depth, occurrences) %>%
    ungroup()
  
  return(depth_info)
  
}

# Run it!
gbif_depth_ranges <- extract_gbif_depth_info(swebif_data)

# Check it!
gbif_depth_ranges

# Write it!
write.csv(gbif_depth_ranges, "Path/to/where/I/want/it.csv", row.names = FALSE, quote = FALSE)



# Now, for the temperature assignments, load the global GBIF data into R.

gbif_data_path <- "Path/to/chp_3_GBIF_Globaldata.csv"
gbif_data <- read.csv(gbif_data_path, sep = "\t")

# Filter out the "Koster Historical Biodiversity Assessment" dataset.

data_identifier <- "51d0bd32-e215-45ea-a04d-47a474336125"

gbif_data <- gbif_data[gbif_data$datasetKey != data_identifier,]


# As Bio-ORACLE has two separate temperature layers (decadal mean 2000s & 2010s) date is 
# important to get a correct value. Dates only including year or year/month become 
# NAs when converted to date, therefore use this function to include incomplete dates.
restructure_dates <- function(date) {
  
  # If the date is just a year (e.g., "2004")
  if (grepl("^\\d{4}$", date)) {
    # Assign the closest date to half the year (June 30th)
    result <- paste0(date, "-06-30")
    return(result)
  }
  # If the date has a year and month (e.g., "2011-05")
  else if (grepl("^\\d{4}-\\d{2}$", date)) {
    # Assign the 15th of that month
    result <- paste0(date, "-15")
    return(result)
  }
  # If the date spans two months (e.g., "2006-05/2006-06")
  else if (grepl("^\\d{4}-\\d{2}/\\d{4}-\\d{2}$", date)) {
    range <- strsplit(date, "/")[[1]]
    second_month <- range[2]
    return(paste0(second_month, "-01")) # First day of the second month
  }
  # If the date is already in full date format (e.g., "2011-05-12")
  else {
    result <- date
    return(result)
  }
}

# Apply this function while converting eventDate to Date format.
gbif_data$eventDate <- as.Date(sapply(gbif_data$eventDate, restructure_dates))



# OPTIONAL If you want to group by different taxonomic ranks, this will add a taxon column 
# that allows multiple species to be grouped under e.g. genus or family in the same column
# as species.

gbif_data <- gbif_data %>%
  mutate(taxon = case_when(
    species == "Protanthea simplex" ~ "Protanthea simplex",
    species == "Sabella pavonina" ~ "Sabella pavonina",
    species == "Geodia barretti" ~ "Geodia barretti",
    species == "Mycale lingua" ~ "Mycale lingua",
    species == "Polycarpa pomaria" ~ "Polycarpa pomaria",
    family == "Serpulidae" ~ "Serpulidae",
    genus == "Molgula" ~ "Molgula",
    genus == "Ascidia" ~ "Ascidia",
    species == "Ciona intestinalis" ~ "Ciona intestinalis",
    species == "Porania pulvillus" ~ "Porania pulvillus",
    species == "Acesta excavata" ~ "Acesta excavata",
    family == "Actiniidae" ~ "Actiniidae",
    genus == "Munida" ~ "Munida",
    species == "Phakellia ventilabrum" ~ "Phakellia ventilabrum",
    species == "Lithodes maja" ~ "Lithodes maja",
    species == "Caryophyllia smithii" ~ "Caryophyllia smithii",
    species == "Alcyonium digitatum" ~ "Alcyonium digitatum",
    TRUE ~ NA_character_  # Default if no condition is met
  )) %>%
  filter( #Remove uncertain positions
    coordinateUncertaintyInMeters <= 100 | is.na(coordinateUncertaintyInMeters)
  ) %>%
  filter( #Remove NA positions
    !is.na(decimalLatitude)
  )



# Separate data by decade to match Bio-Oracle layers.

gbif_data_2000 <- gbif_data %>%
  filter(eventDate >= "2000-01-01" & eventDate <= "2009-12-31")

gbif_data_2010 <- gbif_data %>%
  filter(eventDate >= "2010-01-01" & eventDate <= "2019-12-31")


# Set up parameters for Bio-ORACLE data download.


# Write out the identifier of your desired variable at max/min/mean depth of each 
# Bio-ORACLE grid cell. Here temperature is chosen.

dataset_id_shallow <- "thetao_baseline_2000_2019_depthmin"
dataset_id_mean <- "thetao_baseline_2000_2019_depthmean"
dataset_id_deep <- "thetao_baseline_2000_2019_depthmax"


# Identifier for the two different decade layers.

time = c('2000-01-01T00:00:00Z', "2010-01-01T00:00:00Z")

# Check the geographic bounds of GBIF data
print(paste("Latitude - min:", min(gbif_data$decimalLatitude), "max:", max(gbif_data$decimalLatitude)))
print(paste("Longitude - min:", min(gbif_data$decimalLongitude), "max:", max(gbif_data$decimalLongitude)))

# Set the geographic bounds. Values must be rounded to 0.025 precision and never beyond
# +- 179.975
latitude = c(-43.35, 81.35)
longitude = c(-179.975, 179.975)

# Combine these parameters as a list.
constraints = list(time, latitude, longitude)
names(constraints) = c("time", "latitude", "longitude")

# Select the type of temperature value you want from each decade (e.g. mean, max, etc.)
variables = c("thetao_mean")

# Download Bio-ORACLE data.
shallowtemps <- download_layers(dataset_id_shallow, variables, constraints)
meantemps <- download_layers(dataset_id_mean, variables, constraints)
deeptemps <- download_layers(dataset_id_deep, variables, constraints)


# Extract a dataframe of temperature values separated by decade for each depth type.

#2000-2010 decadal mean temp
shallowtemp_2000 <- as.data.frame(shallowtemps$thetao_mean_1, xy = TRUE, na.rm = FALSE)
meantemp_2000 <- as.data.frame(meantemps$thetao_mean_1, xy = TRUE, na.rm = FALSE)
deeptemp_2000 <- as.data.frame(deeptemps$thetao_mean_1, xy = TRUE, na.rm = FALSE)

#2010-2020 decadal mean temps
shallowtemp_2010 <- as.data.frame(shallowtemps$thetao_mean_2, xy = TRUE, na.rm = FALSE)
meantemp_2010 <- as.data.frame(meantemps$thetao_mean_2, xy = TRUE, na.rm = FALSE)
deeptemp_2010 <- as.data.frame(deeptemps$thetao_mean_2, xy = TRUE, na.rm = FALSE)



# As most GBIF occurrences do not list the depth at which the occurrence was made, we do
# not know if the temperature at mean, max, or min depth is the most accurate. Therefore,
# this next part will remove all occurrences from grid cells that have a temperature
# difference between max/min depth that exceeds a predetermined threshold.

temp_diff_threshold <- 2

# Apply this cleanup to 2000s dataset
seafloor_temps_2000 <- shallowtemp_2000 %>%
  rename(
    temp_mindepth = thetao_mean_1
  ) %>%
  mutate(
    temp_meandepth = meantemp_2000$thetao_mean_1,
    temp_maxdepth = deeptemp_2000$thetao_mean_1,
    temp_diff = temp_mindepth - temp_maxdepth
  ) %>%
  mutate(temp_meandepth = ifelse(abs(temp_diff) > temp_diff_threshold, NA, temp_meandepth),
         temp_diff = ifelse(abs(temp_diff) > temp_diff_threshold, NA, temp_diff)
  )

# Apply this cleanup to 2010s dataset
seafloor_temps_2010 <- shallowtemp_2010 %>%
  rename(
    temp_mindepth = thetao_mean_2
  ) %>%
  mutate(
    temp_meandepth = meantemp_2010$thetao_mean_2,
    temp_maxdepth = deeptemp_2010$thetao_mean_2,
    temp_diff = temp_mindepth - temp_maxdepth
  ) %>%
  mutate(temp_meandepth = ifelse(abs(temp_diff) > temp_diff_threshold, NA, temp_meandepth),
         temp_diff = ifelse(abs(temp_diff) > temp_diff_threshold, NA, temp_diff)
  )


# OPTIONAL. This script usually processes some big data, here you can free up some memory
# by removing these variables.
rm(constraints,shallowtemp_2000,shallowtemp_2010,shallowtemps,meantemp_2000,meantemp_2010,meantemps,
   deeptemp_2000,deeptemp_2010,deeptemps)


# Create a Spatraster object of your cleaned temperature grids allowing for a geographic
# connection between GBIF occurrences and corresponding temperature values

# This creates an empty grid set to the same dimensions as the temperature grids
ext <- ext(min(seafloor_temps_2000$x), max(seafloor_temps_2000$x), min(seafloor_temps_2000$y), max(seafloor_temps_2000$y))
raster_grid <- rast(ext, resolution = 0.025)

# Convert the seafloor temp df to a vector and integrate with the empty grid
temp_vector_2000 <- vect(seafloor_temps_2000, geom = c("x", "y"))
temp_raster_2000 <- rasterize(temp_vector_2000, raster_grid, field = "temp_meandepth")

# OPTIONAL. Free up more space.
rm(temp_vector_2000, seafloor_temps_2000)

temp_vector_2010 <- vect(seafloor_temps_2010, geom = c("x", "y"))
temp_raster_2010 <- rasterize(temp_vector_2010, raster_grid, field = "temp_meandepth")

# OPTIONAL. Free up more space.
rm(temp_vector_2010, seafloor_temps_2010)

# Extract temperature from the SpatRaster temperature grid and assign a corresponing 
# temperature to each occurrence.

# 2000s
temperatures_2000 <- terra::extract(temp_raster_2000, gbif_data_2000[, c("decimalLongitude", "decimalLatitude")], ID = FALSE)
result_2000 <- cbind(gbif_data_2000, seafloor_temperature = temperatures_2000$last)

# 2010s
temperatures_2010 <- terra::extract(temp_raster_2010, gbif_data_2010[, c("decimalLongitude", "decimalLatitude")], ID = FALSE)
result_2010 <- cbind(gbif_data_2010, seafloor_temperature = temperatures_2010$last)

#Combine both dfs and remove NA temperatures
combined_results <- rbind(result_2000[which(!is.na(result_2000$seafloor_temperature)),],result_2010[which(!is.na(result_2010$seafloor_temperature)),])


# CHECKPOINT. You have now added temperature values to all GBIF occurrences that lie on 
# Bio-ORACLE grid-cells that don't have internal temperature differences of more than
# the previously set threshold. Export these results and continue with post-processing 
# below.

write.table(combined_results, "Path/to/where/you/want/your/result.csv", row.names = FALSE, quote = FALSE, sep = "\t")




# To avoid sampling bias, we want to remove occurrences that are too close to each other.
# This following section will create a function that randomly selects an occurrence and
# removes all occurrences within a 1km range.

# NOTE, as this function operates randomly output may vary

# Function to calculate distance between two sets of lat/lon coordinates
calculate_distance <- function(lat1, lon1, lat2, lon2) {
  distVincentySphere(cbind(lon1, lat1), cbind(lon2, lat2))  # Using Vincenty distance formula
}

# Function that selects a random row, adds it to a list, removes all nearby rows, and moves
# on to the next row until all are a set distance apart. This function will take some time. 
#Have a cup of tea, take a break and come back later :)


# ONLY TO BE USED WITH group_by() to ensure that this is operated within each taxon, 
# rather than for all occurrences.

filter_distances <- function(group, taxon_name, threshold) {
  selected_rows <- list()  # List to store selected rows
  total_rows <- nrow(group)
  print(paste("Processing taxon:", taxon_name, "- Total rows:", total_rows))
  
  iteration <- 1
  while (nrow(group) > 0) {
    # Select a random row from the group
    random_row <- group[sample(1:nrow(group), 1), ]
    
    # Add the selected row to the list
    selected_rows <- append(selected_rows, list(random_row))
    print(paste("Iteration", iteration, ":", taxon_name, "- Random row selected at",
                random_row$decimalLatitude, random_row$decimalLongitude))
    
    # Calculate the distance from the selected row to all other rows
    distances <- mapply(calculate_distance, 
                        lat1 = random_row$decimalLatitude, lon1 = random_row$decimalLongitude,
                        lat2 = group$decimalLatitude, lon2 = group$decimalLongitude)
    
    # Remove rows within 1 km distance
    group <- group[distances > threshold, ]
    print(paste("Rows remaining after removing nearby rows:", nrow(group)))
    
    # Increment the iteration counter
    iteration <- iteration + 1
  }
  
  print(paste("Completed taxon:", taxon_name, "- Selected rows:", length(selected_rows)))
  
  # Combine all selected rows for this taxon group
  return(do.call(rbind, selected_rows))
}


# Set a distance threshold for the minimum permitted distance for occurrences to be from
# one another

distance_threshold <- 1000

# Apply the process to each taxon

# # NOTE, as this function operates randomly output may vary

filtered_result <- combined_results %>%
  group_by(taxon) %>%
  group_modify(~ filter_distances(.x, .y, threshold = distance_threshold)) %>%
  ungroup()


# We also want to remove highly isolated occurrences as they likely are incorrectly labeled.

# This function removes occurrences far from all others
remove_isolated_rows <- function(result, taxon_col, max_distance) {
  
  # Function to process each taxon
  process_taxon_isolation <- function(group) {
    total_rows <- nrow(group)
    if (total_rows <= 1) return(group)  # If only one row, nothing to remove
    
    print(paste("Processing taxon:", unique(group[[taxon_col]]), "-", total_rows, "rows"))
    
    # Calculate pairwise distances
    pairwise_distances <- outer(1:total_rows, 1:total_rows, 
                                Vectorize(function(i, j) 
                                  calculate_distance(group$decimalLatitude[i], group$decimalLongitude[i], 
                                                     group$decimalLatitude[j], group$decimalLongitude[j])))
    
    # Calculate the minimum non-zero distance for each row
    nearest_distances <- apply(pairwise_distances, 1, function(x) min(x[x > 0], na.rm = TRUE))
    
    # Logical vector indicating whether to keep each row
    keep_indices <- rowSums(pairwise_distances <= max_distance & pairwise_distances > 0, na.rm = TRUE) > 0
    
    # Print messages for rows that are removed
    removed_rows <- which(!keep_indices)
    if (length(removed_rows) > 0) {
      messages <- paste0(
        "Removed row at Latitude: ", group$decimalLatitude[removed_rows],
        ", Longitude: ", group$decimalLongitude[removed_rows],
        " - Distance to nearest occurrence: ", round(nearest_distances[removed_rows], 2), " m"
      )
      cat(messages, sep = "\n")
    }
    
    # Return filtered group
    return(group[keep_indices, ])
  }
  
  # Apply the function grouped by the taxon column
  result_cleaned <- result %>%
    group_by(!!sym(taxon_col)) %>%
    group_split() %>%
    lapply(process_taxon_isolation) %>%
    bind_rows()
  
  return(result_cleaned)
}


# Set a maximum permitted distance in meters between an occurrence and its nearest neighbor.

max_dist <- 3000000


# Apply the function.

result_cleaned <- remove_isolated_rows(filtered_result, taxon_col = "taxon", max_distance = max_dist)

# OPTIONAL. Check each taxon's number of occurrences

for (dude in unique(result_cleaned$taxon)) {
  print(paste0(dude," no. of occurrences: ",length(result_cleaned$taxon[which(result_cleaned$taxon==dude)])))
}


# Now you have a result filtered to remove potentially mislabeled GBIF occurrences and
# reduced sampling bias by setting a maximum/minimum permitted distance between occurrences.

# Export this result.
write.table(result_cleaned, "Path/to/where/you/want/the/filtered_result.csv", row.names = FALSE, quote = FALSE, sep = "\t")


# OPTIONAL - This is where you can load the extracted occurrences with temperatures used
# by Nilsson et al. (2025). 
# old_temp_occurrences <- read.csv('Path/to/chp_3_publication_temp_occurrences.csv')


# To have a cleaner result, apply the following function to get descriptive statistics 
# (including median preferred temperature) of each taxon's temperature ranges.

# Find median as well as upper/lower limit of central 68% of temperatures for each taxon.
# Here you can swap 'result cleaned' for 'table_s3'
central_range <- result_cleaned %>%
  group_by(taxon) %>%
  summarise(
    lower_16th = quantile(seafloor_temperature, probs = 0.16, na.rm = TRUE),
    median = median(seafloor_temperature),
    upper_84th = quantile(seafloor_temperature, probs = 0.84, na.rm = TRUE),
    range_size = upper_84th-lower_16th,
    standard_deviation = sd(seafloor_temperature),
    n_occurrences = length(seafloor_temperature),
    .groups = "drop"
  ) %>%
  mutate(across(
    c(lower_16th, median, upper_84th, range_size, standard_deviation),
    ~ sprintf("%.5f", .)
  )) %>%
  arrange(median)


# Export the output
write.csv(central_range, "Path/to/where/you/want/central_ranges.csv", row.names = FALSE, quote = FALSE)
