

# Hello, and welcome to this R script in which you will investigate abundance trends
# of different taxa and functional traits. Then you will through statistical tests 
# investigate the impact of the different traits on the abundance change coefficients
# of the investigated taxa.


# Load dplyr.

library(dplyr)


# Read the "Koster Historical Biodiversity Assessment" df
my_data_path <- "Path/to/koster_historical_biodiversity_assessment_occurrence.txt"
my_data <- read.csv(my_data_path, sep = "\t")

# Use this function to assign each abundance score (grouped by depth and transect) with 
# values for year, depth, square root of depth, days of the year elapsed, and days of the 
# year elapsed squared. It then runs a GLM on the effect of these parameters on abundance 
# grouped by taxon.

create_taxonmodel <- function(taxon,df){
  
  taxonframe <- data.frame(abundance = scale(df$organismQuantity[which(df$scientificName==taxon)]),
                           year = as.numeric(substr(df$eventDate[which(df$scientificName==taxon)], 1,4)),
                           depth = scale(df$verbatimDepth[which(df$scientificName==taxon)]),
                           sqrt_depth = scale(sqrt(df$verbatimDepth[which(df$scientificName==taxon)]))
  )
  
  first_day_of_year <- as.Date(paste0(taxonframe$year, "-01-01"))
  
  days_elapsed <- as.numeric(difftime(as.Date(df$eventDate[which(df$scientificName==taxon)]), first_day_of_year, units = "days")) + 1
  
  taxonframe$day <- scale(days_elapsed)
  
  taxonframe$daysq <- scale(days_elapsed**2)
  
  taxonframe$year <- scale(taxonframe$year)
  
  print(length(taxonframe$year))
  
  return(summary(lm(abundance~year+depth+sqrt_depth+day+daysq, data = taxonframe)))
  
}

# Extract each taxon's name
unique_taxa <- unique(my_data$scientificName)

# Create an empty list to store results
taxonmodel_results <- list()

# Loop through each taxon and store results
for (taxon in unique_taxa) {
  cat("\nRunning model for:", taxon, "\n")  # Print progress
  taxonmodel_results[[taxon]] <- create_taxonmodel(taxon, my_data)
}

# Check results and make sure to manually copy them over to a .txt doc!
taxonmodel_results


# Traits part!

# Read the traits df
traits<- read.csv("Path/to/chp_4_traits_df.csv")

# Add each trait as a column to the original occurrence dataset, joining by taxon name.
my_data <- my_data %>%
  left_join(traits, by = "scientificName")


# Run the same type of model as for taxon-specific abundance change but grouped by trait
# type

create_traitmodel <- function(trait,type,df){
  
  traitframe <- data.frame(abundance = scale(df$organismQuantity[which(trait==type)]),
                           year = as.numeric(substr(df$eventDate[which(trait==type)], 1,4)),
                           depth = scale(df$verbatimDepth[which(trait==type)]),
                           sqrt_depth = scale(sqrt(df$verbatimDepth[which(trait==type)]))
  )
  
  first_day_of_year <- as.Date(paste0(traitframe$year, "-01-01"))
  
  days_elapsed <- as.numeric(difftime(as.Date(df$eventDate[which(trait==type)]), first_day_of_year, units = "days")) + 1
  
  traitframe$day <- scale(days_elapsed)
  
  traitframe$daysq <- scale(days_elapsed**2)
  
  traitframe$year <- scale(traitframe$year)
  
  print(paste0("Number of occurrences: ",length(traitframe$year)))
  
  return(summary(lm(abundance~year+depth+sqrt_depth+day+daysq, data = traitframe)))
  
}

# Run it, change column for a new trait, change value for abundance change of specific
# trait class.

# Extract each unique trait of a category
unique_traits <- unique(my_data$temp_category)

# Create an empty list to store results
traitmodel_results <- list()

# !Change trait column in this when running this for a new trait type. 

# Loop through each trait and store results
for (trait in unique_traits) {
  cat("\nRunning model for:", trait, "\n")  # Print progress
  traitmodel_results[[trait]] <- create_traitmodel(my_data$temp_category,trait, my_data)
}                                             # Change trait here!

# Check results and make sure to manually copy them over to a .txt doc!
traitmodel_results

# Repeat until you have all trait lm results! Remember that manual copying!



# Now, we can use the taxon-specific abundance-change coefficients to investigate drivers 
# behind their change.

# List taxon coefficients in the same order as unique(traits$scientificName)
taxon_coefficients = c(0.158, 0.030, -0.441, -0.018, 0.020, -0.227, 0.028, 
                 0.043, 0.116, -0.057, -0.092, 0.031, 0.067, 0.005,
                 -0.050, 0.0850, 0.040)

traits$coefficients <- taxon_coefficients

# Run linear regressions and wilcox-tests of how traits of preferred temperature, size, and 
# lifestyle impact the taxon-specific abundance change coefficient.
summary(lm(scale(coefficients)~scale(median_temp), data = traits))
wilcox.test(coefficients ~ size, data = traits)
wilcox.test(coefficients ~ lifestyle, data = traits)

# Done! Don't forget to copy these to .txt files either :)
