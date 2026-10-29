# Load the required packages for your data cleaning, filtering and integration

library(readr)
library(dplyr)
library(tidyverse)
library(worrms)
library(CoordinateCleaner)
library(data.table)

#-----------------------------------------------------------------

# Load your data

raw_data <- read_delim(
  "test_data/raw_data_gbif.csv",
  delim = "\t",
  quote = "",
  escape_double = FALSE,
  trim_ws = FALSE,
  show_col_types = FALSE
)

# General filtering to discard iNaturalist and observation.org records <2008 and <2004
# Discard taxa different from species

# 1- Create a subset with the iNat and obs.org records with years <2008 & <2004
# If you're using a simple format GBIF data file, use the publishingOrgKey column to
# be able to filter Observation.org data

CS_outlier_filter <- raw_data |>
  filter(
    (publishingOrgKey == "28eb1a3f-1c15-4a95-931a-4af90ecb574d" & year < 2008) | # this is iNaturalist pubOrgKey
    (publishingOrgKey == "c8d737e0-2ff8-42e8-b8fc-6b805d26fc5f" & year < 2004) # this is Observation.org pubOrgKey
    )

# 2- Create a vector of the GBIFID assignated to iNat and obs.org records with years <2008 & <2004

CS_outlier_filt_vect <- CS_outlier_filter$gbifID

# 3- Filter the filtered CS outlier vect and non-species taxa from the raw data
# Keep only the iNaturalist and observation.org records for HUMAN_OBSERVATION

filtered_raw_data <- raw_data |>
  filter(!gbifID %in% CS_outlier_filt_vect,
         basisOfRecord != "HUMAN_OBSERVATION" |
           publishingOrgKey %in% c("28eb1a3f-1c15-4a95-931a-4af90ecb574d", 
                                   "c8d737e0-2ff8-42e8-b8fc-6b805d26fc5f")) |>
  drop_na(species)

# Data Scope

range(filtered_raw_data$year, na.rm = TRUE)

# The next step is to know from what year 90% of the data are from

years_90_onward <- filtered_raw_data |>
  filter(!is.na(year)) |>
  count(year, name = "n_records") |>
  arrange(desc(year)) |>
  mutate(cum_prop = cumsum(n_records) / sum(n_records)) |>
  filter(cum_prop >= 0.90) |>
  slice(1) |>
  pull(year)

#-----------------------------------------------------------------

# Taxonomic harmonization and habitat classification

# 1- Create a species vector to run through the WoRMS API

species_vector_worms <- filtered_raw_data |>
  distinct(species) |>
  pull(species)

# 2- Create a loop that sends a request for 100 records at a time, then takes a 1 minute 
# break before sending the next request. Each export will be saved in the current working directory.
# This allows the WoRMS API to not clash.
# This is a safety in case one batch fails.

# 2.1- Split the vector into batches

species_batches <- split(species_vector_worms, ceiling(seq_along(species_vector_worms) / 100))

# 2.2- Create a "batches" folder in your directory before running the next lines of code

for(i in seq_along(species_batches)) {
  
  cat("Running batch", i, "of", length(species_batches), "\n")
  
  batch_result <- tryCatch({
    wm_records_names(
      name = species_batches[[i]],
      marine_only = FALSE
    )
  }, error = function(e) {
    cat("Batch", i, "failed:", conditionMessage(e), "\n")
    return(NULL)
  })
  
  # Save each batch
  saveRDS(batch_result, file = paste0("batches/wm_batch_", i, ".rds"))
  
  # Wait between batches
  if(i < length(species_batches)) Sys.sleep(60)
}

  # Save the batches into one variable

all_batches <- list.files(path = "batches", pattern = "wm_batch_.*\\.rds$", full.name=TRUE) |>
  lapply(readRDS) |>
  bind_rows(.id = "batch")

# 2.3- Save the results as separate files for further data filtering and integration

write_csv(all_batches, "tax_harm_data.csv")
saveRDS(all_batches, "tax_harm_data.RDS")

# 3- Load the dataset resulting from the taxonomic harmonization

tax_harm_results <- read_csv("tax_harm_data.csv")

#-----------------------------------------------------------------

# Taxonomic harmonization

# 1- Filter the fossil species and non-valid names

# Check if there's any fossil species names

tax_harm_results |>
  filter(isExtinct == 1)

# Drop the species that don't have valid names associated with them

tax_harm_results_filt <- tax_harm_results |>
  drop_na(valid_name) |>
  distinct(valid_name)

# 2- Taxonomic harmonization integration 

# Match the valid species names with the species from the original data set

filtered_data_taxa_clean <- filtered_raw_data |>
  inner_join(tax_harm_results_filt,
    by = c("species" = "valid_name")) |>
  mutate(basisOfRecord = recode(basisOfRecord,
                     "HUMAN_OBSERVATION" = "CS",
                     "PRESERVED_SPECIMEN" = "NHC")) |>
  rename("source" = "basisOfRecord") # re-classification of the data sources

# 3- Save the taxonomic harmonized data set; this one will be used for further analyses
# that don't contemplate double counting of records based on the habitat classification

write_csv(filtered_data_taxa_clean, "filtered_data_taxa_clean.csv")

#-----------------------------------------------------------------

filtered_data_taxa_clean <- read_csv("filtered_data_taxa_clean.csv")

# Clean dataset overview

filtered_data_taxa_clean |> # This is for knowing how many records per source
  count(source, name = "n_records")

filtered_data_taxa_clean |> # This is for knowing how many unique species
  summarise(n_distinct(species))

filtered_data_taxa_clean |> # This is to know how many records per each CS project
  filter(source == "CS") |>
  count(publishingOrgKey, name = "n_records")

filtered_data_taxa_clean |> # This is to know how many NHC projects
  filter(source == "NHC") |>
  drop_na(publishingOrgKey) |>
  count(publishingOrgKey, name = "n_records")

#-----------------------------------------------------------------

# Spatial filtering and function for CoordinateCleaner pipeline

# 1- Filter the records with uncertainty > 5000 and without coordinates
# Make sure that the coordinates are set between lat= -90, 90 & lon = -180, 180

filtered_spatial_data <- filtered_data_taxa_clean |>
  filter(is.na(coordinateUncertaintyInMeters) |
           coordinateUncertaintyInMeters <= 5000) |>
  drop_na(decimalLatitude, decimalLongitude) |>
  filter( between(decimalLatitude, -90, 90),
          between(decimalLongitude, -180, 180))

# 2- Create a folder for batches

dir.create("clean_coordinate_batches", showWarnings = FALSE)

# 3- Extract a vector of the species from the previously filtered data set

species <- unique(filtered_spatial_data$species)

# 4- Set the size of your batches; this can be done based on the amount of records 
# you have. We decided for 3000

batch_size <- 3000

# 5- Set the number of batches to run

species_batches <- split(
  species,
  ceiling(seq_along(species) / batch_size)
)

# 6- Process one batch at a time; this function is going to detect the problematic
# spatial records and get rid of them. If you have a manageable amount of records
# you can test for outliers with the clean_coordinates function

for (i in seq_along(species_batches)) {
  
  cat("Running batch", i, "of", length(species_batches), "\n")
  
  df <- filtered_spatial_data |>
    filter(species %in% species_batches[[i]])
  
  cleaned <- clean_coordinates(
    x = df,
    lon = "decimalLongitude",
    lat = "decimalLatitude",
    species = "species",
    tests = c("centroids", "equal", "gbif", "institutions", "zeros")
  ) |>
    filter(.summary) # This will keep the records that are not flagged
  
  saveRDS(
    cleaned,
    file = sprintf("clean_coordinate_batches/batch_%03d.rds", i)
  )
  
  rm(df, cleaned)
  gc()
}

# 7- Upload and join the RDS files for the obtained batches and turn them into a vector

files <- list.files(
  "clean_coordinate_batches",
  pattern = "\\.rds$",
  full.names = TRUE
)

# 8- Map the data frame from the vector and assign it to a new variable

filtered_data_spatial_clean <- map_dfr(files, readRDS)

# 9- Save the filtered spatial data set; this one will be used for further spatial analyses

write_csv(filtered_data_spatial_clean, "filtered_data_spatial_clean.csv")

#-----------------------------------------------------------------

# Habitat classification

# 1- Mutate the N/A values for a comprehensible classification

habitats_class <- tax_harm_results |>
  mutate(across(c(isMarine, isFreshwater, isTerrestrial, isBrackish),
                ~replace_na(., 0))) ## This replaces the N/A in these columns for 0

# 2- Classify into the four main categories
# Separate those species that have multiple habitats representation
# Get rid of the species that the habitat is unknown

habitats_class_filt <- habitats_class |>
  rowwise() |>
  mutate(habitat = paste(
    c("marine"[isMarine == 1],
      "brackish"[isBrackish == 1],
      "freshwater"[isFreshwater == 1],
      "terrestrial"[isTerrestrial == 1]),
    collapse = "; "
  )) |>
  ungroup() |>
  mutate(habitat = ifelse(habitat == "", "unknown", habitat)) |>
  separate_rows(habitat, sep = "; ") |>
  filter(habitat != "unknown")

# 3- Keep the columns of interest

habitats_class_final <- habitats_class_filt |>
  select(scientificname, valid_name, habitat) |>
  drop_na(valid_name)

# 4- Join the habitat classification to the taxonomic harmonized data frame
# Join the two data frames using the species column

filtered_data_taxa_habitat <- filtered_data_taxa_clean |>
  inner_join(
    habitats_class_final |>
      select(valid_name, habitat) |>
      distinct(),
    by = c("species" = "valid_name"),
    relationship = "many-to-many"
  )

# Save the resulting data set
# This data set is going to be used on source complementarity and taxonomic bias analyses

write_csv(filtered_data_taxa_habitat, "filtered_data_taxa_habitat.csv")

# 5- Join the habitat classification to the cleaned spatial data frame
# Join the two data frames using the species column

filtered_data_taxa_habitat_spatial <- filtered_data_spatial_clean |>
  inner_join(
    habitats_class_final |>
      select(valid_name, habitat) |>
      distinct(),
    by = c("species" = "valid_name"),
    relationship = "many-to-many"
  )

# Save the resulting data set
# This data set is going to be used on spatial bias per habitat analysis

write_csv(filtered_data_taxa_habitat_spatial, "filtered_data_taxa_habitat_spatial.csv")
