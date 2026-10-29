# ------------------------------------------------------------------------------

# GENERAL DATA OVERVIEW

# ------------------------------------------------------------------------------
# 1 - Load required packages

library(dplyr)
library(readr)
library(tidyverse)
library(httr)
library(jsonlite)

# ------------------------------------------------------------------------------
# 2 - Load data

filtered_data_taxa_clean <- read.csv(
  "filtered_data_taxa_clean.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------------------
# 3 - Clean institution codes and save original names

filtered_data <- filtered_data_taxa_clean |>
  mutate(
    institutionCode_original = institutionCode,
    institutionCode = trimws(as.character(institutionCode)),
    institutionCode = na_if(institutionCode, ""),
    institutionCode_key = tolower(institutionCode),
    publishingOrgKey = trimws(as.character(publishingOrgKey)),
    publishingOrgKey = na_if(publishingOrgKey, ""))

# ------------------------------------------------------------------------------
# 4 - Identify institutions associated with publishing organizations

publishing_org_datasets <- filtered_data |>
  distinct(institutionCode, publishingOrgKey) |>
  arrange(institutionCode, publishingOrgKey)

# ------------------------------------------------------------------------------
# 5 - Identify institution name variations by publishing organization

institution_variations <- filtered_data |>
  distinct(publishingOrgKey, institutionCode) |>
  group_by(publishingOrgKey) |>
  summarise(n_names = n_distinct(institutionCode, na.rm = TRUE),
    institution_names = paste(sort(unique(na.omit(institutionCode))), collapse = "; "),
    .groups = "drop") |>
  filter(n_names > 1) |>
  arrange(desc(n_names))

# ------------------------------------------------------------------------------
# 6 - Create dictionary of institution names

institution_dictionary <- c(
  "csiro" = "CSIRO",
  "csiro, australia" = "CSIRO",
  "csiro oceans and atmosphere" = "CSIRO Oceans and Atmosphere",
  "csiro oceans and atmosphere, australia" = "CSIRO Oceans and Atmosphere",
  
  "arc" = "ARC",
  "the atlantic reference centre (arc)" = "ARC",
  "the atlantic reference centre" = "ARC",
  
  "rom" = "ROM",
  "royal ontario museum" = "ROM",
  
  "invemar" = "INVEMAR",
  "instituto de investigaciones marinas y costeras (invemar)" = "INVEMAR",
  
  "mizpan" = "MIZPAN",
  "muzeum i instytut zoologii polskiej akademii nauk" = "MIZPAN",
  
  "cnr-ismar" = "CNR-ISMAR",
  "cnr-ismar-ve" = "CNR-ISMAR",
  
  "rbins" = "RBINS",
  "rbins-scientific heritage" = "RBINS",
  
  "zmb" = "ZMB",
  "zmb/moll" = "ZMB",
  
  "zsm" = "ZSM",
  "zsm/mol" = "ZSM",
  
  "io pan" = "IO PAN",
  "instytut oceanologii polskiej akademii nauk" = "IO PAN",
  
  "shirshov institute of oceanology ran" = "Shirshov Institute of Oceanology",
  "p.p. shirshov institute of oceanology" = "Shirshov Institute of Oceanology"
)

# ------------------------------------------------------------------------------
# 7 - Standardize institution codes

# Use the dictionary for known alternative names.
# For capitalization differences, use the first available spelling.

institution_name_lookup <- filtered_data |>
  filter(!is.na(institutionCode_key)) |>
  group_by(institutionCode_key) |>
  summarise(preferred_name = first(institutionCode), .groups = "drop")

filtered_data <- filtered_data |>
  left_join(institution_name_lookup, by = "institutionCode_key") |>
  mutate(institutionCode_standardized = coalesce(
      unname(institution_dictionary[institutionCode_key]),
      preferred_name)) |>
  select(-institutionCode_key, -preferred_name)

# ------------------------------------------------------------------------------
# 8 - Identify missing institution codes

missing_institutions <- filtered_data |>
  filter(is.na(institutionCode_standardized)) |>
  distinct(publishingOrgKey) |>
  filter(!is.na(publishingOrgKey))

# ------------------------------------------------------------------------------
# 9 - Function to retrieve publishing organization names from GBIF

get_publisher_name <- function(key) {
  
  tryCatch({
    
    url <- paste0("https://api.gbif.org/v1/organization/", key)
    
    response <- httr::GET(url, httr::timeout(30))
    
    httr::stop_for_status(response)
    
    result <- httr::content(response, as = "parsed", type = "application/json")
    
    if (is.null(result$title)) {
      return(NA_character_)
    }
    
    as.character(result$title)
    
  }, error = function(e) {
    
    message("Could not retrieve ", key, ": ", conditionMessage(e))
    
    NA_character_
    
  })
}

# ------------------------------------------------------------------------------
# 10 - Retrieve missing publishing organization names directly from GBIF

publisher_lookup <- missing_institutions |>
  mutate(publishingOrganization = vapply(publishingOrgKey, get_publisher_name,
      character(1)))

# ------------------------------------------------------------------------------
# 11 - Add publishing organization names to dataset

filtered_data <- filtered_data |>
  left_join(publisher_lookup, by = "publishingOrgKey")

# ------------------------------------------------------------------------------
# 12 - Create final standardized institution labels

# Original institution codes take priority.
# Publisher names are used only when institution codes are missing.

filtered_data <- filtered_data |>
  mutate(institution_label_type = case_when(
      !is.na(institutionCode_standardized) ~ "Institution",
      !is.na(publishingOrganization) &
        publishingOrganization != "" ~ "Publisher proxy",
      TRUE ~ "Unidentified"),
    institutionCode_standardized = coalesce(
      institutionCode_standardized,
      publishingOrganization))

# ------------------------------------------------------------------------------
# 13 - Check standardized institution names

institution_changes <- filtered_data |>
  distinct(institutionCode_original,
    institutionCode_standardized,
    institution_label_type,
    publishingOrgKey,
    publishingOrganization) |>
  arrange(institutionCode_standardized, institutionCode_original)

# ------------------------------------------------------------------------------
# 14 - Identify recovered institution labels

recovered_institutions <- filtered_data |>
  filter(is.na(institutionCode_original) |
      trimws(as.character(institutionCode_original)) == "") |>
  distinct(
    publishingOrgKey,
    publishingOrganization,
    institutionCode_standardized,
    institution_label_type) |>
  arrange(institutionCode_standardized)

# ------------------------------------------------------------------------------
# 15 - Identify NHC institutions and CS projects

summary_data <- filtered_data |>
  filter(source %in% c("CS", "NHC"),
    !is.na(species),
    trimws(as.character(species)) != "") |>
  mutate(species = trimws(as.character(species)),
    project_institution = trimws(as.character(institutionCode_standardized)),
    project_institution = na_if(project_institution, "")) |>
  filter(!is.na(project_institution)) |>
  mutate(dataset_id = paste(source, project_institution, sep = " - "))

# ------------------------------------------------------------------------------
# 16 - Identify unique and shared species across all datasets

species_datasets <- summary_data |>
  distinct(dataset_id, species) |>
  group_by(species) |>
  mutate(n_datasets = n()) |>
  ungroup()

# ------------------------------------------------------------------------------
# 17 - Calculate number of occurrences and total species

dataset_counts <- summary_data |>
  group_by(source, project_institution, dataset_id) |>
  summarise(n_occurrences = n(), n_species = n_distinct(species), .groups = "drop")

# ------------------------------------------------------------------------------
# 18 - Calculate unique and shared species

species_counts <- species_datasets |>
  group_by(dataset_id) |>
  summarise(n_unique_species = sum(n_datasets == 1),
    n_shared_species = sum(n_datasets > 1),
    .groups = "drop")

# ------------------------------------------------------------------------------
# 19 - Identify most recorded species in each dataset

most_recorded <- summary_data |>
  count(dataset_id, species, name = "n_records") |>
  group_by(dataset_id) |>
  arrange(desc(n_records), species, .by_group = TRUE) |>
  slice_head(n = 1) |>
  ungroup() |>
  mutate(most_recorded_species = paste0(species, " (", n_records, ")")) |>
  select(dataset_id, most_recorded_species)

# ------------------------------------------------------------------------------
# 20 - Combine all summaries

dataset_summary <- dataset_counts |>
  left_join(species_counts, by = "dataset_id") |>
  left_join(most_recorded, by = "dataset_id") |>
  select(
    source,
    project_institution,
    n_occurrences,
    n_species,
    n_unique_species,
    n_shared_species,
    most_recorded_species) |>
  arrange(source, desc(n_occurrences))

# ------------------------------------------------------------------------------
# 21 - Verify species counts

stopifnot(all(dataset_summary$n_species ==
      dataset_summary$n_unique_species +
      dataset_summary$n_shared_species))

# ------------------------------------------------------------------------------
# 22 - Export summary

write.csv(dataset_summary, "CS_NHC_dataset_summary.csv", row.names = FALSE)

# ------------------------------------------------------------------------------

# SPATIAL DATA OVERVIEW

# ------------------------------------------------------------------------------

# 1 - Load your data

data <- read.csv("filtered_data_spatial_clean.csv", header = TRUE)

data |>
  summarise(n_species = n_distinct(species, na.rm = TRUE))

data |>
  filter(source == "NHC") |>
  summarise(n_species = n_distinct(species, na.rm = TRUE))

# ------------------------------------------------------------------------------
# 2 - Count occurrence records

record_counts <- data |>
  filter(source %in% c("NHC", "CS")) |>
  summarise(total_records = n(), NHC_records = sum(source == "NHC"),
    CS_records = sum(source == "CS"))

# ------------------------------------------------------------------------------
# 3 - Identify unique and shared species

NHC_species <- data |>
  filter(source == "NHC", !is.na(species), species != "") |>
  pull(species) |>
  unique()

CS_species <- data |>
  filter(source == "CS", !is.na(species), species != "") |>
  pull(species) |>
  unique()

NHC_only <- setdiff(NHC_species, CS_species)

CS_only <- setdiff(CS_species, NHC_species)

shared_species <- intersect(NHC_species, CS_species)

# ------------------------------------------------------------------------------
# 4 - Count species

species_counts <- data.frame(total_species = length(union(NHC_species, CS_species)),
  NHC_unique = length(NHC_only),
  CS_unique = length(CS_only),
  shared = length(shared_species))

# ------------------------------------------------------------------------------
# 5 - Combine results

summary_counts <- bind_cols(record_counts, species_counts)

print(summary_counts)

# ------------------------------------------------------------------------------
# 6 - Print results as a sentence

cat(
  format(summary_counts$total_records, big.mark = ","),
  " records (",
  format(summary_counts$NHC_records, big.mark = ","),
  " NHC, ",
  format(summary_counts$CS_records, big.mark = ","),
  " CS) and ",
  format(summary_counts$total_species, big.mark = ","),
  " species (",
  format(summary_counts$NHC_unique, big.mark = ","),
  " NHC unique, ",
  format(summary_counts$CS_unique, big.mark = ","),
  " CS unique, ",
  format(summary_counts$shared, big.mark = ","),
  " shared).",
  sep = ""
)
