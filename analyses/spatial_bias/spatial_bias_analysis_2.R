library(dplyr)
library(ggplot2)
library(rstatix)
library(tidyverse)
library(sf)
library(terra)
library(rnaturalearth)
library(cowplot)
library(MASS)
library(performance)
library(car)

# Upload the filtered habitat data frame obtained from the data integration 

filtered_data_taxa_habitat_spatial <- read.csv("filtered_data_taxa_habitat_spatial.csv", 
                                 header = TRUE)

## ---------------------------------------------------------------------
# Setting coordinates as a sf object
# Creating the grid for the analysis

# 1- Set lon/lat as sf object; set the crs as WGS84
occ_sf_hab <- st_as_sf(filtered_data_taxa_habitat_spatial,
  coords = c("decimalLongitude", "decimalLatitude"),
  crs = 4326,
  remove = FALSE)

# 2- Get the world country boundaries; set it as sf
world <- ne_countries(scale = "medium", returnclass = "sf")

# 3- Re-project the country boundaries and occ coordinates to Equal Earth proj
world_eq <- st_transform(world, 8857)
occ_eq_hab <- st_transform(occ_sf_hab, 8857)

# 4- Create global equal-area hexagonal grid; set the cell size to 150 Km
grid <- st_make_grid(world_eq, # Use the re-projected (crs = 8857) world country boundaries
  cellsize = 150000,
  square = FALSE) |>
  st_sf() |>
  mutate(hex_id = row_number()) # Assign an id for each cell 

# 5- Join the occurrence records to the grid
sf_data_joined_hab <- st_join(occ_eq_hab, grid, join = st_within)

# --------------------------------------------------------------------

# BOXPLOT ANALYSIS

habitat_summary_records <- sf_data_joined_hab |>
  st_drop_geometry() |>
  filter(!is.na(hex_id), !is.na(source), !is.na(habitat)) |>
  group_by(hex_id, source, habitat) |>
  summarise(record_abundance = n(), .groups = "drop")

habitat_summary_spp <- sf_data_joined_hab |>
  st_drop_geometry() |>
  filter(!is.na(hex_id), !is.na(source), !is.na(habitat), !is.na(species)) |>
  group_by(hex_id, source, habitat) |>
  summarise(species_richness = n_distinct(species), .groups = "drop")

## ------------------------------------------------------------
## I f you want the habitat boxplots for each source next to each other; instead of facet

source_cols <- c("NHC" = "#6BAED6", "CS" = "#83f28f")

plot_habt_rec <- ggplot(habitat_summary_records,
       aes(x = habitat,
           y = log1p(record_abundance),
           fill = source)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.alpha = 0.3) +
  scale_fill_manual(values = source_cols) +
  labs(x = "Habitat",
    y = "Record count per grid cell",
    fill = "Source") +
  theme_bw()
plot_habt_rec

plot_habt_spp <- ggplot(habitat_summary_spp,
       aes(x = habitat,
           y = log1p(species_richness),
           fill = source)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.alpha = 0.3) +
  scale_fill_manual(values = source_cols) +
  labs(x = "Habitat",
    y = "Species count per grid cell",
    fill = "Source") +
  theme_bw()
plot_habt_spp

## -------------------------------------------------------------------
## I f you want the habitat boxplots for each source not next to each other; and with facet

hab_cols <- c(
  "terrestrial" = "#439b21",
  "marine"      = "#1f78b4",
  "freshwater"  = "#e4d609",
  "brackish"    = "#bc177d")

plot_facet <- ggplot(habitat_summary,
       aes(x = habitat,
           y = log(record_abundance),
           fill = habitat)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
  facet_wrap(~source) +
  scale_fill_manual(values = hab_cols) +
  labs(x = "Habitat",
    y = "Record count per grid cell") +
  theme_bw() +
  theme(legend.position = "none")

plot_facet

## -------------------------------------------------------------------

# Kruskal-Wallis tests

# NHC
kruskal.test(species_richness ~ habitat, data = filter(habitat_summary_spp, source == "NHC"))

# CS
kruskal.test(record_abundance ~ habitat, data = filter(habitat_summary_records, source == "CS"))

## If the Kruskal–Wallis test is significant, identify which habitats differ:
  
habitat_summary |>
  filter(source == "NHC") |>
  dunn_test(record_abundance ~ habitat, p.adjust.method = "holm")

habitat_summary |>
  filter(source == "CS") |>
  dunn_test(record_abundance ~ habitat, p.adjust.method = "holm")
