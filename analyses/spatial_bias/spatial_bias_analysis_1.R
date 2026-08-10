# Load required packages

library(dplyr)
library(tidyverse)
library(ggplot2)
library(sf)
library(terra)
library(rnaturalearth)
library(cowplot)
library(MASS)
library(ggeffects)


# Upload the filtered spatial data frame obtained from the data cleaning 

filtered_data_spatial_clean <- read.csv("filtered_data_spatial_clean.csv", 
                                 header = TRUE)

## ---------------------------------------------------------------------
# Setting coordinates as a sf object
# Creating the grid for the analysis

# 1- Set lon/lat as sf object; set the crs as WGS84
occ_sf <- st_as_sf(filtered_data_spatial_clean,
  coords = c("decimalLongitude", "decimalLatitude"),
  crs = 4326,
  remove = FALSE)

# 2- Get the world country boundaries; set it as sf
world <- ne_countries(scale = "medium", returnclass = "sf")

# 3- Re-project the country boundaries and occ coordinates to Equal Earth proj
world_eq <- st_transform(world, 8857)
occ_eq   <- st_transform(occ_sf, 8857)

# 4- Create global equal-area hexagonal grid; set the cell size to 150 Km
grid <- st_make_grid(
  world_eq, # Use the re-projected (crs = 8857) world country boundaries
  cellsize = 150000,
  square = FALSE) |>
  st_sf() |>
  mutate(hex_id = row_number()) # Assign an id for each cell 

# 5- Join the occurrence records to the grid
sf_data_joined <- st_join(occ_eq, grid, join = st_within)

# 6- Summarize the data per Hex ID and Source; 

grid_source_summary <- sf_data_joined |>
  st_drop_geometry() |>
  filter(!is.na(hex_id), !is.na(source), !is.na(species)) |>
  mutate(abs_latitude = abs(decimalLatitude),
    source = factor(source, levels = c("NHC", "CS"))) |> # set your data sources as factors
  group_by(hex_id, source) |>
  summarise( # the next values are summarized by the hex_id and the source
    record_abundance = n(), # number of records
    species_count = n_distinct(species), # species count
    mean_abs_latitude = mean(abs_latitude, na.rm = TRUE), # absolute latitude
    .groups = "drop")

# Make sure the number of records is close (after NA drops) to the dataset used for the analysis
sum(grid_source_summary$record_abundance)

# 7- Fit your model; the glm.nb function will log transform your data
# Use a maxit of 100 if your algorithm is not converging

nb_model <- glm.nb(record_abundance ~ mean_abs_latitude * source + species_count,
  data = grid_source_summary,
  control = glm.control(maxit = 100))

summary(nb_model) # This will give you the results of the model

# 8- Graph your predicted model if necessary

# Firstly, predict your results based on your model

pred <- ggpredict(nb_model, terms = c("mean_abs_latitude", "source"))

# Graph based on the predictions for each data source

ggplot(pred, aes(x = x, y = predicted, colour = group, fill = group)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.2,
    colour = NA) +
  labs(x = "Mean absolute latitude (°)",
    y = "Predicted record abundance",
    colour = "Source",
    fill = "Source") +
  theme_bw()
