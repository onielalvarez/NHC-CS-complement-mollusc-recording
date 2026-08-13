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
library(mgcv)
library(DHARMa)
library(mgcViz)
library(gratia)
library(ggplot2)


# Upload the filtered spatial data frame obtained from the data cleaning 

filtered_data_spatial_clean <- read.csv("Data/filtered_data_spatial_clean.csv", 
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

# examine the distribution of the response
summary(grid_source_summary$record_abundance)
summary(grid_source_summary$mean_abs_latitude)
summary(grid_source_summary$species_count)

# 7- Fit your model; the glm.nb function will log transform your data
# Use a maxit of 100 if your algorithm is not converging

nb_model <- glm.nb(record_abundance ~ mean_abs_latitude * source + log10(species_count),
  data = grid_source_summary,
  control = glm.control(maxit = 100))

summary(nb_model) # This will give you the results of the model

# check model assumptions
plot(nb_model)

# check dispersion
sum(residuals(nb_model, type = "pearson")^2) /
  df.residual(nb_model)

# let's try a GAM

gam_nb_test <- gam(
  record_abundance ~ source +
    s(mean_abs_latitude, by = source, k = 10) +
    s(species_count, k = 10),
  family = nb(),
  data = grid_source_summary,
  method = "REML"
)

summary(gam_nb_test)

# check if this is the right k value to use
# if edf is close to k then we need to try larger k
summary(gam_nb_test)$s.table

plot(gam_nb_test, pages = 1, shade = TRUE)


sim_gam <- simulateResiduals(gam_nb_test)

plot(sim_gam)

testDispersion(sim_gam)
testZeroInflation(sim_gam)
testOutliers(
  sim_gam,
  type = "bootstrap"
)
# the diagnostics are not bad!

# I also tested different k values for the smooth terms, and I found that the best diagnoistics came from 
# both smooth terms having a k of 10. If we go with the GAM model, we will have to explain that 
# different k values were tested and we landed on 10 for both smooth terms based on model fit.


# plot the data
draw(
  gam_nb_test,
  select = c(
    "s(mean_abs_latitude):sourceNHC",
    "s(mean_abs_latitude):sourceCS"
  )
) &
  labs(
    x = "Absolute latitude (°)",
    y = "Partial effect on record abundance"
  ) &
  theme_light()



