# Load required packages

library(dplyr)
library(tidyverse)
library(ggplot2)
library(sf)
library(terra)
library(rnaturalearth)
library(biscale)
library(cowplot)


# Upload the filtered spatial data frame obtained from the data cleaning 

filtered_data_spatial_clean <- read.csv("filtered_data_spatial_clean.csv", 
                                 header = TRUE)

## -----------------------------------------------------------------------------
# Setting coordinates as a sf object for further maps
# Creating the grid for the maps and analyses
# Reassignment of coordinate and map projection to Equal Earth

# 1- Set lon/lat as sf object; set the crs as WGS84
occ_sf <- st_as_sf(
  filtered_data_spatial_clean,
  coords = c("decimalLongitude", "decimalLatitude"),
  crs = 4326
)

# 2- Get the world country boundaries; set it as sf
world <- ne_countries(scale = "medium", returnclass = "sf")

# 3- Re-project the country boundaries and occ coordinates to Equal Earth proj
world_eq <- st_transform(world, 8857)
occ_eq   <- st_transform(occ_sf, 8857)

# 4- Create global equal-area hexagonal grid; set the cell size to 150 Km
grid <- st_make_grid(
  world_eq, # Use the re-projected (crs = 8857) world country boundaries
  cellsize = 150000,
  square = FALSE
) |>
  st_sf() |>
  mutate(hex_id = row_number()) # Assign an id for each cell

# 5- Join the occurrence records to the grid
sf_data_joined <- st_join(
  occ_eq,
  grid,
  join = st_within
)

## -----------------------------------------------------------------------------
# SOURCE COMPLEMENTARITY BIVARIATE MAP

# 1- Summarize NHC records and CS records per cell

source_comp_hex_summary <- sf_data_joined |>
  st_drop_geometry() |>
  filter(!is.na(hex_id)) |>
  group_by(hex_id) |>
  summarise(
    n_NHC = sum(source == "NHC", na.rm = TRUE),
    n_CS  = sum(source == "CS", na.rm = TRUE),
    .groups = "drop"
  ) |>
  right_join(grid, by = "hex_id") |>
  st_as_sf() |>
  mutate(
    n_NHC = replace_na(n_NHC, 0),
    n_CS  = replace_na(n_CS, 0),
    log_n_NHC = log10(n_NHC + 1),
    log_n_CS = log10(n_CS + 1)
  )

# 2- Create a bivariate classification

source_comp_hex_bi <- source_comp_hex_summary |>
  filter(n_NHC > 0, n_CS > 0) |>
  bi_class(
    x = log_n_NHC,
    y = log_n_CS,
    style = "quantile",
    dim = 3
  )

# 3- Assign a legend to the bivariate classification

source_comp_legend <- bi_legend(
  pal = "PurpleOr",
  dim = 3,
  xlab = "More records NHC",
  ylab = "More records CS",
  size = 8) +
  theme(axis.title = element_text(size = 29))

# 4- Create the bivariate map

graticules <- st_graticule(
  crs = st_crs(8857),
  datum = st_crs(4326),
  lon = seq(-180, 180, by = 30),
  lat = seq(-60, 60, by = 30))

# 5- Create right-side latitude labels

lat_labels <- data.frame(lon = rep(175, 5),
  lat = c(60, 30, 0, -30, -60),
  label = c("60°N", "30°N", "0°", "30°S", "60°S"))

lat_labels_sf <- st_as_sf(lat_labels, coords = c("lon", "lat"), crs = 4326) |>
  st_transform(8857)

# 6- Create bivariate map

source_comp_map <- ggplot() +
  # Graticules first — behind everything
  geom_sf(data = graticules, color = "grey70", linewidth = 0.25, fill = NA) +
  # Countries boundaries
  geom_sf(data = world_eq, fill = "white", color = "black", linewidth = 0.2) +
  # Bivariate hexagons
  geom_sf(data = source_comp_hex_bi, aes(fill = bi_class), alpha = 0.8, color = NA) +
  bi_scale_fill(pal = "PurpleOr", dim = 3) +
  coord_sf(crs = st_crs(8857), datum = NA, expand = FALSE, clip = "off") +
  geom_sf_text(data = lat_labels_sf, aes(label = label), size = 3, hjust = -0.2,
    color = "black") +
  bi_theme() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    plot.margin = margin(t = 5.5, r = 35, b = 5.5, l = 5.5))

source_comp_map

# 6.1- Save the bivariate legend of the map

ggsave(filename = "Figures/Legends/source_comp_bivariate_legend.png", plot = source_comp_legend,
       width = 5, height = 5, dpi = 300)

# 6.2- Do this if you want to join the bivariate map and legend in one template; 
# adjust the position of the legend to your liking

source_comp_joined_map <- ggdraw() +
  draw_plot(source_comp_map) +
  draw_plot(source_comp_legend,
            x = -0.025,
            y = 0.7,
            width = 0.22,
            height = 0.22)

source_comp_joined_map

# 6.3- Do this If you want to add the grid cells that only have CS or NHC records

NHC_only_hex <- source_comp_hex_summary |>
  filter(n_NHC > 0, n_CS == 0)

CS_only_hex <- source_comp_hex_summary |>
  filter(n_CS > 0, n_NHC == 0)

source_comp_map_2 <- ggplot() +
  geom_sf(data = NHC_only_hex,
    fill = "#6BAED6", # NHC-only cells
    color = NA,
    alpha = 0.8) +
  geom_sf(data = CS_only_hex,
    fill = "#83f28f", # CS-only cells
    color = NA,
    alpha = 0.8) +
  geom_sf(data = source_comp_hex_bi,
    aes(fill = bi_class), # Shared NHC + CS cells (bivariate)
    color = NA,
    alpha = 0.8) +
  bi_scale_fill(pal = "PurpleOr", dim = 3) +
  geom_sf(data = world_eq,
    fill = NA,
    color = "black",
    linewidth = 0.2) +
  coord_sf(crs = st_crs(8857)) +
  bi_theme() +
  theme(legend.position = "none")

source_comp_map_2

## -----------------------------------------------------------------------------
# SAMPLING EFFORT BIVARIATE MAP
# Filter the sf_data_joined for each of the sources being analyzed

# 1- Summarize records and species per cell

smp_eff_hex_summary_CS <- sf_data_joined |>
  st_drop_geometry() |>
  filter(source == "CS") |> # Source filtering here
  filter(!is.na(hex_id)) |>
  group_by(hex_id) |>
  summarise(
    n_obs = n(),
    n_species = n_distinct(species),
    .groups = "drop"
  ) |>
  right_join(grid, by = "hex_id") |>
  st_as_sf() |>
  mutate(
    n_obs = replace_na(n_obs, 0),
    n_species = replace_na(n_species, 0),
    log_n_obs = log10(n_obs + 1),
    log_n_species = log10(n_species + 1)
  )

smp_eff_hex_summary_NHC <- sf_data_joined |>
  st_drop_geometry() |>
  filter(source == "NHC") |> # Source filtering here
  filter(!is.na(hex_id)) |>
  group_by(hex_id) |>
  summarise(
    n_obs = n(),
    n_species = n_distinct(species),
    .groups = "drop"
  ) |>
  right_join(grid, by = "hex_id") |>
  st_as_sf() |>
  mutate(
    n_obs = replace_na(n_obs, 0),
    n_species = replace_na(n_species, 0),
    log_n_obs = log10(n_obs + 1),
    log_n_species = log10(n_species + 1)
  )

# 2- Create a bivariate classification based on one of the hex summary

smp_eff_hex_bi_CS <- smp_eff_hex_summary_CS |>
  filter(n_obs > 0, n_species > 0) |>
  bi_class(
    x = log_n_obs,
    y = log_n_species,
    style = "quantile",
    dim = 3 # This can be adjusted to be a 3 X 3 classification and so on
  )

smp_eff_hex_bi_NHC <- smp_eff_hex_summary_NHC |>
  filter(n_obs > 0, n_species > 0) |>
  bi_class(
    x = log_n_obs,
    y = log_n_species,
    style = "quantile",
    dim = 3 # This can be adjusted to be a 3 X 3 classification and so on
  )

# 3- Assign a legend to the bivariate classification

smp_eff_legend <- bi_legend(
  pal = "PinkGrn",
  dim = 3,
  xlab = "More records",
  ylab = "More species",
  size = 8) +
  theme(axis.title = element_text(size = 34))

# 4- Create the bivariate maps

smp_eff_map_CS <- ggplot() +
  # Graticules — behind continents
  geom_sf(data = graticules, color = "grey70", linewidth = 0.25) +
  # Countries
  geom_sf(data = world_eq, fill = "white", color = "black", linewidth = 0.2) +
  # Bivariate hexagons
  geom_sf(data = smp_eff_hex_bi_CS, aes(fill = bi_class), alpha = 0.8, color = NA) +
  # Bivariate colors
  bi_scale_fill(pal = "PinkGrn", dim = 3) +
  coord_sf(crs = st_crs(8857), datum = NA, expand = FALSE, clip = "off") +
  # Latitude labels
  geom_sf_text(data = lat_labels_sf, aes(label = label), size = 3, hjust = -0.2,
    color = "black") +
  theme_void() +
  theme(legend.position = "none", panel.border = element_blank(),
    plot.margin = margin(
      t = 5.5,
      r = 35,
      b = 5.5,
      l = 5.5))

smp_eff_map_CS

smp_eff_map_NHC <- ggplot() +
  # Graticules — behind continents
  geom_sf(data = graticules, color = "grey70", linewidth = 0.25) +
  # Countries
  geom_sf(data = world_eq, fill = "white", color = "black", linewidth = 0.2) +
  # Bivariate hexagons
  geom_sf(data = smp_eff_hex_bi_NHC, aes(fill = bi_class), alpha = 0.8, color = NA) +
  # Bivariate colors
  bi_scale_fill(pal = "PinkGrn", dim = 3) +
  coord_sf(crs = st_crs(8857), datum = NA, expand = FALSE, clip = "off") +
  # Latitude labels
  geom_sf_text(data = lat_labels_sf, aes(label = label), size = 3, hjust = -0.2,
    color = "black") +
  theme_void() +
  theme(legend.position = "none",
    panel.border = element_blank(),
    plot.margin = margin(
      t = 5.5,
      r = 35,
      b = 5.5,
      l = 5.5))

smp_eff_map_NHC

# 4.1- Save the bivariate legend of the map

ggsave(filename = "Figures/Legends/smp_eff_bivariate_legend.png", plot = smp_eff_legend,
  width = 5, height = 5, dpi = 300)

# 4.2- Do this if you want to Join the bivariate map and legend in one template; 
# adjust the position of the legend to your liking

final_smp_eff_map_CS <- ggdraw() +
  draw_plot(smp_eff_map_CS) +
  draw_plot(smp_eff_legend, x = -0.025, y = 0.78, width = 0.22, height = 0.22)

final_smp_eff_map_CS

final_smp_eff_map_NHC <- ggdraw() +
  draw_plot(smp_eff_map_NHC) +
  draw_plot(smp_eff_legend, x = -0.025, y = 0.78, width = 0.22, height = 0.22)

final_smp_eff_map_NHC

# if you want to do a multi-map figure, you can export the maps without the 
# bivariate legend attached to them; export the bivariate legend by itself
# and compose your figure in Adobe Photoshop or Illustrator

