# Load the next packages:

library(readr)
library(dplyr)
library(tidyr)
library(sf)
library(terra)
library(rnaturalearth)
library(ggplot2)
library(scales)

##-----------------------------------------------------------------------

# 1- Upload your data

clean_data <- read_csv("filtered_data_taxa_clean.csv")

# 2 - Remove records without country code

clean_data <- clean_data |>
  filter(!is.na(countryCode))

# 3 - Add continent using country code

world <- ne_countries(scale = "medium", returnclass = "sf")

country_lookup <- world |>
  st_drop_geometry() |>
  select(iso_a2, continent) |>
  distinct()

clean_data <- clean_data |>
  left_join(country_lookup, by = c("countryCode" = "iso_a2"))

# 4 - Assign continent to records with missing continent using lat/long

world_sv <- ne_countries(scale = "medium", returnclass = "sv")

missing_ids <- which(is.na(clean_data$continent))

if (length(missing_ids) > 0) {
  pts_missing <- vect(clean_data[missing_ids, ],
    geom = c("decimalLongitude", "decimalLatitude"),
    crs = "EPSG:4326")
  
  coord_continent <- terra::extract(world_sv["continent"], pts_missing)
  
  clean_data$continent_coord <- NA_character_
  
  clean_data$continent_coord[missing_ids] <-
    coord_continent$continent
  
  clean_data <- clean_data |>
    mutate(continent = coalesce(continent, continent_coord)) |>
    select(-continent_coord)
}

# 5 - Calculate annual and cumulative records

# For CS

df_continent_CS_clean <- clean_data |>
  filter( source == "CS", year >= 2004, !is.na(year), !is.na(continent)) |>
  count(continent, year, name = "annual_records") |> # Number of records per continent per year
  group_by(continent) |>
  complete(year = seq(min(year), max(year), by = 1),
    fill = list(annual_records = 0)) |> # Add years with zero records
  arrange(continent, year) |>
  mutate(cumulative_records = cumsum(annual_records)) |> # Calculate cumulative records
  ungroup()

# For NHC

df_continent_NHC_clean <- clean_data |>
  filter( source == "NHC", !is.na(year), !is.na(continent)) |>
  count(continent, year, name = "annual_records") |> # Number of records per continent per year
  group_by(continent) |>
  complete(year = seq(min(year), max(year), by = 1),
           fill = list(annual_records = 0)) |> # Add years with zero records
  arrange(continent, year) |>
  mutate(cumulative_records = cumsum(annual_records)) |> # Calculate cumulative records
  ungroup()

# 6 - Calculate CAGR and create labels

# For CS

growth_labels_CS <- df_continent_CS_clean |>
  filter(cumulative_records > 0) |>
  arrange(continent, year) |>
  group_by(continent) |>
  summarise(first_year = first(year), last_year = last(year),
    first_records = first(cumulative_records),
    last_records = last(cumulative_records),
    annual_growth = ifelse(last_year > first_year,
      ((last_records / first_records) ^ (1 / (last_year - first_year)) - 1) * 100, #CAGR
      NA_real_),
    .groups = "drop") |>
  distinct(continent, .keep_all = TRUE) |>
  mutate(label_continent = recode(continent, "Seven seas (open ocean)" = "Open Ocean"),
    label = paste0(label_continent, ": ", round(annual_growth, 1), " %"))

# For MHC

growth_labels_NHC <- df_continent_NHC_clean |>
  filter(cumulative_records > 0) |>
  arrange(continent, year) |>
  group_by(continent) |>
  summarise(first_year = first(year), last_year = last(year),
            first_records = first(cumulative_records),
            last_records = last(cumulative_records),
            annual_growth = ifelse(last_year > first_year,
                                   ((last_records / first_records) ^ (1 / (last_year - first_year)) - 1) * 100, #CAGR
                                   NA_real_),
            .groups = "drop") |>
  distinct(continent, .keep_all = TRUE) |>
  mutate(label_continent = recode(continent, "Seven seas (open ocean)" = "Open Ocean"),
         label = paste0(label_continent, ": ", round(annual_growth, 1), " %"))

# 7 - Create label position

# For CS

label_positions_CS <- df_continent_CS_clean |>
  group_by(continent) |>
  slice_max(order_by = year, n = 1, with_ties = FALSE) |>
  ungroup() |>
  distinct(continent, .keep_all = TRUE) |>
  left_join(growth_labels |>
              select(continent, label), by = "continent")

# For NHC

label_positions_NHC <- df_continent_NHC_clean |>
  group_by(continent) |>
  slice_max(order_by = year, n = 1, with_ties = FALSE) |>
  ungroup() |>
  distinct(continent, .keep_all = TRUE) |>
  left_join(growth_labels_NHC |>
              select(continent, label), by = "continent")

# 8 - Plot

# For CS

temp_plot_CS <- ggplot(df_continent_CS_clean, aes(x = year, 
                                                  y = cumulative_records, 
                                                  color = continent)) +
  geom_line(linewidth = 1.2) +
  geom_text(data = label_positions_CS, aes(x = year, 
                                        y = cumulative_records, 
                                        label = label),
    color = "black",
    hjust = -0.05,
    show.legend = FALSE) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(breaks = pretty_breaks(), expand = expansion(mult = c(0.02, 0.18))) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "Number of CS observations", color = "Continent") +
  theme_classic() +
  theme(legend.position = "none", plot.margin = margin(10, 130, 10, 10))

temp_plot_CS

# For NHC

temp_plot_NHC <- ggplot(df_continent_NHC_clean, aes(x = year, 
                                                  y = cumulative_records, 
                                                  color = continent)) +
  geom_line(linewidth = 1.2) +
  geom_text(data = label_positions_NHC, aes(x = year, 
                                        y = cumulative_records, 
                                        label = label),
            color = "black",
            hjust = -0.05,
            show.legend = FALSE) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(breaks = pretty_breaks(), expand = expansion(mult = c(0.02, 0.18))) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "Number of NHC records", color = "Continent") +
  theme_classic() +
  theme(legend.position = "none", plot.margin = margin(10, 130, 10, 10))

temp_plot_NHC

# You can create your figure with the resulting plots on Adobe Photoshop or Illustrator
