# For the complementary analysis you have to use two data sets, 
# one without the duplicates by habitats and one with them,
# run your models for each one (global and habitats) and graph based on both data sets

## Load the required packages

library(readr)
library(dplyr)
library(tidyverse)
library(performance)
library(see)
library(ggplot2)
library(forcats)
library(binom)

#-----------------------------------------------------------------------

# Load your data

filtered_data_taxa_clean <- read_csv("filtered_data_taxa_clean.csv")
filtered_data_taxa_habitat <- read_csv("filtered_data_taxa_habitat.csv")

# Pivot your data sets based on the data sources; this will allow a comparison
# for the recording of species between sources

# 1- Summary for the overall trend

summary_colms <- filtered_data_taxa_clean |>
  group_by(class, species, source) |>
  summarise(n_records = n(), .groups = "drop") |> # Counting number of records per species per source
  pivot_wider(names_from = source, values_from = n_records, values_fill = 0)

summary(lm(log10(CS + 1) ~ log10(NHC + 1), data = summary_colms))

# 1.1- Overall R2 label

r2_label_overall <- summary(lm(log10(CS + 1) ~ log10(NHC + 1), data = summary_colms))$r.squared

r2_label_overall <- tibble::tibble(label = paste0("R² = ", round(r2_label_overall, 2)))

# 2- Summary for each habitat trend

summary_colms_hab <- filtered_data_taxa_habitat |>
  group_by(class, species, source, habitat) |> 
  summarise(n_records = n(), .groups = "drop") |> # Counting number of records 
                                                  #per species per source per habitat
  pivot_wider(names_from = source, values_from = n_records, values_fill = 0)

# Marine

summary_colms_hab_mar <- filtered_data_taxa_habitat |>
  group_by(class, species, source, habitat) |> 
  filter(habitat == "marine") |>
  summarise(n_records = n(), .groups = "drop") |> # Counting number of records 
  #per species per source per habitat
  pivot_wider(names_from = source, values_from = n_records, values_fill = 0)

summary(lm(log10(CS + 1) ~ log10(NHC + 1), data = summary_colms_hab_mar))

# 2.1- R2 label for each habitat

r2_labels_hab <- summary_colms_hab |>
  group_by(habitat) |>
  summarise(r2 = summary(lm(log10(CS + 1) ~ log10(NHC + 1)))$r.squared,
    label = paste0("R² = ", round(r2, 2)),
    .groups = "drop")

#-----------------------------------------------------------------------

# Relationship graphs; one for the overall trend and a facet_wrap graph for
# the habitats; add a dashed abline for the 1 X 1 expected trend

# 3- OVERALL TREND

# 3.1- Graph plot for overall trend

overall_trend_graph <- ggplot(summary_colms,
       aes(x = log10(NHC + 1),
           y = log10(CS + 1))) +
  geom_point(
    colour = "darkgrey",
    alpha = 0.1,
    size = 1.5
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    colour = "red"
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    colour = "darkgreen",
    linewidth = 0.8
  ) +
  geom_text(
    data = r2_label_overall,
    aes(x = 0.2, y = 5, label = label), # This is the X and Y label positions
    inherit.aes = FALSE,
    hjust = 0,
    size = 3.5
  ) +
  labs(
    x = "Natural History Collection records",
    y = "Citizen Science observations"
  ) +
  theme_light()

overall_trend_graph

#-----------------------------------------------------------------------

# 4- TRENDS PER HABITAT

# 4.1- Colors for each habitat trend

hab_cols <- c(
  "terrestrial" = "#439b21",
  "marine"      = "#1f78b4",
  "freshwater"  = "#e4d609",
  "brackish"    = "#bc177d"
)

# 4.2- Graph plot for habitat trends

habitat_trends_graph <- ggplot(summary_colms_hab,
       aes(x = log10(NHC + 1),
           y = log10(CS + 1))) +
  geom_point(
    aes(colour = 'darkgrey'),
    alpha = 0.08,
    size = 0.5
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    colour = "red"
  ) +
  geom_smooth(
    aes(colour = habitat),
    method = "lm", # these are the trends for each habitat
    se = TRUE,
    linewidth = 0.8
  ) +
  geom_smooth(
    data = summary_colms,
    aes(x = log10(NHC + 1),
        y = log10(CS + 1)),
    inherit.aes = FALSE,
    method = "lm", # this is the overall trend
    se = FALSE,
    colour = "black",
    linewidth = 0.3,
    linetype = "dashed"
  ) +
  geom_text(
    data = r2_labels_hab,
    aes(
      x = 0.2,
      y = 5.0,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3
  ) +
  facet_wrap(~habitat, ncol = 2) +
  scale_colour_manual(values = hab_cols) +
  labs(
    x = "Natural History Collection records",
    y = "Citizen Science observations"
  ) +
  theme_light() +
  theme(legend.position = "none")

habitat_trends_graph

##---------------------------------------------------------------------

# STATISTICAL TESTS

## Linear Models

model_global <- lm(log10(CS + 1) ~ log10(NHC + 1), data = summary_colms)

summary(model_global)

## Test of assumptions for models

check_heteroscedasticity(model_global)

