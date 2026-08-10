# Load required packages

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(stringr)
library(forcats)
library(scales)

#--------------------------------------------------------------------------

# 1 - taxonomic complementarity between data sources

# Upload your data

filtered_data_taxa_clean <- read_csv("filtered_data_taxa_clean.csv")

# 1.1 - One row per species, class and source
species_source <- filtered_data_taxa_clean |>
  drop_na(class, source, species) |>
  distinct(class, source, species)

# 1.2- Presence of each species in each source
species_presence <- species_source |>
  group_by(class, species) |>
  summarise(
    in_NHC = any(source == "NHC"),
    in_CS  = any(source == "CS"),
    .groups = "drop") |>
  mutate(status = case_when(
      in_NHC & in_CS ~ "Both",
      in_NHC ~ "NHC only",
      TRUE ~ "CS only"))

# 1.3 - Count species by class and status
plot_data <- species_presence |>
  count(class, status)

plot_data$status <- factor(plot_data$status,
  levels = c("NHC only", "Both", "CS only"))

# 1.4 - Plot

ggplot(plot_data,
       aes(x = log(n),
           y = class,
           fill = status)) +
  geom_col(width = 0.9) +
  scale_fill_manual(values = c(
    "NHC only" = "#33638DFF",
    "Both"     = "#5AB4AC",
    "CS only"  = "#70DB65")) +
  labs(
    x = "Species count",
    y = NULL,
    fill = "Species representation:") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = c(0.80, 0.90),
    legend.background = element_rect(fill = "white", colour = "black"),
    legend.title = element_text(face = "bold")
  )

#--------------------------------------------------------------------------

# 2 - Proportional difference in species count between data sources 
# (taxonomic dominance by source)

# Upload your data

filtered_data_taxa_habitat <- read_csv("filtered_data_taxa_habitat.csv")

# 2.1 - Get the number of species per class, habitat and source

filtered_data_taxa_habitat <- filtered_data_taxa_habitat |>
  distinct(species, class, habitat, source)

# 2.2 - Set new labels to the class and habitat columns

filtered_data_taxa_habitat <- filtered_data_taxa_habitat |>
  mutate(class = str_to_title(class), habitat = str_to_title(habitat))

# 2.3- Convert to presence/absence per source

filtered_data_taxa_habitat_counts <- filtered_data_taxa_habitat |>
  mutate(present = 1) |>
  pivot_wider(names_from = source, values_from = present, values_fill = 0)

# 2.4- Summarise per class and habitat

summary_df <- filtered_data_taxa_habitat_counts |>
  group_by(class, habitat) |>
  summarise(CS = sum(CS > 0), NHC = sum(NHC > 0), total = CS + NHC, .groups = "drop") |>
  mutate(diff = (CS - NHC) / total) # This is the proportional difference

# 2.5- Keep only dominant direction

summary_df <- summary_df |>
  filter(diff != 0) |>
  mutate(direction = ifelse(diff > 0, "CS", "NHC"))

# 2.6- Order classes by strongest effect

summary_df <- summary_df |>
  mutate(diff = (CS - NHC) / total, abs_diff = abs(diff),
    class = forcats::fct_reorder(class, abs_diff, .fun = max))

# 2.7- Plot

# Prepare class factor
summary_df$class <- factor(summary_df$class)

# Assign small y offsets for habitats to avoid overlap; and a color for each habitat
habitat_levels <- unique(summary_df$habitat)
offsets <- seq(-0.2, 0.2, length.out = length(habitat_levels))
names(offsets) <- habitat_levels

summary_df <- summary_df |>
  mutate(y_pos = as.numeric(class) + offsets[habitat])

habitat_colors <- c(
  "Brackish" = "#bc177d",
  "Freshwater" = "#e4d609",
  "Marine" = "#1f78b4",
  "Terrestrial" = "#439b21")

# 2.8 - Plot
ggplot(summary_df, aes(x = diff, y = y_pos, color = habitat)) +
  geom_segment(aes(x = 0, xend = diff, y = y_pos, yend = y_pos), linewidth = 1.1) + # horizontal lines for each habitat
  geom_point(aes(y = y_pos), size = 2.5) + # Points for each proportional value by class
  geom_vline(xintercept = 0, color = "darkgray", linetype = "solid", linewidth = 1) + # vertical line at 0
  scale_y_continuous(
    breaks = 1:length(levels(summary_df$class)), # numeric y-axis with class labels
    labels = levels(summary_df$class)
  ) +
  scale_x_continuous(limits = c(-1, 1)) + # Proportions
  scale_color_manual(values = habitat_colors) +
  labs(x = "Proportional Difference in Species Count",
    y = "Mollusk Class",
    color = "Habitats") +
  annotate("text", x = -1, y = 1.8, label = "NHC Dominant", hjust = 0, vjust = 1.5, size = 4) + # annotate left/right labels
  annotate("text", x = 0.9, y = 1.8, label = "CS Dominant", hjust = 1, vjust = 1.5, size = 4) +
  theme_light() +
  theme(legend.position = "bottom",
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 11),
    axis.text = element_text(color = "black"))

