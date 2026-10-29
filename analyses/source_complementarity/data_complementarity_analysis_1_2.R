# ------------------------------------------------------------------------------
# Load required packages

library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(scales)
library(eulerr)

# ------------------------------------------------------------------------------
# 1 - Load your data

filtered_data_taxa_clean <- read.csv("filtered_data_taxa_clean.csv", header = TRUE)

# ------------------------------------------------------------------------------

# VENN DIAGRAM

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# 2 - Create species lists

NHC_species <- filtered_data_taxa_clean %>%
  filter(source == "NHC") %>%
  pull(species) %>%
  unique() %>%
  na.omit()

CS_species <- filtered_data_taxa_clean %>%
  filter(source == "CS") %>%
  pull(species) %>%
  unique() %>%
  na.omit()

# ------------------------------------------------------------------------------
# 3 - Calculate unique and shared species

NHC_only <- setdiff(NHC_species, CS_species)

CS_only <- setdiff(CS_species, NHC_species)

shared_species <- intersect(NHC_species, CS_species)

# Count species
n_NHC <- length(NHC_only)
n_CS <- length(CS_only)
n_shared <- length(shared_species)

# Total unique species across both datasets
total_species <- n_NHC + n_CS + n_shared

# Calculate percentages
pct_NHC <- 100 * n_NHC / total_species
pct_CS <- 100 * n_CS / total_species
pct_shared <- 100 * n_shared / total_species

# ------------------------------------------------------------------------------
# 4 - Create area-proportional Venn diagram

venn_data <- c(
  "NHC" = n_NHC,
  "CS" = n_CS,
  "NHC&CS" = n_shared
)

venn_fit <- euler(venn_data)

# ------------------------------------------------------------------------------
# 5 - Create labels with counts and percentages

region_labels <- c(
  paste0(
    "NHC species\n",
    scales::comma(n_NHC),
    " (", round(pct_NHC, 1), "%)"
  ),
  paste0(
    "CS species\n",
    scales::comma(n_CS),
    " (", round(pct_CS, 1), "%)"
  ),
  paste0(
    "Shared species\n",
    scales::comma(n_shared),
    " (", round(pct_shared, 1), "%)"
  )
)

# ------------------------------------------------------------------------------
# 6 - Plot diagram

plot(
  venn_fit,
  
  fills = list(
    fill = c("#6BAED6", "#83f28f"),
    alpha = 0.8
  ),
  
  edges = list(
    col = "#092D3A",
    lwd = 2
  ),
  
  quantities = list(
    labels = region_labels,
    fontface = "bold",
    fontsize = 14
  ),
  
  labels = FALSE,
  
  main = NULL
)

# ------------------------------------------------------------------------------

# SPECIES RANK BASED ON NUMBER OF RECORDS

# ------------------------------------------------------------------------------
# 1 - Count records per species for NHC and CS

species_records <- filtered_data_taxa_clean %>%
  filter(
    source %in% c("NHC", "CS"),
    !is.na(species),
    trimws(as.character(species)) != ""
  ) %>%
  count(source, species, name = "n_records") %>%
  group_by(source) %>%
  arrange(desc(n_records), species, .by_group = TRUE) %>%
  mutate(
    species_rank = row_number(),
    cumulative_records = cumsum(n_records),
    cumulative_percent = 100 * cumulative_records / sum(n_records)
  ) %>%
  ungroup()

# ------------------------------------------------------------------------------
# 2 - Calculate the number of species accounting for 50% of records

half_records <- species_records %>%
  group_by(source) %>%
  summarise(
    total_species = n(),
    total_records = sum(n_records),
    n_species_50 = min(species_rank[cumulative_percent >= 50]),
    .groups = "drop"
  ) %>%
  mutate(
    legend_label = paste0(
      "50% of records: ",
      scales::comma(n_species_50),
      " species"
    ),
    percent_species_50 = 100 * n_species_50 / total_species
  )

# ------------------------------------------------------------------------------
# 3 - Identify the 10 species with the most records

top_species <- species_records %>%
  group_by(source) %>%
  slice_min(
    order_by = species_rank,
    n = 10,
    with_ties = FALSE
  ) %>%
  ungroup()

# ------------------------------------------------------------------------------
# 4 - Identify 10 random species with 1 to 100 records

set.seed(123)

random_species <- species_records %>%
  filter(
    n_records >= 10,
    n_records <= 100
  ) %>%
  group_by(source) %>%
  group_modify(
    ~ slice_sample(.x, n = min(10, nrow(.x)))
  ) %>%
  ungroup()

# ------------------------------------------------------------------------------
# 5 - Create a function to plot each source

plot_rank_records <- function(source_name) {
  
  # Filter data for selected source
  plot_data <- species_records %>%
    filter(source == source_name)
  
  # Get the 50% threshold
  threshold <- half_records %>%
    filter(source == source_name)
  
  # Select top 10 species
  labels_data <- top_species %>%
    filter(source == source_name)
  
  # Select 10 random species with 1 to 100 records
  random_labels_data <- random_species %>%
    filter(source == source_name)
  
  # Create graph
  ggplot(
    plot_data,
    aes(x = species_rank, y = n_records)
  ) +
    
    # Gray shading representing the first 50% of records
    annotate(
      "rect",
      xmin = 1,
      xmax = threshold$n_species_50,
      ymin = -Inf,
      ymax = Inf,
      fill = "gray75",
      alpha = 0.4
    ) +
    
    # Vertical line at 50% threshold
    geom_vline(
      xintercept = threshold$n_species_50,
      colour = "gray40",
      linewidth = 0.5,
      linetype = "dashed"
    ) +
    
    # Rank-record curve
    geom_line(
      colour = "black",
      linewidth = 0.7
    ) +
    
    # Points for top 10 species
    geom_point(
      data = labels_data,
      colour = "black",
      size = 1.8
    ) +
    
    # Labels for top 10 species
    geom_text_repel(
      data = labels_data,
      aes(label = species),
      size = 3.2,
      fontface = "italic",
      colour = "black",
      direction = "both",
      box.padding = 0.5,
      point.padding = 0.3,
      max.overlaps = Inf,
      min.segment.length = 0,
      seed = 123
    ) +
    
    # Points for 10 random low-record species
    geom_point(
      data = random_labels_data,
      colour = "red",
      size = 1.8
    ) +
    
    # Labels for 10 random species with 1 to 100 records
    geom_text_repel(
      data = random_labels_data,
      aes(label = species),
      size = 3.2,
      fontface = "italic",
      colour = "red",
      direction = "both",
      box.padding = 0.5,
      point.padding = 0.3,
      max.overlaps = Inf,
      min.segment.length = 0,
      seed = 123
    ) +
    
    # Label showing the 50% threshold
    annotate(
      "label",
      x = threshold$n_species_50,
      y = max(plot_data$n_records) * 0.65,
      label = threshold$legend_label,
      size = 3.5,
      fill = "white",
      linewidth = 0.2
    ) +
    
    # Log10 species rank
    scale_x_log10(
      breaks = c(1, 10, 100, 1000, 10000, 100000),
      labels = scales::label_comma(),
      expand = expansion(mult = c(0.02, 0.05))
    ) +
    
    # Log10 number of records
    scale_y_log10(
      breaks = scales::breaks_log(n = 6),
      labels = scales::label_comma(),
      expand = expansion(mult = c(0.02, 0.15))
    ) +
    
    labs(
      title = source_name,
      x = "Species rank (log10 scale)",
      y = "Number of occurrence records (log10 scale)"
    ) +
    
    coord_cartesian(clip = "off") +
    
    theme_classic() +
    
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "plain",
        size = 14
      ),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      plot.margin = margin(15, 30, 15, 15)
    )
}

# ------------------------------------------------------------------------------
# 6 - Create separate figures for NHC and CS

NHC_rank_plot <- plot_rank_records("NHC")

CS_rank_plot <- plot_rank_records("CS")

# Display figures
NHC_rank_plot

CS_rank_plot

# ------------------------------------------------------------------------------
# 7 - OPTIONAL: Combine figures

combined_plot <- NHC_rank_plot + CS_rank_plot +
  plot_layout(ncol = 2)

combined_plot

# ------------------------------------------------------------------------------
# 8 - Save figures

ggsave(
  "NHC_rank_records.png",
  NHC_rank_plot,
  width = 9,
  height = 6,
  dpi = 300
)

ggsave(
  "CS_rank_records.png",
  CS_rank_plot,
  width = 9,
  height = 6,
  dpi = 300
)

ggsave(
  "NHC_CS_rank_records.png",
  combined_plot,
  width = 16,
  height = 7,
  dpi = 300
)