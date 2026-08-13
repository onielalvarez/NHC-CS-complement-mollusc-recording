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

library(smatr)
library(ggpmisc)
library(lmodel2)

#-----------------------------------------------------------------------

# Load your data

filtered_data_taxa_clean <- read_csv("Data/filtered_data_taxa_clean.csv")
filtered_data_taxa_habitat <- read_csv("Data/filtered_data_taxa_habitat.csv")

# Pivot your data sets based on the data sources; this will allow a comparison
# for the recording of species between sources

# 1- Summary for the overall trend

summary_colms <- filtered_data_taxa_clean |>
  group_by(class, species, source) |>
  summarise(n_records = n(), .groups = "drop") |> # Counting number of records per species per source
  pivot_wider(names_from = source, values_from = n_records, values_fill = 0)

summary(summary_colms$CS)

model_overall <- lm(log10(CS + 1) ~ log10(NHC + 1), data = summary_colms)

# let's see if the model meets model assumptions
plot(model_overall)
# it does not, so we have to think about alternative ways to model the data.

summary(lm(log10(CS + 1) ~ log10(NHC + 1), data = summary_colms))
confint(lm(log10(CS + 1) ~ log10(NHC + 1), data = summary_colms))

# Spearman's rank sum test accommodates the strongly positively skewed observation counts and reduces sensitivity 
# to extreme values and differences in absolute reporting. It also does not assume linearity
# this tells us whether there is a correlation between species counts across platforms
cor.test(
  summary_colms$NHC,
  summary_colms$CS,
  method = "spearman"
)

# the paired wilcoxon test accounts for the fact that observations are paired (measurements of observations per species)
# this tells us if one source generally has higher counts per species than the other
wilcox.test(
  summary_colms$NHC,
  summary_colms$CS,
  paired = TRUE
)

# this model obtains a slope, intercept, and 95% confidence interval without requiring one variable to explain the other 
# (which is what is happening in the lm model), so regardless of which varaiable is response and which is predictor, the slope will 
# be the same. This is not the case with the linear model. 
# this tells us how the CS counts scale relative to NHC
sma_overall <- sma(
  log10(CS + 1) ~ log10(NHC + 1),
  data = summary_colms
)

summary(sma_overall)

summary(sma_overall)$Slope

# check model assumptions
plot(sma_overall, which="qq")
plot(sma_overall, which="residual")
plot(sma_overall, which="default")
# looks better than the linear model

sma_slope <- summary(sma_overall)$Slope
sma_intercept <- summary(sma_overall)$Int

# 1.1- Overall R2 label

r2_label_overall <- summary(sma_overall)$r2

r2_label_overall <- tibble::tibble(label = paste0("R² = ", round(r2_label_overall, 2)))

# 2- Summary for each habitat trend

summary_colms_hab <- filtered_data_taxa_habitat |>
  group_by(class, species, source, habitat) |> 
  summarise(n_records = n(), .groups = "drop") |> # Counting number of records 
                                                  #per species per source per habitat
  pivot_wider(names_from = source, values_from = n_records, values_fill = 0)

# Run SMA models for each habitat

r2_labels_hab <- summary_colms_hab |>
  group_by(habitat) |>
  summarise(r2 = summary(lm(log10(CS + 1) ~ log10(NHC + 1)))$r.squared,
    label = paste0("R² = ", round(r2, 2)),
    .groups = "drop")


sma_habitat <- summary_colms_hab |>
  group_by(habitat) |>
  group_modify(~ {
    model <- sma(
      log10(CS + 1) ~ log10(NHC + 1),
      data = .x
    )
    
    data.frame(
      slope = coef(model)[2],
      intercept = coef(model)[1]
    )
  })

sma_habitat <- summary_colms_hab |>
  group_by(habitat) |>
  group_modify(~ {
    
    model <- sma(
      log10(CS + 1) ~ log10(NHC + 1),
      data = .x
    )
    
    data.frame(
      slope = summary(model)$Slope,
      intercept = summary(model)$Int
    )
  })

sma_habitat

# check model assumptions
sma_models_habitat <- summary_colms_hab |>
  group_split(habitat) |>
  setNames(unique(summary_colms_hab$habitat)) |>
  lapply(function(dat) {
    
    sma(
      log10(CS + 1) ~ log10(NHC + 1),
      data = dat
    )
  })

# Terrestrial
plot(sma_models_habitat[["terrestrial"]], which = "qq")
plot(sma_models_habitat[["terrestrial"]], which = "residual")
plot(sma_models_habitat[["terrestrial"]], which = "default")

# Marine
plot(sma_models_habitat[["marine"]], which = "qq")
plot(sma_models_habitat[["marine"]], which = "residual")
plot(sma_models_habitat[["marine"]], which = "default")

# Freshwater
plot(sma_models_habitat[["freshwater"]], which = "qq")
plot(sma_models_habitat[["freshwater"]], which = "residual")
plot(sma_models_habitat[["freshwater"]], which = "default")

# Brackish
plot(sma_models_habitat[["brackish"]], which = "qq")
plot(sma_models_habitat[["brackish"]], which = "residual")
plot(sma_models_habitat[["brackish"]], which = "default")

# it looks okay

#-----------------------------------------------------------------------

# Relationship graphs; one for the overall trend and a facet_wrap graph for
# the habitats; add a dashed abline for the 1 X 1 expected trend

# 3- OVERALL TREND

# 3.1- Graph plot for overall trend

overall_trend_graph <- ggplot(
  summary_colms,
  aes(
    x = log10(NHC + 1),
    y = log10(CS + 1)
  )
) +
  geom_point(
    colour = "darkgrey",
    alpha = 0.1,
    size = 1.5
  ) +
  
  # 1:1 line
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    colour = "red"
  ) +
  
  # SMA line
  stat_ma_line(
    se = TRUE,
    linewidth = 0.8
  ) +
  
  geom_text(
    data = r2_label_overall,
    aes(
      x = 0.2,
      y = 5,
      label = label
    ),
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

habitat_trends_graph <- ggplot(
  summary_colms_hab,
  aes(
    x = log10(NHC + 1),
    y = log10(CS + 1)
  )
) +
  
  geom_point(
    colour = "darkgrey",
    alpha = 0.08,
    size = 0.5
  ) +
  
  # 1:1 line
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    colour = "red"
  ) +
  
  # Habitat-specific SMA
  stat_ma_line(
    aes(colour = habitat),
    se = TRUE,
    linewidth = 0.8
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


