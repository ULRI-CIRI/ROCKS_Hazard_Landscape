# ===========================================================
# fig03ab_coverage_density.R
# ===========================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(here)

#========================================================
# Output directory
#========================================================

out_dir <- here::here(
  "outputs",
  "figures"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

#========================================================
# Load data
#========================================================

rocks <- readRDS(
  here::here("data", "derived", "rocks.rds")
)

rocks_ids <- rocks %>%
  distinct(cichem_id)

hed_p_by_chem_ebmdhed <- readRDS(
  here::here("data", "derived", "hed_p_by_chem_ebmdhed.rds")
)

b_mc7_aed <- readRDS(
  here::here("data", "derived", "b_mc7_aed.rds")
)

bp_aed <- readr::read_csv(
  here::here(
    "data",
    "raw",
    "rocks_profile_pod_table_aeds_08sept26.csv"
  )
)
#========================================================
# Create final chemical-level minimum AED
#
# For each chemical, use the minimum positive AED from:
#   1. b_mc7_aed$p5_aed_ac50_sumcl50
#   2. bp_aed$min_pod_aed_sumcl50
#========================================================

mc7_aed_by_chem <- b_mc7_aed %>%
  transmute(
    cichem_id,
    aed_p05 = suppressWarnings(
      as.numeric(p5_aed_ac50_sumcl50)
    )
  ) %>%
  filter(
    !is.na(cichem_id),
    cichem_id != "",
    !is.na(aed_p05),
    is.finite(aed_p05),
    aed_p05 > 0
  ) %>%
  group_by(cichem_id) %>%
  summarise(
    aed_p05 = min(aed_p05),
    .groups = "drop"
  )


bp_aed_by_chem <- bp_aed %>%
  transmute(
    cichem_id,
    aed_min_pod = suppressWarnings(
      as.numeric(min_pod_aed_sumcl50)
    )
  ) %>%
  filter(
    !is.na(cichem_id),
    cichem_id != "",
    !is.na(aed_min_pod),
    is.finite(aed_min_pod),
    aed_min_pod > 0
  ) %>%
  group_by(cichem_id) %>%
  summarise(
    aed_min_pod = min(aed_min_pod),
    .groups = "drop"
  )


aed_final_by_chem <- full_join(
  mc7_aed_by_chem,
  bp_aed_by_chem,
  by = "cichem_id"
) %>%
  mutate(
    minimum_aed = case_when(
      !is.na(aed_p05) & !is.na(aed_min_pod) ~
        pmin(aed_p05, aed_min_pod),

      !is.na(aed_p05) ~
        aed_p05,

      !is.na(aed_min_pod) ~
        aed_min_pod,

      TRUE ~
        NA_real_
    ),

    minimum_aed_source = case_when(
      !is.na(aed_p05) &
        !is.na(aed_min_pod) &
        aed_p05 < aed_min_pod ~
        "Minimum AED",

      !is.na(aed_p05) &
        !is.na(aed_min_pod) &
        aed_min_pod < aed_p05 ~
        "Minimum broad-profiling POD AED",

      !is.na(aed_p05) &
        !is.na(aed_min_pod) &
        aed_p05 == aed_min_pod ~
        "Equal values",

      !is.na(aed_p05) ~
        "Minimum AED only",

      !is.na(aed_min_pod) ~
        "Minimum broad-profiling POD AED only",

      TRUE ~
        NA_character_
    )
  ) %>%
  filter(
    !is.na(minimum_aed),
    is.finite(minimum_aed),
    minimum_aed > 0
  ) %>%
  semi_join(
    rocks_ids,
    by = "cichem_id"
  )


#========================================================
# Prepare eBMDHED values
#========================================================

hed_plot <- hed_p_by_chem_ebmdhed %>%
  filter(
    !is.na(hed_p05),
    hed_p05 > 0,
    n_hed_values > 0
  ) %>%
  semi_join(
    rocks_ids,
    by = "cichem_id"
  ) %>%
  select(
    cichem_id,
    ebmdhed_p05 = hed_p05
  )

#========================================================
# Prepare final AED values
#========================================================

aed_plot <- aed_final_by_chem %>%
  select(
    cichem_id,
    minimum_aed,
    aed_p05,
    aed_min_pod,
    minimum_aed_source
  )

#========================================================
# Restrict to chemicals with both AED and eBMDHED values
#========================================================

potency_overlap <- inner_join(
  hed_plot,
  aed_plot,
  by = "cichem_id"
)

cat(
  "ROCKS chemicals with both AED and eBMDHED values:",
  scales::comma(n_distinct(potency_overlap$cichem_id)),
  "\n"
)

#========================================================
# Convert to long format
#========================================================

plot_df <- potency_overlap %>%
  pivot_longer(
    cols = c(
      ebmdhed_p05,
      minimum_aed
    ),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(
      metric,
      "ebmdhed_p05" = "p05 eBMDHED",
      "minimum_aed" = "Minimum AED"
    ),
    metric = factor(
      metric,
      levels = c(
        "Minimum AED",
        "p05 eBMDHED"
      )
    ),
    log10_value = log10(value)
  )

#========================================================
# Summary statistics for medians and labels
#========================================================

potency_summary <- plot_df %>%
  group_by(metric) %>%
  summarise(
    n_chem = n_distinct(cichem_id),
    median_value = median(value, na.rm = TRUE),
    median_log10_value = median(log10_value, na.rm = TRUE),
    p25_log10_value = quantile(log10_value, probs = 0.25, na.rm = TRUE),
    p75_log10_value = quantile(log10_value, probs = 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    metric_label = paste0(
      metric,
      " (n=",
      scales::comma(n_chem),
      ")"
    )
  )

print(potency_summary)

plot_df <- plot_df %>%
  left_join(
    potency_summary %>%
      select(
        metric,
        metric_label
      ),
    by = "metric"
  ) %>%
  mutate(
    metric_label = factor(
      metric_label,
      levels = potency_summary$metric_label
    )
  )

#========================================================
# Histogram + density overlay
#========================================================

col_hed <- "#BCE4F7"
col_aed <- "#0A32FF"

legend_labels <- c(
  "Minimum AED" =
    paste0(
      "Minimum AED (n=",
      comma(
        potency_summary$n_chem[
          potency_summary$metric == "Minimum AED"
        ]
      ),
      ")"
    ),

  "p05 eBMDHED" =
    paste0(
      "p05 eBMDHED (n=",
      comma(
        potency_summary$n_chem[
          potency_summary$metric == "p05 eBMDHED"
        ]
      ),
      ")"
    )
)

p_potency_overlap <- ggplot(
  plot_df,
  aes(
    x = log10_value,
    color = metric
  )
) +

geom_histogram(
  aes(
    fill = metric,
    y = after_stat(density)
  ),
  position = "identity",
  bins = 40
) +

geom_density(
  linewidth = 1.5
) +

  # median lines
  geom_vline(
    data = potency_summary,
    aes(
      xintercept = median_log10_value,
      color = metric
    ),
    linetype = "dashed",
    linewidth = 1
  ) +

scale_fill_manual(
  values = c(
    "Minimum AED" = col_aed,
    "p05 eBMDHED" = col_hed
  ),
  labels = legend_labels
) +

scale_color_manual(
  values = c(
    "Minimum AED" = "#1d4169",
    "p05 eBMDHED" = "#3bb3d8"
  ),
  labels = legend_labels
)+

  scale_x_continuous(
    breaks = seq(-8, 6, by = 1),
    minor_breaks = seq(-8, 6, by = 0.5)
  ) +

  labs(
    x = expression(
      log[10] *
        "(chemical-level potency, mg/kg-day)"
    ),
    y = "Density",
    fill = NULL,
    color = NULL
  ) +

  theme_minimal(base_size = 16) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )

p_potency_overlap

#========================================================
# Save figure 
#========================================================

ggsave(
  file.path(
    out_dir,
    "fig03b_potency_distribution_overlap.png"
  ),
  plot = p_potency_overlap,
  width = 6,
  height = 4,
  dpi = 300
)

#========================================================
# Panel A: Potency coverage of ROCKS inventory
#========================================================

# Colors matched to potency figure
potency_cols <- c(
  "Minimum AED\n(oral)" = "#0A32FF",
  "p05 eBMDHED\n(oral)" = "#BCE4F7",
  "Both" = "#122C49"
)

panel_a_chem <- rocks_ids %>%
  left_join(
    hed_plot %>%
      distinct(cichem_id) %>%
      mutate(has_ebmdhed = TRUE),
    by = "cichem_id"
  ) %>%
  left_join(
    aed_final_by_chem %>%
      distinct(cichem_id) %>%
      mutate(has_aed = TRUE),
    by = "cichem_id"
  ) %>%
  mutate(
    has_ebmdhed = coalesce(has_ebmdhed, FALSE),
    has_aed = coalesce(has_aed, FALSE)
  )

n_rocks <- n_distinct(rocks_ids$cichem_id)

panel_a_dat <- bind_rows(
  panel_a_chem %>%
    filter(has_aed) %>%
    summarise(
      group = "Minimum AED\n(oral)",
      n = n_distinct(cichem_id)
    ),

  panel_a_chem %>%
    filter(has_ebmdhed) %>%
    summarise(
      group = "p05 eBMDHED\n(oral)",
      n = n_distinct(cichem_id)
    ),

  panel_a_chem %>%
    filter(has_aed, has_ebmdhed) %>%
    summarise(
      group = "Both",
      n = n_distinct(cichem_id)
    )
) %>%
  mutate(
    pct_inventory = 100 * n / n_rocks
  ) %>%
  mutate(
    group = factor(
      group,
      levels = rev(c(
        "Minimum AED\n(oral)",
        "p05 eBMDHED\n(oral)",
        "Both"
      ))
    ),
    text_color = case_when(
      group == "p05 eBMDHED\n(oral)" ~ "black",
      TRUE ~ "white"
    )
  )


# Plot
panel_a <- ggplot(
  panel_a_dat,
  aes(
    x = pct_inventory,
    y = group,
    fill = group
  )
) +
  geom_col(
    width = 0.7
  ) +

  geom_text(
    aes(
      x = pct_inventory / 2,
      label = paste0(
        round(pct_inventory, 1),
        "%\n(n=",
        scales::comma(n),
        ")"
      ),
      color = text_color
    ),
    size = 5,
    lineheight = 0.9
  ) +

  scale_fill_manual(
    values = potency_cols
  ) +

  scale_color_identity() +

  guides(
    fill = "none"
  ) +

  scale_x_continuous(
    limits = c(
      0,
      max(panel_a_dat$pct_inventory) * 1.10
    ),
    expand = expansion(
      mult = c(0, 0.02)
    )
  ) +

  labs(
    x = "% of ROCKS inventory",
    y = NULL
  ) +

  theme_minimal(
    base_size = 16
  ) +

  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(
      size = 14
    ),
    axis.text.x = element_text(
      size = 13
    ),
    axis.title.x = element_text(
      size = 16
    )
  ) +
  theme(
    plot.margin = margin(
      t = 0,
      r = 5,
      b = 0,
      l = 5
    )
  )

panel_a

ggsave(
  file.path(
    out_dir,
    "fig03a_data_coverage.png"
  ),
  plot = panel_a,
  width = 6,
  height = 4,
  dpi = 300
)

#############################
# PANEL C #
#############################

chem_class <- readRDS(
  here("data", "derived", "chem_class.rds")
)


#========================================================
# Create final chemical-level AED
#
# Final AED is the minimum positive value across:
#   1. p5_aed_ac50_sumcl50 from b_mc7_aed
#   2. min_pod_aed_sumcl50 from bp_aed
#========================================================

mc7_aed_by_chem <- b_mc7_aed %>%
  transmute(
    cichem_id,
    aed_mc7 = as.numeric(p5_aed_ac50_sumcl50)
  ) %>%
  filter(
    !is.na(cichem_id),
    cichem_id != "",
    !is.na(aed_mc7),
    is.finite(aed_mc7),
    aed_mc7 > 0
  ) %>%
  group_by(cichem_id) %>%
  summarise(
    aed_mc7 = min(aed_mc7),
    .groups = "drop"
  )


bp_aed_by_chem <- bp_aed %>%
  transmute(
    cichem_id,
    aed_bp = as.numeric(min_pod_aed_sumcl50)
  ) %>%
  filter(
    !is.na(cichem_id),
    cichem_id != "",
    !is.na(aed_bp),
    is.finite(aed_bp),
    aed_bp > 0
  ) %>%
  group_by(cichem_id) %>%
  summarise(
    aed_bp = min(aed_bp),
    .groups = "drop"
  )

aed_by_chem <- full_join(
  mc7_aed_by_chem,
  bp_aed_by_chem,
  by = "cichem_id"
) %>%
  mutate(
    final_aed = case_when(
      !is.na(aed_mc7) & !is.na(aed_bp) ~
        pmin(aed_mc7, aed_bp),

      !is.na(aed_mc7) ~
        aed_mc7,

      !is.na(aed_bp) ~
        aed_bp,

      TRUE ~
        NA_real_
    ),

    final_aed_source = case_when(
      !is.na(aed_mc7) &
        !is.na(aed_bp) &
        aed_mc7 < aed_bp ~
        "B_MC7 p5",

      !is.na(aed_mc7) &
        !is.na(aed_bp) &
        aed_bp < aed_mc7 ~
        "Broad profiling minimum POD",

      !is.na(aed_mc7) &
        !is.na(aed_bp) &
        aed_mc7 == aed_bp ~
        "Equal in both sources",

      !is.na(aed_mc7) ~
        "B_MC7 p5 only",

      !is.na(aed_bp) ~
        "Broad profiling minimum POD only",

      TRUE ~
        NA_character_
    )
  ) %>%
  filter(
    !is.na(final_aed),
    is.finite(final_aed),
    final_aed > 0
  )

#========================================================
# Overlapping AED + HED chemicals
#========================================================

# USE eBMDHED VALUES

hed_overlap <- hed_p_by_chem_ebmdhed %>%
  filter(
    !is.na(hed_p05),
    hed_p05 > 0,
    n_hed_values > 0
  ) %>%
  select(
    cichem_id,
    hed_p05
  )

aed_overlap <- aed_by_chem %>%
  select(
    cichem_id,
    aed_mc7,
    aed_bp,
    final_aed,
    final_aed_source
  )

overlap_chems <- inner_join(
  hed_overlap,
  aed_overlap,
  by = "cichem_id"
) %>%
  semi_join(
    rocks_ids,
    by = "cichem_id"
  )

cat(
  "ROCKS overlap chemicals:",
  scales::comma(
    n_distinct(overlap_chems$cichem_id)
  ),
  "\n"
)

#========================================================
# HED
#========================================================

hed_plot_df <- overlap_chems %>%
  transmute(
    cichem_id,
    metric = "p05 eBMDHED",
    value = hed_p05,
    log10_value = log10(hed_p05)
  )

#========================================================
# AED
#========================================================

aed_plot_df <- overlap_chems %>%
  transmute(
    cichem_id,
    metric = "Minimum AED",
    value = final_aed,
    log10_value = log10(final_aed)
  )

#========================================================
# Combine
#========================================================

plot_df <- bind_rows(
  hed_plot_df,
  aed_plot_df
)

#========================================================
# Join ClassyFire class
#========================================================

plot_df <- plot_df %>%
  left_join(
    chem_class %>%
      select(
        cichem_id,
        class
      ) %>%
      distinct(),
    by = "cichem_id"
  ) %>%
  mutate(
    class = case_when(
      is.na(class) |
        str_trim(class) == "" ~ "Unclassified",
      TRUE ~ class
    )
  )

#========================================================
# Keep well-populated classes
#========================================================

# Classes retained if n >= 20 overlapping AED/HED chemicals
min_n_class <- 20

class_counts <- plot_df %>%
  distinct(
    cichem_id,
    class
  ) %>%
  count(
    class,
    name = "n_chem"
  ) %>%
  filter(
    class != "Unclassified",
    n_chem >= min_n_class
  )

plot_df_top <- plot_df %>%
  semi_join(
    class_counts,
    by = "class"
  )

#========================================================
# Labels
#========================================================

plot_df_top <- plot_df_top %>%
  left_join(
    class_counts,
    by = "class"
  ) %>%
  mutate(
    class_label = paste0(
  str_wrap(class, 20),
  "\n(n=",
  comma(n_chem),
  ")"
))


#========================================================
# Assessment values
#========================================================

tam_assess <- readRDS(
  here::here("data", "derived", "tam_assess.rds")
)

assess_class <- tam_assess %>%
  filter(
    population == "General Population"
  ) %>%
  filter(
    toxval_units=="mg/kg-day"
  )%>%
  filter(
    !is.na(toxval_numeric),
    toxval_numeric > 0
  ) %>%
  select(
    cichem_id,
    toxval_numeric
  ) %>%
  semi_join(
    rocks_ids,
    by = "cichem_id"
  ) %>%

  left_join(
    chem_class %>%
      select(
        cichem_id,
        class
      ) %>%
      distinct(),
    by = "cichem_id"
  ) %>%
  mutate(
    class = case_when(
      is.na(class) |
        str_trim(class) == "" ~ "Unclassified",
      TRUE ~ class
    ),
    log10_value = log10(toxval_numeric)
  ) %>%
  semi_join(
    class_counts,
    by = "class"
  )

#========================================================
# Assessment interval summary
#========================================================

assess_summary <- assess_class %>%
  group_by(class) %>%
  summarise(
    assess_n = n_distinct(cichem_id),
    assess_p25 = quantile(
      log10_value,
      probs = 0.25,
      na.rm = TRUE
    ),
    assess_median = median(
      log10_value,
      na.rm = TRUE
    ),
    assess_p75 = quantile(
      log10_value,
      probs = 0.75,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  filter(
    assess_n >= 3
  )

#========================================================
# Order classes by assessment median
#========================================================

class_order <- assess_summary %>%
  arrange(assess_median) %>%
  pull(class)

# append classes lacking sufficient assessment values
class_order <- c(
  class_order,
  setdiff(
    unique(plot_df_top$class),
    class_order
  )
)

#========================================================
# Apply ordering
#========================================================

plot_df_top <- plot_df_top %>%
  mutate(
    class = factor(
      class,
      levels = class_order
    )
  )

assess_summary <- assess_summary %>%
  mutate(
    class = factor(
      class,
      levels = class_order
    )
  )

#========================================================
# Create labels AFTER ordering
#========================================================

plot_df_top <- plot_df_top %>%
  mutate(
    class_label = forcats::fct_reorder(
      paste0(
        class,
        "\n(n=",
        scales::comma(n_chem),
        ")"
      ),
      as.numeric(class)
    )
  )

assess_summary <- assess_summary %>%
  left_join(
    class_counts,
    by = "class"
  ) %>%
  mutate(
    class_label = factor(
      paste0(
        class,
        "\n(n=",
        scales::comma(n_chem),
        ")"
      ),
      levels = levels(plot_df_top$class_label)
    )
  )

#========================================================
# Plot
#========================================================

p_class_box <- ggplot(
  plot_df_top,
  aes(
    x = log10_value,
    y = class_label,
    fill = metric
  )
) +

  geom_boxplot(
    position = position_dodge(width = 0.85),
    width = 0.65,
    outlier.alpha = 0.1
  ) +

  # assessment IQR
  geom_errorbarh(
    data = assess_summary,
    aes(
      xmin = assess_p25,
      xmax = assess_p75,
      y = class_label
    ),
    inherit.aes = FALSE,
    height = 0,
    linewidth = 1,
    color = "black",
    position = position_nudge(y = -0.35)
  ) +

  # assessment median
  geom_point(
    data = assess_summary,
    aes(
      x = assess_median,
      y = class_label
    ),
    inherit.aes = FALSE,
    shape = 21,
    fill = "white",
    color = "black",
    stroke = 0.8,
    size = 2.5,
    position = position_nudge(y = -0.35)
  ) +

  scale_fill_manual(
    values = c(
      "p05 eBMDHED" = "#BCE4F7",
      "Minimum AED" = "#0A32FF"
    )
  ) +

  labs(
    x = expression(
      log[10] *
      "(chemical-level potency)"
    ),
    y = NULL,
    fill = NULL
  ) +
scale_x_continuous(
  breaks = seq(-6, 5, by = 1),
  minor_breaks = seq(-6, 5, by = 0.5)
)+

  theme_minimal(base_size=16) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 14)
  )+
coord_flip(clip = "off") +

theme(
  legend.position = "bottom",
  legend.text = element_text(size = 16),
  legend.key.size = unit(1.1,"cm"),
  panel.grid.minor = element_blank(),
  axis.text.y = element_text(size = 14),
  axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    size = 14
  ),
  plot.margin = margin(10, 20, 20, 80)
)

p_class_box

ggsave(
    file.path(out_dir,
    "fig03_potency_by_class.png"
  ),
  plot = p_class_box,
  width = 14,
  height = 9
)

#========================================================
# Compare class-level potency rankings across data streams
#========================================================
#
# Purpose:
#   Quantify whether chemical classes have similar relative potency rankings
#   across Minimum AED, p05 eBMDHED, and assessment value distributions.
#
# Notes:
#   - Rankings are based on class-level median log10 potency.
#   - Lower values indicate more potent/lower dose values.
#   - Spearman rho is used because the question is whether classes are ordered
#     similarly across data streams, not whether the medians are numerically equal.
#========================================================

# AED and eBMDHED class medians
class_summary_aed_hed <- plot_df_top %>%
  group_by(class, metric) %>%
  summarise(
    n_chem = n_distinct(cichem_id),
    median_log10_value = median(log10_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = metric,
    values_from = c(n_chem, median_log10_value),
    names_sep = "_"
  ) %>%
  rename(
    n_aed = `n_chem_Minimum AED`,
    n_hed = `n_chem_p05 eBMDHED`,
    median_aed = `median_log10_value_Minimum AED`,
    median_hed = `median_log10_value_p05 eBMDHED`
  )

# Assessment value class medians
class_summary_assess <- assess_summary %>%
  transmute(
    class = as.character(class),
    n_assess = assess_n,
    median_assess = assess_median
  )

# Combine class-level medians
class_rank_compare <- class_summary_aed_hed %>%
  mutate(
    class = as.character(class)
  ) %>%
  left_join(
    class_summary_assess,
    by = "class"
  ) %>%
  mutate(
    rank_aed = rank(median_aed, ties.method = "average"),
    rank_hed = rank(median_hed, ties.method = "average"),
    rank_assess = rank(median_assess, ties.method = "average", na.last = "keep")
  )

# Inspect
print(class_rank_compare)

#--------------------------------------------------------
# Spearman rank correlations
#--------------------------------------------------------

class_rank_correlations <- tibble(
  comparison = c(
    "Minimum AED vs p05 eBMDHED",
    "Minimum AED vs assessment values",
    "p05 eBMDHED vs assessment values"
  ),
  n_classes = c(
    sum(complete.cases(class_rank_compare$median_aed,
                       class_rank_compare$median_hed)),
    sum(complete.cases(class_rank_compare$median_aed,
                       class_rank_compare$median_assess)),
    sum(complete.cases(class_rank_compare$median_hed,
                       class_rank_compare$median_assess))
  ),
  spearman_rho = c(
    cor(
      class_rank_compare$median_aed,
      class_rank_compare$median_hed,
      method = "spearman",
      use = "complete.obs"
    ),
    cor(
      class_rank_compare$median_aed,
      class_rank_compare$median_assess,
      method = "spearman",
      use = "complete.obs"
    ),
    cor(
      class_rank_compare$median_hed,
      class_rank_compare$median_assess,
      method = "spearman",
      use = "complete.obs"
    )
  )
)

print(class_rank_correlations)

## COMBINE

library(png)
library(grid)
library(cowplot)

#--------------------------------------------------
# Read images
#--------------------------------------------------

img_a <- rasterGrob(
  readPNG(
    here::here(
      "outputs",
      "figures",
      "fig03a_data_coverage.png"
    )
  ),
  interpolate = TRUE
)

img_b <- rasterGrob(
  readPNG(
    here::here(
      "outputs",
      "figures",
      "fig03b_potency_distribution_overlap.png"
    )
  ),
  interpolate = TRUE
)

img_c <- rasterGrob(
  readPNG(
    here::here(
      "outputs",
      "figures",
      "fig03_potency_by_class.png"
    )
  ),
  interpolate = TRUE
)


#--------------------------------------------------
# Convert to ggplot objects
#--------------------------------------------------

p_a <- ggdraw() +
  draw_plot(img_a)

p_b <- ggdraw() +
  draw_plot(img_b)

p_c <- ggdraw() +
  draw_plot(img_c)

top_row <- plot_grid(
  p_a,
  p_b,
  ncol = 2,
  rel_widths = c(1, 1),
  labels = c("A", "B"),
  label_fontface = "bold",
  label_size = 14,
  label_x = 0.02,
  label_y = 0.98,
  hjust = 0,
  vjust = 1
)

fig3 <- plot_grid(
  top_row,
  p_c,
  ncol = 1,
  rel_heights = c(1, 2),
  labels = c("", "C"),
  label_fontface = "bold",
  label_size = 14,
  label_x = c(0, 0.02),
  label_y = c(1, 0.98),
  hjust = 0,
  vjust = 1
)


fig3

ggsave(
  file.path(
    out_dir,
    "fig03_combined.png"
  ),
  plot = fig3,
  width = 9,
  height = 9,
  dpi = 600
)

#========================================================
# Paired comparison of Minimum AED vs p05 eBMDHED
#
# Lower dose = greater potency
#
# Primary question:
# Are Minimum AEDs generally lower than p05 eBMDHEDs
# for the same chemicals?
#========================================================

potency_test_df <- potency_overlap %>%
  transmute(
    cichem_id,
    minimum_aed,
    ebmdhed_p05,

    log10_aed = log10(minimum_aed),
    log10_ebmdhed = log10(ebmdhed_p05),

    # Negative values mean AED is more potent
    log10_ratio = log10(minimum_aed / ebmdhed_p05),

    # Ratio < 1 means AED is more potent
    potency_ratio = minimum_aed / ebmdhed_p05,

    aed_more_potent = minimum_aed < ebmdhed_p05,
    equal_potency = minimum_aed == ebmdhed_p05
  ) %>%
  filter(
    is.finite(log10_aed),
    is.finite(log10_ebmdhed),
    is.finite(log10_ratio)
  )

# Confirm exactly one paired observation per chemical
stopifnot(
  nrow(potency_test_df) ==
    n_distinct(potency_test_df$cichem_id)
)

#--------------------------------------------------------
# Primary test: paired Wilcoxon signed-rank test
#--------------------------------------------------------

wilcox_potency <- wilcox.test(
  x = potency_test_df$log10_aed,
  y = potency_test_df$log10_ebmdhed,
  paired = TRUE,
  alternative = "less",
  exact = FALSE,
  conf.int = TRUE,
  conf.level = 0.95
)

print(wilcox_potency)

#========================================================
# Descriptive effect-size summary
#========================================================

potency_effect_summary <- potency_test_df %>%
  summarise(
    n_chem = n(),

    median_minimum_aed =
      median(minimum_aed),

    median_p05_ebmdhed =
      median(ebmdhed_p05),

    median_log10_ratio =
      median(log10_ratio),

    median_aed_to_ebmdhed_ratio =
      median(potency_ratio),

    median_fold_lower =
      1 / median(potency_ratio),

    n_aed_more_potent =
      sum(aed_more_potent),

    n_ebmdhed_more_potent =
      sum(minimum_aed > ebmdhed_p05),

    n_equal =
      sum(equal_potency),

    pct_aed_more_potent =
      100 * mean(aed_more_potent)
  )

print(potency_effect_summary)
