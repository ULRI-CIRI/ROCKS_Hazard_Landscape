# ===========================================================
# fig01_inventory_coverage.R
# ===========================================================

library(dplyr)
library(ggplot2)
library(here)
library(scales)
library(patchwork)
library(tidyr)
library(dplyr)
library(stringr)

out_dir <- here::here(
  "outputs",
  "figures"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# INPUT

tam <- readRDS(
  here("data", "derived", "tam.rds")
)

tam_assess <- readRDS(
  here("data", "derived", "tam_assess.rds")
)

rocks <- readRDS(
  here("data", "derived", "rocks.rds")
)

b_exec_sum <- readRDS(
  here("data", "derived", "b_exec_sum.rds")
)

bp <- readRDS(
  here("data", "derived", "bp.rds")
)

source(here("R", "ulri_theme.R"))

#----------------------------
# 1. Create unique chemical-level flags by cichem_id
#----------------------------

tam_mrsv_flag <- tam %>%
  filter(toxval_type_supercategory == "Mortality Response Summary Value") %>%
  distinct(cichem_id) %>%
  mutate(has_in_vivo_mrsv = TRUE)

tam_drsv_flag <- tam %>%
  filter(toxval_type_supercategory == "Dose Response Summary Value") %>%
  distinct(cichem_id) %>%
  mutate(has_in_vivo_drsv = TRUE)

ivt_flag <- b_exec_sum %>%
  distinct(cichem_id) %>%
  mutate(has_in_vitro = TRUE)

bp_flag <- bp %>%
  distinct(cichem_id) %>%
  mutate(has_broad_profiling = TRUE)

assess_flag <- tam_assess %>%
  filter(population=="General Population")%>%
  filter(toxval_type!="CTV") %>%
  distinct(cichem_id) %>%
  mutate(has_assessment = TRUE)

## Table of assessment values ##

toxval_by_source <- tam_assess %>%
  filter(
    population == "General Population",
    toxval_type != "CTV"
  ) %>%
  group_by(source) %>%
  summarise(
    toxval_types = str_c(sort(unique(toxval_type)), collapse = ", "),
    n_types = n_distinct(toxval_type),
    .groups = "drop"
  ) %>%
  arrange(source)

#----------------------------
# 2. Join flags onto one-row-per-chemical ROCKS inventory
#----------------------------

rocks_panelA <- rocks %>%
  distinct(cichem_id, .keep_all = TRUE) %>%
  left_join(tam_mrsv_flag, by = "cichem_id") %>%
  left_join(tam_drsv_flag, by = "cichem_id") %>%
  left_join(ivt_flag, by = "cichem_id") %>%
  left_join(bp_flag, by = "cichem_id") %>%
  left_join(assess_flag, by = "cichem_id") %>%
  mutate(
    has_in_vivo_mrsv = if_else(is.na(has_in_vivo_mrsv), FALSE, has_in_vivo_mrsv),
    has_in_vivo_drsv = if_else(is.na(has_in_vivo_drsv), FALSE, has_in_vivo_drsv),
    has_in_vitro = if_else(is.na(has_in_vitro), FALSE, has_in_vitro),
    has_broad_profiling = if_else(is.na(has_broad_profiling), FALSE, has_broad_profiling),
    has_assessment = if_else(is.na(has_assessment), FALSE, has_assessment)
  )

#----------------------------
# 3. Collapse to mutually exclusive inventory-wide coverage groups
#----------------------------

rocks_panelA <- rocks_panelA %>%
  mutate(
    has_in_vivo =  has_in_vivo_drsv, # dose response summary values only - this still allows acute duration studies
    has_any_data = has_in_vivo | has_in_vitro | has_broad_profiling,
    coverage_group = case_when(
      has_assessment ~ "Assessment value",
      has_any_data & !has_assessment ~ "Toxicity data only",
      TRUE ~ "No toxicity data or assessment value"
    )
  )

#----------------------------
# 4. Summarize counts and percentages
#----------------------------

panelA_df <- rocks_panelA %>%
  count(coverage_group, name = "n") %>%
  mutate(
    coverage_group = factor(
      coverage_group,
      levels = c(
        "Assessment value",
        "Toxicity data only",
        "No toxicity data or assessment value"
      )
    )
  ) %>%
  arrange(coverage_group) %>%
  mutate(
    pct = n / sum(n),
    label = paste0(comma(n), "\n(", percent(pct, accuracy = 0.1), ")")
  )

print(panelA_df)

# Save stats for manuscript

assessment_pct <- panelA_df %>%
  filter(coverage_group == "Assessment value") %>%
  pull(pct)

tox_pct <- panelA_df %>%
  filter(coverage_group == "Toxicity data only") %>%
  pull(pct)

tox_n <- panelA_df %>%
  filter(coverage_group == "Toxicity data only") %>%
  pull(n)

no_tox_n <- panelA_df %>%
  filter(coverage_group == "No toxicity data or assessment value") %>%
  pull(n)

no_tox_pct <- panelA_df %>%
  filter(coverage_group == "No toxicity data or assessment value") %>%
  pull(pct)

assessment_n <- panelA_df %>%
  filter(coverage_group == "Assessment value") %>%
  pull(n)

inventory_n <- sum(panelA_df$n)

manuscript_stats <- list(
  inventory_n = inventory_n,
  assessment_n = assessment_n,
  assessment_pct = assessment_pct,
  tox_n = tox_n,
  tox_pct = tox_pct,
  no_tox_n = no_tox_n,
  no_tox_pct = no_tox_pct
)

# compute cumulative top positions so labels can be placed above each segment
panelA_df <- panelA_df %>%
  arrange(coverage_group) %>%
  mutate(pos_top = cumsum(pct)) 

n_inventory <- sum(panelA_df$n)

print(n_inventory)

#----------------------------
# 5. Plot single horizontal stacked bar
#----------------------------

panelA_labels <- panelA_df %>%
  mutate(
    coverage_group = factor(
      coverage_group,
      levels = c(
        "No toxicity data or assessment value",
        "Toxicity data only",
        "Assessment value"
      )
    )
  ) %>%
  arrange(coverage_group) %>%
  mutate(
    cumulative = cumsum(pct),
    mid_point = cumulative - pct/2,
    outside = pct < 0.05
  )
p_panelA <- ggplot(panelA_df, aes(x = 1, y = pct, fill = coverage_group)) +
  geom_col(width = 0.6, color = "white", linewidth = 0.5) +
  # labels inside larger segments
  geom_text(
    data = panelA_labels %>% filter(!outside),
    aes(
      y = mid_point,
      label = percent(pct, accuracy = 0.1)
    ),
    size = 6
  ) +

  # labels outside for tiny segments
  geom_text(
    data = panelA_labels %>% filter(outside),
    aes(
      y = 1.01,   # just beyond 100%
      label = percent(pct, accuracy = 0.1)
    ),
    hjust = 0,
    size = 6
  ) +
  coord_flip() +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.15))
  ) +
  scale_fill_manual(
    values = c(
      "Assessment value" = "#0A32FF",               
      "Toxicity data only" = "#BCE4F7",     
      "No toxicity data or assessment value" = "#E5DDCA"),
    drop = FALSE,
    guide = guide_legend(reverse=TRUE)
  ) +
  labs(
    title = paste0("ROCKS Chemical Collection (n = ", comma(n_inventory), ")"),
    x = NULL,
    y = "Percent of ROCKS chemicals",
    fill = NULL
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(face = "bold", size = 20),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  panel.grid.major.x = element_blank(),
  legend.position = "bottom",
  legend.text = element_text(size = 16),
  legend.key.size = unit(1.2, "cm")
)
p_panelA

#----------------------------
# 6. Save figure
#----------------------------

ggsave(
  file.path(out_dir, 
    "fig01_panelA_inventory_coverage_summary.png"
  ),
  plot = p_panelA,
  width = 12,
  height = 3,
  units = "in",
  dpi = 300

)

saveRDS(
  panelA_df,
  here(
    "data",
    "derived",
    "fig01_panelA_df.rds"
  )
)

message(
  "fig01_panelA_inventory_coverage.R completed successfully"
)

message(
  "Inventory chemicals: ",
  comma(n_inventory)
)

message(
  "Figure saved to: figures/fig01_panelA_inventory_coverage_summary.png"
)

# Quantify multiple inventories #

# ============================================================
# Coverage among chemicals in Multiple Inventories
# ============================================================

multiple_inventory_df <- rocks_panelA %>%
  filter(
    sector_class == "Multiple Inventories"
  ) %>%
  mutate(
    coverage_group = factor(
      coverage_group,
      levels = c(
        "Assessment value",
        "Toxicity data only",
        "No toxicity data or assessment value"
      )
    )
  )

n_multiple_inventory <- n_distinct(
  multiple_inventory_df$cichem_id
)

multiple_inventory_summary <- multiple_inventory_df %>%
  count(
    coverage_group,
    name = "n",
    .drop = FALSE
  ) %>%
  mutate(
    n_inventory = n_multiple_inventory,
    pct = n / n_inventory,
    pct_label = percent(
      pct,
      accuracy = 0.1
    ),
    label = paste0(
      comma(n),
      "\n(",
      pct_label,
      ")"
    )
  ) %>%
  arrange(coverage_group)

print(multiple_inventory_summary)


######################################
# PANEL B
# EVIDENCE STREAM OVERLAP
######################################

#----------------------------
# 1. Prep flags
#----------------------------
df <- rocks_panelA %>%
  mutate(
    has_in_vivo = has_in_vivo_drsv,   # key change
    assessment_status = if_else(has_assessment, "Assessment", "No assessment")
  )

#----------------------------
# 2. Evidence combinations
#----------------------------
df_combo <- df %>%
  mutate(
    has_in_vivo = coalesce(has_in_vivo, FALSE),
    has_in_vitro = coalesce(has_in_vitro, FALSE),
    has_broad_profiling = coalesce(has_broad_profiling, FALSE),

    combo = paste0(
      ifelse(has_in_vivo, "in_vivo + ", ""),
      ifelse(has_in_vitro, "in_vitro + ", ""),
      ifelse(has_broad_profiling, "broad", "")
    ),
    combo = gsub(" \\+ $", "", combo),
    combo = ifelse(combo == "", "no data", combo)
  ) %>%
  filter(combo != "no data")

df_combo <- df_combo %>%
  mutate(combo = recode(combo,
    in_vivo                     = "In vivo",
    in_vitro                    = "Targeted bioactivity",
    broad                       = "Broad profiling",
    `in_vivo + in_vitro`        = "In vivo + Targeted bioactivity",
    `in_vivo + broad`           = "In vivo + Broad profiling",
    `in_vitro + broad`          = "Targeted bioactivity + Broad profiling",
    `in_vivo + in_vitro + broad`= "In vivo + Targeted bioactivity + Broad profiling"
  ))

df_combo <- df_combo %>%
  add_count(combo, name = "combo_n") %>%
  filter(combo_n >= 20)


#----------------------------
# 3. Counts for bars
#----------------------------
bar_df <- df_combo %>%
  count(combo, assessment_status, name = "n")

totals_df <- bar_df %>%
  group_by(combo) %>%
  summarise(
    total = sum(n),
    .groups = "drop"
  ) %>%
  mutate(
    pct_inventory = total / inventory_n,
    pct_label = percent(pct_inventory, accuracy = 0.1)
  ) %>%
  arrange(desc(total))

label_df <- totals_df %>%
  left_join(
    bar_df %>%
      group_by(combo) %>%
      summarise(
        y_pos = sum(n) / 2,
        .groups = "drop"
      ),
    by = "combo"
  )

bar_df$combo <- factor(bar_df$combo, levels = totals_df$combo)
totals_df$combo <- factor(totals_df$combo, levels = totals_df$combo)

#----------------------------
# 4. Matrix (dot panel)
#----------------------------
matrix_df <- df_combo %>%
  distinct(combo, has_in_vivo, has_in_vitro, has_broad_profiling) %>%
  pivot_longer(
    cols = c(has_in_vivo, has_in_vitro, has_broad_profiling),
    names_to = "stream",
    values_to = "present"
  )

matrix_df$combo <- factor(matrix_df$combo, levels = levels(bar_df$combo))

#----------------------------
# 5. Top bar plot
#----------------------------
p_bar <- ggplot(
  bar_df,
  aes(x = combo, y = n, fill = assessment_status)
) +
  geom_col(
    color = "white",
    linewidth = 0.4
  ) +
geom_text(
  data = label_df,
  aes(
    x = combo,
    y = y_pos,
    label = percent(
      pct_inventory,
      accuracy = 0.1
    )
  ),
  inherit.aes = FALSE,
  color = "black",
  size = 5) +
  scale_fill_manual(
    values = c(
      "Assessment" = "#0A32FF",
      "No assessment" = "#BCE4F7"
    )
  ) +
  labs(
    y = "Number of chemicals",
    x = NULL,
    fill = NULL
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

#----------------------------
# 6. Dot matrix plot
#----------------------------
p_matrix <- ggplot(matrix_df,
                   aes(x = combo, y = stream)) +
  geom_line(
    data = subset(matrix_df, present),
    aes(group = combo),
    color = "black",
    linewidth = 1
  ) +
  geom_point(
    aes(color = present),
    size = 4
  ) +
  scale_color_manual(
    values = c(
      "TRUE" = "black",
      "FALSE" = "grey80"
    ),
    guide = "none"
  ) +
  scale_y_discrete(
    labels = c(
      has_in_vivo = "In vivo",
      has_in_vitro = "Targeted bioactivity",
      has_broad_profiling = "Broad profiling"
    )
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 18) +
theme(
  axis.text.x = element_blank(),
  panel.grid.minor = element_blank(),
  legend.position = "none"
)

#######################
# PANEL C #
#######################

#----------------------------
# 1. Define sector indicator columns in rocks
#    
#----------------------------

sector_cols <- c(
  "consumer",
  "industrial",
  "agricultural",
  "food",
  "pharmaceutical",
  "biomonitoring",
  "env_monitoring"
)

#----------------------------
# 2. Create unique chemical-level flags by cichem_id
#----------------------------

tam_flag <- tam %>%
  distinct(cichem_id) %>%
  mutate(has_in_vivo = TRUE)

ivt_flag <- b_exec_sum %>%
  distinct(cichem_id) %>%
  mutate(has_in_vitro = TRUE)

bp_flag <- bp %>%
  distinct(cichem_id) %>%
  mutate(has_broad_profiling = TRUE)

assess_flag <- tam_assess %>%
  filter(population=="General Population")%>%
  distinct(cichem_id) %>%
  mutate(has_assessment = TRUE)

#----------------------------
# 3. Join flags onto one-row-per-chemical ROCKS inventory
#----------------------------

rocks_sector_panel <- rocks %>%
  distinct(cichem_id, .keep_all = TRUE) %>%
  left_join(tam_flag, by = "cichem_id") %>%
  left_join(ivt_flag, by = "cichem_id") %>%
  left_join(bp_flag, by = "cichem_id") %>%
  left_join(assess_flag, by = "cichem_id") %>%
  mutate(
    has_in_vivo = if_else(is.na(has_in_vivo), FALSE, has_in_vivo),
    has_in_vitro = if_else(is.na(has_in_vitro), FALSE, has_in_vitro),
    has_broad_profiling = if_else(is.na(has_broad_profiling), FALSE, has_broad_profiling),
    has_assessment = if_else(is.na(has_assessment), FALSE, has_assessment),
    has_any_data = has_in_vivo | has_in_vitro | has_broad_profiling,
    coverage_group = case_when(
      has_assessment ~ "Assessment value present",
      has_any_data & !has_assessment ~ "Toxicity data, no assessment value",
      TRUE ~ "No toxicity data"
    )
  )

#----------------------------
# 4. Pivot sector columns longer and coerce membership to logical
#----------------------------

primary_sector_cols <- c(
  "consumer",
  "industrial",
  "agricultural",
  "food",
  "pharmaceutical"
)

fallback_sector_cols <- c(
  "biomonitoring",
  "env_monitoring"
)

sector_long <- rocks_sector_panel %>%
  select(cichem_id, coverage_group, all_of(sector_cols)) %>%
  mutate(
    # Convert the sector fields to explicit logical flags.
    # coalesce() prevents NA values from propagating.
    across(
      all_of(sector_cols),
      ~ coalesce(.x == TRUE | .x == 1, FALSE)
    ),

    has_primary_sector = if_any(
      all_of(primary_sector_cols),
      ~ .x
    ),

    has_biomonitoring = biomonitoring,
    has_env_monitoring = env_monitoring,

    # Fallback categories apply only when there is no primary sector.
    biomonitoring_unknown =
      !has_primary_sector & has_biomonitoring,

    env_monitoring_unknown =
      !has_primary_sector & has_env_monitoring,

    # True unknown means no primary or fallback annotation.
    unknown =
      !has_primary_sector &
      !has_biomonitoring &
      !has_env_monitoring
  ) %>%
  pivot_longer(
    cols = c(
      all_of(primary_sector_cols),
      biomonitoring_unknown,
      env_monitoring_unknown,
      unknown
    ),
    names_to = "sector",
    values_to = "in_sector_raw"
  ) %>%
  filter(in_sector_raw)


#----------------------------
# 5. Optional prettier sector labels
#----------------------------

sector_long <- sector_long %>%
  mutate(
sector = recode(
  sector,
  consumer = "Consumer",
  industrial = "Industrial",
  agricultural = "Agricultural",
  food = "Food",
  pharmaceutical = "Pharmaceutical",
  biomonitoring_unknown = "Unknown\nBiomonitoring",
  env_monitoring_unknown = "Unknown\nEnv. monitoring",
  unknown = "Unknown"
)
  )


#----------------------------
# 6. Summarize counts by sector and coverage group
#----------------------------

sector_counts <- sector_long %>%
  count(sector, coverage_group, name = "n")

sector_order <- c(
  "Pharmaceutical",
  "Agricultural",
  "Food",
  "Consumer",
  "Industrial",
  "Unknown\nBiomonitoring",
  "Unknown\nEnv. monitoring",
  "Unknown"
)

sector_totals <- sector_long %>%
  count(sector, name = "total_n")

sector_totals <- sector_totals %>%
  mutate(
    sector = factor(sector, levels = sector_order)
  ) %>%
  arrange(sector) %>%
  mutate(
    sector_label = paste0(sector, " (n=", comma(total_n), ")")
  )

sector_counts <- sector_counts %>%
  left_join(sector_totals, by = "sector") %>%
  mutate(
    sector_label = factor(sector_label, levels = rev(sector_totals$sector_label)),
    coverage_group = factor(
      coverage_group,
      levels = c(
        "Assessment value present",
        "Toxicity data, no assessment value",
        "No toxicity data"
      )
    )
  )

#----------------------------
# 7. Compute within-sector percentages
#----------------------------

sector_plot_df <- sector_counts %>%
  group_by(sector_label) %>%
  mutate(
    pct = n / sum(n),

    # Light blue: always label
    # Beige: label only when >= 10%
    # Dark blue: handled separately outside the bar
    pct_label = case_when(
      coverage_group == "Toxicity data, no assessment value" ~
        percent(pct, accuracy = 1),

      coverage_group == "No toxicity data" & pct >= 0.1 ~
        percent(pct, accuracy = 1),

      TRUE ~ ""
    )
  ) %>%
  ungroup()


# Labels for ALL assessment values, placed to the right
assessment_labels <- sector_plot_df %>%
  filter(coverage_group == "Assessment value present") %>%
  mutate(
    label = percent(pct, accuracy = 1)
  )

#test
assessment_labels <- sector_plot_df %>%
  filter(coverage_group == "Assessment value present") %>%
  mutate(
    label_y = c(
      1.03,  # Ag
      1.03,  # 
      1.03,  # Food
      1.03,  # Consumer
      1.03,  # Pharm
      1.07,  # Biomonitoring
      1.06   # Env monitoring
    )
  )%>%
  mutate(
    label = percent(pct, accuracy = 1)
  )

p_sector <- ggplot(
  sector_plot_df,
  aes(x = sector_label, y = pct, fill = coverage_group)
) +
  geom_col(
    width = 0.75,
    color = "white",
    linewidth = 0.4
  ) +

  # Labels inside the bars:
  # all light-blue segments and beige segments >= 10%
  geom_text(
    aes(
      label = pct_label,
      color = coverage_group
    ),
    position = position_stack(vjust = 0.5),
    size = 5,
    lineheight = 0.95
  ) +

geom_text(
  data = assessment_labels,
  aes(
    x = sector_label,
    y = label_y,
    label = label
  ),
  inherit.aes = FALSE,
  hjust = 0,
  size = 5,
  color = "black"
)+

  coord_flip(clip = "off") +

  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1.15),
    expand = expansion(mult = c(0, 0))
  ) +

  scale_fill_manual(
    values = c(
      "Assessment value present" = "#0A32FF",
      "Toxicity data, no assessment value" = "#BCE4F7",
      "No toxicity data" = "#E5DDCA"
    ),
    drop = FALSE
  ) +

  scale_color_manual(
    values = c(
      "Assessment value present" = "white",
      "Toxicity data, no assessment value" = "black",
      "No toxicity data" = "black"
    ),
    guide = "none"
  ) +

  labs(
    x = NULL,
    y = "Percent of chemicals within sector"
  ) +

  theme_minimal(base_size = 16) +
  theme(
    axis.text.y = element_text(size = 14),
    legend.position = "none",
    panel.grid.minor = element_blank(),
      plot.margin = margin(5.5, 45, 5.5, 5.5)
  )

print(p_sector)

#----------------------------
# 9. Save figure
#----------------------------

ggsave(
  filename = here("outputs","figures", "fig1_panel_sector_coverage.png"),
  plot = p_sector,
  width = 6,
  height = 6,
  dpi = 300
)


#----------------------------
# 7. Combine (patchwork)
#----------------------------

p_upset <- p_bar / p_matrix +
  plot_layout(heights = c(2, 1))

print(p_upset)

ggsave(
  file.path(out_dir, 
    "fig01_panelB_evidence_overlap.png"
  ),
  plot = p_upset,
  width = 6,
  height = 6,
  dpi = 300
)

saveRDS(
  bar_df,
  here(
    "data",
    "derived",
    "fig01_panelB_bar_df.rds"
  )
)

saveRDS(
  matrix_df,
  here(
    "data",
    "derived",
    "fig01_panelB_matrix_df.rds"
  )
)

#################################
# Combine all Figure 1 panels #
#################################

library(cowplot)
library(png)
library(grid)

img_top <- rasterGrob(
  readPNG(
    here(
      "outputs",
      "figures",
      "fig01_panelA_inventory_coverage_summary.png"
    )
  ),
  interpolate = TRUE
)

img_left <- rasterGrob(
  readPNG(
    here(
      "outputs",
      "figures",
      "fig01_panelB_evidence_overlap.png"
    )
  ),
  interpolate = TRUE
)

img_right <- rasterGrob(
  readPNG(
    here(
      "outputs",
      "figures",
      "fig1_panel_sector_coverage.png"
    )
  ),
  interpolate = TRUE
)

p_top <- ggdraw() + draw_grob(img_top)
p_left <- ggdraw() + draw_grob(img_left)
p_right <- ggdraw() + draw_grob(img_right)

bottom_row <- plot_grid(
  p_left,
  p_right,
  ncol = 2,
  rel_widths = c(1, 1),
  labels = c("B", "C"),
  label_fontface = "bold",
  label_size = 16
)

final_fig <- plot_grid(
  p_top,
  bottom_row,
  ncol = 1,
  rel_heights = c(0.4, 1),
  labels = "A",
  label_fontface = "bold",
  label_size = 16
)

final_fig

ggsave(
  file.path(out_dir, 
  "combined.png"),
  final_fig,
  width = 10,
  height = 7,
  dpi = 600
)
