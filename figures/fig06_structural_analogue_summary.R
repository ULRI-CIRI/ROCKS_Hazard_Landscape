# ===========================================================
# fig06_structural_analogue_summary.R
#
# Purpose:
#   Create a 3-panel summary figure for structural analogue
#   coverage among data-poor ROCKS chemicals.
#
# Panels:
#   A. ROCKS inventory fingerprint/data denominator
#   B. Coverage across Tanimoto thresholds
#   C. Number of source analogues per covered target at primary threshold
# ===========================================================


# ===========================================================
# 0) Libraries
# ===========================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(here)
library(scales)
library(patchwork)
library(readr)
library(stringr)
library(ggrepel)


# ===========================================================
# 1) Paths
# ===========================================================

derived_dir <- here::here("data", "derived")
fig_dir <- here::here("outputs", "figures")

dir.create(
  fig_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ===========================================================
# 2) Load inputs
# ===========================================================

fp_all <- readRDS(
  file.path(derived_dir, "fp_all_anydata_rocks_flag.rds")
)

coverage_by_threshold <- readRDS(
  file.path(derived_dir, "coverage_by_threshold.rds")
)

target_match_counts <- readRDS(
  file.path(derived_dir, "target_match_counts_thr070.rds")
)

summary_now <- readRDS(
  file.path(derived_dir, "summary_now_thr070.rds")
)


# ===========================================================
# 3) Settings
# ===========================================================

primary_threshold <- 0.7

base_family <- "sans"

panel_label_size <- 12
axis_title_size <- 10
axis_text_size <- 9
legend_text_size <- 9
legend_title_size <- 10

col_fingerprint_missing <- "#E5DDCA"
col_has_data <- "#0A32FF"
col_target <- "#BCE4F7"
col_coverage <- "#123A28"
col_match_dist <- "#BCE4F7"

# ===========================================================
# 4) Build Panel A data
# ===========================================================
# Panel A summarizes the ROCKS denominator hierarchy:
#   total ROCKS
#   ROCKS without fingerprints
#   fingerprinted ROCKS with selected data
#   fingerprinted ROCKS without selected data, i.e., targets

rocks_total <- fp_all %>%
  filter(on_rocks == 1L) %>%
  summarise(
    n_rocks_with_fp = n_distinct(cichem_id),
    n_rocks_with_any_data = n_distinct(cichem_id[has_any_data == 1L]),
    n_rocks_without_any_data = n_distinct(cichem_id[has_any_data == 0L])
  )

fingerprint_coverage <- tibble::tibble(
  n_rocks_total = 238117L,
  n_rocks_with_fp = rocks_total$n_rocks_with_fp,
  n_rocks_without_fp = 238117L - rocks_total$n_rocks_with_fp
) %>%
  mutate(
    pct_rocks_with_fp = 100 * n_rocks_with_fp / n_rocks_total,
    pct_rocks_without_fp = 100 * n_rocks_without_fp / n_rocks_total
  )

panel_a_dat <- tibble::tibble(
  status = c(
    "No Morgan fingerprint",
    "Fingerprint + toxicity data",
    "Fingerprint + no toxicity data"
  ),
  n_chemicals = c(
    fingerprint_coverage$n_rocks_without_fp,
    rocks_total$n_rocks_with_any_data,
    rocks_total$n_rocks_without_any_data
  )
) %>%
  mutate(
    pct_total_rocks = 100 * n_chemicals / fingerprint_coverage$n_rocks_total,
    status = factor(
      status,
      levels = c(
        "Fingerprint + no toxicity data",
        "Fingerprint + toxicity data",
        "No Morgan fingerprint"
      )
    ),
    label = paste0(
      comma(n_chemicals),
      "\n",
      sprintf("%.1f%%", pct_total_rocks)
    )
  )

# ===========================================================
# 5) Panel A
# ===========================================================

panel_a_dat <- panel_a_dat %>%
  mutate(
    ymax = cumsum(n_chemicals),
    ymin = lag(ymax, default = 0),
    midpoint = (ymin + ymax) / 2
  )
panel_a <- ggplot(
  panel_a_dat,
  aes(
    x = "ROCKS inventory",
    y = n_chemicals,
    fill = status
  )
) +
  geom_col(
    width = 0.55,
    color = "white",
    linewidth = 0.3
  ) +
geom_text(
  data = filter(panel_a_dat,
                status != "Fingerprint + toxicity data"),
  aes(
    label = label,
    color = status
  ),
  position = position_stack(vjust = 0.5),
  size = 6,
  lineheight = 0.9
) +
geom_text_repel(
  data = filter(panel_a_dat,
                status == "Fingerprint + toxicity data"),
  aes(
    x = "ROCKS inventory",
    y = midpoint,
    label = label
  ),
  color = "black",
  size = 6,
  direction = "x",
  nudge_x = -1,
  segment.color = "black",
  segment.size = 0.4,
  box.padding = 0.5,
  point.padding = 0.2,
  min.segment.length = 0,
  force = 5
)+
  scale_color_manual(
    values = c(
      "No Morgan fingerprint" = "black",
      "Fingerprint + toxicity data" = "white",
      "Fingerprint + no toxicity data" = "black"
    ),
    guide = "none"
  ) +
  scale_y_continuous(
    labels = comma,
    expand = expansion(mult = c(0, 0.03))
  ) +
  scale_fill_manual(
    values = c(
      "No Morgan fingerprint" = col_fingerprint_missing,
      "Fingerprint + toxicity data" = col_has_data,
      "Fingerprint + no toxicity data" = col_target
    ),
    name = NULL,
    guide = guide_legend(reverse= TRUE)
  ) +
  labs(
    x = NULL,
    y = "Number of chemicals" )+
  coord_flip()+
  theme_minimal(base_size = 16) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 16),
  legend.key.size = unit(1.2, "cm")
  )
panel_a

ggsave(
  file.path(fig_dir, 
    "fig06_panelA_denominator.png"
  ),
  plot = panel_a,
  width = 12,
  height = 3,
  units = "in",
  dpi = 300

)

# ===========================================================
# 6) Build Panel B data
# ===========================================================

panel_b_dat <- coverage_by_threshold %>%
  mutate(
    threshold = as.numeric(threshold),
    threshold_label = ifelse(
      threshold == primary_threshold,
      paste0(threshold, " primary"),
      as.character(threshold)
    )
  )


# ===========================================================
# 7) Panel B
# ===========================================================

panel_b <- ggplot(
  panel_b_dat,
  aes(
    x = threshold,
    y = pct_targets_with_match
  )
) +
  geom_line(
    color = col_coverage,
    linewidth = 0.8
  ) +
  geom_point(
    color = col_coverage,
    size = 6
  ) +
  geom_vline(
    xintercept = primary_threshold,
    linetype = "dashed",
    color = "gray45",
    linewidth = 0.4
  ) +
#  geom_text(
#    data = panel_b_dat %>%
#      filter(threshold == primary_threshold),
#    aes(
#      label = paste0(
#        "Primary threshold\n",
#        sprintf("%.1f%%", pct_targets_with_match)
#      )
#    ),
#    nudge_x = 0.035,
#    hjust = 0,
#    size = 5,
#    color = "gray20"
#  ) +
  scale_x_continuous(
    breaks = sort(unique(panel_b_dat$threshold))
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    x = "Jaccard similarity threshold",
    y = "Targets with at least one source analogue" ) +
  theme_minimal(base_size = 16)

panel_b

ggsave(
  file.path(fig_dir, 
    "fig06_panelB_threshold_coverage.png"
  ),
  plot = panel_b,
  width = 6,
  height = 6,
  dpi = 300
)

# ===========================================================
# 8) Build Panel C data
# ===========================================================
# Cumulative analogue availability among all data-poor ROCKS targets.
# Denominator is total targets, not just covered targets.

n_total_targets <- summary_now$n_total_targets[1]

panel_c_dat <- target_match_counts %>%
  summarise(
    `≥1` = sum(n_matches >= 1),
    `≥2` = sum(n_matches >= 2),
    `≥5` = sum(n_matches >= 5),
    `≥10` = sum(n_matches >= 10),
    `≥20` = sum(n_matches >= 20),
    `≥50` = sum(n_matches >= 50)
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "minimum_analogues",
    values_to = "n_targets"
  ) %>%
  mutate(
    pct_targets = 100 * n_targets / n_total_targets,
    minimum_analogues = factor(
      minimum_analogues,
      levels = c("≥1", "≥2", "≥5", "≥10", "≥20", "≥50")
    ),
    label = paste0(
      comma(n_targets))
    )


# ===========================================================
# 9) Panel C
# ===========================================================

panel_c <- ggplot(
  panel_c_dat,
  aes(
    x = minimum_analogues,
    y = pct_targets
  )
) +
  geom_col(
    fill = col_match_dist,
    width = 0.65
  ) +
  geom_text(
    aes(label = label),
    vjust = -0.25,
    size = 5,
    lineheight = 0.9
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    limits = c(0, max(panel_c_dat$pct_targets, na.rm = TRUE) * 1.18),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "Min. number source analogues at 0.7",
    y = "Targets") +
  theme_minimal(base_size = 16)

panel_c

ggsave(
  file.path(fig_dir, 
    "fig06_panelC_analogue_depth.png"
  ),
  plot = panel_c,
  width = 6,
  height = 6,
  dpi = 300
)

# Combine all panels
library(cowplot)
library(png)
library(grid)

img_top <- rasterGrob(
  readPNG(
    here(
      "outputs",
      "figures",
      "fig06_panelA_denominator.png"
    )
  ),
  interpolate = TRUE
)

img_left <- rasterGrob(
  readPNG(
    here(
      "outputs",
      "figures",
      "fig06_panelB_threshold_coverage.png"
    )
  ),
  interpolate = TRUE
)

img_right <- rasterGrob(
  readPNG(
    here(
      "outputs",
      "figures",
      "fig06_panelC_analogue_depth.png"
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
  file.path(fig_dir, 
  "fig06_combined.png"),
  final_fig,
  width = 10,
  height = 7,
  dpi = 600
)




