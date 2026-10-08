# ===========================================================
# fig02_depth_by_sector.R
# ===========================================================

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(scales)
library(readr)
library(here)
library(ggnewscale)

# INPUT

broad_profiling <- readr::read_csv(
  here::here(
    "data",
    "raw",
    "rocks_profile_pod_table_aeds_08sept26.csv"
  )
)

# ------------------------------------------------------------
# 0. Load derived objects
# ------------------------------------------------------------

rocks <- readRDS(here::here("data", "derived", "rocks.rds"))
tam <- readRDS(here::here("data", "derived", "tam.rds"))
tam_assess <- readRDS(here::here("data", "derived", "tam_assess.rds"))
b_assay_data_sum <- readRDS(here::here("data", "derived", "b_assay_data_sum.rds"))
b_mc7_assay <- readRDS(here::here("data", "derived", "b_mc7_assay.rds"))
b_exec_sum <- readRDS(
  here("data", "derived", "b_exec_sum.rds")
)

out_dir <- here::here("outputs", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source(
  here(
    "R",
    "ulri_theme.R"
  )
)

# ------------------------------------------------------------
# 1. Helpers
# ------------------------------------------------------------

pick_col <- function(df, candidates, required = TRUE) {
  hit <- candidates[candidates %in% names(df)]

  if (length(hit) == 0 && required) {
    stop(
      "None of these columns were found: ",
      paste(candidates, collapse = ", "),
      "\nAvailable columns are:\n",
      paste(names(df), collapse = ", ")
    )
  }

  if (length(hit) == 0) return(NA_character_)

  hit[[1]]
}

clean_sector_name <- function(x) {
  x %>%
    str_replace_all("_", " ") %>%
    str_replace_all("\\.", " ") %>%
    str_squish() %>%
    str_to_title()
}

# ------------------------------------------------------------
# 2. Create long sector table from ROCKS
# ------------------------------------------------------------
id_col_rocks <- "cichem_id"

sector_cols <- c(
  "consumer",
  "industrial",
  "agricultural",
  "food",
  "pharmaceutical",
  "biomonitoring",
  "env_monitoring"
)

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

rocks2 <- rocks %>%
  mutate(
    any_primary_sector =
      rowSums(across(all_of(primary_sector_cols),
                     ~ replace_na(.x, 0))) > 0
  )

sector_long <- rocks2 %>%
  select(
    chem_id = all_of(id_col_rocks),
    all_of(primary_sector_cols)
  ) %>%
  pivot_longer(
    cols = all_of(primary_sector_cols),
    names_to = "sector",
    values_to = "in_sector"
  ) %>%
  mutate(
    in_sector = in_sector == 1,
    sector = clean_sector_name(sector)
  ) %>%
  filter(
    in_sector,
    !is.na(chem_id),
    chem_id != ""
  ) %>%
  select(chem_id, sector)

fallback_sector_long <- rocks2 %>%
  filter(!any_primary_sector) %>%
  select(
    chem_id = all_of(id_col_rocks),
    all_of(fallback_sector_cols)
  ) %>%
  pivot_longer(
    cols = all_of(fallback_sector_cols),
    names_to = "sector",
    values_to = "in_sector"
  ) %>%
mutate(
  in_sector = in_sector == 1,
  sector = clean_sector_name(sector),
  sector = dplyr::recode(
    sector,
    "Biomonitoring" = "Unknown\nBiomonitoring",
    "Env Monitoring" = "Unknown\nEnv. Monitoring"
  )
) %>%
  filter(
    in_sector,
    !is.na(chem_id),
    chem_id != ""
  ) %>%
  select(chem_id, sector)

unknown_sector_long <- rocks2 %>%
  filter(
    !any_primary_sector,
    rowSums(across(all_of(fallback_sector_cols),
                   ~ replace_na(.x, 0))) == 0
  ) %>%
  transmute(
    chem_id = .data[[id_col_rocks]],
    sector = "Unknown"
  )

sector_long <- bind_rows(
  sector_long,
  fallback_sector_long,
  unknown_sector_long
) %>%
  distinct()

# ------------------------------------------------------------
# 3. In vivo DRSV depth per chemical
# ------------------------------------------------------------

id_col_tam <- "cichem_id"

study_col_tam <- "source_hash"

if (is.na(study_col_tam)) {
  message("No study-level identifier found. Using toxval_id/row-level DRSV record counts instead.")
  study_col_tam <- pick_col(tam, c("toxval_id"), required = TRUE)
}

toxval_id_col <- pick_col(tam, c("toxval_id"), required = FALSE)

drsv_depth <- tam %>%
  filter(toxval_type_supercategory == "Dose Response Summary Value") %>%
  filter(!is.na(.data[[id_col_tam]]), .data[[id_col_tam]] != "") %>%
  group_by(chem_id = .data[[id_col_tam]]) %>%
  summarise(
    n_in_vivo_drsv_studies = n_distinct(.data[[study_col_tam]], na.rm = TRUE),
    n_in_vivo_drsv_records = if (!is.na(toxval_id_col)) {
      n_distinct(.data[[toxval_id_col]], na.rm = TRUE)
    } else {
      n()
    },
    .groups = "drop"
  )

# ------------------------------------------------------------
# 4. In vitro targeted assay endpoint depth per chemical
# ------------------------------------------------------------

targeted_depth <- b_exec_sum %>%
  transmute(
    chem_id = .data[["cichem_id"]],
    n_in_vitro_targeted_assay_endpoints =
      as.numeric(n_assay_endpoints_tested)
  ) %>%
  distinct()

# ------------------------------------------------------------
# 4b. Broad-profiling POD depth per chemical
# ------------------------------------------------------------

broad_depth <- broad_profiling %>%
  transmute(
    chem_id = cichem_id,
    n_broad_profiling_unique_pods = as.numeric(n_pods_uniq)
  ) %>%
  filter(
    !is.na(chem_id),
    chem_id != ""
  ) %>%
  group_by(chem_id) %>%
  summarise(
    # Protect against accidental duplicate chemical rows
    n_broad_profiling_unique_pods = max(
      n_broad_profiling_unique_pods,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    # max(..., na.rm = TRUE) returns -Inf if all values were NA
    n_broad_profiling_unique_pods = if_else(
      is.infinite(n_broad_profiling_unique_pods),
      NA_real_,
      n_broad_profiling_unique_pods
    )
  )

# ------------------------------------------------------------
# 5. Combine chemical-level depth with sectors
# ------------------------------------------------------------

depth_by_sector <- sector_long %>%
  left_join(
    drsv_depth,
    by = "chem_id"
  ) %>%
  left_join(
    targeted_depth,
    by = "chem_id"
  ) %>%
  left_join(
    broad_depth,
    by = "chem_id"
  ) %>%
  mutate(
    n_in_vivo_drsv_studies = replace_na(
      n_in_vivo_drsv_studies,
      0L
    ),
    n_in_vivo_drsv_records = replace_na(
      n_in_vivo_drsv_records,
      0L
    ),
    n_in_vitro_targeted_assay_endpoints = replace_na(
      n_in_vitro_targeted_assay_endpoints,
      0
    ),
    n_broad_profiling_unique_pods = replace_na(
      n_broad_profiling_unique_pods,
      0
    )
  )

# ------------------------------------------------------------
# 6. Build threshold data
# ------------------------------------------------------------

threshold_df <- depth_by_sector %>%
  select(
    chem_id,
    sector,
    n_in_vivo_drsv_studies,
    n_in_vitro_targeted_assay_endpoints,
    n_broad_profiling_unique_pods
  ) %>%
  distinct() %>%
  pivot_longer(
    cols = c(
      n_in_vivo_drsv_studies,
      n_in_vitro_targeted_assay_endpoints,
      n_broad_profiling_unique_pods
    ),
    names_to = "stream",
    values_to = "count"
  ) %>%
  mutate(
    stream = recode(
      stream,
      n_in_vivo_drsv_studies = "In vivo",
      n_in_vitro_targeted_assay_endpoints = "Targeted",
      n_broad_profiling_unique_pods = "Broad profiling"
    ),
    count = replace_na(
      as.numeric(count),
      0
    )
  )

threshold_lookup <- tibble::tribble(
  ~stream,             ~threshold_label, ~threshold_value, ~threshold_order,
  "In vivo",           "≥1",             1,                1,
  "In vivo",           "≥2",             2,                2,
  "In vivo",           "≥5",             5,                3,
  "Targeted",          "≥1",             1,                1,
  "Targeted",          "≥100",           100,              2,
  "Targeted",          "≥300",           300,              3,
  "Broad profiling",   "≥1",             1,                1,
  "Broad profiling",   "≥2",             2,                2,
  "Broad profiling",   "≥4",             4,                3
)

# ------------------------------------------------------------
# 7. Calculate percent of each sector inventory meeting thresholds
# ------------------------------------------------------------

sector_totals <- threshold_df %>%
  distinct(sector, chem_id) %>%
  count(sector, name = "n_sector_total")

threshold_counts_raw <- threshold_df %>%
  left_join(threshold_lookup, by = "stream") %>%
  filter(count >= threshold_value) %>%
  group_by(
    sector,
    stream,
    threshold_label,
    threshold_value,
    threshold_order
  ) %>%
  summarise(
    n_chems = n_distinct(chem_id),
    .groups = "drop"
  )

threshold_counts <- threshold_df %>%
  distinct(sector, stream) %>%
  left_join(threshold_lookup, by = "stream") %>%
  left_join(
    threshold_counts_raw,
    by = c(
      "sector",
      "stream",
      "threshold_label",
      "threshold_value",
      "threshold_order"
    )
  ) %>%
  left_join(sector_totals, by = "sector") %>%
  mutate(
    n_chems = replace_na(n_chems, 0L),
    pct_sector = 100 * n_chems / n_sector_total
  )

# ------------------------------------------------------------
# 8. Order sectors
# ------------------------------------------------------------

sector_order <- c(
  "Industrial",
  "Consumer",
  "Agricultural",
  "Food",
  "Pharmaceutical",
  "Unknown\nBiomonitoring",
  "Unknown\nEnv. Monitoring"
)

stream_order <- c(
  "In vivo",
  "Targeted bioactivity",
  "Broad profiling"
)

stream_order <- c(
  "In vivo",
  "Targeted bioactivity",
  "Broad profiling"
)

threshold_counts <- threshold_counts %>%
  mutate(
    stream = recode(
      stream,
      "Targeted" = "Targeted bioactivity"
    ),
    sector = factor(sector, levels = sector_order),
    stream = factor(stream, levels = stream_order)
  )

# ------------------------------------------------------------
# 9. Split data by stream for separate legends
# ------------------------------------------------------------

threshold_targeted <- threshold_counts %>%
  filter(stream == "Targeted bioactivity") %>%
  mutate(
    threshold_label = factor(
      threshold_label,
      levels = c("≥1", "≥100", "≥300")
    )
  )

threshold_invivo <- threshold_counts %>%
  filter(stream == "In vivo") %>%
  mutate(
    threshold_label = factor(
      threshold_label,
      levels = c("≥1", "≥2", "≥5")
    )
  )
threshold_broad <- threshold_counts %>%
  filter(stream == "Broad profiling") %>%
  mutate(
    threshold_label = factor(
      threshold_label,
      levels = c("≥1", "≥2", "≥4")
    )
  )

# Add assessment annotation----


assessment_chems <- tam_assess %>%
  filter(
    population %in% c(
      "General Population"
    )
  ) %>%
  distinct(
    chem_id = .data[["cichem_id"]]
  )

assessment_by_sector <- sector_long %>%
  mutate(
    has_assessment = chem_id %in% assessment_chems$chem_id
  ) %>%
  group_by(sector) %>%
  summarise(
    n_sector_total = n_distinct(chem_id),
    n_assessment = n_distinct(chem_id[has_assessment]),
    pct_assessment = 100 * n_assessment / n_sector_total,
    .groups = "drop"
  ) %>%
  mutate(
    sector = factor(
      sector,
      levels = levels(threshold_counts$sector)
    )
  )

assessment_by_sector <- assessment_by_sector %>%
  mutate(
    stream = factor(
      "In vivo",
      levels = stream_order
    )
  )

# ------------------------------------------------------------
# 10. Final figure
# ------------------------------------------------------------

p_threshold_pct <- ggplot() +

  # Targeted assay endpoints
  geom_col(
    data = threshold_targeted,
    aes(
      x = sector,
      y = pct_sector,
      fill = threshold_label
    ),
    position = position_dodge(width = 0.82),
    width = 0.72,
    color = "white",
    linewidth = 0.2
  ) +
  scale_fill_manual(
    name = "Assay endpoints",
    values = c(
      "≥1" = "#f3968fd5",
      "≥100" = "#FF6256",
      "≥300" = "#c05c54"
    ),
    breaks = c("≥1", "≥100", "≥300"),
    guide = guide_legend(order =3,
      nrow = 1,
      title.position = "top"
    )
  ) +

  ggnewscale::new_scale_fill() +

  # In vivo DRSVs
  geom_col(
    data = threshold_invivo,
    aes(
      x = sector,
      y = pct_sector,
      fill = threshold_label
    ),
    position = position_dodge(width = 0.82),
    width = 0.72,
    color = "white",
    linewidth = 0.2
  ) +
  scale_fill_manual(
    name = "DRSVs",
    values = c(
      "≥1" = "#1fd0d67c",
      "≥2" = "#1FD0D6",
      "≥5" = "#38a0a3"
    ),
    breaks = c("≥1", "≥2", "≥5"),
    guide = guide_legend(order =2,
      nrow = 1,
      title.position = "top"
    )
  ) +

  ggnewscale::new_scale_fill() +

  # Broad-profiling unique PODs
  geom_col(
    data = threshold_broad,
    aes(
      x = sector,
      y = pct_sector,
      fill = threshold_label
    ),
    position = position_dodge(width = 0.82),
    width = 0.72,
    color = "white",
    linewidth = 0.2
  ) +
  scale_fill_manual(
    name = "Unique PODs",
    values = c(
      "≥1" = "#2792848f",
      "≥2" = "#279284",
      "≥4" = "#2b6e65"
    ),
    breaks = c("≥1", "≥2", "≥4"),
    guide = guide_legend(order=4,
      nrow = 1,
      title.position = "top"
    )
  ) +

  # Assessment-value annotation
scale_linetype_manual(
  name = NULL,
  values = c("Assessment values" = 1),
  guide = guide_legend(order = 1)
) +
scale_shape_manual(
  name = NULL,
  values = c("Assessment values" = 21),
  guide = guide_legend(order = 1)
) +
  geom_line(
    data = assessment_by_sector,
    aes(
      x = sector,
      y = pct_assessment,
      group = 1,
      linetype = "Assessment values"
    ),
    inherit.aes = FALSE,
    linewidth = 0.8,
    color = "black"
  ) +
  geom_point(
    data = assessment_by_sector,
    aes(
      x = sector,
      y = pct_assessment,
      shape = "Assessment values"
    ),
    inherit.aes = FALSE,
    size = 2.8,
    stroke = 1,
    color = "black",
    fill = "white"
  ) +

facet_wrap(
  ~ stream,
  nrow = 1,
  drop = FALSE
) +
  scale_y_continuous(
    labels = label_percent(scale = 1),
    expand = expansion(mult = c(0, 0.06))
  ) +
  labs(
    x = NULL,
    y = "Chemicals within sector (%)"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    strip.text = element_text(face = "bold",size=18),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      angle = 35,
      hjust = 1,
      vjust = 1
    ),
    legend.position = "bottom",
    legend.box = "horizontal"
  )

p_threshold_pct

ggsave(
  file.path(out_dir, "rocks_depth_threshold_pct_by_sector.png"),
  p_threshold_pct,
  width = 14,
  height = 7,
  dpi = 300
)

message(
  "fig02_depth_by_sector.R completed successfully"
)

message(
  "Figure saved:"
)

message(
  file.path(
    out_dir,
    "rocks_depth_threshold_pct_by_sector.png"
  )
)

