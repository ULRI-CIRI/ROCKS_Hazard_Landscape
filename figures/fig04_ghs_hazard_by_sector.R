# =============================================================================
# figS3_ghs_hazard_by_sector.R
# =============================================================================

library(tidyverse)
library(here)
library(stringr)
library(scales)

# -----------------------------------------------------------------------------
# 1. File paths
# -----------------------------------------------------------------------------

ghs_file <- here("data", "derived", "ghs_results.csv")
cid_file <- here("data", "raw", "cids.txt")
rocks_file <- here("data", "derived", "rocks.rds")

fig_dir <- here("outputs", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

out_heatmap_png <- here(
  "outputs",
  "figures",
  "ghs_selected_hazards_by_sector_heatmap.png"
)
# -----------------------------------------------------------------------------
# 2. Selected GHS categories
# -----------------------------------------------------------------------------
# Keep the same order as Figure 4 Panel A/B.

selected_ghs_categories <- c(
  "Reproductive toxicity (D)",
  "Reproductive toxicity (W)",
  "Germ cell mutagenicity (D)",
  "Germ cell mutagenicity (W)",
  "Carcinogenicity (D)",
  "Carcinogenicity (W)",
  "STOT-RE (D)",
  "STOT-RE (W)"
)

# -----------------------------------------------------------------------------
# 3. Read input data
# -----------------------------------------------------------------------------

ghs_raw <- readr::read_csv(
  ghs_file,
  show_col_types = FALSE
)

cid_xwalk <- readr::read_tsv(
  cid_file,
  col_names = c("inchikey", "cid"),
  show_col_types = FALSE
) %>%
  mutate(
    cid = as.character(cid),
    inchikey = as.character(inchikey)
  ) %>%
  distinct(inchikey, cid)

rocks <- readRDS(rocks_file) %>%
  mutate(
    across(
      any_of(c("cichem_id", "inchikey", "sector")),
      as.character
    )
  )

# -----------------------------------------------------------------------------
# 4. Build sector annotation table
# -----------------------------------------------------------------------------
# Primary sectors may overlap.
#
# Biomonitoring and environmental monitoring are treated as fallback categories:
# they are assigned only when a chemical has no primary-sector assignment.
#
# Chemicals with no primary or fallback assignment are labeled "Unknown".
# -----------------------------------------------------------------------------

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

sector_cols <- c(
  primary_sector_cols,
  fallback_sector_cols
)

# Check that all required sector columns are present
missing_sector_cols <- setdiff(
  sector_cols,
  names(rocks)
)

if (length(missing_sector_cols) > 0) {
  stop(
    "The following required sector columns are missing from rocks: ",
    paste(missing_sector_cols, collapse = ", ")
  )
}

sector_annot <- rocks %>%
  distinct(
    cichem_id,
    .keep_all = TRUE
  ) %>%
  select(
    cichem_id,
    all_of(sector_cols)
  ) %>%
  mutate(
    # Convert TRUE/FALSE, 1/0, and missing values to logical indicators
    across(
      all_of(sector_cols),
      ~ coalesce(
        .x == TRUE | .x == 1,
        FALSE
      )
    ),

    # A chemical has a primary sector if at least one primary flag is TRUE
    has_primary_sector = if_any(
      all_of(primary_sector_cols),
      ~ .x
    ),

    # Fallback categories are used only when there is no primary sector
    biomonitoring_unknown =
      !has_primary_sector &
      biomonitoring,

    env_monitoring_unknown =
      !has_primary_sector &
      env_monitoring,

    # Fully unknown means no primary sector and no fallback sector
    unknown =
      !has_primary_sector &
      !biomonitoring &
      !env_monitoring
  ) %>%
  pivot_longer(
    cols = c(
      all_of(primary_sector_cols),
      biomonitoring_unknown,
      env_monitoring_unknown,
      unknown
    ),
    names_to = "sector",
    values_to = "in_sector"
  ) %>%
  filter(
    in_sector
  ) %>%
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
  ) %>%
  distinct(
    cichem_id,
    sector
  )

# -----------------------------------------------------------------------------
# 5. Attach InChIKey and bridge PubChem results back to cichem_id
# -----------------------------------------------------------------------------

ghs_raw_joined <- ghs_raw %>%
  mutate(
    cid = as.character(cid)
  ) %>%
  left_join(
    cid_xwalk,
    by = "cid",
    suffix = c("", "_xwalk")
  )

# Handle possible duplicate InChIKey columns after the join.
if ("inchikey_xwalk" %in% names(ghs_raw_joined)) {
  if ("inchikey" %in% names(ghs_raw_joined)) {
    ghs_raw_joined <- ghs_raw_joined %>%
      mutate(
        inchikey = coalesce(
          as.character(inchikey),
          as.character(inchikey_xwalk)
        )
      ) %>%
      select(-inchikey_xwalk)
  } else {
    ghs_raw_joined <- ghs_raw_joined %>%
      rename(
        inchikey = inchikey_xwalk
      )
  }
}

if (!"inchikey" %in% names(ghs_raw_joined)) {
  stop("No inchikey column was available after joining ghs_results.csv to cids.txt.")
}

ghs_base <- ghs_raw_joined %>%
  mutate(
    cid = as.character(cid),
    inchikey = as.character(inchikey)
  ) %>%
  left_join(
    rocks %>%
      select(cichem_id, inchikey) %>%
      distinct(),
    by = "inchikey"
  ) %>%
  filter(
    !is.na(cichem_id)
  )

message(
  "GHS records mapped to ROCKS cichem_id: ",
  scales::comma(n_distinct(ghs_base$cichem_id))
)

# -----------------------------------------------------------------------------
# 6. Map GHS H-codes to selected health hazard categories
# -----------------------------------------------------------------------------

map_ghs_health_category <- function(h_code) {

  case_when(

    # Acute Toxicity
    h_code %in% c(
      "H300", "H301", "H302",
      "H310", "H311", "H312",
      "H330", "H331", "H332", "H333"
    ) ~ "Acute toxicity",

    # Skin Corrosion / Irritation
    h_code %in% c(
      "H314", "H315", "H316"
    ) ~ "Skin corrosion / irritation",

    # Serious Eye Damage / Eye Irritation
    h_code %in% c(
      "H318", "H319", "H320"
    ) ~ "Serious eye damage / eye irritation",

    # Respiratory Sensitization
    h_code == "H334" ~ "Respiratory sensitization",

    # Skin Sensitization
    h_code == "H317" ~ "Skin sensitization",

    # Germ Cell Mutagenicity
    h_code %in% c(
      "H340", "H341"
    ) ~ "Germ cell mutagenicity",

    # Carcinogenicity
    str_detect(h_code, "^H350[A-Za-z]*$") |
      h_code == "H351" ~ "Carcinogenicity",

    # Reproductive Toxicity
    str_detect(h_code, "^H360[A-Za-z]*$") |
      str_detect(h_code, "^H361[A-Za-z]*$") |
      h_code == "H362" ~ "Reproductive toxicity",

    # STOT Single Exposure
    h_code %in% c(
      "H335", "H336",
      "H370", "H371"
    ) ~ "STOT single exposure",

    # STOT Repeated Exposure
    h_code %in% c(
      "H372", "H373"
    ) ~ "STOT-RE",

    # Aspiration Hazard
    h_code == "H304" ~ "Aspiration hazard",

    TRUE ~ "Other / non-health or unmapped"
  )
}

map_ghs_signal_word <- function(h_code) {

  case_when(

    # Acute toxicity
    h_code %in% c("H300", "H301", "H310", "H311", "H330", "H331") ~ "Danger",
    h_code %in% c("H302", "H312", "H332") ~ "Warning",

    # Aspiration
    h_code == "H304" ~ "Danger",

    # Skin / eye
    h_code %in% c("H314", "H318") ~ "Danger",
    h_code %in% c("H315", "H316", "H319", "H320") ~ "Warning",

    # Sensitization
    h_code == "H334" ~ "Danger",
    h_code == "H317" ~ "Warning",

    # STOT single exposure
    h_code == "H370" ~ "Danger",
    h_code %in% c("H335", "H336", "H371") ~ "Warning",

    # Mutagenicity
    h_code == "H340" ~ "Danger",
    h_code == "H341" ~ "Warning",

    # Carcinogenicity
    str_detect(h_code, "^H350[A-Za-z]*$") ~ "Danger",
    h_code == "H351" ~ "Warning",

    # Reproductive toxicity
    str_detect(h_code, "^H360[A-Za-z]*$") ~ "Danger",
    str_detect(h_code, "^H361[A-Za-z]*$") ~ "Warning",
    h_code == "H362" ~ "Warning",

    # STOT repeated exposure
    h_code == "H372" ~ "Danger",
    h_code == "H373" ~ "Warning",

    TRUE ~ NA_character_
  )
}

ghs_long <- ghs_base %>%
  mutate(
    cid = as.character(cid),
    hazard_statements = as.character(hazard_statements)
  ) %>%
  filter(
    !is.na(hazard_statements),
    hazard_statements != ""
  ) %>%
  mutate(
    h_code_list = str_extract_all(
      hazard_statements,
      "\\bH[0-9]{3}[A-Za-z]*\\b"
    )
  ) %>%
  select(
    cid,
    inchikey,
    cichem_id,
    any_of(c("dtxsid")),
    any_of(c("compound_name")),
    any_of(c("signal_word")),
    hazard_statements,
    h_code_list
  ) %>%
  unnest(h_code_list) %>%
  rename(
    h_code = h_code_list
  ) %>%
  distinct(
    cichem_id,
    cid,
    h_code,
    .keep_all = TRUE
  )

ghs_health <- ghs_long %>%
  mutate(
    ghs_health_category = map_ghs_health_category(h_code),
    ghs_signal_word = map_ghs_signal_word(h_code)
  ) %>%
  filter(
    ghs_health_category != "Other / non-health or unmapped",
    !is.na(ghs_signal_word)
  ) %>%
  mutate(
    hazard_severity = paste0(
      ghs_health_category,
      " (",
      substr(ghs_signal_word, 1, 1),
      ")"
    )
  )

# -----------------------------------------------------------------------------
# 7. Join sector annotations to GHS health data
# -----------------------------------------------------------------------------

ghs_sector <- ghs_health %>%
  left_join(
    sector_annot,
    by = "cichem_id"
  ) %>%
  filter(
    !is.na(sector),
    sector != ""
  )

message(
  "GHS health-classified chemicals with sector annotations: ",
  scales::comma(n_distinct(ghs_sector$cichem_id))
)

# -----------------------------------------------------------------------------
# 8. Sector denominators
# -----------------------------------------------------------------------------
# Denominator:
# ROCKS chemicals represented in the GHS retrieval results and assigned to
# each sector.
#
# Primary sectors may overlap. Each chemical is counted only once within
# a given sector.
#
# Unknown Biomonitoring, Unknown Env. monitoring, and Unknown are mutually
# exclusive fallback categories.
# -----------------------------------------------------------------------------

sector_denominators <- ghs_base %>%
  distinct(
    cichem_id
  ) %>%
  inner_join(
    sector_annot,
    by = "cichem_id"
  ) %>%
  filter(
    !is.na(sector),
    sector != ""
  ) %>%
  distinct(
    cichem_id,
    sector
  ) %>%
  count(
    sector,
    name = "n_chemicals_in_sector"
  ) %>%
  arrange(
    desc(n_chemicals_in_sector),
    sector
  )

print(
  sector_denominators,
  n = Inf
)

# -----------------------------------------------------------------------------
# 9. Sector-by-hazard summary
# -----------------------------------------------------------------------------

sector_hazard_summary <- ghs_sector %>%
  filter(
    hazard_severity %in% selected_ghs_categories
  ) %>%
  distinct(
    cichem_id,
    sector,
    hazard_severity
  ) %>%
  count(
    sector,
    hazard_severity,
    name = "n_chemicals"
  ) %>%
  left_join(
    sector_denominators,
    by = "sector"
  ) %>%
  mutate(
    pct_of_sector = if_else(
      n_chemicals_in_sector > 0,
      n_chemicals / n_chemicals_in_sector,
      NA_real_
    )
  )

# -----------------------------------------------------------------------------
# 10. Complete sector x hazard matrix
# -----------------------------------------------------------------------------

all_sectors <- sector_denominators %>%
  arrange(
    desc(n_chemicals_in_sector),
    sector
  ) %>%
  pull(sector)

heatmap_sector <- tidyr::expand_grid(
  sector = all_sectors,
  hazard_severity = selected_ghs_categories
) %>%
  left_join(
    sector_hazard_summary,
    by = c(
      "sector",
      "hazard_severity"
    )
  ) %>%
  left_join(
    sector_denominators,
    by = "sector",
    suffix = c("", "_denom")
  ) %>%
  mutate(
    n_chemicals = replace_na(
      n_chemicals,
      0L
    ),

    n_chemicals_in_sector = coalesce(
      n_chemicals_in_sector,
      n_chemicals_in_sector_denom
    ),

    pct_of_sector = if_else(
      n_chemicals_in_sector > 0,
      n_chemicals / n_chemicals_in_sector,
      NA_real_
    ),

    hazard_severity = factor(
      hazard_severity,
      levels = selected_ghs_categories
    )
  ) %>%
  select(
    sector,
    hazard_severity,
    n_chemicals,
    n_chemicals_in_sector,
    pct_of_sector
  )

heatmap_sector_plot <- heatmap_sector

# -----------------------------------------------------------------------------
# 11. Order sectors by denominator size and create plot labels
# -----------------------------------------------------------------------------

sector_order <- sector_denominators %>%
  semi_join(
    heatmap_sector_plot %>%
      distinct(sector),
    by = "sector"
  ) %>%
  arrange(
    desc(n_chemicals_in_sector),
    sector
  ) %>%
  pull(sector)

sector_labels <- sector_denominators %>%
  filter(
    sector %in% sector_order
  ) %>%
  mutate(
    sector = factor(
      sector,
      levels = sector_order
    )
  ) %>%
  arrange(
    sector
  ) %>%
  mutate(
    sector_label = paste0(
      as.character(sector),
      "\n(n=",
      scales::comma(n_chemicals_in_sector),
      ")"
    )
  )

sector_label_order <- sector_labels %>%
  pull(sector_label)

heatmap_sector_plot <- heatmap_sector_plot %>%
  mutate(
    sector = as.character(sector)
  ) %>%
  left_join(
    sector_labels %>%
      mutate(
        sector = as.character(sector)
      ) %>%
      select(
        sector,
        sector_label
      ),
    by = "sector"
  ) %>%
  mutate(
    sector = factor(
      sector,
      levels = sector_order
    ),

    sector_label = factor(
      sector_label,
      levels = sector_label_order
    ),

    hazard_severity = factor(
      as.character(hazard_severity),
      levels = selected_ghs_categories
    )
  )

# QC: confirm final displayed order
tibble(
  displayed_order = seq_along(sector_order),
  sector = sector_order
) %>%
  left_join(
    sector_denominators,
    by = "sector"
  ) %>%
  print(n = Inf)
# -----------------------------------------------------------------------------
# 12. Plot sector heatmap
# -----------------------------------------------------------------------------
max_pct <- max(
  heatmap_sector_plot$pct_of_sector,
  na.rm = TRUE
)

p_sector_heatmap <- ggplot(
  heatmap_sector_plot,
  aes(
    x = hazard_severity,
    y = sector_label,
    fill = pct_of_sector
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.3
  ) +
  scale_x_discrete(
    drop = FALSE
  ) +
  scale_y_discrete(
    drop = FALSE
  ) +

scale_fill_gradientn(
  colours = c(
    "#F7FBFF",
    "#0a91ff52",
    "#0eaffab4",
    "#0A32FF"
  ),
  name = "Percent of sector",
  breaks = pretty(c(0, max_pct), n = 3),
  labels = scales::label_percent(accuracy = 1)
)+
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_minimal(
    base_size = 14
  ) +
  theme(
    panel.grid = element_blank(),

    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),

    axis.text.y = element_text(
      size = 12
    ),

    plot.title = element_text(
      face = "bold"
    ),

    legend.position = "bottom"
  ) +
  coord_flip()

print(p_sector_heatmap)

ggsave(
  filename = out_heatmap_png,
  plot = p_sector_heatmap,
  width = 7.5,
  height = 6,
  dpi = 300
)
