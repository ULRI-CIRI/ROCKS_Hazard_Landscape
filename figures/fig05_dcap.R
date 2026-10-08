out_dir <- here::here(
  "outputs",
  "figures"
)

tam_assess <- readRDS(
  here::here("data", "derived", "tam_assess.rds")
)

# =============================================================================
# Assessment value coverage by sector
# Traditional assessment values + DCAP CTVs
# =============================================================================

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

# -------------------------------------------------------------------------
# Traditional assessment values
# -------------------------------------------------------------------------

assess_flag <- tam_assess %>%
  filter(
    population == "General Population",
    toxval_type != "CTV"
  ) %>%
  distinct(cichem_id) %>%
  mutate(
    has_assessment = TRUE
  )

# -------------------------------------------------------------------------
# DCAP values
# -------------------------------------------------------------------------

dcap_flag <- tam_assess %>%
  filter(
    population == "General Population",
    toxval_type == "CTV"
  ) %>%
  distinct(cichem_id) %>%
  mutate(
    has_dcap = TRUE
  )

# -------------------------------------------------------------------------
# Inventory with coverage flags
# -------------------------------------------------------------------------

rocks_assessment_panel <- rocks %>%
  distinct(cichem_id, .keep_all = TRUE) %>%
  left_join(
    assess_flag,
    by = "cichem_id"
  ) %>%
  left_join(
    dcap_flag,
    by = "cichem_id"
  ) %>%
  mutate(
    has_assessment =
      replace_na(has_assessment, FALSE),

    has_dcap =
      replace_na(has_dcap, FALSE),

    coverage_group = case_when(
      has_assessment ~ "Traditional assessment value",

      has_dcap & !has_assessment ~
        "DCAP value only",

      TRUE ~
        "No assessment value"
    )
  )

# -------------------------------------------------------------------------
# Sector assignments
# -------------------------------------------------------------------------

sector_long <- rocks_assessment_panel %>%
  select(
    cichem_id,
    coverage_group,
    all_of(sector_cols)
  ) %>%
  mutate(
    across(
      all_of(sector_cols),
      ~ coalesce(.x == TRUE | .x == 1, FALSE)
    ),

    has_primary_sector =
      if_any(
        all_of(primary_sector_cols),
        ~ .x
      ),

    biomonitoring_unknown =
      !has_primary_sector & biomonitoring,

    env_monitoring_unknown =
      !has_primary_sector & env_monitoring,

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
    values_to = "in_sector_raw"
  ) %>%

  filter(in_sector_raw) %>%

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
# -------------------------------------------------------------------------
# Sector totals
# -------------------------------------------------------------------------

sector_totals <- sector_long %>%
  distinct(
    cichem_id,
    sector
  ) %>%
  count(
    sector,
    name = "total_n"
  )

# -------------------------------------------------------------------------
# Coverage counts
# -------------------------------------------------------------------------

sector_counts <- sector_long %>%
  filter(
    coverage_group != "No assessment value"
  ) %>%
  distinct(
    cichem_id,
    sector,
    coverage_group
  ) %>%
  count(
    sector,
    coverage_group,
    name = "n"
  ) %>%
  left_join(
    sector_totals,
    by = "sector"
  ) %>%
  mutate(
    pct = n / total_n
  )

# -------------------------------------------------------------------------
# Order sectors by total assessment coverage
# -------------------------------------------------------------------------

sector_order <- sector_counts %>%
  group_by(sector) %>%
  summarise(
    pct_total = sum(pct),
    .groups = "drop"
  ) %>%
  arrange(desc(pct_total)) %>%
  pull(sector)

# -------------------------------------------------------------------------
# Create ordered labels
# -------------------------------------------------------------------------

sector_labels <- sector_totals %>%
  mutate(
    sector = factor(
      sector,
      levels = sector_order
    )
  ) %>%
  arrange(sector) %>%
  mutate(
    sector_label = paste0(
      sector,
      "\n(n=",
      scales::comma(total_n),
      ")"
    )
  )

sector_counts <- sector_counts %>%
  left_join(
    sector_labels %>%
      select(
        sector,
        sector_label
      ),
    by = "sector"
  ) %>%
  mutate(
    sector_label = factor(
      sector_label,
      levels = sector_labels$sector_label
    ),
coverage_group = factor(
  coverage_group,
  levels = c(
    "DCAP value only",
    "Traditional assessment value"
  )
    )
  )

dcap_labels <- sector_counts %>%
  select(
    sector_label,
    coverage_group,
    pct
  ) %>%
  pivot_wider(
    names_from = coverage_group,
    values_from = pct,
    values_fill = 0
  ) %>%
  mutate(
    pct_blue =
      `DCAP value only` /
      (`DCAP value only` + `Traditional assessment value`),

    label = scales::percent(
      pct_blue,
      accuracy = 1
    ),

    y =
      `Traditional assessment value` +
      (`DCAP value only` / 2)
  ) %>%
  filter(`DCAP value only` > 0.003)

# -------------------------------------------------------------------------
# Labels for totals at top of bars
# -------------------------------------------------------------------------

label_dat <- sector_counts %>%
  group_by(sector_label) %>%
  summarise(
    pct_total = sum(pct),
    .groups = "drop"
  )

# -------------------------------------------------------------------------
# Plot
# -------------------------------------------------------------------------

p_assessment_sector <- ggplot(
  sector_counts,
  aes(
    x = sector_label,
    y = pct,
    fill = coverage_group
  )
) +
  geom_col() +
geom_text(
  data = dcap_labels,
  aes(
    x = sector_label,
    y = y,
    label = label
  ),
  inherit.aes = FALSE,
  color = "white",
  size = 4,
  fontface = "bold"
)+
  scale_fill_manual(
    values = c(
      "Traditional assessment value" = "#BCE4F7",
      "DCAP value only" = "#0A32FF"
    )
  ) +
  scale_y_continuous(
    labels = scales::percent_format()
  ) +
  labs(
    x = NULL,
    y = "Chemicals with oral assessment value (%)",
    fill = NULL,
  ) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "bottom",
  axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    size=14
  )
)

p_assessment_sector

ggsave(
  file.path(out_dir, 
    "fig05_dcap.png"
  ),
  plot = p_assessment_sector,
  width = 7,
  height = 7,
  dpi = 300
)

n_inventory <- n_distinct(rocks$cichem_id)

n_assessment <- tam_assess %>%
  filter(
    population == "General Population",
    toxval_type != "CTV"
  ) %>%
  distinct(cichem_id) %>%
  nrow()

n_dcap <- tam_assess %>%
  filter(
    population == "General Population",
    toxval_type == "CTV"
  ) %>%
  distinct(cichem_id) %>%
  nrow()

n_combined <- tam_assess %>%
  filter(
    population == "General Population"
  ) %>%
  distinct(cichem_id) %>%
  nrow()

dcap_only <- tam_assess %>%
  filter(
    population == "General Population",
    toxval_type == "CTV"
  ) %>%
  distinct(cichem_id) %>%
  anti_join(
    tam_assess %>%
      filter(
        population == "General Population",
        toxval_type != "CTV"
      ) %>%
      distinct(cichem_id),
    by = "cichem_id"
  ) %>%
  nrow()
dcap_only / n_inventory
n_combined / n_assessment
remaining_gap <- 1 - (n_combined / n_inventory)
remaining_gap
