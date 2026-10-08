library(readxl)
library(dplyr)
library(here)
library(janitor)

source(here("R", "ulri_theme.R"))

# =============================================================================
# Read EU REACH registered substances
# =============================================================================

eu_reg_subs <- read_xlsx(
  here("data", "raw", "reach-reg-substances_eu_2023_FILTERED.xlsx"),
  sheet = "results"
) %>%
  clean_names() 

# Quick checks
head(eu_reg_subs)

n_missing_cas <- eu_reg_subs %>%
  filter(cas_number == "-") %>%
  nrow()

n_unique_cas <- n_distinct(eu_reg_subs$cas_number)

cat("Rows with missing CAS ('-'):", n_missing_cas, "\n")
cat("Unique CAS RNs:", n_unique_cas, "\n")

# =============================================================================
# Create unique CAS RN list
# =============================================================================

unique_casrn <- eu_reg_subs %>%
  filter(cas_number != "-") %>%
  distinct(cas_number)

# Optional: build SQL IN clause
escaped_vals <- gsub("'", "''", unique_casrn$cas_number)

values_clause <- paste0(
  "'", escaped_vals, "'",
  collapse = ", "
)

# =============================================================================
# Assign REACH Annex based on tonnage band
# =============================================================================

annex_vii <- c(
  "Intermediate use only",
  "≥ 1 to < 10 tonnes",
  "≥ 0 to < 10 tonnes",
  "≥ 0 tonnes",
  "≥ 1 tonnes"
)

annex_viii <- c(
  "≥ 10 to < 100 tonnes",
  "≥ 10 tonnes"
)

annex_ix <- c(
  "≥ 100 to < 1000 tonnes",
  "≥ 100 tonnes"
)

annex_x <- c(
  "≥ 1000 to < 10000 tonnes",
  "≥ 10000 to < 100000 tonnes",
  "≥ 1000 tonnes",
  "≥ 100000 to < 1000000 tonnes",
  "≥ 10000 tonnes",
  "≥ 1000000 to < 10000000 tonnes",
  "≥ 100000 tonnes",
  "≥ 1000000 tonnes",
  "≥ 10000000 to < 100000000 tonnes",
  "≥ 10000000 tonnes",
  "≥ 100000000 to < 1000000000 tonnes"
)

eu_reg_subs <- eu_reg_subs %>%
  mutate(
    Annex = case_when(
      total_tonnage_band %in% annex_vii  ~ "VII",
      total_tonnage_band %in% annex_viii ~ "VIII",
      total_tonnage_band %in% annex_ix   ~ "IX",
      total_tonnage_band %in% annex_x    ~ "X",
      TRUE ~ NA_character_
    )
  )

# Check assignments
eu_reg_subs %>%
  count(Annex, sort = TRUE)

## TEST

library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(here)
library(readr)

# =============================================================================
# Purpose:
# Summarize hazard data availability by REACH Annex
#   - number of in vivo study records per chemical
#   - number of in vitro assay records per chemical
#
# Output:
#   Boxplot of record/assay counts by REACH Annex
# =============================================================================


# =============================================================================
# 1. Create chemical-level REACH Annex table
# =============================================================================

reach_annex_chem <- eu_reg_subs %>%
  filter(!is.na(Annex), cas_number != "-") %>%
  distinct(cas_number, Annex) %>%
  mutate(
    Annex = factor(
      Annex,
      levels = c("VII", "VIII", "IX", "X"),
      ordered = TRUE
    )
  )

# Optional check
reach_annex_chem %>%
  count(Annex, name = "n_reach_registered_chemicals")


# =============================================================================
# 2. Join REACH Annex chemicals to ROCKS identifiers
# =============================================================================
# Assumes rocks_ids contains cas_number and cichem_id.
# If your CAS column is named differently, update the join below.
rocks_ids<-readRDS(here::here("data","derived","rocks.rds"))

reach_annex_ids <- reach_annex_chem %>%
  rename(casrn=cas_number)%>%
  left_join(
    rocks_ids %>%
      distinct(casrn, cichem_id, dtxsid),
    by = "casrn"
  )

# Quick QC
reach_annex_ids %>%
  summarise(
    n_reach_cas = n_distinct(casrn),
    n_with_cichem_id = n_distinct(cichem_id[!is.na(cichem_id)]),
    n_with_dtxsid = n_distinct(dtxsid[!is.na(dtxsid)])
  )


# =============================================================================
# 3. Count in vivo studies per chemical
# =============================================================================

# dose response values only
tam <- readRDS(
  here("data", "derived", "tam.rds")
) %>%
  filter(toxval_type_supercategory=="Dose Response Summary Value")

invivo_counts <- tam %>%
  filter(!is.na(cichem_id)) %>%
  group_by(cichem_id) %>%
  summarise(
    n_invivo_studies = n(),
    .groups = "drop"
  )


# =============================================================================
# 4. Count in vitro assays per chemical
# =============================================================================
#
b_exec_sum <- readRDS(
  here("data", "derived", "b_exec_sum.rds")
)

invitro_counts <- b_exec_sum %>%
  filter(!is.na(cichem_id)) %>%
  select(
    cichem_id,
    n_invitro_assays = n_assay_endpoints_tested
  )


# =============================================================================
# 5. Build chemical-level Annex evidence table
# =============================================================================

annex_evidence_counts <- reach_annex_ids %>%
  filter(!is.na(cichem_id)) %>%
  distinct(cichem_id, Annex) %>%
  left_join(invivo_counts, by = "cichem_id") %>%
  left_join(invitro_counts, by = "cichem_id") %>%
  mutate(
    n_invivo_studies = replace_na(n_invivo_studies, 0L),
    n_invitro_assays = replace_na(n_invitro_assays, 0L)
  )

# QC: one row per REACH chemical/identifier/Annex
annex_evidence_counts %>%
  count(Annex, name = "n_chemicals")

# =============================================================================
# 6B. Data availability summary by Annex
# =============================================================================

annex_availability <- annex_evidence_counts %>%
  mutate(
    has_invivo = n_invivo_studies > 0,
    has_invitro = n_invitro_assays > 0,
    has_either = has_invivo | has_invitro,
    has_both = has_invivo & has_invitro
  ) %>%
  group_by(Annex) %>%
  summarise(
    n_chemicals = n_distinct(cichem_id),

    n_with_invivo = sum(has_invivo),
    pct_with_invivo = n_with_invivo / n_chemicals * 100,

    n_with_invitro = sum(has_invitro),
    pct_with_invitro = n_with_invitro / n_chemicals * 100,

    n_with_either = sum(has_either),
    pct_with_either = n_with_either / n_chemicals * 100,

    n_with_both = sum(has_both),
    pct_with_both = n_with_both / n_chemicals * 100,

    .groups = "drop"
  )

print(annex_availability)


# =============================================================================
# 6. Summary statistics by Annex
# =============================================================================

annex_summary <- annex_evidence_counts %>%
  group_by(Annex) %>%
  summarise(
    n_chemicals = n_distinct(cichem_id),

    median_invivo_studies = median(n_invivo_studies),
    mean_invivo_studies = mean(n_invivo_studies),
    pct_with_invivo = mean(n_invivo_studies > 0) * 100,

    median_invitro_assays = median(n_invitro_assays),
    mean_invitro_assays = mean(n_invitro_assays),
    pct_with_invitro = mean(n_invitro_assays > 0) * 100,

    median_invitro_nonzero =
      median(
        n_invitro_assays[n_invitro_assays > 0],
        na.rm = TRUE
      ),

    p90_invitro =
      quantile(
        n_invitro_assays,
        0.90,
        na.rm = TRUE
      ),

    .groups = "drop"
  )

print(annex_summary)


# =============================================================================
# 7. Reshape for plotting
# =============================================================================

plot_dat <- annex_evidence_counts %>%
  select(cichem_id, Annex, n_invivo_studies, n_invitro_assays) %>%
  pivot_longer(
    cols = c(n_invivo_studies, n_invitro_assays),
    names_to = "evidence_stream",
    values_to = "n_records"
  ) %>%
  mutate(
    evidence_stream = recode(
      evidence_stream,
      n_invivo_studies = "In vivo studies",
      n_invitro_assays = "In vitro assays"
    ),
    evidence_stream = factor(
      evidence_stream,
      levels = c("In vivo studies", "In vitro assays")
    )
  )

# =============================================================================
# 7B. Plot percent of chemicals with data by Annex
# =============================================================================

## Version with annotations

availability_plot_pct <- annex_summary %>%
  select(
    Annex,

    pct_with_invivo,
    pct_with_invitro,

    mean_invivo_studies,
    median_invivo_studies,

    mean_invitro_assays,
    median_invitro_assays
  ) %>%
  pivot_longer(
    cols = c(pct_with_invivo, pct_with_invitro),
    names_to = "evidence_stream",
    values_to = "percent_with_data"
  ) %>%
  mutate(
    evidence_stream = recode(
      evidence_stream,
      pct_with_invivo = "In vivo DRSV",
      pct_with_invitro = "Targeted"
    ),
    evidence_stream = factor(
      evidence_stream,
      levels = c(
        "In vivo DRSV",
        "Targeted"
      )
    ),

    mean_count = case_when(
      evidence_stream == "In vivo DRSV" ~ mean_invivo_studies,
      evidence_stream == "Targeted" ~ mean_invitro_assays,
      TRUE ~ NA_real_
    ),

    median_count = case_when(
      evidence_stream == "In vivo DRSV" ~ median_invivo_studies,
      evidence_stream == "Targeted" ~ median_invitro_assays,
      TRUE ~ NA_real_
    ),

    count_label = paste0(
      "mean=", round(mean_count, 1),
      "\nmedian=", round(median_count, 1)
    )
  )


evidence_stream_colors <- c(
  "In vivo DRSV" = ULRIcolors[3],
  "Targeted" = ULRIcolors[6]
)

p_annex_availability_pct <- ggplot(
  availability_plot_pct,
  aes(
    x = Annex,
    y = percent_with_data,
    fill = evidence_stream
  )
) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65,
    color = "grey25",
    linewidth = 0.2
  ) +

  # Mean and median count labels on bars
  geom_text(
    aes(
      label = count_label,
      group = evidence_stream
    ),
    position = position_dodge(width = 0.75),
    vjust = -0.25,
    size = 5.0,
    lineheight = 0.9,
    color = "black"
  ) +

  scale_fill_manual(
    values = evidence_stream_colors
  ) +

  scale_y_continuous(
    limits = c(0, 110),
    breaks = seq(0, 100, 20),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    title = "Hazard data availability by REACH Annex",
    x = "REACH Annex",
    y = "Chemicals with data (%)",
    fill = "Evidence stream"
  ) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

p_annex_availability_pct

# =============================================================================
# 8. Boxplot by REACH Annex
# =============================================================================

annex_colors <- c(
  "VII" = ULRIcolors[7],
  "VIII" = ULRIcolors[5],
  "IX" = ULRIcolors[4],
  "X" = ULRIcolors[1]
)

p_annex_boxplot <- ggplot(
  plot_dat,
  aes(
    x = Annex,
    y = n_records,
    fill = Annex
  )
) +
  geom_boxplot(
    width = 0.65,
    outlier.shape = NA,
    alpha = 0.85
  ) +
  geom_jitter(
    aes(color = Annex),
    width = 0.15,
    height = 0,
    alpha = 0.12,
    size = 0.8,
    show.legend = FALSE
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 21,
    size = 2.8,
    fill = "white",
    color = "black"
  ) +
  facet_wrap(
    ~ evidence_stream,
    scales = "free_y"
  ) +

  scale_fill_manual(
    values = annex_colors
  ) +
  scale_color_manual(
    values = annex_colors
  ) +

  scale_y_continuous(
    trans = "log1p",
    breaks = c(0, 1, 5, 10, 25, 50, 100, 250, 500, 1000)
  ) +
  labs(
    title = "Hazard data availability by REACH Annex",
    subtitle = "Distribution of in vivo study and in vitro assay counts per chemical",
    x = "REACH Annex",
    y = "Number of records per chemical, log1p scale",
    fill = "REACH Annex"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )


p_annex_boxplot


# =============================================================================
# 9. Save figure
# =============================================================================

ggsave(
  filename = here("outputs", "figures", "reach_annex_invivo_invitro_boxplot.png"),
  plot = p_annex_boxplot,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  filename = here("outputs", "figures", "reach_annex_invivo_invitro_boxplot.pdf"),
  plot = p_annex_boxplot,
  width = 8,
  height = 5
)


# =============================================================================
# TSCA Active Inventory
# =============================================================================

tsca_active <- read_xlsx(
  here("data", "raw", "Chemical List TSCA_ACTIVE_NCTI_0125-2026-07-27.xlsx")
) %>%
  clean_names() %>%
  transmute(
    cas_number = cas_rn,   # update to actual column name
    inventory = "TSCA Active"
  )

reach_inventory <- eu_reg_subs %>%
  filter(cas_number != "-") %>%
  distinct(cas_number) %>%
  transmute(
    cas_number,
    inventory = "EU REACH"
  )

inventory_membership <- bind_rows(
  reach_inventory,
  tsca_active
) %>%
  distinct()