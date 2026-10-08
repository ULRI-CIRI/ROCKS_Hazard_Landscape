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

library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(here)
library(readr)

# =============================================================================
# Summarize hazard data availability by REACH Annex
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
# 4b. Count broad-profiling unique PODs per chemical
# =============================================================================

broad_profiling <- readxl::read_xlsx(
  here::here(
    "data",
    "raw",
    "rocks_profile_pod_table.xlsx"
  )
)

broad_depth <- broad_profiling %>%
  transmute(
    cichem_id = cichem_id,
    n_broad_profiling_unique_pods = as.numeric(n_pods_uniq)
  ) %>%
  filter(
    !is.na(cichem_id),
    cichem_id != ""
  ) %>%
  group_by(cichem_id) %>%
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

# =============================================================================
# 5. Build chemical-level Annex evidence table
# =============================================================================

annex_evidence_counts <- reach_annex_ids %>%
  filter(!is.na(cichem_id)) %>%
  distinct(cichem_id, Annex) %>%
  left_join(
    invivo_counts,
    by = "cichem_id"
  ) %>%
  left_join(
    invitro_counts,
    by = "cichem_id"
  ) %>%
  left_join(
    broad_depth,
    by = "cichem_id"
  ) %>%
  mutate(
    n_invivo_studies = replace_na(
      n_invivo_studies,
      0L
    ),
    n_invitro_assays = replace_na(
      n_invitro_assays,
      0
    ),
    n_broad_profiling_unique_pods = replace_na(
      n_broad_profiling_unique_pods,
      0
    )
  )


# =============================================================================
# 6B. Data availability summary by Annex
# =============================================================================

annex_availability <- annex_evidence_counts %>%
  mutate(
    has_invivo = n_invivo_studies > 0,
    has_invitro = n_invitro_assays > 0,
    has_broad = n_broad_profiling_unique_pods > 0,

    has_either =
      has_invivo |
      has_invitro |
      has_broad,

    has_all_three =
      has_invivo &
      has_invitro &
      has_broad
  ) %>%
  group_by(Annex) %>%
  summarise(
    n_chemicals = n_distinct(cichem_id),

    n_with_invivo = sum(has_invivo),
    pct_with_invivo =
      100 * n_with_invivo / n_chemicals,

    n_with_invitro = sum(has_invitro),
    pct_with_invitro =
      100 * n_with_invitro / n_chemicals,

    n_with_broad = sum(has_broad),
    pct_with_broad =
      100 * n_with_broad / n_chemicals,

    n_with_any = sum(has_either),
    pct_with_any =
      100 * n_with_any / n_chemicals,

    n_with_all_three = sum(has_all_three),
    pct_with_all_three =
      100 * n_with_all_three / n_chemicals,

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

    median_invivo_studies =
      median(n_invivo_studies),
    mean_invivo_studies =
      mean(n_invivo_studies),
    pct_with_invivo =
      mean(n_invivo_studies > 0) * 100,

    median_invitro_assays =
      median(n_invitro_assays),
    mean_invitro_assays =
      mean(n_invitro_assays),
    pct_with_invitro =
      mean(n_invitro_assays > 0) * 100,

    median_invitro_nonzero =
      median(
        n_invitro_assays[
          n_invitro_assays > 0
        ],
        na.rm = TRUE
      ),

    p90_invitro =
      quantile(
        n_invitro_assays,
        0.90,
        na.rm = TRUE
      ),

    median_broad_pods =
      median(n_broad_profiling_unique_pods),
    mean_broad_pods =
      mean(n_broad_profiling_unique_pods),
    pct_with_broad =
      mean(n_broad_profiling_unique_pods > 0) * 100,

    median_broad_nonzero =
      median(
        n_broad_profiling_unique_pods[
          n_broad_profiling_unique_pods > 0
        ],
        na.rm = TRUE
      ),

    p90_broad =
      quantile(
        n_broad_profiling_unique_pods,
        0.90,
        na.rm = TRUE
      ),

    .groups = "drop"
  )

print(annex_summary)

# =============================================================================
# 6C. Build threshold table
# =============================================================================

depth_by_annex <- annex_evidence_counts %>%
  transmute(
    chem_id = cichem_id,
    Annex,
    n_in_vivo_drsv_studies =
      n_invivo_studies,
    n_in_vitro_targeted_assay_endpoints =
      n_invitro_assays,
    n_broad_profiling_unique_pods =
      n_broad_profiling_unique_pods
  )

threshold_df <- depth_by_annex %>%
  select(
    chem_id,
    Annex,
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
      n_in_vivo_drsv_studies =
        "In vivo",
      n_in_vitro_targeted_assay_endpoints =
        "Targeted bioactivity",
      n_broad_profiling_unique_pods =
        "Broad profiling"
    ),
    count = replace_na(
      as.numeric(count),
      0
    )
  )

threshold_lookup <- tibble::tribble(
  ~stream,                 ~threshold_label, ~threshold_value, ~threshold_order,
  "In vivo",               "≥1",             1,                1,
  "In vivo",               "≥2",             2,                2,
  "In vivo",               "≥5",             5,                3,
  "Targeted bioactivity",  "≥1",             1,                1,
  "Targeted bioactivity",  "≥100",           100,              2,
  "Targeted bioactivity",  "≥300",           300,              3,
  "Broad profiling",       "≥1",             1,                1,
  "Broad profiling",       "≥2",             2,                2,
  "Broad profiling",       "≥4",             4,                3
)

# =============================================================================
# 7. Percent of chemicals meeting evidence thresholds
# =============================================================================

annex_totals <- threshold_df %>%
  distinct(Annex, chem_id) %>%
  count(Annex, name = "n_annex_total")

annex_x_labels <- annex_totals %>%
  mutate(
    x_label = paste0(
      Annex,
      "\n(",
      scales::comma(n_annex_total),
      ")"
    )
  ) %>%
  select(
    Annex,
    x_label
  ) %>%
  tibble::deframe()

threshold_counts_raw <- threshold_df %>%
  left_join(
    threshold_lookup,
    by = "stream"
  ) %>%
  filter(
    count >= threshold_value
  ) %>%
  group_by(
    Annex,
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
  distinct(Annex, stream) %>%
  left_join(
    threshold_lookup,
    by = "stream"
  ) %>%
  left_join(
    threshold_counts_raw,
    by = c(
      "Annex",
      "stream",
      "threshold_label",
      "threshold_value",
      "threshold_order"
    )
  ) %>%
  left_join(
    annex_totals,
    by = "Annex"
  ) %>%
  mutate(
    n_chems = replace_na(n_chems, 0),
    pct_annex = 100 * n_chems / n_annex_total
  )

# =============================================================================
# 8. Assessment values by Annex
# =============================================================================

tam_assess <- readRDS(
  here("data", "derived", "tam_assess.rds")
)

assessment_chems <- tam_assess %>%
  filter(
    population == "General Population"
  ) %>%
  distinct(
    chem_id = dtxsid
  )

assessment_by_annex <- reach_annex_ids %>%
  filter(!is.na(dtxsid)) %>%
  distinct(dtxsid, Annex) %>%
  mutate(
    has_assessment =
      dtxsid %in% assessment_chems$chem_id
  ) %>%
  group_by(Annex) %>%
  summarise(
    n_total = n_distinct(dtxsid),
    n_assessment =
      n_distinct(dtxsid[has_assessment]),
    pct_assessment =
      100 * n_assessment / n_total,
    .groups = "drop"
  ) %>%
  mutate(
    stream = "In vivo"
  )

threshold_invivo <- threshold_counts %>%
  filter(stream == "In vivo") %>%
  mutate(
    threshold_label = factor(
      threshold_label,
      levels = c("≥1", "≥2", "≥5")
    )
  )

threshold_targeted <- threshold_counts %>%
  filter(stream == "Targeted") %>%
  mutate(
    threshold_label = factor(
      threshold_label,
      levels = c("≥1", "≥100", "≥300")
    )
  )
# =============================================================================
# 8B. Order panels and prepare stream-specific plotting data
# =============================================================================

stream_order <- c(
  "In vivo",
  "Targeted bioactivity",
  "Broad profiling"
)

threshold_counts <- threshold_counts %>%
  mutate(
    Annex = factor(
      Annex,
      levels = c(
        "VII",
        "VIII",
        "IX",
        "X"
      ),
      ordered = TRUE
    ),
    stream = factor(
      stream,
      levels = stream_order
    )
  )

threshold_invivo <- threshold_counts %>%
  filter(
    stream == "In vivo"
  ) %>%
  mutate(
    threshold_label = factor(
      threshold_label,
      levels = c(
        "≥1",
        "≥2",
        "≥5"
      )
    )
  )

threshold_targeted <- threshold_counts %>%
  filter(
    stream == "Targeted bioactivity"
  ) %>%
  mutate(
    threshold_label = factor(
      threshold_label,
      levels = c(
        "≥1",
        "≥100",
        "≥300"
      )
    )
  )

threshold_broad <- threshold_counts %>%
  filter(
    stream == "Broad profiling"
  ) %>%
  mutate(
    threshold_label = factor(
      threshold_label,
      levels = c(
        "≥1",
        "≥2",
        "≥4"
      )
    )
  )
# =============================================================================
# 9. Plot
# =============================================================================

p_annex_thresholds <- ggplot() +

  # Targeted bioactivity assay endpoints
  geom_col(
    data = threshold_targeted,
    aes(
      x = Annex,
      y = pct_annex,
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
    breaks = c(
      "≥1",
      "≥100",
      "≥300"
    ),
    guide = guide_legend(
      order = 2,
      nrow = 1,
      title.position = "top"
    )
  ) +

  ggnewscale::new_scale_fill() +

  # In vivo DRSVs
  geom_col(
    data = threshold_invivo,
    aes(
      x = Annex,
      y = pct_annex,
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
    breaks = c(
      "≥1",
      "≥2",
      "≥5"
    ),
    guide = guide_legend(
      order = 1,
      nrow = 1,
      title.position = "top"
    )
  ) +

  ggnewscale::new_scale_fill() +

  # Broad-profiling unique PODs
  geom_col(
    data = threshold_broad,
    aes(
      x = Annex,
      y = pct_annex,
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
    breaks = c(
      "≥1",
      "≥2",
      "≥4"
    ),
    guide = guide_legend(
      order = 3,
      nrow = 1,
      title.position = "top"
    )
  ) +

  scale_x_discrete(
    labels = annex_x_labels
  ) +

  facet_wrap(
    ~stream,
    nrow = 1,
    drop = FALSE
  ) +

  scale_y_continuous(
    labels = scales::label_percent(
      scale = 1
    ),
    expand = expansion(
      mult = c(0, 0.06)
    )
  ) +

  labs(
    x = "REACH Annex",
    y = "Chemicals within Annex (%)"
  ) +

  theme_minimal(base_size = 16) +

  theme(
    strip.text = element_text(
      face = "bold"
    ),
    axis.text.x = element_text(
      lineheight = 0.9
    ),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal"
  )

p_annex_thresholds

out_dir <- here::here(
  "outputs",
  "figures"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  file.path(
    out_dir,
    "depth_by_reach_annex.png"
  ),
  p_annex_thresholds,
  width = 13,
  height = 5.5,
  dpi = 300
)

### COUNTS

# Total ROCKS inventory size
n_rocks <- rocks_ids %>%
  summarise(
    n_rocks = n_distinct(cichem_id)
  ) %>%
  pull(n_rocks)

# Count ROCKS chemicals in each Annex
rocks_annex_summary <- reach_annex_ids %>%
  filter(
    !is.na(cichem_id),
    !is.na(Annex)
  ) %>%
  distinct(
    cichem_id,
    Annex
  ) %>%
  count(
    Annex,
    name = "n_rocks_in_annex"
  ) %>%
  mutate(
    pct_rocks_inventory =
      100 * n_rocks_in_annex / n_rocks
  )

rocks_annex_summary

# =============================================================================
# Count and percent of chemicals in each REACH Annex that are in ROCKS
# =============================================================================

# Create unique ROCKS CAS-to-ID lookup
rocks_cas_lookup <- rocks_ids %>%
  filter(
    !is.na(casrn),
    !is.na(cichem_id)
  ) %>%
  distinct(
    casrn,
    cichem_id
  )

# Join REACH Annex chemicals to ROCKS by CAS RN
reach_annex_rocks_overlap <- reach_annex_chem %>%
  rename(
    casrn = cas_number
  ) %>%
  left_join(
    rocks_cas_lookup,
    by = "casrn"
  )

# Summarize ROCKS representation within each Annex
rocks_overlap_by_annex <- reach_annex_rocks_overlap %>%
  group_by(Annex) %>%
  summarise(
    n_reach_annex_chemicals = n_distinct(casrn),
    n_in_rocks = n_distinct(cichem_id[!is.na(cichem_id)]),
    pct_in_rocks = 100 * n_in_rocks / n_reach_annex_chemicals,
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0(
      scales::comma(n_in_rocks),
      " / ",
      scales::comma(n_reach_annex_chemicals),
      " (",
      round(pct_in_rocks, 1),
      "%)"
    )
  )

print(rocks_overlap_by_annex)

rocks_overlap_overall <- reach_annex_rocks_overlap %>%
  summarise(
    n_reach_annex_chemicals = n_distinct(casrn),
    n_in_rocks = n_distinct(cichem_id[!is.na(cichem_id)]),
    pct_in_rocks = 100 * n_in_rocks / n_reach_annex_chemicals
  ) %>%
  mutate(
    label = paste0(
      scales::comma(n_in_rocks),
      " / ",
      scales::comma(n_reach_annex_chemicals),
      " (",
      round(pct_in_rocks, 1),
      "%)"
    )
  )

rocks_overlap_overall

# =============================================================================
# QC: Is REACH-to-ROCKS missingness driven by missing identifiers?
# =============================================================================

reach_annex_mapping_qc <- reach_annex_chem %>%
  rename(
    casrn = cas_number
  ) %>%
  mutate(
    casrn_clean = trimws(casrn),
    has_usable_cas = !is.na(casrn_clean) &
      casrn_clean != "" &
      casrn_clean != "-"
  ) %>%
  left_join(
    rocks_ids %>%
      mutate(
        casrn_clean = trimws(casrn)
      ) %>%
      filter(
        !is.na(casrn_clean),
        casrn_clean != "",
        casrn_clean != "-"
      ) %>%
      distinct(
        casrn_clean,
        cichem_id,
        dtxsid
      ),
    by = "casrn_clean"
  ) %>%
  mutate(
    mapping_status = case_when(
      !has_usable_cas ~ "Missing/unusable CAS",
      has_usable_cas & !is.na(cichem_id) ~ "Mapped to ROCKS",
      has_usable_cas & is.na(cichem_id) ~ "Usable CAS, not mapped to ROCKS",
      TRUE ~ "Other"
    )
  )

mapping_qc_by_annex <- reach_annex_mapping_qc %>%
  distinct(
    Annex,
    casrn_clean,
    mapping_status,
    cichem_id
  ) %>%
  group_by(Annex, mapping_status) %>%
  summarise(
    n_chemicals = n_distinct(casrn_clean),
    .groups = "drop"
  ) %>%
  group_by(Annex) %>%
  mutate(
    n_annex_total = sum(n_chemicals),
    pct_annex = 100 * n_chemicals / n_annex_total
  ) %>%
  ungroup() %>%
  arrange(
    Annex,
    desc(n_chemicals)
  )

print(mapping_qc_by_annex)
