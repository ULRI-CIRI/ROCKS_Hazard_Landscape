# ============================================================
# CITox-TAM v1 source/source-table summary
#
# Input:
#   data/derived/tam.rds
#
# Output:
#   output/citox_tam_v1_source_summary.csv
# ============================================================

library(dplyr)
library(tidyr)
library(readr)
library(here)


# ------------------------------------------------------------
# 1. Load CITox-TAM data
# ------------------------------------------------------------

tam <- readRDS(
  here::here("data", "derived", "tam.rds")
)

tam_all_renamed <- tam %>%
  dplyr::mutate(
    type = dplyr::coalesce(
      response_type,
      assessment_value_name_abbreviated
    ),

    numeric = dplyr::coalesce(
      response_value,
      assessment_value
    ),

    units = dplyr::coalesce(
      response_units,
      assessment_units
    ),

    # Common high-level classification for release summary
    type_supercategory = dplyr::coalesce(
      response_type_supercategory,
      assessment_value_type
    )
  )


# ------------------------------------------------------------
# 2. Summarize records by source/source table/type
# ------------------------------------------------------------

source_summary_long <- tam_all_renamed %>%
  dplyr::mutate(
    type_supercategory = dplyr::coalesce(
      type_supercategory,
      "Unclassified"
    )
  ) %>%
  dplyr::group_by(
    source,
    source_table,
    type_supercategory
  ) %>%
  dplyr::summarise(
    n_records = dplyr::n_distinct(tamdb_id),
    n_chemicals = dplyr::n_distinct(chemical_id),
    .groups = "drop"
  )


# ------------------------------------------------------------
# Add totals for each source table
# ------------------------------------------------------------

source_table_totals <- tam_all_renamed %>%
  dplyr::filter(
    !is.na(source),
    !is.na(source_table)
  ) %>%
  dplyr::group_by(
    source,
    source_table
  ) %>%
  dplyr::summarise(
    n_records = dplyr::n_distinct(tamdb_id),
    n_chemicals = dplyr::n_distinct(chemical_id),
    .groups = "drop"
  )


# ------------------------------------------------------------
# Convert supercategory counts to wide format
# ------------------------------------------------------------

source_summary_wide <- source_summary_long %>%
  tidyr::pivot_wider(
    id_cols = c(
      source,
      source_table
    ),
    names_from = type_supercategory,
    values_from = c(
      n_records,
      n_chemicals
    ),
    values_fill = 0,
    names_glue = "{type_supercategory}_{.value}"
  ) %>%
  dplyr::left_join(
    source_table_totals,
    by = c(
      "source",
      "source_table"
    )
  ) %>%
  dplyr::relocate(
    source,
    source_table,
    n_records,
    n_chemicals
  ) %>%
  dplyr::arrange(
    source,
    source_table
  )


# ------------------------------------------------------------
# 5. QC
# ------------------------------------------------------------

cat(
  "Sources:",
  n_distinct(source_summary_wide$source),
  "\n"
)

cat(
  "Source tables:",
  nrow(source_summary_wide),
  "\n"
)

cat(
  "Expected for CITox-TAM v1: 64 sources and 85 source tables\n"
)


# ------------------------------------------------------------
# 6. Write release table
# ------------------------------------------------------------

write_csv(
  source_summary_wide,
  here::here(
    "manuscript",
    "citox_tam_v1_source_summary.csv"
  )
)
