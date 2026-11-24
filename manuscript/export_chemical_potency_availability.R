# ==============================================================================
# export_chemical_potency_availability.R
#
# Purpose:
#   Create a one-row-per-chemical workbook summarizing availability of:
#
#     1. InvitroDB data
#     2. InvitroDB AEDs
#     3. Broad Profiling data
#     4. Broad Profiling AEDs
#     5. All TAM DRSVs
#     6. Selected oral TAM DRSVs
#
# Identifier policy:
#   The ROCKS inventory is defined at the cichem_id level.
#   Therefore, cichem_id is the primary key for joins and summaries.
#
# TAM DRSV definitions:
#   drsvs_all.rds:
#     All qualifying TAM DRSVs before route/unit restrictions and
#     study-group selection.
#
#   hazard_landscape_drsvs.rds:
#     Selected oral DRSVs after the oral/unit restrictions and
#     one-DRSV-per-chemical-study-group selection.
#
# Output:
#   manuscript/chemical_potency_availability.xlsx
# ==============================================================================

# ==============================================================================
# PACKAGES
# ==============================================================================

library(dplyr)
library(readr)
library(stringr)
library(janitor)
library(here)
library(openxlsx)


# ==============================================================================
# SETUP
# ==============================================================================

dir.create(
  here::here("manuscript"),
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

normalize_id_types <- function(df) {
  df %>%
    mutate(
      across(
        any_of(
          c(
            "cichem_id",
            "dtxsid",
            "cicomp_id",
            "inchi",
            "inchikey",
            "smiles"
          )
        ),
        as.character
      )
    )
}


standardize_dtxsid_column <- function(df) {
  df <- df %>%
    janitor::clean_names()

  dtxsid_candidates <- c(
    "dtxsid",
    "dsstox_substance_id",
    "dsstox_sid",
    "substance_id"
  )

  existing_dtxsid_col <- intersect(
    dtxsid_candidates,
    names(df)
  )

  if (
    !"dtxsid" %in% names(df) &&
      length(existing_dtxsid_col) > 0
  ) {
    names(df)[
      names(df) == existing_dtxsid_col[[1]]
    ] <- "dtxsid"
  }

  df %>%
    normalize_id_types()
}


add_cichem_id_from_dtxsid <- function(df, id_crosswalk) {
  df <- df %>%
    standardize_dtxsid_column()

  if (!"dtxsid" %in% names(df)) {
    warning(
      "No dtxsid-like column found. cichem_id was not added."
    )

    return(df)
  }

  dtxsid_map <- id_crosswalk %>%
    select(
      cichem_id,
      dtxsid
    ) %>%
    filter(
      !is.na(dtxsid),
      dtxsid != ""
    ) %>%
    distinct()

  if ("cichem_id" %in% names(df)) {
    df <- df %>%
      left_join(
        dtxsid_map %>%
          rename(
            cichem_id_from_dtxsid = cichem_id
          ),
        by = "dtxsid"
      ) %>%
      mutate(
        cichem_id = coalesce(
          as.character(cichem_id),
          as.character(cichem_id_from_dtxsid)
        )
      ) %>%
      select(-cichem_id_from_dtxsid)
  } else {
    df <- df %>%
      left_join(
        dtxsid_map,
        by = "dtxsid"
      )
  }

  df %>%
    normalize_id_types()
}


first_nonmissing <- function(x) {
  x <- as.character(x)

  x <- x[
    !is.na(x) &
      stringr::str_squish(x) != ""
  ]

  if (length(x) == 0) {
    return(NA_character_)
  }

  x[[1]]
}


collapse_values <- function(x) {
  x <- as.character(x)

  x <- x[
    !is.na(x) &
      stringr::str_squish(x) != ""
  ]

  x <- sort(unique(x))

  if (length(x) == 0) {
    return(NA_character_)
  }

  paste(
    x,
    collapse = " | "
  )
}


collapse_numeric_values <- function(x, digits = 6) {
  x <- suppressWarnings(
    as.numeric(x)
  )

  x <- x[
    !is.na(x) &
      is.finite(x) &
      x > 0
  ]

  x <- sort(unique(x))

  if (length(x) == 0) {
    return(NA_character_)
  }

  formatted_values <- format(
    signif(
      x,
      digits = digits
    ),
    scientific = TRUE,
    trim = TRUE
  )

  paste(
    formatted_values,
    collapse = " | "
  )
}

require_columns <- function(
  df,
  required_columns,
  dataset_name
) {
  missing_columns <- setdiff(
    required_columns,
    names(df)
  )

  if (length(missing_columns) > 0) {
    stop(
      paste0(
        dataset_name,
        " is missing the following required columns: ",
        paste(
          missing_columns,
          collapse = ", "
        )
      )
    )
  }

  invisible(TRUE)
}


# ==============================================================================
# READ CORE DATASETS
# ==============================================================================

message("\nReading core datasets...")


rocks <- readRDS(
  here::here(
    "data",
    "derived",
    "rocks.rds"
  )
) %>%
  janitor::clean_names() %>%
  normalize_id_types()


id_crosswalk <- readRDS(
  here::here(
    "data",
    "derived",
    "id_crosswalk.rds"
  )
) %>%
  janitor::clean_names() %>%
  normalize_id_types()


b_exec_sum <- readRDS(
  here::here(
    "data",
    "derived",
    "b_exec_sum.rds"
  )
) %>%
  janitor::clean_names() %>%
  normalize_id_types()


b_mc7_aed <- readRDS(
  here::here(
    "data",
    "derived",
    "b_mc7_aed.rds"
  )
) %>%
  janitor::clean_names() %>%
  normalize_id_types()


bp <- readRDS(
  here::here(
    "data",
    "derived",
    "bp.rds"
  )
) %>%
  janitor::clean_names() %>%
  normalize_id_types()


bp_aed <- readr::read_csv(
  here::here(
    "data",
    "raw",
    "rocks_profile_pod_table_aeds_08sept26.csv"
  ),
  show_col_types = FALSE
) %>%
  janitor::clean_names() %>%
  standardize_dtxsid_column()


drsvs_all <- readRDS(
  here::here(
    "data",
    "derived",
    "drsvs_all.rds"
  )
) %>%
  janitor::clean_names() %>%
  normalize_id_types()


tam_drsv <- readRDS(
  here::here(
    "data",
    "derived",
    "hazard_landscape_drsvs.rds"
  )
) %>%
  janitor::clean_names() %>%
  normalize_id_types()


# ==============================================================================
# ADD OR COMPLETE cichem_id
# ==============================================================================

b_exec_sum <- add_cichem_id_from_dtxsid(
  b_exec_sum,
  id_crosswalk
)


b_mc7_aed <- add_cichem_id_from_dtxsid(
  b_mc7_aed,
  id_crosswalk
)


bp <- add_cichem_id_from_dtxsid(
  bp,
  id_crosswalk
)


bp_aed <- add_cichem_id_from_dtxsid(
  bp_aed,
  id_crosswalk
)


drsvs_all <- add_cichem_id_from_dtxsid(
  drsvs_all,
  id_crosswalk
)


tam_drsv <- add_cichem_id_from_dtxsid(
  tam_drsv,
  id_crosswalk
)


# ==============================================================================
# CHECK REQUIRED IDENTIFIER COLUMNS
# ==============================================================================

require_columns(
  rocks,
  "cichem_id",
  "rocks.rds"
)


require_columns(
  b_exec_sum,
  "cichem_id",
  "b_exec_sum.rds"
)


require_columns(
  b_mc7_aed,
  "cichem_id",
  "b_mc7_aed.rds"
)


require_columns(
  bp,
  "cichem_id",
  "bp.rds"
)


require_columns(
  bp_aed,
  "cichem_id",
  "rocks_profile_pod_table_aeds_08sept26.csv"
)


require_columns(
  drsvs_all,
  c(
    "cichem_id",
    "toxval_numeric"
  ),
  "drsvs_all.rds"
)


require_columns(
  tam_drsv,
  c(
    "cichem_id",
    "toxval_numeric"
  ),
  "hazard_landscape_drsvs.rds"
)

# ==============================================================================
# DATA AVAILABILITY FLAGS
# ==============================================================================

message("\nBuilding data availability flags...")


# ------------------------------------------------------------------------------
# InvitroDB: any data
# ------------------------------------------------------------------------------

invitrodb_data_flag <- b_exec_sum %>%
  filter(
    !is.na(cichem_id),
    cichem_id != ""
  ) %>%
  distinct(cichem_id) %>%
  mutate(
    has_invitrodb_data = "Yes"
  )


# ------------------------------------------------------------------------------
# InvitroDB: AED available
#
# Availability is based specifically on p5_aed_ac50_sumcl50.
# ------------------------------------------------------------------------------

invitrodb_aed_flag <- b_mc7_aed %>%
  filter(
    !is.na(cichem_id),
    cichem_id != "",
    !is.na(p5_aed_ac50_sumcl50),
    p5_aed_ac50_sumcl50 > 0
  ) %>%
  distinct(cichem_id) %>%
  mutate(
    has_invitrodb_aed = "Yes"
  )


# ------------------------------------------------------------------------------
# Broad Profiling: any data
# ------------------------------------------------------------------------------

broad_profiling_data_flag <- bp %>%
  filter(
    !is.na(cichem_id),
    cichem_id != ""
  ) %>%
  distinct(cichem_id) %>%
  mutate(
    has_broad_profiling_data = "Yes"
  )


# ------------------------------------------------------------------------------
# Broad Profiling: AED available
#
# Availability is based specifically on min_pod_aed_sumcl50.
# ------------------------------------------------------------------------------

broad_profiling_aed_flag <- bp_aed %>%
  filter(
    !is.na(cichem_id),
    cichem_id != "",
    !is.na(min_pod_aed_sumcl50),
    min_pod_aed_sumcl50 > 0
  ) %>%
  distinct(cichem_id) %>%
  mutate(
    has_broad_profiling_aed = "Yes"
  )


# ------------------------------------------------------------------------------
# TAM: any DRSV
#
# Based on the complete DRSV dataset.
# ------------------------------------------------------------------------------

tam_drsv_all_flag <- drsvs_all %>%
  filter(
    !is.na(cichem_id),
    cichem_id != "",
    !is.na(toxval_numeric),
    toxval_numeric > 0
  ) %>%
  distinct(cichem_id) %>%
  mutate(
    has_any_tam_drsv = "Yes"
  )


# ------------------------------------------------------------------------------
# TAM: selected oral DRSV
#
# Based on the final hazard-landscape DRSV dataset.
# ------------------------------------------------------------------------------

tam_selected_drsv_flag <- tam_drsv %>%
  filter(
    !is.na(cichem_id),
    cichem_id != "",
    !is.na(toxval_numeric),
    toxval_numeric > 0
  ) %>%
  distinct(cichem_id) %>%
  mutate(
    has_selected_oral_tam_drsv = "Yes"
  )


# ==============================================================================
# BUILD ONE-ROW-PER-CHEMICAL IDENTIFIER TABLE
# ==============================================================================

message("\nBuilding chemical identifier table...")


tam_chemical_names <- drsvs_all %>%
  filter(
    !is.na(cichem_id),
    cichem_id != ""
  ) %>%
  group_by(cichem_id) %>%
  summarise(
    chemical_name_tam = if ("chemical_name" %in% names(drsvs_all)) {
      first_nonmissing(chemical_name)
    } else {
      NA_character_
    },
    .groups = "drop"
  )


chemical_identifier_table <- rocks %>%
  filter(
    !is.na(cichem_id),
    cichem_id != ""
  ) %>%
  group_by(cichem_id) %>%
  summarise(
    dtxsid = if ("dtxsid" %in% names(rocks)) {
      first_nonmissing(dtxsid)
    } else {
      NA_character_
    },

    casrn = if ("casrn" %in% names(rocks)) {
      first_nonmissing(casrn)
    } else {
      NA_character_
    },

    chemical_name_rocks = if ("chemical_name" %in% names(rocks)) {
      first_nonmissing(chemical_name)
    } else {
      NA_character_
    },

    smiles = if ("smiles" %in% names(rocks)) {
      first_nonmissing(smiles)
    } else {
      NA_character_
    },

    inchi = if ("inchi" %in% names(rocks)) {
      first_nonmissing(inchi)
    } else {
      NA_character_
    },

    .groups = "drop"
  ) %>%
  left_join(
    tam_chemical_names,
    by = "cichem_id"
  ) %>%
  mutate(
    chemical_name = coalesce(
      chemical_name_rocks,
      chemical_name_tam
    )
  ) %>%
  select(
    cichem_id,
    dtxsid,
    casrn,
    chemical_name,
    smiles,
    inchi
  )


# ==============================================================================
# BUILD FINAL ONE-ROW-PER-CHEMICAL TABLE
# ==============================================================================

message("\nBuilding final chemical availability table...")


chemical_potency_availability <- chemical_identifier_table %>%

  left_join(
    invitrodb_data_flag,
    by = "cichem_id"
  ) %>%

  left_join(
    invitrodb_aed_flag,
    by = "cichem_id"
  ) %>%

  left_join(
    broad_profiling_data_flag,
    by = "cichem_id"
  ) %>%

  left_join(
    broad_profiling_aed_flag,
    by = "cichem_id"
  ) %>%

  left_join(
    tam_drsv_all_flag,
    by = "cichem_id"
  ) %>%

  left_join(
    tam_selected_drsv_flag,
    by = "cichem_id"
  ) %>%

  mutate(
    across(
      starts_with("has_"),
      ~ coalesce(.x, "No")
    ),

    has_any_potency = if_else(
      has_invitrodb_aed == "Yes" |
        has_broad_profiling_aed == "Yes" |
        has_any_tam_drsv == "Yes",
      "Yes",
      "No"
    )
  ) %>%

  select(
    cichem_id,
    dtxsid,
    casrn,
    chemical_name,

    has_invitrodb_data,
    has_invitrodb_aed,

    has_broad_profiling_data,
    has_broad_profiling_aed,

    has_any_tam_drsv,
    has_selected_oral_tam_drsv,

    has_any_potency,

    smiles,
    inchi
  ) %>%

  arrange(
    chemical_name,
    cichem_id
  )

chemical_potency_summary <- tibble::tibble(
  metric = c(
    "Total ROCKS chemicals",
    "Chemicals with any InvitroDB data",
    "Chemicals with InvitroDB AED",
    "Chemicals with any Broad Profiling data",
    "Chemicals with Broad Profiling AED",
    "Chemicals with any TAM DRSV",
    "Chemicals with selected oral TAM DRSV",
    "Chemicals with any potency information",
    "Chemicals with no potency information"
  ),

  n_chemicals = c(
    nrow(chemical_potency_availability),

    sum(
      chemical_potency_availability$has_invitrodb_data == "Yes"
    ),

    sum(
      chemical_potency_availability$has_invitrodb_aed == "Yes"
    ),

    sum(
      chemical_potency_availability$has_broad_profiling_data == "Yes"
    ),

    sum(
      chemical_potency_availability$has_broad_profiling_aed == "Yes"
    ),

    sum(
      chemical_potency_availability$has_any_tam_drsv == "Yes"
    ),

    sum(
      chemical_potency_availability$has_selected_oral_tam_drsv == "Yes"
    ),

    sum(
      chemical_potency_availability$has_any_potency == "Yes"
    ),

    sum(
      chemical_potency_availability$has_any_potency == "No"
    )
  )
) %>%

  mutate(
    percent_of_rocks = round(
      100 * n_chemicals / nrow(chemical_potency_availability),
      1
    )
  )

# ==============================================================================
# WRITE EXCEL WORKBOOK
# ==============================================================================

message("\nWriting Excel workbook...")


potency_export_path <- here::here(
  "manuscript",
  "chemical_potency_availability.xlsx"
)


workbook <- openxlsx::createWorkbook()


# ------------------------------------------------------------------------------
# Chemical availability sheet
# ------------------------------------------------------------------------------

openxlsx::addWorksheet(
  workbook,
  "chemical_availability"
)

openxlsx::writeDataTable(
  workbook,
  sheet = "chemical_availability",
  x = chemical_potency_availability,
  tableStyle = "TableStyleMedium2"
)

openxlsx::freezePane(
  workbook,
  sheet = "chemical_availability",
  firstRow = TRUE
)

openxlsx::setColWidths(
  workbook,
  sheet = "chemical_availability",
  cols = 1:ncol(chemical_potency_availability),
  widths = "auto"
)


# ------------------------------------------------------------------------------
# Summary sheet
# ------------------------------------------------------------------------------

openxlsx::addWorksheet(
  workbook,
  "summary"
)

openxlsx::writeDataTable(
  workbook,
  sheet = "summary",
  x = chemical_potency_summary,
  tableStyle = "TableStyleMedium2"
)

openxlsx::freezePane(
  workbook,
  sheet = "summary",
  firstRow = TRUE
)

openxlsx::setColWidths(
  workbook,
  sheet = "summary",
  cols = 1:ncol(chemical_potency_summary),
  widths = "auto"
)


# ------------------------------------------------------------------------------
# Save workbook
# ------------------------------------------------------------------------------

openxlsx::saveWorkbook(
  workbook,
  file = potency_export_path,
  overwrite = TRUE
)


message(
  "\nChemical potency availability workbook saved to: ",
  potency_export_path
)
