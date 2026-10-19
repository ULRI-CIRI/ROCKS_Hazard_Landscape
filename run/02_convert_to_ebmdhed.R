# run/convert_oral_drsvs_to_eBMDHED_hazard_landscape.R

# Standalone hazard landscape workflow:
# 1. Map toxicological_effect to toxicological_effect_category
# 2. Map effect categories to conceptual models
# 3. Convert oral DRSVs to eBMD.HED
#
# Cancer is retained
# base_class is derived directly from toxval_type

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(tibble)
  library(here)
})

#---------------------------
# output folders
#---------------------------
dir.create(
  here::here("outputs", "effect_mapping", "tables"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  here::here("outputs", "ebmdhed_conversion", "tables"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  here::here("data", "derived"),
  recursive = TRUE,
  showWarnings = FALSE
)

#===============================================================================
# USER-EDITABLE SETTINGS
#===============================================================================

drop_missing_ebmd_hed <- TRUE

# Oral allometric scaling body weights, kg.
# CF2 is calculated as (human_bw / animal_bw)^(1/4).
# eBMD.HED = oral_value / CF2 for non-HED animal values.
body_weight_kg <- c(
  Human = 80,
  Rat = 0.25,
  Mouse = 0.03,
  Rabbit = 2.5,
  Dog = 10
)

CF2_lookup <- (body_weight_kg["Human"] / body_weight_kg)^(1 / 4)
CF2_lookup["Human"] <- 1

#===============================================================================
# HELPER FUNCTIONS
#===============================================================================

lookup_named_numeric <- function(x, lookup) {
  x_chr <- as.character(x)
  out <- unname(lookup[x_chr])
  as.numeric(out)
}

has_term <- function(x, pattern) {
  stringr::str_detect(
    dplyr::coalesce(as.character(x), ""),
    stringr::regex(pattern, ignore_case = TRUE)
  )
}

derive_base_class <- function(toxval_type) {
  x <- toupper(as.character(toxval_type))
  x <- stringr::str_squish(x)

  dplyr::case_when(
    stringr::str_detect(x, "BMDL|BMCL") ~ "BMDL",
    stringr::str_detect(x, "\\bBMD\\b|\\bBMC\\b") ~ "BMD",

    stringr::str_detect(x, "NOAEL|NOAEC") ~ "NOA",
    stringr::str_detect(x, "LOAEL|LOAEC") ~ "LOA",

    stringr::str_detect(x, "\\bNOEL\\b|\\bNOEC\\b|\\bNEL\\b|\\bNEC\\b") ~ "NE",
    stringr::str_detect(x, "\\bLOEL\\b|\\bLOEC\\b|\\bLEL\\b|\\bLEC\\b") ~ "LE",

    TRUE ~ NA_character_
  )
}

derive_toxval_type_standard <- function(base_class) {
  dplyr::case_when(
    base_class %in% c("NOA", "NOAEL") ~ "NOAEL",
    base_class %in% c("LOA", "LOAEL") ~ "LOAEL",
    base_class %in% c("NE", "NOEL", "NEL") ~ "NOEL",
    base_class %in% c("LE", "LOEL", "LEL") ~ "LOEL",
    base_class == "BMDL" ~ "BMDL",
    base_class == "BMD" ~ "BMD",
    TRUE ~ NA_character_
  )
}

derive_pod_class <- function(toxval_type_standard) {
  dplyr::case_when(
    toxval_type_standard == "BMDL" ~ "BMDL",
    toxval_type_standard == "BMD" ~ "BMD",
    toxval_type_standard %in% c("NOAEL", "LOAEL", "NOEL", "LOEL") ~ "LEVEL",
    TRUE ~ NA_character_
  )
}

derive_effect_type <- function(study_type, toxicological_effect_category) {
  study_type_clean <- study_type %>%
    as.character() %>%
    stringr::str_squish() %>%
    stringr::str_to_lower()

  effect_category_clean <- toxicological_effect_category %>%
    as.character() %>%
    stringr::str_squish() %>%
    stringr::str_to_lower()

  dplyr::case_when(
    stringr::str_detect(study_type_clean, "development|repro") ~
      "reproductive developmental",

    stringr::str_detect(effect_category_clean, "development|repro") ~
      "reproductive developmental",

    TRUE ~ "repeat dose"
  )
}

derive_duration_adjustment <- function(
  study_type,
  duration_class,
  toxicological_effect_category
) {
  effect_type <- derive_effect_type(
    study_type = study_type,
    toxicological_effect_category = toxicological_effect_category
  )

  duration_clean <- duration_class %>%
    as.character() %>%
    stringr::str_squish() %>%
    stringr::str_to_lower()

  dplyr::case_when(
    effect_type == "reproductive developmental" ~ "no adjustment",

    stringr::str_detect(duration_clean, "chronic") ~ "chronic",

    stringr::str_detect(duration_clean, "subchronic|sub-chronic") ~
      "subchronic",

    stringr::str_detect(duration_clean, "short") ~
      "short-term",

    is.na(duration_clean) | duration_clean == "" | duration_clean == "na" ~
      "unspecified",

    TRUE ~ "unspecified"
  )
}

#===============================================================================
# 1. READ HAZARD LANDSCAPE DRSVS
#===============================================================================

drsv <- readRDS(
  here::here("data", "derived", "hazard_landscape_drsvs.rds")
)

if (!"row_id" %in% names(drsv)) {
  drsv <- drsv %>%
    dplyr::mutate(row_id = dplyr::row_number())
}

required_drsv_cols <- c(
  "row_id",
  "dtxsid",
  "toxval_type",
  "toxval_numeric",
  "toxicological_effect",
  "study_type",
  "duration_class",
  "common_name"
)

missing_drsv_cols <- setdiff(required_drsv_cols, names(drsv))

if (length(missing_drsv_cols) > 0) {
  stop(
    "hazard_landscape_drsvs.rds is missing required columns: ",
    paste(missing_drsv_cols, collapse = ", ")
  )
}

#===============================================================================
# 2. TOXICOLOGICAL EFFECT TERM -> EFFECT CATEGORY
#===============================================================================

mapping_file <- here::here(
  "config",
  "toxicological_effect_categories_mapped.csv"
)

if (!file.exists(mapping_file)) {
  stop(
    "Missing mapping file: ",
    mapping_file
  )
}

effect_map_raw <- readr::read_csv(
  mapping_file,
  show_col_types = FALSE
)

required_effect_map_cols <- c(
  "toxicological_effect",
  "toxicological_effect_category"
)

missing_effect_map_cols <- setdiff(
  required_effect_map_cols,
  names(effect_map_raw)
)

if (length(missing_effect_map_cols) > 0) {
  stop(
    "toxicological_effect_categories_mapped.csv is missing required columns: ",
    paste(missing_effect_map_cols, collapse = ", ")
  )
}

effect_map <- effect_map_raw %>%
  dplyr::select(
    toxicological_effect,
    toxicological_effect_category
  ) %>%
  dplyr::mutate(
    toxicological_effect = toxicological_effect %>%
      as.character() %>%
      stringr::str_squish() %>%
      stringr::str_to_lower(),

    toxicological_effect_category = toxicological_effect_category %>%
      as.character() %>%
      stringr::str_squish(),

    toxicological_effect_category =
      dplyr::na_if(toxicological_effect_category, "")
  ) %>%
  dplyr::filter(
    !is.na(toxicological_effect),
    toxicological_effect != ""
  ) %>%
  dplyr::distinct(
    toxicological_effect,
    .keep_all = TRUE
  )

effect_long <- drsv %>%
  dplyr::select(
    row_id,
    toxicological_effect
  ) %>%
  tidyr::separate_rows(
    toxicological_effect,
    sep = "\\|"
  ) %>%
  dplyr::mutate(
    toxicological_effect = toxicological_effect %>%
      as.character() %>%
      stringr::str_squish() %>%
      stringr::str_to_lower()
  ) %>%
  dplyr::filter(
    !is.na(toxicological_effect),
    toxicological_effect != ""
  )

effect_long_mapped <- effect_long %>%
  dplyr::left_join(
    effect_map,
    by = "toxicological_effect"
  )

missing_toxicological_effect_terms <- effect_long_mapped %>%
  dplyr::filter(
    is.na(toxicological_effect_category)
  ) %>%
  dplyr::count(
    toxicological_effect,
    sort = TRUE,
    name = "n_rows"
  )

readr::write_csv(
  missing_toxicological_effect_terms,
  here::here(
    "outputs",
    "effect_mapping",
    "tables",
    "missing_toxicological_effect_terms.csv"
  )
)

effect_category_by_row <- effect_long_mapped %>%
  dplyr::group_by(row_id) %>%
  dplyr::summarise(
    toxicological_effect_category_mapped = paste(
      unique(stats::na.omit(toxicological_effect_category)),
      collapse = "|"
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    toxicological_effect_category_mapped =
      dplyr::na_if(toxicological_effect_category_mapped, "")
  )

hazard_effect_category_mapped <- drsv %>%
  dplyr::left_join(
    effect_category_by_row,
    by = "row_id"
  ) %>%
  dplyr::mutate(
    effect_category_mapping_missing =
      is.na(toxicological_effect_category_mapped)
  )

saveRDS(
  hazard_effect_category_mapped,
  here::here(
    "data",
    "derived",
    "hazard_effect_category_mapped.rds"
  )
)

#===============================================================================
# 3. EFFECT CATEGORY -> CONCEPTUAL MODEL
#===============================================================================

model_map_file <- here::here(
  "config",
  "conceptual_model_map.csv"
)

if (!file.exists(model_map_file)) {
  stop(
    "Missing conceptual model mapping file: ",
    model_map_file
  )
}

model_map <- readr::read_csv(
  model_map_file,
  show_col_types = FALSE
)

required_model_map_cols <- c(
  "standard",
  "model1",
  "model2"
)

missing_model_map_cols <- setdiff(
  required_model_map_cols,
  names(model_map)
)

if (length(missing_model_map_cols) > 0) {
  stop(
    "conceptual_model_map.csv is missing required columns: ",
    paste(missing_model_map_cols, collapse = ", ")
  )
}

model_map_clean <- model_map %>%
  dplyr::select(
    standard,
    model1,
    model2
  ) %>%
  dplyr::mutate(
    standard = standard %>%
      as.character() %>%
      stringr::str_squish() %>%
      stringr::str_to_lower(),

    model1 = model1 %>%
      as.character() %>%
      stringr::str_squish(),

    model2 = model2 %>%
      as.character() %>%
      stringr::str_squish(),

    model1 = dplyr::na_if(model1, ""),
    model2 = dplyr::na_if(model2, "")
  ) %>%
  dplyr::distinct(
    standard,
    .keep_all = TRUE
  ) %>%
  dplyr::rename(
    model1_all = model1,
    model2_all = model2
  )

mapping_base <- hazard_effect_category_mapped %>%
  dplyr::mutate(
    study_type_clean = study_type %>%
      as.character() %>%
      stringr::str_squish() %>%
      stringr::str_to_lower(),

    effect_type_mapped = dplyr::case_when(
      stringr::str_detect(study_type_clean, "development|repro") ~ "repro dev",
      TRUE ~ "repeat dose"
    )
  )

mapping_long <- mapping_base %>%
  dplyr::select(
    row_id,
    study_type,
    effect_type_mapped,
    toxicological_effect_category_mapped
  ) %>%
  tidyr::separate_rows(
    toxicological_effect_category_mapped,
    sep = "\\|"
  ) %>%
  dplyr::mutate(
    toxicological_effect_category_mapped =
      toxicological_effect_category_mapped %>%
      as.character() %>%
      stringr::str_squish() %>%
      stringr::str_to_lower()
  ) %>%
  dplyr::filter(
    !is.na(toxicological_effect_category_mapped),
    toxicological_effect_category_mapped != ""
  )

mapping_long_joined <- mapping_long %>%
  dplyr::left_join(
    model_map_clean,
    by = c("toxicological_effect_category_mapped" = "standard")
  )

unmapped_toxicological_effect_categories <- mapping_long_joined %>%
  dplyr::filter(
    is.na(model1_all),
    is.na(model2_all)
  ) %>%
  dplyr::count(
    toxicological_effect_category_mapped,
    study_type,
    effect_type_mapped,
    sort = TRUE,
    name = "n_rows"
  )

readr::write_csv(
  unmapped_toxicological_effect_categories,
  here::here(
    "outputs",
    "effect_mapping",
    "tables",
    "unmapped_toxicological_effect_categories.csv"
  )
)

mapping_collapsed <- mapping_long_joined %>%
  dplyr::group_by(row_id) %>%
  dplyr::summarise(
    toxicological_effect_category_conceptual = paste(
      unique(stats::na.omit(toxicological_effect_category_mapped)),
      collapse = "|"
    ),

    n_effect_categories_mapped = dplyr::n_distinct(
      toxicological_effect_category_mapped,
      na.rm = TRUE
    ),

    model1_all = paste(
      unique(stats::na.omit(model1_all)),
      collapse = "|"
    ),

    model2_all = paste(
      unique(stats::na.omit(model2_all)),
      collapse = "|"
    ),

    .groups = "drop"
  ) %>%
  dplyr::mutate(
    toxicological_effect_category_conceptual =
      dplyr::na_if(toxicological_effect_category_conceptual, ""),

    model1_all =
      dplyr::na_if(model1_all, ""),

    model2_all =
      dplyr::na_if(model2_all, "")
  )

mapping_row_level <- mapping_base %>%
  dplyr::select(
    row_id,
    effect_type_mapped
  ) %>%
  dplyr::left_join(
    mapping_collapsed,
    by = "row_id"
  ) %>%
  dplyr::mutate(
    multiple_flag = dplyr::case_when(
      n_effect_categories_mapped > 1 ~ "multiple",

      has_term(model1_all, "cont") &
        has_term(model1_all, "det") ~ "multiple",

      has_term(model1_all, "cont") &
        has_term(model1_all, "stoch") ~ "multiple",

      has_term(model1_all, "det") &
        has_term(model1_all, "stoch") ~ "multiple",

      TRUE ~ NA_character_
    )
  )

mapping_final <- mapping_row_level %>%
  dplyr::mutate(
    final_model1 = dplyr::case_when(
      has_term(model1_all, "cont") &
        has_term(model1_all, "det") ~ "continuous",

      has_term(model1_all, "cont") &
        has_term(model1_all, "stoch") ~ "continuous",

      has_term(model1_all, "det") &
        has_term(model1_all, "stoch") ~ "continuous",

      has_term(model1_all, "det") ~ "quantal-deterministic",

      has_term(model1_all, "stoch") ~ "quantal-stochastic",

      has_term(model1_all, "cont") ~ "continuous",

      TRUE ~ NA_character_
    ),

    model2 = dplyr::case_when(
      has_term(model2_all, "det") &
        has_term(model2_all, "stoch") &
        effect_type_mapped == "repeat dose" ~ "quantal-deterministic",

      has_term(model2_all, "det") &
        has_term(model2_all, "stoch") &
        effect_type_mapped == "repro dev" ~ "quantal-stochastic",

      has_term(model2_all, "det") ~ "quantal-deterministic",

      has_term(model2_all, "stoch") ~ "quantal-stochastic",

      TRUE ~ NA_character_
    ),

    final_model2 = dplyr::case_when(
      is.na(model2) ~ NA_character_,

      model2 == "quantal-stochastic" ~ "quantal-stochastic",

      model2 == "quantal-deterministic" ~ "quantal-deterministic",

      multiple_flag == "multiple" &
        effect_type_mapped == "repeat dose" ~ "quantal-deterministic",

      multiple_flag == "multiple" &
        effect_type_mapped == "repro dev" ~ "quantal-stochastic",

      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::select(
    row_id,
    toxicological_effect_category_conceptual,
    n_effect_categories_mapped,
    effect_type_mapped,
    model1_all,
    model2_all,
    multiple_flag,
    final_model1,
    final_model2
  )

hazard_effect_mapped <- hazard_effect_category_mapped %>%
  dplyr::left_join(
    mapping_final,
    by = "row_id"
  ) %>%
  dplyr::mutate(
    conceptual_model_missing = is.na(final_model1)
  )

model_assignment_summary <- hazard_effect_mapped %>%
  dplyr::count(
    final_model1,
    final_model2,
    sort = TRUE,
    name = "n_rows"
  )

mapping_coverage_summary <- tibble::tibble(
  metric = c(
    "n_input_rows",
    "n_output_rows",
    "n_rows_with_effect_category_mapped",
    "n_rows_missing_effect_category_mapped",
    "n_unique_missing_toxicological_effect_terms",
    "n_rows_with_final_model1",
    "n_rows_with_final_model2",
    "n_rows_missing_final_model1",
    "n_unmapped_toxicological_effect_categories"
  ),
  value = c(
    nrow(drsv),
    nrow(hazard_effect_mapped),
    sum(!is.na(hazard_effect_mapped$toxicological_effect_category_mapped)),
    sum(is.na(hazard_effect_mapped$toxicological_effect_category_mapped)),
    nrow(missing_toxicological_effect_terms),
    sum(!is.na(hazard_effect_mapped$final_model1)),
    sum(!is.na(hazard_effect_mapped$final_model2)),
    sum(is.na(hazard_effect_mapped$final_model1)),
    nrow(unmapped_toxicological_effect_categories)
  )
)

saveRDS(
  hazard_effect_mapped,
  here::here(
    "data",
    "derived",
    "hazard_effect_mapped.rds"
  )
)

readr::write_csv(
  model_assignment_summary,
  here::here(
    "outputs",
    "effect_mapping",
    "tables",
    "model_assignment_summary.csv"
  )
)

readr::write_csv(
  mapping_coverage_summary,
  here::here(
    "outputs",
    "effect_mapping",
    "tables",
    "mapping_coverage_summary.csv"
  )
)

#===============================================================================
# 4. ORAL DRSV -> eBMD.HED CONVERSION
#===============================================================================

dat <- hazard_effect_mapped

required_conversion_cols <- c(
  "row_id",
  "dtxsid",
  "toxval_type",
  "toxval_numeric",
  "final_model1",
  "final_model2",
  "study_type",
  "duration_class",
  "common_name"
)

missing_conversion_cols <- setdiff(
  required_conversion_cols,
  names(dat)
)

if (length(missing_conversion_cols) > 0) {
  stop(
    "hazard_effect_mapped is missing required columns for eBMD.HED conversion: ",
    paste(missing_conversion_cols, collapse = ", ")
  )
}

effect_category_candidates <- c(
  "toxicological_effect_category_mapped",
  "toxicological_effect_category_conceptual",
  "toxicological_effect_category"
)

effect_category_col <- effect_category_candidates[
  effect_category_candidates %in% names(dat)
][1]

if (is.na(effect_category_col)) {
  warning(
    "No toxicological effect category column found. ",
    "Duration/type derivation will rely on study_type only."
  )

  dat$toxicological_effect_category_for_conversion <- NA_character_
} else {
  dat <- dat %>%
    dplyr::mutate(
      toxicological_effect_category_for_conversion =
        .data[[effect_category_col]]
    )
}

dat <- dat %>%
  dplyr::mutate(
    selected_value_oral = toxval_numeric
  )

if (!"toxval_numeric_hed" %in% names(dat)) {
  dat$toxval_numeric_hed <- 0
}

dat <- dat %>%
  dplyr::mutate(
    toxval_numeric_hed = dplyr::case_when(
      as.character(toxval_numeric_hed) %in% c(
        "1", "TRUE", "True", "true", "yes", "YES"
      ) ~ 1,
      TRUE ~ 0
    )
  )

dat <- dat %>%
  dplyr::mutate(
    toxval_numeric_hed = dplyr::case_when(
      stringr::str_detect(
        as.character(toxval_type),
        stringr::regex("HED|human equivalent", ignore_case = TRUE)
      ) ~ 1,
      TRUE ~ toxval_numeric_hed
    )
  )

if ("selected_value_basis" %in% names(dat)) {
  dat <- dat %>%
    dplyr::mutate(
      toxval_numeric_hed = dplyr::case_when(
        stringr::str_detect(
          as.character(selected_value_basis),
          stringr::regex("HED|human equivalent", ignore_case = TRUE)
        ) ~ 1,
        TRUE ~ toxval_numeric_hed
      )
    )
}

dat <- dat %>%
  dplyr::mutate(
    final_model1 = final_model1 %>%
      as.character() %>%
      stringr::str_squish() %>%
      dplyr::na_if("") %>%
      dplyr::na_if("-"),

    final_model2 = final_model2 %>%
      as.character() %>%
      stringr::str_squish() %>%
      dplyr::na_if("") %>%
      dplyr::na_if("-"),

    common_name = common_name %>%
      as.character() %>%
      stringr::str_squish(),

    common_name = dplyr::case_when(
      stringr::str_to_lower(common_name) == "rat" ~ "Rat",
      stringr::str_to_lower(common_name) == "mouse" ~ "Mouse",
      stringr::str_to_lower(common_name) == "rabbit" ~ "Rabbit",
      stringr::str_to_lower(common_name) == "dog" ~ "Dog",
      stringr::str_to_lower(common_name) == "human" ~ "Human",
      TRUE ~ common_name
    ),

    base_class = derive_base_class(toxval_type),

    toxval_type_standard =
      derive_toxval_type_standard(base_class),

    pod_class =
      derive_pod_class(toxval_type_standard),

    type = derive_effect_type(
      study_type = study_type,
      toxicological_effect_category =
        toxicological_effect_category_for_conversion
    ),

    duration_adjustment = derive_duration_adjustment(
      study_type = study_type,
      duration_class = duration_class,
      toxicological_effect_category =
        toxicological_effect_category_for_conversion
    )
  )

#---------------------------
# conversion factors
#---------------------------
dat <- dat %>%
  dplyr::mutate(
    CF1 = dplyr::case_when(
      duration_adjustment == "subchronic" ~ 2,
      duration_adjustment == "short-term" ~ 5,
      duration_adjustment == "chronic" ~ 1,
      duration_adjustment == "no adjustment" ~ 1,
      duration_adjustment == "unspecified" ~ 1,
      TRUE ~ 1
    ),

    CF2 = dplyr::case_when(
      is.na(final_model1) ~ NA_real_,
      toxval_numeric_hed == 1 ~ 1,
      common_name %in% names(CF2_lookup) ~
        lookup_named_numeric(common_name, CF2_lookup),
      TRUE ~ NA_real_
    ),

    CF3 = dplyr::case_when(
      toxval_type_standard %in% c("LOAEL", "LOEL") ~ 3,
      TRUE ~ 1
    ),

    CF41 = dplyr::case_when(
      is.na(final_model1) ~ NA_real_,

      pod_class == "BMD" ~ 1,

      pod_class == "BMDL" &
        final_model1 == "continuous" ~ 1 / 3,

      pod_class == "BMDL" &
        final_model1 == "quantal-deterministic" ~ 1 / 9,

      pod_class == "BMDL" &
        final_model1 == "quantal-stochastic" ~ 1 / 3,

      pod_class == "LEVEL" &
        final_model1 == "continuous" ~ 1 / 3,

      pod_class == "LEVEL" &
        final_model1 == "quantal-deterministic" ~ 2 / 9,

      pod_class == "LEVEL" &
        final_model1 == "quantal-stochastic" ~ 2 / 3,

      TRUE ~ 1
    ),

    CF42 = dplyr::case_when(
      is.na(final_model1) | is.na(final_model2) ~ NA_real_,

      pod_class == "BMD" ~ 1,

      pod_class == "BMDL" &
        final_model2 == "continuous" ~ 1 / 3,

      pod_class == "BMDL" &
        final_model2 == "quantal-deterministic" ~ 1 / 9,

      pod_class == "BMDL" &
        final_model2 == "quantal-stochastic" ~ 1 / 3,

      pod_class == "LEVEL" &
        final_model2 == "continuous" ~ 1 / 3,

      pod_class == "LEVEL" &
        final_model2 == "quantal-deterministic" ~ 2 / 9,

      pod_class == "LEVEL" &
        final_model2 == "quantal-stochastic" ~ 2 / 3,

      TRUE ~ 1
    ),

    CF5 = 1
  )

#---------------------------
# derive eBMD.HED 
#---------------------------
dat <- dat  %>%
  dplyr::mutate(
    denom1 = CF1 * CF2 * CF3 * CF41 * CF5,

    denom2 = dplyr::case_when(
      is.na(final_model2) ~ NA_real_,
      TRUE ~ CF1 * CF2 * CF3 * CF42 * CF5
    ),

    eBMD.HED1 = selected_value_oral / denom1,

    eBMD.HED2 = dplyr::case_when(
      is.na(final_model2) ~ NA_real_,
      TRUE ~ selected_value_oral / denom2
    ),

    eBMD.HED = dplyr::case_when(
      !is.na(final_model2) &
        !is.na(eBMD.HED1) &
        !is.na(eBMD.HED2) &
        eBMD.HED1 > 0 &
        eBMD.HED2 > 0 ~
        10^(0.5 * (log10(eBMD.HED1) + log10(eBMD.HED2))),

      TRUE ~ eBMD.HED1
    ),

    eBMD.HED2 = dplyr::case_when(
      is.na(final_model2) ~ NA_real_,
      TRUE ~ eBMD.HED2
    )
  )


#
hed_p_by_chem_ebmdhed <- dat %>%
  filter(common_name != "Human") %>%
  filter(!is.na(cichem_id), !is.na(eBMD.HED)) %>%
  group_by(cichem_id) %>%
  summarise(
    n_hed_values = n(),

    hed_p05 = as.numeric(
      quantile(
        eBMD.HED,
        probs = 0.05,
        na.rm = TRUE,
        type = 7
      )
    ),

        hed_p18 = as.numeric(
      quantile(
        eBMD.HED,
        probs = 0.1877,
        na.rm = TRUE,
        type = 7
      )
    ),

    hed_p25 = as.numeric(
      quantile(
        eBMD.HED,
        probs = 0.25,
        na.rm = TRUE,
        type = 7
      )
    ),

    .groups = "drop"
  )

saveRDS(
  hed_p_by_chem_ebmdhed, here::here(
  "data",
  "derived",
  "hed_p_by_chem_ebmdhed.rds")
)


#===============================================================================
# 5. DIAGNOSTICS AND OUTPUTS
#===============================================================================

diagnostic_summary <- dat %>%
  dplyr::summarise(
    n_rows = dplyr::n(),
    n_missing_selected_value_oral = sum(is.na(selected_value_oral)),
    n_nonpositive_selected_value_oral = sum(
      !is.na(selected_value_oral) & selected_value_oral <= 0
    ),
    n_missing_toxicological_effect_category_mapped =
      sum(is.na(toxicological_effect_category_mapped)),
    n_missing_final_model1 = sum(is.na(final_model1)),
    n_missing_final_model2 = sum(is.na(final_model2)),
    n_missing_toxval_type = sum(is.na(toxval_type)),
    n_missing_base_class = sum(is.na(base_class)),
    n_missing_toxval_type_standard = sum(is.na(toxval_type_standard)),
    n_missing_pod_class = sum(is.na(pod_class)),
    n_missing_common_name = sum(is.na(common_name)),
    n_missing_CF2 = sum(is.na(CF2)),
    n_missing_CF41 = sum(is.na(CF41)),
    n_missing_denom1 = sum(is.na(denom1)),
    n_missing_eBMD_HED1 = sum(is.na(eBMD.HED1)),
    n_missing_eBMD_HED = sum(is.na(eBMD.HED)),
    n_nonpositive_eBMD_HED1 = sum(!is.na(eBMD.HED1) & eBMD.HED1 <= 0),
    n_nonpositive_eBMD_HED2 = sum(!is.na(eBMD.HED2) & eBMD.HED2 <= 0)
  )

print(diagnostic_summary)

dat_unfiltered <- dat

if (drop_missing_ebmd_hed) {
  dat <- dat %>%
    dplyr::filter(!is.na(eBMD.HED))
}

conversion_summary <- tibble::tibble(
  metric = c(
    "n_input_rows",
    "n_output_rows",
    "n_removed_missing_eBMD_HED",
    "n_unique_missing_toxicological_effect_terms",
    "n_unmapped_toxicological_effect_categories",
    "n_missing_final_model1",
    "n_with_final_model2",
    "n_bmdl_rows",
    "n_bmd_rows",
    "n_level_rows",
    "n_human_rows",
    "n_already_hed_rows",
    "n_missing_CF2",
    "n_missing_base_class"
  ),
  value = c(
    as.character(nrow(dat_unfiltered)),
    as.character(nrow(dat)),
    as.character(nrow(dat_unfiltered) - nrow(dat)),
    as.character(nrow(missing_toxicological_effect_terms)),
    as.character(nrow(unmapped_toxicological_effect_categories)),
    as.character(sum(is.na(dat_unfiltered$final_model1))),
    as.character(sum(!is.na(dat_unfiltered$final_model2))),
    as.character(sum(dat_unfiltered$pod_class == "BMDL", na.rm = TRUE)),
    as.character(sum(dat_unfiltered$pod_class == "BMD", na.rm = TRUE)),
    as.character(sum(dat_unfiltered$pod_class == "LEVEL", na.rm = TRUE)),
    as.character(sum(dat_unfiltered$common_name == "Human", na.rm = TRUE)),
    as.character(sum(dat_unfiltered$toxval_numeric_hed == 1, na.rm = TRUE)),
    as.character(sum(is.na(dat_unfiltered$CF2))),
    as.character(sum(is.na(dat_unfiltered$base_class)))
  )
)

cf_summary <- dat_unfiltered %>%
  dplyr::count(
    duration_adjustment,
    type,
    common_name,
    toxval_numeric_hed,
    toxval_type_standard,
    pod_class,
    CF1,
    CF2,
    CF3,
    CF5,
    sort = TRUE,
    name = "n_rows"
  )

missing_cf2_summary <- dat_unfiltered %>%
  dplyr::filter(is.na(CF2)) %>%
  dplyr::count(
    common_name,
    toxval_numeric_hed,
    final_model1,
    sort = TRUE,
    name = "n_rows"
  )

missing_model_summary <- dat_unfiltered %>%
  dplyr::filter(is.na(final_model1)) %>%
  dplyr::count(
    toxicological_effect_category_for_conversion,
    toxicological_effect_category_mapped,
    toxval_type,
    study_type,
    sort = TRUE,
    name = "n_rows"
  )

missing_base_class_summary <- dat_unfiltered %>%
  dplyr::filter(is.na(base_class)) %>%
  dplyr::count(
    toxval_type,
    sort = TRUE,
    name = "n_rows"
  )

#---------------------------
# save final outputs
#---------------------------
saveRDS(
  dat,
  here::here(
    "data",
    "derived",
    "hazard_eBMD_HED.rds"
  )
)

saveRDS(
  dat_unfiltered,
  here::here(
    "data",
    "derived",
    "hazard_eBMD_HED_unfiltered.rds"
  )
)

readr::write_csv(
  conversion_summary,
  here::here(
    "outputs",
    "ebmdhed_conversion",
    "tables",
    "ebmdhed_conversion_summary.csv"
  )
)

readr::write_csv(
  cf_summary,
  here::here(
    "outputs",
    "ebmdhed_conversion",
    "tables",
    "ebmdhed_cf_summary.csv"
  )
)

readr::write_csv(
  missing_cf2_summary,
  here::here(
    "outputs",
    "ebmdhed_conversion",
    "tables",
    "ebmdhed_missing_cf2_summary.csv"
  )
)

readr::write_csv(
  missing_model_summary,
  here::here(
    "outputs",
    "ebmdhed_conversion",
    "tables",
    "ebmdhed_missing_model_summary.csv"
  )
)

readr::write_csv(
  missing_base_class_summary,
  here::here(
    "outputs",
    "ebmdhed_conversion",
    "tables",
    "ebmdhed_missing_base_class_summary.csv"
  )
)

message(
  "Hazard landscape oral DRSV to eBMD.HED workflow complete: ",
  "n_input_rows = ", nrow(dat_unfiltered),
  "; n_output_rows = ", nrow(dat),
  "; n_removed_missing_eBMD_HED = ", nrow(dat_unfiltered) - nrow(dat),
  "; n_unique_missing_toxicological_effect_terms = ",
  nrow(missing_toxicological_effect_terms),
  "; n_unmapped_toxicological_effect_categories = ",
  nrow(unmapped_toxicological_effect_categories),
  "; n_with_final_model2 = ", sum(!is.na(dat_unfiltered$final_model2)),
  "; n_bmdl_rows = ", sum(dat_unfiltered$pod_class == "BMDL", na.rm = TRUE),
  "; n_bmd_rows = ", sum(dat_unfiltered$pod_class == "BMD", na.rm = TRUE),
  "; n_level_rows = ", sum(dat_unfiltered$pod_class == "LEVEL", na.rm = TRUE),
  "; n_missing_CF2 = ", sum(is.na(dat_unfiltered$CF2)),
  "; n_missing_base_class = ", sum(is.na(dat_unfiltered$base_class)),
  "; n_missing_final_model1 = ", sum(is.na(dat_unfiltered$final_model1))
)
