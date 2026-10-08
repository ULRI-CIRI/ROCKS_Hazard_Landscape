# ===========================================================
# 02_build_potency_objects.R
#
# Purpose:
#   Build chemical-level potency objects from TAM,
#   ToxCast AEDs, and human dose-response data.
# ===========================================================

library(here)
library(dplyr)
library(stringr)

# ===========================================================
# 0) Paths and input checks
# ===========================================================

derived_dir <- here("data", "derived")

required_files <- c(
  here("data", "derived", "tam.rds"),
  here("data", "derived", "rocks.rds"),
  here("data", "derived", "b_mc7_aed.rds")
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing required input files:\n",
    paste(missing_files, collapse = "\n")
  )
}

# ===========================================================
# 1) Load input objects
# ===========================================================

tam <- readRDS(here("data", "derived", "tam.rds"))

rocks <- readRDS(here("data", "derived", "rocks.rds"))

b_mc7_aed <- readRDS(
  here("data", "derived", "b_mc7_aed.rds")
)

# ===========================================================
# 2) EPA BW3/4 conversion function
#    HED = AnimalDose * (BW_animal/BW_human)^0.25
# ===========================================================

hed_bw34 <- function(
    dose_mgkgd,
    bw_animal_kg,
    bw_human_kg = 80
) {

  dose_mgkgd * (bw_animal_kg / bw_human_kg)^0.25

}

# ===========================================================
# 3) Settings
# ===========================================================

BW_HUMAN <- 80

bw_defaults <- c(
  "mouse" = 0.03,
  "rat" = 0.25,
  "hamster" = 0.10,
  "guinea pig" = 0.40,
  "rabbit" = 2.5,
  "dog" = 10.0,
  "cat" = 2.5
)

#===========================================================
# 4) Add HED column using species default body weights
# ===========================================================

tam_hed <- tam %>%
  filter(!toxval_type %in% c("UCDL (10)", "FEL")) %>%
  filter(study_duration_class != "acute") %>%
  filter(
    toxval_type_supercategory ==
      "Dose Response Summary Value"
  ) %>%
  filter(toxval_units == "mg/kg-day") %>%
  filter(exposure_route == "oral") %>%
mutate(
  already_hed =
    grepl("HED", toxval_type, ignore.case = TRUE),

  species_lc =
    tolower(trimws(common_name)),

  bw_animal_kg =
    unname(bw_defaults[species_lc]),

  hed = case_when(
    already_hed ~ toxval_numeric,
    is.na(toxval_numeric) ~ NA_real_,
    is.na(bw_animal_kg) ~ NA_real_,
    TRUE ~ hed_bw34(
      toxval_numeric,
      bw_animal_kg,
      BW_HUMAN
    )
  )
)

# ===========================================================
# 5) Quick QA summaries
# ===========================================================

cat(
  "\nRows already marked HED (kept as-is): ",
  sum(tam_hed$already_hed, na.rm = TRUE),
  "\n"
)

cat(
  "Rows converted to HED (computed): ",
  sum(
    !tam_hed$already_hed &
      !is.na(tam_hed$hed),
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Rows missing BW (could not convert): ",
  sum(
    !tam_hed$already_hed &
      is.na(tam_hed$bw_animal_kg),
    na.rm = TRUE
  ),
  "\n\n"
)

missing_species <- tam_hed %>%
  filter(
    !already_hed,
    is.na(bw_animal_kg)
  ) %>%
  count(species_lc, sort = TRUE)

if (nrow(missing_species) > 0) {

  cat(
    "Species with missing BW defaults ",
    "(edit bw_defaults or provide BW column):\n"
  )

  print(missing_species)

}

# ===========================================================
# 6) Chemical-level HED p05 + p25
# ===========================================================

hed_p_by_chem <- tam_hed %>%
  filter(common_name != "Human") %>%
  filter(!is.na(cichem_id), !is.na(hed)) %>%
  group_by(cichem_id) %>%
  summarise(
    n_hed_values = n(),

    hed_p05 = as.numeric(
      quantile(
        hed,
        probs = 0.05,
        na.rm = TRUE,
        type = 7
      )
    ),

    hed_p25 = as.numeric(
      quantile(
        hed,
        probs = 0.25,
        na.rm = TRUE,
        type = 7
      )
    ),

    .groups = "drop"
  )

print(
  hed_p_by_chem %>%
    arrange(hed_p05) %>%
    head(10)
)


# ===========================================================
# 7) Chemical-level Human p05 + p25
# ===========================================================

human_p_by_chem <- tam_hed %>%
  filter(common_name == "Human") %>%
  filter(!is.na(cichem_id), !is.na(toxval_numeric)) %>%
  group_by(cichem_id) %>%
  summarise(
    n_human_values = n(),

    human_p05 = as.numeric(
      quantile(
        toxval_numeric,
        probs = 0.05,
        na.rm = TRUE,
        type = 7
      )
    ),

    human_p25 = as.numeric(
      quantile(
        toxval_numeric,
        probs = 0.25,
        na.rm = TRUE,
        type = 7
      )
    ),

    .groups = "drop"
  )

# ===========================================================
# 8) Build combined chemical-level table
# ===========================================================

hed_aed_hum <- rocks %>%
  distinct(cichem_id) %>%
  left_join(
    hed_p_by_chem,
    by = "cichem_id"
  ) %>%
  left_join(
    b_mc7_aed,
    by = "cichem_id"
  ) %>%
  left_join(
    human_p_by_chem,
    by = "cichem_id"
  )

# ===========================================================
# 9) Remove temporary QA columns
# ===========================================================

tam_hed <- tam_hed %>%
  select(
    -already_hed,
    -species_lc
  )

# ===========================================================
# 10) Save outputs
# ===========================================================

saveRDS(
  tam_hed,
  file.path(derived_dir, "tam_hed.rds")
)

saveRDS(
  hed_p_by_chem,
  file.path(derived_dir, "hed_p_by_chem.rds")
)

saveRDS(
  human_p_by_chem,
  file.path(derived_dir, "human_p_by_chem.rds")
)

saveRDS(
  hed_aed_hum,
  file.path(derived_dir, "hed_aed_hum.rds")
)

# ===========================================================
# 11) Completion summary
# ===========================================================

message("\n02_build_potency_objects.R completed successfully")

message("Rows in tam_hed: ", nrow(tam_hed))

message(
  "Chemicals in hed_aed_hum: ",
  nrow(hed_aed_hum)
)

message("\nSaved:")

message(
  " - ",
  file.path(derived_dir, "tam_hed.rds")
)

message(
  " - ",
  file.path(derived_dir, "hed_p_by_chem.rds")
)

message(
  " - ",
  file.path(derived_dir, "human_p_by_chem.rds")
)

message(
  " - ",
  file.path(derived_dir, "hed_aed_hum.rds")
)
