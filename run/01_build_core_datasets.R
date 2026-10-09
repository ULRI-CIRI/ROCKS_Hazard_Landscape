# ===========================================================
# 01_build_core_datasets.R
#
# Purpose:
#   Build core datasets required for the Hazard Landscape
#   analysis workflow from source databases and raw input files.
#
# Identifier policy:
#   The ROCKS inventory is defined at the cichem_id level.
#   Therefore, cichem_id is the primary key for all inventory-level
#   joins, coverage summaries, manuscript statistics, and figures.
#
#   Some evidence sources, including Broad Profiling and InvitroDB,
#   report chemicals by dtxsid only. These sources are mapped to
#   cichem_id during this build step using id_crosswalk.rds, which
#   is created directly from the built ROCKS inventory.
#
#   dtxsid is retained as a source/reporting identifier, but downstream
#   analyses use cichem_id as the primary join and grouping key.

# ===========================================================

library(DBI)
library(RPostgres)
library(dplyr)
library(stringr)
library(janitor)
library(here)

# =====================================================
# SETUP
# =====================================================

dir.create(
  here::here("data", "derived"),
  recursive = TRUE,
  showWarnings = FALSE
)

normalize_id_types <- function(df) {
  df %>%
    mutate(
      across(
        any_of(c(
          "cichem_id",
          "dtxsid",
          "cicomp_id",
          "inchi",
          "inchikey",
          "smiles"
        )),
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

  existing_dtxsid_col <- intersect(dtxsid_candidates, names(df))

  if (!"dtxsid" %in% names(df) && length(existing_dtxsid_col) > 0) {
    names(df)[names(df) == existing_dtxsid_col[1]] <- "dtxsid"
  }

  df %>%
    normalize_id_types()
}

add_cichem_id_from_dtxsid <- function(df, id_crosswalk) {
  df <- df %>%
    standardize_dtxsid_column()

  if (!"dtxsid" %in% names(df)) {
    warning("No dtxsid-like column found. cichem_id was not added.")
    return(df)
  }

  dtxsid_map <- id_crosswalk %>%
    select(cichem_id, dtxsid) %>%
    filter(!is.na(dtxsid), dtxsid != "") %>%
    distinct()

  if ("cichem_id" %in% names(df)) {
    df <- df %>%
      left_join(
        dtxsid_map %>%
          rename(cichem_id_from_dtxsid = cichem_id),
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
      left_join(dtxsid_map, by = "dtxsid")
  }

  df %>%
    normalize_id_types()
}

coalesce_suffix_column <- function(df, base_col, suffix_col) {
  if (suffix_col %in% names(df)) {
    if (!base_col %in% names(df)) {
      df[[base_col]] <- NA_character_
    }

    df[[base_col]] <- coalesce(
      as.character(df[[base_col]]),
      as.character(df[[suffix_col]])
    )

    df[[suffix_col]] <- NULL
  }

  df
}

# =====================================================
# TAM DATASETS
# =====================================================

message("\nBuilding TAM datasets...")

con_tam <- dbConnect(
  RPostgres::Postgres(),
  dbname   = Sys.getenv("DB_NAME"),
  host     = Sys.getenv("DB_HOST"),
  port     = Sys.getenv("DB_PORT"),
  user     = Sys.getenv("DB_USER"),
  password = Sys.getenv("DB_PASSWORD")
)

on.exit({
  try(DBI::dbDisconnect(con_tam), silent = TRUE)
}, add = TRUE)

tamdb.schema <- "invivodb_v01"

base_query <- paste0(
  "SELECT DISTINCT * ",
  " FROM ", tamdb.schema, ".tamdb b ",
  "LEFT JOIN ", tamdb.schema, ".source_chemical a on a.chemical_id=b.chemical_id ",
  "LEFT JOIN ", tamdb.schema, ".species d on b.species_id=d.species_id ",
  "LEFT JOIN ", tamdb.schema, ".type_dictionary e on b.toxval_type=e.toxval_type ",
  "LEFT JOIN ", tamdb.schema, ".record_source f on b.toxval_id=f.toxval_id ",
  "WHERE f.record_source_level not in ('origin')"
) %>%
  gsub("source_chemical..", "a.", ., fixed = TRUE) %>%
  gsub("tamdb..", "b.", ., fixed = TRUE) %>%
  gsub("species..", "d.", ., fixed = TRUE) %>%
  gsub("tamdb_type_dictionary..", "e.", ., fixed = TRUE) %>%
  gsub("record_source..", "f.", ., fixed = TRUE)

tam_all <- dbGetQuery(con_tam, base_query)

tam <- tam_all %>%
  clean_names() %>%
  filter(
    toxval_type_supercategory %in% c(
      "Dose Response Summary Value",
      "Mortality Response Summary Value"
    )
  ) %>%
  distinct()

saveRDS(
  tam,
  here::here("data", "derived", "tam.rds")
)
drsvs<- tam %>% # oral only
  filter(toxval_type_supercategory=="Dose Response Summary Value")%>%
  filter(exposure_route=="oral")%>%
  filter(toxval_units=="mg/kg-day")

saveRDS(
  drsvs,
  here::here("data", "derived", "hazard_landscape_drsvs.rds")
)

tam_assess <- tam_all %>%
  clean_names() %>%
  filter(
    !toxval_type_supercategory %in% c(
      "Dose Response Summary Value",
      "Mortality Response Summary Value"
    )
  ) %>%
  distinct()

tam_assess <- tam_assess %>%
  mutate(
    source_group = case_when(
      str_detect(source, "IRIS") ~ "US EPA IRIS",
      str_detect(source, "PPRTV")~ "US EPA PPRTV",
      str_detect(source, "HEAST") ~ "US EPA HEAST",
      str_detect(source, "HHTV") ~ "US EPA HHTV",
      str_detect(source, "ATSDR") ~ "US ATSDR",
      str_detect(source, "DCAP") ~ "DCAP",
      str_detect(source, "Health Canada") ~ "Health Canada",
      str_detect(source, "JECFA") ~ "WHO JECFA",
      str_detect(source, "EFSA") ~ "EFSA",
      str_detect(source, "TCEQ") ~ "TCEQ",
      str_detect(source, "OEHHA") ~ "Cal OEHHA",
      TRUE ~ NA_character_
    )
  )

priority_sources <- c(
  "US EPA IRIS",
  "US EPA PPRTV",
  "US EPA HEAST",
  "US EPA HHTV",
  "US ATSDR",
  "Health Canada",
  "WHO JECFA",
  "EFSA",
  "TCEQ",
  "Cal OEHHA",
  "DCAP"
)

tam_assess_filt <- tam_assess %>%
  filter(source_group %in% priority_sources) %>%
  filter(!toxval_type %in% c("ESL","Odor threshold"))

saveRDS(
  tam_assess_filt,
  here::here("data", "derived", "tam_assess.rds")
)

# =====================================================
# ROCKS INVENTORY INCLUDING FINGERPRINTS
# =====================================================

message("\nBuilding ROCKS inventory...")

connect_pipelines <- function() {
  DBI::dbConnect(
    RPostgres::Postgres(),
    host     = Sys.getenv("DB_HOST"),
    port     = Sys.getenv("DB_PORT"),
    dbname   = "pipelines",
    user     = Sys.getenv("USEUSER"),
    password = Sys.getenv("USEPASSWORD"),
    sslmode  = Sys.getenv("PGSSLMODE", unset = "require")
  )
}

pipe_con <- connect_pipelines()

on.exit({
  try(DBI::dbDisconnect(pipe_con), silent = TRUE)
}, add = TRUE)

chemical_ids <- dbGetQuery(
  pipe_con,
  "
  SELECT
    cichem,
    dtxsid,
    smiles,
    inchi,
    inchikey,
    casrn
  FROM public.chemical_identifier
  "
) %>%
  clean_names() %>%
  rename(cichem_id = cichem) %>%
  normalize_id_types() %>%
  distinct()

rocks_sector <- read.csv(
  here::here("data", "raw", "rocks_sector_table_2026-08-27.csv")
) %>%
  clean_names() %>%
  normalize_id_types()

if ("cichem" %in% names(rocks_sector) && !"cichem_id" %in% names(rocks_sector)) {
  rocks_sector <- rocks_sector %>%
    rename(cichem_id = cichem)
}

chemical_ids_for_join <- chemical_ids %>%
  rename(
    id_dtxsid   = dtxsid,
    id_smiles   = smiles,
    id_inchi    = inchi,
    id_inchikey = inchikey,
    id_casrn = casrn
  )

rocks_ids <- rocks_sector %>%
  left_join(
    chemical_ids_for_join,
    by = "cichem_id"
  )

for (id_col in c("dtxsid", "smiles", "inchi", "inchikey", "casrn")) {
  if (!id_col %in% names(rocks_ids)) {
    rocks_ids[[id_col]] <- NA_character_
  }
}

rocks_ids <- rocks_ids %>%
  mutate(
    dtxsid   = coalesce(as.character(dtxsid),   as.character(id_dtxsid)),
    smiles   = coalesce(as.character(smiles),   as.character(id_smiles)),
    inchi    = coalesce(as.character(inchi),    as.character(id_inchi)),
    inchikey = coalesce(as.character(inchikey), as.character(id_inchikey)),
    casrn = coalesce(as.character(casrn), as.character(id_casrn))

  ) %>%
  select(-starts_with("id_")) %>%
  normalize_id_types()

# read complete fingerprint file
fp_raw <- read.csv(
  here::here("data", "raw", "chemical_fingerprint_densmap_8.18.2026.csv")
) %>%
  clean_names() %>%
  normalize_id_types()

# reduced version used for rocks
fp <- fp_raw %>%
  select(
    any_of(
      c(
        "cichem_id",
        "cicomp_id",
        "inchi",
        "inchikey",
        "smiles",
        "densmap_1",
        "densmap_2"
      )
    )
  ) %>%
  distinct()

if ("cichem_id" %in% names(fp) && "cichem_id" %in% names(rocks_ids)) {
  rocks <- rocks_ids %>%
    left_join(fp, by = "cichem_id", suffix = c("", "_fp"))
} else if ("cicomp_id" %in% names(fp) && "cicomp_id" %in% names(rocks_ids)) {
  rocks <- rocks_ids %>%
    left_join(fp, by = "cicomp_id", suffix = c("", "_fp"))
} else if ("inchi" %in% names(fp) && "inchi" %in% names(rocks_ids)) {
  rocks <- rocks_ids %>%
    left_join(fp, by = "inchi", suffix = c("", "_fp"))
} else if ("inchikey" %in% names(fp) && "inchikey" %in% names(rocks_ids)) {
  rocks <- rocks_ids %>%
    left_join(fp, by = "inchikey", suffix = c("", "_fp"))
} else {
  warning("No clear fingerprint join key found. rocks.rds will be saved without fingerprint coordinates.")
  rocks <- rocks_ids
}

rocks <- rocks %>%
  coalesce_suffix_column("smiles", "smiles_fp") %>%
  coalesce_suffix_column("inchi", "inchi_fp") %>%
  coalesce_suffix_column("inchikey", "inchikey_fp") %>%
    coalesce_suffix_column("casrn", "casrn_fp") %>%
  normalize_id_types() %>%
  distinct()

saveRDS(
  rocks,
  here::here("data", "derived", "rocks.rds")
)

# =====================================================
# ID CROSSWALK BUILT FROM ROCKS
# =====================================================

message("\nBuilding id_crosswalk from rocks...")

id_crosswalk <- rocks %>%
  select(
    any_of(c(
      "cichem_id",
      "dtxsid",
      "cicomp_id",
      "smiles",
      "inchi",
      "inchikey"
    ))
  ) %>%
  normalize_id_types() %>%
  distinct()

saveRDS(
  id_crosswalk,
  here::here("data", "derived", "id_crosswalk.rds")
)

fp_complete <- fp_raw %>%
  left_join(
    id_crosswalk,
    by = c("cichem" = "cichem_id")
  ) %>%
  normalize_id_types() %>%
  distinct()

saveRDS(
  fp_complete,
  here::here("data", "derived", "fp_complete.rds")
)

# ===========================================================
# Get classyfire class for each chemical
#============================================================

chem_class <- dbGetQuery(
  pipe_con,
  "
  SELECT
      ci.cichem,
      cf.kingdom,
      cf.superclass,
      cf.class,
      cf.subclass
  FROM chemical_identifier ci
  LEFT OUTER JOIN descriptors_classyfire cf
    ON ci.chemical_id = cf.chemical_id
  "
) %>%
  janitor::clean_names() %>%
  rename(
    cichem_id = cichem
  ) %>%
  distinct()

chem_class <- rocks %>%
  distinct(cichem_id) %>%
  left_join(
    chem_class,
    by = "cichem_id"
  )

saveRDS(
  chem_class,
  here::here(
    "data",
    "derived",
    "chem_class.rds"
  )
)

# =====================================================
# BROAD PROFILING
# =====================================================

message("\nBuilding Broad Profiling dataset...")

bp <- read.csv(
  here::here("data", "raw", "httr_htpp_screen_table.csv")
) %>%
  standardize_dtxsid_column() %>%
  add_cichem_id_from_dtxsid(id_crosswalk)

saveRDS(
  bp,
  here::here("data", "derived", "bp.rds")
)

# =====================================================
# INVITRODB
# =====================================================

message("\nBuilding InvitroDB/ToxCast datasets...")

load(
  here::here("data", "raw", "rocks_toxcast_summary_tables_05252026.RData")
)

b_exec_sum <- tables_to_save[["bioactivity_exec_summ_final"]] %>%
  add_cichem_id_from_dtxsid(id_crosswalk)

b_assay_sum <- tables_to_save[["assay_list_summary"]] %>%
  add_cichem_id_from_dtxsid(id_crosswalk)

b_assay_data_sum <- tables_to_save[["assay_data_summary"]] %>%
  add_cichem_id_from_dtxsid(id_crosswalk)

b_mc7_aed <- tables_to_save[["mc7_aed_summary"]] %>%
  add_cichem_id_from_dtxsid(id_crosswalk)

b_mc7_assay <- tables_to_save[["mc7_assay_summary"]] %>%
  add_cichem_id_from_dtxsid(id_crosswalk)

b_mc7_aed_exec <- b_mc7_aed %>%
  left_join(
    b_exec_sum,
    by = "dtxsid",
    suffix = c("_mc7", "_exec")
  )

saveRDS(
  b_mc7_aed_exec,
  here::here("data", "derived", "b_mc7_aed_exec.rds")
)

saveRDS(
  b_exec_sum,
  here::here("data", "derived", "b_exec_sum.rds")
)

saveRDS(
  b_assay_sum,
  here::here("data", "derived", "b_assay_sum.rds")
)

saveRDS(
  b_assay_data_sum,
  here::here("data", "derived", "b_assay_data_sum.rds")
)

saveRDS(
  b_mc7_aed,
  here::here("data", "derived", "b_mc7_aed.rds")
)

saveRDS(
  b_mc7_assay,
  here::here("data", "derived", "b_mc7_assay.rds")
)

# =====================================================
# COMPLETION MESSAGE
# =====================================================

message("\n01_build_core_datasets.R completed successfully")
