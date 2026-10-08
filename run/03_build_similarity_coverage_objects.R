# ===========================================================
# 03_build_similarity_coverage_objects.R
#
# Purpose:
#   Build structural analogue coverage objects for ROCKS chemicals
#   lacking selected toxicity or in vitro data, using an expanded
#   source universe of chemicals with data in TAM DRSV

# How many data-poor ROCKS chemicals have evidence-bearing analogues above fixed Tanimoto thresholds?
#
# Key design:
#   Targets are restricted to ROCKS chemicals without data.
#   Sources are not restricted to ROCKS and include any fingerprinted
#   chemical with data in TAM DRSV
# Notes:
#
#   Target chemicals are ROCKS chemicals lacking both selected data streams.
#
#   Structural similarity is calculated using binary Morgan fingerprint
#   Tanimoto/Jaccard similarity.
#
#   This is a build script only. Plotting in
#   a separate figure script.
# ===========================================================


# ===========================================================
# 0) Libraries
# ===========================================================

library(dplyr)
library(here)
library(Matrix)
library(janitor)
library(DBI)
library(RPostgres)

# ===========================================================
# 1) Paths and input checks
# ===========================================================

derived_dir <- here::here("data", "derived")

dir.create(
  derived_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

### PULL IN MORGAN FINGERPRINTS FOR CHEMICALS NOT ON ROCKS ###
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

fingerprints <- dbGetQuery(
  pipe_con,
  "
  SELECT
    ci.cichem,
    ci.dtxsid,
    ci.smiles,
    ci.inchi,
    ci.inchikey,
    ci.casrn,
    fp.*
  FROM public.descriptors_morgan_fingerprint fp
  LEFT JOIN public.chemical_identifier ci
    ON fp.chemical_id = ci.chemical_id
  "
) %>%
  clean_names() %>%
  rename(cichem_id = cichem) %>%
  normalize_id_types()


# ===========================================================
# 2) Load input objects
# ===========================================================

rocks<- readRDS(
  here::here("data", "derived", "rocks.rds")
)

tam_drsvs<- readRDS(
  here::here("data", "derived", "tam.rds"))%>%
  filter(toxval_type_supercategory=="Dose Response Summary Value")%>%
  filter(study_duration_class!="acute")

b_exec_sum <- readRDS(
  here("data", "derived", "b_exec_sum.rds")
)


# ===========================================================
# 3) Settings
# ===========================================================

primary_threshold <- 0.7

thresholds <- c(0.6, 0.7, 0.8)

chunk_size <- 500


# ===========================================================
# 5) Build broad data availability flag
# ===========================================================

tam_data_ids <- tam_drsvs %>%
  clean_names() %>%
  normalize_id_types() %>%
  filter(!is.na(cichem_id)) %>%
  distinct(cichem_id) %>%
  mutate(has_tam_drsv = TRUE)

invitro_data_ids <- b_exec_sum %>%
  clean_names() %>%
  normalize_id_types() %>%
  filter(!is.na(cichem_id)) %>%
  distinct(cichem_id) %>%
  mutate(has_invitro = TRUE)

#data_available_flag <- full_join(
#  tam_data_ids,
#  invitro_data_ids,
#  by = "cichem_id"
#) %>%
#  mutate(
#    has_tam_drsv = coalesce(has_tam_drsv, FALSE),
#    has_invitro = coalesce(has_invitro, FALSE),
#    has_any_data = has_tam_drsv | has_invitro
#  )

# latest version considers repeat dose DRSVs only for analogues
data_available_flag <- tam_data_ids %>%
  mutate(has_any_data = TRUE)

# Rebuild rocks IDs from the current rocks object

rocks_ids <- rocks %>%
  clean_names() %>%
  normalize_id_types() %>%
  filter(!is.na(cichem_id)) %>%
  distinct(cichem_id) %>%
  mutate(on_rocks = TRUE)


# Rebuild fp_all from the corrected fingerprints object

fp_all <- fingerprints %>%
  clean_names() %>%
  normalize_id_types() %>%
  left_join(
    data_available_flag,
    by = "cichem_id"
  ) %>%
  left_join(
    rocks_ids,
    by = "cichem_id"
  ) %>%
  mutate(
  #  has_tam_drsv = coalesce(has_tam_drsv, FALSE),
  #  has_invitro = coalesce(has_invitro, FALSE),
    has_any_data = coalesce(has_any_data, FALSE),
    on_rocks = coalesce(on_rocks, FALSE),
    has_any_data = as.integer(has_any_data),
    on_rocks = as.integer(on_rocks)
  ) %>%
  filter(
    has_any_data == 1L |
      on_rocks == 1L
  )

saveRDS(
  fp_all,
  file.path(derived_dir, "fp_all_anydata_rocks_flag.rds")
)


# ===========================================================
# 7) QA summaries
# ===========================================================

availability_summary <- fp_all %>%
  summarise(
    n_total_fingerprinted = n(),
    n_with_any_data = sum(has_any_data == 1L, na.rm = TRUE),
    n_without_any_data = sum(has_any_data == 0L, na.rm = TRUE),
    n_on_rocks = sum(on_rocks == 1L, na.rm = TRUE),
    n_rocks_with_any_data = sum(on_rocks == 1L & has_any_data == 1L, na.rm = TRUE),
    n_rocks_without_any_data = sum(on_rocks == 1L & has_any_data == 0L, na.rm = TRUE),
    pct_rocks_without_any_data = 100 * n_rocks_without_any_data / n_on_rocks
  )

source_target_summary <- fp_all %>%
  mutate(
    analysis_group = case_when(
      on_rocks == 1L & has_any_data == 0L ~ "target_rocks_without_any_data",
      has_any_data == 1L ~ "source_with_any_data",
      TRUE ~ "neither_source_nor_target"
    )
  ) %>%
  count(analysis_group, name = "n_chemicals") %>%
  mutate(
    pct_chemicals = 100 * n_chemicals / sum(n_chemicals)
  )

source_rocks_status <- fp_all %>%
  filter(has_any_data == 1L) %>%
  count(on_rocks, name = "n_sources") %>%
  mutate(
    rocks_status = ifelse(
      on_rocks == 1L,
      "source_on_rocks",
      "source_not_on_rocks"
    ),
    pct_sources = 100 * n_sources / sum(n_sources)
  ) %>%
  select(
    rocks_status,
    on_rocks,
    n_sources,
    pct_sources
  )

saveRDS(
  availability_summary,
  file.path(derived_dir, "availability_summary_anydata_tamdrsv_invitro.rds")
)

saveRDS(
  source_target_summary,
  file.path(derived_dir, "source_target_summary_rocks_targets_anydata_sources.rds")
)

saveRDS(
  source_rocks_status,
  file.path(derived_dir, "source_rocks_status_anydata_tamdrsv_invitro.rds")
)

print(availability_summary)
print(source_target_summary)
print(source_rocks_status)


# ===========================================================
# 8) Build sparse fingerprint matrix
# ===========================================================

build_sparse_from_fp <- function(fp_strings) {

  split_bits <- strsplit(fp_strings, "")

  i <- integer()
  j <- integer()

  for (row_idx in seq_along(split_bits)) {

    on_bits <- which(
      split_bits[[row_idx]] == "1"
    )

    if (length(on_bits) > 0) {

      i <- c(
        i,
        rep.int(row_idx, length(on_bits))
      )

      j <- c(
        j,
        on_bits
      )
    }
  }

  sparseMatrix(
    i = i,
    j = j,
    x = 1L,
    dims = c(length(fp_strings), 2048)
  )
}

fp_sparse <- build_sparse_from_fp(
  fp_all$fp_bit
)

gc()


# ===========================================================
# 9) Define source and target sets
# ===========================================================

source_idx <- which(fp_all$has_any_data == 1L)

target_idx <- which(
  fp_all$on_rocks == 1L &
    fp_all$has_any_data == 0L
)

source_ids <- fp_all$cichem_id[source_idx]
target_ids <- fp_all$cichem_id[target_idx]

source_mat <- fp_sparse[source_idx, , drop = FALSE]
target_mat <- fp_sparse[target_idx, , drop = FALSE]

source_pop <- Matrix::rowSums(source_mat)
target_pop <- Matrix::rowSums(target_mat)

message("Number of source chemicals with any data: ", length(source_ids))
message("Number of ROCKS target chemicals without any data: ", length(target_ids))
message(
  "Number of non-ROCKS source chemicals: ",
  sum(fp_all$has_any_data == 1L & fp_all$on_rocks == 0L, na.rm = TRUE)
)


# ===========================================================
# 10) Chunked Tanimoto/Jaccard matching function
# ===========================================================

get_tanimoto_matches_chunked <- function(
    target_mat,
    source_mat,
    target_ids,
    source_ids,
    target_pop,
    source_pop,
    threshold = 0.7,
    chunk_size = 500
) {

  n_targets <- NROW(target_mat)

  n_sources <- NROW(source_mat)

  if (n_targets == 0 || n_sources == 0) {
    return(
      data.frame(
        target_id = character(),
        source_id = character(),
        tanimoto = numeric()
      )
    )
  }

  starts <- seq.int(
    1L,
    n_targets,
    by = chunk_size
  )

  out <- vector(
    "list",
    length(starts)
  )

  for (k in seq_along(starts)) {

    start <- starts[k]

    end <- min(
      start + chunk_size - 1L,
      n_targets
    )

    chunk <- target_mat[start:end, , drop = FALSE]

    chunk_ids <- target_ids[start:end]

    chunk_pop <- target_pop[start:end]

    inter <- chunk %*% Matrix::t(source_mat)

    inter <- as.matrix(inter)

    union <- outer(
      chunk_pop,
      source_pop,
      "+"
    ) - inter

    sim <- inter / union

    sim[union == 0] <- NA_real_

    idx <- which(
      sim >= threshold,
      arr.ind = TRUE
    )

    if (nrow(idx) > 0) {
      out[[k]] <- data.frame(
        target_id = chunk_ids[idx[, 1]],
        source_id = source_ids[idx[, 2]],
        tanimoto = sim[idx]
      )
    }

    if (k %% 10 == 0 || k == length(starts)) {
      cat(
        "Processed chunk",
        k,
        "of",
        length(starts),
        "for threshold",
        threshold,
        "\n"
      )

      gc()
    }
  }

  bind_rows(out)
}


# ===========================================================
# 11) Run threshold based matching at primary threshold
# ===========================================================

tanimoto_matches <- get_tanimoto_matches_chunked(
  target_mat = target_mat,
  source_mat = source_mat,
  target_ids = target_ids,
  source_ids = source_ids,
  target_pop = target_pop,
  source_pop = source_pop,
  threshold = primary_threshold,
  chunk_size = chunk_size
)

tanimoto_matches <- tanimoto_matches %>%
  filter(target_id != source_id)

saveRDS(
  tanimoto_matches,
  file.path(derived_dir, "tanimoto_matches_thr070.rds")
)


# ===========================================================
# 12) Coverage summary at primary threshold
# ===========================================================

coverage_summary <- tanimoto_matches %>%
  distinct(target_id) %>%
  summarise(
    n_targets_with_match = n(),
    n_total_targets = length(unique(target_ids)),
    pct_targets_with_match = 100 * n_targets_with_match / n_total_targets
  )

saveRDS(
  coverage_summary,
  file.path(derived_dir, "coverage_summary_thr070.rds")
)

print(coverage_summary)


# ===========================================================
# 13) Number of source analogues per target
# ===========================================================

target_match_counts <- tanimoto_matches %>%
  count(
    target_id,
    name = "n_matches"
  )

if (nrow(target_match_counts) > 0) {

  target_match_summary <- target_match_counts %>%
    summarise(
      n_targets_with_match = n(),
      min_matches = min(n_matches, na.rm = TRUE),
      p25_matches = as.numeric(quantile(n_matches, 0.25, na.rm = TRUE)),
      median_matches = median(n_matches, na.rm = TRUE),
      mean_matches = mean(n_matches, na.rm = TRUE),
      p75_matches = as.numeric(quantile(n_matches, 0.75, na.rm = TRUE)),
      max_matches = max(n_matches, na.rm = TRUE),
      n_with_1plus = sum(n_matches >= 1),
      n_with_2plus = sum(n_matches >= 2),
      n_with_5plus = sum(n_matches >= 5),
      n_with_10plus = sum(n_matches >= 10),
      n_with_20plus = sum(n_matches >= 20),
      n_with_50plus = sum(n_matches >= 50)
    )

} else {

  target_match_summary <- data.frame(
    n_targets_with_match = 0,
    min_matches = NA_real_,
    p25_matches = NA_real_,
    median_matches = NA_real_,
    mean_matches = NA_real_,
    p75_matches = NA_real_,
    max_matches = NA_real_,
    n_with_1plus = 0,
    n_with_2plus = 0,
    n_with_5plus = 0,
    n_with_10plus = 0,
    n_with_20plus = 0,
    n_with_50plus = 0
  )

}

saveRDS(
  target_match_counts,
  file.path(derived_dir, "target_match_counts_thr070.rds")
)

saveRDS(
  target_match_summary,
  file.path(derived_dir, "target_match_summary_thr070.rds")
)

# get test list to run genra
targets_20plus_dtxsid <- target_match_counts %>%
  filter(n_matches >= 20) %>%
  left_join(
    fp_all %>%
      select(cichem_id, dtxsid) %>%
      distinct(),
    by = c("target_id" = "cichem_id")
  ) %>%
  pull(dtxsid) %>%
  unique() %>%
  sort()

cat(
  paste(
    targets_20plus_dtxsid,
    collapse = ","
  )
)

# ===========================================================
# 14) Coverage across multiple thresholds
# ===========================================================

coverage_by_threshold <- lapply(
  thresholds,
  function(thr) {

    if (thr == primary_threshold) {

      matches <- tanimoto_matches

    } else {

      matches <- get_tanimoto_matches_chunked(
        target_mat = target_mat,
        source_mat = source_mat,
        target_ids = target_ids,
        source_ids = source_ids,
        target_pop = target_pop,
        source_pop = source_pop,
        threshold = thr,
        chunk_size = chunk_size
      )

    }

    data.frame(
      threshold = thr,
      n_targets_with_match = n_distinct(matches$target_id),
      n_total_targets = length(unique(target_ids))
    ) %>%
      mutate(
        pct_targets_with_match = 100 * n_targets_with_match / n_total_targets
      )

  }
) %>%
  bind_rows()

saveRDS(
  coverage_by_threshold,
  file.path(derived_dir, "coverage_by_threshold.rds")
)

print(coverage_by_threshold)


# ===========================================================
# 18) Residual chemicals without similar data rich analogues
# ===========================================================

covered_now <- unique(
  tanimoto_matches$target_id
)

residual_ids <- setdiff(
  unique(target_ids),
  covered_now
)

summary_now <- data.frame(
  n_total_targets = length(unique(target_ids)),
  n_covered_now = length(covered_now),
  n_residual_uncovered = length(residual_ids),
  pct_covered_now = 100 * length(covered_now) / length(unique(target_ids)),
  pct_residual_uncovered = 100 * length(residual_ids) / length(unique(target_ids))
)

saveRDS(
  residual_ids,
  file.path(derived_dir, "residual_uncoveredthr070.rds")
)

saveRDS(
  summary_now,
  file.path(derived_dir, "summary_now_thr070.rds")
)

print(summary_now)


# ===========================================================
# 19) Attach source and target metadata
# ===========================================================

target_metadata <- fp_all %>%
  filter(cichem_id %in% unique(target_ids)) %>%
  select(
    cichem_id,
    has_any_data,
    on_rocks,
    everything(),
    -fp_bit
  )

source_metadata <- fp_all %>%
  filter(cichem_id %in% unique(source_ids)) %>%
  select(
    cichem_id,
    has_any_data,
    on_rocks,
    everything(),
    -fp_bit
  )

saveRDS(
  target_metadata,
  file.path(derived_dir, "target_metadata_without.rds")
)

saveRDS(
  source_metadata,
  file.path(derived_dir, "source_metadata_with.rds")
)

##

best_neighbor <- tanimoto_matches %>%
  group_by(target_id) %>%
  summarise(
    best_similarity = max(tanimoto)
  )
summary(best_neighbor$best_similarity)

# ===========================================================
# 20) Completion summary
# ===========================================================

message("\n03_build_similarity_coverage_objects.R completed successfully")

message("\nAvailability summary:")
print(availability_summary)

message("\nSource and target summary:")
print(source_target_summary)

message("\nPrimary threshold coverage summary:")
print(coverage_summary)

message("\nTarget match summary:")
print(target_match_summary)

message("\nCoverage by threshold:")
print(coverage_by_threshold)

message("\nResidual target summary:")
print(summary_now)

message("\nSaved derived objects to:")

message(" - ", file.path(derived_dir, "fp_all_hed_aed_hum_flag.rds"))
message(" - ", file.path(derived_dir, "availability_summary_hed_aed_hum.rds"))
message(" - ", file.path(derived_dir, "source_target_summary_hed_aed_hum.rds"))
message(" - ", file.path(derived_dir, "tanimoto_matches_hed_aed_hum_thr070.rds"))
message(" - ", file.path(derived_dir, "coverage_summary_hed_aed_hum_thr070.rds"))
message(" - ", file.path(derived_dir, "target_match_counts_hed_aed_hum_thr070.rds"))
message(" - ", file.path(derived_dir, "target_match_summary_hed_aed_hum_thr070.rds"))
message(" - ", file.path(derived_dir, "coverage_by_threshold_hed_aed_hum.rds"))
message(" - ", file.path(derived_dir, "residual_uncovered_hed_aed_hum_thr070.rds"))
message(" - ", file.path(derived_dir, "summary_now_hed_aed_hum_thr070.rds"))
message(" - ", file.path(derived_dir, "target_metadata_without_hed_aed_hum.rds"))
message(" - ", file.path(derived_dir, "source_metadata_with_hed_aed_hum.rds"))