# ===========================================================
# 00_build_all.R
#
# Purpose:
#   Run all build scripts required for the Hazard Landscape
#   workflow.
# ===========================================================

message("Starting full build pipeline")

source(here::here("run", "01_build_core_datasets.R"))

source(here::here("run", "02_build_potency_objects.R"))

source(here::here("run", "02b_convert_to_ebmdhed.R"))

source(here::here("run", "03_build_similarity_coverage_objects.R"))

source(here::here("run", "04_map_annex_categories.R"))

message("Build pipeline completed successfully")
