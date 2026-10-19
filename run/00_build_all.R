# ===========================================================
# 00_build_all.R
#
# Purpose:
#   Run all build scripts required for the Hazard Landscape
#   workflow.
# ===========================================================

message("Starting full build pipeline")

source(here::here("run", "01_build_core_datasets.R"))

source(here::here("run", "02_convert_to_ebmdhed.R"))

source(here::here("run", "03_build_similarity_coverage_objects.R"))

message("Build pipeline completed successfully")
