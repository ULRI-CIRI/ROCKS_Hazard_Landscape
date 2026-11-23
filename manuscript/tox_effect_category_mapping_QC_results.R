library(readxl)
library(dplyr)
library(tidyr)

qc <- read_excel(here::here("data", "raw", "qc_of_mapped_tox_effects.xlsx")) %>%
  clean_names()

# 1. Rows where human and model categories both exist
human_model_comp <- qc %>%
  filter(
    !is.na(human_category),
    !is.na(model_category)
  )

# How many human vs model disagreements?
human_model_disagreements <- human_model_comp %>%
  filter(human_category != model_category)

n_human_model_disagreements <- nrow(human_model_disagreements)

n_human_model_disagreements

disagreement_summary <- human_model_disagreements %>%
  filter(!is.na(final_category)) %>%
  summarise(
    n_disagreements = n(),

    n_final_agrees_model = sum(final_category == model_category),
    pct_final_agrees_model = 100 * mean(final_category == model_category),

    n_final_agrees_human = sum(final_category == human_category),
    pct_final_agrees_human = 100 * mean(final_category == human_category),

    n_final_agrees_neither = sum(
      final_category != model_category &
        final_category != human_category
    ),
    pct_final_agrees_neither = 100 *
      mean(
        final_category != model_category &
          final_category != human_category
      )
  )

disagreement_summary

model_only_summary <- qc %>%
  filter(
    is.na(human_category),
    !is.na(model_category),
    !is.na(final_category)
  ) %>%
  summarise(
    n_model_only = n(),

    n_agreement = sum(final_category == model_category),
    pct_agreement = 100 * mean(final_category == model_category),

    n_disagreement = sum(final_category != model_category),
    pct_disagreement = 100 * mean(final_category != model_category)
  )

model_only_summary

n_human_model_disagreements
print(disagreement_summary, width = Inf)
model_only_summary
