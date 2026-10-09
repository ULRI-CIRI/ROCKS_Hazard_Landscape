#Uterotrophic, Hershberger, any multi-generation reproductive or EOGRT studies 
# in TAMs to quantify

rocks <- readRDS(
  here::here("data", "derived", "rocks.rds")
)

tam <- readRDS(
  here::here("data", "derived", "tam.rds")
)

drsvs<- tam %>%
  filter(toxval_type_supercategory=="Dose Response Summary Value")

hershutero<- tam %>% filter(study_type=="in vivo endocrine screen")

mgtg <- drsvs %>%
  filter(
    grepl(
      "two-generation|multi-generation|416|443|extended one|870.3800",
      guideline,
      ignore.case = TRUE
    ) |
    grepl(
      "multi-gen|extended one-gen",
      study_duration_class,
      ignore.case = TRUE
    ) |
    grepl(
      "multi-gen|extended one-gen",
      study_duration_class_original,
      ignore.case = TRUE
  ) 
)%>%
  filter(!is.na(cichem_id))%>%
  select(c(source,cichem_id,dtxsid,name, toxval_type,toxval_numeric,toxval_units,study_type,duration_days,study_duration_class,study_duration_class_original,guideline, generation, lifestage, sex, toxicological_effect))

# number of unique cichem_id values
mgtg %>%
  summarise(n = n_distinct(cichem_id)) #1984

# Unique chemical counts
#
# Denominators
n_rocks <- n_distinct(rocks$cichem_id)

n_drsv <- drsvs %>%
  distinct(cichem_id) %>%
  nrow()

# Unique chemicals in each category
endocrine_ids <- hershutero %>%
  distinct(cichem_id) %>%
  pull(cichem_id)

mgtg_ids <- mgtg %>%
  distinct(cichem_id) %>%
  pull(cichem_id)

combined_ids <- union(
  endocrine_ids,
  mgtg_ids
)

summary_tbl <- tibble(
  data_stream = c(
    "Endocrine screens (Hershberger/Uterotrophic)",
    "Multi-generation reproductive or EOGRT studies",
    "Either endocrine screen or multi-generation/EOGRT"
  ),
  n_chemicals = c(
    length(endocrine_ids),
    length(mgtg_ids),
    length(combined_ids)
  )
) %>%
  mutate(
    pct_rocks = 100 * n_chemicals / n_rocks,
    pct_drsv = 100 * n_chemicals / n_drsv
  )

summary_tbl
