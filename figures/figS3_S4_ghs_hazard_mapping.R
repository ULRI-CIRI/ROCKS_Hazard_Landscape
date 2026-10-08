# run 01 to get rocks_ids

rocks_ids_inchi<-rocks_ids %>% select(inchikey)

writeLines(rocks_ids_inchi$inchikey, "inchikeys.txt")

# run through PubChem ID exchange to get CIDs
# https://pubchem.ncbi.nlm.nih.gov/idexchange/idexchange.cgi

# output from PubChem
cids <- tibble::tibble(
  cid = readLines(here::here("data", "raw", "cids.txt"))
)

library(tidyr)

cids <- cids %>%
  separate(
    cid,
    into = c("inchikey", "cid"),
    sep = "\t"
  )

## Run "pulling ghs.py" to get ghs_results.csv output ##

# =============================================================================
# Figure: GHS health hazard categories by chemical class
# =============================================================================

library(tidyverse)
library(here)
library(stringr)
library(scales)

# -----------------------------------------------------------------------------
# 1. File paths
# -----------------------------------------------------------------------------

ghs_file <- here("data", "derived", "ghs_results.csv")
cid_file <- here("data", "raw", "cids.txt")
class_file <- here("data", "derived", "chem_class.rds")

fig_dir <- here("outputs", "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

out_summary <- here(
  "data",
  "derived",
  "ghs_health_category_by_chemical_class_summary.csv"
)

out_heatmap_pct <- here(
  "figures",
  "ghs_health_category_by_chemical_class_heatmap_pct.png"
)

out_heatmap_pct_pdf <- here(
  "figures",
  "ghs_health_category_by_chemical_class_heatmap_pct.pdf"
)

out_heatmap_n <- here(
  "figures",
  "ghs_health_category_by_chemical_class_heatmap_counts.png"
)

out_heatmap_n_pdf <- here(
  "figures",
  "ghs_health_category_by_chemical_class_heatmap_counts.pdf"
)

# -----------------------------------------------------------------------------
# 2. Read data
# -----------------------------------------------------------------------------

ghs_raw <- readr::read_csv(
  ghs_file,
  show_col_types = FALSE
) 

cid_xwalk <- readr::read_tsv(
  cid_file,
  col_names = c("inchikey", "cid"),
  show_col_types = FALSE
) %>%
  mutate(
    cid = as.character(cid),
    inchikey = as.character(inchikey)
  ) %>%
  distinct(inchikey, cid)

chem_class <- readRDS(class_file)

# Standardize ID/class columns
chem_class <- chem_class %>%
  mutate(
    across(
      any_of(c("cid", "inchikey", "dtxsid", "cichem_id", "class", "superclass")),
      as.character
    )
  )


# -----------------------------------------------------------------------------
# 3. Basic checks
# -----------------------------------------------------------------------------

required_ghs_cols <- c("cid", "hazard_statements")

missing_ghs_cols <- setdiff(required_ghs_cols, names(ghs_raw))

if (length(missing_ghs_cols) > 0) {
  stop(
    "Missing required columns in ghs_results.csv: ",
    paste(missing_ghs_cols, collapse = ", ")
  )
}

if (!"cid" %in% names(cid_xwalk) || !"inchikey" %in% names(cid_xwalk)) {
  stop("cid_xwalk must contain columns named 'inchikey' and 'cid'.")
}

if (!any(c("class", "superclass") %in% names(chem_class))) {
  stop("chem_class must contain at least one of: 'class' or 'superclass'.")
}

# Choose class variable
class_var <- "superclass"

if (!class_var %in% names(chem_class)) {
  class_var <- "superclass"
}

message("Using chemical class variable: ", class_var)

# -----------------------------------------------------------------------------
# 4. Attach InChIKey to GHS results
# -----------------------------------------------------------------------------

ghs_raw <- ghs_raw %>%
  mutate(cid = as.character(cid)) %>%
  left_join(
    cid_xwalk,
    by = "cid"
  ) %>%
  rename(inchikey=inchikey.y)


# -----------------------------------------------------------------------------
# 5. Join chemical class
# -----------------------------------------------------------------------------
chem_class <- readRDS(
  here("data", "derived", "chem_class.rds")
)
rocks_ids <- readRDS(
  here("data", "derived", "rocks.rds")
)

# bridge PubChem back to cichem_id
ghs_class <- ghs_raw %>%
  mutate(cid = as.character(cid)) %>%
  left_join(
    rocks_ids %>%
      select(cichem_id, inchikey) %>%
      distinct(),
    by = "inchikey"
  ) %>%
  left_join(
    chem_class,
    by = "cichem_id"
  )

# Remove unclassified/missing class labels
ghs_class <- ghs_class %>%
  mutate(
    chemical_class = .data[[class_var]],
    chemical_class = na_if(chemical_class, ""),
    chemical_class = na_if(chemical_class, "NA"),
    chemical_class = na_if(chemical_class, "Unclassified"),
    chemical_class = na_if(chemical_class, "unclassified")
  ) %>%
  filter(!is.na(chemical_class))

# -----------------------------------------------------------------------------
# 6. Extract H-codes
# -----------------------------------------------------------------------------

ghs_long <- ghs_class %>%
  mutate(
    cid = as.character(cid),
    hazard_statements = as.character(hazard_statements)
  ) %>%
  filter(!is.na(hazard_statements), hazard_statements != "") %>%
  mutate(
    h_code_list = str_extract_all(
      hazard_statements,
      "\\bH[0-9]{3}[A-Za-z]*\\b"
    )
  ) %>%
  select(
    cid,
    inchikey,
    any_of(c("cichem_id", "dtxsid")),
    compound_name,
    signal_word,
    chemical_class,
    hazard_statements,
    h_code_list
  ) %>%
  unnest(h_code_list) %>%
  rename(h_code = h_code_list) %>%
  distinct(cid, chemical_class, h_code, .keep_all = TRUE)

# -----------------------------------------------------------------------------
# 7. Map H-codes to broad GHS health hazard categories
# -----------------------------------------------------------------------------

map_ghs_health_category <- function(h_code) {

  case_when(

    # Acute Toxicity
    h_code %in% c(
      "H300","H301","H302",
      "H310","H311","H312",
      "H330","H331","H332","H333"
    ) ~ "Acute toxicity",

    # Skin Corrosion / Irritation
    h_code %in% c(
      "H314","H315","H316"
    ) ~ "Skin corrosion / irritation",

    # Serious Eye Damage / Eye Irritation
    h_code %in% c(
      "H318","H319","H320"
    ) ~ "Serious eye damage / eye irritation",

    # Respiratory Sensitization
    h_code == "H334" ~ "Respiratory sensitization",

    # Skin Sensitization
    h_code == "H317" ~ "Skin sensitization",

    # Germ Cell Mutagenicity
    h_code %in% c(
      "H340","H341"
    ) ~ "Germ cell mutagenicity",

    # Carcinogenicity
    str_detect(
      h_code,
      "^H350[A-Za-z]*$"
    ) |
      h_code == "H351" ~ "Carcinogenicity",

    # Reproductive Toxicity
    str_detect(
      h_code,
      "^H360[A-Za-z]*$"
    ) |
      str_detect(
        h_code,
        "^H361[A-Za-z]*$"
      ) |
      h_code == "H362" ~ "Reproductive toxicity",

    # STOT Single Exposure
    h_code %in% c(
      "H335","H336",
      "H370","H371"
    ) ~ "STOT-SE",

    # STOT Repeated Exposure
    h_code %in% c(
      "H372","H373"
    ) ~ "STOT-RE",

    # Aspiration Hazard
    h_code == "H304" ~ "Aspiration hazard",

    TRUE ~ "Other / non-health or unmapped"
  )
}

map_ghs_signal_word <- function(h_code) {
  case_when(
    # Acute toxicity
    h_code %in% c("H300", "H301", "H310", "H311", "H330", "H331") ~ "Danger",
    h_code %in% c("H302", "H312", "H332") ~ "Warning",

    # Aspiration
    h_code == "H304" ~ "Danger",

    # Skin / eye
    h_code %in% c("H314", "H318") ~ "Danger",
    h_code %in% c("H315", "H316", "H319", "H320") ~ "Warning",

    # Sensitization
    h_code == "H334" ~ "Danger",
    h_code == "H317" ~ "Warning",

    # STOT single exposure
    h_code == "H370" ~ "Danger",
    h_code %in% c("H335", "H336", "H371") ~ "Warning",

    # Mutagenicity
    h_code == "H340" ~ "Danger",
    h_code == "H341" ~ "Warning",

    # Carcinogenicity
    str_detect(h_code, "^H350[A-Za-z]*$") ~ "Danger",
    h_code == "H351" ~ "Warning",

    # Reproductive toxicity
    str_detect(h_code, "^H360[A-Za-z]*$") ~ "Danger",
    str_detect(h_code, "^H361[A-Za-z]*$") ~ "Warning",
    h_code == "H362" ~ "Warning",

    # STOT repeated exposure
    h_code == "H372" ~ "Danger",
    h_code == "H373" ~ "Warning",

    TRUE ~ NA_character_
  )
}

ghs_health <- ghs_long %>%
  mutate(
    ghs_health_category = map_ghs_health_category(h_code),
    ghs_signal_word = map_ghs_signal_word(h_code)
  ) %>%
  filter(
    ghs_health_category != "Other / non-health or unmapped",
    !is.na(ghs_signal_word)
  ) %>%
  mutate(
    hazard_severity = paste0(
      ghs_health_category,
      " (",
      substr(ghs_signal_word, 1, 1),
      ")"
    )
  )

# -----------------------------------------------------------------------------
# 8. Class-by-hazard-severity summary
# -----------------------------------------------------------------------------

class_denominators <- ghs_class %>%
  mutate(
    chemical_class = as.character(chemical_class)
  ) %>%
  distinct(cid, chemical_class) %>%
  count(
    chemical_class,
    name = "n_chemicals_in_class"
  )

class_hazard_summary <- ghs_health %>%
  mutate(
    chemical_class = as.character(chemical_class),
    hazard_severity = as.character(hazard_severity)
  ) %>%
  distinct(
    cid,
    chemical_class,
    hazard_severity,
    h_code
  ) %>%
  group_by(
    chemical_class,
    hazard_severity
  ) %>%
  summarise(
    n_chemicals = n_distinct(cid),
    n_distinct_h_codes = n_distinct(h_code),
    h_codes = paste(sort(unique(h_code)), collapse = "; "),
    .groups = "drop"
  ) %>%
  left_join(
    class_denominators,
    by = "chemical_class"
  ) %>%
  mutate(
    pct_of_class = n_chemicals / n_chemicals_in_class
  ) %>%
  arrange(
    desc(n_chemicals_in_class),
    chemical_class,
    hazard_severity
  )

# -----------------------------------------------------------------------------
# 8A. Superclass label counts for heatmap axis
# -----------------------------------------------------------------------------
# X = number of unique chemicals in the superclass with at least one mapped
#     health GHS hazard statement
# Y = total number of unique chemicals in the superclass
selected_ghs_categories <- c(
  "Reproductive toxicity (D)",
  "Reproductive toxicity (W)",
  "Germ cell mutagenicity (D)",
  "Germ cell mutagenicity (W)",
  "Carcinogenicity (D)",
  "Carcinogenicity (W)",
  "STOT-RE (D)",
  "STOT-RE (W)"
)

class_ghs_counts <- ghs_health %>%
  filter(
    hazard_severity %in% selected_ghs_categories
  ) %>%
  mutate(
    chemical_class = as.character(chemical_class)
  ) %>%
  distinct(
    cid,
    chemical_class
  ) %>%
  count(
    chemical_class,
    name = "n_with_ghs_health_statement"
  )

class_label_lookup <- class_denominators %>%
  mutate(
    chemical_class = as.character(chemical_class)
  ) %>%
  left_join(
    class_ghs_counts,
    by = "chemical_class"
  ) %>%
  mutate(
    n_with_ghs_health_statement = replace_na(
      n_with_ghs_health_statement,
      0L
    ),
chemical_class_label = paste0(
  chemical_class,
  " (",
  scales::comma(n_chemicals_in_class),
  ")"
)
  )

# -----------------------------------------------------------------------------
# 9. Prepare complete matrix for heatmap
# -----------------------------------------------------------------------------

all_classes <- class_denominators %>%
  mutate(
    chemical_class = as.character(chemical_class)
  ) %>%
  arrange(desc(n_chemicals_in_class)) %>%
  pull(chemical_class)

all_hazard_severity <- c(
  "Acute toxicity (D)",
  "Acute toxicity (W)",

  "Skin corrosion / irritation (D)",
  "Skin corrosion / irritation (W)",

  "Serious eye damage / eye irritation (D)",
  "Serious eye damage / eye irritation (W)",

  "Respiratory sensitization (D)",
  "Skin sensitization (W)",

  "Germ cell mutagenicity (D)",
  "Germ cell mutagenicity (W)",

  "Carcinogenicity (D)",
  "Carcinogenicity (W)",

  "Reproductive toxicity (D)",
  "Reproductive toxicity (W)",

  "STOT-SE (D)",
  "STOT-SE (W)",

  "STOT-RE (D)",
  "STOT-RE (W)",

  "Aspiration hazard (D)"
)

heatmap_dat <- tidyr::expand_grid(
  chemical_class = all_classes,
  hazard_severity = all_hazard_severity
) %>%
  left_join(
    class_hazard_summary,
    by = c(
      "chemical_class",
      "hazard_severity"
    )
  ) %>%
  left_join(
    class_denominators,
    by = "chemical_class",
    suffix = c("", "_denom")
  ) %>%
  mutate(
    n_chemicals = replace_na(n_chemicals, 0L),
    n_distinct_h_codes = replace_na(n_distinct_h_codes, 0L),
    h_codes = replace_na(h_codes, ""),
    n_chemicals_in_class = coalesce(
      n_chemicals_in_class,
      n_chemicals_in_class_denom
    ),
    pct_of_class = n_chemicals / n_chemicals_in_class,
    chemical_class = factor(
      chemical_class,
      levels = rev(all_classes)
    ),
    hazard_severity = factor(
      hazard_severity,
      levels = all_hazard_severity
    )
  )

min_class_n <- 10

heatmap_dat_plot <- heatmap_dat %>%
  filter(n_chemicals_in_class >= min_class_n)

# -----------------------------------------------------------------------------
# Restrict to selected GHS categories for main figure
# -----------------------------------------------------------------------------

selected_ghs_categories <- c(
  "Reproductive toxicity (D)",
  "Reproductive toxicity (W)",
  "Germ cell mutagenicity (D)",
  "Germ cell mutagenicity (W)",
  "Carcinogenicity (D)",
  "Carcinogenicity (W)",
  "STOT-RE (D)",
  "STOT-RE (W)"
)

heatmap_dat_plot <- heatmap_dat_plot %>%
  filter(
    hazard_severity %in% selected_ghs_categories
  ) %>%
  mutate(
    hazard_severity = factor(
      as.character(hazard_severity),
      levels = selected_ghs_categories
    )
  )

# -----------------------------------------------------------------------------
# 10A. Cluster chemical classes AND GHS hazard-severity categories
# -----------------------------------------------------------------------------
# Cluster rows:
#   chemical classes by their GHS hazard-severity profiles.
#
# Cluster columns:
#   GHS hazard-severity categories by their distribution across chemical classes.
#
# This lets both axes be ordered by similarity rather than preserving the
# original GHS category order.
# -----------------------------------------------------------------------------

cluster_input <- heatmap_dat_plot %>%
  select(
    chemical_class,
    hazard_severity,
    pct_of_class
  ) %>%
  distinct() %>%
  mutate(
    chemical_class = as.character(chemical_class),
    hazard_severity = as.character(hazard_severity)
  ) %>%
  pivot_wider(
    names_from = hazard_severity,
    values_from = pct_of_class,
    values_fill = 0
  )

# Convert to class x hazard matrix
cluster_mat <- cluster_input %>%
  tibble::column_to_rownames("chemical_class") %>%
  as.matrix()

# QC checks
stopifnot(!is.null(rownames(cluster_mat)))
stopifnot(!any(rownames(cluster_mat) == ""))
stopifnot(!is.null(colnames(cluster_mat)))
stopifnot(!any(colnames(cluster_mat) == ""))

# -----------------------------
# Row clustering: chemical class
# -----------------------------

row_hclust <- hclust(
  dist(cluster_mat, method = "euclidean"),
  method = "ward.D2"
)

clustered_classes <- row_hclust$labels[row_hclust$order]

# -----------------------------------------------------------------------------
# 10B. Heatmap: percent of class with each GHS hazard-severity category
# -----------------------------------------------------------------------------
ulri_blue_seq <- c(
  "#D9D9D9", # light neutral
  "#BCE4F7", # pale blue
  "#00CAFB", # cyan
  "#00689C", # medium blue
  "#122C49"  # dark navy
)

# -----------------------------------------------------------------------------
# Add superclass labels with GHS statement coverage: X/Y
# -----------------------------------------------------------------------------

plot_class_levels <- rev(clustered_classes)

plot_class_label_levels <- class_label_lookup$chemical_class_label[
  match(
    plot_class_levels,
    class_label_lookup$chemical_class
  )
]

heatmap_dat_plot <- heatmap_dat_plot %>%
  mutate(
    chemical_class = factor(
      as.character(chemical_class),
      levels = plot_class_levels
    )
  ) %>%
  left_join(
    class_label_lookup %>%
      select(
        chemical_class,
        chemical_class_label,
        n_with_ghs_health_statement,
        n_chemicals_in_class
      ) %>%
      mutate(
        chemical_class = as.character(chemical_class)
      ),
    by = c("chemical_class" = "chemical_class")
  ) %>%
  mutate(
    chemical_class_label = factor(
      chemical_class_label,
      levels = plot_class_label_levels
    )
  )

p_heatmap_pct_signal <- ggplot(
  heatmap_dat_plot,
  aes(
    x = hazard_severity,
    y = chemical_class_label,
    fill = pct_of_class
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.3
  ) +
  scale_x_discrete(drop = FALSE) +
scale_fill_gradientn(
  colours = c(
    "#f7f7f7",  # almost white
    "#D9D9D9",
    "#BCE4F7",
    "#00CAFB",
    "#00689C",
    "#122C49"
  ),
  values = scales::rescale(
    c(0, 0.01, 0.03, 0.07, 0.11, 0.15)
  ),
  limits = c(0, 0.15),
  breaks = c(0, 0.05, 0.10, 0.15),
  labels = scales::label_percent(accuracy = 1)
) +
  labs(
    x = NULL,
    y = "Chemical superclass"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size=14
    ),
    axis.text.y = element_text(size = 14),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

print(p_heatmap_pct_signal)
ghs_a<-p_heatmap_pct_signal

ggsave(
  filename = here(
    "outputs",
    "figures",
    "ghs_hazard_severity_by_chemical_superclass_heatmap_pct.png"
  ),
  plot = p_heatmap_pct_signal,
  width = 9,
  height = 12,
  dpi = 300
)

# =============================================================================
# 11. Enrichment analysis: GHS hazard categories by chemical class
# =============================================================================
#
# Question:
#   Which chemical classes are enriched for each GHS health hazard-severity
#   category relative to the rest of the GHS-mapped inventory?
#
# Test:
#   Fisher's exact test for each chemical_class x hazard_severity pair.
#
# 2 x 2 table:
#
#                         Has hazard     Does not have hazard
#   In class                  a                  b
#   Not in class              c                  d
#
# Interpretation:
#   odds_ratio > 1 and FDR < threshold = enriched in class
#   odds_ratio < 1 and FDR < threshold = depleted in class
#
# Notes:
#   - This analysis depends on PubChem GHS labels and should be interpreted as
#     enrichment among chemicals with retrieved GHS classifications, not as true
#     prevalence across the full chemical inventory.
# =============================================================================

library(tidyverse)
library(here)
library(scales)

# -----------------------------------------------------------------------------
# 11A. Create chemical-level hazard presence/absence table
# -----------------------------------------------------------------------------
class_var <- "class"
ghs_class <- ghs_raw %>%
  mutate(cid = as.character(cid)) %>%
  left_join(
    rocks_ids %>%
      select(cichem_id, inchikey) %>%
      distinct(),
    by = "inchikey"
  ) %>%
  left_join(
    chem_class,
    by = "cichem_id"
  )
ghs_class <- ghs_class %>%
  mutate(
    chemical_class = .data[[class_var]],
    chemical_class = na_if(chemical_class, ""),
    chemical_class = na_if(chemical_class, "NA"),
    chemical_class = na_if(chemical_class, "Unclassified"),
    chemical_class = na_if(chemical_class, "unclassified")
  ) %>%
  filter(!is.na(chemical_class))
# Class-level denominators for enrichment analysis
class_denominators_enrich <- ghs_class %>%
  distinct(cid, chemical_class) %>%
  count(
    chemical_class,
    name = "n_chemicals_in_class"
  )
ghs_long <- ghs_class %>%
  mutate(
    cid = as.character(cid),
    hazard_statements = as.character(hazard_statements)
  ) %>%
  filter(!is.na(hazard_statements), hazard_statements != "") %>%
  mutate(
    h_code_list = str_extract_all(
      hazard_statements,
      "\\bH[0-9]{3}[A-Za-z]*\\b"
    )
  ) %>%
  select(
    cid,
    inchikey,
    any_of(c("cichem_id", "dtxsid")),
    compound_name,
    signal_word,
    chemical_class,
    hazard_statements,
    h_code_list
  ) %>%
  unnest(h_code_list) %>%
  rename(h_code = h_code_list) %>%
  distinct(cid, chemical_class, h_code, .keep_all = TRUE)
ghs_health <- ghs_long %>%
  mutate(
    ghs_health_category = map_ghs_health_category(h_code),
    ghs_signal_word = map_ghs_signal_word(h_code)
  ) %>%
  filter(
    ghs_health_category != "Other / non-health or unmapped",
    !is.na(ghs_signal_word)
  ) %>%
  mutate(
    hazard_severity = paste0(
      ghs_health_category,
      " (",
      substr(ghs_signal_word, 1, 1),
      ")"
    )
  )

min_class_n_enrich <- 10

# Chemical-level universe for enrichment analysis
# This represents chemicals with retrieved GHS results and assigned chemical class
ghs_universe <- ghs_class %>%
  mutate(
    cid = as.character(cid),
    chemical_class = as.character(chemical_class)
  ) %>%
  filter(
    !is.na(cid),
    !is.na(chemical_class)
  ) %>%
  distinct(
    cid,
    chemical_class
  )

ghs_universe_plot <- ghs_universe %>%
  semi_join(
    class_denominators_enrich %>%
      filter(n_chemicals_in_class >= min_class_n_enrich) %>%
      select(chemical_class),
    by = "chemical_class"
  )

# Chemical x hazard severity presence
chemical_hazard_presence <- ghs_health %>%
  mutate(
    cid = as.character(cid),
    chemical_class = as.character(chemical_class),
    hazard_severity = as.character(hazard_severity)
  ) %>%
  distinct(
    cid,
    chemical_class,
    hazard_severity
  ) %>%
  semi_join(
    ghs_universe_plot,
    by = c("cid", "chemical_class")
  ) %>%
  mutate(
    has_hazard = TRUE
  )

# Complete chemical x hazard matrix
chemical_hazard_matrix <- tidyr::expand_grid(
  ghs_universe_plot,
  hazard_severity = all_hazard_severity
) %>%
  left_join(
    chemical_hazard_presence,
    by = c(
      "cid",
      "chemical_class",
      "hazard_severity"
    )
  ) %>%
  mutate(
    has_hazard = replace_na(has_hazard, FALSE)
  )

# Optional QC
message(
  "Chemicals in enrichment universe: ",
  scales::comma(n_distinct(chemical_hazard_matrix$cid))
)

message(
  "Classes in enrichment universe: ",
  scales::comma(n_distinct(chemical_hazard_matrix$chemical_class))
)

# -----------------------------------------------------------------------------
# 11B. Fisher enrichment test for each class x hazard category
# -----------------------------------------------------------------------------

run_ghs_enrichment <- function(dat, class_i, hazard_i) {

  test_dat <- dat %>%
    mutate(
      in_class = chemical_class == class_i,
      has_hazard = hazard_severity == hazard_i & has_hazard
    ) %>%
    group_by(cid, chemical_class, in_class) %>%
    summarise(
      has_hazard = any(has_hazard),
      .groups = "drop"
    )

  a <- test_dat %>%
    filter(in_class, has_hazard) %>%
    summarise(n = n_distinct(cid), .groups = "drop") %>%
    pull(n)

  b <- test_dat %>%
    filter(in_class, !has_hazard) %>%
    summarise(n = n_distinct(cid), .groups = "drop") %>%
    pull(n)

  c <- test_dat %>%
    filter(!in_class, has_hazard) %>%
    summarise(n = n_distinct(cid), .groups = "drop") %>%
    pull(n)

  d <- test_dat %>%
    filter(!in_class, !has_hazard) %>%
    summarise(n = n_distinct(cid), .groups = "drop") %>%
    pull(n)

  # Replace length-zero outputs with zero
  a <- ifelse(length(a) == 0, 0, a)
  b <- ifelse(length(b) == 0, 0, b)
  c <- ifelse(length(c) == 0, 0, c)
  d <- ifelse(length(d) == 0, 0, d)

  fisher_mat <- matrix(
    c(a, b, c, d),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      class = c("in_class", "not_in_class"),
      hazard = c("has_hazard", "no_hazard")
    )
  )

  fisher_res <- fisher.test(fisher_mat)

  n_class <- a + b
  n_not_class <- c + d
  n_total <- a + b + c + d

  pct_class <- ifelse(n_class > 0, a / n_class, NA_real_)
  pct_not_class <- ifelse(n_not_class > 0, c / n_not_class, NA_real_)

  expected_in_class <- ifelse(
    n_total > 0,
    n_class * ((a + c) / n_total),
    NA_real_
  )

  tibble(
    chemical_class = class_i,
    hazard_severity = hazard_i,

    n_class = n_class,
    n_not_class = n_not_class,

    n_class_with_hazard = a,
    n_not_class_with_hazard = c,

    pct_class_with_hazard = pct_class,
    pct_not_class_with_hazard = pct_not_class,

    expected_class_with_hazard = expected_in_class,
    observed_expected_ratio = ifelse(
      expected_in_class > 0,
      a / expected_in_class,
      NA_real_
    ),

    odds_ratio = unname(fisher_res$estimate),
    p_value = fisher_res$p.value
  )
}

classes_for_enrichment <- ghs_universe_plot %>%
  distinct(chemical_class) %>%
  arrange(chemical_class) %>%
  pull(chemical_class)

hazards_for_enrichment <- all_hazard_severity

ghs_enrichment <- tidyr::expand_grid(
  chemical_class = classes_for_enrichment,
  hazard_severity = hazards_for_enrichment
) %>%
  mutate(
    res = purrr::map2(
      chemical_class,
      hazard_severity,
      ~ run_ghs_enrichment(
        dat = chemical_hazard_matrix,
        class_i = .x,
        hazard_i = .y
      )
    )
  ) %>%
  select(res) %>%
  unnest(res) %>%
  mutate(
    fdr = p.adjust(p_value, method = "BH"),
    log2_odds_ratio = log2(odds_ratio),
    enrichment_direction = case_when(
      fdr < 0.05 & odds_ratio > 1 ~ "Enriched",
      fdr < 0.05 & odds_ratio < 1 ~ "Depleted",
      TRUE ~ "Not significant"
    )
  ) %>%
  arrange(
    fdr,
    desc(abs(log2_odds_ratio))
  )


# Print top enriched combinations
ghs_enrichment %>%
  filter(
    enrichment_direction == "Enriched"
  ) %>%
  select(
    chemical_class,
    hazard_severity,
    n_class,
    n_class_with_hazard,
    pct_class_with_hazard,
    pct_not_class_with_hazard,
    odds_ratio,
    observed_expected_ratio,
    p_value,
    fdr
  ) %>%
  arrange(fdr, desc(odds_ratio)) %>%
  print(n = 50)

# -----------------------------------------------------------------------------
# 11C. Plot enriched GHS category-class combinations
# -----------------------------------------------------------------------------

ghs_enrichment_plot_dat <- ghs_enrichment %>%
  filter(
    enrichment_direction == "Enriched",
    n_class_with_hazard >= 3
  ) %>%
  mutate(
    chemical_class = forcats::fct_reorder(
      chemical_class,
      log2_odds_ratio,
      .fun = median,
      .desc = FALSE
    ),
    hazard_severity = forcats::fct_reorder(
      hazard_severity,
      log2_odds_ratio,
      .fun = median,
      .desc = TRUE
    ),
    neg_log10_fdr = -log10(fdr)
  )

p_ghs_enrichment <- ggplot(
  ghs_enrichment_plot_dat,
  aes(
    x = hazard_severity,
    y = chemical_class
  )
) +
  geom_point(
    aes(
      size = pct_class_with_hazard,
      color = log2_odds_ratio
    ),
    alpha = 0.9
  ) +
  scale_color_gradientn(
    colors = c(
      "#BCE4F7",
      "#00CAFB",
      "#00689C",
      "#122C49"
    ),
    name = expression(log[2]~"odds ratio")
  ) +
  scale_size_continuous(
    labels = scales::label_percent(accuracy = 1),
    name = "Percent of\nclass"
  ) +
  labs(
    title = "Chemical classes enriched for GHS health hazard categories",
    subtitle = "Fisher's exact test, BH-adjusted FDR < 0.05; points require at least 3 chemicals in class with hazard",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

print(p_ghs_enrichment)

ggsave(
  filename = here(
    "outputs",
    "figures",
    "ghs_class_hazard_enrichment_dotplot.png"
  ),
  plot = p_ghs_enrichment,
  width = 13,
  height = 8,
  dpi = 300
)

ggsave(
  filename = here(
    "outputs",
    "figures",
    "ghs_class_hazard_enrichment_dotplot.pdf"
  ),
  plot = p_ghs_enrichment,
  width = 13,
  height = 8
)

# =============================================================================
# Supplemental GHS enrichment dot plot:
# selected chronic/systemic hazard categories only
# =============================================================================

library(tidyverse)
library(here)
library(scales)

selected_ghs_categories <- c(
  "Reproductive toxicity (D)",
  "Reproductive toxicity (W)",
  "Germ cell mutagenicity (D)",
  "Germ cell mutagenicity (W)",
  "Carcinogenicity (D)",
  "Carcinogenicity (W)",
  "STOT-RE (D)",
  "STOT-RE (W)"
)

ghs_enrichment_selected <- ghs_enrichment %>%
  filter(
    enrichment_direction == "Enriched",
    hazard_severity %in% selected_ghs_categories,
    fdr < 0.01,
n_class_with_hazard >= 5
  ) %>%
  mutate(
    hazard_severity = factor(
      hazard_severity,
      levels = selected_ghs_categories
    )
  )

# Optional: keep only classes that have at least one selected enrichment
classes_to_show <- ghs_enrichment_selected %>%
  distinct(chemical_class) %>%
  pull(chemical_class)

ghs_enrichment_selected <- ghs_enrichment_selected %>%
  filter(
    chemical_class %in% classes_to_show
  ) %>%
  mutate(
    chemical_class = forcats::fct_reorder(
      chemical_class,
      log2_odds_ratio,
      .fun = median,
      .desc = FALSE
    ),
    neg_log10_fdr = -log10(fdr)
  )

p_ghs_enrichment_selected <- ggplot(
  ghs_enrichment_selected,
  aes(
    x = hazard_severity,
    y = chemical_class
  )
) +
  geom_point(
    aes(
      size = pct_class_with_hazard,
      color = log2_odds_ratio
    ),
    alpha = 0.9
  ) +
  scale_color_gradientn(
    colors = c(
      "#BCE4F7",
      "#00CAFB",
      "#00689C",
      "#122C49"
    ),
    name = expression(log[2]~"odds ratio")
  ) +
scale_size_continuous(
  breaks = c(0.05, 0.25),
  labels = scales::label_percent(accuracy = 1),
  name = NULL,
  range = c(1.5, 6)
)+
  labs(
    x = NULL,
    y = "Chemical class"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size=12
    ),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

p_ghs_enrichment_selected
ghs_b<-p_ghs_enrichment_selected

ggsave(
  filename = here(
    "outputs",
    "figures",
    "ghs_class_hazard_enrichment_selected_chronic_categories.png"
  ),
  plot = p_ghs_enrichment_selected,
  width = 9,
  height = 15,
  dpi = 300
)

library(cowplot)

ghs_ab<-plot_grid(
  ghs_a,
  ghs_b,
  labels = c("A", "B"),
  ncol = 2,
  align = "v",
  axis = "lr"
)

ggsave(
  filename = here(
    "outputs",
    "figures",
    "ghs_combined_figure.png"
  ),
  plot = ghs_ab,
  width = 12,
  height = 8,
  dpi = 300
)

# =============================================================================
# GHS coverage summary for ROCKS inventory
# =============================================================================
# Purpose:
#   Quantify the fraction of chemicals on the ROCKS/global inventory with:
#     1. any retrieved PubChem GHS hazard statement
#     2. any mapped GHS health hazard classification
#     3. any selected chronic/systemic GHS health hazard classification
#
# Notes:
#   - Counts are chemical-level counts based on cichem_id.
#   - This avoids double-counting chemicals with multiple PubChem CIDs or
#     multiple GHS hazard statements.
# =============================================================================


# Inventory denominator
inventory_ids <- rocks_ids %>%
  distinct(cichem_id)

n_inventory <- n_distinct(inventory_ids$cichem_id)

# -------------------------------------------------------------------------
# Chemicals with any retrieved PubChem GHS hazard statement
# -------------------------------------------------------------------------

ghs_any_statement_ids <- ghs_raw %>%
  mutate(
    cid = as.character(cid),
    hazard_statements = as.character(hazard_statements)
  ) %>%
  filter(
    !is.na(hazard_statements),
    hazard_statements != ""
  ) %>%
  left_join(
    rocks_ids %>%
      select(cichem_id, inchikey) %>%
      distinct(),
    by = "inchikey"
  ) %>%
  filter(
    !is.na(cichem_id)
  ) %>%
  distinct(cichem_id)

# -------------------------------------------------------------------------
# Chemicals with any mapped GHS health hazard classification
# -------------------------------------------------------------------------

ghs_any_health_ids <- ghs_health %>%
  filter(
    !is.na(cichem_id)
  ) %>%
  distinct(cichem_id)

# -------------------------------------------------------------------------
# Chemicals with one of the selected chronic/systemic hazard categories
# -------------------------------------------------------------------------

ghs_selected_health_ids <- ghs_health %>%
  filter(
    !is.na(cichem_id),
    hazard_severity %in% selected_ghs_categories
  ) %>%
  distinct(cichem_id)

# -------------------------------------------------------------------------
# Summary table
# -------------------------------------------------------------------------

ghs_inventory_coverage <- tibble(
  category = c(
    "Inventory chemicals",
    "Any PubChem GHS hazard statement",
    "Any mapped GHS health hazard classification",
    "Any selected chronic/systemic GHS classification"
  ),
  n_chemicals = c(
    n_inventory,
    n_distinct(ghs_any_statement_ids$cichem_id),
    n_distinct(ghs_any_health_ids$cichem_id),
    n_distinct(ghs_selected_health_ids$cichem_id)
  )
) %>%
  mutate(
    pct_inventory = n_chemicals / n_inventory,
    pct_inventory_label = scales::percent(
      pct_inventory,
      accuracy = 0.1
    )
  )

print(ghs_inventory_coverage)

# coverage from all mapped categories

ghs_all_category_coverage <- ghs_health %>%
  filter(!is.na(cichem_id)) %>%
  distinct(
    cichem_id,
    ghs_health_category
  ) %>%
  count(
    ghs_health_category,
    name = "n_chemicals"
  ) %>%
  mutate(
    pct_inventory = n_chemicals / n_inventory,
    pct_inventory_label = scales::percent(
      pct_inventory,
      accuracy = 0.1
    )
  ) %>%
  arrange(desc(n_chemicals))

print(ghs_all_category_coverage)

# coverage from selected categories
ghs_selected_category_coverage_collapsed <- ghs_health %>%
  filter(
    !is.na(cichem_id),
    hazard_severity %in% selected_ghs_categories
  ) %>%
  distinct(
    cichem_id,
    ghs_health_category
  ) %>%
  count(
    ghs_health_category,
    name = "n_chemicals"
  ) %>%
  mutate(
    pct_inventory = n_chemicals / n_inventory,
    pct_inventory_label = scales::percent(
      pct_inventory,
      accuracy = 0.1
    )
  ) %>%
  arrange(desc(n_chemicals))

print(ghs_selected_category_coverage_collapsed)

