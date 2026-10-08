# Figure 7

# ---- Packages ----
library(DBI)
library(data.table)
library(fastcluster)
library(flextable)
library(ggplot2)
library(ggrepel)
library(ggupset)
library(grid)
library(gridExtra)
library(httr)
library(httr2)
library(jsonlite)
library(dplyr)
library(Matrix)
library(openxlsx)
library(RcppAnnoy)
library(tibble)
library(tidyr)
library(usethis)
library(uwot)
library(viridis)
library(here)

# ---- Paths and options ----
input_file <- here::here("data", "derived", "er_ar_summary_input_06sept26.RData")
output_dir <- here::here(
  "outputs",
  "figures"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}
load(input_file)

# ---- ULRI colors ----
ULcolors <- c(
  "#122C49", "#000095", "#0A32FF", "#00689C", "#577E9E", "#00CAFB",
  "#BCE4F7", "#123A28", "#00948F", "#00A451", "#58595B", "#B38874",
  "#E5DDCA", "#5B0428", "#98443D", "#CA0123", "#FF6255", "#FF9D55"
)

# ---- OPERA v2.9 property coverage ----
message("Chemicals without a molecular formula: ",
        length(unique(rocks_chemistry[is.na(mol_formula), cichem])))

cols_to_check <- names(rocks_amenability)[grepl("_ad$", names(rocks_amenability))]
ad_counts <- rbindlist(lapply(cols_to_check, function(col) {
  tbl <- table(rocks_amenability[[col]], useNA = "always")
  data.table(column = col, value = names(tbl), count = as.integer(tbl))
}))
ad_counts[, value := fifelse(is.na(value) | value == "<NA>", "No Prediction", value)]
ad_counts[, value := factor(value, levels = c("TRUE", "FALSE", "No Prediction"))]

# fig7 panel A
fig_property_coverage <- ggplot(ad_counts, aes(column, count, fill = value)) +
  geom_col() +
  geom_text(aes(label = count, color = value),
            position = position_stack(vjust = 0.5), size = 4, fontface = "bold") +
  scale_fill_manual(
    values = c("TRUE" = "#0A32FF", "FALSE" = "#BCE4F7", "No Prediction" = "#E5DDCA"),
    labels = c("TRUE" = "Inside AD", "FALSE" = "Outside AD", "No Prediction" = "No Prediction"),
    name = "Applicability Domain (AD)"
  ) +
  scale_color_manual(values = c("TRUE" = "white", "FALSE" = "black", "No Prediction" = "black"), guide = "none") +
  scale_x_discrete(labels = c(logp_ad = "LogP", vp_ad = "VP", hl_ad = "HLC", ws_ad = "Water Sol")) +
  labs(x = "Predicted Physicochemical Properties", y = "Chemical Count") +
  theme_minimal(base_size=14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
fig_property_coverage

# Water-solubility exploration: less than 10 uM, assuming log10-M.
water_sol_limit <- log10(0.000010)
ws_limit <- rocks_amenability[!is.na(ws) & ws < water_sol_limit]
message("Chemicals with predicted water solubility below 10 uM: ", nrow(ws_limit))

# ---- In vitro screening amenability ----
rocks_amenability[, in_vitro_screening_amenable := fcase(
  is.na(logp) & is.na(vp) & is.na(hl), "UNKNOWN",
  !is.na(logp_ad) & !logp_ad & !is.na(vp_ad) & !vp_ad & !is.na(hl_ad) & !hl_ad, "UNKNOWN",
  is.na(logp) & is.na(vp) & !is.na(hl),
  fifelse(!is.na(hl_ad) & hl_ad & hl > -3.0, "FALSE", "TRUE"),
  !is.na(average_mass) & average_mass > 100 & !is.na(logp) & logp < 7 &
    !is.na(logp_ad) & logp_ad & !is.na(vp) & vp < 2 & !is.na(vp_ad) & vp_ad, "TRUE",
  !is.na(average_mass) & average_mass > 100 & !is.na(logp) & logp < 7 &
    !is.na(logp_ad) & logp_ad & (is.na(vp_ad) | !vp_ad), "TRUE",
  default = "FALSE"
)]

# ---- Screening-candidate UpSet data ----
combined_endo_summary_classy[, any_positive :=
  as.integer(any(c(agonist, antagonist, binding) == 1, na.rm = TRUE)),
  by = .(cichem_id, model)]
combined_toxcast_summary[, any_positive :=
  as.integer(any(c(agonist, antagonist) >= 0.1, na.rm = TRUE)),
  by = .(cichem_id, model)]

endo_collapsed <- combined_endo_summary_classy[, .(
  endo_any_positive = as.integer(any(any_positive == 1, na.rm = TRUE)),
  endo_model_available = as.integer(any(!is.na(agonist) | !is.na(antagonist) | !is.na(binding)))
), by = cichem_id]

toxcast_collapsed <- combined_toxcast_summary[, .(
  toxcast_any_positive = as.integer(any(any_positive == 1, na.rm = TRUE)),
  toxcast_model_available = as.integer(any(!is.na(agonist) | !is.na(antagonist)))
), by = cichem_id]

upset_data <- merge(
  rocks_amenability[, .(cichem, in_vitro_screening_amenable)], endo_collapsed,
  by.x = "cichem", by.y = "cichem_id", all.x = TRUE
)
upset_data <- merge(
  upset_data, toxcast_collapsed,
  by.x = "cichem", by.y = "cichem_id", all.x = TRUE
)

upset_data[, datasets := Map(
  function(
    endo,
    toxcast,
    endo_avail,
    toxcast_avail,
    amenable
  ) {
    sets <- character()

    if (!is.na(amenable) && amenable == "TRUE") {
      sets <- c(sets, "Predicted Screening Amenable")
    }

    if (!is.na(amenable) && amenable == "UNKNOWN") {
      sets <- c(sets, "Cannot Determine Amenability")
    }

    if (!is.na(endo) && endo == 1) {
      sets <- c(sets, "In Silico Model Positive")
    } else if (!is.na(endo_avail) && endo_avail == 1) {
      sets <- c(sets, "In Silico Model Negative")
    } else {
      sets <- c(sets, "In Silico Model Not Available")
    }

    if (!is.na(toxcast) && toxcast == 1) {
      sets <- c(sets, "ToxCast Model Positive")
    } else if (!is.na(toxcast_avail) && toxcast_avail == 1) {
      sets <- c(sets, "ToxCast Model Negative")
    }

    # Return the character vector directly.
    sets
  },
  endo_any_positive,
  toxcast_any_positive,
  endo_model_available,
  toxcast_model_available,
  in_vitro_screening_amenable
)]


upset_data[, `Screening Candidates` := fifelse(
  vapply(datasets, function(x) {
    "Predicted Screening Amenable" %in% x &&
      "In Silico Model Positive" %in% x &&
      !any(c("Cannot Determine Amenability", "In Silico Model Negative",
             "ToxCast Model Positive", "ToxCast Model Negative") %in% x)
  }, logical(1)), "Screening Candidates", "Other"
)]

upset_data2 <- copy(upset_data)

upset_data2[, datasets := lapply(datasets, function(x) {

  x <- setdiff(x, "Cannot Determine Amenability")

  if (any(c(
    "ToxCast Model Positive",
    "ToxCast Model Negative"
  ) %in% x)) {

    x <- c(
      setdiff(
        x,
        c(
          "ToxCast Model Positive",
          "ToxCast Model Negative"
        )
      ),
      "ToxCast Model Available"
    )
  }

  unique(x)
})]

# fig7 panel B
fig_upset <- ggplot(upset_data2, aes(datasets, fill = `Screening Candidates`)) +
  geom_bar() +
  geom_text(stat = "count", aes(label = after_stat(count)),
            vjust = -0.5, hjust = -0.1, angle = 45, size = 4) +
  scale_x_upset() +
  scale_y_log10(
    breaks = c(1, 5, 10, 50, 100, 500, 1000, 5000, 10000, 50000, 100000),
    labels = scales::comma, expand = expansion(mult = c(0, 0.3))
  ) +
  scale_fill_manual(values = c("Screening Candidates" = "#0A32FF", "Other" = "#BCE4F7"), name = NULL) +
  labs(x = NULL, y = "Chemical Count") +
  theme_minimal(base_size=14) +
  theme_combmatrix(combmatrix.label.text = element_text(size = 12))+
  theme(legend.position = "bottom")
fig_upset

# ---- Extended CERAPP/COMPARA UpSet plot ----
endo_collapsed_by_model <- combined_endo_summary_classy[, .(
  endo_any_positive = as.integer(any(any_positive == 1, na.rm = TRUE)),
  endo_model_available = as.integer(any(!is.na(agonist) | !is.na(antagonist) | !is.na(binding)))
), by = .(cichem_id, model)]

endo_wide <- dcast(endo_collapsed_by_model, cichem_id ~ model,
                   value.var = c("endo_any_positive", "endo_model_available"), fill = NA)
upset_data3 <- merge(rocks_amenability[, .(cichem, in_vitro_screening_amenable)],
                     endo_wide, by.x = "cichem", by.y = "cichem_id", all.x = TRUE)
upset_data3 <- merge(upset_data3, toxcast_collapsed,
                     by.x = "cichem", by.y = "cichem_id", all.x = TRUE)

required_endo_cols <- c("endo_any_positive_CERAPP", "endo_any_positive_COMPARA",
                        "endo_model_available_CERAPP", "endo_model_available_COMPARA")
missing_endo_cols <- setdiff(required_endo_cols, names(upset_data3))
if (length(missing_endo_cols)) stop("Missing expected model columns: ", paste(missing_endo_cols, collapse = ", "))

upset_data3[, datasets := Map(function(cerapp, compara, toxcast, cerapp_avail,
                                        compara_avail, toxcast_avail, amenable) {
  sets <- character()
  if (!is.na(amenable) && amenable == "TRUE") sets <- c(sets, "Predicted Screening Amenable")
  if (!is.na(amenable) && amenable == "UNKNOWN") sets <- c(sets, "Cannot Determine Amenability")
  if (!is.na(cerapp) && cerapp == 1) sets <- c(sets, "CERAPP Positive")
  else if (!is.na(cerapp_avail) && cerapp_avail == 1) sets <- c(sets, "CERAPP Negative")
  if (!is.na(compara) && compara == 1) sets <- c(sets, "COMPARA Positive")
  else if (!is.na(compara_avail) && compara_avail == 1) sets <- c(sets, "COMPARA Negative")
  if (!is.na(toxcast) && toxcast == 1) sets <- c(sets, "ToxCast Model Positive")
  else if (!is.na(toxcast_avail) && toxcast_avail == 1) sets <- c(sets, "ToxCast Model Negative")
  sets
}, endo_any_positive_CERAPP, endo_any_positive_COMPARA, toxcast_any_positive,
   endo_model_available_CERAPP, endo_model_available_COMPARA,
   toxcast_model_available, in_vitro_screening_amenable)]

upset_data3[, datasets := lapply(datasets, function(x) {
  if (any(c("ToxCast Model Positive", "ToxCast Model Negative") %in% x)) {
    x <- c(setdiff(x, c("ToxCast Model Positive", "ToxCast Model Negative")),
           "ToxCast Model Available")
  }
  x
})]
upset_data3[, `Screening Candidates` := fifelse(
  vapply(datasets, function(x) {
    "Predicted Screening Amenable" %in% x &&
      any(c("CERAPP Positive", "COMPARA Positive") %in% x) &&
      !"ToxCast Model Available" %in% x
  }, logical(1)), "Screening Candidates", "Other"
)]

fig_upset3 <- ggplot(upset_data3, aes(datasets, fill = `Screening Candidates`)) +
  geom_bar() +
  geom_text(stat = "count", aes(label = after_stat(count)),
            vjust = -0.5, hjust = -0.1, angle = 45, size = 3) +
  scale_x_upset() +
  scale_y_log10(breaks = c(1, 5, 10, 50, 100, 500, 1000, 5000, 10000, 50000, 100000),
                labels = scales::comma, expand = expansion(mult = c(0, 0.3))) +
  scale_fill_manual(values = c("Screening Candidates" = "#0A32FF", "Other" = "#BCE4F7"), name = NULL) +
  labs(x = "Dataset Intersection", y = "Chemical Count") +
  theme_minimal() +
  theme_combmatrix(combmatrix.label.text = element_text(size = 10))

# ---- Save property-coverage and UpSet figures ----
fig_combined_ab <- arrangeGrob(
  arrangeGrob(fig_property_coverage, top = textGrob("A", x = 0, hjust = -0.5, gp = gpar(fontsize = 16, fontface = "bold"))),
  arrangeGrob(fig_upset, top = textGrob("B", x = 0, hjust = -0.5, gp = gpar(fontsize = 16, fontface = "bold"))),
  ncol = 1
)
tiff(file.path(output_dir, "fig7_06sept2026.tif"), width = 7, height = 9, units = "in", res = 300)
grid.draw(fig_combined_ab)
dev.off()
tiff(file.path(output_dir, "supp_fig_extended_upset_06sept2026.tif"), width = 7, height = 7, units = "in", res = 300)
grid.draw(fig_upset3)
dev.off()

# ---- Screening candidates ----
screening_candidates <- upset_data2[vapply(datasets, function(x) {
  "Predicted Screening Amenable" %in% x &&
    "In Silico Model Positive" %in% x &&
    !"ToxCast Model Available" %in% x
}, logical(1))]
screening_candidates_cichem <- screening_candidates$cichem
fwrite(screening_candidates, file.path(output_dir, "screening_candidates_06sept2026.csv"))

# ---- Morgan fingerprints and saved UMAP ----
rocks_screening <- rocks_morganfp[cichem %in% screening_candidates_cichem]
rocks_screening <- rocks_screening[!is.na(fp) & fp != ""]
message("Screening candidates with fingerprints: ", nrow(rocks_screening))

umap_file <- here::here("data", "derived", "umap_result_06sept2026.rds")
if (!file.exists(umap_file)) stop("UMAP file not found: ", umap_file)
umap_result <- readRDS(umap_file)
if (nrow(umap_result) != nrow(rocks_screening)) {
  stop("UMAP rows (", nrow(umap_result), ") do not match fingerprint rows (", nrow(rocks_screening), ").")
}

umap_dt <- data.table(cichem = rocks_screening$cichem,
                      primary_name = rocks_screening$primary_name,
                      UMAP1 = umap_result[, 1], UMAP2 = umap_result[, 2])
sector_lookup <- unique(combined_toxcast_summary[, .(cichem = cichem_id, sector = sector_class)])
umap_dt <- merge(umap_dt, sector_lookup, by = "cichem", all.x = TRUE)
umap_dt_unique <- unique(umap_dt, by = "cichem")

# ---- k-means and representatives ----
umap_matrix <- as.matrix(umap_dt_unique[, .(UMAP1, UMAP2)])
rownames(umap_matrix) <- umap_dt_unique$cichem
set.seed(42)
wss <- vapply(2:50, function(n_clusters) {
  kmeans(umap_matrix, centers = n_clusters, nstart = 5, iter.max = 100)$tot.withinss
}, numeric(1))
png(file.path(output_dir, "kmeans_elbow_plot.png"), width = 7, height = 5, units = "in", res = 300)
plot(2:50, wss, type = "b", pch = 19, xlab = "Number of clusters k",
     ylab = "Total within-cluster sum of squares", main = "Elbow plot: choose k at the bend")
dev.off()

k <- 20
set.seed(42)
km_result <- kmeans(umap_matrix, centers = k, nstart = 25, iter.max = 200)
umap_dt_unique[, cluster := factor(km_result$cluster)]
centroids <- as.data.table(km_result$centers)
setnames(centroids, c("centroid_UMAP1", "centroid_UMAP2"))
centroids[, cluster := factor(seq_len(.N))]
umap_dt_unique[, dist_to_centroid := {
  center <- centroids[cluster == .BY$cluster]
  sqrt((UMAP1 - center$centroid_UMAP1)^2 + (UMAP2 - center$centroid_UMAP2)^2)
}, by = cluster]
medoid_reps <- umap_dt_unique[, .SD[which.min(dist_to_centroid)], by = cluster]

sector_colors <- c(
  "Consumer Only" = "#000095", "Food Only" = "#00689C", "Industrial Only" = "#E5DDCA",
  "Agricultural Only" = "#00A451", "Multiple Inventories" = "#BCE4F7",
  "Pharmaceutical Only" = "#CA0123", "Unknown" = "#FF9D55"
)
fig_medoid <- ggplot(umap_dt_unique, aes(UMAP1, UMAP2, color = sector)) +
  geom_point(alpha = 0.6, size = 0.3) +
  geom_point(data = medoid_reps, color = "black", shape = 8, size = 3) +
  geom_label_repel(data = medoid_reps, aes(label = cluster), color = "black",
                   fontface = "bold", size = 3, fill = "white", max.overlaps = Inf) +
  scale_color_manual(values = sector_colors) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  theme_minimal() +
  labs(subtitle = paste0("k = ", k, " clusters, representatives marked with a star"),
       x = "UMAP 1", y = "UMAP 2", color = "Sector")
ggsave(file.path(output_dir, "fig7_supplemental_umap_medoids_06sept2026.tif"),
       fig_medoid, width = 8, height = 6, dpi = 300, compression = "lzw")

fig_umap_clusters <- ggplot(umap_dt_unique, aes(UMAP1, UMAP2, color = cluster)) +
  geom_point(alpha = 0.5, size = 0.6) +
  geom_label_repel(data = medoid_reps, aes(label = cluster), color = "black",
                   fontface = "bold", size = 3, fill = "white", max.overlaps = Inf) +
  guides(color = guide_legend(override.aes = list(size = 5, alpha = 1))) +
  theme_minimal() +
  labs(title = paste0("UMAP colored by k-means cluster (k = ", k, ")"),
       x = "UMAP 1", y = "UMAP 2", color = "Cluster")
ggsave(file.path(output_dir, "fig_supplemental_umap_clusters_06sept2026.tif"),
       fig_umap_clusters, width = 8, height = 6, dpi = 300, compression = "lzw")

# Select approximately 500 compounds, balanced across clusters.
compounds_per_cluster <- ceiling(500 / k)
diverse_selection <- umap_dt_unique[, .SD[order(dist_to_centroid)][seq_len(min(.N, compounds_per_cluster))],
                                    by = cluster]
stopifnot(uniqueN(diverse_selection$cichem) == nrow(diverse_selection))
fwrite(diverse_selection[, .(cluster, cichem, primary_name, sector, UMAP1, UMAP2, dist_to_centroid)],
       file.path(output_dir, "diverse_ROCKS_endocrine_screening_selection_500_06sept2026.csv"))

# ---- Overlay candidates on full ROCKS densMAP ----
rocks_embed_file <- file.path("data", "raw", "chemical_fingerprint_densmap_8.18.2026.csv")
if (!file.exists(rocks_embed_file)) stop("densMAP file not found: ", rocks_embed_file)
rocks_embed <- fread(rocks_embed_file)

cichem_to_cicomp <- unique(rocks_chemistry[, .(cichem, cicomp_id = cicomp)])
screening_with_cicomp <- merge(
  umap_dt_unique[, .(cichem, cluster, dist_to_centroid, sector)], cichem_to_cicomp, by = "cichem")
screening_densmap <- merge(screening_with_cicomp,
                           rocks_embed[, .(cicomp, densmap_1, densmap_2)],
                           by.x = "cicomp_id", by.y = "cicomp")
screening_densmap[, selected_500 := cichem %in% diverse_selection$cichem]
screening_densmap[, is_medoid := cichem %in% medoid_reps$cichem]

# Prefer an existing cichem column in rocks_embed. Otherwise map cicomp to cichem.
if (!"cichem" %in% names(rocks_embed)) {
  rocks_embed <- merge(rocks_embed, cichem_to_cicomp, by.x = "cicomp", by.y = "cicomp_id", all.x = TRUE)
}
rocks_embed_sector <- merge(rocks_embed, sector_lookup, by = "cichem", all.x = TRUE)
medoid_densmap <- screening_densmap[is_medoid == TRUE]

# fig7 panel C
fig_densmap_candidates <- ggplot() +
  geom_point(data = rocks_embed_sector, aes(densmap_1, densmap_2, color = sector), alpha = 0.15, size = 0.2) +
  geom_point(data = screening_densmap, aes(densmap_1, densmap_2, color = sector), alpha = 0.4, size = 0.4) +
  geom_point(data = screening_densmap[selected_500 == TRUE], aes(densmap_1, densmap_2, fill = "500 diverse picks"),
             shape = 24, size = 3, stroke = 0) +
  geom_point(data = medoid_densmap, aes(densmap_1, densmap_2), color = "black", shape = 8, size = 4) +
  geom_label_repel(data = medoid_densmap, aes(densmap_1, densmap_2, label = cluster),
                   color = "black", fontface = "bold", size = 3, fill = "white", max.overlaps = Inf) +
  scale_color_manual(values = sector_colors) +
  scale_fill_manual(name = "Selection", values = c("500 diverse picks" = "#0A32FF")) +
  coord_cartesian(xlim = c(-8, 17), ylim = c(-1, 20)) +
  theme_minimal(base_size=14) +
  theme(panel.grid.minor = element_blank(), plot.margin = margin(20, 20, 20, 20)) +
  labs(x = "densMAP 1", y = "densMAP 2", color = "Sector")+
  guides(
  color = guide_legend(
    override.aes = list(size = 4, alpha = 1)
  ),
  fill  = guide_legend(
    override.aes = list(
      size = 4,
      shape = 24,
      fill="#0A32FF"
    )
  )
)
fig_densmap_candidates

fig_combined_abc <- arrangeGrob(
  arrangeGrob(fig_property_coverage, top = textGrob("A", x = 0, hjust = -0.5, gp = gpar(fontsize = 16, fontface = "bold"))),
  arrangeGrob(fig_upset, top = textGrob("B", x = 0, hjust = -0.5, gp = gpar(fontsize = 16, fontface = "bold"))),
  arrangeGrob(fig_densmap_candidates, top = textGrob("C", x = 0, hjust = -0.5, gp = gpar(fontsize = 16, fontface = "bold"))),
  layout_matrix = rbind(c(1, 2), c(3, 3)), widths = c(1, 1), heights = c(1, 1)
)
png(
  file.path(output_dir, "fig7_allpanels_06sept2026.png"),
  width = 13,
  height = 10,
  units = "in",
  res = 600
)
grid.draw(fig_combined_abc)
dev.off()

# ---- ER/AR activity and prediction coverage ----
collapse_mode <- function(x) {
  if (all(is.na(x))) NA_integer_ else as.integer(any(x == 1, na.rm = TRUE))
}
endo_unique <- combined_endo_summary_classy[model %in% c("CERAPP", "COMPARA"), .(
  binding = collapse_mode(binding), agonist = collapse_mode(agonist),
  antagonist = collapse_mode(antagonist), is_qsar_ready = first(is_qsar_ready),
  discard_error = first(discard_error), dtxsid = first(dtxsid), kingdom = first(kingdom),
  superclass = first(superclass), class = first(class), subclass = first(subclass)
), by = .(cichem_id, model)]

cerapp_compara_long <- melt(endo_unique,
  id.vars = c("cichem_id", "dtxsid", "model", "is_qsar_ready", "discard_error"),
  measure.vars = c("binding", "agonist", "antagonist"),
  variable.name = "Mode", value.name = "Prediction")
cerapp_compara_long[, Prediction := fcase(
  is.na(Prediction), "NA", Prediction == 0, "0", Prediction == 1, "1")]
cerapp_compara_long[, Prediction := factor(Prediction, levels = c("NA", "0", "1"))]

fig_in_silico_endo_coverage <- ggplot(cerapp_compara_long, aes(model, fill = Prediction)) +
  geom_bar() + facet_wrap(~Mode) +
  geom_text(stat = "count", aes(label = after_stat(count), color = Prediction),
            position = position_stack(vjust = 0.5), size = 3, fontface = "bold") +
  scale_fill_manual(values = c("NA" = "#E5DDCA", "0" = "#BCE4F7", "1" = "#0A32FF"),
                    labels = c("No Prediction", "Inactive", "Active"), name = "Prediction") +
  scale_color_manual(values = c("NA" = "black", "0" = "black", "1" = "white"), guide = "none") +
  labs(x = "Model", y = "Chemical Count") + theme_minimal()

toxcast_long <- melt(combined_toxcast_summary,
  id.vars = c("cichem_id", "dtxsid", "model"), measure.vars = c("agonist", "antagonist"),
  variable.name = "Mode", value.name = "AUC")
toxcast_long <- toxcast_long[!is.na(dtxsid) & !(model == "ToxCast_ER_2017" & Mode == "antagonist")]
toxcast_long[, Prediction := fcase(
  is.na(AUC), "NA", AUC > 0.1, "Active", AUC > 0.001, "Equivocal", default = "Inactive")]
toxcast_long[, Prediction := factor(Prediction, levels = c("Inactive", "Equivocal", "Active", "NA"))]

fig_toxcast_endo_coverage <- ggplot(toxcast_long[Prediction != "NA"], aes(model, fill = Prediction)) +
  geom_bar() + facet_wrap(~Mode, scales = "free_x") +
  geom_text(stat = "count", aes(label = after_stat(count), color = Prediction),
            position = position_stack(vjust = 0.5), size = 3, fontface = "bold") +
  scale_fill_manual(values = c("Inactive" = "#BCE4F7", "Equivocal" = ULcolors[13], "Active" = "#0A32FF"),
                    name = "Prediction") +
  scale_color_manual(values = c("Inactive" = "black", "Equivocal" = "black", "Active" = "white"), guide = "none") +
  labs(x = "Model", y = "Chemical Count") + theme_minimal()

png(file.path(output_dir, "fig_supplement_barchart_endo_coverage.png"),
    width = 8, height = 8, units = "in", res = 300)
grid.arrange(
  arrangeGrob(fig_in_silico_endo_coverage, top = textGrob("A", x = 0, hjust = -0.5, gp = gpar(fontsize = 16, fontface = "bold"))),
  arrangeGrob(fig_toxcast_endo_coverage, top = textGrob("B", x = 0, hjust = -0.5, gp = gpar(fontsize = 16, fontface = "bold"))),
  ncol = 1
)
dev.off()

# ---- ClassyFire heatmaps ----
heatmap_data <- combined_endo_summary_classy[
  model %in% c("CERAPP", "COMPARA") & !is.na(superclass) & superclass != "",
  .(any_positive = as.integer(any(c(agonist, antagonist, binding) == 1, na.rm = TRUE)),
    any_predicted = as.integer(any(!is.na(c(agonist, antagonist, binding))))),
  by = .(cichem_id, superclass, model)]
heatmap_summary <- heatmap_data[, .(
  n_chemicals = .N,
  percent_positive = 100 * sum(any_positive) / .N,
  percent_predicted = 100 * sum(any_predicted) / .N,
  count_positive = sum(any_positive), count_predicted = sum(any_predicted)
), by = .(superclass, model)]
heatmap_counts_long <- melt(heatmap_summary,
  id.vars = c("superclass", "model", "n_chemicals"),
  measure.vars = c("count_positive", "count_predicted"),
  variable.name = "metric", value.name = "count")
heatmap_summary_long <- melt(heatmap_summary,
  id.vars = c("superclass", "model", "n_chemicals"),
  measure.vars = c("percent_positive", "percent_predicted"),
  variable.name = "metric", value.name = "percent")

superclass_order_p <- heatmap_summary_long[metric == "percent_positive"][order(-percent), unique(superclass)]
heatmap_summary_long[, superclass := factor(superclass, levels = superclass_order_p)]
superclass_order_n <- heatmap_counts_long[metric == "count_positive"][order(-count), unique(superclass)]
heatmap_counts_long[, superclass := factor(superclass, levels = superclass_order_n)]

fig_heatmap_percent <- ggplot(heatmap_summary_long, aes(model, superclass, fill = percent)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(percent, 1), color = ifelse(percent > 50, "white", "black")), size = 4) +
  scale_fill_gradient(low = "white", high = "#0A32FF", limits = c(0, 100), name = "%") +
  scale_color_identity() +
  facet_wrap(~metric, labeller = labeller(metric = c(percent_positive = "% Positive", percent_predicted = "% Predicted"))) +
  labs(title = "Chemical Predictions by Superclass and Model", x = "Model", y = "Chemical Superclass") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12), plot.title = element_text(hjust = 0.5, size = 14),
        aspect.ratio = 4)

fig_heatmap_counts <- ggplot(heatmap_counts_long, aes(model, superclass, fill = count)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = count, color = ifelse(count > max(count) / 2, "white", "black")), size = 4) +
  scale_fill_gradient(low = "white", high = "#0A32FF", name = "Count") +
  scale_color_identity() +
  facet_wrap(~metric, labeller = labeller(metric = c(count_positive = "Count Positive", count_predicted = "Count Predicted"))) +
  labs(title = "Chemical Prediction Counts by Superclass and Model", x = "Model", y = "Chemical Superclass") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 10), plot.title = element_text(hjust = 0.5, size = 14),
        aspect.ratio = 4)

ggsave(file.path(output_dir, "supplement_classyfire_heatmap_percentage.png"),
       fig_heatmap_percent, width = 12, height = 20, dpi = 300)
ggsave(file.path(output_dir, "supplement_classyfire_heatmap_counts.png"),
       fig_heatmap_counts, width = 12, height = 20, dpi = 300)

message("Analysis complete. Outputs written to: ", normalizePath(output_dir))
