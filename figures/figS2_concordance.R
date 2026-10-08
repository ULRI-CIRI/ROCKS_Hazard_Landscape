
# =============================================================================
# Analyze concordance between HED p05 and AED p05
# =============================================================================

library(tidyverse)
library(here)
library(broom)
library(scales)

# -----------------------------------------------------------------------------
# 1. Load and assemble analysis data
# -----------------------------------------------------------------------------

message("\nLoading input datasets...")

rocks <- readRDS(here("data", "derived", "rocks.rds"))

aed_p <- readRDS(here::here("data", "derived", "b_mc7_aed_exec.rds")) %>%
  select(
    cichem_id_exec,
    n_aed_ac50_sumcl50,
    p5_aed_ac50_sumcl50,
    pct_assay_endpoints_active
  ) %>%
  rename(
    cichem_id = cichem_id_exec
  )

hed_p_by_chem_ebmdhed <- readRDS(
  here::here("data", "derived", "hed_p_by_chem_ebmdhed.rds")
) %>%
  select(
    cichem_id,
    hed_p05,
    n_hed_values
  ) %>%
  filter(
    !is.na(hed_p05),
    hed_p05 > 0,
    n_hed_values > 0
  )

p5_hed_aed <- rocks %>%
  distinct(cichem_id) %>%
  left_join(
    hed_p_by_chem_ebmdhed,
    by = "cichem_id"
  ) %>%
  left_join(
    aed_p,
    by = "cichem_id"
  )

# -----------------------------------------------------------------------------
# 2. Clean paired analysis dataset
# -----------------------------------------------------------------------------

message("\nPreparing paired HED-AED analysis dataset...")

p5_hed_aed_analysis <- p5_hed_aed %>%
  filter(
    !is.na(cichem_id),
    !is.na(hed_p05),
    !is.na(p5_aed_ac50_sumcl50),
    hed_p05 > 0,
    p5_aed_ac50_sumcl50 > 0,
    !is.na(n_aed_ac50_sumcl50),
    !is.na(pct_assay_endpoints_active),
    !is.na(n_hed_values)
  ) %>%
  mutate(
    log10_hed_p05 = log10(hed_p05),
    log10_aed_p05 = log10(p5_aed_ac50_sumcl50),

    # Difference on log10 scale.
    # Positive = HED p05 is greater than AED p05.
    # Negative = HED p05 is lower than AED p05.
    log10_diff_hed_minus_aed = log10_hed_p05 - log10_aed_p05,

    # Absolute discordance in log10 units.
    abs_log10_diff = abs(log10_diff_hed_minus_aed),

    # Fold difference, expressed symmetrically as the larger/smaller value.
    fold_difference = 10^abs_log10_diff,

    # Log-transform depth variables for modeling because these are positive,
    # right-skewed count-like variables.
    log10_n_aed_ac50 = log10(n_aed_ac50_sumcl50 + 1),
    log10_n_hed_values = log10(n_hed_values + 1),

    # Optional categorical bins for plots.
    aed_depth_bin = case_when(
      n_aed_ac50_sumcl50 <= quantile(n_aed_ac50_sumcl50, 1 / 3, na.rm = TRUE) ~ "Low AED endpoint depth",
      n_aed_ac50_sumcl50 <= quantile(n_aed_ac50_sumcl50, 2 / 3, na.rm = TRUE) ~ "Medium AED endpoint depth",
      TRUE ~ "High AED endpoint depth"
    ),

    active_pct_bin = case_when(
      pct_assay_endpoints_active <= quantile(pct_assay_endpoints_active, 1 / 3, na.rm = TRUE) ~ "Low % active endpoints",
      pct_assay_endpoints_active <= quantile(pct_assay_endpoints_active, 2 / 3, na.rm = TRUE) ~ "Medium % active endpoints",
      TRUE ~ "High % active endpoints"
    ),

    hed_depth_bin = case_when(
      n_hed_values <= quantile(n_hed_values, 1 / 3, na.rm = TRUE) ~ "Low HED evidence depth",
      n_hed_values <= quantile(n_hed_values, 2 / 3, na.rm = TRUE) ~ "Medium HED evidence depth",
      TRUE ~ "High HED evidence depth"
    )
  ) %>%
  mutate(
    aed_depth_bin = factor(
      aed_depth_bin,
      levels = c(
        "Low AED endpoint depth",
        "Medium AED endpoint depth",
        "High AED endpoint depth"
      )
    ),
    active_pct_bin = factor(
      active_pct_bin,
      levels = c(
        "Low % active endpoints",
        "Medium % active endpoints",
        "High % active endpoints"
      )
    ),
    hed_depth_bin = factor(
      hed_depth_bin,
      levels = c(
        "Low HED evidence depth",
        "Medium HED evidence depth",
        "High HED evidence depth"
      )
    )
  )

message("Number of paired chemicals: ", nrow(p5_hed_aed_analysis))

saveRDS(
  p5_hed_aed_analysis,
  here("data", "derived", "p5_hed_aed_concordance_dat.rds")
)



# =============================================================================
# Threshold analysis:
# How does HED-AED concordance change with increasing in vitro signal/coverage?
#
# Assumes p5_hed_aed_analysis already contains:
#   cichem_id
#   hed_p05
#   p5_aed_ac50_sumcl50
#   n_aed_ac50_sumcl50
#   pct_assay_endpoints_active
#   log10_hed_p05
#   log10_aed_p05
#
# If log10 columns do not exist, they are created below.
# =============================================================================

library(tidyverse)
library(here)
library(scales)

# -----------------------------------------------------------------------------
# 1. Prepare analysis dataset
# -----------------------------------------------------------------------------

hed_aed_dat <- p5_hed_aed_analysis %>%
  filter(
    !is.na(hed_p05),
    !is.na(p5_aed_ac50_sumcl50),
    hed_p05 > 0,
    p5_aed_ac50_sumcl50 > 0,
    !is.na(n_aed_ac50_sumcl50),
    !is.na(pct_assay_endpoints_active)
  )

message("Number of paired chemicals: ", nrow(hed_aed_dat))

# -----------------------------------------------------------------------------
# 2. Helper function for threshold-based concordance
# -----------------------------------------------------------------------------

calc_threshold_concordance <- function(dat, threshold_var, thresholds, label) {

  map_dfr(thresholds, function(thresh) {

    dat_thresh <- dat %>%
      filter(.data[[threshold_var]] >= thresh)

    n_thresh <- nrow(dat_thresh)

    if (n_thresh < 5) {
      return(
        tibble(
          threshold_type = label,
          threshold = thresh,
          n = n_thresh,
          pearson_r = NA_real_,
          pearson_r2 = NA_real_,
          spearman_rho = NA_real_,
          rmse_log10 = NA_real_,
          median_abs_log10_diff = NA_real_
        )
      )
    }

    pearson_r <- cor(
      dat_thresh$log10_aed_p05,
      dat_thresh$log10_hed_p05,
      method = "pearson",
      use = "complete.obs"
    )

    spearman_rho <- cor(
      dat_thresh$log10_aed_p05,
      dat_thresh$log10_hed_p05,
      method = "spearman",
      use = "complete.obs"
    )

    lm_fit <- lm(
      log10_hed_p05 ~ log10_aed_p05,
      data = dat_thresh
    )

    pred <- predict(lm_fit)

    tibble(
      threshold_type = label,
      threshold = thresh,
      n = n_thresh,
      pearson_r = pearson_r,
      pearson_r2 = pearson_r^2,
      spearman_rho = spearman_rho,
      rmse_log10 = sqrt(mean((dat_thresh$log10_hed_p05 - pred)^2, na.rm = TRUE)),
      median_abs_log10_diff = median(
        abs(dat_thresh$log10_hed_p05 - dat_thresh$log10_aed_p05),
        na.rm = TRUE
      )
    )
  })
}

# -----------------------------------------------------------------------------
# 3. Define thresholds
# -----------------------------------------------------------------------------

# Use interpretable thresholds rather than every unique value.
# Adjust these if your distributions suggest better cut points.

endpoint_thresholds <- c(
  0, 5, 10, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250
)

endpoint_thresholds <- endpoint_thresholds[
  endpoint_thresholds <= max(hed_aed_dat$n_aed_ac50_sumcl50, na.rm = TRUE)
]

pct_active_thresholds <- c(
  0, 1, 2.5, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50
)

pct_active_thresholds <- pct_active_thresholds[
  pct_active_thresholds <= max(hed_aed_dat$pct_assay_endpoints_active, na.rm = TRUE)
]

# -----------------------------------------------------------------------------
# 4. Calculate threshold-based metrics
# -----------------------------------------------------------------------------

threshold_concordance <- bind_rows(
  calc_threshold_concordance(
    dat = hed_aed_dat,
    threshold_var = "n_aed_ac50_sumcl50",
    thresholds = endpoint_thresholds,
    label = "Minimum AED endpoint count"
  ),
  calc_threshold_concordance(
    dat = hed_aed_dat,
    threshold_var = "pct_assay_endpoints_active",
    thresholds = pct_active_thresholds,
    label = "Minimum % active endpoints"
  )
)

print(threshold_concordance)

saveRDS(
  threshold_concordance,
  here("data", "derived", "hed_aed_threshold_concordance.rds")
)

# -----------------------------------------------------------------------------
# 5. Simple plot: R2 and Spearman rho by threshold
# -----------------------------------------------------------------------------

plot_dat <- threshold_concordance %>%
  pivot_longer(
    cols = c(pearson_r2, spearman_rho),
    names_to = "metric",
    values_to = "value"
  ) %>%
mutate(
  metric = recode(
    metric,
    pearson_r2 = "Pearson R²",
    spearman_rho = "Spearman rho"
  )
)

p_threshold_concordance <- ggplot(
  plot_dat,
  aes(
    x = threshold,
    y = value,
    color = metric,
    group = metric
  )
) +
  geom_line(
    linewidth = 0.9
  ) +
  geom_point(
    aes(size = n),
    alpha = 0.85
  ) +
  geom_text(
    data = plot_dat %>%
      filter(metric == expression("Pearson " * R^2)),
    aes(
      label = paste0("N=", n)
    ),
    color = "grey25",
    size = 3,
    vjust = -1,
    show.legend = FALSE
  ) +
  facet_wrap(
    ~ threshold_type,
    scales = "free_x"
  ) +
  scale_color_manual(
    values = c(
      "Pearson R²" = "#122C49",
      "Spearman rho" = "#0A32FF"
    ),
    labels = c(
      "Pearson R²",
      "Spearman rho"
    )
  ) +
  scale_size_continuous(
    name = "Number of chemicals",
    range = c(2, 6)
  ) +
  coord_cartesian(
    ylim = c(0, 1)
  ) +
  labs(
    x = "Minimum threshold",
    y = "Correlation / concordance metric",
    color = NULL,
    title = "HED-AED concordance by in vitro signal and assay coverage",
    subtitle = "Each point restricts the analysis to chemicals meeting or exceeding the threshold"
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey85", color = "grey40"),
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )

p_threshold_concordance

library(ggrepel)

r2plot <- ggplot(
  threshold_concordance,
  aes(
    x = threshold,
    y = pearson_r2
  )
) +
  geom_line(
    color = "#122C49",
    linewidth = 1
  ) +
  geom_point(
    aes(size = n),
    color = "#0A32FF"
  ) +
  geom_text_repel(
    aes(label = paste0("N=", n)),
    box.padding = 0.4,
    point.padding = 0.3,
    min.segment.length = 0,
    max.overlaps = Inf,
    seed = 123
  ) +
  facet_wrap(
    ~ threshold_type,
    scales = "free_x"
  ) +
  labs(
    x = "Minimum threshold",
    y = expression(R^2)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(size = 16, face = "bold")
  )


r2plot

ggsave(
  here("outputs","figures", "hed_aed_threshold_concordance.png"),
  r2plot,
  width = 11,
  height = 7,
  dpi = 300
)

high_conf <- p5_hed_aed_analysis %>%
  filter(
    pct_assay_endpoints_active >= 25,
    n_aed_ac50_sumcl50 >= 150
  )

nrow(high_conf)

m_high_conf <- lm(
  log10(hed_p05) ~ log10(p5_aed_ac50_sumcl50),
  data = high_conf
)

summary(m_high_conf)
summary(m_high_conf)$r.squared


cor.test(
  log10(high_conf$hed_p05),
  log10(high_conf$p5_aed_ac50_sumcl50),
  method = "pearson"
)

cor.test(
  log10(high_conf$hed_p05),
  log10(high_conf$p5_aed_ac50_sumcl50),
  method = "spearman"
)

# all for comparison

nrow(p5_hed_aed_analysis)

m_all <- lm(
  log10(hed_p05) ~ log10(p5_aed_ac50_sumcl50),
  data = p5_hed_aed_analysis
)

summary(m_all)
summary(m_all)$r.squared


cor.test(
  log10(p5_hed_aed_analysis$hed_p05),
  log10(p5_hed_aed_analysis$p5_aed_ac50_sumcl50),
  method = "pearson"
)

cor.test(
  log10(p5_hed_aed_analysis$hed_p05),
  log10(p5_hed_aed_analysis$p5_aed_ac50_sumcl50),
  method = "spearman"
)
