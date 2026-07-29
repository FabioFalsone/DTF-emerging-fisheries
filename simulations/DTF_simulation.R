# =====================================================================================
# DTF simulation script aligned with the final manuscript
# =====================================================================================
# Purpose
# -------
# This script reproduces the simulation-based validation and sensitivity analyses
# described in Section 2.6 of the manuscript:
#
# "A Decision-Tree Framework for the Sustainable Management of Emerging Fisheries
# and Fishing Innovations" by Falsone et al.
#
# The script includes:
#   - a two-phase workflow comprising calibration and full-factorial performance
#     evaluation;
#   - four composite-score weighting schemes;
#   - the automated landing-based screening procedure implemented in the Shiny
#     application;
#   - Supplementary Figures S1-S5; and
#   - main-text Figure 2, including panels A-C.
#
# The workflow, analytical settings and outputs are organised to match the final
# manuscript and its supplementary material.
# =====================================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
})

# -------------------------------- SETTINGS -------------------------------------------

set.seed(123)

# Output directory
OUT_DIR <- getwd()
FIG_DIR <- file.path(OUT_DIR, "figures")
TAB_DIR <- file.path(OUT_DIR, "tables")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)

# Simulation design exactly as described in the manuscript
SCENARIOS <- c(
  "no_change",
  "abrupt_moderate",
  "abrupt_strong",
  "gradual_expansion",
  "expansion_peak_collapse"
)

N_YEARS_GRID  <- c(8, 10, 15, 20, 30)
PHI_GRID      <- c(0, 0.3, 0.6)
CV_GRID       <- c(0.10, 0.25)
EFFECT_GRID   <- c(0.40, 0.80, 1.25)
BREAK_POS_GRID <- c(0.30, 0.50, 0.70)

# Phase 1: calibration
CALIBRATION_N <- 20000

# Phase 2: full factorial performance evaluation
N_REPS_FULL <- 100

# Breakpoint-testing parameters
MIN_POINTS <- 3
ALPHA <- 0.05

# Calibration utility settings
UTILITY_LAMBDA <- 2  # useful detection within +/-2 years - 2 * false positive rate
FPR_CONSTRAINTS <- c(0.05, 0.10, 0.15)
SELECTED_FPR_CAP <- 0.10

# Threshold grid screened in phase 1
THRESHOLD_GRID <- expand.grid(
  min_percent_increase = c(25, 50, 75, 100),
  min_abs_d            = c(0.5, 0.8, 1.2),
  require_p            = c(FALSE, TRUE),
  min_post_pos_diffs   = c(0, 1, 2),
  min_score            = c(0, 0.25, 0.50, 1.00),
  stringsAsFactors = FALSE
) %>%
  mutate(threshold_id = dplyr::row_number())

# Four weighting schemes exactly as reported in the final paper
WEIGHT_GRID <- tibble::tribble(
  ~weight_id,            ~w_p, ~w_d, ~w_inc, ~weight_label,
  "equal",               1.0,  1.0,   1.0,   "Equal",
  "effect_priority",     0.5,  1.5,   1.0,   "Effect-priority",
  "increase_priority",   0.5,  1.0,   1.5,   "Increase-priority",
  "no_p",                0.0,  1.0,   1.0,   "No-p-value"
)

# -------------------------------- HELPERS --------------------------------------------

safe_p <- function(x) {
  ifelse(is.na(x) | !is.finite(x), 1, pmin(pmax(x, 0), 1))
}

safe_var <- function(x) {
  if (length(x) < 2) return(NA_real_)
  stats::var(x, na.rm = TRUE)
}

cohens_d <- function(g1, g2) {
  n1 <- length(g1); n2 <- length(g2)
  v1 <- safe_var(g1); v2 <- safe_var(g2)
  if (!is.finite(v1) || !is.finite(v2) || n1 < 2 || n2 < 2) return(NA_real_)
  pooled <- sqrt(((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2))
  if (!is.finite(pooled) || pooled == 0) return(NA_real_)
  (mean(g2) - mean(g1)) / pooled
}

safe_shapiro_p <- function(x) {
  if (length(x) < 3 || length(x) > 5000 || length(unique(x)) < 3) return(NA_real_)
  tryCatch(stats::shapiro.test(x)$p.value, error = function(e) NA_real_)
}

safe_fligner_p <- function(g1, g2) {
  if (length(g1) < 2 || length(g2) < 2) return(NA_real_)
  tryCatch(stats::fligner.test(list(g1, g2))$p.value, error = function(e) NA_real_)
}

safe_ljung_p <- function(x, lag = NULL) {
  if (length(x) < 8) return(NA_real_)
  if (is.null(lag)) lag <- min(3, length(x) - 2)
  tryCatch(stats::Box.test(x, lag = lag, type = "Ljung-Box")$p.value,
           error = function(e) NA_real_)
}

lag1_acf <- function(x) {
  if (length(x) < 4 || stats::sd(x) == 0) return(NA_real_)
  tryCatch(as.numeric(stats::acf(x, plot = FALSE, lag.max = 1)$acf[2]),
           error = function(e) NA_real_)
}

ar1_errors <- function(n, phi = 0, sd = 1) {
  e <- numeric(n)
  e[1] <- rnorm(1, 0, sd)
  if (n > 1) {
    innov_sd <- ifelse(abs(phi) < 1, sd * sqrt(1 - phi^2), sd)
    for (i in 2:n) e[i] <- phi * e[i - 1] + rnorm(1, 0, innov_sd)
  }
  e
}

wilson_ci <- function(x, n, conf = 0.95) {
  if (is.na(x) || is.na(n) || n == 0) return(c(NA_real_, NA_real_))
  z <- qnorm(1 - (1 - conf) / 2)
  p <- x / n
  denom <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denom
  half <- (z / denom) * sqrt((p * (1 - p) / n) + (z^2 / (4 * n^2)))
  c(max(0, centre - half), min(1, centre + half))
}

# -------------------------------- SIMULATION MODEL -----------------------------------

simulate_fishery_series <- function(n_years,
                                    scenario,
                                    phi = 0,
                                    cv = 0.10,
                                    effect = 0.80,
                                    break_pos = 0.50,
                                    base = 100,
                                    start_year = 2000) {
  years <- start_year + seq_len(n_years) - 1
  t <- seq_len(n_years)

  true_idx <- round(n_years * break_pos)
  true_idx <- max(MIN_POINTS + 1, min(true_idx, n_years - MIN_POINTS))

  mu <- rep(base, n_years)
  has_change <- scenario != "no_change"

  if (scenario == "no_change") {
    true_idx <- NA_integer_
    mu <- rep(base, n_years)
  }

  if (scenario == "abrupt_moderate") {
    mu[(true_idx + 1):n_years] <- base * (1 + effect)
  }

  if (scenario == "abrupt_strong") {
    strong_effect <- max(effect, 1.25)
    mu[(true_idx + 1):n_years] <- base * (1 + strong_effect)
  }

  if (scenario == "gradual_expansion") {
    after <- (true_idx + 1):n_years
    if (length(after) > 0) {
      ramp <- seq(0, effect, length.out = length(after))
      mu[after] <- base * (1 + ramp)
    }
  }

  if (scenario == "expansion_peak_collapse") {
    after <- (true_idx + 1):n_years
    if (length(after) > 0) {
      peak_effect <- max(effect, 0.80)
      peak_len <- max(2, ceiling(length(after) * 0.45))
      rise <- seq(0, peak_effect, length.out = peak_len)
      fall <- seq(peak_effect * 0.90, -0.20, length.out = length(after) - peak_len)
      path <- c(rise, fall)
      mu[after] <- base * (1 + path)
      mu <- pmax(mu, base * 0.20)
    }
  }

  eps <- ar1_errors(n_years, phi = phi, sd = cv)
  landing <- mu * exp(eps - 0.5 * cv^2)

  data.frame(
    year = years,
    t = t,
    landing = landing,
    true_idx = true_idx,
    true_year = ifelse(is.na(true_idx), NA_integer_, years[true_idx]),
    has_change = has_change,
    scenario = scenario,
    phi = phi,
    cv = cv,
    effect = effect,
    break_pos = break_pos,
    n_years = n_years,
    stringsAsFactors = FALSE
  )
}

# ------------------------------ DTF SCREENING LOGIC ----------------------------------

candidate_table <- function(df, min_points = MIN_POINTS, alpha = ALPHA) {
  n <- nrow(df)
  if (n < min_points * 2) return(tibble())

  out <- vector("list", n - 2 * min_points + 1)
  counter <- 1

  for (i in min_points:(n - min_points)) {
    g1 <- df$landing[1:i]
    g2 <- df$landing[(i + 1):n]

    m1 <- mean(g1)
    m2 <- mean(g2)
    inc <- ifelse(m1 > 0, 100 * (m2 - m1) / m1, NA_real_)
    d <- cohens_d(g1, g2)

    p_t <- tryCatch(stats::t.test(g2, g1, alternative = "greater", var.equal = TRUE)$p.value,
                    error = function(e) NA_real_)
    p_w <- tryCatch(stats::wilcox.test(g2, g1, alternative = "greater", exact = FALSE)$p.value,
                    error = function(e) NA_real_)
    p_fligner <- safe_fligner_p(g1, g2)
    p_shap1 <- safe_shapiro_p(g1)
    p_shap2 <- safe_shapiro_p(g2)

    normal_ok <- (!is.na(p_shap1) && !is.na(p_shap2) && p_shap1 >= alpha && p_shap2 >= alpha)
    var_ok <- (!is.na(p_fligner) && p_fligner >= alpha)
    test_used <- ifelse(normal_ok && var_ok, "t-test", "Wilcoxon")
    p_adapt <- ifelse(test_used == "t-test", p_t, p_w)

    after_values <- df$landing[i:n]
    post_diffs <- diff(after_values)
    post_positive_diffs <- ifelse(length(post_diffs) == 0, 0, sum(post_diffs > 0, na.rm = TRUE))

    out[[counter]] <- tibble(
      idx = i,
      year = df$year[i],
      n_before = length(g1),
      n_after = length(g2),
      mean_before = m1,
      mean_after = m2,
      percent_increase = inc,
      cohens_d = d,
      p_t = safe_p(p_t),
      p_wilcox = safe_p(p_w),
      p_adaptive = safe_p(p_adapt),
      test_used = test_used,
      p_fligner = p_fligner,
      p_shapiro_before = p_shap1,
      p_shapiro_after = p_shap2,
      normal_ok = normal_ok,
      var_ok = var_ok,
      post_positive_diffs = post_positive_diffs
    )
    counter <- counter + 1
  }

  bind_rows(out)
}

add_scores <- function(candidates, weight_row) {
  if (nrow(candidates) == 0) return(candidates)
  candidates %>%
    mutate(
      inc_prop = pmax(percent_increase / 100, 0),
      effect_component = pmax(abs(cohens_d), 0),
      p_component = pmax(1 - p_adaptive, 0),
      score = (p_component ^ weight_row$w_p) *
              (effect_component ^ weight_row$w_d) *
              (inc_prop ^ weight_row$w_inc),
      weight_id = weight_row$weight_id,
      weight_label = weight_row$weight_label
    )
}

apply_dtf_threshold <- function(candidates, threshold_row, alpha = ALPHA) {
  if (nrow(candidates) == 0) return(NULL)

  pass <- candidates %>%
    filter(
      is.finite(percent_increase),
      is.finite(cohens_d),
      is.finite(score),
      percent_increase >= threshold_row$min_percent_increase,
      abs(cohens_d) >= threshold_row$min_abs_d,
      post_positive_diffs >= threshold_row$min_post_pos_diffs,
      score >= threshold_row$min_score
    )

  if (isTRUE(threshold_row$require_p)) {
    pass <- pass %>% filter(p_adaptive <= alpha)
  }

  if (nrow(pass) == 0) return(NULL)
  pass %>% arrange(desc(score)) %>% slice(1)
}

calc_metrics <- function(df) {
  tp <- sum(df$tp, na.rm = TRUE)
  fp <- sum(df$fp, na.rm = TRUE)
  tn <- sum(df$tn, na.rm = TRUE)
  fn <- sum(df$fn, na.rm = TRUE)

  precision <- ifelse(tp + fp > 0, tp / (tp + fp), NA_real_)
  recall <- ifelse(tp + fn > 0, tp / (tp + fn), NA_real_)
  f1 <- ifelse(is.finite(precision + recall) && (precision + recall) > 0,
               2 * precision * recall / (precision + recall), NA_real_)

  tibble(
    n = nrow(df),
    TP = tp, FP = fp, TN = tn, FN = fn,
    detection_rate = recall,
    false_positive_rate = ifelse(fp + tn > 0, fp / (fp + tn), NA_real_),
    false_negative_rate = ifelse(fn + tp > 0, fn / (fn + tp), NA_real_),
    precision = precision,
    recall = recall,
    F1 = f1,
    mean_error_years = mean(df$error_years[df$tp], na.rm = TRUE),
    MAE_years = mean(df$abs_error[df$tp], na.rm = TRUE),
    RMSE_years = sqrt(mean(df$error_years[df$tp]^2, na.rm = TRUE)),
    within_1_rate = mean(df$within_1[df$has_change], na.rm = TRUE),
    within_2_rate = mean(df$within_2[df$has_change], na.rm = TRUE),
    autocorr_warning_rate = mean(df$series_autocorr_warning, na.rm = TRUE),
    wilcoxon_selected_rate = mean(df$selected_test_used == "Wilcoxon", na.rm = TRUE)
  )
}

run_one_series_all_settings <- function(sim_id, df) {
  candidates_base <- candidate_table(df)
  series_ljung <- safe_ljung_p(df$landing)
  series_acf1  <- lag1_acf(df$landing)

  out <- list()
  k <- 1
  for (w_i in seq_len(nrow(WEIGHT_GRID))) {
    w <- WEIGHT_GRID[w_i, ]
    candidates <- add_scores(candidates_base, w)

    for (thr_i in seq_len(nrow(THRESHOLD_GRID))) {
      thr <- THRESHOLD_GRID[thr_i, ]
      bp <- apply_dtf_threshold(candidates, thr)

      detected <- !is.null(bp)
      estimated_idx <- if (detected) bp$idx else NA_integer_
      estimated_year <- if (detected) bp$year else NA_integer_
      has_change <- unique(df$has_change)
      true_idx <- unique(df$true_idx)
      true_year <- unique(df$true_year)

      tp <- detected && has_change
      fp <- detected && !has_change
      tn <- !detected && !has_change
      fn <- !detected && has_change

      error_years <- if (tp) estimated_idx - true_idx else NA_real_
      abs_error   <- if (tp) abs(error_years) else NA_real_
      within_1    <- if (tp) abs_error <= 1 else FALSE
      within_2    <- if (tp) abs_error <= 2 else FALSE

      out[[k]] <- tibble(
        sim_id = sim_id,
        scenario = unique(df$scenario),
        n_years = unique(df$n_years),
        phi = unique(df$phi),
        cv = unique(df$cv),
        effect = unique(df$effect),
        break_pos = unique(df$break_pos),
        has_change = has_change,
        true_idx = true_idx,
        true_year = true_year,
        threshold_id = thr$threshold_id,
        min_percent_increase = thr$min_percent_increase,
        min_abs_d = thr$min_abs_d,
        require_p = thr$require_p,
        min_post_pos_diffs = thr$min_post_pos_diffs,
        min_score = thr$min_score,
        weight_id = w$weight_id,
        weight_label = w$weight_label,
        detected = detected,
        estimated_idx = estimated_idx,
        estimated_year = estimated_year,
        tp = tp, fp = fp, tn = tn, fn = fn,
        error_years = error_years,
        abs_error = abs_error,
        within_1 = within_1,
        within_2 = within_2,
        series_ljung_p = series_ljung,
        series_acf1 = series_acf1,
        series_autocorr_warning = !is.na(series_ljung) && series_ljung < ALPHA,
        selected_score = if (detected) bp$score else NA_real_,
        selected_percent_increase = if (detected) bp$percent_increase else NA_real_,
        selected_cohens_d = if (detected) bp$cohens_d else NA_real_,
        selected_p_adaptive = if (detected) bp$p_adaptive else NA_real_,
        selected_test_used = if (detected) bp$test_used else NA_character_,
        selected_post_positive_diffs = if (detected) bp$post_positive_diffs else NA_integer_
      )
      k <- k + 1
    }
  }
  bind_rows(out)
}

run_one_series_selected_setting <- function(sim_id, df, selected_threshold, selected_weight) {
  candidates_base <- candidate_table(df)
  series_ljung <- safe_ljung_p(df$landing)
  series_acf1  <- lag1_acf(df$landing)

  candidates <- add_scores(candidates_base, selected_weight)
  bp <- apply_dtf_threshold(candidates, selected_threshold)

  detected <- !is.null(bp)
  estimated_idx <- if (detected) bp$idx else NA_integer_
  estimated_year <- if (detected) bp$year else NA_integer_
  has_change <- unique(df$has_change)
  true_idx <- unique(df$true_idx)
  true_year <- unique(df$true_year)

  tp <- detected && has_change
  fp <- detected && !has_change
  tn <- !detected && !has_change
  fn <- !detected && has_change

  error_years <- if (tp) estimated_idx - true_idx else NA_real_
  abs_error   <- if (tp) abs(error_years) else NA_real_
  within_1    <- if (tp) abs_error <= 1 else FALSE
  within_2    <- if (tp) abs_error <= 2 else FALSE

  tibble(
    sim_id = sim_id,
    scenario = unique(df$scenario),
    n_years = unique(df$n_years),
    phi = unique(df$phi),
    cv = unique(df$cv),
    effect = unique(df$effect),
    break_pos = unique(df$break_pos),
    has_change = has_change,
    true_idx = true_idx,
    true_year = true_year,
    detected = detected,
    estimated_idx = estimated_idx,
    estimated_year = estimated_year,
    tp = tp, fp = fp, tn = tn, fn = fn,
    error_years = error_years,
    abs_error = abs_error,
    within_1 = within_1,
    within_2 = within_2,
    series_ljung_p = series_ljung,
    series_acf1 = series_acf1,
    series_autocorr_warning = !is.na(series_ljung) && series_ljung < ALPHA,
    selected_score = if (detected) bp$score else NA_real_,
    selected_percent_increase = if (detected) bp$percent_increase else NA_real_,
    selected_cohens_d = if (detected) bp$cohens_d else NA_real_,
    selected_p_adaptive = if (detected) bp$p_adaptive else NA_real_,
    selected_test_used = if (detected) bp$test_used else NA_character_,
    selected_post_positive_diffs = if (detected) bp$post_positive_diffs else NA_integer_,
    threshold_id = selected_threshold$threshold_id,
    weight_id = selected_weight$weight_id,
    weight_label = selected_weight$weight_label
  )
}

# ------------------------------ PHASE 1: CALIBRATION ---------------------------------

message("Phase 1/2 - building calibration sample (", CALIBRATION_N, " series)...")

calibration_design <- tibble(
  sim_id = seq_len(CALIBRATION_N),
  scenario = sample(SCENARIOS, CALIBRATION_N, replace = TRUE),
  n_years = sample(N_YEARS_GRID, CALIBRATION_N, replace = TRUE),
  phi = sample(PHI_GRID, CALIBRATION_N, replace = TRUE),
  cv = sample(CV_GRID, CALIBRATION_N, replace = TRUE),
  effect = sample(EFFECT_GRID, CALIBRATION_N, replace = TRUE),
  break_pos = sample(BREAK_POS_GRID, CALIBRATION_N, replace = TRUE)
)

calibration_results_list <- vector("list", CALIBRATION_N)
for (i in seq_len(CALIBRATION_N)) {
  if (i %% 250 == 0) message("  Calibration series ", i, " / ", CALIBRATION_N)
  row <- calibration_design[i, ]
  df <- simulate_fishery_series(
    n_years = row$n_years,
    scenario = row$scenario,
    phi = row$phi,
    cv = row$cv,
    effect = row$effect,
    break_pos = row$break_pos
  )
  calibration_results_list[[i]] <- run_one_series_all_settings(row$sim_id, df)
}

calibration_results <- bind_rows(calibration_results_list)
write_csv(calibration_results, file.path(TAB_DIR, "calibration_results_all.csv"))

summary_by_threshold <- calibration_results %>%
  group_by(threshold_id, weight_id, weight_label,
           min_percent_increase, min_abs_d, require_p,
           min_post_pos_diffs, min_score) %>%
  group_modify(~calc_metrics(.x)) %>%
  ungroup() %>%
  mutate(utility = within_2_rate - UTILITY_LAMBDA * false_positive_rate)

write_csv(summary_by_threshold, file.path(TAB_DIR, "calibration_summary_by_threshold.csv"))

weight_sensitivity <- summary_by_threshold %>%
  group_by(weight_id, weight_label) %>%
  summarise(
    median_utility = median(utility, na.rm = TRUE),
    median_detection = median(within_2_rate, na.rm = TRUE),
    median_fpr = median(false_positive_rate, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(weight_sensitivity, file.path(TAB_DIR, "weight_sensitivity_summary.csv"))

recommended_settings <- bind_rows(lapply(FPR_CONSTRAINTS, function(fpr_cap) {
  summary_by_threshold %>%
    filter(false_positive_rate <= fpr_cap) %>%
    arrange(desc(utility), desc(F1), MAE_years) %>%
    slice(1) %>%
    mutate(fpr_cap = fpr_cap)
}))

write_csv(recommended_settings, file.path(TAB_DIR, "recommended_settings_by_fpr_cap.csv"))

selected_row <- recommended_settings %>%
  filter(abs(fpr_cap - SELECTED_FPR_CAP) < 1e-12) %>%
  slice(1)

selected_threshold <- THRESHOLD_GRID %>% filter(threshold_id == selected_row$threshold_id)
selected_weight    <- WEIGHT_GRID %>% filter(weight_id == selected_row$weight_id)

message("Selected setting for phase 2:")
print(selected_row)

# ------------------------------ PHASE 2: FULL FACTORIAL ------------------------------

message("Phase 2/2 - full factorial performance evaluation...")

performance_design <- expand.grid(
  rep = seq_len(N_REPS_FULL),
  scenario = SCENARIOS,
  n_years = N_YEARS_GRID,
  phi = PHI_GRID,
  cv = CV_GRID,
  effect = EFFECT_GRID,
  break_pos = BREAK_POS_GRID,
  stringsAsFactors = FALSE
) %>%
  mutate(sim_id = row_number())

message("Total full-factorial series: ", nrow(performance_design))

performance_results_list <- vector("list", nrow(performance_design))
for (i in seq_len(nrow(performance_design))) {
  if (i %% 500 == 0) message("  Performance series ", i, " / ", nrow(performance_design))
  row <- performance_design[i, ]
  df <- simulate_fishery_series(
    n_years = row$n_years,
    scenario = row$scenario,
    phi = row$phi,
    cv = row$cv,
    effect = row$effect,
    break_pos = row$break_pos
  )
  performance_results_list[[i]] <- run_one_series_selected_setting(row$sim_id, df, selected_threshold, selected_weight)
}

performance_results <- bind_rows(performance_results_list)
write_csv(performance_results, file.path(TAB_DIR, "performance_results_selected_setting.csv"))

# -------------------------------- TABLE SUMMARIES ------------------------------------

summary_by_condition <- performance_results %>%
  group_by(scenario, n_years, phi, cv, effect, break_pos) %>%
  group_modify(~calc_metrics(.x)) %>%
  ungroup()

write_csv(summary_by_condition, file.path(TAB_DIR, "performance_summary_by_condition.csv"))

summary_by_scenario <- performance_results %>%
  filter(has_change) %>%
  group_by(scenario) %>%
  summarise(
    useful_detection = mean(within_2, na.rm = TRUE),
    detection_rate = mean(detected, na.rm = TRUE),
    mean_abs_error = mean(abs_error[tp], na.rm = TRUE),
    .groups = "drop"
  )
write_csv(summary_by_scenario, file.path(TAB_DIR, "performance_summary_by_scenario.csv"))

# -------------------------------- SUPPLEMENTARY FIGURES ------------------------------

# S1. Illustrative trajectories
set.seed(321)
illustrative_df <- bind_rows(lapply(c(0.10, 0.25), function(cv_i) {
  bind_rows(lapply(c("no_change", "abrupt_moderate", "abrupt_strong", "expansion_peak_collapse", "gradual_expansion"), function(sc) {
    bind_rows(lapply(c(0, 0.3, 0.6), function(phi_i) {
      simulate_fishery_series(
        n_years = 20,
        scenario = sc,
        phi = phi_i,
        cv = cv_i,
        effect = 0.80,
        break_pos = 0.50
      )
    }))
  }))
})) %>%
  mutate(
    scenario_label = recode(scenario,
      no_change = "No change",
      abrupt_moderate = "Abrupt moderate",
      abrupt_strong = "Abrupt strong",
      expansion_peak_collapse = "Expansion-peak-collapse",
      gradual_expansion = "Gradual expansion"
    ),
    cv_label = paste0("CV = ", cv),
    phi_label = factor(phi, levels = c(0, 0.3, 0.6))
  )

p_s1 <- ggplot(illustrative_df, aes(x = year, y = landing, colour = phi_label, group = interaction(phi_label, cv, scenario))) +
  geom_line(linewidth = 0.8) +
  geom_vline(aes(xintercept = true_year), linetype = "dashed", colour = "black", data = distinct(illustrative_df, scenario, cv, phi, true_year, scenario_label, cv_label) %>% filter(!is.na(true_year))) +
  facet_grid(scenario_label ~ cv_label, scales = "free_y", switch = "y") +
  scale_colour_manual(values = c("0" = "#2C7FB8", "0.3" = "#E69F00", "0.6" = "#1FA187"), name = "AR(1) phi") +
  labs(
    title = "Illustrative simulated landing trajectories",
    subtitle = "Rows show scenarios, columns show observation noise, colours show AR(1) autocorrelation",
    x = "Year",
    y = "Simulated landings"
  ) +
  theme_bw(base_size = 13) +
  theme(strip.background = element_rect(fill = "grey90"),
        legend.position = "bottom")

ggsave(file.path(FIG_DIR, "Supplementary_Figure_S1.png"), p_s1, width = 14, height = 14, dpi = 300)

# S2. Useful detection heatmap (averaged over CV and breakpoint position)
s2_df <- performance_results %>%
  filter(has_change) %>%
  group_by(scenario, effect, phi, n_years) %>%
  summarise(useful_detection = mean(within_2, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    scenario_label = recode(scenario,
      abrupt_moderate = "Abrupt moderate",
      abrupt_strong = "Abrupt strong",
      expansion_peak_collapse = "Expansion-peak-collapse",
      gradual_expansion = "Gradual expansion"
    )
  ) %>%
  filter(!is.na(useful_detection), scenario != "no_change")

p_s2 <- ggplot(s2_df, aes(x = factor(n_years), y = factor(scenario_label,
                                                         levels = c("Gradual expansion", "Expansion-peak-collapse", "Abrupt strong", "Abrupt moderate")),
                          fill = useful_detection)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = round(100 * useful_detection)), size = 4) +
  facet_grid(effect ~ phi, labeller = label_both, switch = "y") +
  scale_fill_gradient(low = "#f7f7f7", high = "#333333", limits = c(0, 1), labels = function(x) paste0(100 * x, "%"), name = "Useful\ndetection") +
  labs(
    title = "Useful detection across simulation conditions",
    subtitle = "Cell values show percentage of detections within +/-2 years",
    x = "Time-series length (years)",
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "grey90"),
        legend.position = "right")

ggsave(file.path(FIG_DIR, "Supplementary_Figure_S2.png"), p_s2, width = 14, height = 10, dpi = 300)

# S3. Detailed temporal error across conditions
s3_df <- performance_results %>%
  filter(tp, scenario != "no_change") %>%
  mutate(
    scenario_label = recode(scenario,
      abrupt_moderate = "Abrupt moderate",
      abrupt_strong = "Abrupt strong",
      expansion_peak_collapse = "Expansion-peak-collapse",
      gradual_expansion = "Gradual expansion"
    )
  )

p_s3 <- ggplot(s3_df, aes(x = factor(n_years), y = error_years)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(outlier.alpha = 0.15, width = 0.6) +
  facet_grid(scenario_label ~ phi, labeller = label_both, switch = "y") +
  labs(
    title = "Detailed temporal error across simulation conditions",
    subtitle = "Estimated transition year minus true simulated transition boundary",
    x = "Time-series length (years)",
    y = "Temporal error (years)"
  ) +
  theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "grey90"))

ggsave(file.path(FIG_DIR, "Supplementary_Figure_S3.png"), p_s3, width = 14, height = 14, dpi = 300)

# S4. Weight sensitivity of utility
p_s4 <- summary_by_threshold %>%
  mutate(weight_label = factor(weight_label,
                               levels = c("Equal", "Effect-priority", "Increase-priority", "No-p-value"))) %>%
  ggplot(aes(x = weight_label, y = utility)) +
  geom_boxplot(fill = "grey80", width = 0.55) +
  stat_summary(fun = mean, geom = "point", size = 3, colour = "black") +
  labs(
    title = "Weight-sensitivity analysis",
    subtitle = "Utility = useful detection within +/-2 years - 2 x false-positive rate",
    x = NULL,
    y = "Utility"
  ) +
  theme_bw(base_size = 13)

ggsave(file.path(FIG_DIR, "Supplementary_Figure_S4.png"), p_s4, width = 11, height = 7, dpi = 300)

# S5. Trade-off between useful detection and false positives
tradeoff_df <- summary_by_threshold %>%
  mutate(weight_label = factor(weight_label,
                               levels = c("Equal", "Effect-priority", "Increase-priority", "No-p-value")))

p_s5 <- ggplot(tradeoff_df, aes(x = false_positive_rate, y = within_2_rate, colour = weight_label)) +
  geom_point(alpha = 0.65, size = 2.8) +
  geom_vline(xintercept = c(0.05, 0.10, 0.15), linetype = "dashed", alpha = 0.7) +
  scale_x_continuous(labels = function(x) paste0(round(100 * x), "%")) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
  labs(
    title = "Trade-off between detection and false positives",
    subtitle = "Dashed lines indicate false-positive constraints of 5%, 10% and 15%",
    x = "False-positive rate",
    y = "Useful detection rate within +/-2 years",
    colour = "Weighting scheme"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

ggsave(file.path(FIG_DIR, "Supplementary_Figure_S5.png"), p_s5, width = 12, height = 8, dpi = 300)

# -------------------------------- MAIN FIGURE 2 --------------------------------------

# Figure 2A. Useful detection by scenario
fig2a_df <- performance_results %>%
  filter(has_change, scenario != "no_change") %>%
  group_by(scenario) %>%
  summarise(useful_detection = mean(within_2, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    scenario_label = factor(recode(scenario,
      abrupt_moderate = "Abrupt moderate",
      abrupt_strong = "Abrupt strong",
      expansion_peak_collapse = "Expansion-peak-collapse",
      gradual_expansion = "Gradual expansion"
    ), levels = c("Abrupt moderate", "Abrupt strong", "Expansion-peak-collapse", "Gradual expansion"))
  )

p2a <- ggplot(fig2a_df, aes(x = scenario_label, y = useful_detection)) +
  geom_col(fill = "grey80", colour = "grey30", width = 0.62) +
  geom_text(aes(label = paste0(round(100 * useful_detection), "%")), vjust = -0.4, size = 4.5) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, 1.08)) +
  labs(
    title = "A. Useful detection",
    subtitle = "Detections within +/-2 years of the true simulated transition boundary",
    x = NULL,
    y = "Useful detection rate"
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

# Figure 2B. False-positive rate by noise, autocorrelation and time-series length
fig2b_df <- performance_results %>%
  filter(!has_change) %>%
  group_by(cv, n_years, phi) %>%
  summarise(
    fp_n = sum(fp),
    n = n(),
    false_positive_rate = mean(detected, na.rm = TRUE),
    .groups = "drop"
  )

fig2b_ci <- t(mapply(function(x, n) wilson_ci(x, n), fig2b_df$fp_n, fig2b_df$n))
fig2b_df$ci_low  <- fig2b_ci[, 1]
fig2b_df$ci_high <- fig2b_ci[, 2]
fig2b_df <- fig2b_df %>%
  mutate(cv_label = paste0("CV = ", cv))

p2b <- ggplot(fig2b_df, aes(x = n_years, y = false_positive_rate, colour = factor(phi), group = factor(phi))) +
  geom_hline(yintercept = 0.10, linetype = "dashed") +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.35, alpha = 0.7) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.8) +
  facet_wrap(~cv_label, nrow = 1) +
  scale_colour_manual(values = c("0" = "#1F77B4", "0.3" = "#E69F00", "0.6" = "#009E73"), name = "AR(1) autocorrelation, phi") +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
  labs(
    title = "B. False-positive rate",
    subtitle = "Points show estimates; vertical bars show 95% Wilson confidence intervals",
    x = "Time-series length (years)",
    y = "False-positive rate"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# Figure 2C. Temporal accuracy by scenario
fig2c_df <- performance_results %>%
  filter(tp, scenario != "no_change") %>%
  group_by(scenario) %>%
  summarise(
    med = median(error_years, na.rm = TRUE),
    q1 = quantile(error_years, 0.25, na.rm = TRUE),
    q3 = quantile(error_years, 0.75, na.rm = TRUE),
    mae = mean(abs_error, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    scenario_label = factor(recode(scenario,
      abrupt_moderate = "Abrupt moderate",
      abrupt_strong = "Abrupt strong",
      expansion_peak_collapse = "Expansion-peak-collapse",
      gradual_expansion = "Gradual expansion"
    ), levels = c("Abrupt moderate", "Abrupt strong", "Expansion-peak-collapse", "Gradual expansion"))
  )

p2c <- ggplot(fig2c_df, aes(x = scenario_label, y = med)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = q1, ymax = q3), width = 0.15, linewidth = 1) +
  geom_point(size = 3.2) +
  geom_text(aes(label = paste0("MAE = ", round(mae, 1), " y")), nudge_y = 0.25, size = 4.2) +
  labs(
    title = "C. Temporal accuracy",
    subtitle = "Points show median temporal error; bars show interquartile range",
    x = NULL,
    y = "Temporal error (years)"
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

# Save panels separately
fig2_dir <- file.path(FIG_DIR, "Figure_2_panels")
dir.create(fig2_dir, showWarnings = FALSE)
ggsave(file.path(fig2_dir, "Figure_2A_useful_detection.png"), p2a, width = 8, height = 5, dpi = 300)
ggsave(file.path(fig2_dir, "Figure_2B_false_positive_rate.png"), p2b, width = 10, height = 6, dpi = 300)
ggsave(file.path(fig2_dir, "Figure_2C_temporal_accuracy.png"), p2c, width = 8, height = 5, dpi = 300)

# Try to save a combined Figure 2 if patchwork is available
if (requireNamespace("patchwork", quietly = TRUE)) {
  combined_fig2 <- p2a / p2b / p2c + patchwork::plot_layout(heights = c(1, 1.25, 1))
  ggsave(file.path(FIG_DIR, "Figure_2_combined.png"), combined_fig2, width = 10, height = 16, dpi = 300)
}

# -------------------------------- FINAL EXPORTS --------------------------------------

write_csv(summary_by_condition, file.path(TAB_DIR, "summary_by_condition_selected_setting.csv"))
write_csv(fig2a_df, file.path(TAB_DIR, "figure2A_values.csv"))
write_csv(fig2b_df, file.path(TAB_DIR, "figure2B_values.csv"))
write_csv(fig2c_df, file.path(TAB_DIR, "figure2C_values.csv"))

message("Simulation workflow complete.")
message("Outputs written to:")
message("  Figures: ", FIG_DIR)
message("  Tables : ", TAB_DIR)
message("Selected threshold/weight combination used in phase 2:")
print(selected_row)
