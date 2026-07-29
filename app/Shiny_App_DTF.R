# =======================================================================================
# DTF – New Fisheries Decision Tool with Automated and Expert-Supported Pathways (ENGLISH)
# Rewritten for clarity and usability for external users (reviewer-friendly)
#
# Key usability improvements:
# - "How to use" is the FIRST tab and acts as the entry point
# - Clear, step-by-step guidance + troubleshooting inside the app
# - A "Demo: completed example run" that matches the provided example pathway
# - A "Load example data" + "Download example CSV" to practice data format
# - Decision Tree defaults to "Path only" view for readability
# - Path is shown as: (1) highlighted graph, (2) readable step list, (3) bullet list
# - Outputs use plain language, not internal IDs
# - Two initial pathways: automated landing screening for eligible series and an expert-supported pathway for shorter series
# =======================================================================================

library(shiny)
library(bslib)
library(visNetwork)
library(DT)
library(ggplot2)
library(gridExtra)
library(grid)
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

# ------------------------------ ASSUMPTION-AWARE BREAKPOINT ANALYSIS -------------------------
#
# The implementation follows the analytical procedure described in the manuscript:
# 1) candidate splits require at least three annual observations before and after;
# 2) normality is assessed separately in the pre- and post-breakpoint groups using Shapiro-Wilk;
# 3) homogeneity of variance is assessed using Fligner-Killeen;
# 4) a one-sided Student t-test is used when both assumptions are supported;
# 5) otherwise, a one-sided Wilcoxon rank-sum / Mann-Whitney U test is used;
# 6) the selected p-value, Cohen's d and percentage increase are combined in the DTF score;
# 7) Ljung-Box and lag-1 autocorrelation are reported as time-series diagnostics.
#
# These outputs are diagnostic early-warning evidence and not confirmatory stock-assessment results.

safe_p_value <- function(expr) {
  tryCatch({
    out <- suppressWarnings(expr)
    if (length(out) == 0 || !is.finite(out)) NA_real_ else as.numeric(out)
  }, error = function(e) NA_real_)
}

format_p_value <- function(x, digits = 3) {
  if (length(x) == 0 || is.na(x) || !is.finite(x)) return("not evaluated")
  format.pval(x, digits = digits, eps = 10^-digits)
}

# Check whether the uploaded series is eligible for automated screening.
# The quantitative pathway requires at least six unique, consecutive annual observations.
landing_series_status <- function(df, min_years = 6) {
  if (is.null(df) || !all(c("year", "landing") %in% names(df))) {
    return(list(
      eligible = FALSE,
      n_years = 0L,
      consecutive = FALSE,
      duplicated_years = FALSE,
      non_positive_landings = FALSE,
      message = "No valid landing series is currently available."
    ))
  }

  years <- suppressWarnings(as.numeric(df$year))
  landings <- suppressWarnings(as.numeric(df$landing))
  valid <- is.finite(years) & is.finite(landings)
  years <- years[valid]
  landings <- landings[valid]

  n_years <- length(years)
  duplicated_years <- anyDuplicated(years) > 0
  consecutive <- n_years >= 2 && !duplicated_years && all(diff(sort(years)) == 1)
  non_positive_landings <- any(landings <= 0, na.rm = TRUE)
  eligible <- n_years >= min_years && consecutive

  if (duplicated_years) {
    message <- "The uploaded series contains duplicated years and is not eligible for automated screening."
  } else if (n_years < min_years) {
    message <- paste0(
      "The uploaded series contains ", n_years, " annual observations; at least ",
      min_years, " consecutive observations are required for automated screening."
    )
  } else if (!consecutive) {
    message <- paste0(
      "The uploaded series contains gaps between years. At least ", min_years,
      " consecutive annual observations are required for automated screening."
    )
  } else {
    message <- paste0(
      "The uploaded series is eligible for automated screening (", n_years,
      " consecutive annual observations)."
    )
  }

  if (eligible && non_positive_landings) {
    message <- paste(
      message,
      "Warning: zero or negative landing values are present; percentage increases may be unstable."
    )
  }

  list(
    eligible = eligible,
    n_years = n_years,
    consecutive = consecutive,
    duplicated_years = duplicated_years,
    non_positive_landings = non_positive_landings,
    message = message
  )
}

series_temporal_diagnostics <- function(df) {
  x <- as.numeric(df$landing)
  n <- length(x)

  lb_lag <- if (n >= 4) min(3L, max(1L, floor(n / 4)), n - 1L) else NA_integer_
  ljung_box_p <- if (!is.na(lb_lag)) {
    safe_p_value(stats::Box.test(x, lag = lb_lag, type = "Ljung-Box")$p.value)
  } else {
    NA_real_
  }

  acf_lag1 <- if (n >= 3) {
    acf_obj <- tryCatch(
      stats::acf(x, plot = FALSE, na.action = na.pass),
      error = function(e) NULL
    )
    if (is.null(acf_obj) || length(acf_obj$acf) < 2) NA_real_ else as.numeric(acf_obj$acf[2])
  } else {
    NA_real_
  }

  list(
    n_years = n,
    ljung_box_lag = lb_lag,
    ljung_box_p = ljung_box_p,
    acf_lag1 = acf_lag1,
    autocorrelation_warning =
      (!is.na(ljung_box_p) && ljung_box_p < 0.05) ||
      (!is.na(acf_lag1) && abs(acf_lag1) >= 0.40)
  )
}

candidate_breakpoint_result <- function(df, idx, alpha = 0.05) {
  g_before <- as.numeric(df$landing[1:idx])
  g_after <- as.numeric(df$landing[(idx + 1):nrow(df)])

  shapiro_before_p <- if (
    length(g_before) >= 3 && length(g_before) <= 5000 &&
    length(unique(g_before)) >= 3
  ) {
    safe_p_value(stats::shapiro.test(g_before)$p.value)
  } else {
    NA_real_
  }

  shapiro_after_p <- if (
    length(g_after) >= 3 && length(g_after) <= 5000 &&
    length(unique(g_after)) >= 3
  ) {
    safe_p_value(stats::shapiro.test(g_after)$p.value)
  } else {
    NA_real_
  }

  group <- factor(c(rep("Before", length(g_before)), rep("After", length(g_after))))
  values <- c(g_before, g_after)
  fligner_p <- safe_p_value(stats::fligner.test(values ~ group)$p.value)

  # "Reasonably supported" is interpreted conservatively:
  # both group-normality tests and the homogeneity test must be evaluable and non-significant.
  normality_supported <-
    is.finite(shapiro_before_p) && shapiro_before_p >= alpha &&
    is.finite(shapiro_after_p) && shapiro_after_p >= alpha
  homoscedasticity_supported <- is.finite(fligner_p) && fligner_p >= alpha

  student_t_p <- safe_p_value(
    stats::t.test(g_after, g_before, alternative = "greater", var.equal = TRUE)$p.value
  )
  wilcoxon_p <- safe_p_value(
    stats::wilcox.test(g_after, g_before, alternative = "greater", exact = FALSE)$p.value
  )

  if (normality_supported && homoscedasticity_supported && is.finite(student_t_p)) {
    selected_test <- "One-sided Student t-test"
    selected_p <- student_t_p
  } else {
    selected_test <- "One-sided Wilcoxon rank-sum test"
    selected_p <- wilcoxon_p
  }

  # Defensive fallback in the rare case where the selected test cannot be evaluated.
  if (!is.finite(selected_p) && is.finite(student_t_p)) {
    selected_test <- "One-sided Student t-test (fallback)"
    selected_p <- student_t_p
  }

  mean_before <- mean(g_before, na.rm = TRUE)
  mean_after <- mean(g_after, na.rm = TRUE)
  var_before <- stats::var(g_before, na.rm = TRUE)
  var_after <- stats::var(g_after, na.rm = TRUE)

  pooled_sd <- sqrt(
    ((length(g_before) - 1) * var_before + (length(g_after) - 1) * var_after) /
      (length(g_before) + length(g_after) - 2)
  )

  cohens_d <- if (is.finite(pooled_sd) && pooled_sd > 0) {
    (mean_after - mean_before) / pooled_sd
  } else {
    NA_real_
  }

  percent_increase <- if (is.finite(mean_before) && mean_before != 0) {
    ((mean_after - mean_before) / mean_before) * 100
  } else {
    NA_real_
  }

  score <- if (
    is.finite(selected_p) && is.finite(cohens_d) &&
    is.finite(percent_increase) && percent_increase > 0
  ) {
    (1 - selected_p) * abs(cohens_d) * (percent_increase / 100)
  } else {
    NA_real_
  }

  supported <- isTRUE(
    is.finite(selected_p) && selected_p < alpha &&
    is.finite(percent_increase) && percent_increase > 0 &&
    is.finite(score)
  )

  data.frame(
    breakpoint_year = df$year[idx],
    idx = idx,
    n_before = length(g_before),
    n_after = length(g_after),
    mean_before = mean_before,
    mean_after = mean_after,
    percent_increase = percent_increase,
    cohens_d = cohens_d,
    shapiro_before_p = shapiro_before_p,
    shapiro_after_p = shapiro_after_p,
    fligner_p = fligner_p,
    normality_supported = normality_supported,
    homoscedasticity_supported = homoscedasticity_supported,
    student_t_p = student_t_p,
    wilcoxon_p = wilcoxon_p,
    selected_test = selected_test,
    selected_p = selected_p,
    score = score,
    supported = supported,
    stringsAsFactors = FALSE
  )
}

run_breakpoint_pipeline <- function(df, alpha = 0.05, min_points = 3) {
  df <- df[order(df$year), , drop = FALSE]
  status <- landing_series_status(df, min_years = min_points * 2)

  if (!isTRUE(status$eligible)) {
    return(list(
      status = status,
      temporal_diagnostics = NULL,
      candidates = data.frame(),
      selected_breakpoint = NULL,
      alpha = alpha,
      min_points = min_points
    ))
  }

  candidate_rows <- lapply(
    min_points:(nrow(df) - min_points),
    function(idx) candidate_breakpoint_result(df, idx, alpha = alpha)
  )
  candidates <- do.call(rbind, candidate_rows)
  rownames(candidates) <- NULL

  supported_candidates <- candidates[
    candidates$supported & is.finite(candidates$score),
    ,
    drop = FALSE
  ]

  # Retain the highest-ranked candidate for transparent diagnostic display even when
  # it does not meet the one-sided statistical-support criterion used to route the DTF.
  ranked_candidates <- candidates[
    is.finite(candidates$score) & is.finite(candidates$percent_increase) &
      candidates$percent_increase > 0,
    ,
    drop = FALSE
  ]

  highest_ranked <- NULL
  if (nrow(ranked_candidates) > 0) {
    highest_row <- ranked_candidates[which.max(ranked_candidates$score), , drop = FALSE]
    highest_ranked <- as.list(highest_row[1, , drop = FALSE])
    highest_ranked$year <- highest_ranked$breakpoint_year
  }

  selected <- NULL
  if (nrow(supported_candidates) > 0) {
    selected_row <- supported_candidates[which.max(supported_candidates$score), , drop = FALSE]
    selected <- as.list(selected_row[1, , drop = FALSE])
    selected$year <- selected$breakpoint_year
  }

  list(
    status = status,
    temporal_diagnostics = series_temporal_diagnostics(df),
    candidates = candidates,
    highest_ranked_candidate = highest_ranked,
    selected_breakpoint = selected,
    alpha = alpha,
    min_points = min_points
  )
}

# Backward-compatible function names retained for older code and saved app sessions.
statistical_breakpoint_analysis <- function(df, alpha = 0.05, min_points = 3, test_type = "both") {
  run_breakpoint_pipeline(df, alpha = alpha, min_points = min_points)
}

select_best_breakpoint <- function(res) {
  if (is.null(res)) return(NULL)
  if (!is.null(res$selected_breakpoint)) return(res$selected_breakpoint)
  NULL
}

plot_breakpoint <- function(df, bp) {
  if (is.null(bp)) {
    return(
      ggplot(df, aes(year, landing)) +
        geom_line(linewidth = 0.9, alpha = 0.85) +
        geom_point(size = 2.4) +
        labs(
          title = "No supported breakpoint detected",
          subtitle = "No candidate split met the selected one-sided statistical-support criterion.",
          x = "Year", y = "Landing"
        ) +
        theme_minimal(base_size = 14)
    )
  }

  idx <- as.integer(bp$idx)
  df_means <- data.frame(
    year = c(df$year[1], df$year[idx], df$year[idx + 1], df$year[nrow(df)]),
    mean = c(bp$mean_before, bp$mean_before, bp$mean_after, bp$mean_after),
    group = c("Before", "Before", "After", "After")
  )

  ggplot(df, aes(year, landing)) +
    geom_line(linewidth = 0.9, alpha = 0.80) +
    geom_point(size = 2.4) +
    geom_vline(xintercept = bp$year, color = "#E74C3C", linewidth = 1.2) +
    geom_line(data = df_means, aes(year, mean, linetype = group), linewidth = 1) +
    scale_linetype_manual(values = c("Before" = "dashed", "After" = "dashed")) +
    labs(
      title = paste("Selected diagnostic breakpoint – year:", bp$year),
      subtitle = paste0(
        bp$selected_test, "; p = ", format_p_value(bp$selected_p),
        "; mean increase = ", round(bp$percent_increase, 1), "%"
      ),
      x = "Year", y = "Landing"
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none")
}

build_analysis_grob <- function(df, bp, pipeline = NULL) {
  p <- plot_breakpoint(df, bp)

  if (is.null(bp)) {
    stats_df <- data.frame(
      Item = c("Result", "Meaning"),
      Value = c(
        "No supported breakpoint",
        "No candidate split met the selected one-sided statistical-support criterion."
      )
    )
  } else {
    temporal <- if (!is.null(pipeline)) pipeline$temporal_diagnostics else NULL
    stats_df <- data.frame(
      Item = c(
        "Breakpoint year", "Selected test", "Selected p-value",
        "Increase (%)", "Cohen's d", "Composite score",
        "Shapiro p (before)", "Shapiro p (after)", "Fligner-Killeen p",
        "Ljung-Box p", "Lag-1 autocorrelation"
      ),
      Value = c(
        bp$year,
        bp$selected_test,
        format_p_value(bp$selected_p),
        round(bp$percent_increase, 1),
        round(bp$cohens_d, 3),
        round(bp$score, 3),
        format_p_value(bp$shapiro_before_p),
        format_p_value(bp$shapiro_after_p),
        format_p_value(bp$fligner_p),
        if (is.null(temporal)) "not evaluated" else format_p_value(temporal$ljung_box_p),
        if (is.null(temporal) || is.na(temporal$acf_lag1)) "not evaluated" else round(temporal$acf_lag1, 3)
      ),
      stringsAsFactors = FALSE
    )
  }

  tab <- tableGrob(
    stats_df,
    rows = NULL,
    theme = ttheme_minimal(base_size = 11, padding = unit(c(4, 4), "mm"))
  )
  gridExtra::arrangeGrob(p, tab, ncol = 2, widths = c(2.2, 1.15))
}

candidate_table_for_display <- function(pipeline) {
  if (is.null(pipeline) || is.null(pipeline$candidates) || nrow(pipeline$candidates) == 0) {
    return(data.frame(Message = "No candidate breakpoint results are available."))
  }

  x <- pipeline$candidates
  score_rank <- rep(NA_integer_, nrow(x))
  finite_score <- is.finite(x$score)
  if (any(finite_score)) {
    score_rank[finite_score] <- rank(-x$score[finite_score], ties.method = "min")
  }

  selected_year <- if (is.null(pipeline$selected_breakpoint)) {
    NA_real_
  } else {
    as.numeric(pipeline$selected_breakpoint$year)
  }

  data.frame(
    `Score rank` = score_rank,
    `Breakpoint year` = x$breakpoint_year,
    `n before` = x$n_before,
    `n after` = x$n_after,
    `Mean before` = round(x$mean_before, 3),
    `Mean after` = round(x$mean_after, 3),
    `Increase (%)` = round(x$percent_increase, 2),
    `Cohen's d` = round(x$cohens_d, 3),
    `Shapiro p (before)` = signif(x$shapiro_before_p, 3),
    `Shapiro p (after)` = signif(x$shapiro_after_p, 3),
    `Normality supported` = ifelse(x$normality_supported, "Yes", "No"),
    `Fligner-Killeen p` = signif(x$fligner_p, 3),
    `Homoscedasticity supported` = ifelse(x$homoscedasticity_supported, "Yes", "No"),
    `Selected test` = x$selected_test,
    `Student t p` = signif(x$student_t_p, 3),
    `Wilcoxon p` = signif(x$wilcoxon_p, 3),
    `Selected p` = signif(x$selected_p, 3),
    `Composite score` = round(x$score, 4),
    `Statistical support` = ifelse(x$supported, "Yes", "No"),
    `DTF-selected signal` = ifelse(
      is.finite(selected_year) & x$breakpoint_year == selected_year,
      "Yes",
      "No"
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

selected_candidate_table <- function(pipeline) {
  if (is.null(pipeline)) {
    return(data.frame(
      Diagnostic = "Analysis status",
      Result = "No analysis has been run.",
      Interpretation = "Upload an eligible landing series or use the Run / refresh analysis button.",
      stringsAsFactors = FALSE
    ))
  }

  candidate <- pipeline$selected_breakpoint
  candidate_status <- "DTF-selected diagnostic signal"
  if (is.null(candidate)) {
    candidate <- pipeline$highest_ranked_candidate
    candidate_status <- "Highest-ranked candidate; not selected as a DTF signal"
  }

  if (is.null(candidate)) {
    return(data.frame(
      Diagnostic = "Candidate result",
      Result = "No positive candidate split was available.",
      Interpretation = "The series did not provide a positive finite breakpoint score.",
      stringsAsFactors = FALSE
    ))
  }

  normality_text <- if (isTRUE(candidate$normality_supported)) "Supported" else "Not supported"
  variance_text <- if (isTRUE(candidate$homoscedasticity_supported)) "Supported" else "Not supported"
  support_text <- if (isTRUE(candidate$supported)) "Yes" else "No"

  data.frame(
    Diagnostic = c(
      "Candidate status",
      "Breakpoint year",
      "Observations before / after",
      "Shapiro-Wilk p-value: before",
      "Shapiro-Wilk p-value: after",
      "Normality assumption",
      "Fligner-Killeen p-value",
      "Homogeneity-of-variance assumption",
      "Selected one-sided comparison",
      "Student t-test p-value",
      "Wilcoxon rank-sum p-value",
      "Selected p-value",
      "Mean before / after",
      "Increase in mean landings",
      "Cohen's d",
      "Composite score",
      "Statistical support at alpha = 0.05"
    ),
    Result = c(
      candidate_status,
      as.character(candidate$year),
      paste0(candidate$n_before, " / ", candidate$n_after),
      format_p_value(candidate$shapiro_before_p),
      format_p_value(candidate$shapiro_after_p),
      normality_text,
      format_p_value(candidate$fligner_p),
      variance_text,
      candidate$selected_test,
      format_p_value(candidate$student_t_p),
      format_p_value(candidate$wilcoxon_p),
      format_p_value(candidate$selected_p),
      paste0(round(candidate$mean_before, 3), " / ", round(candidate$mean_after, 3)),
      paste0(round(candidate$percent_increase, 2), "%"),
      if (is.finite(candidate$cohens_d)) round(candidate$cohens_d, 3) else "not evaluated",
      if (is.finite(candidate$score)) round(candidate$score, 4) else "not evaluated",
      support_text
    ),
    Interpretation = c(
      if (isTRUE(candidate$supported)) {
        "This candidate is used to route the automated DTF pathway."
      } else {
        "Shown for transparency; it does not route the DTF to the positive-signal branch."
      },
      "Last year in the pre-breakpoint segment.",
      "Candidate splits require at least three observations in each segment.",
      if (is.finite(candidate$shapiro_before_p) && candidate$shapiro_before_p >= 0.05) "No evidence against normality in the pre-breakpoint segment." else "Normality is not supported or could not be evaluated.",
      if (is.finite(candidate$shapiro_after_p) && candidate$shapiro_after_p >= 0.05) "No evidence against normality in the post-breakpoint segment." else "Normality is not supported or could not be evaluated.",
      "Both segment-specific Shapiro-Wilk tests must be evaluable and p >= 0.05.",
      if (is.finite(candidate$fligner_p) && candidate$fligner_p >= 0.05) "No evidence of unequal variances." else "Homogeneity of variance is not supported or could not be evaluated.",
      "Fligner-Killeen p >= 0.05 is treated as support for homoscedasticity.",
      "Student's t-test is selected only when normality and homoscedasticity are supported; otherwise Wilcoxon is used.",
      "Reported for transparency even when it is not selected.",
      "Reported for transparency even when it is not selected.",
      "This p-value enters the composite score and the statistical-support criterion.",
      "Segment means used to quantify the change.",
      "Relative difference between post- and pre-breakpoint means.",
      "Standardised magnitude of the difference.",
      "Operational ranking metric: (1 - p) × |d| × proportional increase.",
      "A positive increase and selected one-sided p < 0.05 are required."
    ),
    stringsAsFactors = FALSE
  )
}

# ------------------------------- DTF STRUCTURE --------------------------------------------------

dtf_original <- list(
  "start" = list(
    type = "question",
    text = "Are you aware that fishers are using new fishing gear?",
    choices = list("Yes" = "q6", "No" = "q_years")
  ),
  
  "q_years" = list(
    type = "question",
    text = "Are at least six consecutive annual observations of total landings available for the target species?",
    choices = list("Yes" = "q2", "No" = "q4a")
  ),
  
  "q2" = list(
    type = "question",
    text = "Does the automated landing-based screening identify a sudden and persistent increase in total landings across all fishing gears targeting the species?",
    requires_data = TRUE,
    choices = list("Yes" = "q3", "No" = "a_continue_monitoring")
  ),
  
  "q3" = list(type = "question", text = "Did this peak occur in the past?", choices = list("Yes" = "q11", "No" = "q4a")),
  
  "q4a" = list(
    type = "question",
    text = "Are fishing-effort data or reliable effort proxies available?",
    choices = list(
      "Yes" = "q4",
      "No" = "a_collect_effort_info"
    )
  ),
  
  "a_collect_effort_info" = list(
    type = "question",
    text = paste(
      "Collect any available information on fishing effort before proceeding.",
      "Effort proxies can include fishing capacity (number of vessels) and fishing activity (number of fishing days or trips).",
      "How would you like to proceed?"
    ),
    choices = list(
      "Proceed based on expert knowledge that the fishing pattern may have changed" = "q5",
      "Stop here and collect all needed information" = "a_stop_rerun_effort"
    )
  ),
  "a_stop_rerun_effort" = list(
    type = "action",
    text = paste(
      "Stop here and collect fishing-effort information or reliable effort proxies before continuing.",
      "Relevant proxies can include fishing capacity (number of vessels) and fishing activity (number of fishing days or trips).",
      "Once the information is available, re-run the tool."
    )
  ),
  
  "q4" = list(
    type = "question",
    text = "Do the available fishing-effort data, together with expert knowledge, suggest that the fishing pattern has changed?",
    choices = list(
      "Yes" = "q5",
      "No" = "q17"
    )
  ),
  
  "q5" = list(
    type = "question",
    text = "What are the reasons for the change in the fishing pattern?",
    choices = list(
      "New fishing gear" = "q6",
      "New fishing grounds" = "q17",
      "The reasons for the change in fishing pattern are unknown" = "a_collect_info"
    )
  ),
  
  "a_collect_info" = list(
    type = "action",
    text = paste(
      "Conduct a frame survey to collect the information needed to identify the reason behind the change in fishing pattern. Re-run the tool when the information becomes available."
    )
  ),
  
  "q6" = list(
    type = "question",
    text = "Do you know the features of the new fishing gear?",
    choices = list("Yes" = "q7", "No, Collect any information on the new gear" = "q7_bis")
  ),
  
  "q7" = list(
    type = "question",
    text = "Could you provide an assessment of the target stock and evaluate the impact of the new fishing gear?",
    choices = list("Yes" = "q8", "No, Collect more data and consider the stock as overexploited" = "q10")
  ),
  
  "q7_bis" = list(
    type = "question",
    text = "Could you provide an assessment of the target stock and evaluate the impact of the new fishing gear?",
    choices = list("Yes" = "q8", "No, Collect more data and consider the stock as overexploited" = "q10")
  ),
  
  "q8" = list(type = "question", text = "What is the stock status?", choices = list("Sustainable" = "q9", "Overexploited" = "q10")),
  
  "q9" = list(type = "question", text = "Are any management measures in place?",
              choices = list("Yes" = "a_continue_monitoring_assess", "No" = "a_establish_management_measures")),
  
  "q10" = list(type = "question", text = "Are any management measures in place?",
               choices = list("Yes" = "a_modify_management_measures", "No" = "a_define_management_measures")),
  
  "q11" = list(type = "question", text = "Do you have biomass or catch trends?", choices = list("Yes" = "q12", "No" = "q4a")),
  
  "q12" = list(
    type = "question",
    text = "Does the biomass or catch trend suggest the sudden increase is due to natural stock fluctuations?",
    choices = list("Yes" = "q13", "No" = "q4a")
  ),
  
  "q13" = list(
    type = "question",
    text = "Could you provide an assessment of the stock status?",
    choices = list("Yes" = "q14", "No, Improve the available data and provide an assessment of the stock status for the target species" = "q14")
  ),
  
  "q14" = list(type = "question", text = "What is the stock status?", choices = list("Sustainable" = "q15", "Overexploited" = "q16")),
  
  "q15" = list(type = "question", text = "Are any management measures in place?",
               choices = list("Yes" = "a_continue_monitoring", "No" = "a_define_appropriate_management")),
  
  "q16" = list(type = "question", text = "Are any management measures in place?",
               choices = list("Yes" = "a_assess_modify_measures", "No" = "a_assess_define_measures")),
  
  "q17" = list(type = "question", text = "Do you have a trend of the biomass index or catch rates for the target species?",
               choices = list("Yes" = "q18", "No" = "q13")),
  
  "q18" = list(
    type = "question",
    text = "Does the biomass index or catch rates increase during the same period when landings rise?",
    choices = list("Yes" = "q13", "No, There has probably been an increase in fishing effort on the target species" = "q13")
  ),
  
  "q19" = list(type = "question", text = "What are the reasons for the change in the fishing pattern?",
               choices = list("New fishing gear" = "q6", "New fishing grounds" = "q17")),
  
  "a_continue_monitoring" = list(type = "action", text = "Continue monitoring."),
  
  "a_continue_monitoring_assess" = list(type = "action", text = "Continue monitoring and assess the species annually."),
  
  "a_establish_management_measures" = list(
    type = "action",
    text = "Establish management measures to limit the catch from the new fishing gear: establish licensing, effort controls, and quota regulations and technical measures before the fishery expands uncontrollably. Continue monitoring and assess the species annually."
  ),
  
  "a_modify_management_measures" = list(
    type = "action",
    text = "Modify the current management measures and implement a precautionary LIMIT (i.e. reduce fishing days) or STOP the fishery using the new gear. Continue monitoring and assess the species annually."
  ),
  
  "a_define_management_measures" = list(
    type = "action",
    text = "Define appropriate management measures and implement a precautionary LIMIT (i.e. reduce fishing days) or STOP the fishery using the new gear. Continue monitoring and assess the species annually."
  ),
  
  "a_define_appropriate_management" = list(type = "action", text = "Define appropriate management measure and continue monitoring."),
  
  "a_assess_modify_measures" = list(
    type = "action",
    text = "Assess the stock status annually and modify the management measures, such as licensing, quotas, effort limits, and closures, as needed."
  ),
  
  "a_assess_define_measures" = list(
    type = "action",
    text = "Assess the stock status annually and define the management measures, such as licensing, quotas, effort limits, and closures, as needed."
  )
)

# ------------------------------- READABLE LABEL HELPERS -----------------------------------------

create_wrapped_label <- function(text, max_chars_per_line = 50) {
  words <- strsplit(text, "\\s+")[[1]]
  if (nchar(text) <= max_chars_per_line) return(text)
  lines <- character(); current <- ""
  for (w in words) {
    test <- if (current == "") w else paste(current, w)
    if (nchar(test) <= max_chars_per_line) current <- test
    else { if (nzchar(current)) lines <- c(lines, current); current <- w }
  }
  if (nzchar(current)) lines <- c(lines, current)
  paste(lines, collapse = "\n")
}

gv_label <- function(x, max_chars_per_line = 42) {
  txt <- create_wrapped_label(x, max_chars_per_line)
  txt <- gsub("\\\\", "\\\\\\\\", txt)
  txt <- gsub("\"", "\\\\\"", txt)
  gsub("\n", "\\n", txt, fixed = TRUE)
}

filter_dtf_by_path <- function(original_dtf, path_taken) {
  filtered <- list()
  for (id in path_taken) if (id %in% names(original_dtf)) filtered[[id]] <- original_dtf[[id]]
  filtered
}

build_graphviz_dot <- function(dtf_data, path_taken) {
  if (!length(dtf_data)) return("digraph G { rankdir=TB; }")
  nodes_ids <- names(dtf_data)
  
  dot <- 'digraph G {
    graph [rankdir=TB, nodesep=0.6, ranksep=1.0, bgcolor="white"];
    node  [shape=box, style=filled, color="#123A63", fillcolor="#E8F2FF",
           fontname="Helvetica", fontsize=26, margin="0.18,0.12"];
    edge  [fontname="Helvetica", fontsize=22, color="#5B6B7A", arrowsize=0.8];
  '
  
  for (id in nodes_ids) {
    nd <- dtf_data[[id]]
    is_action <- nd$type == "action"
    fill <- if (is_action) "#E6FFF2" else "#E8F2FF"
    border <- if (is_action) "#1F7A4D" else "#123A63"
    lbl <- if (is_action) paste0("RECOMMENDATION: ", nd$text) else nd$text
    lbl <- gv_label(lbl, if (is_action) 46 else 42)
    
    if (id %in% path_taken) {
      fill <- if (is_action) "#1F7A4D" else "#123A63"
      border <- fill
      dot <- paste0(dot, sprintf('"%s" [label="%s", fillcolor="%s", color="%s", fontcolor="white"];\n',
                                 id, lbl, fill, border))
    } else {
      dot <- paste0(dot, sprintf('"%s" [label="%s", fillcolor="%s", color="%s"];\n',
                                 id, lbl, fill, border))
    }
  }
  
  if (length(path_taken) > 1) {
    for (i in 1:(length(path_taken)-1)) {
      from <- path_taken[i]; to <- path_taken[i+1]
      ch_node <- dtf_original[[from]]
      elab <- ""
      if (!is.null(ch_node$choices)) {
        nm <- names(ch_node$choices)
        vl <- unlist(ch_node$choices, use.names = FALSE)
        if (to %in% vl) elab <- nm[which(vl == to)[1]]
      }
      elab <- gv_label(elab, 24)
      dot <- paste0(dot, sprintf('"%s" -> "%s" [color="#1E6BB8", fontcolor="#1E6BB8", penwidth=4, label="%s"];\n',
                                 from, to, elab))
    }
  }
  
  paste0(dot, "}")
}

# Convert internal path IDs to a user-friendly step list
path_to_step_list <- function(dtf, path_taken) {
  if (length(path_taken) <= 1) return(data.frame(Step = integer(), Type = character(), Text = character(), Choice = character()))
  out <- data.frame(Step = integer(), Type = character(), Text = character(), Choice = character(), stringsAsFactors = FALSE)
  
  step_num <- 1
  for (i in 1:(length(path_taken)-1)) {
    from_id <- path_taken[i]
    to_id   <- path_taken[i+1]
    from <- dtf[[from_id]]
    to   <- dtf[[to_id]]
    if (is.null(from)) next
    
    choice_label <- ""
    if (!is.null(from$choices)) {
      nm <- names(from$choices)
      vl <- unlist(from$choices, use.names = FALSE)
      if (to_id %in% vl) choice_label <- nm[which(vl == to_id)[1]]
    }
    
    # Record the question step
    if (!is.null(from$type) && from$type == "question") {
      out <- rbind(out, data.frame(
        Step = step_num,
        Type = "Question",
        Text = from$text,
        Choice = choice_label,
        stringsAsFactors = FALSE
      ))
      step_num <- step_num + 1
    }
  }
  
  # Append final recommendation if last node is action
  last_id <- tail(path_taken, 1)
  last <- dtf[[last_id]]
  if (!is.null(last) && last$type == "action") {
    out <- rbind(out, data.frame(
      Step = step_num,
      Type = "Recommendation",
      Text = last$text,
      Choice = "",
      stringsAsFactors = FALSE
    ))
  }
  
  out
}

generate_bullet_list_path <- function(dtf, path_taken) {
  steps <- path_to_step_list(dtf, path_taken)
  if (nrow(steps) == 0) return("No answers have been selected yet.\nTip: click “Demo: completed example run” to see a full example.")
  
  out <- character()
  for (i in seq_len(nrow(steps))) {
    s <- steps[i, ]
    if (s$Type == "Question") {
      out <- c(out, paste0(i, ". ", s$Text, "  →  ", s$Choice))
    } else {
      out <- c(out, paste0("\nFINAL RECOMMENDATION:\n", s$Text))
    }
  }
  paste(out, collapse = "\n")
}

# ------------------------------ EXAMPLE DATA (self-contained) -----------------------------------

# The app will use example_data.csv if present in the app folder; otherwise it uses embedded data.
load_example_dataset <- function() {
  if (file.exists("example_data.csv")) {
    df <- read.csv("example_data.csv", stringsAsFactors = FALSE)
  } else {
    df <- data.frame(
      year = c(2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018),
      landing = c(129,237,471,509,671,816,1111,1154,958,646,710,431,430,532,169)
    )
  }
  
  # Standardize columns
  if (!("year" %in% names(df))) {
    ycol <- if ("years" %in% names(df)) "years" else names(df)[1]
    df$year <- suppressWarnings(as.numeric(df[[ycol]]))
  }
  land_cols <- names(df)[grepl("landing|landings|catch", names(df), ignore.case = TRUE)]
  if (!("landing" %in% names(df)) && length(land_cols) > 0) {
    df$landing <- suppressWarnings(as.numeric(df[[land_cols[1]]]))
  }
  
  df <- df[, c("year", "landing")]
  df <- df[stats::complete.cases(df), , drop = FALSE]
  df <- df[order(df$year), , drop = FALSE]
  df
}

example_df <- load_example_dataset()

# Demo path MUST match the provided example flow (first answer is NO)
demo_path <- c(
  "start", "q_years", "q2", "q3", "q4a", "q4", "q5", "q6", "q7", "q8", "q10", "a_define_management_measures"
)

# ------------------------------ UI ---------------------------------------------------------------

theme_blue <- bs_theme(
  version = 5,
  bootswatch = "flatly",
  primary = "#1E6BB8",
  base_font = font_google("Inter")
)

ui <- fluidPage(
  theme = theme_blue,
  
  tags$head(tags$style(HTML("
    body { background: #F6FAFF; }
    .panel-card { background: white; border-radius: 14px; border: 1px solid rgba(30,107,184,0.16); box-shadow: 0 6px 18px rgba(30,107,184,0.08); padding: 14px; }
    .note-box { background: #E8F2FF; border: 1px solid rgba(30,107,184,0.22); border-radius: 12px; padding: 12px; }
    .success-box { background: #E6FFF2; border: 1px solid rgba(31,122,77,0.30); border-radius: 12px; padding: 12px; }
    .danger-box { background: #FFF2F2; border: 1px solid rgba(192,57,43,0.25); border-radius: 12px; padding: 12px; }
    h2, h3, h4, h5 { color: #123A63; }
    .btn-primary { background-color: #1E6BB8; border-color: #1E6BB8; }
    .btn-warning { border-radius: 10px; }
    .shiny-input-container label { color: #123A63; font-weight: 600; }
    .small-muted { color: rgba(18,58,99,0.75); font-size: 12px; }
    .badge-step { display:inline-block; padding:4px 8px; border-radius: 999px; background:#123A63; color:white; font-weight:600; font-size:12px; }
  "))),
  
  titlePanel("DTF – New Fisheries Decision Tool"),
  
  sidebarLayout(
    sidebarPanel(
      div(class = "panel-card",
          
          div(class = "note-box",
              tags$b("What this tool does: "),
              "Answer the questions on the left. Your pathway is highlighted in the Decision Tree. ",
              "The initial pathway is selected according to whether a new gear is already known and whether an eligible landing series is available.",
              tags$br(),
              tags$span(class = "small-muted", "New users: start with “Demo: completed example run”.")
          ),
          br(),
          
          h3("1) Answer the questionnaire"),
          uiOutput("question_ui"),
          br(),
          
          fluidRow(
            column(6, actionButton("restart", "Restart", class = "btn btn-warning", style = "width:100%;")),
            column(6, actionButton("demo_run", "Demo: completed example run", class = "btn btn-primary", style = "width:100%;"))
          ),
          br(), br(),
          
          h3("2) Landing data (automated pathway)"),
          div(class = "note-box",
              tags$ol(
                tags$li("Use example data to learn the required format."),
                tags$li("Or upload your own CSV with columns: year + landing/catch."),
                tags$li("Automated screening requires at least 6 consecutive annual observations.")
              )
          ),
          br(),
          fluidRow(
            column(6, actionButton("load_example", "Load example data", class = "btn btn-primary", style = "width:100%;")),
            column(6, downloadButton("download_example_csv", "Download example CSV", style = "width:100%;"))
          ),
          br(),
          fileInput("data_file", "Upload your CSV", accept = c(".csv", ".txt")),
          tags$div(class = "small-muted", "Expected columns: 'year' and a landing/catch column (e.g., 'landing' or 'catch')."),
          br(),
          checkboxInput("show_data_preview", "Show data preview", value = FALSE),
          conditionalPanel(condition = "input.show_data_preview",
                           DT::dataTableOutput("data_preview")),
          br(),
          
          h3("3) Exports"),
          div(class = "small-muted", "Exports work at any time. The Network export uses the highlighted pathway."),
          br(),
          h5("Network"),
          downloadButton("download_network_png", "Download Network (PNG, HD)"),
          downloadButton("download_network_pdf", "Download Network (PDF)"),
          br(), br(),
          h5("Analysis figure"),
          downloadButton("export_analysis_png", "Export Analysis (PNG)"),
          downloadButton("export_analysis_pdf", "Export Analysis (PDF)"),
          br(), br(),
          h5("Decision summary (text)"),
          downloadButton("download_path", "Download Bullet Summary (TXT)"),
          downloadButton("download_bullet_pdf", "Download Bullet Summary (PDF)"),
          br(), br(),
          
          checkboxInput("show_table", "Show full decision table (advanced)", value = FALSE)
      )
    ),
    
    mainPanel(
      tabsetPanel(
        
        # ---------------------- HOW TO USE FIRST ----------------------
        tabPanel(
          "How to use",
          h2("How to use this tool"),
          div(class = "note-box",
              tags$b("In short: "),
              "Answer questions on the left. Your pathway is highlighted step-by-step in the Decision Tree.",
              tags$br(),
              tags$span(class = "small-muted", "Tip: click “Demo: completed example run” to see a complete example immediately.")
          ),
          br(),
          h4("Workflow (4 steps)"),
          tags$ol(
            tags$li(tags$b("Answer the questionnaire."), " One question at a time (left panel)."),
            tags$li(tags$b("Watch the Decision Tree."), " The path highlights as you answer."),
            tags$li(tags$b("Initial pathway selection."), " Eligible landing series enter automated screening; shorter or unavailable series enter the expert-supported pathway."),
            tags$li(tags$b("Read the recommendation."), " A final management recommendation is shown at the end.")
          ),
          br(),
          h4("What you will see in the outputs"),
          tags$ul(
            tags$li(tags$b("Decision Tree tab:"), " the pathway is highlighted (default view shows only your pathway for readability)."),
            tags$li(tags$b("Bullet summary tab:"), " a concise summary you can export."),
            tags$li(tags$b("Data Analysis tab:"), " results of the breakpoint test if data are provided.")
          ),
          br(),
          h4("Troubleshooting"),
          div(class = "danger-box",
              tags$ul(
                tags$li(tags$b("I don't know what to click:"), " use the Demo button to see a full example."),
                tags$li(tags$b("Tree is hard to read:"), " keep 'Path only' view in the Decision Tree tab."),
                tags$li(tags$b("Automated analysis does not run:"), " upload at least 6 consecutive annual observations; otherwise use the expert-supported pathway."),
                tags$li(tags$b("CSV not accepted:"), " check it contains a year column and a landing/catch column.")
              )
          )
        ),
        
        # ---------------------- DECISION TREE ----------------------
        tabPanel(
          "Decision Tree",
          fluidRow(
            column(
              8,
              div(class = "note-box",
                  tags$b("Decision Tree view: "),
                  "Use this graph to see your pathway. The default is “Path only” to improve readability.",
                  tags$br(),
                  tags$span(class = "small-muted", "If you want the full framework, switch to “Full tree”.")
              )
            ),
            column(
              4,
              radioButtons(
                "tree_view_mode",
                label = NULL,
                choices = c("Path only (recommended)" = "path", "Full tree" = "full"),
                selected = "path",
                inline = TRUE
              )
            )
          ),
          br(),
          visNetworkOutput("network", height = "720px")
        ),
        
 
        # ---------------------- BULLET SUMMARY ----------------------
        tabPanel(
          "Bullet summary",
          h3("Concise text summary (exportable)"),
          verbatimTextOutput("bullet_list_path"),
          br(),
          downloadButton("download_path2", "Download bullet summary (TXT)")
        ),
        
        # ---------------------- DATA ANALYSIS ----------------------
        tabPanel(
          "Data Analysis",
          conditionalPanel(
            condition = "output.show_analysis_tab",
            h3("Assumption-aware landing breakpoint analysis"),
            div(class = "note-box",
                tags$b("What is reported: "),
                "For every eligible candidate split, the app displays separate Shapiro-Wilk tests ",
                "for the pre- and post-breakpoint segments, the Fligner-Killeen variance test, ",
                "both one-sided Student t-test and Wilcoxon rank-sum p-values, the comparison selected ",
                "by the assumption rule, Cohen's d, percentage increase and the composite DTF score. ",
                "Ljung-Box and lag-1 autocorrelation are reported at series level. ",
                "Minimum requirement: at least 6 consecutive annual observations."
            ),
            br(),
            actionButton(
              "run_analysis_btn",
              "Run / refresh full statistical analysis",
              class = "btn btn-primary"
            ),
            br(), br(),
            verbatimTextOutput("detailed_analysis"),
            br(),
            h4("Selected or highest-ranked candidate: full assumption diagnostics"),
            tags$p(
              class = "small-muted",
              "This table remains visible even when no breakpoint meets the DTF statistical-support criterion."
            ),
            DT::dataTableOutput("selected_breakpoint_diagnostics"),
            br(),
            plotOutput("detailed_plot", height = "420px"),
            br(),
            h4("Series-level temporal diagnostics"),
            DT::dataTableOutput("series_diagnostics"),
            br(),
            h4("All candidate breakpoint diagnostics"),
            tags$p(
              class = "small-muted",
              "Scroll horizontally to inspect normality, variance, Student t-test, Wilcoxon, selected p-value and score for every candidate split."
            ),
            DT::dataTableOutput("candidate_breakpoints"),
            br(),
            h4("Data summary"),
            DT::dataTableOutput("data_summary")
          ),
          conditionalPanel(
            condition = "!output.show_analysis_tab",
            h3("Automated analysis not available"),
            textOutput("analysis_unavailable_message")
          )
        )
      ),
      
      conditionalPanel(
        condition = "input.show_table",
        br(),
        h3("Full decision table (advanced)"),
        DT::dataTableOutput("decision_table")
      )
    )
  )
)

# ------------------------------ SERVER -----------------------------------------------------------

server <- function(input, output, session) {
  
  values <- reactiveValues(
    current_node = "start",
    path_taken = c("start"),
    user_data = NULL,
    analysis_results = NULL,
    breakpoint_result = NULL
  )
  
  reset_all <- function() {
    values$current_node <- "start"
    values$path_taken <- c("start")
    values$user_data <- NULL
    values$analysis_results <- NULL
    values$breakpoint_result <- NULL
  }
  
  observeEvent(input$restart, {
    reset_all()
  })
  
  # Example data load + download
  observeEvent(input$load_example, {
    values$user_data <- example_df
    if (isTRUE(landing_series_status(values$user_data)$eligible)) {
      pipeline <- run_breakpoint_pipeline(values$user_data, alpha = 0.05, min_points = 3)
      values$analysis_results <- pipeline
      values$breakpoint_result <- pipeline$selected_breakpoint
      showNotification(
        "Example landing data loaded and the full assumption-aware analysis was run.",
        type = "message"
      )
    } else {
      values$analysis_results <- NULL
      values$breakpoint_result <- NULL
      showNotification("Example landing data loaded.", type = "message")
    }
  })
  
  output$download_example_csv <- downloadHandler(
    filename = function() paste0("example_landing_data_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) write.csv(example_df, file, row.names = FALSE)
  )
  
  # Demo run: completed example pathway that matches the provided figure
  observeEvent(input$demo_run, {
    values$user_data <- example_df
    values$path_taken <- demo_path
    values$current_node <- tail(demo_path, 1)
    
    # Compute the full manuscript-aligned analysis for the demonstration.
    if (isTRUE(landing_series_status(values$user_data)$eligible)) {
      pipeline <- run_breakpoint_pipeline(values$user_data, alpha = 0.05, min_points = 3)
      values$analysis_results <- pipeline
      values$breakpoint_result <- pipeline$selected_breakpoint
    } else {
      values$analysis_results <- NULL
      values$breakpoint_result <- NULL
    }
    
    showNotification("Demo loaded: full example pathway (starting with 'No' at the first question).", type = "message")
  })
  
  # Data upload
  observeEvent(input$data_file, {
    req(input$data_file)
    tryCatch({
      df <- read.csv(input$data_file$datapath, stringsAsFactors = FALSE)
      
      if (!("year" %in% names(df) || "years" %in% names(df))) {
        showNotification("Upload error: the file must contain a 'year' column (or 'years').", type = "error")
        return()
      }
      ycol <- if ("year" %in% names(df)) "year" else "years"
      
      land_cols <- names(df)[grepl("landing|landings|catch", names(df), ignore.case = TRUE)]
      if (!length(land_cols)) {
        showNotification("Upload error: no landing/catch column found. Try naming it 'landing' or 'catch'.", type = "error")
        return()
      }
      lcol <- land_cols[1]
      
      dfc <- data.frame(
        year = suppressWarnings(as.numeric(df[[ycol]])),
        landing = suppressWarnings(as.numeric(df[[lcol]]))
      )
      dfc <- dfc[stats::complete.cases(dfc), , drop = FALSE]
      dfc <- dfc[order(dfc$year), , drop = FALSE]
      
      if (nrow(dfc) < 2) {
        showNotification("Upload error: at least 2 valid (year, landing) rows are required.", type = "error")
        return()
      }
      
      values$user_data <- dfc
      status <- landing_series_status(dfc)

      if (isTRUE(status$eligible)) {
        pipeline <- run_breakpoint_pipeline(dfc, alpha = 0.05, min_points = 3)
        values$analysis_results <- pipeline
        values$breakpoint_result <- pipeline$selected_breakpoint
      } else {
        values$analysis_results <- NULL
        values$breakpoint_result <- NULL
      }

      showNotification(
        paste(
          "Landing data uploaded successfully.",
          status$message,
          if (isTRUE(status$eligible)) "Full assumption-aware diagnostics are now available in the Data Analysis tab." else ""
        ),
        type = if (isTRUE(status$eligible)) "message" else "warning",
        duration = 10
      )
    }, error = function(e) {
      showNotification(paste("Upload error:", e$message), type = "error")
    })
  })
  
  # Data preview
  output$data_preview <- DT::renderDataTable({
    req(values$user_data)
    DT::datatable(values$user_data, options = list(pageLength = 10, scrollX = TRUE))
  })
  
  # Run the manuscript-aligned assumption-aware analysis only for an eligible landing series.
  perform_breakpoint_analysis <- function() {
    req(values$user_data)
    df <- values$user_data
    status <- landing_series_status(df)

    if (!isTRUE(status$eligible)) {
      values$analysis_results <- NULL
      values$breakpoint_result <- NULL
      showNotification(
        paste("Automated analysis not performed.", status$message,
              "Proceeding through the expert-supported pathway is recommended."),
        type = "warning",
        duration = 8
      )
      return(FALSE)
    }

    pipeline <- run_breakpoint_pipeline(df, alpha = 0.05, min_points = 3)
    values$analysis_results <- pipeline
    values$breakpoint_result <- pipeline$selected_breakpoint
    TRUE
  }

  observeEvent(input$run_analysis_btn, {
    if (is.null(values$user_data)) {
      showNotification("Upload landing data before running the analysis.", type = "warning")
      return()
    }

    if (perform_breakpoint_analysis()) {
      showNotification(
        "Full statistical analysis completed. Normality, variance, Student t-test, Wilcoxon and temporal diagnostics have been updated.",
        type = "message",
        duration = 8
      )
    }
  }, ignoreInit = TRUE)

  # Questionnaire UI (clear wording + guidance)
  output$question_ui <- renderUI({
    node <- dtf_original[[values$current_node]]
    if (is.null(node)) return(div(h4("End of decision tree"), style = "color: #C0392B;"))
    
    # Special case: select the initial pathway from the uploaded landing series when available.
    if (values$current_node == "q_years" && !is.null(values$user_data)) {
      status <- landing_series_status(values$user_data)
      next_node <- if (isTRUE(status$eligible)) "q2" else "q4a"
      selected_answer <- if (isTRUE(status$eligible)) "Yes" else "No"
      values$current_node <- next_node
      values$path_taken <- c(values$path_taken, next_node)
      showNotification(
        paste0("Landing-series check: ", status$message,
               " The app selected the '", selected_answer, "' branch."),
        type = if (isTRUE(status$eligible)) "message" else "warning",
        duration = 8
      )
      return(div(
        div(class = if (isTRUE(status$eligible)) "success-box" else "note-box",
            tags$b("Initial pathway selected. "),
            status$message)
      ))
    }
    
    # Special case: Q2 requires eligible landing data and is answered automatically.
    if (values$current_node == "q2") {
      if (is.null(values$user_data)) {
        return(div(
          div(class = "danger-box",
              tags$b("Automated pathway requires landing data: "),
              "Upload an eligible landing series (or load the example data) to run the automated screening.",
              tags$br(),
              tags$span(class = "small-muted", "If an eligible series is unavailable, continue through the expert-supported pathway.")
          ),
          br(),
          h4(node$text),
          actionButton("expert_path_q2_btn", "Use expert-supported pathway", class = "btn btn-primary", style = "width: 100%;")
        ))
      } else {
        ok <- perform_breakpoint_analysis()
        if (!ok) {
          return(div(
            div(class = "danger-box",
                tags$b("Landing series not eligible for automated screening: "),
                landing_series_status(values$user_data)$message,
                tags$br(),
                tags$span(class = "small-muted", "Continue through the expert-supported pathway or upload an eligible series.")
            ),
            br(),
            h4(node$text),
            actionButton("expert_path_q2_btn", "Use expert-supported pathway", class = "btn btn-primary", style = "width: 100%;")
          ))
        }
        
        # Auto-route based on breakpoint existence
        if (!is.null(values$breakpoint_result)) {
          values$current_node <- "q3"
          values$path_taken <- c(values$path_taken, "q3")
          showNotification("Automated screening identified a diagnostic increase → proceeding on the 'Yes' branch.", type = "message")
        } else {
          values$current_node <- "a_continue_monitoring"
          values$path_taken <- c(values$path_taken, "a_continue_monitoring")
          showNotification("Automated screening did not identify a clear increase → proceeding on the 'No' branch.", type = "message")
        }
        
        return(div(
          div(class = "success-box",
              tags$b("Automated screening completed. "),
              "The tool selected the next step based on the diagnostic result."
          )
        ))
      }
    }
    
    if (node$type == "question") {
      btns <- lapply(names(node$choices), function(ch) {
        actionButton(
          paste0("choice_", ch),
          ch,
          class = "btn btn-primary",
          style = "margin: 4px 0; width: 100%;"
        )
      })
      
      return(div(
        tags$span(class = "badge-step", "Current question"),
        br(), br(),
        h4(node$text),
        do.call(div, btns)
      ))
    }
    
    if (node$type == "action") {
      return(div(
        div(class = "success-box",
            tags$b("Final recommendation"),
            tags$p(node$text, style = "font-weight: 700; margin-top: 8px;")
        ),
        br(),
        tags$div(class = "note-box",
                 tags$b("Next: "),
                 "You can export the highlighted network and the bullet summary from the left panel."
        )
      ))
    }
    
    div(p("Unknown node type."))
  })
  
  # Choice handling (Yes/No)
  observeEvent(input$choice_Yes, {
    node <- dtf_original[[values$current_node]]
    req(node, node$type == "question", "Yes" %in% names(node$choices))
    to <- node$choices[["Yes"]]
    values$current_node <- to
    values$path_taken <- c(values$path_taken, to)
  }, ignoreInit = TRUE)
  
  observeEvent(input$choice_No, {
    node <- dtf_original[[values$current_node]]
    req(node, node$type == "question", "No" %in% names(node$choices))
    to <- node$choices[["No"]]
    values$current_node <- to
    values$path_taken <- c(values$path_taken, to)
  }, ignoreInit = TRUE)
  
  # Q2: when eligible landing data are unavailable, use the expert-supported pathway.
  observeEvent(input$expert_path_q2_btn, {
    if (values$current_node == "q2") {
      values$current_node <- "q4a"
      values$path_taken <- c(values$path_taken, "q4a")
      showNotification(
        "Proceeding through the expert-supported pathway because automated screening is unavailable.",
        type = "warning"
      )
    }
  }, ignoreInit = TRUE)
  
  # Other choices (non-Yes/No) – dynamic observers
  observe({
    node <- dtf_original[[values$current_node]]
    req(node, node$type == "question", !is.null(node$choices))
    lapply(names(node$choices), function(ch) {
      if (!ch %in% c("Yes", "No")) {
        tgt <- node$choices[[ch]]
        observeEvent(input[[paste0("choice_", ch)]], {
          values$current_node <- tgt
          values$path_taken <- c(values$path_taken, tgt)
        }, ignoreInit = TRUE, once = TRUE)
      }
    })
  })
  
  # Analysis tab visibility: quantitative outputs are shown only for eligible series.
  output$show_analysis_tab <- reactive({
    !is.null(values$user_data) && isTRUE(landing_series_status(values$user_data)$eligible)
  })
  outputOptions(output, "show_analysis_tab", suspendWhenHidden = FALSE)
  
  output$analysis_unavailable_message <- renderText({
    if (is.null(values$user_data)) {
      return("Upload landing data (or load the example) to assess eligibility for automated screening. The DTF can still be used through the expert-supported pathway when landing data are unavailable.")
    }
    paste0(
      landing_series_status(values$user_data)$message,
      " The DTF remains operational through the expert-supported pathway based on fishing-effort data, reliable proxies and expert knowledge."
    )
  })
  
  # ---------------------- Decision Tree (readable; default path-only) ----------------------
  
  output$network <- renderVisNetwork({
    mode <- input$tree_view_mode
    dtf_to_show <- if (identical(mode, "path")) {
      filter_dtf_by_path(dtf_original, values$path_taken)
    } else {
      dtf_original
    }
    
    nodes <- data.frame(id = names(dtf_to_show), stringsAsFactors = FALSE)
    nodes$label <- sapply(names(dtf_to_show), function(x) {
      nd <- dtf_to_show[[x]]
      txt <- if (nd$type == "question") nd$text else paste("RECOMMENDATION:", nd$text)
      create_wrapped_label(txt, 34)
    })
    nodes$shape <- "box"
    
    is_action <- sapply(dtf_to_show, function(x) x$type == "action")
    nodes$color <- ifelse(is_action, "#33C27F", "#1E6BB8")
    
    # Highlight nodes on the taken path
    in_path <- nodes$id %in% values$path_taken
    nodes$color[in_path & !is_action] <- "#123A63"
    nodes$color[in_path & is_action]  <- "#1F7A4D"
    
    nodes$font.size <- if (identical(mode, "path")) 22 else 18
    nodes$margin <- if (identical(mode, "path")) 18 else 12
    
    # Build edges (only those present in dtf_to_show)
    # Build edges
    if (identical(mode, "path")) {
      
      # ONLY edges along the selected path
      edges <- data.frame(from = character(), to = character(), label = character(),
                          color = character(), width = numeric(),
                          stringsAsFactors = FALSE)
      
      if (length(values$path_taken) > 1) {
        for (i in 1:(length(values$path_taken) - 1)) {
          from_id <- values$path_taken[i]
          to_id   <- values$path_taken[i + 1]
          
          # label = the chosen answer text
          nd <- dtf_original[[from_id]]
          lab <- ""
          if (!is.null(nd$choices)) {
            nm <- names(nd$choices)
            vl <- unlist(nd$choices, use.names = FALSE)
            if (to_id %in% vl) lab <- nm[which(vl == to_id)[1]]
          }
          
          edges <- rbind(edges, data.frame(
            from = from_id,
            to   = to_id,
            label = create_wrapped_label(lab, 22),
            color = "#000000",
            width = 3.5,
            stringsAsFactors = FALSE
          ))
        }
      }
      
    } else {
      
      # FULL tree: show all edges
      edges <- data.frame(from = character(), to = character(), label = character(),
                          color = character(), width = numeric(),
                          stringsAsFactors = FALSE)
      
      for (node_id in names(dtf_to_show)) {
        nd <- dtf_to_show[[node_id]]
        if (nd$type == "question" && !is.null(nd$choices)) {
          for (choice_text in names(nd$choices)) {
            to_id <- nd$choices[[choice_text]]
            if (!to_id %in% names(dtf_to_show)) next
            
            edges <- rbind(edges, data.frame(
              from = node_id,
              to   = to_id,
              label = create_wrapped_label(choice_text, 22),
              color = "rgba(90,110,130,0.30)",
              width = 1,
              stringsAsFactors = FALSE
            ))
          }
        }
      }
    }
    
    
    visNetwork(nodes, edges) %>%
      visOptions(
        highlightNearest = TRUE,
        nodesIdSelection = list(
          enabled = TRUE,
          useLabels = TRUE,
          main = "Jump to node"
        )
      ) %>%
      visLayout(
        randomSeed = 123,
        improvedLayout = TRUE,
        hierarchical = list(
          enabled = TRUE,
          levelSeparation = if (identical(mode, "path")) 320 else 230,
          nodeSpacing = if (identical(mode, "path")) 260 else 220,
          treeSpacing = if (identical(mode, "path")) 300 else 260,
          direction = "UD",
          sortMethod = "directed"
        )
      ) %>%
      visPhysics(enabled = FALSE) %>%
      visInteraction(dragNodes = TRUE, dragView = TRUE, zoomView = TRUE) %>%
      visNodes(
        shadow = FALSE,
        borderWidth = 2,
        font = list(multi = "md", color = "white"),
        color = list(highlight = "#000000")
      ) %>%
      visEdges(
        arrows = "to",
        smooth = list(enabled = TRUE, type = "cubicBezier", roundness = 0.35),
        font = list(
          size = 18,
          align = "horizontal",
          vadjust = -20
        ),
        chosen = FALSE,
        labelHighlightBold = FALSE
      )
    
    
  })
  

  # Bullet summary text
  clean_text <- function(x) {
    x <- gsub("\u2192", " -> ", x)      # Unicode arrow
    x <- gsub("<U\\+2192>", " -> ", x)  # Literal <U+2192>
    x
  }
  
  output$bullet_list_path <- renderText({
    txt <- clean_text(generate_bullet_list_path(dtf_original, values$path_taken))
    clean_text(txt)
  })
  
  
  # Bullet downloads
  output$download_path <- downloadHandler(
    filename = function() paste0("decision_summary_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"),
    content = function(file) {
      txt <- generate_bullet_list_path(dtf_original, values$path_taken)
      writeLines(clean_text(txt), file)
    }
  )
  
  output$download_path2 <- downloadHandler(
    filename = function() paste0(bullet_filename(), ".txt"),
    content = function(file){
      txt <- generate_bullet_list_path(dtf_original, values$path_taken)
      writeLines(clean_text(txt), file)
    }
  )
  output$download_bullet_pdf <- downloadHandler(
    filename = function() paste0(bullet_filename(), ".pdf"),
    content = function(file) {
      txt <- generate_bullet_list_path(dtf_original, values$path_taken)
      pdf(file, width = 8.27, height = 11.69) # A4
      grid::grid.newpage()
      g <- grid::textGrob(txt, x = 0.02, y = 0.98, just = c("left","top"),
                          gp = grid::gpar(fontsize = 12), default.units = "npc")
      grid::grid.draw(g)
      dev.off()
    }
  )
  
  # ---------------------- Data Analysis outputs ----------------------

  output$detailed_analysis <- renderText({
    req(values$user_data)

    status <- landing_series_status(values$user_data)
    if (!isTRUE(status$eligible)) {
      return(paste0(
        "Analysis not performed.\n",
        "Reason: ", status$message, "\n",
        "Use the expert-supported pathway to continue the DTF."
      ))
    }

    pipeline <- values$analysis_results
    if (is.null(pipeline)) {
      return(
        "No analysis has been run yet. Click 'Run / refresh full statistical analysis' or reach the automated landing-screening node."
      )
    }

    temporal <- pipeline$temporal_diagnostics
    selected <- pipeline$selected_breakpoint
    displayed <- if (is.null(selected)) pipeline$highest_ranked_candidate else selected
    displayed_label <- if (is.null(selected)) {
      "Highest-ranked candidate (not supported as a DTF signal)"
    } else {
      "DTF-selected diagnostic signal"
    }

    temporal_text <- paste0(
      "SERIES-LEVEL TEMPORAL DIAGNOSTICS\n",
      "Ljung-Box lag: ",
      if (is.null(temporal) || is.na(temporal$ljung_box_lag)) "not evaluated" else temporal$ljung_box_lag,
      "\n",
      "Ljung-Box p-value: ",
      if (is.null(temporal)) "not evaluated" else format_p_value(temporal$ljung_box_p),
      "\n",
      "Lag-1 autocorrelation: ",
      if (is.null(temporal) || is.na(temporal$acf_lag1)) "not evaluated" else round(temporal$acf_lag1, 3),
      "\n",
      "Temporal-autocorrelation warning: ",
      if (is.null(temporal)) "not evaluated" else ifelse(temporal$autocorrelation_warning, "YES", "NO"),
      "\n\n"
    )

    small_sample_note <- paste0(
      "CAUTION\n",
      "Candidate splits may contain only three observations per segment. ",
      "Shapiro-Wilk and Fligner-Killeen tests have low power in very small samples; ",
      "the tests and the selected p-value are diagnostic early-warning evidence, not confirmatory inference.\n\n"
    )

    if (is.null(displayed)) {
      return(paste0(
        "Analysis method: assumption-aware DTF breakpoint screening\n\n",
        temporal_text,
        small_sample_note,
        "RESULT\n",
        "No positive finite candidate breakpoint score was available."
      ))
    }

    paste0(
      "Analysis method: assumption-aware DTF breakpoint screening\n\n",
      temporal_text,
      small_sample_note,
      "CANDIDATE DISPLAYED\n",
      displayed_label, "\n",
      "Breakpoint year: ", displayed$year, "\n",
      "n before / n after: ", displayed$n_before, " / ", displayed$n_after, "\n\n",
      "ASSUMPTION TESTS\n",
      "Shapiro-Wilk p-value before breakpoint: ", format_p_value(displayed$shapiro_before_p), "\n",
      "Shapiro-Wilk p-value after breakpoint: ", format_p_value(displayed$shapiro_after_p), "\n",
      "Normality supported: ", ifelse(displayed$normality_supported, "YES", "NO"), "\n",
      "Fligner-Killeen p-value: ", format_p_value(displayed$fligner_p), "\n",
      "Homoscedasticity supported: ", ifelse(displayed$homoscedasticity_supported, "YES", "NO"), "\n\n",
      "ONE-SIDED COMPARISONS\n",
      "Student t-test p-value: ", format_p_value(displayed$student_t_p), "\n",
      "Wilcoxon rank-sum p-value: ", format_p_value(displayed$wilcoxon_p), "\n",
      "Selected comparison: ", displayed$selected_test, "\n",
      "Selected p-value: ", format_p_value(displayed$selected_p), "\n\n",
      "MAGNITUDE AND DTF RANKING\n",
      "Mean before: ", round(displayed$mean_before, 3), "\n",
      "Mean after: ", round(displayed$mean_after, 3), "\n",
      "Increase in mean landings: ", round(displayed$percent_increase, 2), "%\n",
      "Cohen's d: ", if (is.finite(displayed$cohens_d)) round(displayed$cohens_d, 3) else "not evaluated", "\n",
      "Composite score: ", if (is.finite(displayed$score)) round(displayed$score, 4) else "not evaluated", "\n",
      "Statistical support at alpha = ", pipeline$alpha, ": ",
      ifelse(isTRUE(displayed$supported), "YES", "NO"), "\n\n",
      "DTF ROUTING\n",
      if (isTRUE(displayed$supported)) {
        "The candidate is used as the positive automated-screening signal."
      } else {
        "The candidate is displayed for transparency but does not route the DTF to the positive-signal branch."
      }
    )
  })

  output$detailed_plot <- renderPlot({
    req(values$user_data)
    req(isTRUE(landing_series_status(values$user_data)$eligible))
    plot_breakpoint(values$user_data, values$breakpoint_result)
  })

  output$selected_breakpoint_diagnostics <- DT::renderDataTable({
    req(values$user_data)
    req(isTRUE(landing_series_status(values$user_data)$eligible))

    DT::datatable(
      selected_candidate_table(values$analysis_results),
      options = list(pageLength = 25, dom = "t", scrollX = TRUE),
      rownames = FALSE
    )
  })

  output$series_diagnostics <- DT::renderDataTable({
    req(values$user_data)
    req(isTRUE(landing_series_status(values$user_data)$eligible))

    temporal <- if (is.null(values$analysis_results)) {
      series_temporal_diagnostics(values$user_data)
    } else {
      values$analysis_results$temporal_diagnostics
    }

    series_df <- data.frame(
      Diagnostic = c(
        "Number of annual observations",
        "Ljung-Box lag",
        "Ljung-Box p-value",
        "Lag-1 autocorrelation",
        "Temporal-autocorrelation warning"
      ),
      Result = c(
        nrow(values$user_data),
        if (is.na(temporal$ljung_box_lag)) "not evaluated" else temporal$ljung_box_lag,
        format_p_value(temporal$ljung_box_p),
        if (is.na(temporal$acf_lag1)) "not evaluated" else round(temporal$acf_lag1, 3),
        ifelse(temporal$autocorrelation_warning, "YES", "NO")
      ),
      Interpretation = c(
        "Automated screening requires at least six consecutive annual observations.",
        "Lag used in the Ljung-Box test.",
        if (!is.na(temporal$ljung_box_p) && temporal$ljung_box_p < 0.05) {
          "Evidence of temporal dependence."
        } else {
          "No significant Ljung-Box evidence of temporal dependence at alpha = 0.05."
        },
        "Values with larger absolute magnitude indicate stronger lag-1 dependence.",
        "Triggered when Ljung-Box p < 0.05 or |lag-1 autocorrelation| >= 0.40."
      ),
      stringsAsFactors = FALSE
    )

    DT::datatable(
      series_df,
      options = list(pageLength = 10, dom = "t", scrollX = TRUE),
      rownames = FALSE
    )
  })

  output$candidate_breakpoints <- DT::renderDataTable({
    req(values$user_data)
    req(isTRUE(landing_series_status(values$user_data)$eligible))

    if (is.null(values$analysis_results)) {
      return(DT::datatable(
        data.frame(Message = "Run the automated landing-screening node to generate candidate diagnostics."),
        options = list(dom = "t"),
        rownames = FALSE
      ))
    }

    DT::datatable(
      candidate_table_for_display(values$analysis_results),
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = FALSE
    )
  })

  output$data_summary <- DT::renderDataTable({
    req(values$user_data)
    req(isTRUE(landing_series_status(values$user_data)$eligible))
    df <- values$user_data
    summary_df <- data.frame(
      Item = c(
        "Number of years", "First year", "Last year",
        "Years consecutive", "Duplicated years",
        "Min landing", "Max landing", "Mean landing", "Median landing"
      ),
      Value = c(
        nrow(df), min(df$year), max(df$year),
        ifelse(all(diff(sort(df$year)) == 1), "Yes", "No"),
        ifelse(anyDuplicated(df$year) > 0, "Yes", "No"),
        round(min(df$landing), 2), round(max(df$landing), 2),
        round(mean(df$landing), 2), round(median(df$landing), 2)
      ),
      stringsAsFactors = FALSE
    )
    DT::datatable(summary_df, options = list(pageLength = 20, dom = "t"), rownames = FALSE)
  })

  output$decision_table <- DT::renderDataTable({
    dt <- data.frame(
      ID   = names(dtf_original),
      Type = sapply(dtf_original, function(x) x$type),
      Text = sapply(dtf_original, function(x) x$text),
      stringsAsFactors = FALSE
    )
    DT::datatable(dt, options = list(pageLength = 15, scrollX = TRUE))
  })
  
  # ---------------------- Exports: Network (Graphviz path-only) ----------------------
  
  export_graphviz_svg <- function(path_taken) {
    filtered <- filter_dtf_by_path(dtf_original, path_taken)
    dot <- build_graphviz_dot(filtered, path_taken)
    gr <- DiagrammeR::grViz(dot)
    DiagrammeRsvg::export_svg(gr)
  }
  
  output$download_network_png <- downloadHandler(
    filename = function() paste0("decision_tree_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) {
      svg_txt <- export_graphviz_svg(values$path_taken)
      rsvg::rsvg_png(charToRaw(svg_txt), file = file, width = 7000)
    }
  )
  
  output$download_network_pdf <- downloadHandler(
    filename = function() paste0("decision_tree_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
    content = function(file) {
      svg_txt <- export_graphviz_svg(values$path_taken)
      rsvg::rsvg_pdf(charToRaw(svg_txt), file = file)
    }
  )
  
  # Export analysis figure
  output$export_analysis_png <- downloadHandler(
    filename = function() paste0("breakpoint_analysis_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) {
      req(values$user_data)
      req(isTRUE(landing_series_status(values$user_data)$eligible))
      g <- build_analysis_grob(values$user_data, values$breakpoint_result, values$analysis_results)
      png(file, width = 4800, height = 2400, res = 300)
      grid::grid.newpage(); grid::grid.draw(g); dev.off()
    }
  )
  
  output$export_analysis_pdf <- downloadHandler(
    filename = function() paste0("breakpoint_analysis_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
    content = function(file) {
      req(values$user_data)
      req(isTRUE(landing_series_status(values$user_data)$eligible))
      g <- build_analysis_grob(values$user_data, values$breakpoint_result, values$analysis_results)
      pdf(file, width = 16, height = 8, family = "Helvetica")
      grid::grid.newpage(); grid::grid.draw(g); dev.off()
    }
  )
}

shinyApp(ui = ui, server = server)
