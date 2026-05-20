#' Run random-effects meta-analysis on CPMM slope-change estimates
#'
#' Performs random-effects meta-analysis for each protein using cohort-specific
#' CPMM slope-change estimates.
#'
#' The input should be a combined CPMM results data frame containing one row per
#' protein per cohort. At minimum, the data frame must contain columns for
#' protein name, Delta, and SE Delta.
#'
#' @param all_cohort_results_df Combined CPMM results data frame.
#' @param protein_col Name of the protein column. Default is "Protein".
#' @param delta_col Name of the CPMM slope-change estimate column. Default is "Delta".
#' @param se_delta_col Name of the standard error column for Delta. Default is "SE Delta".
#' @param min_cohorts Minimum number of valid cohorts required per protein.
#' Default is 3.
#' @param fdr_method Multiple-testing correction method passed to p.adjust().
#' Default is "BH".
#' @param alpha FDR significance threshold. Default is 0.05.
#'
#' @return A tibble with random-effects meta-analysis results per protein.
#' @export
run_cpmm_meta_analysis <- function(
    all_cohort_results_df,
    protein_col = "Protein",
    delta_col = "Delta",
    se_delta_col = "SE Delta",
    min_cohorts = 3,
    fdr_method = "BH",
    alpha = 0.05
) {

  ## ---- dependencies ----
  pkgs <- c("metafor", "tibble", "dplyr")

  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required.", p), call. = FALSE)
    }
  }

  ## ---- input checks ----
  if (!is.data.frame(all_cohort_results_df)) {
    stop("all_cohort_results_df must be a data frame.", call. = FALSE)
  }

  required_cols <- c(protein_col, delta_col, se_delta_col)

  missing_cols <- setdiff(required_cols, names(all_cohort_results_df))

  if (length(missing_cols) > 0) {
    stop(
      sprintf(
        "The following required columns were not found: %s",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (!is.numeric(all_cohort_results_df[[delta_col]])) {
    stop(sprintf("Column '%s' must be numeric.", delta_col), call. = FALSE)
  }

  if (!is.numeric(all_cohort_results_df[[se_delta_col]])) {
    stop(sprintf("Column '%s' must be numeric.", se_delta_col), call. = FALSE)
  }

  proteins <- unique(all_cohort_results_df[[protein_col]])

  meta_results <- vector("list", length(proteins))

  ## ---- protein-wise random-effects meta-analysis ----
  for (i in seq_along(proteins)) {

    protein <- proteins[[i]]

    subdf <- all_cohort_results_df[
      all_cohort_results_df[[protein_col]] == protein,
      ,
      drop = FALSE
    ]

    subdf <- subdf[
      !is.na(subdf[[delta_col]]) &
        !is.na(subdf[[se_delta_col]]) &
        subdf[[se_delta_col]] > 0,
      ,
      drop = FALSE
    ]

    if (nrow(subdf) < min_cohorts) {
      meta_results[[i]] <- NULL
      next
    }

    yi <- subdf[[delta_col]]
    vi <- subdf[[se_delta_col]]^2

    fit <- tryCatch(
      metafor::rma.uni(
        yi = yi,
        vi = vi,
        method = "PM",
        test = "z"
      ),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      meta_results[[i]] <- NULL
      next
    }

    meta_delta <- as.numeric(fit$b[1])
    meta_se <- as.numeric(fit$se[1])

    meta_z <- if (!is.na(meta_se) && meta_se > 0) {
      meta_delta / meta_se
    } else {
      NA_real_
    }

    meta_p <- if (!is.na(meta_z)) {
      2 * stats::pnorm(abs(meta_z), lower.tail = FALSE)
    } else {
      NA_real_
    }

    ci_low <- meta_delta - stats::qnorm(0.975) * meta_se
    ci_high <- meta_delta + stats::qnorm(0.975) * meta_se

    meta_results[[i]] <- tibble::tibble(
      Protein = protein,
      `N Cohorts` = nrow(subdf),

      `Meta Delta` = meta_delta,
      `Meta SE` = meta_se,
      `Meta 95% CI Lower` = ci_low,
      `Meta 95% CI Upper` = ci_high,
      `Meta Z` = meta_z,
      `Meta P-value` = meta_p,

      Q = as.numeric(fit$QE),
      df = as.integer(nrow(subdf) - 1),
      `Tau^2` = as.numeric(fit$tau2),
      `I^2` = as.numeric(fit$I2)
    )
  }

  meta_results_df <- dplyr::bind_rows(meta_results)

  if (nrow(meta_results_df) == 0) {
    return(
      tibble::tibble(
        Protein = character(),
        `N Cohorts` = integer(),
        `Meta Delta` = numeric(),
        `Meta SE` = numeric(),
        `Meta 95% CI Lower` = numeric(),
        `Meta 95% CI Upper` = numeric(),
        `Meta Z` = numeric(),
        `Meta P-value` = numeric(),
        Q = numeric(),
        df = integer(),
        `Tau^2` = numeric(),
        `I^2` = numeric(),
        `Meta Adjusted P-value (FDR)` = numeric(),
        `Meta Significant` = logical()
      )
    )
  }

  meta_results_df <- meta_results_df[
    !is.na(meta_results_df$`Meta P-value`),
    ,
    drop = FALSE
  ]

  meta_results_df$`Meta Adjusted P-value (FDR)` <- stats::p.adjust(
    meta_results_df$`Meta P-value`,
    method = fdr_method
  )

  meta_results_df$`Meta Significant` <-
    meta_results_df$`Meta Adjusted P-value (FDR)` < alpha

  meta_results_df <- meta_results_df |>
    dplyr::arrange(
      .data$`Meta Adjusted P-value (FDR)`,
      .data$`Meta P-value`
    )

  tibble::as_tibble(meta_results_df)
}
