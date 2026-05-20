#' Fit change-point mixed models (CPMM) across proteins
#'
#' Fits a random-intercept CPMM for each protein in status-change subjects:
#'
#'   protein ~ before_onset + after_onset + covariates + (1 | subject_id)
#'
#' The change point is fixed at 0 on the years-since-onset scale.
#'
#' @details
#' Piecewise predictors are defined as:
#'   before_onset = pmax(0, 0 - years_since_onset)
#'   after_onset  = pmax(0, years_since_onset - 0)
#'
#' The raw coefficient for before_onset is sign-flipped to produce the
#' interpretable pre-onset slope:
#'
#'   Beta 1 = -Beta Before Raw
#'   Beta 2 =  Beta After Raw
#'
#' The slope-change contrast is:
#'
#'   Delta = Beta After Raw + Beta Before Raw
#'
#' which is equivalent to:
#'
#'   Delta = Beta 2 - Beta 1
#'
#' @param df_status_change Data frame containing status-change subjects only.
#' @param protein_list Character vector of protein column names.
#' @param covariates Fixed-effect covariates. Default is c("SEX", "BASELINE_AGE").
#' @param subject_id_col Subject ID column. Default is "SUBID".
#' @param years_since_onset_col Years-since-onset column. Default is "years_since_onset".
#' @param change_point Numeric change point. Default is 0.
#' @param reml Logical; whether to fit models using REML. Default is TRUE.
#'
#' @return A tibble with one row per protein containing CPMM coefficients,
#' slope-change contrast estimates, confidence intervals, and model-fit metrics.
#' @export
fit_cpmm_all_proteins <- function(
    df_status_change,
    protein_list,
    covariates = c("SEX", "BASELINE_AGE"),
    subject_id_col = "SUBID",
    years_since_onset_col = "years_since_onset",
    change_point = 0,
    reml = TRUE
) {

  ## ---- dependencies ----
  pkgs <- c("lme4", "lmerTest", "tibble", "dplyr")

  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required.", p), call. = FALSE)
    }
  }

  ## ---- input checks ----
  if (!is.data.frame(df_status_change)) {
    stop("df_status_change must be a data frame.", call. = FALSE)
  }

  if (!(subject_id_col %in% names(df_status_change))) {
    stop(sprintf("Subject ID column '%s' was not found.", subject_id_col), call. = FALSE)
  }

  if (!(years_since_onset_col %in% names(df_status_change))) {
    stop(sprintf("Years-since-onset column '%s' was not found.", years_since_onset_col), call. = FALSE)
  }

  missing_covariates <- setdiff(covariates, names(df_status_change))
  if (length(missing_covariates) > 0) {
    stop(
      sprintf(
        "The following covariates were not found in the data: %s",
        paste(missing_covariates, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  ## ---- helper to safely quote variable names ----
  quote_var <- function(x) {
    paste0("`", gsub("`", "\\\\`", x), "`")
  }

  ## ---- add piecewise terms ----
  df <- df_status_change

  yso <- df[[years_since_onset_col]]

  df$before_onset <- pmax(0, change_point - yso)
  df$after_onset  <- pmax(0, yso - change_point)

  out <- vector("list", length(protein_list))

  ## ---- model loop ----
  for (i in seq_along(protein_list)) {

    protein <- protein_list[[i]]

    empty_row <- tibble::tibble(
      Protein = protein,

      `Beta 1` = NA_real_,
      `SE Beta 1` = NA_real_,
      `Beta 2` = NA_real_,
      `SE Beta 2` = NA_real_,

      `Beta Before Raw` = NA_real_,
      `Beta After Raw` = NA_real_,
      `Var Before Raw` = NA_real_,
      `Var After Raw` = NA_real_,
      `Cov Before After` = NA_real_,

      Delta = NA_real_,
      `SE Delta` = NA_real_,
      `Delta 95% CI Lower` = NA_real_,
      `Delta 95% CI Upper` = NA_real_,

      Intercept = NA_real_,
      AIC = NA_real_,
      BIC = NA_real_,
      MSE = NA_real_,
      `N Obs` = NA_integer_
    )

    if (!(protein %in% names(df))) {
      out[[i]] <- empty_row
      next
    }

    rhs <- paste(
      c("before_onset", "after_onset", covariates),
      collapse = " + "
    )

    fml <- stats::as.formula(
      paste0(
        quote_var(protein),
        " ~ ",
        rhs,
        " + (1 | ",
        quote_var(subject_id_col),
        ")"
      )
    )

    fit <- tryCatch(
      suppressWarnings(
        lmerTest::lmer(
          formula = fml,
          data = df,
          REML = reml,
          na.action = stats::na.omit
        )
      ),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      out[[i]] <- empty_row
      next
    }

    coefs <- summary(fit)$coefficients
    vc <- as.matrix(stats::vcov(fit))

    ## ---- raw fixed-effect coefficients ----
    beta_before_raw <- if ("before_onset" %in% rownames(coefs)) {
      coefs["before_onset", "Estimate"]
    } else {
      NA_real_
    }

    beta_after_raw <- if ("after_onset" %in% rownames(coefs)) {
      coefs["after_onset", "Estimate"]
    } else {
      NA_real_
    }

    ## ---- interpretable slopes ----
    slope_before <- -beta_before_raw
    slope_after <- beta_after_raw

    se_before_raw <- if ("before_onset" %in% rownames(coefs)) {
      coefs["before_onset", "Std. Error"]
    } else {
      NA_real_
    }

    se_after_raw <- if ("after_onset" %in% rownames(coefs)) {
      coefs["after_onset", "Std. Error"]
    } else {
      NA_real_
    }

    intercept <- if ("(Intercept)" %in% rownames(coefs)) {
      coefs["(Intercept)", "Estimate"]
    } else {
      NA_real_
    }

    ## ---- variance-covariance components ----
    var_before_raw <- if ("before_onset" %in% rownames(vc)) {
      vc["before_onset", "before_onset"]
    } else {
      NA_real_
    }

    var_after_raw <- if ("after_onset" %in% rownames(vc)) {
      vc["after_onset", "after_onset"]
    } else {
      NA_real_
    }

    cov_before_after <- if (
      "before_onset" %in% rownames(vc) &&
      "after_onset" %in% colnames(vc)
    ) {
      vc["before_onset", "after_onset"]
    } else {
      NA_real_
    }

    ## ---- slope-change contrast ----
    delta <- beta_after_raw + beta_before_raw

    var_delta <- var_after_raw + var_before_raw + 2 * cov_before_after

    if (is.na(var_delta) || var_delta < 0) {
      se_delta <- NA_real_
      ci_low <- NA_real_
      ci_high <- NA_real_
    } else {
      se_delta <- sqrt(var_delta)
      ci_low <- delta - 1.96 * se_delta
      ci_high <- delta + 1.96 * se_delta
    }

    ## ---- model-fit metrics ----
    model_df <- stats::model.frame(fit)

    y <- stats::model.response(model_df)

    yhat <- tryCatch(
      stats::predict(fit, newdata = model_df, re.form = NA, allow.new.levels = TRUE),
      error = function(e) rep(NA_real_, length(y))
    )

    mse <- mean((y - yhat)^2, na.rm = TRUE)

    out[[i]] <- tibble::tibble(
      Protein = protein,

      `Beta 1` = slope_before,
      `SE Beta 1` = se_before_raw,
      `Beta 2` = slope_after,
      `SE Beta 2` = se_after_raw,

      `Beta Before Raw` = beta_before_raw,
      `Beta After Raw` = beta_after_raw,
      `Var Before Raw` = var_before_raw,
      `Var After Raw` = var_after_raw,
      `Cov Before After` = cov_before_after,

      Delta = delta,
      `SE Delta` = se_delta,
      `Delta 95% CI Lower` = ci_low,
      `Delta 95% CI Upper` = ci_high,

      Intercept = intercept,
      AIC = stats::AIC(fit),
      BIC = stats::BIC(fit),
      MSE = mse,
      `N Obs` = stats::nobs(fit)
    )
  }

  dplyr::bind_rows(out)
}
