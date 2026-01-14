#' Fit change-point mixed models (CPMM) across proteins
#'
#' Fits a random-intercept CPMM per protein
#'
#'   protein ~ before_onset + after_onset + covariates + (1 | subject_id)
#'

#'
#' @details
#' The change point is fixed at 0 on `years_since_onset_col`.
#' Piecewise terms are defined as:
#'   - before_onset = pmax(0, -years_since_onset)
#'   - after_onset  = pmax(0,  years_since_onset)

#' @param df_status_change Data frame for status-change subjects only.
#' @param protein_list Character vector of protein column names.
#' @param covariates Fixed-effect covariates. Default c("SEX","BASELINE_AGE").
#' @param subject_id_col Subject ID column. Default "SUBID".
#' @param years_since_onset_col Years-since-onset column. Default "years_since_onset".
#'
#' @return A tibble with one row per protein containing Beta 1–2, SEs,
#' intercept, and model fit metrics.
#' @export
fit_cpmm_all_proteins <- function(
    df_status_change,
    protein_list,
    covariates = c("SEX", "BASELINE_AGE"),
    subject_id_col = "SUBID",
    years_since_onset_col = "years_since_onset"
) {

  ## ---- dependencies ----
  pkgs <- c("lme4", "lmerTest", "tibble", "dplyr")
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required.", p), call. = FALSE)
    }
  }

  ## ---- add piecewise terms ----
  df <- df_status_change
  yso <- df[[years_since_onset_col]]
  df$before_onset <- pmax(0, -yso)
  df$after_onset  <- pmax(0,  yso)

  out <- vector("list", length(protein_list))

  for (i in seq_along(protein_list)) {
    protein <- protein_list[[i]]

    if (!(protein %in% names(df))) {
      out[[i]] <- tibble::tibble(
        Protein = protein,
        `Beta 1` = NA_real_, `SE Beta 1` = NA_real_,
        `Beta 2` = NA_real_, `SE Beta 2` = NA_real_,
        Intercept = NA_real_,
        AIC = NA_real_, BIC = NA_real_, MSE = NA_real_
      )
      next
    }

    rhs <- paste(c("before_onset", "after_onset", covariates), collapse = " + ")

    fml <- as.formula(
      paste0(protein, " ~ ", rhs, " + (1|", subject_id_col, ")")
    )

    fit <- tryCatch(
      lmerTest::lmer(fml, data = df, REML = TRUE),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      out[[i]] <- tibble::tibble(
        Protein = protein,
        `Beta 1` = NA_real_, `SE Beta 1` = NA_real_,
        `Beta 2` = NA_real_, `SE Beta 2` = NA_real_,
        Intercept = NA_real_,
        AIC = NA_real_, BIC = NA_real_, MSE = NA_real_
      )
      next
    }

    coefs <- summary(fit)$coefficients

    beta_before <- if ("before_onset" %in% rownames(coefs))
      -coefs["before_onset", "Estimate"] else NA_real_

    se_before <- if ("before_onset" %in% rownames(coefs))
      coefs["before_onset", "Std. Error"] else NA_real_

    beta_after <- if ("after_onset" %in% rownames(coefs))
      coefs["after_onset", "Estimate"] else NA_real_

    se_after <- if ("after_onset" %in% rownames(coefs))
      coefs["after_onset", "Std. Error"] else NA_real_

    intercept <- if ("(Intercept)" %in% rownames(coefs))
      coefs["(Intercept)", "Estimate"] else NA_real_

    y <- df[[protein]]
    yhat <- predict(fit, newdata = df, allow.new.levels = TRUE)
    mse <- mean((y - yhat)^2, na.rm = TRUE)

    out[[i]] <- tibble::tibble(
      Protein = protein,
      `Beta 1` = beta_before, `SE Beta 1` = se_before,
      `Beta 2` = beta_after,  `SE Beta 2` = se_after,
      Intercept = intercept,
      AIC = AIC(fit),
      BIC = BIC(fit),
      MSE = mse
    )
  }

  dplyr::bind_rows(out)
}
