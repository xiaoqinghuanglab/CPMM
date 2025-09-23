#' Fit change-point MMs across multiple proteins and strata
#'
#' Fits three random-intercept MMs per protein:
#'  - Status-change:  protein ~ before_onset + after_onset + covariates + (1 | subject_id)
#'  - Normal-only:    protein ~ before_onset + covariates + (1 | subject_id)
#'  - Abnormal-only:  protein ~ after_onset  + covariates + (1 | subject_id)
#'
#' The change point is fixed at 0 on `years_since_onset_col`. We define
#' `before_onset = pmax(0, -years_since_onset)` and
#' `after_onset  = pmax(0,  years_since_onset)`.
#'
#' @param df_status_change Data frame for status-change subjects.
#' @param df_normal Data frame for always-normal subjects.
#' @param df_abnormal Data frame for always-abnormal subjects.
#' @param protein_list Character vector of protein column names to model.
#' @param covariates Character vector of additional fixed-effect covariates.
#'   Defaults to c("SEX","BASELINE_AGE").
#' @param subject_id_col Column name for subject ID (random intercept). Default "SUBID".
#' @param years_since_onset_col Column name for years since onset. Default "years_since_onset".
#'
#' @return A tibble with one row per protein containing Beta 1–4, their SEs,
#'   the status-model intercept, and AIC/BIC/MSE per stratum.
#' @export
#'
#' @examples
#' # See tests for a minimal synthetic example.
fit_cpmm_all_proteins <- function(
    df_status_change,
    df_normal,
    df_abnormal,
    protein_list,
    covariates = c("SEX", "BASELINE_AGE"),
    subject_id_col = "SUBID",
    years_since_onset_col = "years_since_onset"
) {
  # deps
  requireNamespace("lme4", quietly = TRUE)
  requireNamespace("lmerTest", quietly = TRUE)
  requireNamespace("tibble", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)

  # helper to add piecewise terms
  add_piecewise_terms <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(df)
    yso <- df[[years_since_onset_col]]
    df$before_onset <- pmax(0, -yso)
    df$after_onset  <- pmax(0,  yso)
    df
  }

  # helper: fit a single LMM and extract slopes/SE
  fit_lmm_and_get_slopes <- function(df, rhs_terms) {
    if (is.null(df) || nrow(df) == 0) {
      return(list(
        slope_before = NA_real_, slope_after = NA_real_,
        se_before = NA_real_, se_after = NA_real_,
        fit = NULL
      ))
    }

    # build formula: <protein> ~ <rhs> + (1|id) later per-protein
    list(slope_before = NA_real_, slope_after = NA_real_,
         se_before = NA_real_, se_after = NA_real_, fit = NULL,
         rhs_terms = rhs_terms)
  }

  # prepare data with piecewise terms
  dfs <- list(
    status = add_piecewise_terms(df_status_change),
    normal = add_piecewise_terms(df_normal),
    abnormal = add_piecewise_terms(df_abnormal)
  )

  out_rows <- vector("list", length(protein_list))

  for (i in seq_along(protein_list)) {
    protein <- protein_list[[i]]

    # Build RHS per group
    rhs_status  <- paste(c("before_onset", "after_onset", covariates), collapse = " + ")
    rhs_normal  <- paste(c("before_onset", covariates), collapse = " + ")
    rhs_abnorm  <- paste(c("after_onset",  covariates), collapse = " + ")

    # Fitting sub-function (per protein & group)
    .fit_group <- function(df, rhs) {
      if (is.null(df) || nrow(df) == 0 || !(protein %in% names(df))) {
        return(list(
          slope_before = NA_real_, slope_after = NA_real_,
          se_before = NA_real_, se_after = NA_real_,
          fit = NULL
        ))
      }
      # assemble formula with random intercept
      fml <- as.formula(
        paste0(protein, " ~ ", rhs, " + (1|", subject_id_col, ")")
      )
      # robust fit with tryCatch
      fit <- tryCatch(
        lmerTest::lmer(formula = fml, data = df, REML = TRUE),
        error = function(e) NULL,
        warning = function(w) suppressWarnings(
          lmerTest::lmer(formula = fml, data = df, REML = TRUE)
        )
      )
      if (is.null(fit)) {
        return(list(
          slope_before = NA_real_, slope_after = NA_real_,
          se_before = NA_real_, se_after = NA_real_,
          fit = NULL
        ))
      }
      coefs <- summary(fit)$coefficients
      # pull estimates and SEs if present
      est_before <- if ("before_onset" %in% rownames(coefs)) coefs["before_onset","Estimate"] else NA_real_
      se_before  <- if ("before_onset" %in% rownames(coefs)) coefs["before_onset","Std. Error"] else NA_real_
      est_after  <- if ("after_onset"  %in% rownames(coefs)) coefs["after_onset","Estimate"]  else NA_real_
      se_after   <- if ("after_onset"  %in% rownames(coefs)) coefs["after_onset","Std. Error"] else NA_real_
      list(
        slope_before = if (is.na(est_before)) NA_real_ else -est_before, # match Python sign
        slope_after  = est_after,
        se_before = se_before,
        se_after  = se_after,
        fit = fit
      )
    }

    res_status <- .fit_group(dfs$status,  rhs_status)
    res_normal <- .fit_group(dfs$normal,  rhs_normal)
    res_abnorm <- .fit_group(dfs$abnormal, rhs_abnorm)

    # metrics helper
    .metrics <- function(fit, df) {
      if (is.null(fit) || is.null(df)) return(c(AIC = NA_real_, BIC = NA_real_, MSE = NA_real_))
      y  <- df[[protein]]
      yhat <- tryCatch(
        predict(fit, newdata = df, allow.new.levels = TRUE, re.form = NULL),
        error = function(e) rep(NA_real_, length(y))
      )
      mse <- if (all(is.na(yhat))) NA_real_ else mean((y - yhat)^2, na.rm = TRUE)
      c(AIC = stats::AIC(fit), BIC = stats::BIC(fit), MSE = mse)
    }

    met_status <- .metrics(res_status$fit, dfs$status)
    met_normal <- .metrics(res_normal$fit, dfs$normal)
    met_abnorm <- .metrics(res_abnorm$fit, dfs$abnormal)

    # intercept from status model (if available)
    intercept <- NA_real_
    if (!is.null(res_status$fit)) {
      sc <- summary(res_status$fit)$coefficients
      if ("(Intercept)" %in% rownames(sc)) intercept <- sc["(Intercept)","Estimate"]
    }

    out_rows[[i]] <- tibble::tibble(
      Protein = protein,
      `Beta 1` = res_status$slope_before, `SE Beta 1` = res_status$se_before,
      `Beta 2` = res_normal$slope_before, `SE Beta 2` = res_normal$se_before,
      `Beta 3` = res_status$slope_after,  `SE Beta 3` = res_status$se_after,
      `Beta 4` = res_abnorm$slope_after,  `SE Beta 4` = res_abnorm$se_after,
      Intercept = intercept,
      `AIC Status` = met_status["AIC"], `BIC Status` = met_status["BIC"], `MSE Status` = met_status["MSE"],
      `AIC Normal` = met_normal["AIC"], `BIC Normal` = met_normal["BIC"], `MSE Normal` = met_normal["MSE"],
      `AIC Abnormal` = met_abnorm["AIC"], `BIC Abnormal` = met_abnorm["BIC"], `MSE Abnormal` = met_abnorm["MSE"]
    )
  }

  dplyr::bind_rows(out_rows)
}
