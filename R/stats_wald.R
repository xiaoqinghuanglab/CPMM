#' Wald test for pre- vs post-onset slope difference (CPMM)
#'
#' Performs a Wald test comparing the pre-onset and post-onset slopes
#' from a CPMM:
#'
#' \deqn{W = (\beta_1 - \beta_2)^2 / (SE_1^2 + SE_2^2), \quad
#' p = 1 - F_{\chi^2_1}(W)}
#'
#' Optionally applies Benjamini–Hochberg FDR correction, flags significance,
#' and ranks proteins by (adjusted) p-value.
#'
#' @param results_df A data.frame/tibble from \code{fit_cpmm_all_proteins()}
#'   containing columns: \code{Protein}, \code{Beta 1}, \code{SE Beta 1},
#'   \code{Beta 2}, \code{SE Beta 2}.
#' @param adjust_p Logical; apply BH-FDR? Default \code{TRUE}.
#' @param alpha Numeric; significance threshold. Default \code{0.05}.
#'
#' @return A tibble with columns:
#'   \itemize{
#'     \item \code{Protein}
#'     \item \code{Beta 1}, \code{SE Beta 1}
#'     \item \code{Beta 2}, \code{SE Beta 2}
#'     \item \code{Wald Statistic}
#'     \item \code{P-value}
#'     \item \code{Adjusted P-value} (if \code{adjust_p = TRUE})
#'     \item \code{Significant} (logical; if \code{adjust_p = TRUE})
#'     \item \code{Rank}
#'   }
#' @export
compute_wald_test <- function(
    results_df,
    adjust_p = TRUE,
    alpha = 0.05
) {
  # deps
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("tibble", quietly = TRUE)

  # expected columns
  need <- c("Protein", "Beta 1", "SE Beta 1", "Beta 2", "SE Beta 2")
  miss <- setdiff(need, names(results_df))
  if (length(miss)) {
    stop(
      "compute_wald_test(): missing columns in results_df: ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }

  df <- tibble::as_tibble(results_df)

  # drop rows with missing betas or SEs (Python: continue)
  keep <- !apply(df[, need[-1], drop = FALSE], 1, function(x) any(is.na(x)))
  df_use <- df[keep, , drop = FALSE]

  if (nrow(df_use) == 0) {
    out <- tibble::tibble(
      Protein = character(),
      `Beta 1` = numeric(), `SE Beta 1` = numeric(),
      `Beta 2` = numeric(), `SE Beta 2` = numeric(),
      `Wald Statistic` = numeric(),
      `P-value` = numeric(),
      Rank = integer()
    )
    if (adjust_p) {
      out$`Adjusted P-value` <- numeric()
      out$Significant <- logical()
    }
    return(out)
  }

  # Wald statistic
  wald_stat <- (df_use$`Beta 1` - df_use$`Beta 2`)^2 /
    (df_use$`SE Beta 1`^2 + df_use$`SE Beta 2`^2)

  pval <- stats::pchisq(wald_stat, df = 1, lower.tail = FALSE)

  out <- df_use |>
    dplyr::transmute(
      Protein,
      `Beta 1`, `SE Beta 1`,
      `Beta 2`, `SE Beta 2`,
      `Wald Statistic` = wald_stat,
      `P-value` = pval
    )

  if (isTRUE(adjust_p)) {
    out$`Adjusted P-value` <- stats::p.adjust(out$`P-value`, method = "BH")
    out$Significant <- out$`Adjusted P-value` < alpha
    rank_col <- "Adjusted P-value"
  } else {
    rank_col <- "P-value"
  }

  out |>
    dplyr::arrange(.data[[rank_col]]) |>
    dplyr::mutate(Rank = dplyr::row_number())
}
