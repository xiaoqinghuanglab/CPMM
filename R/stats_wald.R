#' Wald tests for pre- vs post-onset slope differences
#'
#' Compares (Beta 1 vs Beta 3) and (Beta 2 vs Beta 4) using 1-df Wald tests:
#' \deqn{W = (b_i - b_j)^2 / (se_i^2 + se_j^2),\quad p = 1 - F_{\chi^2_1}(W).}
#' Optionally applies BH-FDR to each set of p-values, ranks results by the chosen
#' (adjusted) p-value, and flags significance at \code{alpha}.
#'
#' @param results_df A data.frame/tibble from \code{fit_cpmm_all_proteins()}
#'   containing columns: \code{Protein}, \code{Beta 1}, \code{SE Beta 1},
#'   \code{Beta 2}, \code{SE Beta 2}, \code{Beta 3}, \code{SE Beta 3},
#'   \code{Beta 4}, \code{SE Beta 4}.
#' @param adjust_p Logical; apply BH-FDR? Default \code{TRUE}.
#' @param rank_by Integer {1,2}; rank by test 1 (status-change slopes) or test 2
#'   (normal vs abnormal slopes). Default \code{1}.
#' @param alpha Numeric; significance threshold. Default \code{0.05}.
#'
#' @return A tibble with columns:
#'   \itemize{
#'     \item \code{Protein}
#'     \item \code{Beta 1}, \code{SE Beta 1}, \code{Beta 3}, \code{SE Beta 3},
#'           \code{Wald Statistic 1}, \code{P-value 1}, and (if adjusted)
#'           \code{Adjusted P-value 1}, \code{Significant 1}
#'     \item \code{Beta 2}, \code{SE Beta 2}, \code{Beta 4}, \code{SE Beta 4},
#'           \code{Wald Statistic 2}, \code{P-value 2}, and (if adjusted)
#'           \code{Adjusted P-value 2}, \code{Significant 2}
#'     \item \code{Rank} (1 = most significant according to \code{rank_by})
#'   }
#' @export
compute_wald_test <- function(
    results_df,
    adjust_p = TRUE,
    rank_by = 1,
    alpha = 0.05
) {
  # deps
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("tibble", quietly = TRUE)

  # expected columns
  need <- c("Protein",
            "Beta 1","SE Beta 1","Beta 3","SE Beta 3",
            "Beta 2","SE Beta 2","Beta 4","SE Beta 4")
  miss <- setdiff(need, names(results_df))
  if (length(miss)) {
    stop("compute_wald_test(): missing columns in results_df: ",
         paste(miss, collapse = ", "), call. = FALSE)
  }

  df <- tibble::as_tibble(results_df)

  # filter out rows with any NA in required columns (matches Python "continue")
  req_num_cols <- need[-1]                   # drop "Protein"
  keep <- !apply(df[, req_num_cols, drop = FALSE], 1, function(r) any(is.na(r)))
  df_use <- df[keep, , drop = FALSE]

  # If nothing to compute, return an empty tibble with proper columns
  if (nrow(df_use) == 0) {
    out <- tibble::tibble(
      Protein = character(),
      `Beta 1` = numeric(), `SE Beta 1` = numeric(),
      `Beta 3` = numeric(), `SE Beta 3` = numeric(),
      `Wald Statistic 1` = numeric(), `P-value 1` = numeric(),
      `Beta 2` = numeric(), `SE Beta 2` = numeric(),
      `Beta 4` = numeric(), `SE Beta 4` = numeric(),
      `Wald Statistic 2` = numeric(), `P-value 2` = numeric()
    )
    if (adjust_p) {
      out$`Adjusted P-value 1` <- numeric()
      out$`Significant 1` <- logical()
      out$`Adjusted P-value 2` <- numeric()
      out$`Significant 2` <- logical()
    }
    out$Rank <- integer()
    return(out)
  }

  # compute Wald stats & p-values
  w1 <- ( (df_use$`Beta 1` - df_use$`Beta 3`)^2 ) /
    ( df_use$`SE Beta 1`^2 + df_use$`SE Beta 3`^2 )
  w2 <- ( (df_use$`Beta 2` - df_use$`Beta 4`)^2 ) /
    ( df_use$`SE Beta 2`^2 + df_use$`SE Beta 4`^2 )

  p1 <- stats::pchisq(w1, df = 1, lower.tail = FALSE)
  p2 <- stats::pchisq(w2, df = 1, lower.tail = FALSE)

  out <- df_use |>
    dplyr::transmute(
      Protein,
      `Beta 1`, `SE Beta 1`,
      `Beta 3`, `SE Beta 3`,
      `Wald Statistic 1` = w1, `P-value 1` = p1,
      `Beta 2`, `SE Beta 2`,
      `Beta 4`, `SE Beta 4`,
      `Wald Statistic 2` = w2, `P-value 2` = p2
    )

  # BH-FDR adjustment (per-test family, same as your loop over i in Python)
  if (isTRUE(adjust_p)) {
    adj1 <- stats::p.adjust(out$`P-value 1`, method = "BH")
    adj2 <- stats::p.adjust(out$`P-value 2`, method = "BH")
    out$`Adjusted P-value 1` <- adj1
    out$`Significant 1` <- adj1 < alpha
    out$`Adjusted P-value 2` <- adj2
    out$`Significant 2` <- adj2 < alpha
  }

  # ranking
  if (rank_by %in% c(1L, 2L)) {
    rank_col <- if (isTRUE(adjust_p)) paste0("Adjusted P-value ", rank_by) else paste0("P-value ", rank_by)
  } else {
    stop("rank_by must be 1 or 2.", call. = FALSE)
  }

  out <- out |>
    dplyr::arrange(.data[[rank_col]]) |>
    dplyr::mutate(Rank = dplyr::row_number())

  out
}
