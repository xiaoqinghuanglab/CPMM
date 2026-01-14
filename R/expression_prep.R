#' Prepare long-format expression data across multiple strata
#'
#' Melts selected genes from one or more data frames into a single
#' long-format expression table. Each input can optionally be filtered
#' by diagnosis/category and assigned a custom source label.
#'
#' @param inputs A list of lists. Each element must contain:
#'   - df: data frame
#'   - label: source label to assign
#'   - category_col (optional): column name for filtering
#'   - categories (optional): vector of category values to keep
#' @param subset_genes Character vector of gene/protein columns to include.
#'
#' @return A tibble with columns Gene, Expression, Source.
#'   Gene is an ordered factor matching subset_genes.
#' @export
prepare_combined_expression <- function(
    inputs,
    subset_genes
) {
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("tidyr", quietly = TRUE)
  requireNamespace("tibble", quietly = TRUE)

  melt_one <- function(df, label) {
    if (is.null(df)) return(NULL)
    if (!all(subset_genes %in% names(df))) return(NULL)

    df |>
      dplyr::select(dplyr::all_of(subset_genes)) |>
      tidyr::pivot_longer(
        cols = dplyr::everything(),
        names_to = "Gene",
        values_to = "Expression"
      ) |>
      dplyr::mutate(
        Expression = suppressWarnings(as.numeric(.data$Expression)),
        Source = label
      )
  }

  out <- list()

  for (i in seq_along(inputs)) {
    x <- inputs[[i]]
    df <- x$df
    label <- x$label

    if (is.null(df) || is.null(label)) next

    # optional category filtering
    if (!is.null(x$category_col) && !is.null(x$categories)) {
      if (x$category_col %in% names(df)) {
        df <- df[df[[x$category_col]] %in% x$categories, , drop = FALSE]
      } else {
        next
      }
    }

    melted <- melt_one(df, label)
    if (!is.null(melted)) out[[length(out) + 1]] <- melted
  }

  combined <- dplyr::bind_rows(out)

  if (!nrow(combined)) {
    return(tibble::tibble(
      Gene = factor(levels = subset_genes),
      Expression = numeric(),
      Source = character()
    ))
  }

  combined |>
    dplyr::mutate(
      Gene = factor(.data$Gene, levels = subset_genes, ordered = TRUE)
    ) |>
    tibble::as_tibble()
}
