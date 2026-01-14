#' Volcano plot for CPMM Wald test results
#'
#' Builds a volcano plot using the signed slope difference
#' \eqn{\Delta \beta = \beta_{after} - \beta_{before}} on the x-axis
#' and \eqn{-\log_{10}(p)} on the y-axis.
#'
#' @param wald_df Data frame from \code{compute_wald_test()}.
#' @param beta_before_col Column for pre-onset slope. Default \code{"Beta 1"}.
#' @param beta_after_col Column for post-onset slope. Default \code{"Beta 2"}.
#' @param pval_col Column with raw p-values. Default \code{"P-value"}.
#' @param fdr_col Column with adjusted p-values. Default \code{"Adjusted P-value"}.
#' @param protein_col Column with protein/gene names. Default \code{"Protein"}.
#' @param pval_threshold Raw p-value cutoff. Default \code{0.05}.
#' @param fdr_threshold FDR cutoff. Default \code{0.05}.
#' @param annotate Logical; annotate significant proteins? Default \code{TRUE}.
#' @param annotate_list Optional vector of protein names to annotate.
#' @param style_config Optional list of ggplot2 theme overrides.
#' @param export Logical; save plot? Default \code{FALSE}.
#' @param export_dir Directory to save files. Default \code{"Figures"}.
#' @param export_name Base filename. Default \code{"wald_volcano"}.
#' @param export_formats File formats to export. Default \code{c("pdf","svg")}.
#'
#' @return A ggplot object.
#' @export
plot_wald_volcano <- function(
    wald_df,
    beta_before_col = "Beta 1",
    beta_after_col  = "Beta 2",
    pval_col = "P-value",
    fdr_col  = "Adjusted P-value",
    protein_col = "Protein",
    pval_threshold = 0.05,
    fdr_threshold  = 0.05,
    annotate = TRUE,
    annotate_list = NULL,
    style_config = NULL,
    export = FALSE,
    export_dir = "Figures",
    export_name = "wald_volcano",
    export_formats = c("pdf","svg")
) {
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("tibble", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)

  df <- tibble::as_tibble(wald_df)

  # Required columns
  need <- c(beta_before_col, beta_after_col, pval_col, protein_col)
  miss <- setdiff(need, names(df))
  if (length(miss)) {
    stop("plot_wald_volcano(): missing columns: ",
         paste(miss, collapse = ", "), call. = FALSE)
  }

  # Derived quantities
  df <- df |>
    dplyr::mutate(
      delta_beta = .data[[beta_after_col]] - .data[[beta_before_col]],
      neglog10p  = -log10(.data[[pval_col]]),
      Significance = dplyr::if_else(
        .data[[pval_col]] < pval_threshold &
          (!fdr_col %in% names(df) | .data[[fdr_col]] < fdr_threshold),
        "Significant", "Not Significant"
      )
    ) |>
    dplyr::filter(is.finite(delta_beta), is.finite(neglog10p))

  # Base plot
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$delta_beta, y = .data$neglog10p, fill = .data$Significance)
  ) +
    ggplot2::geom_hline(
      yintercept = -log10(pval_threshold),
      linetype = "dashed",
      linewidth = 0.5,
      color = "black"
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.5,
      color = "steelblue"
    ) +
    ggplot2::geom_point(
      shape = 21,
      color = "black",
      stroke = 0.3,
      alpha = 0.75,
      size = 2
    ) +
    ggplot2::scale_fill_manual(
      values = c("Not Significant" = "gray70", "Significant" = "#B24745")
    ) +
    ggplot2::labs(
      x = expression(Delta * " slope (Post − Pre)"),
      y = expression(-log[10](p)),
      title = "CPMM Wald Test Volcano Plot",
      fill = "Significance"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(color = "#DDDDDD", linewidth = 0.4),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = 11),
      axis.text  = ggplot2::element_text(size = 10),
      plot.title = ggplot2::element_text(size = 13, face = "bold"),
      legend.position = "right"
    )

  # Annotations
  if (isTRUE(annotate)) {
    lab_df <- if (!is.null(annotate_list)) {
      df[df[[protein_col]] %in% annotate_list, , drop = FALSE]
    } else {
      df[df$Significance == "Significant", , drop = FALSE]
    }

    if (nrow(lab_df)) {
      p <- p + ggplot2::geom_text(
        data = lab_df,
        ggplot2::aes(label = .data[[protein_col]]),
        size = 3,
        vjust = -0.3,
        hjust = 0.5
      )
    }
  }

  # Style overrides
  if (!is.null(style_config)) {
    p <- p + do.call(ggplot2::theme, style_config)
  }

  # Export
  if (isTRUE(export)) {
    if (!dir.exists(export_dir)) dir.create(export_dir, recursive = TRUE)
    for (fmt in export_formats) {
      ggplot2::ggsave(
        filename = file.path(export_dir, paste0(export_name, ".", fmt)),
        plot = p,
        width = 7,
        height = 5,
        dpi = 300
      )
    }
  }

  p
}
