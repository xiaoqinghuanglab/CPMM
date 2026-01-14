#' Quadrant plot of CPMM beta coefficients
#'
#' Plots pre-onset slope (Beta 1) vs post-onset slope (Beta 2) with axes at 0.
#' Points are highlighted by FDR significance and quadrant counts are displayed.
#'
#' @param wald_df Data frame from compute_wald_test().
#' @param beta_x_col Column for pre-onset slope. Default "Beta 1".
#' @param beta_y_col Column for post-onset slope. Default "Beta 2".
#' @param fdr_col Column for adjusted p-values. Default "Adjusted P-value".
#' @param protein_col Column with protein identifiers. Default "Protein".
#' @param fdr_threshold Numeric FDR cutoff. Default 0.05.
#' @param annotate Logical; annotate significant proteins? Default TRUE.
#' @param style_config Optional ggplot2 theme overrides.
#' @param export Logical; save plot? Default FALSE.
#' @param export_dir Output directory. Default "Figures".
#' @param export_name Base filename. Default "quadrant_plot".
#' @param export_formats File formats. Default c("pdf","svg").
#'
#' @return A ggplot object.
#' @export
plot_quadrant_beta <- function(
    wald_df,
    beta_x_col = "Beta 1",
    beta_y_col = "Beta 2",
    fdr_col = "Adjusted P-value",
    protein_col = "Protein",
    fdr_threshold = 0.05,
    annotate = TRUE,
    style_config = NULL,
    export = FALSE,
    export_dir = "Figures",
    export_name = "quadrant_plot",
    export_formats = c("pdf", "svg")
) {
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("tibble", quietly = TRUE)

  df <- tibble::as_tibble(wald_df)

  # Basic column checks
  need <- c(beta_x_col, beta_y_col, fdr_col, protein_col)
  miss <- setdiff(need, names(df))
  if (length(miss)) {
    stop("plot_quadrant_beta(): missing columns: ",
         paste(miss, collapse = ", "), call. = FALSE)
  }

  # Numeric vectors
  x <- suppressWarnings(as.numeric(df[[beta_x_col]]))
  y <- suppressWarnings(as.numeric(df[[beta_y_col]]))
  fdr <- suppressWarnings(as.numeric(df[[fdr_col]]))

  keep <- is.finite(x) & is.finite(y) & is.finite(fdr)
  dfp <- tibble::tibble(
    x = x[keep],
    y = y[keep],
    Protein = df[[protein_col]][keep],
    Significant = fdr[keep] < fdr_threshold
  )

  if (!nrow(dfp)) {
    stop("plot_quadrant_beta(): no finite values to plot.", call. = FALSE)
  }

  # Quadrant counts
  q1 <- sum(dfp$x > 0 & dfp$y > 0)
  q2 <- sum(dfp$x < 0 & dfp$y > 0)
  q3 <- sum(dfp$x < 0 & dfp$y < 0)
  q4 <- sum(dfp$x > 0 & dfp$y < 0)

  # Axis limits (10% padding)
  xlim <- range(dfp$x) * 1.1
  ylim <- range(dfp$y) * 1.1

  p <- ggplot2::ggplot(dfp, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(
      ggplot2::aes(color = Significant),
      alpha = 0.8,
      size = 2
    ) +
    ggplot2::scale_color_manual(
      values = c("FALSE" = "grey50", "TRUE" = "#B24745")
    ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.5) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.5) +
    ggplot2::coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE) +
    ggplot2::labs(
      x = "Pre-onset slope (Beta 1)",
      y = "Post-onset slope (Beta 2)",
      title = "Quadrant Plot of CPMM Slopes"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      legend.position = "none",
      panel.border = ggplot2::element_rect(color = "black", fill = NA)
    )

  # Annotate significant proteins
  if (isTRUE(annotate) && any(dfp$Significant)) {
    p <- p + ggplot2::geom_text(
      data = dfp[dfp$Significant, ],
      ggplot2::aes(label = Protein),
      size = 2.8,
      hjust = 1,
      vjust = 0
    )
  }

  # Quadrant labels
  p <- p +
    ggplot2::annotate("text", x = xlim[2] * 0.7, y = ylim[2] * 0.9, label = paste0("n=", q1)) +
    ggplot2::annotate("text", x = xlim[1] * 0.7, y = ylim[2] * 0.9, label = paste0("n=", q2)) +
    ggplot2::annotate("text", x = xlim[1] * 0.7, y = ylim[1] * 0.9, label = paste0("n=", q3)) +
    ggplot2::annotate("text", x = xlim[2] * 0.7, y = ylim[1] * 0.9, label = paste0("n=", q4))

  if (!is.null(style_config)) {
    p <- p + do.call(ggplot2::theme, style_config)
  }

  if (isTRUE(export)) {
    if (!dir.exists(export_dir)) dir.create(export_dir, recursive = TRUE)
    for (fmt in export_formats) {
      ggplot2::ggsave(
        file.path(export_dir, paste0(export_name, ".", fmt)),
        p, width = 6, height = 6, dpi = 300
      )
    }
  }

  p
}
