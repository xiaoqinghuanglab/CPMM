#' Boxplots of expression across diagnostic categories for selected genes
#'
#' @param combined_expr Long-format data frame with columns Gene, Expression, Source (defaults).
#' @param gene_order Character vector giving the x-axis order of genes.
#' @param hue_col Column for grouping (e.g., diagnostic categories). Default "Source".
#' @param expression_col Column with expression values. Default "Expression".
#' @param gene_col Column with gene names. Default "Gene".
#' @param palette Optional named character vector mapping categories -> colors.
#'   If NULL, a CPMM default diagnostic palette is used.
#' @param title Plot title.
#' @param y_limit Optional numeric; upper y-axis limit.
#' @param figsize Numeric length-2 vector of inches \code{c(width, height)} for export.
#' @param rotation Integer angle for x tick labels. Default 45.
#' @param export Logical; save plot to disk? Default FALSE.
#' @param export_dir Directory to save into if \code{export=TRUE}. Default "Figures".
#' @param export_name Base filename (no extension). Default "gene_expression_boxplot".
#' @param export_formats Character vector of extensions (e.g., c("pdf","svg")).
#'
#' @return A ggplot object.
#' @export
plot_expression_boxplot <- function(
    combined_expr,
    gene_order,
    hue_col        = "Source",
    expression_col = "Expression",
    gene_col       = "Gene",
    palette        = NULL,
    title          = "Expression of Top Genes Across Categories",
    y_limit        = NULL,
    figsize        = c(32, 8),
    rotation       = 45,
    export         = FALSE,
    export_dir     = "Figures",
    export_name    = "gene_expression_boxplot",
    export_formats = c("pdf", "svg")
) {
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("tibble", quietly = TRUE)
  requireNamespace("rlang", quietly = TRUE)

  # ---- checks ----
  stopifnot(is.data.frame(combined_expr))
  need <- c(hue_col, expression_col, gene_col)
  miss <- setdiff(need, names(combined_expr))
  if (length(miss)) {
    stop("plot_expression_boxplot(): missing columns: ",
         paste(miss, collapse = ", "), call. = FALSE)
  }
  if (length(figsize) != 2) {
    stop("figsize must be length-2 numeric: c(width, height).", call. = FALSE)
  }

  df <- tibble::as_tibble(combined_expr)

  # ensure gene order
  df[[gene_col]] <- factor(df[[gene_col]], levels = gene_order, ordered = TRUE)

  # ---- default CPMM palette ----
  if (is.null(palette)) {
    palette <- c(
      Normal_only  = "#00A1D5",
      SCD          = "#DF8F44",
      MCI          = "#B24745",
      AD_Dementia  = "#79AF97",
      FTD_Dementia = "#6A6599"
    )
  }

  # ensure palette covers all groups
  groups <- sort(unique(df[[hue_col]]))
  if (is.null(names(palette))) {
    names(palette) <- groups
  }
  if (!all(groups %in% names(palette))) {
    pal_tmp <- rep(palette, length.out = length(groups))
    names(pal_tmp) <- groups
    palette <- pal_tmp
  }
  pal_used <- unname(palette[groups])

  # ---- plot ----
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data[[gene_col]],
      y = .data[[expression_col]],
      fill = .data[[hue_col]]
    )
  ) +
    ggplot2::geom_boxplot(
      width = 0.9,
      linewidth = 0.8,
      position = ggplot2::position_dodge2(preserve = "single"),
      outlier.shape = 21,
      outlier.fill  = "white",
      outlier.color = "black",
      outlier.size  = 1.8,
      outlier.stroke = 0.5
    ) +
    ggplot2::scale_fill_manual(values = pal_used, name = "Category") +
    ggplot2::labs(
      x = "Gene",
      y = "Expression",
      title = title
    ) +
    ggplot2::theme_minimal(base_size = 20) +
    ggplot2::theme(
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(linewidth = 0.4),
      panel.border       = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1.5),
      axis.title         = ggplot2::element_text(size = 20),
      axis.text          = ggplot2::element_text(size = 18),
      axis.text.x        = ggplot2::element_text(angle = rotation, hjust = 1),
      legend.title       = ggplot2::element_text(size = 18),
      legend.text        = ggplot2::element_text(size = 18),
      plot.title         = ggplot2::element_text(size = 22, face = "bold"),
      legend.position    = "right"
    )

  if (!is.null(y_limit) && is.finite(y_limit)) {
    p <- p + ggplot2::coord_cartesian(ylim = c(NA, y_limit))
  }

  # ---- export ----
  if (isTRUE(export)) {
    if (!dir.exists(export_dir)) {
      dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
    }
    for (fmt in export_formats) {
      ggplot2::ggsave(
        filename = file.path(export_dir, paste0(export_name, ".", fmt)),
        plot = p,
        width = figsize[1],
        height = figsize[2],
        dpi = 300
      )
    }
  }

  p
}
