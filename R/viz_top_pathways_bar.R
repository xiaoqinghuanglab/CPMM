#' Top pathways by gene count (bar plot)
#'
#' Ranks pathways by unique gene count, colors by biological category, and
#' annotates each bar with a label derived from the max of \code{logq_col}
#' per pathway (default: "LogQ = {val:.2f}").
#'
#' @param df Data frame with pathway enrichment results.
#' @param pathway_col Column for pathway names. Default `"Cleaned_Pathway"`.
#' @param gene_col Column for gene identifiers. Default `"Gene"`.
#' @param logq_col Column for log-transformed q-values (e.g., -log10 FDR). Default `"LogQValue"`.
#' @param category_col Column for biological categories. Default `"BioCategory_Manual"`.
#' @param top_n Number of pathways to display. Default `15`.
#' @param palette Optional vector of hex colors for categories. If `NULL`, a default
#' palette is used and truncated to the number of categories.
#' @param annotate Annotate bars? Default `TRUE`.
#' @param annotation_format Template for annotations. Supports `{val:.2f}`, `{val:.1e}`, or `{val}`.
#'   Default `"LogQ = {val:.2f}"`.
#' @param title Plot title. Default `"Top Pathways by Gene Count"`.
#' @param export Save figure to disk? Default `FALSE`.
#' @param export_dir Output directory. Default `"Figures"`.
#' @param export_name Base filename (no extension). Default `"top_pathways_bar"`.
#' @param export_formats Vector of formats (e.g., `c("pdf","svg")`).
#'
#' @return A `ggplot` object.
#' @export
plot_top_pathways_bar <- function(
    df,
    pathway_col   = "Cleaned_Pathway",
    gene_col      = "Gene",
    logq_col      = "LogQValue",
    category_col  = "BioCategory_Manual",
    top_n         = 15,
    palette       = NULL,
    annotate      = TRUE,
    annotation_format = "LogQ = {val:.2f}",
    title         = "Top Pathways by Gene Count",
    export        = FALSE,
    export_dir    = "Figures",
    export_name   = "top_pathways_bar",
    export_formats = c("pdf","svg")
) {
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("tibble", quietly = TRUE)
  requireNamespace("rlang", quietly = TRUE)

  # ---- checks ----
  need <- c(pathway_col, gene_col, logq_col, category_col)
  miss <- setdiff(need, names(df))
  if (length(miss)) {
    stop("plot_top_pathways_bar(): missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }

  # ---- summarise & rank ----
  # unique gene count per pathway
  top_pathways <- df |>
    dplyr::group_by(.data[[pathway_col]]) |>
    dplyr::summarise(
      GeneCount = dplyr::n_distinct(.data[[gene_col]]),
      !!logq_col := max(.data[[logq_col]], na.rm = TRUE),
      !!category_col := dplyr::first(.data[[category_col]]),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(.data$GeneCount)) |>
    dplyr::slice_head(n = top_n)

  if (!nrow(top_pathways)) stop("No pathways available to plot.", call. = FALSE)

  # order y by GeneCount descending
  top_pathways[[pathway_col]] <- factor(
    top_pathways[[pathway_col]],
    levels = rev(top_pathways[[pathway_col]][order(top_pathways$GeneCount, decreasing = TRUE)])
  )

  # ---- palette ----
  all_cats <- sort(unique(df[[category_col]]))

  if (is.null(palette)) {
    palette <- c(
      "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
      "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF"
    )
  }

  # recycle safely if needed
  palette <- rep(palette, length.out = length(all_cats))
  names(palette) <- all_cats

  # colors used in the plotted subset
  cats_used <- levels(factor(top_pathways[[category_col]]))
  pal_used <- unname(palette[cats_used])

  # ---- small helper to render annotation_format ----
  fmt_label <- function(val, template) {
    if (grepl("\\{val:", template)) {
      spec <- sub("^.*\\{val:([^}]*)\\}.*$", "\\1", template)
      if (grepl("^\\.\\d+f$", spec)) {
        digits <- as.integer(sub("^\\.(\\d+)f$", "\\1", spec))
        rep <- formatC(val, format = "f", digits = digits)
      } else if (grepl("^\\.\\d+e$", spec)) {
        digits <- as.integer(sub("^\\.(\\d+)e$", "\\1", spec))
        rep <- formatC(val, format = "e", digits = digits)
      } else {
        rep <- as.character(val)
      }
      sub("\\{val:[^}]*\\}", rep, template)
    } else if (grepl("\\{val\\}", template)) {
      sub("\\{val\\}", as.character(val), template)
    } else {
      # fallback: append value
      paste0(template, " ", val)
    }
  }

  ann_labs <- vapply(top_pathways[[logq_col]], fmt_label, character(1), template = annotation_format)

  # ---- build plot ----
  p <- ggplot2::ggplot(
    top_pathways,
    ggplot2::aes(x = .data$GeneCount, y = .data[[pathway_col]], fill = .data[[category_col]])
  ) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::scale_fill_manual(values = pal_used, name = "Biological Category") +
    ggplot2::labs(x = "Gene Count", y = "Pathway", title = title) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 10),
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      legend.position = "right"
    )

  if (isTRUE(annotate)) {
    x_max <- max(top_pathways$GeneCount, na.rm = TRUE)
    p <- p +
      ggplot2::geom_text(
        ggplot2::aes(label = ann_labs),
        hjust = 0, nudge_x = 0.3, size = 3
      ) +
      ggplot2::coord_cartesian(xlim = c(0, x_max + max(1, 0.15 * x_max)), clip = "off")
  }

  # ---- export ----
  if (isTRUE(export)) {
    if (!dir.exists(export_dir)) dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
    for (fmt in export_formats) {
      ggplot2::ggsave(
        filename = file.path(export_dir, paste0(export_name, ".", fmt)),
        plot = p, width = 9, height = 6, dpi = 300
      )
    }
  }

  p
}
