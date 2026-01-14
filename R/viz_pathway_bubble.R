#' Pathway enrichment bubble plot (compact by category)
#'
#' Creates a compact bubble plot of pathway enrichment colored by -log10(FDR)
#' and sized by the unique gene count per pathway/source.
#'
#' @param df Data frame with pathway results (already cleaned/categorized).
#' @param pathway_col Column with cleaned pathway names. Default `"Cleaned_Pathway"`.
#' @param category_col Column with biological categories. Default `"BioCategory_Manual"`.
#' @param source_col Column indicating enrichment source (e.g., DAVID/Reactome). Default `"Source"`.
#' @param logq_col Column with -log10(FDR) or similar. Default `"LogQValue"`.
#' @param gene_col Column with gene identifiers (for unique counts). Default `"Gene"`.
#' @param size_scale Max bubble size in mm (area-scaled). Default `15`.
#' @param cmap Brewer palette name (e.g., `"RdBu"` or `"RdBu_r"`). Default `"RdBu_r"`.
#' @param title Plot title. Default `"Pathway Enrichment by Source"`.
#' @param style_config Optional list of `ggplot2::theme()` overrides.
#' @param export Save the plot? Default `FALSE`.
#' @param export_dir Output directory. Default `"Figures"`.
#' @param export_name Base filename (no extension). Default `"pathway_bubble_plot"`.
#' @param export_formats Vector of formats, e.g. `c("pdf","svg")`.
#'
#' @return A `ggplot` object.
#' @export
plot_pathway_bubble <- function(
    df,
    pathway_col  = "Cleaned_Pathway",
    category_col = "BioCategory_Manual",
    source_col   = "Source",
    logq_col     = "LogQValue",
    gene_col     = "Gene",
    size_scale   = 15,
    cmap         = "RdBu_r",
    title        = "Pathway Enrichment by Source",
    style_config = NULL,
    export       = FALSE,
    export_dir   = "Figures",
    export_name  = "pathway_bubble_plot",
    export_formats = c("pdf","svg")
) {
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("tibble", quietly = TRUE)
  requireNamespace("rlang", quietly = TRUE)

  need <- c(pathway_col, category_col, source_col, logq_col, gene_col)
  miss <- setdiff(need, names(df))
  if (length(miss)) {
    stop("plot_pathway_bubble(): missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }

  # -------- Filter valid pathways --------
  df_valid <- df[!is.na(df[[pathway_col]]), , drop = FALSE]

  # -------- Compute compact y positions grouped by category --------
  pathway_order_df <- df_valid |>
    dplyr::select(dplyr::all_of(c(pathway_col, category_col))) |>
    dplyr::distinct() |>
    dplyr::arrange(.data[[category_col]], .data[[pathway_col]])

  # counts per category and running offsets with a gap of 1
  cat_counts <- pathway_order_df |>
    dplyr::count(.data[[category_col]], name = "n") |>
    dplyr::arrange(.data[[category_col]]) |>
    dplyr::mutate(offset = dplyr::lag(cumsum(n + 1L), default = 0L))

  pathway_order_df <- pathway_order_df |>
    dplyr::left_join(cat_counts, by = rlang::set_names(category_col, category_col)) |>
    dplyr::group_by(.data[[category_col]]) |>
    dplyr::mutate(compact_y_pos = (dplyr::row_number() - 1L) + .data$offset) |>
    dplyr::ungroup()

  # map y pos to main df
  ypos_map <- stats::setNames(pathway_order_df$compact_y_pos, pathway_order_df[[pathway_col]])
  df_valid$compact_y_pos <- unname(ypos_map[df_valid[[pathway_col]]])

  # -------- Aggregate for plotting: max(-log10 FDR), unique gene count --------
  plot_df <- df_valid |>
    dplyr::group_by(.data[[pathway_col]], .data[[source_col]], .data$compact_y_pos) |>
    dplyr::summarise(
      LogQValue = max(.data[[logq_col]], na.rm = TRUE),
      GeneCount = dplyr::n_distinct(.data[[gene_col]]),
      .groups = "drop"
    )

  if (!nrow(plot_df)) stop("plot_pathway_bubble(): nothing to plot after aggregation.", call. = FALSE)

  # x positions / order for sources
  sources <- unique(plot_df[[source_col]])
  plot_df[[source_col]] <- factor(plot_df[[source_col]], levels = sources)

  # category separators (end y per category)
  cat_ends <- pathway_order_df |>
    dplyr::group_by(.data[[category_col]]) |>
    dplyr::summarise(end_y = max(.data$compact_y_pos), .groups = "drop")

  # legend size breaks = quartiles of unique GeneCount
  uniq_sizes <- sort(unique(plot_df$GeneCount))
  if (length(uniq_sizes) >= 3) {
    legend_sizes <- unique(round(stats::quantile(uniq_sizes, probs = c(0.25, 0.5, 0.75), names = FALSE)))
  } else {
    legend_sizes <- uniq_sizes
  }

  # palette handling ("RdBu" vs "RdBu_r")
  pal_name <- sub("_r$", "", cmap)
  pal_dir  <- if (grepl("_r$", cmap)) -1 else 1

  # dynamic height roughly: (max_y + 2) * 0.5 inches
  max_y <- max(pathway_order_df$compact_y_pos, na.rm = TRUE)
  export_height <- (max_y + 2) * 0.5
  export_width  <- 7.5

  # -------- Build plot --------
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[[source_col]], y = .data$compact_y_pos)) +
    ggplot2::geom_point(
      ggplot2::aes(size = .data$GeneCount, fill = .data$LogQValue),
      shape = 21, color = "black", stroke = 0.3
    ) +
    # Y axis: pathways (reverse so first at top)
    ggplot2::scale_y_reverse(
      breaks = pathway_order_df$compact_y_pos,
      labels = pathway_order_df[[pathway_col]]
    ) +
    # X axis: sources
    ggplot2::scale_x_discrete() +
    # Colorbar for -log10(FDR)
    ggplot2::scale_fill_distiller(palette = pal_name, direction = pal_dir, name = "-log10(FDR)") +
    # Bubble size (area-scaled) with quartile legend
    ggplot2::scale_size_area(
      max_size = size_scale,
      breaks = legend_sizes,
      labels = paste0(legend_sizes, " genes"),
      name = "Gene Count"
    ) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "grey90"),
      plot.title = ggplot2::element_text(size = 16, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 10)
    )

  # category separators
  for (yy in cat_ends$end_y) {
    p <- p + ggplot2::geom_hline(yintercept = yy + 0.5, linetype = "dashed", linewidth = 0.6, color = "grey50")
  }

  # optional theme overrides
  if (!is.null(style_config) && length(style_config)) {
    p <- p + do.call(ggplot2::theme, style_config)
  }

  # -------- Export --------
  if (isTRUE(export)) {
    if (!dir.exists(export_dir)) dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
    for (fmt in export_formats) {
      ggplot2::ggsave(
        filename = file.path(export_dir, paste0(export_name, ".", fmt)),
        plot = p, width = export_width, height = export_height, dpi = 300
      )
    }
  }

  p
}
