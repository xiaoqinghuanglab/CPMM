#' Pathway–Gene membership heatmap (grouped by category)
#'
#' Creates a compact heatmap showing which genes belong to which pathways.
#' Y-axis pathways are grouped by biological category with dashed separators.
#' Each tile is filled by the pathway's category color.
#'
#' @param df Data frame containing pathways, categories, and genes.
#' @param pathway_col Column name for pathways. Default `"Cleaned_Pathway"`.
#' @param category_col Column name for pathway categories. Default `"BioCategory_Manual"`.
#' @param gene_col Column name for genes. Default `"Gene"`.
#' @param palette Vector of hex colors for categories (names optional). If `NULL`,
#'   a default palette is used and truncated to the number of categories.
#' @param style_config Optional list of `ggplot2::theme()` overrides.
#' @param title Plot title. Default `"Pathway-Gene Membership Heatmap"`.
#' @param export Save the plot? Default `FALSE`.
#' @param export_dir Directory to save into. Default `"Figures"`.
#' @param export_name Base filename (no extension). Default `"pathway_gene_heatmap"`.
#' @param export_formats Vector of extensions, e.g. `c("pdf","svg")`.
#'
#' @return A `ggplot` object.
#' @export
plot_pathway_gene_heatmap <- function(
    df,
    pathway_col  = "Cleaned_Pathway",
    category_col = "BioCategory_Manual",
    gene_col     = "Gene",
    palette      = NULL,
    style_config = NULL,
    title        = "Pathway-Gene Membership Heatmap",
    export       = FALSE,
    export_dir   = "Figures",
    export_name  = "pathway_gene_heatmap",
    export_formats = c("pdf","svg")
) {
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("tibble", quietly = TRUE)
  requireNamespace("rlang", quietly = TRUE)

  need <- c(pathway_col, category_col, gene_col)
  miss <- setdiff(need, names(df))
  if (length(miss)) {
    stop("plot_pathway_gene_heatmap(): missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }

  # --------------- Prepare data ---------------
  df_valid <- df[!is.na(df[[pathway_col]]), , drop = FALSE]

  # distinct pathway/category; sort by category then pathway
  pathway_order_df <- df_valid |>
    dplyr::select(dplyr::all_of(c(pathway_col, category_col))) |>
    dplyr::distinct() |>
    dplyr::arrange(.data[[category_col]], .data[[pathway_col]])

  # category offsets with a gap of 1 line between categories
  cat_counts <- pathway_order_df |>
    dplyr::count(.data[[category_col]], name = "n") |>
    dplyr::arrange(.data[[category_col]]) |>
    dplyr::mutate(offset = dplyr::lag(cumsum(n + 1L), default = 0L))

  pathway_order_df <- pathway_order_df |>
    dplyr::left_join(cat_counts, by = rlang::set_names(category_col, category_col)) |>
    dplyr::group_by(.data[[category_col]]) |>
    dplyr::mutate(compact_y_pos = (dplyr::row_number() - 1L) + .data$offset) |>
    dplyr::ungroup()

  # membership pairs (present = 1), keep distinct and map y position + category
  heat_pairs <- df_valid |>
    dplyr::select(dplyr::all_of(c(pathway_col, gene_col))) |>
    dplyr::distinct() |>
    dplyr::left_join(pathway_order_df, by = rlang::set_names(pathway_col, pathway_col))

  # gene order (sorted) and numeric x positions
  gene_order <- sort(unique(heat_pairs[[gene_col]]))
  gene_x_map <- stats::setNames(seq_along(gene_order) - 1L, gene_order)  # 0..N-1 for easy tick math
  heat_pairs$x_pos <- unname(gene_x_map[heat_pairs[[gene_col]]])

  # --------------- Palette (by category) ---------------
  bio_categories <- sort(unique(pathway_order_df[[category_col]]))
  if (is.null(palette)) {
    default_palette <- c(
      "#DF8F44", "#00A1D5", "#B24745", "#79AF97", "#6A6599",
      "#374E55", "#80796B", "#AA4499", "#117733", "#999933", "#882255"
    )
    palette <- default_palette[seq_len(min(length(default_palette), length(bio_categories)))]
  }
  if (is.null(names(palette))) names(palette) <- bio_categories
  cat_pal <- palette
  # ensure palette has entries for all categories
  if (!all(bio_categories %in% names(cat_pal))) {
    # recycle if unnamed/short
    extra_cols <- rep(cat_pal, length.out = length(bio_categories))
    names(extra_cols) <- bio_categories
    cat_pal <- extra_cols
  }

  # --------------- Scales and separators ---------------
  # y breaks/labels from pathway_order_df (reverse so first at top)
  y_breaks <- pathway_order_df$compact_y_pos
  y_labels <- pathway_order_df[[pathway_col]]

  # category divider positions: end row index per category
  cat_ends <- pathway_order_df |>
    dplyr::group_by(.data[[category_col]]) |>
    dplyr::summarise(end_y = max(.data$compact_y_pos), .groups = "drop")

  max_y <- max(pathway_order_df$compact_y_pos, na.rm = TRUE)
  export_height <- (max_y + 2) * 0.5
  export_width  <- length(gene_order) * 0.6

  # --------------- Build plot ---------------
  # Use geom_tile with y as numeric positions and x as numeric 0..N-1 (so we can draw gridlines easily)
  p <- ggplot2::ggplot(heat_pairs, ggplot2::aes(x = .data$x_pos, y = .data$compact_y_pos)) +
    ggplot2::geom_tile(
      ggplot2::aes(fill = .data[[category_col]]),
      width = 0.98, height = 0.98, color = "black", linewidth = 0.3
    ) +
    # Reverse y so first is at the top; label with pathway names
    ggplot2::scale_y_reverse(breaks = y_breaks, labels = y_labels, expand = c(0, 0)) +
    # x axis with ticks at 0..N-1 and gene labels rotated
    ggplot2::scale_x_continuous(
      breaks = seq_along(gene_order) - 1L,
      labels = gene_order,
      expand = c(0, 0)
    ) +
    ggplot2::scale_fill_manual(values = cat_pal, name = "Category") +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 16, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 10),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.6)
    )

  # gridlines
  # vertical lines at boundaries
  vgrid <- tibble::tibble(x = seq(0, length(gene_order), by = 1))
  hgrid <- tibble::tibble(y = seq(-0.5, max_y + 1.5, by = 1))
  p <- p +
    ggplot2::geom_vline(data = vgrid, ggplot2::aes(xintercept = .data$x), color = "grey70", linewidth = 0.3) +
    ggplot2::geom_hline(data = hgrid, ggplot2::aes(yintercept = .data$y), color = "grey70", linewidth = 0.3)

  # dashed separators between categories (at end_y + 0.5)
  if (nrow(cat_ends)) {
    p <- p + ggplot2::geom_hline(
      data = cat_ends,
      ggplot2::aes(yintercept = .data$end_y + 0.5),
      linetype = "dashed", linewidth = 0.6, color = "grey50"
    )
  }

  # optional theme overrides
  if (!is.null(style_config) && length(style_config)) {
    p <- p + do.call(ggplot2::theme, style_config)
  }

  # --------------- Export ---------------
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
