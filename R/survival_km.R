#' Flexible KM plot for biomarker threshold crossing
#'
#' Allows user-defined 2- or 3-group Kaplan–Meier comparisons using
#' arbitrary cohorts or filtered subsets. Performs a multivariate
#' log-rank test across all supplied groups and shows an at-risk table.
#'
#' @param biomarker_name Character; biomarker column name.
#' @param threshold Numeric; threshold defining event occurrence.
#' @param groups Named list defining groups. Each element must contain:
#'   - df: data.frame
#'   - filter: (optional) logical expression evaluated within df
#' @param save Logical; save SVG output?
#' @param save_path Directory to save plot.
#' @param palette Optional named vector of colors.
#' @param time_points Numeric vector of x-axis ticks.
#'
#' @return A patchwork/ggplot object.
#' @export
plot_km_with_threshold <- function(
    biomarker_name,
    threshold,
    groups,
    save = FALSE,
    save_path = "Figures",
    palette = NULL,
    time_points = seq(-10, 4, by = 2)
) {
  requireNamespace("survival", quietly = TRUE)
  requireNamespace("survminer", quietly = TRUE)
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("patchwork", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("rlang", quietly = TRUE)

  if (length(groups) < 2 || length(groups) > 3) {
    stop("Provide exactly 2 or 3 groups.", call. = FALSE)
  }

  # ---- Build event data ----
  event_list <- lapply(names(groups), function(g) {
    spec <- groups[[g]]
    df <- spec$df

    if (!is.null(spec$filter)) {
      df <- dplyr::filter(df, !!spec$filter)
    }

    compute_event_df(df, threshold, biomarker_name, g)
  })

  event_df <- dplyr::bind_rows(event_list)
  if (!nrow(event_df)) stop("No event data generated.", call. = FALSE)

  event_df$group <- factor(event_df$group, levels = names(groups))

  # ---- Time shifting ----
  min_time <- suppressWarnings(min(c(event_df$time, time_points), na.rm = TRUE))
  offset <- if (min_time < 0) -min_time else 0
  event_df$time_shifted <- event_df$time + offset
  ticks_shifted <- time_points + offset

  # ---- KM fit ----
  sf <- survival::survfit(
    survival::Surv(time_shifted, event) ~ group,
    data = event_df
  )

  # ---- Log-rank (multigroup) ----
  sd <- survival::survdiff(
    survival::Surv(time_shifted, event) ~ group,
    data = event_df
  )
  p_val <- stats::pchisq(sd$chisq, df = length(sd$n) - 1, lower.tail = FALSE)
  p_text <- sprintf("p = %.3g", p_val)

  # ---- Palette ----
  if (is.null(palette)) {
    palette <- c(
      "#00A1D5", "#B24745", "#DF8F44"
    )[seq_along(groups)]
    names(palette) <- names(groups)
  }

  strata_levels <- sub("^group=", "", names(sf$strata))
  pal_vec <- unname(palette[strata_levels])

  # ---- Plot ----
  g <- survminer::ggsurvplot(
    sf,
    data = event_df,
    conf.int = TRUE,
    risk.table = TRUE,
    risk.table.height = 0.25,
    break.time.by = diff(time_points)[1],
    xlim = range(ticks_shifted),
    palette = pal_vec,
    legend.title = "Group",
    ggtheme = ggplot2::theme_minimal(base_size = 14)
  )

  g$plot <- g$plot +
    ggplot2::geom_vline(xintercept = offset, linetype = "dashed", color = "grey60") +
    ggplot2::scale_x_continuous(
      breaks = ticks_shifted,
      labels = time_points
    ) +
    ggplot2::annotate(
      "text",
      x = max(ticks_shifted) * 0.6,
      y = 0.1,
      label = p_text,
      size = 5
    ) +
    ggplot2::labs(
      title = sprintf("%s ≥ %.2f", biomarker_name, threshold),
      x = "Years to Onset",
      y = "Proportion Not Yet Crossed"
    )

  if (!is.null(g$table)) {
    g$table <- g$table +
      ggplot2::scale_x_continuous(
        breaks = ticks_shifted,
        labels = time_points
      )
  }

  combined <- g$plot / g$table + patchwork::plot_layout(heights = c(3, 1))

  if (save) {
    dir.create(save_path, showWarnings = FALSE, recursive = TRUE)
    ggplot2::ggsave(
      file.path(save_path, paste0("KM_", biomarker_name, ".svg")),
      combined,
      width = 7,
      height = 9,
      dpi = 300
    )
  }

  combined
}
