#' Overlay age-based Kaplan-Meier curves for two onset definitions
#'
#' Overlays High versus Low baseline-expression Kaplan-Meier curves from two
#' survival data frames. The background data frame typically represents the
#' original onset definition, while the main data frame represents the
#' onset-minus-one-year sensitivity analysis.
#'
#' This function uses age as the only supported time axis.
#'
#' @param surv_df_main Subject-level survival data frame for the main analysis,
#' usually onset minus 1 year.
#' @param surv_df_bg Subject-level survival data frame for the background
#' comparison, usually original onset.
#' @param protein_name Character; protein name shown in the plot title.
#' @param main_label Label for the main survival curves.
#' Default is "Onset - 1 year".
#' @param bg_label Label for the background survival curves.
#' Default is "Original onset".
#' @param cohort_label Optional cohort label shown in the title. Default is NULL.
#' @param time_col Age-based survival time column. Default is "age_time".
#' @param event_col Event indicator column. Default is "event".
#' @param group_col Expression group column. Default is "expr_group".
#' @param save Logical; whether to save the plot. Default is FALSE.
#' @param save_path Directory to save the plot. Default is "Figures/survival".
#' @param file_prefix File prefix. Default is "KM_overlay".
#' @param palette_main Named colors for the main curves.
#' @param palette_bg Named colors for the background curves.
#' @param show_main_ci Logical; whether to show confidence interval ribbons for
#' the main curves. Default is TRUE.
#'
#' @return A ggplot object.
#' @export
plot_expression_km_onset_overlay <- function(
    surv_df_main,
    surv_df_bg,
    protein_name,
    main_label = "Onset - 1 year",
    bg_label = "Original onset",
    cohort_label = NULL,
    time_col = "age_time",
    event_col = "event",
    group_col = "expr_group",
    save = FALSE,
    save_path = "Figures/survival",
    file_prefix = "KM_overlay",
    palette_main = c(
      Low = "#00A1D5",
      High = "#B24745"
    ),
    palette_bg = c(
      Low = "#A9D9EC",
      High = "#E8B8B7"
    ),
    show_main_ci = TRUE
) {

  pkgs <- c("survival", "ggplot2", "dplyr", "tibble")

  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required.", p), call. = FALSE)
    }
  }

  if (!is.data.frame(surv_df_main)) {
    stop("surv_df_main must be a data frame.", call. = FALSE)
  }

  if (!is.data.frame(surv_df_bg)) {
    stop("surv_df_bg must be a data frame.", call. = FALSE)
  }

  required_cols <- c(time_col, event_col, group_col)

  missing_main <- setdiff(required_cols, names(surv_df_main))
  if (length(missing_main) > 0) {
    stop(
      sprintf(
        "surv_df_main is missing required columns: %s",
        paste(missing_main, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  missing_bg <- setdiff(required_cols, names(surv_df_bg))
  if (length(missing_bg) > 0) {
    stop(
      sprintf(
        "surv_df_bg is missing required columns: %s",
        paste(missing_bg, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  prepare_survival_df <- function(df, timeline_label) {
    df |>
      dplyr::filter(
        !is.na(.data[[time_col]]),
        !is.na(.data[[event_col]]),
        !is.na(.data[[group_col]])
      ) |>
      dplyr::mutate(
        .time = as.numeric(.data[[time_col]]),
        .event = as.integer(.data[[event_col]]),
        .group = as.character(.data[[group_col]]),
        .timeline = timeline_label
      ) |>
      dplyr::filter(.data$.group %in% c("Low", "High"))
  }

  plot_main <- prepare_survival_df(surv_df_main, main_label)
  plot_bg <- prepare_survival_df(surv_df_bg, bg_label)

  if (nrow(plot_main) == 0) {
    stop(
      sprintf("No data available for %s using age_time.", protein_name),
      call. = FALSE
    )
  }

  main_groups <- unique(plot_main$.group)

  if (length(main_groups) < 2) {
    stop(
      "Both Low and High groups are required in surv_df_main for Kaplan-Meier comparison.",
      call. = FALSE
    )
  }

  plot_main$.group <- factor(plot_main$.group, levels = c("Low", "High"))

  if (nrow(plot_bg) > 0) {
    plot_bg$.group <- factor(plot_bg$.group, levels = c("Low", "High"))
  }

  make_km_df <- function(df, group_name, timeline_label) {
    group_df <- df[df$.group == group_name, , drop = FALSE]

    if (nrow(group_df) == 0) {
      return(NULL)
    }

    fit <- survival::survfit(
      survival::Surv(.time, .event) ~ 1,
      data = group_df
    )

    fit_summary <- summary(fit)

    if (length(fit_summary$time) == 0) {
      return(NULL)
    }

    curve_label <- sprintf("%s (%s)", group_name, timeline_label)

    km_df <- tibble::tibble(
      time = fit_summary$time,
      survival = fit_summary$surv,
      lower = fit_summary$lower,
      upper = fit_summary$upper,
      expr_group = group_name,
      timeline = timeline_label,
      curve_label = curve_label
    )

    km_df |>
      dplyr::mutate(
        lower = ifelse(is.na(.data$lower), .data$survival, .data$lower),
        upper = ifelse(is.na(.data$upper), .data$survival, .data$upper),
        lower = pmax(.data$lower, 0),
        upper = pmin(.data$upper, 1),
        expr_group = factor(.data$expr_group, levels = c("Low", "High"))
      ) |>
      dplyr::filter(
        !is.na(.data$time),
        !is.na(.data$survival),
        !is.na(.data$lower),
        !is.na(.data$upper)
      )
  }

  main_km <- dplyr::bind_rows(
    make_km_df(plot_main, "Low", main_label),
    make_km_df(plot_main, "High", main_label)
  )

  bg_km <- dplyr::bind_rows(
    make_km_df(plot_bg, "Low", bg_label),
    make_km_df(plot_bg, "High", bg_label)
  )

  if (nrow(main_km) == 0) {
    stop("No Kaplan-Meier estimates could be generated for surv_df_main.", call. = FALSE)
  }

  lr <- survival::survdiff(
    survival::Surv(.time, .event) ~ .group,
    data = plot_main
  )

  p_val <- stats::pchisq(
    lr$chisq,
    df = length(lr$n) - 1,
    lower.tail = FALSE
  )

  p_text <- if (is.na(p_val)) {
    "p = NA"
  } else {
    sprintf("p = %.4f", p_val)
  }

  combined_time <- c(plot_main$.time, plot_bg$.time)
  combined_time <- combined_time[!is.na(combined_time)]

  x_min <- min(combined_time, na.rm = TRUE)
  x_max <- max(combined_time, na.rm = TRUE)

  x_limits <- c(x_min - 1, x_max + 1)

  plot_title <- if (is.null(cohort_label)) {
    protein_name
  } else {
    sprintf("%s (%s)", protein_name, cohort_label)
  }

  low_bg_label <- sprintf("Low (%s)", bg_label)
  high_bg_label <- sprintf("High (%s)", bg_label)
  low_main_label <- sprintf("Low (%s)", main_label)
  high_main_label <- sprintf("High (%s)", main_label)

  line_colors <- c(
    setNames(palette_bg[["Low"]], low_bg_label),
    setNames(palette_bg[["High"]], high_bg_label),
    setNames(palette_main[["Low"]], low_main_label),
    setNames(palette_main[["High"]], high_main_label)
  )

  legend_order <- c(
    low_bg_label,
    high_bg_label,
    low_main_label,
    high_main_label
  )

  if (nrow(bg_km) > 0) {
    bg_km$curve_label <- factor(bg_km$curve_label, levels = legend_order)
  }
  main_km$curve_label <- factor(main_km$curve_label, levels = legend_order)

  p <- ggplot2::ggplot()

  # Background/original onset curves: lighter colors, no CI
  if (nrow(bg_km) > 0) {
    p <- p +
      ggplot2::geom_step(
        data = bg_km,
        ggplot2::aes(
          x = .data$time,
          y = .data$survival,
          color = .data$curve_label
        ),
        linewidth = 1.2,
        alpha = 1
      )
  }

  # Main/onset-minus-1 CI ribbon: colored by Low/High only, no legend entry
  if (isTRUE(show_main_ci)) {
    p <- p +
      ggplot2::geom_ribbon(
        data = main_km,
        ggplot2::aes(
          x = .data$time,
          ymin = .data$lower,
          ymax = .data$upper,
          fill = .data$expr_group
        ),
        alpha = 0.18,
        color = NA,
        show.legend = FALSE
      )
  }

  # Main/onset-minus-1 curves: darker colors
  p <- p +
    ggplot2::geom_step(
      data = main_km,
      ggplot2::aes(
        x = .data$time,
        y = .data$survival,
        color = .data$curve_label
      ),
      linewidth = 1.35,
      alpha = 1
    ) +
    ggplot2::annotate(
      "text",
      x = x_limits[1] + 0.60 * diff(x_limits),
      y = 0.08,
      label = p_text,
      size = 5
    ) +
    ggplot2::scale_color_manual(
      values = line_colors,
      breaks = legend_order,
      drop = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        Low = palette_main[["Low"]],
        High = palette_main[["High"]]
      ),
      drop = FALSE
    ) +
    ggplot2::coord_cartesian(
      xlim = x_limits,
      ylim = c(0, 1)
    ) +
    ggplot2::labs(
      title = plot_title,
      x = "Age",
      y = "Survival Probability",
      color = NULL,
      fill = NULL
    ) +
    ggplot2::theme_classic(
      base_size = 14,
      base_family = "serif"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = 16,
        margin = ggplot2::margin(b = 10)
      ),
      axis.title = ggplot2::element_text(size = 15),
      axis.text = ggplot2::element_text(size = 13),
      axis.line = ggplot2::element_line(linewidth = 0.8),
      axis.ticks = ggplot2::element_line(linewidth = 0.7),
      panel.grid.major.y = ggplot2::element_line(
        linetype = "dashed",
        color = "#DDDDDD",
        linewidth = 0.5
      ),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = c(0.05, 0.10),
      legend.justification = c(0, 0),
      legend.background = ggplot2::element_rect(
        fill = "white",
        color = "#D0D0D0",
        linewidth = 0.4
      ),
      legend.text = ggplot2::element_text(size = 11),
      legend.key = ggplot2::element_blank()
    )

  if (isTRUE(save)) {
    if (!requireNamespace("svglite", quietly = TRUE)) {
      stop("Package 'svglite' is required to save SVG files.", call. = FALSE)
    }

    dir.create(save_path, recursive = TRUE, showWarnings = FALSE)

    outfile <- file.path(
      save_path,
      paste0(file_prefix, "_", protein_name, "_age_time.svg")
    )

    ggplot2::ggsave(
      filename = outfile,
      plot = p,
      width = 8,
      height = 7,
      dpi = 300,
      device = "svg"
    )
  }

  p
}
