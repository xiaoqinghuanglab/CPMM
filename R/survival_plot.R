#' Plot Kaplan-Meier curves by baseline expression group
#'
#' Plots Kaplan-Meier survival curves comparing High versus Low baseline
#' expression groups using age as the time axis.
#'
#' This function expects the subject-level output from
#' \code{make_expression_survival_df()}.
#'
#' @param surv_df Subject-level survival data frame.
#' @param protein_name Character; protein name to show in the plot title.
#' @param time_col Age-based survival time column. Default is "age_time".
#' @param event_col Event indicator column. Default is "event".
#' @param group_col Expression group column. Default is "expr_group".
#' @param save Logical; whether to save the plot. Default is FALSE.
#' @param save_path Directory to save the plot. Default is "Figures/survival".
#' @param file_prefix File prefix. Default is "KM".
#' @param palette Named color vector for Low and High groups.
#'
#' @return A ggplot object.
#' @export
plot_expression_km <- function(
    surv_df,
    protein_name,
    time_col = "age_time",
    event_col = "event",
    group_col = "expr_group",
    save = FALSE,
    save_path = "Figures/survival",
    file_prefix = "KM",
    palette = c(
      Low = "#00A1D5",
      High = "#B24745"
    )
) {

  pkgs <- c("survival", "ggplot2", "dplyr", "tibble")

  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required.", p), call. = FALSE)
    }
  }

  if (!is.data.frame(surv_df)) {
    stop("surv_df must be a data frame.", call. = FALSE)
  }

  required_cols <- c(time_col, event_col, group_col)
  missing_cols <- setdiff(required_cols, names(surv_df))

  if (length(missing_cols) > 0) {
    stop(
      sprintf(
        "surv_df is missing required columns: %s",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  plot_df <- surv_df |>
    dplyr::filter(
      !is.na(.data[[time_col]]),
      !is.na(.data[[event_col]]),
      !is.na(.data[[group_col]])
    ) |>
    dplyr::mutate(
      .time = as.numeric(.data[[time_col]]),
      .event = as.integer(.data[[event_col]]),
      .group = as.character(.data[[group_col]])
    ) |>
    dplyr::filter(.data$.group %in% c("Low", "High"))

  if (nrow(plot_df) == 0) {
    stop(
      sprintf("No data available for %s using age_time.", protein_name),
      call. = FALSE
    )
  }

  available_groups <- unique(plot_df$.group)

  if (length(available_groups) < 2) {
    stop(
      "Both Low and High groups are required for Kaplan-Meier comparison.",
      call. = FALSE
    )
  }

  plot_df$.group <- factor(
    plot_df$.group,
    levels = c("Low", "High")
  )

  surv_formula <- stats::as.formula("survival::Surv(.time, .event) ~ .group")

  sf <- survival::survfit(
    formula = surv_formula,
    data = plot_df
  )

  lr <- survival::survdiff(
    formula = surv_formula,
    data = plot_df
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

  ## ---- Convert survfit object to plotting data ----
  sf_summary <- summary(sf)

  km_df <- tibble::tibble(
    time = sf_summary$time,
    survival = sf_summary$surv,
    lower = sf_summary$lower,
    upper = sf_summary$upper,
    strata = sf_summary$strata
  )

  km_df <- km_df |>
    dplyr::mutate(
      lower = ifelse(is.na(.data$lower), .data$survival, .data$lower),
      upper = ifelse(is.na(.data$upper), .data$survival, .data$upper),
      lower = pmax(.data$lower, 0),
      upper = pmin(.data$upper, 1),
      expr_group = sub("^\\.group=", "", .data$strata),
      expr_group = factor(.data$expr_group, levels = c("Low", "High"))
    ) |>
    dplyr::filter(
      !is.na(.data$time),
      !is.na(.data$survival),
      !is.na(.data$lower),
      !is.na(.data$upper)
    )

  x_min <- min(plot_df$.time, na.rm = TRUE)
  x_max <- max(plot_df$.time, na.rm = TRUE)
  x_limits <- c(x_min - 1, x_max + 1)

  ## ---- Plot ----
  p <- ggplot2::ggplot(
    km_df,
    ggplot2::aes(
      x = .data$time,
      y = .data$survival,
      color = .data$expr_group,
      fill = .data$expr_group
    )
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = .data$lower,
        ymax = .data$upper
      ),
      alpha = 0.18,
      color = NA
    ) +
    ggplot2::geom_step(
      linewidth = 1.2
    ) +
    ggplot2::annotate(
      "text",
      x = x_limits[1] + 0.68 * diff(x_limits),
      y = 0.08,
      label = p_text,
      size = 5
    ) +
    ggplot2::scale_color_manual(
      values = palette[c("Low", "High")],
      drop = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = palette[c("Low", "High")],
      drop = FALSE
    ) +
    ggplot2::coord_cartesian(
      xlim = x_limits,
      ylim = c(0, 1)
    ) +
    ggplot2::labs(
      title = protein_name,
      x = "Age",
      y = "Survival Probability",
      color = "Expression group",
      fill = "Expression group"
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
      legend.position = "right",
      legend.title = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 12),
      legend.key = ggplot2::element_blank()
    )

  if (isTRUE(save)) {
    if (!requireNamespace("svglite", quietly = TRUE)) {
      stop("Package 'svglite' is required to save SVG files.", call. = FALSE)
    }

    dir.create(save_path, recursive = TRUE, showWarnings = FALSE)

    outfile <- file.path(
      save_path,
      paste0(file_prefix, "_", protein_name, ".svg")
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
