#' Plot CPMM trajectory
#'
#' Fits a change-point mixed model with a fixed changepoint at 0:
#'   protein ~ before_onset + after_onset + covariates + (1 | SUBID)
#'
#' Plots:
#'   • Global fitted CPMM trajectory (population-level prediction)
#'   • Raw subject-level scatter points (outline only)
#'   • Vertical onset line at 0
#'
#' @param df_status_change Data frame of status-change (converter) subjects.
#' @param protein Character; protein column name.
#' @param covariates Character vector of covariates. Default c("SEX","BASELINE_AGE").
#' @param years_since_onset_col Column name. Default "years_since_onset".
#' @param subject_id_col Subject ID column. Default "SUBID".
#' @param main_color Line/point color. Default muted red.
#' @param export Logical; save plot?
#' @param export_dir Output directory.
#' @param export_formats File formats.
#'
#' @return A ggplot object.
#' @export
plot_cpmm <- function(
    df_status_change,
    protein,
    covariates = c("SEX", "BASELINE_AGE"),
    years_since_onset_col = "years_since_onset",
    subject_id_col = "SUBID",
    main_color = "#B24745",
    export = FALSE,
    export_dir = "Figures",
    export_formats = c("pdf", "svg")
) {
  requireNamespace("lmerTest", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("tibble", quietly = TRUE)

  # ---- Add piecewise terms ----
  df <- df_status_change |>
    dplyr::mutate(
      before_onset = pmax(0, - .data[[years_since_onset_col]]),
      after_onset  = pmax(0,   .data[[years_since_onset_col]])
    )

  # ---- Mean covariates (status-change cohort) ----
  mean_cov <- lapply(covariates, function(cov) {
    v <- df[[cov]]
    if (is.numeric(v)) mean(v, na.rm = TRUE)
    else names(sort(table(v), decreasing = TRUE))[1]
  })
  names(mean_cov) <- covariates

  # ---- Model ----
  rhs <- paste(c("before_onset", "after_onset", covariates), collapse = " + ")
  fml <- as.formula(
    paste0(protein, " ~ ", rhs, " + (1 | ", subject_id_col, ")")
  )

  fit <- lmerTest::lmer(fml, data = df, REML = TRUE)

  # ---- Predict population-level trajectory ----
  df_pred <- df
  for (cov in names(mean_cov)) {
    df_pred[[cov]] <- mean_cov[[cov]]
  }

  df_pred$.pred <- predict(fit, newdata = df_pred, re.form = NA)

  trend <- df_pred |>
    dplyr::group_by(.data[[years_since_onset_col]]) |>
    dplyr::summarise(y = mean(.data$.pred, na.rm = TRUE), .groups = "drop")

  # ---- X ticks every 2 years (GTEx style) ----
  min_x <- floor(min(df[[years_since_onset_col]]) / 2) * 2
  max_x <- ceiling(max(df[[years_since_onset_col]]) / 2) * 2

  # ---- Plot ----
  p <- ggplot2::ggplot() +
    # Scatter (outline only)
    ggplot2::geom_point(
      data = df,
      ggplot2::aes(
        x = .data[[years_since_onset_col]],
        y = .data[[protein]]
      ),
      shape = 21,
      fill = "white",
      color = main_color,
      alpha = 0.6,
      size = 2
    ) +
    # Fitted CPMM trajectory
    ggplot2::geom_line(
      data = trend,
      ggplot2::aes(
        x = .data[[years_since_onset_col]],
        y = .data$y
      ),
      color = main_color,
      linewidth = 1.5
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.8) +
    ggplot2::scale_x_continuous(breaks = seq(min_x, max_x, by = 2)) +
    ggplot2::labs(
      x = "Years Since Onset",
      y = "Abundance level",
      title = protein
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(
        linetype = "dashed",
        color = "#E0E0E0",
        linewidth = 0.8
      ),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(linewidth = 1.2),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
    )

  # ---- Export ----
  if (isTRUE(export)) {
    dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
    for (fmt in export_formats) {
      ggplot2::ggsave(
        file.path(export_dir, paste0(protein, "_CPMM.", fmt)),
        p,
        width = 6,
        height = 4,
        dpi = 300
      )
    }
  }

  p
}
