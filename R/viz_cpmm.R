#' Plot CPMM trajectories (status-change, normal-only, abnormal-only)
#'
#' Fits 3 random-intercept LMMs for a given protein with a fixed changepoint at 0:
#' \itemize{
#'   \item Status-change:  \code{protein ~ before_onset + after_onset + covariates + (1|SUBID)}
#'   \item Normal-only:    \code{protein ~ before_onset + covariates + (1|SUBID)}
#'   \item Abnormal-only:  \code{protein ~ after_onset  + covariates + (1|SUBID)}
#' }
#' It predicts each model at the mean of the provided covariates (from the status-change
#' cohort), averages predictions within each \code{years_since_onset} to get a smooth
#' "global" trend, and plots all three curves. A dashed line marks onset (0).
#'
#' @param df_status_change Data frame of status-change patients.
#' @param df_normal Data frame of normal-only patients.
#' @param df_abnormal Data frame of abnormal-only patients.
#' @param protein Character; protein column to analyze.
#' @param covariates Character vector of covariate columns. Default \code{c("SEX","BASELINE_AGE")}.
#' @param years_since_onset_col Column name for years since onset. Default \code{"years_since_onset"}.
#' @param subject_id_col Subject ID column (random intercept). Default \code{"SUBID"}.
#' @param group_colors Named vector for colors with keys \code{"status_change","normal","abnormal"}.
#' @param group_labels Named vector for legend labels with same keys as \code{group_colors}.
#' @param style_config Optional list of ggplot2 theme overrides (e.g. \code{list(panel.grid = element_blank())}).
#' @param export Logical; save the plot to \code{export_dir} in \code{export_formats}.
#' @param export_dir Directory to save into. Default \code{"Figures"}.
#' @param export_formats Character vector of file extensions (e.g. \code{c("pdf","svg")}). Default \code{c("pdf","svg")}.
#'
#' @return A ggplot object.
#' @export
plot_cpmm <- function(
    df_status_change,
    df_normal,
    df_abnormal,
    protein,
    covariates = c("SEX", "BASELINE_AGE"),
    years_since_onset_col = "years_since_onset",
    subject_id_col = "SUBID",
    group_colors = NULL,
    group_labels = NULL,
    style_config = NULL,
    export = FALSE,
    export_dir = "Figures",
    export_formats = c("pdf", "svg")
) {
  # deps
  requireNamespace("lme4", quietly = TRUE)
  requireNamespace("lmerTest", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("tibble", quietly = TRUE)

  # Defaults
  if (is.null(group_colors)) {
    group_colors <- c(
      status_change = "#DF8F44",
      normal        = "#00A1D5",
      abnormal      = "#B24745"
    )
  }
  if (is.null(group_labels)) {
    group_labels <- c(
      status_change = "Status Change",
      normal        = "Normal",
      abnormal      = "Abnormal"
    )
  }

  # mean covariates computed from status-change cohort (only those present)
  mean_covariates <- lapply(covariates, function(cov) {
    if (cov %in% names(df_status_change)) {
      # numeric mean for numeric; for factor/char use the most frequent level
      v <- df_status_change[[cov]]
      if (is.numeric(v)) {
        mean(v, na.rm = TRUE)
      } else {
        tbl <- sort(table(v), decreasing = TRUE)
        names(tbl)[1]
      }
    } else {
      NULL
    }
  })
  names(mean_covariates) <- covariates
  mean_covariates <- Filter(Negate(is.null), mean_covariates)

  # helper: add piecewise terms
  add_pw <- function(df) {
    if (is.null(df) || !nrow(df)) return(df)
    yso <- df[[years_since_onset_col]]
    df$before_onset <- pmax(0, -yso)
    df$after_onset  <- pmax(0,  yso)
    df
  }

  # helper: fit one group, predict at mean covariates, average by years_since_onset
  fit_and_trend <- function(df, rhs_terms, key) {
    if (is.null(df) || !nrow(df)) return(NULL)
    if (!all(c(subject_id_col, protein, years_since_onset_col) %in% names(df))) return(NULL)

    df <- add_pw(df)

    fml <- stats::as.formula(
      paste0(protein, " ~ ", rhs_terms, " + (1|", subject_id_col, ")")
    )

    fit <- tryCatch(
      lmerTest::lmer(fml, data = df, REML = TRUE),
      error = function(e) NULL,
      warning = function(w) suppressWarnings(lmerTest::lmer(fml, data = df, REML = TRUE))
    )
    if (is.null(fit)) return(NULL)

    df_pred <- df
    # set covariates to mean values where available
    for (cov in names(mean_covariates)) {
      if (cov %in% names(df_pred)) {
        val <- mean_covariates[[cov]]
        if (is.numeric(df_pred[[cov]]) && is.numeric(val)) {
          df_pred[[cov]] <- val
        } else {
          df_pred[[cov]] <- as.character(val)
        }
      }
    }

    # population-level prediction (exclude random effects): re.form = NA
    pred <- tryCatch(
      predict(fit, newdata = df_pred, allow.new.levels = TRUE, re.form = NA),
      error = function(e) rep(NA_real_, nrow(df_pred))
    )

    tr <- dplyr::tibble(
      !!years_since_onset_col := df_pred[[years_since_onset_col]],
      .pred = pred
    ) |>
      dplyr::group_by(.data[[years_since_onset_col]]) |>
      dplyr::summarise(y = mean(.data$.pred, na.rm = TRUE), .groups = "drop") |>
      dplyr::mutate(group = key)

    tr
  }

  rhs_status  <- paste(c("before_onset", "after_onset", covariates), collapse = " + ")
  rhs_normal  <- paste(c("before_onset", covariates), collapse = " + ")
  rhs_abnorm  <- paste(c("after_onset",  covariates), collapse = " + ")

  trend_status <- fit_and_trend(df_status_change, rhs_status, "status_change")
  trend_normal <- fit_and_trend(df_normal,        rhs_normal, "normal")
  trend_abnorm <- fit_and_trend(df_abnormal,      rhs_abnorm, "abnormal")

  trends <- dplyr::bind_rows(trend_status, trend_normal, trend_abnorm)
  if (!nrow(trends)) {
    stop("plot_cpmm(): no trends could be computed — check inputs/columns.", call. = FALSE)
  }

  # dynamic x ticks (rounded to step 5)
  all_x <- c(
    df_status_change[[years_since_onset_col]],
    df_normal[[years_since_onset_col]],
    df_abnormal[[years_since_onset_col]]
  )
  all_x <- all_x[is.finite(all_x)]
  if (length(all_x)) {
    min_x <- floor(min(all_x, na.rm = TRUE) / 5) * 5
    max_x <- ceiling(max(all_x, na.rm = TRUE) / 5) * 5
    breaks_x <- seq(min_x, max_x, by = 5)
  } else {
    breaks_x <- ggplot2::waiver() # let ggplot choose
  }

  # factor order for legend & color mapping
  trends$group <- factor(trends$group, levels = c("status_change","normal","abnormal"))

  # build plot
  p <- ggplot2::ggplot(
    trends,
    ggplot2::aes(x = .data[[years_since_onset_col]], y = .data$y, color = .data$group)
  ) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.6, color = "black") +
    ggplot2::scale_color_manual(
      values = group_colors[levels(trends$group)],
      breaks = levels(trends$group),
      labels = group_labels[levels(trends$group)]
    ) +
    ggplot2::labs(
      x = "Years Since Onset",
      y = paste0(protein, " Expression"),
      title = paste0(protein, " Trajectory Aligned to Onset"),
      color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(linetype = "dotted", linewidth = 0.3),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(linewidth = 0.6, color = "black"),
      legend.position = "right",
      plot.title = ggplot2::element_text(face = "bold")
    )

  if (!is.null(style_config) && length(style_config)) {
    # allow simple theme overrides via list of elements, e.g., list(legend.position = "bottom")
    p <- p + do.call(ggplot2::theme, style_config)
  }

  if (!is.null(breaks_x) && !identical(breaks_x, ggplot2::waiver())) {
    p <- p + ggplot2::scale_x_continuous(breaks = breaks_x)
  }

  if (isTRUE(export)) {
    if (!dir.exists(export_dir)) dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
    for (fmt in export_formats) {
      outfile <- file.path(export_dir, paste0(protein, "_cpmm.", fmt))
      ggplot2::ggsave(outfile, p, width = 6, height = 4, dpi = 300)
    }
  }

  p
}
