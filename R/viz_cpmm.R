#' Plot CPMM trajectories across one or more cohorts
#'
#' Fits and overlays cohort-specific CPMM trajectories for one or more cohorts.
#' Each cohort is modeled separately using:
#'
#'   protein ~ before_onset + after_onset + covariates + (1 | subject_id)
#'
#' The fitted trajectory represents the population-level fixed-effect prediction.
#'
#' @param cohort_dfs Named list of cohort data frames. Each list element should
#' contain the protein column, subject ID column, years-since-onset column, and
#' covariates.
#' @param protein Character; protein column name to plot.
#' @param covariates Character vector of fixed-effect covariates.
#' Default is c("SEX", "BASELINE_AGE").
#' @param subject_id_col Subject ID column name. Default is "SUBID".
#' @param years_since_onset_col Years-since-onset column name.
#' Default is "years_since_onset".
#' @param cohort_colors Optional named character vector of colors. Names should
#' match names(cohort_dfs). If NULL, ggplot default colors are used.
#' @param export Logical; whether to save the plot. Default is FALSE.
#' @param outpath Optional output path. Default is NULL.
#' @param remove_outliers_for_plot Logical; remove outliers from scatter and
#' spaghetti layers only. Model fitting is always done on the full available data.
#' Default is FALSE.
#' @param outlier_method Outlier method. Currently only "iqr" is supported.
#' @param show_spaghetti Logical; whether to show subject-level trajectories.
#' Default is TRUE.
#' @param n_points Number of prediction points per cohort. Default is 100.
#'
#' @return A ggplot object.
#' @export
plot_cpmm <- function(
    cohort_dfs,
    protein,
    covariates = c("SEX", "BASELINE_AGE"),
    subject_id_col = "SUBID",
    years_since_onset_col = "years_since_onset",
    cohort_colors = NULL,
    export = FALSE,
    outpath = NULL,
    remove_outliers_for_plot = FALSE,
    outlier_method = "iqr",
    show_spaghetti = TRUE,
    n_points = 100
) {

  pkgs <- c("lme4", "lmerTest", "ggplot2", "dplyr", "tibble")

  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required.", p), call. = FALSE)
    }
  }

  if (!is.list(cohort_dfs)) {
    stop("cohort_dfs must be a named list of data frames.", call. = FALSE)
  }

  if (is.null(names(cohort_dfs)) || any(names(cohort_dfs) == "")) {
    names(cohort_dfs) <- paste0("Cohort_", seq_along(cohort_dfs))
  }

  if (!outlier_method %in% c("iqr")) {
    stop("Currently, only outlier_method = 'iqr' is supported.", call. = FALSE)
  }

  quote_var <- function(x) {
    paste0("`", gsub("`", "\\\\`", x), "`")
  }

  get_plot_df <- function(df, protein) {
    plot_df <- df

    if (!isTRUE(remove_outliers_for_plot)) {
      return(plot_df)
    }

    q1 <- stats::quantile(plot_df[[protein]], 0.25, na.rm = TRUE)
    q3 <- stats::quantile(plot_df[[protein]], 0.75, na.rm = TRUE)
    iqr <- q3 - q1

    lower_bound <- q1 - 1.5 * iqr
    upper_bound <- q3 + 1.5 * iqr

    plot_df[
      plot_df[[protein]] >= lower_bound &
        plot_df[[protein]] <= upper_bound,
      ,
      drop = FALSE
    ]
  }

  get_covariate_reference <- function(df, covariate) {
    v <- df[[covariate]]

    if (is.numeric(v)) {
      return(mean(v, na.rm = TRUE))
    }

    v_nonmissing <- v[!is.na(v)]

    if (length(v_nonmissing) == 0) {
      return(NA)
    }

    names(sort(table(v_nonmissing), decreasing = TRUE))[1]
  }

  fit_lmm_and_predict <- function(df, protein) {
    if (nrow(df) == 0) {
      return(NULL)
    }

    df <- df |>
      dplyr::mutate(
        before_onset = pmax(0, 0 - .data[[years_since_onset_col]]),
        after_onset = pmax(0, .data[[years_since_onset_col]] - 0)
      )

    rhs <- paste(
      c("before_onset", "after_onset", covariates),
      collapse = " + "
    )

    fml <- stats::as.formula(
      paste0(
        quote_var(protein),
        " ~ ",
        rhs,
        " + (1 | ",
        quote_var(subject_id_col),
        ")"
      )
    )

    fit <- tryCatch(
      suppressWarnings(
        lmerTest::lmer(
          formula = fml,
          data = df,
          REML = TRUE,
          na.action = stats::na.omit
        )
      ),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      return(NULL)
    }

    x_range <- seq(
      min(df[[years_since_onset_col]], na.rm = TRUE),
      max(df[[years_since_onset_col]], na.rm = TRUE),
      length.out = n_points
    )

    df_pred <- data.frame(
      years_since_onset = x_range,
      before_onset = pmax(0, 0 - x_range),
      after_onset = pmax(0, x_range - 0)
    )

    for (cov in covariates) {
      df_pred[[cov]] <- get_covariate_reference(df, cov)

      if (is.factor(df[[cov]])) {
        df_pred[[cov]] <- factor(df_pred[[cov]], levels = levels(df[[cov]]))
      }
    }

    df_pred$prediction <- stats::predict(
      fit,
      newdata = df_pred,
      re.form = NA,
      allow.new.levels = TRUE
    )

    fixed_terms <- stats::delete.response(stats::terms(fit, fixed.only = TRUE))

    x_mat <- stats::model.matrix(fixed_terms, data = df_pred)

    beta_names <- names(lme4::fixef(fit))

    x_mat <- x_mat[, beta_names, drop = FALSE]

    vc <- as.matrix(stats::vcov(fit))
    vc <- vc[beta_names, beta_names, drop = FALSE]

    pred_var <- diag(x_mat %*% vc %*% t(x_mat))
    pred_se <- sqrt(pmax(pred_var, 0))

    df_pred$lower_ci <- df_pred$prediction - 1.96 * pred_se
    df_pred$upper_ci <- df_pred$prediction + 1.96 * pred_se

    tibble::as_tibble(df_pred)
  }

  required_cols <- c(
    subject_id_col,
    years_since_onset_col,
    covariates,
    protein
  )

  plot_layers <- list()
  pred_layers <- list()

  for (cohort_name in names(cohort_dfs)) {
    df_full <- as.data.frame(cohort_dfs[[cohort_name]])

    missing_cols <- setdiff(required_cols, names(df_full))

    if (length(missing_cols) > 0) {
      warning(
        sprintf(
          "Skipping %s because these columns are missing: %s",
          cohort_name,
          paste(missing_cols, collapse = ", ")
        ),
        call. = FALSE
      )
      next
    }

    filter_expr <- stats::complete.cases(df_full[, required_cols, drop = FALSE])
    df_fit <- df_full[filter_expr, , drop = FALSE]

    if (nrow(df_fit) == 0) {
      next
    }

    df_plot <- get_plot_df(df_fit, protein)
    pred_data <- fit_lmm_and_predict(df_fit, protein)

    if (is.null(pred_data)) {
      next
    }

    df_plot$Cohort <- cohort_name
    pred_data$Cohort <- cohort_name

    names(df_plot)[names(df_plot) == years_since_onset_col] <- ".years_since_onset"
    names(df_plot)[names(df_plot) == subject_id_col] <- ".subject_id"

    plot_layers[[cohort_name]] <- tibble::as_tibble(df_plot)
    pred_layers[[cohort_name]] <- tibble::as_tibble(pred_data)
  }

  plot_df <- dplyr::bind_rows(plot_layers)
  pred_df <- dplyr::bind_rows(pred_layers)

  if (nrow(plot_df) == 0 || nrow(pred_df) == 0) {
    stop("No valid cohort data were available for plotting.", call. = FALSE)
  }

  min_x <- floor(min(plot_df$.years_since_onset, na.rm = TRUE) / 2) * 2
  max_x <- ceiling(max(plot_df$.years_since_onset, na.rm = TRUE) / 2) * 2

  p <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = pred_df,
      ggplot2::aes(
        x = .data$years_since_onset,
        ymin = .data$lower_ci,
        ymax = .data$upper_ci,
        fill = .data$Cohort
      ),
      alpha = 0.10,
      color = NA
    )

  if (isTRUE(show_spaghetti)) {
    p <- p +
      ggplot2::geom_line(
        data = plot_df |>
          dplyr::arrange(.data$Cohort, .data$.subject_id, .data$.years_since_onset),
        ggplot2::aes(
          x = .data$.years_since_onset,
          y = .data[[protein]],
          group = interaction(.data$Cohort, .data$.subject_id),
          color = .data$Cohort
        ),
        linewidth = 0.35,
        alpha = 0.08
      )
  }

  p <- p +
    ggplot2::geom_point(
      data = plot_df,
      ggplot2::aes(
        x = .data$.years_since_onset,
        y = .data[[protein]],
        color = .data$Cohort
      ),
      shape = 21,
      fill = "white",
      alpha = 0.75,
      size = 2.2,
      stroke = 0.7
    ) +
    ggplot2::geom_line(
      data = pred_df,
      ggplot2::aes(
        x = .data$years_since_onset,
        y = .data$prediction,
        color = .data$Cohort
      ),
      linewidth = 1.25
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      color = "black",
      linetype = "dashed",
      linewidth = 0.7
    ) +
    ggplot2::scale_x_continuous(breaks = seq(min_x, max_x, by = 2)) +
    ggplot2::labs(
      x = "Years since onset",
      y = "Abundance level",
      title = protein,
      color = NULL,
      fill = NULL
    ) +
    ggplot2::theme_classic(base_size = 14, base_family = "serif") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "plain",
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
        color = "#E0E0E0",
        linewidth = 0.5
      ),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "right",
      legend.text = ggplot2::element_text(size = 12),
      legend.key = ggplot2::element_blank()
    )

  if (!is.null(cohort_colors)) {
    p <- p +
      ggplot2::scale_color_manual(values = cohort_colors, drop = FALSE) +
      ggplot2::scale_fill_manual(values = cohort_colors, drop = FALSE)
  }

  if (isTRUE(export)) {
    if (is.null(outpath)) {
      outpath <- file.path("Figures", "cpmm_overlay", paste0(protein, "_overlay.svg"))
    }

    dir.create(dirname(outpath), recursive = TRUE, showWarnings = FALSE)

    ggplot2::ggsave(
      filename = outpath,
      plot = p,
      width = 12,
      height = 7,
      dpi = 300
    )
  }

  p
}
