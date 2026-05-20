#' Build subject-level survival data from baseline protein expression
#'
#' Creates a subject-level survival data frame for time-to-diagnosis analyses
#' using baseline protein expression groups.
#'
#' Subjects are handled as follows:
#' \itemize{
#'   \item CU subjects are censored at last follow-up.
#'   \item CI subjects are treated as events if include_ci = TRUE.
#'   \item Converter subjects are treated as events at onset age.
#' }
#'
#' Expression groups are defined using baseline protein abundance.
#'
#' @param cu Data frame of cognitively unimpaired or control subjects.
#' @param ci Data frame of cognitively impaired subjects.
#' @param conv Data frame of converter subjects.
#' @param protein Character; protein column name.
#' @param subject_col Subject ID column. Default is "SUBID".
#' @param age_col Visit age column. Default is "PROCEDURE_AGE".
#' @param onset_col Onset age column. Default is "ONSET_AGE".
#' @param expr_group_method Method for defining expression groups. One of
#' "median", "upper25", or "lower25". Default is "median".
#' @param include_ci Logical; whether to include CI subjects. Default is TRUE.
#' @param ci_event_time How to assign event age for CI subjects. One of
#' "baseline" or "onset". Default is "baseline".
#'
#' @return A tibble with one row per subject containing baseline expression,
#' event status, expression group, and age-based event/censoring time.
#' @export
make_expression_survival_df <- function(
    cu,
    ci,
    conv,
    protein,
    subject_col = "SUBID",
    age_col = "PROCEDURE_AGE",
    onset_col = "ONSET_AGE",
    expr_group_method = "median",
    include_ci = TRUE,
    ci_event_time = c("baseline", "onset")
) {

  pkgs <- c("tibble", "dplyr")

  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required.", p), call. = FALSE)
    }
  }

  ci_event_time <- match.arg(ci_event_time)

  if (!expr_group_method %in% c("median", "upper25", "lower25")) {
    stop(
      "expr_group_method must be one of 'median', 'upper25', or 'lower25'.",
      call. = FALSE
    )
  }

  check_required <- function(df, df_name, required_cols) {
    missing_cols <- setdiff(required_cols, names(df))

    if (length(missing_cols) > 0) {
      stop(
        sprintf(
          "%s is missing required columns: %s",
          df_name,
          paste(missing_cols, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  check_required(
    cu,
    "cu",
    c(subject_col, age_col, protein)
  )

  check_required(
    conv,
    "conv",
    c(subject_col, age_col, onset_col, protein)
  )

  if (isTRUE(include_ci)) {
    ci_required <- c(subject_col, age_col, protein)

    if (ci_event_time == "onset") {
      ci_required <- c(ci_required, onset_col)
    }

    check_required(ci, "ci", ci_required)
  }

  make_cu_rows <- function(df) {
    split_df <- split(df, df[[subject_col]])

    rows <- lapply(names(split_df), function(subid) {
      g <- split_df[[subid]]
      g <- g[order(g[[age_col]]), , drop = FALSE]

      if (nrow(g) == 0) {
        return(NULL)
      }

      baseline_age <- g[[age_col]][1]
      last_age <- g[[age_col]][nrow(g)]
      baseline_expr <- g[[protein]][1]

      onset_age <- if (onset_col %in% names(g)) {
        g[[onset_col]][1]
      } else {
        NA_real_
      }

      data.frame(
        SUBID = subid,
        source_group = "CU",
        baseline_age = as.numeric(baseline_age),
        end_age = as.numeric(last_age),
        onset_age = as.numeric(onset_age),
        baseline_expr = as.numeric(baseline_expr),
        event = 0L,
        stringsAsFactors = FALSE
      )
    })

    do.call(rbind, rows)
  }

  make_ci_rows <- function(df) {
    if (!isTRUE(include_ci)) {
      return(NULL)
    }

    split_df <- split(df, df[[subject_col]])

    rows <- lapply(names(split_df), function(subid) {
      g <- split_df[[subid]]
      g <- g[order(g[[age_col]]), , drop = FALSE]

      if (nrow(g) == 0) {
        return(NULL)
      }

      baseline_age <- g[[age_col]][1]
      baseline_expr <- g[[protein]][1]

      onset_age <- if (onset_col %in% names(g)) {
        g[[onset_col]][1]
      } else {
        baseline_age
      }

      end_age <- if (ci_event_time == "baseline") {
        baseline_age
      } else {
        onset_age
      }

      data.frame(
        SUBID = subid,
        source_group = "CI",
        baseline_age = as.numeric(baseline_age),
        end_age = as.numeric(end_age),
        onset_age = as.numeric(onset_age),
        baseline_expr = as.numeric(baseline_expr),
        event = 1L,
        stringsAsFactors = FALSE
      )
    })

    do.call(rbind, rows)
  }

  make_converter_rows <- function(df) {
    split_df <- split(df, df[[subject_col]])

    rows <- lapply(names(split_df), function(subid) {
      g <- split_df[[subid]]
      g <- g[order(g[[age_col]]), , drop = FALSE]

      if (nrow(g) == 0) {
        return(NULL)
      }

      baseline_age <- g[[age_col]][1]
      baseline_expr <- g[[protein]][1]
      onset_age <- g[[onset_col]][1]

      if (is.na(onset_age)) {
        return(NULL)
      }

      data.frame(
        SUBID = subid,
        source_group = "Converters",
        baseline_age = as.numeric(baseline_age),
        end_age = as.numeric(onset_age),
        onset_age = as.numeric(onset_age),
        baseline_expr = as.numeric(baseline_expr),
        event = 1L,
        stringsAsFactors = FALSE
      )
    })

    do.call(rbind, rows)
  }

  rows <- dplyr::bind_rows(
    make_cu_rows(cu),
    make_ci_rows(ci),
    make_converter_rows(conv)
  )

  surv_df <- tibble::as_tibble(rows)

  surv_df <- surv_df |>
    dplyr::filter(
      !is.na(.data$baseline_expr),
      !is.na(.data$baseline_age),
      !is.na(.data$end_age)
    )

  if (nrow(surv_df) == 0) {
    stop("No valid subjects remained after filtering.", call. = FALSE)
  }

  vals <- surv_df$baseline_expr

  cutoff <- switch(
    expr_group_method,
    median = stats::median(vals, na.rm = TRUE),
    upper25 = as.numeric(stats::quantile(vals, 0.75, na.rm = TRUE)),
    lower25 = as.numeric(stats::quantile(vals, 0.25, na.rm = TRUE))
  )

  surv_df <- surv_df |>
    dplyr::mutate(
      expr_group = ifelse(.data$baseline_expr >= cutoff, "High", "Low"),
      cutoff = cutoff,
      Protein = protein,
      age_time = .data$end_age
    )

  surv_df
}
