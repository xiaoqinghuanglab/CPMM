test_that("plot_cpmm returns a ggplot object for multiple cohorts", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("lmerTest")

  set.seed(1)

  make_cohort <- function(
    n_id = 8,
    n_visit = 5,
    y_min = -6,
    y_max = 6,
    cohort_shift = 0
  ) {
    n <- n_id * n_visit

    df <- data.frame(
      SUBID = rep(paste0("S", seq_len(n_id)), each = n_visit),
      years_since_onset = runif(n, y_min, y_max),
      SEX = rbinom(n, 1, 0.5),
      BASELINE_AGE = rnorm(n, 68, 6)
    )

    before <- pmax(0, -df$years_since_onset)
    after <- pmax(0, df$years_since_onset)

    df$P1 <- 1 +
      cohort_shift +
      0.04 * before -
      0.06 * after +
      0.15 * df$SEX +
      0.01 * df$BASELINE_AGE +
      rnorm(n, 0, 0.05)

    df
  }

  cohort_dfs <- list(
    Cohort_A = make_cohort(cohort_shift = 0.00),
    Cohort_B = make_cohort(cohort_shift = 0.15),
    Cohort_C = make_cohort(cohort_shift = -0.10)
  )

  p <- plot_cpmm(
    cohort_dfs = cohort_dfs,
    protein = "P1",
    covariates = c("SEX", "BASELINE_AGE"),
    subject_id_col = "SUBID",
    years_since_onset_col = "years_since_onset",
    show_spaghetti = TRUE
  )

  expect_s3_class(p, "ggplot")
})


test_that("plot_cpmm works with one cohort", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("lmerTest")

  set.seed(2)

  n_id <- 8
  n_visit <- 5
  n <- n_id * n_visit

  df <- data.frame(
    SUBID = rep(paste0("S", seq_len(n_id)), each = n_visit),
    years_since_onset = runif(n, -5, 5),
    SEX = rbinom(n, 1, 0.5),
    BASELINE_AGE = rnorm(n, 67, 5)
  )

  before <- pmax(0, -df$years_since_onset)
  after <- pmax(0, df$years_since_onset)

  df$P1 <- 0.8 +
    0.03 * before -
    0.04 * after +
    0.12 * df$SEX +
    0.01 * df$BASELINE_AGE +
    rnorm(n, 0, 0.05)

  p <- plot_cpmm(
    cohort_dfs = list(Single_Cohort = df),
    protein = "P1",
    show_spaghetti = FALSE
  )

  expect_s3_class(p, "ggplot")
})


test_that("plot_cpmm can use custom cohort colors", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("lmerTest")

  set.seed(3)

  make_cohort <- function(n_id = 8, n_visit = 5) {
    n <- n_id * n_visit

    df <- data.frame(
      SUBID = rep(paste0("S", seq_len(n_id)), each = n_visit),
      years_since_onset = runif(n, -5, 5),
      SEX = rbinom(n, 1, 0.5),
      BASELINE_AGE = rnorm(n, 67, 5)
    )

    before <- pmax(0, -df$years_since_onset)
    after <- pmax(0, df$years_since_onset)

    df$P1 <- 1 +
      0.03 * before -
      0.04 * after +
      0.10 * df$SEX +
      0.01 * df$BASELINE_AGE +
      rnorm(n, 0, 0.05)

    df
  }

  cohort_dfs <- list(
    Discovery = make_cohort(),
    Replication = make_cohort()
  )

  p <- plot_cpmm(
    cohort_dfs = cohort_dfs,
    protein = "P1",
    cohort_colors = c(
      Discovery = "#00565c",
      Replication = "#ba9629"
    )
  )

  expect_s3_class(p, "ggplot")
})


test_that("plot_cpmm removes outliers for plotting only", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("lmerTest")

  set.seed(4)

  make_cohort <- function(n_id = 8, n_visit = 5) {
    n <- n_id * n_visit

    df <- data.frame(
      SUBID = rep(paste0("S", seq_len(n_id)), each = n_visit),
      years_since_onset = runif(n, -5, 5),
      SEX = rbinom(n, 1, 0.5),
      BASELINE_AGE = rnorm(n, 67, 5)
    )

    before <- pmax(0, -df$years_since_onset)
    after <- pmax(0, df$years_since_onset)

    df$P1 <- 1 +
      0.03 * before -
      0.04 * after +
      0.10 * df$SEX +
      0.01 * df$BASELINE_AGE +
      rnorm(n, 0, 0.05)

    df$P1[1] <- df$P1[1] + 10

    df
  }

  cohort_dfs <- list(
    Cohort_A = make_cohort(),
    Cohort_B = make_cohort()
  )

  p <- plot_cpmm(
    cohort_dfs = cohort_dfs,
    protein = "P1",
    remove_outliers_for_plot = TRUE,
    show_spaghetti = TRUE
  )

  expect_s3_class(p, "ggplot")
})


test_that("plot_cpmm errors when no valid cohort data are available", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("lmerTest")

  bad_df <- data.frame(
    SUBID = paste0("S", 1:5),
    years_since_onset = rnorm(5),
    SEX = rbinom(5, 1, 0.5),
    BASELINE_AGE = rnorm(5, 67, 5)
  )

  expect_error(
    suppressWarnings(
      plot_cpmm(
        cohort_dfs = list(Cohort_A = bad_df),
        protein = "P1"
      )
    ),
    "No valid cohort data"
  )
})
