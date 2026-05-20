test_that("fit_cpmm_all_proteins returns updated CPMM output columns", {
  skip_on_cran()
  set.seed(1)

  make_df <- function(ids, yso_min = -5, yso_max = 5) {
    df <- data.frame(
      SUBID = rep(ids, each = 5),
      years_since_onset = runif(length(ids) * 5, yso_min, yso_max),
      SEX = sample(c("M", "F"), length(ids) * 5, replace = TRUE),
      BASELINE_AGE = rnorm(length(ids) * 5, 65, 6)
    )

    before <- pmax(0, -df$years_since_onset)
    after  <- pmax(0,  df$years_since_onset)

    signal <- 0.5 +
      0.03 * before -
      0.05 * after +
      0.2  * (df$SEX == "M") +
      0.01 * df$BASELINE_AGE

    df$P1 <- signal + rnorm(nrow(df), 0, 0.05)
    df$P2 <- signal * 0.5 + rnorm(nrow(df), 0, 0.05)

    df
  }

  ids_status <- paste0("S", 1:10)
  df_status <- make_df(ids_status)

  res <- fit_cpmm_all_proteins(
    df_status_change = df_status,
    protein_list = c("P1", "P2"),
    covariates = c("SEX", "BASELINE_AGE"),
    subject_id_col = "SUBID",
    years_since_onset_col = "years_since_onset"
  )

  expect_s3_class(res, "tbl_df")
  expect_equal(nrow(res), 2L)

  expected_cols <- c(
    "Protein",

    "Beta 1",
    "SE Beta 1",
    "Beta 2",
    "SE Beta 2",

    "Beta Before Raw",
    "Beta After Raw",
    "Var Before Raw",
    "Var After Raw",
    "Cov Before After",

    "Delta",
    "SE Delta",
    "Delta 95% CI Lower",
    "Delta 95% CI Upper",

    "Intercept",
    "AIC",
    "BIC",
    "MSE",
    "N Obs"
  )

  expect_true(all(expected_cols %in% names(res)))

  expect_false("Wald Statistic" %in% names(res))
  expect_false("P-value" %in% names(res))

  numeric_cols <- c(
    "Beta 1",
    "SE Beta 1",
    "Beta 2",
    "SE Beta 2",
    "Beta Before Raw",
    "Beta After Raw",
    "Var Before Raw",
    "Var After Raw",
    "Cov Before After",
    "Delta",
    "SE Delta",
    "Delta 95% CI Lower",
    "Delta 95% CI Upper",
    "Intercept",
    "AIC",
    "BIC",
    "MSE",
    "N Obs"
  )

  expect_true(all(vapply(res[numeric_cols], is.numeric, logical(1))))

  expect_true(all(is.finite(res$`Beta 1`)))
  expect_true(all(is.finite(res$`Beta 2`)))
  expect_true(all(is.finite(res$Delta)))
  expect_true(all(is.finite(res$`SE Delta`)))
  expect_true(all(is.finite(res$AIC)))
  expect_true(all(is.finite(res$BIC)))
  expect_true(all(is.finite(res$MSE)))

  expect_equal(
    res$Delta,
    res$`Beta 2` - res$`Beta 1`,
    tolerance = 1e-8
  )

  expect_equal(
    res$`Beta 1`,
    -res$`Beta Before Raw`,
    tolerance = 1e-8
  )

  expect_equal(
    res$`Beta 2`,
    res$`Beta After Raw`,
    tolerance = 1e-8
  )

  expect_true(all(res$`Delta 95% CI Lower` <= res$Delta))
  expect_true(all(res$`Delta 95% CI Upper` >= res$Delta))
})


test_that("fit_cpmm_all_proteins handles missing protein columns", {
  skip_on_cran()
  set.seed(2)

  df_status <- data.frame(
    SUBID = rep(paste0("S", 1:8), each = 5),
    years_since_onset = runif(40, -5, 5),
    SEX = sample(c("M", "F"), 40, replace = TRUE),
    BASELINE_AGE = rnorm(40, 65, 6)
  )

  before <- pmax(0, -df_status$years_since_onset)
  after  <- pmax(0,  df_status$years_since_onset)

  df_status$P1 <- 1 +
    0.04 * before -
    0.02 * after +
    0.01 * df_status$BASELINE_AGE +
    rnorm(nrow(df_status), 0, 0.05)

  res <- fit_cpmm_all_proteins(
    df_status_change = df_status,
    protein_list = c("P1", "MissingProtein"),
    covariates = c("SEX", "BASELINE_AGE")
  )

  expect_equal(nrow(res), 2L)

  missing_row <- res[res$Protein == "MissingProtein", ]

  expect_true(is.na(missing_row$`Beta 1`))
  expect_true(is.na(missing_row$`Beta 2`))
  expect_true(is.na(missing_row$Delta))
  expect_true(is.na(missing_row$`SE Delta`))
  expect_true(is.na(missing_row$AIC))
  expect_true(is.na(missing_row$BIC))
  expect_true(is.na(missing_row$MSE))
})


test_that("fit_cpmm_all_proteins stops when required columns are missing", {
  skip_on_cran()

  df_missing_subject <- data.frame(
    years_since_onset = runif(20, -5, 5),
    SEX = sample(c("M", "F"), 20, replace = TRUE),
    BASELINE_AGE = rnorm(20, 65, 6),
    P1 = rnorm(20)
  )

  expect_error(
    fit_cpmm_all_proteins(
      df_status_change = df_missing_subject,
      protein_list = "P1"
    ),
    "Subject ID column"
  )

  df_missing_yso <- data.frame(
    SUBID = rep(paste0("S", 1:5), each = 4),
    SEX = sample(c("M", "F"), 20, replace = TRUE),
    BASELINE_AGE = rnorm(20, 65, 6),
    P1 = rnorm(20)
  )

  expect_error(
    fit_cpmm_all_proteins(
      df_status_change = df_missing_yso,
      protein_list = "P1"
    ),
    "Years-since-onset column"
  )

  df_missing_covariate <- data.frame(
    SUBID = rep(paste0("S", 1:5), each = 4),
    years_since_onset = runif(20, -5, 5),
    SEX = sample(c("M", "F"), 20, replace = TRUE),
    P1 = rnorm(20)
  )

  expect_error(
    fit_cpmm_all_proteins(
      df_status_change = df_missing_covariate,
      protein_list = "P1",
      covariates = c("SEX", "BASELINE_AGE")
    ),
    "covariates were not found"
  )
})
