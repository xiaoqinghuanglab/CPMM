test_that("fit_cpmm_all_proteins returns expected columns", {
  skip_on_cran()
  set.seed(1)

  make_df <- function(ids, yso_min = -5, yso_max = 5) {
    df <- data.frame(
      SUBID = rep(ids, each = 4),
      years_since_onset = runif(length(ids) * 4, yso_min, yso_max),
      SEX = sample(c("M","F"), length(ids) * 4, replace = TRUE),
      BASELINE_AGE = rnorm(length(ids) * 4, 65, 6)
    )

    # Generative CPMM signal
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

  ids_status <- paste0("S", 1:6)
  df_status  <- make_df(ids_status)

  res <- fit_cpmm_all_proteins(
    df_status_change = df_status,
    protein_list = c("P1", "P2"),
    covariates = c("SEX", "BASELINE_AGE")
  )

  # structure checks
  expect_s3_class(res, "tbl_df")
  expect_equal(nrow(res), 2L)

  expected_cols <- c(
    "Protein",
    "Beta 1", "SE Beta 1",
    "Beta 2", "SE Beta 2",
    "Intercept",
    "AIC", "BIC", "MSE"
  )

  expect_true(all(expected_cols %in% names(res)))

  # numeric sanity
  expect_true(all(is.numeric(res$`Beta 1`)))
  expect_true(all(is.numeric(res$`Beta 2`)))
  expect_true(all(is.numeric(res$AIC)))
  expect_true(all(is.numeric(res$MSE)))
})
