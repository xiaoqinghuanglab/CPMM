test_that("run_cpmm_meta_analysis returns expected meta-analysis columns", {
  skip_on_cran()
  skip_if_not_installed("metafor")

  all_cohort_results_df <- data.frame(
    Protein = rep(c("P1", "P2", "P3"), each = 3),
    Cohort = rep(c("IADRC", "GNPC_B", "GNPC_P"), times = 3),
    Delta = c(
      0.20, 0.18, 0.22,
      -0.10, -0.12, -0.08,
      0.05, NA, 0.06
    ),
    `SE Delta` = c(
      0.04, 0.05, 0.04,
      0.03, 0.04, 0.03,
      0.02, 0.03, 0.02
    ),
    check.names = FALSE
  )

  res <- run_cpmm_meta_analysis(
    all_cohort_results_df = all_cohort_results_df,
    protein_col = "Protein",
    delta_col = "Delta",
    se_delta_col = "SE Delta",
    min_cohorts = 3
  )

  expect_s3_class(res, "tbl_df")

  expected_cols <- c(
    "Protein",
    "N Cohorts",
    "Meta Delta",
    "Meta SE",
    "Meta 95% CI Lower",
    "Meta 95% CI Upper",
    "Meta Z",
    "Meta P-value",
    "Q",
    "df",
    "Tau^2",
    "I^2",
    "Meta Adjusted P-value (FDR)",
    "Meta Significant"
  )

  expect_true(all(expected_cols %in% names(res)))

  expect_false("Same Direction All" %in% names(res))
  expect_false("N Positive" %in% names(res))
  expect_false("N Negative" %in% names(res))

  expect_true(all(c("P1", "P2") %in% res$Protein))
  expect_false("P3" %in% res$Protein)

  expect_true(all(res$`N Cohorts` == 3))
  expect_true(all(is.finite(res$`Meta Delta`)))
  expect_true(all(is.finite(res$`Meta SE`)))
  expect_true(all(is.finite(res$`Meta P-value`)))
  expect_true(all(is.finite(res$`Meta Adjusted P-value (FDR)`)))

  expect_true(all(res$`Meta 95% CI Lower` <= res$`Meta Delta`))
  expect_true(all(res$`Meta 95% CI Upper` >= res$`Meta Delta`))

  expect_type(res$`Meta Significant`, "logical")
})


test_that("run_cpmm_meta_analysis excludes invalid SE Delta values", {
  skip_on_cran()
  skip_if_not_installed("metafor")

  all_cohort_results_df <- data.frame(
    Protein = c("P1", "P1", "P1", "P2", "P2", "P2"),
    Cohort = c("IADRC", "GNPC_B", "GNPC_P", "IADRC", "GNPC_B", "GNPC_P"),
    Delta = c(0.20, 0.18, 0.22, 0.10, 0.11, 0.12),
    `SE Delta` = c(0.04, 0.05, 0.04, 0.03, 0, NA),
    check.names = FALSE
  )

  res <- run_cpmm_meta_analysis(
    all_cohort_results_df = all_cohort_results_df,
    min_cohorts = 3
  )

  expect_true("P1" %in% res$Protein)
  expect_false("P2" %in% res$Protein)
})


test_that("run_cpmm_meta_analysis stops when required columns are missing", {
  skip_on_cran()

  bad_df <- data.frame(
    Protein = c("P1", "P1", "P1"),
    Delta = c(0.1, 0.2, 0.3)
  )

  expect_error(
    run_cpmm_meta_analysis(bad_df),
    "required columns"
  )
})


test_that("run_cpmm_meta_analysis returns empty tibble when no protein has enough cohorts", {
  skip_on_cran()
  skip_if_not_installed("metafor")

  all_cohort_results_df <- data.frame(
    Protein = c("P1", "P1", "P2", "P2"),
    Cohort = c("IADRC", "GNPC_B", "IADRC", "GNPC_B"),
    Delta = c(0.1, 0.2, -0.1, -0.2),
    `SE Delta` = c(0.04, 0.05, 0.03, 0.04),
    check.names = FALSE
  )

  res <- run_cpmm_meta_analysis(
    all_cohort_results_df = all_cohort_results_df,
    min_cohorts = 3
  )

  expect_s3_class(res, "tbl_df")
  expect_equal(nrow(res), 0L)

  expected_cols <- c(
    "Protein",
    "N Cohorts",
    "Meta Delta",
    "Meta SE",
    "Meta 95% CI Lower",
    "Meta 95% CI Upper",
    "Meta Z",
    "Meta P-value",
    "Q",
    "df",
    "Tau^2",
    "I^2",
    "Meta Adjusted P-value (FDR)",
    "Meta Significant"
  )

  expect_true(all(expected_cols %in% names(res)))
})
