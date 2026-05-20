test_that("make_expression_survival_df creates expected subject-level survival data", {
  cu <- tibble::tibble(
    SUBID = c("CU1", "CU1", "CU2", "CU2"),
    PROCEDURE_AGE = c(60, 64, 62, 67),
    ONSET_AGE = c(NA, NA, NA, NA),
    P1 = c(1.0, 1.2, 2.0, 2.1)
  )

  ci <- tibble::tibble(
    SUBID = c("CI1", "CI1"),
    PROCEDURE_AGE = c(70, 72),
    ONSET_AGE = c(70, 70),
    P1 = c(3.0, 3.2)
  )

  conv <- tibble::tibble(
    SUBID = c("C1", "C1", "C2", "C2"),
    PROCEDURE_AGE = c(65, 68, 66, 69),
    ONSET_AGE = c(68, 68, 69, 69),
    P1 = c(4.0, 4.3, 5.0, 5.3)
  )

  surv_df <- make_expression_survival_df(
    cu = cu,
    ci = ci,
    conv = conv,
    protein = "P1",
    include_ci = TRUE
  )

  expect_s3_class(surv_df, "tbl_df")
  expect_equal(nrow(surv_df), 5L)

  expected_cols <- c(
    "SUBID",
    "source_group",
    "baseline_age",
    "end_age",
    "onset_age",
    "baseline_expr",
    "event",
    "expr_group",
    "cutoff",
    "Protein",
    "age_time"
  )

  expect_true(all(expected_cols %in% names(surv_df)))

  expect_equal(
    surv_df$event[surv_df$SUBID == "CU1"],
    0L
  )

  expect_equal(
    surv_df$event[surv_df$SUBID == "CI1"],
    1L
  )

  expect_equal(
    surv_df$event[surv_df$SUBID == "C1"],
    1L
  )

  expect_equal(
    surv_df$end_age[surv_df$SUBID == "CU1"],
    64
  )

  expect_equal(
    surv_df$end_age[surv_df$SUBID == "CI1"],
    70
  )

  expect_equal(
    surv_df$end_age[surv_df$SUBID == "C1"],
    68
  )

  expect_equal(
    surv_df$age_time[surv_df$SUBID == "C1"],
    68
  )

  expect_true(all(surv_df$expr_group %in% c("High", "Low")))
  expect_true(all(surv_df$Protein == "P1"))
})


test_that("make_expression_survival_df can exclude CI subjects", {
  cu <- tibble::tibble(
    SUBID = c("CU1", "CU1"),
    PROCEDURE_AGE = c(60, 65),
    P1 = c(1.0, 1.2)
  )

  ci <- tibble::tibble(
    SUBID = c("CI1", "CI1"),
    PROCEDURE_AGE = c(70, 72),
    P1 = c(3.0, 3.2)
  )

  conv <- tibble::tibble(
    SUBID = c("C1", "C1"),
    PROCEDURE_AGE = c(65, 68),
    ONSET_AGE = c(68, 68),
    P1 = c(4.0, 4.3)
  )

  surv_df <- make_expression_survival_df(
    cu = cu,
    ci = ci,
    conv = conv,
    protein = "P1",
    include_ci = FALSE
  )

  expect_equal(nrow(surv_df), 2L)
  expect_false("CI" %in% surv_df$source_group)
})


test_that("make_expression_survival_df supports CI event time using onset column", {
  cu <- tibble::tibble(
    SUBID = c("CU1", "CU1"),
    PROCEDURE_AGE = c(60, 65),
    ONSET_AGE_minus = c(NA, NA),
    P1 = c(1.0, 1.2)
  )

  ci <- tibble::tibble(
    SUBID = c("CI1", "CI1"),
    PROCEDURE_AGE = c(70, 72),
    ONSET_AGE_minus = c(69, 69),
    P1 = c(3.0, 3.2)
  )

  conv <- tibble::tibble(
    SUBID = c("C1", "C1"),
    PROCEDURE_AGE = c(65, 68),
    ONSET_AGE_minus = c(67, 67),
    P1 = c(4.0, 4.3)
  )

  surv_df <- make_expression_survival_df(
    cu = cu,
    ci = ci,
    conv = conv,
    protein = "P1",
    onset_col = "ONSET_AGE_minus",
    include_ci = TRUE,
    ci_event_time = "onset"
  )

  expect_equal(
    surv_df$end_age[surv_df$SUBID == "CI1"],
    69
  )

  expect_equal(
    surv_df$onset_age[surv_df$SUBID == "CI1"],
    69
  )

  expect_equal(
    surv_df$end_age[surv_df$SUBID == "C1"],
    67
  )
})


test_that("make_expression_survival_df supports expression grouping methods", {
  cu <- tibble::tibble(
    SUBID = c("CU1", "CU2"),
    PROCEDURE_AGE = c(60, 62),
    P1 = c(1, 2)
  )

  ci <- tibble::tibble(
    SUBID = c("CI1"),
    PROCEDURE_AGE = c(70),
    P1 = c(3)
  )

  conv <- tibble::tibble(
    SUBID = c("C1", "C2"),
    PROCEDURE_AGE = c(65, 66),
    ONSET_AGE = c(68, 69),
    P1 = c(4, 5)
  )

  median_df <- make_expression_survival_df(
    cu = cu,
    ci = ci,
    conv = conv,
    protein = "P1",
    expr_group_method = "median"
  )

  upper_df <- make_expression_survival_df(
    cu = cu,
    ci = ci,
    conv = conv,
    protein = "P1",
    expr_group_method = "upper25"
  )

  lower_df <- make_expression_survival_df(
    cu = cu,
    ci = ci,
    conv = conv,
    protein = "P1",
    expr_group_method = "lower25"
  )

  expect_equal(unique(median_df$cutoff), 3)
  expect_equal(unique(upper_df$cutoff), 4)
  expect_equal(unique(lower_df$cutoff), 2)

  expect_true(all(median_df$expr_group %in% c("High", "Low")))
  expect_true(all(upper_df$expr_group %in% c("High", "Low")))
  expect_true(all(lower_df$expr_group %in% c("High", "Low")))
})


test_that("make_expression_survival_df errors for missing required columns", {
  cu <- tibble::tibble(
    SUBID = c("CU1"),
    PROCEDURE_AGE = c(60)
  )

  ci <- tibble::tibble(
    SUBID = c("CI1"),
    PROCEDURE_AGE = c(70),
    P1 = c(2)
  )

  conv <- tibble::tibble(
    SUBID = c("C1"),
    PROCEDURE_AGE = c(65),
    ONSET_AGE = c(68),
    P1 = c(3)
  )

  expect_error(
    make_expression_survival_df(
      cu = cu,
      ci = ci,
      conv = conv,
      protein = "P1"
    ),
    "cu is missing required columns"
  )
})


test_that("make_expression_survival_df errors for unsupported expression grouping method", {
  cu <- tibble::tibble(
    SUBID = c("CU1"),
    PROCEDURE_AGE = c(60),
    P1 = c(1)
  )

  ci <- tibble::tibble(
    SUBID = c("CI1"),
    PROCEDURE_AGE = c(70),
    P1 = c(2)
  )

  conv <- tibble::tibble(
    SUBID = c("C1"),
    PROCEDURE_AGE = c(65),
    ONSET_AGE = c(68),
    P1 = c(3)
  )

  expect_error(
    make_expression_survival_df(
      cu = cu,
      ci = ci,
      conv = conv,
      protein = "P1",
      expr_group_method = "bad_method"
    ),
    "expr_group_method must be"
  )
})
