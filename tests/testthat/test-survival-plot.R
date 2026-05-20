test_that("plot_expression_km returns a ggplot object", {
  skip_on_cran()
  skip_if_not_installed("survival")
  skip_if_not_installed("ggplot2")

  surv_df <- tibble::tibble(
    SUBID = paste0("S", 1:10),
    source_group = c(rep("CU", 4), rep("CI", 2), rep("Converters", 4)),
    baseline_age = c(60, 61, 62, 63, 70, 71, 64, 65, 66, 67),
    end_age = c(68, 69, 70, 71, 70, 71, 67, 68, 69, 70),
    onset_age = c(NA, NA, NA, NA, 70, 71, 67, 68, 69, 70),
    baseline_expr = c(1, 1.1, 1.2, 1.3, 3, 3.1, 4, 4.1, 4.2, 4.3),
    event = c(0, 0, 0, 0, 1, 1, 1, 1, 1, 1),
    expr_group = c("Low", "Low", "Low", "Low", "High", "High", "High", "High", "High", "High"),
    cutoff = 2.15,
    Protein = "P1",
    age_time = c(68, 69, 70, 71, 70, 71, 67, 68, 69, 70)
  )

  p <- plot_expression_km(
    surv_df = surv_df,
    protein_name = "P1"
  )

  expect_s3_class(p, "ggplot")
})


test_that("plot_expression_km errors when required columns are missing", {
  skip_on_cran()

  bad_df <- tibble::tibble(
    SUBID = paste0("S", 1:4),
    event = c(0, 0, 1, 1),
    expr_group = c("Low", "Low", "High", "High")
  )

  expect_error(
    plot_expression_km(
      surv_df = bad_df,
      protein_name = "P1"
    ),
    "missing required columns"
  )
})


test_that("plot_expression_km requires both Low and High groups", {
  skip_on_cran()
  skip_if_not_installed("survival")
  skip_if_not_installed("ggplot2")

  surv_df <- tibble::tibble(
    SUBID = paste0("S", 1:5),
    event = c(0, 0, 1, 1, 1),
    expr_group = rep("High", 5),
    age_time = c(65, 66, 67, 68, 69)
  )

  expect_error(
    plot_expression_km(
      surv_df = surv_df,
      protein_name = "P1"
    ),
    "Both Low and High groups are required"
  )
})


test_that("plot_expression_km can save an SVG file", {
  skip_on_cran()
  skip_if_not_installed("survival")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("svglite")

  surv_df <- tibble::tibble(
    SUBID = paste0("S", 1:10),
    event = c(0, 0, 0, 0, 1, 1, 1, 1, 1, 1),
    expr_group = c("Low", "Low", "Low", "Low", "High", "High", "High", "High", "High", "High"),
    age_time = c(68, 69, 70, 71, 70, 71, 67, 68, 69, 70)
  )

  tmp_dir <- tempdir()

  p <- plot_expression_km(
    surv_df = surv_df,
    protein_name = "P1",
    save = TRUE,
    save_path = tmp_dir
  )

  expect_s3_class(p, "ggplot")

  expect_true(
    file.exists(file.path(tmp_dir, "KM_P1.svg"))
  )
})
