test_that("plot_expression_km_onset_overlay returns a ggplot object", {
  skip_on_cran()
  skip_if_not_installed("survival")
  skip_if_not_installed("ggplot2")

  surv_df_bg <- tibble::tibble(
    SUBID = paste0("S", 1:10),
    event = c(0, 0, 0, 0, 1, 1, 1, 1, 1, 1),
    expr_group = c(
      "Low", "Low", "Low", "Low", "High",
      "High", "High", "High", "High", "High"
    ),
    age_time = c(68, 69, 70, 71, 70, 71, 67, 68, 69, 70)
  )

  surv_df_main <- tibble::tibble(
    SUBID = paste0("S", 1:10),
    event = c(0, 0, 0, 0, 1, 1, 1, 1, 1, 1),
    expr_group = c(
      "Low", "Low", "Low", "Low", "High",
      "High", "High", "High", "High", "High"
    ),
    age_time = c(68, 69, 70, 71, 69, 70, 66, 67, 68, 69)
  )

  p <- plot_expression_km_onset_overlay(
    surv_df_main = surv_df_main,
    surv_df_bg = surv_df_bg,
    protein_name = "P1",
    main_label = "Early onset",
    bg_label = "Onset",
    cohort_label = "IADRC"
  )

  expect_s3_class(p, "ggplot")
})


test_that("plot_expression_km_onset_overlay creates four legend labels", {
  skip_on_cran()
  skip_if_not_installed("survival")
  skip_if_not_installed("ggplot2")

  surv_df_bg <- tibble::tibble(
    event = c(0, 0, 1, 1),
    expr_group = c("Low", "Low", "High", "High"),
    age_time = c(68, 69, 70, 71)
  )

  surv_df_main <- tibble::tibble(
    event = c(0, 0, 1, 1),
    expr_group = c("Low", "Low", "High", "High"),
    age_time = c(67, 68, 69, 70)
  )

  p <- plot_expression_km_onset_overlay(
    surv_df_main = surv_df_main,
    surv_df_bg = surv_df_bg,
    protein_name = "P1",
    main_label = "Early onset",
    bg_label = "Onset"
  )

  color_scale <- p$scales$get_scales("colour")
  legend_breaks <- color_scale$breaks

  expect_true(all(c(
    "Low (Onset)",
    "High (Onset)",
    "Low (Early onset)",
    "High (Early onset)"
  ) %in% legend_breaks))
})


test_that("plot_expression_km_onset_overlay errors when main required columns are missing", {
  skip_on_cran()

  surv_df_main <- tibble::tibble(
    event = c(0, 1),
    expr_group = c("Low", "High")
  )

  surv_df_bg <- tibble::tibble(
    event = c(0, 1),
    expr_group = c("Low", "High"),
    age_time = c(68, 70)
  )

  expect_error(
    plot_expression_km_onset_overlay(
      surv_df_main = surv_df_main,
      surv_df_bg = surv_df_bg,
      protein_name = "P1"
    ),
    "surv_df_main is missing required columns"
  )
})


test_that("plot_expression_km_onset_overlay errors when background required columns are missing", {
  skip_on_cran()

  surv_df_main <- tibble::tibble(
    event = c(0, 1),
    expr_group = c("Low", "High"),
    age_time = c(67, 69)
  )

  surv_df_bg <- tibble::tibble(
    event = c(0, 1),
    expr_group = c("Low", "High")
  )

  expect_error(
    plot_expression_km_onset_overlay(
      surv_df_main = surv_df_main,
      surv_df_bg = surv_df_bg,
      protein_name = "P1"
    ),
    "surv_df_bg is missing required columns"
  )
})


test_that("plot_expression_km_onset_overlay requires both Low and High groups in main data", {
  skip_on_cran()
  skip_if_not_installed("survival")
  skip_if_not_installed("ggplot2")

  surv_df_main <- tibble::tibble(
    event = c(0, 0, 1, 1),
    expr_group = c("High", "High", "High", "High"),
    age_time = c(67, 68, 69, 70)
  )

  surv_df_bg <- tibble::tibble(
    event = c(0, 0, 1, 1),
    expr_group = c("Low", "Low", "High", "High"),
    age_time = c(68, 69, 70, 71)
  )

  expect_error(
    plot_expression_km_onset_overlay(
      surv_df_main = surv_df_main,
      surv_df_bg = surv_df_bg,
      protein_name = "P1"
    ),
    "Both Low and High groups are required"
  )
})


test_that("plot_expression_km_onset_overlay can save an SVG file", {
  skip_on_cran()
  skip_if_not_installed("survival")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("svglite")

  surv_df_bg <- tibble::tibble(
    event = c(0, 0, 0, 0, 1, 1, 1, 1),
    expr_group = c("Low", "Low", "Low", "Low", "High", "High", "High", "High"),
    age_time = c(68, 69, 70, 71, 70, 71, 67, 68)
  )

  surv_df_main <- tibble::tibble(
    event = c(0, 0, 0, 0, 1, 1, 1, 1),
    expr_group = c("Low", "Low", "Low", "Low", "High", "High", "High", "High"),
    age_time = c(68, 69, 70, 71, 69, 70, 66, 67)
  )

  tmp_dir <- tempdir()

  p <- plot_expression_km_onset_overlay(
    surv_df_main = surv_df_main,
    surv_df_bg = surv_df_bg,
    protein_name = "P1",
    save = TRUE,
    save_path = tmp_dir
  )

  expect_s3_class(p, "ggplot")

  expect_true(
    file.exists(file.path(tmp_dir, "KM_overlay_P1_age_time.svg"))
  )
})
