test_that("plot_quadrant_beta returns a ggplot", {
  df <- tibble::tibble(
    Protein = paste0("P", 1:8),
    `Beta 1` = c(-0.2, 0.1, 0.3, -0.4, 0.2, -0.1, 0.4, -0.3),
    `Beta 2` = c( 0.3, 0.2,-0.1, -0.5, 0.1,  0.2,-0.2,  0.4),
    `Adjusted P-value` = c(0.01, 0.2, 0.03, 0.5, 0.2, 0.01, 0.04, 0.06)
  )

  p <- plot_quadrant_beta(
    wald_df = df,
    beta_x_col = "Beta 1",
    beta_y_col = "Beta 2",
    fdr_col = "Adjusted P-value",
    export = FALSE
  )

  expect_s3_class(p, "ggplot")
})
