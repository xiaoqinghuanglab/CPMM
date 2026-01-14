test_that("plot_wald_volcano returns a ggplot", {
  df <- tibble::tibble(
    Protein = c("A","B","C","D"),
    `Beta 1` = c( 0.20,  0.10, -0.10,  0.05),
    `Beta 2` = c(-0.20,  0.30, -0.05,  0.20),
    `P-value` = c(1e-5, 0.02, 0.5, 0.001),
    `Adjusted P-value` = c(2e-5, 0.04, 0.6, 0.005)
  )

  p <- plot_wald_volcano(
    wald_df = df,
    beta_before_col = "Beta 1",
    beta_after_col  = "Beta 2",
    pval_col = "P-value",
    fdr_col  = "Adjusted P-value",
    protein_col = "Protein",
    export = FALSE
  )

  expect_s3_class(p, "ggplot")
})
