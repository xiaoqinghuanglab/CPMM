test_that("plot_lasso_cv_roc runs and returns expected structure", {

  skip_if_not_installed("glmnet")
  skip_if_not_installed("pROC")
  skip_on_cran()
  set.seed(1)

  # ---- Synthetic dataset ----
  n <- 60
  df <- data.frame(
    P1 = rnorm(n),
    P2 = rnorm(n),
    P3 = rnorm(n),
    PROCEDURE_AGE = rnorm(n, 65, 6),
    Sex_bin = sample(c(0, 1), n, replace = TRUE),
    label = rep(c(0, 1), length.out = n)
  )

  # ---- Run function (allow warning for small folds) ----
  res <- suppressWarnings(
    plot_lasso_cv_roc(
      df = df,
      protein_features = c("P1", "P2", "P3"),
      covariates = c("PROCEDURE_AGE", "Sex_bin"),
      label_col = "label",
      n_folds = 5,
      title = "Test ROC",
      export = FALSE
    )
  )

  # ---- Structural checks ----
  expect_type(res, "list")

  expect_true(all(c(
    "selected_features",
    "final_features",
    "mean_auc",
    "sd_auc",
    "fold_aucs",
    "plot"
  ) %in% names(res)))

  # ---- Feature outputs ----
  expect_type(res$selected_features, "character")
  expect_type(res$final_features, "character")

  # ---- AUC outputs (may be NA in edge cases) ----
  expect_true(is.numeric(res$mean_auc) || is.na(res$mean_auc))
  expect_true(is.numeric(res$sd_auc)   || is.na(res$sd_auc))
  expect_type(res$fold_aucs, "double")

  # ---- Plot ----
  expect_s3_class(res$plot, "ggplot")
})
