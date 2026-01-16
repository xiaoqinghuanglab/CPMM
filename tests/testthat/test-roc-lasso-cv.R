test_that("plot_lasso_cv_roc runs and returns expected structure", {

  skip_if_not_installed("glmnet")
  skip_if_not_installed("pROC")
  skip_on_cran()

  set.seed(123)

  # ---- Synthetic dataset (balanced, CV-safe) ----
  n_per_class <- 50
  n <- n_per_class * 2

  df <- data.frame(
    P1 = c(rnorm(n_per_class,  0.5), rnorm(n_per_class, -0.5)),
    P2 = c(rnorm(n_per_class,  0.3), rnorm(n_per_class, -0.3)),
    P3 = rnorm(n),
    PROCEDURE_AGE = rnorm(n, 65, 6),
    Sex_bin = sample(c(0, 1), n, replace = TRUE),
    label = rep(c(0, 1), each = n_per_class)
  )

  # ---- Run function ----
  res <- plot_lasso_cv_roc(
    df = df,
    protein_features = c("P1", "P2", "P3"),
    covariates = c("PROCEDURE_AGE", "Sex_bin"),
    label_col = "label",
    n_folds = 5,
    title = "Test ROC",
    export = FALSE
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

  # ---- AUC outputs ----
  expect_type(res$mean_auc, "double")
  expect_type(res$sd_auc, "double")
  expect_true(is.finite(res$mean_auc))
  expect_true(res$mean_auc > 0.5)
  expect_length(res$fold_aucs, 5)

  # ---- Plot ----
  expect_s3_class(res$plot, "ggplot")
})
