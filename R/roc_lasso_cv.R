#' Cross-validated ROC analysis with LASSO feature selection
#'
#' Performs LASSO (L1) logistic regression for feature selection, followed by
#' stratified K-fold cross-validation using a standard (L2) logistic regression
#' model. Covariates are optional and are always retained in the final model.
#'
#' @param df Data frame containing predictors and a binary outcome.
#' @param protein_features Character vector of protein / feature column names.
#' @param label_col Character; binary outcome column (0/1). Default "label".
#' @param covariates Optional character vector of covariate column names.
#' @param n_folds Integer; number of CV folds. Default 5.
#' @param seed Integer random seed. Default 42.
#' @param title Plot title.
#' @param export Logical; save plot?
#' @param export_path File path if saving.
#'
#' @return Named list with selected features, AUC stats, and ggplot.
#' @export
plot_lasso_cv_roc <- function(
    df,
    protein_features,
    label_col = "label",
    covariates = NULL,
    n_folds = 5,
    seed = 42,
    title = "ROC Curve",
    export = FALSE,
    export_path = "roc_plot.svg"
) {
  requireNamespace("glmnet", quietly = TRUE)
  requireNamespace("pROC", quietly = TRUE)
  requireNamespace("ggplot2", quietly = TRUE)

  set.seed(seed)
  if (is.null(covariates)) covariates <- character(0)

  features <- unique(c(protein_features, covariates))
  df <- df[, c(features, label_col), drop = FALSE]
  df <- df[stats::complete.cases(df), , drop = FALSE]

  y <- df[[label_col]]
  if (length(unique(y)) < 2) {
    stop("Outcome must contain two classes.", call. = FALSE)
  }

  X_mat <- stats::model.matrix(~ . -1, data = df[, features, drop = FALSE])

  # ---- LASSO ----
  cv_lasso <- glmnet::cv.glmnet(
    X_mat, y, family = "binomial", alpha = 1, nfolds = n_folds
  )

  coef_mat <- as.matrix(glmnet::coef.glmnet(cv_lasso, s = "lambda.1se"))
  selected_features <- rownames(coef_mat)[coef_mat[,1] != 0]
  selected_features <- setdiff(selected_features, "(Intercept)")
  final_features <- unique(c(selected_features, covariates))

  # ---- Stratified folds ----
  fold_id <- split(
    seq_len(nrow(df)),
    interaction(y, sample(seq_len(n_folds), nrow(df), replace = TRUE))
  )

  mean_fpr <- seq(0, 1, length.out = 100)
  tprs <- list()
  aucs <- numeric()

  for (idx in fold_id) {
    train_idx <- setdiff(seq_len(nrow(df)), idx)
    test_idx <- idx

    y_test <- y[test_idx]
    if (length(unique(y_test)) < 2) next

    train_df <- df[train_idx, , drop = FALSE]
    test_df  <- df[test_idx, , drop = FALSE]

    X_train <- stats::model.matrix(
      stats::reformulate(final_features),
      data = train_df
    )
    X_test <- stats::model.matrix(
      stats::reformulate(final_features),
      data = test_df
    )

    mu <- colMeans(X_train)
    sdv <- apply(X_train, 2, stats::sd)
    sdv[sdv == 0] <- 1

    X_train <- scale(X_train, mu, sdv)
    X_test  <- scale(X_test,  mu, sdv)

    fit <- stats::glm(
      y[train_idx] ~ .,
      data = data.frame(y = y[train_idx], X_train),
      family = stats::binomial()
    )

    prob <- as.numeric(stats::predict(
      fit, newdata = data.frame(X_test), type = "response"
    ))

    roc_obj <- pROC::roc(y_test, prob, quiet = TRUE)
    aucs <- c(aucs, as.numeric(pROC::auc(roc_obj)))

    interp <- stats::approx(
      x = 1 - roc_obj$specificities,
      y = roc_obj$sensitivities,
      xout = mean_fpr
    )$y
    interp[1] <- 0
    tprs[[length(tprs) + 1]] <- interp
  }

  # ---- Handle test edge case ----
  if (length(aucs) == 0) {
    warning("All CV folds had a single outcome class; ROC not computed.")
    return(list(
      selected_features = selected_features,
      final_features = final_features,
      mean_auc = NA_real_,
      sd_auc = NA_real_,
      fold_aucs = numeric(),
      plot = ggplot2::ggplot() + ggplot2::theme_void()
    ))
  }

  mean_tpr <- Reduce("+", tprs) / length(tprs)
  mean_tpr[length(mean_tpr)] <- 1

  mean_auc <- mean(aucs)
  sd_auc <- stats::sd(aucs)

  plot_df <- data.frame(FPR = mean_fpr, TPR = mean_tpr)

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(FPR, TPR)) +
    ggplot2::geom_line(color = "#80796B", linewidth = 1.6) +
    ggplot2::geom_abline(linetype = "dashed") +
    ggplot2::labs(
      title = title,
      subtitle = sprintf("Mean AUC = %.3f ± %.3f", mean_auc, sd_auc),
      x = "1 – Specificity",
      y = "Sensitivity"
    ) +
    ggplot2::theme_minimal()

  if (export) {
    ggplot2::ggsave(export_path, p, width = 7, height = 6, dpi = 300)
  }

  list(
    selected_features = selected_features,
    final_features = final_features,
    mean_auc = mean_auc,
    sd_auc = sd_auc,
    fold_aucs = aucs,
    plot = p
  )
}
