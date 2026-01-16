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

  # ------------------ Input prep ------------------
  features <- unique(c(protein_features, covariates))
  df <- df[, c(features, label_col), drop = FALSE]
  df <- df[stats::complete.cases(df), , drop = FALSE]

  y <- df[[label_col]]
  if (length(unique(y)) < 2) {
    stop("Outcome must contain two classes.", call. = FALSE)
  }

  # ------------------ LASSO feature selection ------------------
  X_mat <- stats::model.matrix(~ . -1, data = df[, features, drop = FALSE])

  cv_lasso <- glmnet::cv.glmnet(
    x = X_mat,
    y = y,
    family = "binomial",
    alpha = 1,
    nfolds = n_folds
  )

  coef_mat <- as.matrix(glmnet::coef.glmnet(cv_lasso, s = "lambda.1se"))
  selected_features <- rownames(coef_mat)[coef_mat[, 1] != 0]
  selected_features <- setdiff(selected_features, "(Intercept)")
  final_features <- unique(c(selected_features, covariates))

  if (!any(selected_features %in% protein_features)) {
    warning(
      "No protein features were selected by LASSO; model uses covariates only.",
      call. = FALSE
    )
  }

  # ------------------ Stratified folds ------------------
  make_stratified_folds <- function(y, k, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)

    idx0 <- sample(which(y == 0))
    idx1 <- sample(which(y == 1))

    split0 <- split(idx0, rep(seq_len(k), length.out = length(idx0)))
    split1 <- split(idx1, rep(seq_len(k), length.out = length(idx1)))

    folds <- vector("list", k)
    for (i in seq_len(k)) {
      folds[[i]] <- c(split0[[i]], split1[[i]])
    }
    folds
  }

  folds <- make_stratified_folds(y, n_folds, seed)

  # ------------------ CV ROC ------------------
  mean_fpr <- seq(0, 1, length.out = 200)
  tprs <- list()
  aucs <- numeric()

  for (i in seq_along(folds)) {
    test_idx <- folds[[i]]
    train_idx <- setdiff(seq_len(nrow(df)), test_idx)

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

    prob <- as.numeric(
      stats::predict(fit, newdata = data.frame(X_test), type = "response")
    )

    roc_obj <- pROC::roc(y_test, prob, quiet = TRUE)
    aucs <- c(aucs, as.numeric(pROC::auc(roc_obj)))

    interp <- stats::approx(
      x = 1 - roc_obj$specificities,
      y = roc_obj$sensitivities,
      xout = mean_fpr,
      ties = "ordered",
      rule = 2
    )$y

    interp[1] <- 0
    tprs[[length(tprs) + 1]] <- interp
  }

  # ------------------ Edge case ------------------
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

  # ------------------ Aggregate ROC ------------------
  mean_tpr <- Reduce("+", tprs) / length(tprs)
  mean_tpr[length(mean_tpr)] <- 1

  mean_auc <- mean(aucs)
  sd_auc <- stats::sd(aucs)

  plot_df <- data.frame(
    FPR = mean_fpr,
    TPR = mean_tpr
  )
  plot_df <- plot_df[stats::complete.cases(plot_df), ]

  # ------------------ Plot ------------------
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = FPR, y = TPR)) +
    ggplot2::geom_line(color = "#B24745", linewidth = 1.8, alpha = 0.95) +
    ggplot2::geom_abline(
      slope = 1, intercept = 0,
      linetype = "dashed", color = "black", linewidth = 0.8
    ) +
    ggplot2::labs(
      title = title,
      subtitle = sprintf("Mean AUC = %.3f ± %.3f", mean_auc, sd_auc),
      x = "1 – Specificity",
      y = "Sensitivity"
    ) +
    ggplot2::theme_minimal(base_size = 16) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (isTRUE(export)) {
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
