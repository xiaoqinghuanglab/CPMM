test_that("fit_cpmm_all_proteins returns expected columns", {
  skip_on_cran()
  set.seed(1)
  n_id <- 12
  make_df <- function(ids, yso_min = -5, yso_max = 5, abnormal = FALSE, normal = FALSE) {
    df <- data.frame(
      SUBID = rep(ids, each = 4),
      years_since_onset = runif(length(ids) * 4, yso_min, yso_max),
      SEX = sample(c("M","F"), length(ids) * 4, replace = TRUE),
      BASELINE_AGE = rnorm(length(ids) * 4, 65, 6)
    )
    # Simple generative slopes: before +0.03/yr; after -0.05/yr
    before <- pmax(0, -df$years_since_onset)
    after  <- pmax(0,  df$years_since_onset)
    signal <- 0.5 + 0.03*before - 0.05*after + 0.2*(df$SEX == "M") + 0.01*df$BASELINE_AGE
    noise  <- rnorm(nrow(df), 0, 0.05)
    df$P1 <- signal + noise
    df$P2 <- signal*0.5 + rnorm(nrow(df), 0, 0.05)
    df
  }
  ids_status  <- paste0("S", 1:6)
  ids_normal  <- paste0("N", 1:3)
  ids_abnorm  <- paste0("A", 1:3)

  df_status   <- make_df(ids_status)
  df_normal   <- make_df(ids_normal, yso_min = -6, yso_max = 0)   # mostly before
  df_abnormal <- make_df(ids_abnorm, yso_min = 0,  yso_max = 6)   # mostly after

  res <- fit_cpmm_all_proteins(
    df_status_change = df_status,
    df_normal = df_normal,
    df_abnormal = df_abnormal,
    protein_list = c("P1","P2"),
    covariates = c("SEX","BASELINE_AGE")
  )

  expect_true(all(c("Protein","Beta 1","Beta 2","Beta 3","Beta 4",
                    "SE Beta 1","SE Beta 2","SE Beta 3","SE Beta 4",
                    "Intercept","AIC Status","BIC Status","MSE Status") %in% names(res)))
  expect_equal(nrow(res), 2L)
})
