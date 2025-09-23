test_that("plot_cpmm returns a ggplot object", {
  skip_on_cran()
  set.seed(1)
  mk <- function(n_id, y_min, y_max) {
    n <- n_id * 4
    df <- data.frame(
      SUBID = rep(paste0("S", 1:n_id), each = 4),
      years_since_onset = runif(n, y_min, y_max),
      SEX = sample(c("M","F"), n, TRUE),
      BASELINE_AGE = rnorm(n, 65, 6)
    )
    df$before_onset <- pmax(0, -df$years_since_onset)
    df$after_onset  <- pmax(0,  df$years_since_onset)
    df$P1 <- 0.5 + 0.03*df$before_onset - 0.05*df$after_onset + 0.2*(df$SEX=="M") + 0.01*df$BASELINE_AGE + rnorm(n,0,0.05)
    df
  }
  df_status  <- mk(6, -6, 6)
  df_normal  <- mk(4, -6, 0)
  df_abnorm  <- mk(4,  0, 6)

  p <- plot_cpmm(df_status, df_normal, df_abnorm, protein = "P1")
  expect_s3_class(p, "gg")
})
