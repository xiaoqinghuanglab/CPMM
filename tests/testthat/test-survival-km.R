test_that("plot_km_with_threshold returns a patchwork/ggplot object", {

  set.seed(1)

  make_df <- function(n_id, label, mean_val) {
    n <- n_id * 3
    data.frame(
      SUBID = rep(paste0(label, seq_len(n_id)), each = 3),
      PROCEDURE_AGE = rep(60:62, times = n_id),
      ONSET_AGE = rep(59, n),
      CATEGORY = label,
      NFL = rnorm(n, mean_val, 2)
    )
  }

  df_normal <- make_df(4, "Normal", 8)
  df_status <- rbind(
    make_df(5, "MCI", 15),
    make_df(5, "Other", 9)
  )

  p <- plot_km_with_threshold(
    biomarker_name = "NFL",
    threshold = 12,
    groups = list(
      Normal = list(df = df_normal),
      Status_change = list(df = df_status),
      MCI = list(
        df = df_status,
        filter = rlang::expr(CATEGORY == "MCI")
      )
    ),
    time_points = seq(-2, 4, by = 2),
    save = FALSE
  )

  expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))
})
