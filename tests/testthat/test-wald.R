test_that("compute_wald_test works and ranks correctly (CPMM)", {
  df <- tibble::tibble(
    Protein = c("A", "B", "C"),
    `Beta 1` = c( 0.30,  0.10, 0.00),
    `SE Beta 1` = c(0.10, 0.10, 0.10),
    `Beta 2` = c(-0.10,  0.08, 0.00),
    `SE Beta 2` = c(0.10, 0.10, 0.20)
  )

  res <- compute_wald_test(df, adjust_p = TRUE, alpha = 0.05)

  # Expected columns
  expect_true(all(c(
    "Protein",
    "Beta 1", "SE Beta 1",
    "Beta 2", "SE Beta 2",
    "Wald Statistic",
    "P-value",
    "Adjusted P-value",
    "Significant",
    "Rank"
  ) %in% names(res)))

  # Row count preserved
  expect_equal(nrow(res), 3L)

  # Protein A should rank first (largest slope difference)
  expect_equal(res$Protein[1], "A")

  # Rank should be sequential
  expect_equal(res$Rank, seq_len(3))

  # Significance flag exists and is logical
  expect_type(res$Significant, "logical")

  # Test without p-value adjustment
  res2 <- compute_wald_test(df, adjust_p = FALSE)

  expect_true(all(c(
    "Protein",
    "Wald Statistic",
    "P-value",
    "Rank"
  ) %in% names(res2)))

  expect_equal(nrow(res2), 3L)
  expect_equal(res2$Rank, seq_len(3))
})
