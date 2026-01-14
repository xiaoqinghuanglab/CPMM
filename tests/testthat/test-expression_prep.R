test_that("prepare_combined_expression pivots, coerces, and labels correctly", {

  genes <- c("G1", "G2")

  df_all <- tibble::tibble(
    CATEGORY = c("Normal", "SCD", "MCI", "AD Dementia"),
    G1 = c("1.0", "2.0", "3.0", "4.0"),
    G2 = c("5.0", "6.0", "7.0", "8.0")
  )

  df_normal_only   <- df_all[1, , drop = FALSE]
  df_status_change <- df_all[2:3, , drop = FALSE]
  df_abnormal_only <- df_all[4, , drop = FALSE]

  out <- prepare_combined_expression(
    inputs = list(
      list(
        df = df_normal_only,
        label = "Normal_only"
      ),
      list(
        df = df_status_change,
        label = "Status_change",
        category_col = "CATEGORY",
        categories = c("SCD", "MCI")
      ),
      list(
        df = df_abnormal_only,
        label = "Abnormal_only"
      )
    ),
    subset_genes = genes
  )

  # ---- structural checks ----
  expect_true(all(c("Gene", "Expression", "Source") %in% names(out)))
  expect_s3_class(out$Gene, "factor")
  expect_true(is.ordered(out$Gene))
  expect_equal(levels(out$Gene), genes)
  expect_type(out$Expression, "double")

  # ---- source labels ----
  expect_true(all(
    c("Normal_only", "Status_change", "Abnormal_only") %in% unique(out$Source)
  ))

  # ---- category filtering worked ----
  expect_true(
    nrow(out[out$Source == "Status_change", ]) == length(genes) * 2
  )
})
