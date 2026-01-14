test_that("cpmm_palette returns named hex colors", {
  pal <- cpmm_palette()

  # must be a named character vector
  expect_type(pal, "character")
  expect_true(!is.null(names(pal)))

  # all values must be valid hex colors
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", unname(pal))))

  # key diagnostic groups exist
  expect_true(all(c(
    "Normal_only", "Normal", "CU",
    "SCD", "MCI", "AD_Dementia", "FTD_Dementia"
  ) %in% names(pal)))

  # key cohort groups exist
  expect_true(all(c(
    "Status_change", "Converters", "Abnormal_only"
  ) %in% names(pal)))

  # canonical colors are correct
  expect_equal(pal[["MCI"]], "#B24745")
  expect_equal(pal[["Normal"]], "#00A1D5")
  expect_equal(pal[["Status_change"]], "#DF8F44")
})


test_that("scale_color_cpmm returns a ggplot2 scale", {
  sc <- scale_color_cpmm()
  expect_s3_class(sc, "Scale")
})

test_that("scale_fill_cpmm returns a ggplot2 scale", {
  sc <- scale_fill_cpmm()
  expect_s3_class(sc, "Scale")
})
