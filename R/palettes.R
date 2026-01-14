#' CPMM default color palette for diagnostic and cohort groups
#'
#' Returns a named character vector of hex colors used consistently across
#' CPMM plots (diagnostic categories and cohort-level strata).
#'
#' @return Named character vector of colors.
#' @examples
#' pal <- cpmm_palette()
#' pal["MCI"]
#' pal["Status_change"]
#' @export
cpmm_palette <- function() {
  c(
    # Diagnostic groups
    Normal_only   = "#00A1D5",
    Normal        = "#00A1D5",
    CU            = "#00A1D5",

    SCD           = "#DF8F44",
    MCI           = "#B24745",
    AD_Dementia   = "#79AF97",
    FTD_Dementia  = "#6A6599",

    # Cohort / trajectory groups
    Status_change = "#DF8F44",
    Converters    = "#B24745",
    Abnormal_only = "#B24745"
  )
}

#' ggplot2 color scale using the CPMM palette
#'
#' @inheritParams ggplot2::scale_color_manual
#' @export
scale_color_cpmm <- function(...) {
  ggplot2::scale_color_manual(values = cpmm_palette(), ...)
}

#' @rdname scale_color_cpmm
#' @export
scale_fill_cpmm <- function(...) {
  ggplot2::scale_fill_manual(values = cpmm_palette(), ...)
}
