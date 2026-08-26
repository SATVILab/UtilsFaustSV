test_that("faust_marker_get_usage works", {
  marker_nlevel_vec <- faust_marker_get_usage(
    testthat::test_path()
  )
  marker_vec <- c(
    "CD33", "CD7", "CCR7", "CD8-IgD", "HLA-DR-beads", "CD14", "CD27",
    "CD4", "CD16", "CD20", "TCRgd-CD19", "CD3", "CD45RA", "CXCR5"
  )
  expected <- stats::setNames(rep(2, length(marker_vec)), marker_vec)
  expect_identical( # order doesn't matter, contents does
    marker_nlevel_vec[order(names(marker_nlevel_vec))],
    expected[order(names(expected))]
  )
})
