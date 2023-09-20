test_that("faust_get_marker works", {
  marker_nlevel_vec <- faust_get_marker(usethis::proj_path("/tests/testthat"))
  marker_vec <- c(
    "CD33", "CD7", "CCR7", "CD8-IgD", "HLA-DR-beads", "CD14", "CD27",
    "CD4", "CD16", "CD20", "TCRgd-CD19", "CD3", "CD45RA", "CXCR5"
  )
  expect_identical(
    marker_nlevel_vec, stats::setNames(rep(2, length(marker_vec)), marker_vec)
  )
})
