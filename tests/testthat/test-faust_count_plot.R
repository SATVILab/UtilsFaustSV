test_that("faust_count_plot runs", {
  proj_path <- usethis::proj_path("./tests/testthat")
  # character pop
  faust_count_plot(
    project_path = proj_path,
    pop = c("CD4" = "-")
  )
  path_plot <- file.path(
    proj_path, "faustData", "plotData", "pop_stats", "CD4-.png"
  )
  expect_true(file.exists(path_plot))
  unlink(file.path(proj_path, "faustData", "plotData"), recursive = TRUE)
  unlink(file.path(proj_path, "faustData", "plotData"), recursive = TRUE)
  # list pop
  faust_count_plot(
    project_path = proj_path,
    pop = list(
      c("CD4" = "-"), c("CD4" = "+", "HLA-DR-beads" = "+"),
      c("CD4" = "+", "HLA-DR-beads" = "-")
    )
  )
  path_plot <- file.path(
    proj_path, "faustData", "plotData",
    "pop_stats", "Frequencies of subsets.png"
  )
  expect_true(file.exists(path_plot))

  # NULL pop
  faust_count_plot(
    project_path = proj_path,
    pop = NULL
  )
  path_plot <- file.path(
    proj_path, "faustData", "plotData",
    "pop_stats", "Frequencies of all FAUST subsets.png"
  )
  expect_true(file.exists(path_plot))

  # large plots
  faust_count_plot(
    project_path = proj_path,
    pop = NULL, p_height = 130, p_width = 130
  )
  path_plot <- file.path(
    proj_path, "faustData", "plotData",
    "pop_stats", "Frequencies of all FAUST subsets.png"
  )
  expect_true(file.exists(path_plot))

  unlink(file.path(proj_path, "faustData", "plotData"), recursive = TRUE)
  unlink(file.path(proj_path, "faustData", "plotData"), recursive = TRUE)

  # create exhaustive count matrix
  # mat <- matrix(c(1,2,3,4), byrow = TRUE, ncol = 2)
  # colnames(mat) <- c("CD3+CD4+", "CD3-CD4+")
  # rownames(mat) <- c("sample1", "sample2")
  # saveRDS(
  # mat, file.path(proj_path,
  # "faustData", "exhaustiveFaustCountMatrix.rds")
  # )
  expect_error(faust_count_plot(
    project_path = proj_path,
    pop = c("CD8-IgD" = "-"),
    exhaustive = TRUE
  ))
  faust_count_plot(
    project_path = proj_path,
    pop = c("CD3" = "+"),
    exhaustive = TRUE
  )
  path_plot <- file.path(
    proj_path, "faustData", "plotData", "pop_stats-exhaustive", "CD3+.png"
  )
  expect_true(file.exists(path_plot))
  unlink(file.path(proj_path, "faustData", "plotData"), recursive = TRUE)
})
