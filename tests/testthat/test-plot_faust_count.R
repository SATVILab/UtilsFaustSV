test_that('plot_faust_count runs', {

  proj_path <- usethis::proj_path('/tests/testthat')
  # character pop
  plot_faust_count(project_path = proj_path,
                   pop = c("CD4" = "-"))
  path_plot <- file.path(proj_path, 'faustData', 'plotData', 'pop_stats', 'CD4-.png')
  expect_true(file.exists(path_plot))
  unlink(file.path(proj_path, 'faustData', 'plotData'), recursive = TRUE)
  unlink(file.path(proj_path, 'faustData', 'plotData'), recursive = TRUE)
  # list pop
  plot_faust_count(project_path = proj_path,
                   pop = list(c("CD4" = "-"), c("CD4" = "+", "HLA-DR-beads" = "+"),
                              c("CD4" = "+", "HLA-DR-beads" = "-")))
  path_plot <- file.path(proj_path, 'faustData', 'plotData', 'pop_stats', 'Total counts by subset.png')
  expect_true(file.exists(path_plot))

  unlink(file.path(proj_path, 'faustData', 'plotData'), recursive = TRUE)
  unlink(file.path(proj_path, 'faustData', 'plotData'), recursive = TRUE)
})


