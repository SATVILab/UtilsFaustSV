library(testthat)
test_that("get_pop_counts works", {

  dir_proj <- usethis::proj_path('/tests/testthat')

  full_tbl <- get_pop_counts(dir_proj)
  cd4p_tbl <- get_pop_counts(project_path = dir_proj,
                             pop = c("CD4" = "+"))
  cd4n_tbl <- get_pop_counts(project_path = dir_proj,
                             pop = c("CD4" = "-"))
  expect_identical(ncol(cd4p_tbl), 5L)
  expect_identical(ncol(cd4n_tbl), 8L)

  # check that we have the expected number of columns for cd4p and cd4n
  expect_identical(length(intersect(colnames(cd4p_tbl),
                                    colnames(cd4n_tbl))), 5L)

  expect_identical(length(intersect(colnames(cd4p_tbl),
                                    colnames(cd4n_tbl))), 5L)

  # check that cd4p and cd4n together return all original columns
  expect_identical(ncol(full_tbl),
                   ncol(dplyr::left_join(cd4p_tbl, cd4n_tbl,
                                         by = c('sample', 'exp_unit', 'tot_count',
                                                'tot_count_classified'))))

  # list
  cd4p_tbl <- get_pop_counts(project_path = dir_proj,
                             pop = list(c("CD4" = "+")))
  cd4n_tbl <- get_pop_counts(project_path = dir_proj,
                             pop = list(c("CD4" = "-")))
  cd4_tbl <- get_pop_counts(project_path = dir_proj,
                             pop = list(c("CD4" = "+"),
                                        c("CD4" = "-")))


  expect_identical(colnames(cd4p_tbl)[[5]], "CD4+")
  expect_identical(colnames(cd4n_tbl)[[5]], "CD4-")
  expect_identical(colnames(cd4_tbl)[5:6], c("CD4+", "CD4-"))

  cd4n_hladr_tbl <- get_pop_counts(project_path = dir_proj,
                            pop = list(c("CD4" = "-", "HLA-DR-beads" = "-"),
                                       c("CD4" = "-")))
  expect_identical(cd4n_hladr_tbl$`CD4-HLA-DR-beads-`, c(456, 1578))
  expect_identical(cd4n_hladr_tbl$`CD4-`, c(799, 2286))
})
