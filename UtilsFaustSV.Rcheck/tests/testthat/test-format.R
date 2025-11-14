test_that(".faust_level_convert works", {
  expect_error(.faust_marker_format_level_to_tilde(lvl = "abc"))
  expect_identical(
    .faust_marker_format_level_to_tilde(marker = "CD3", lvl = "~1~2~"), "~1~2~"
  )
  expect_identical(
    .faust_marker_format_level_to_tilde(marker = "CD3", lvl = "Dim", n_lvl = 3),
    "~2~3~"
  )
  expect_error(
    .faust_marker_format_level_to_tilde(marker = "CD3", lvl = "Dim", n_lvl = 2)
  )
  # test common conversion
  expect_error(
    .faust_marker_format_level_to_common(marker = "CD3", lvl = "abc")
  )
  expect_identical(.faust_marker_format_level_to_common(lvl = "+"), "+")
  expect_identical(.faust_marker_format_level_to_common(lvl = "~1~2~"), "-")
  # pop
  # ----------------------

  # to common
  expect_identical(
    .faust_pop_format_level_to_common("CD4+CD3-"), "CD4+CD3-"
  )
  expect_identical(
    .faust_pop_format_level_to_common("CD4+CD3-"), "CD4+CD3-"
  )
  expect_identical(
    .faust_pop_format_level_to_common("CD4~1~3~CD3~4~4~"), "CD4-CD3VeryBright"
  )

  # to tilde
  expect_identical(
    .faust_pop_format_level_to_tilde("CD4~1~3~CD3~2~2~"), "CD4~1~3~CD3~2~2~"
  )
  expect_identical(
    .faust_pop_format_level_to_tilde(
      pop = "CD4DimCD3VeryBrightCD8-",
      marker_to_n_lvl = c("CD4" = 3, "CD8" = 2, "CD3" = 4),
    ),
    "CD4~2~3~CD3~4~4~CD8~1~2~"
  )
  marker_to_n_lvl <- c("CD4" = 3, "CD8-IgD" = 2, "CD3" = 4)
  pop <- "CD4DimCD3VeryBrightCD8-IgD-"
  expect_identical(
    .faust_pop_format_level_to_tilde(
      pop = "CD4DimCD3VeryBrightCD8-IgD-",
      marker_to_n_lvl = c("CD4" = 3, "CD8-IgD" = 2, "CD3" = 4),
    ),
    "CD4~2~3~CD3~4~4~CD8-IgD~1~2~"
  )
})
