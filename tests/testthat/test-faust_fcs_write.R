# ================================
# Create data for tests
# ================================

.create_data_for_tests <- function() {
  # Note that initial GatingSet was gs_cytof_acs

  # directories
  # ----------------------------
  dir_test <- here::here("tests/testthat/faustData")
  if (!dir.exists(dir_test)) dir.create(dir_test)
  dir_inst <- here::here("inst", "extdata", "faustData")

  # analysis map
  # -----------------------------
  analysis_map <- readRDS(file.path(dir_inst, "/metaData/analysisMap.rds"))
  analysis_map <- analysis_map[1:2, , drop = FALSE]
  dir_analysis_map <- file.path(dir_test, "metaData")
  if (!dir.exists(dir_analysis_map)) dir.create(dir_analysis_map)
  saveRDS(analysis_map, file.path(dir_analysis_map, "analysisMap.rds"))

  # flowFrames and FAUST annotations
  # ------------------------------
  # copy original FCS files across
  # (not repeatable off Miguel Rodo's computer,
  # but that does not matter)
  # ------------------------------
  dir_comp_acs <- "C:/Users/migue/Work/Projects/SATVI/CompACS"
  dir_fcs <- file.path(dir_comp_acs, "Output/OutputDataTidyACSCyTOFPreprocess/final/fcs") # nolint
  fn_vec <- c(
    "01-0993 D0 AND 07-1147 DAY0-pid1_ebv-debeaded_2.fcs",
    "01-0993 D0 AND 07-1147 DAY0-pid1_mtbaux-debeaded_2.fcs"
  )
  if (!requireNamespace("here", quietly = TRUE)) {
    utils::install.packages("here")
  }

  if (!dir.exists(here::here("tests/testthat/fcsSource"))) {
    dir.create(here::here("tests/testthat/fcsSource"), recursive = TRUE)
  }

  file.copy(
    file.path(dir_fcs, fn_vec),
    here::here("tests/testthat/fcsSource", fn_vec)
  )
  fs <- flowCore::read.flowSet(
    file.path(dir_fcs, fn_vec)
  )

  for (i in 1:2) {
    # get, subset and save flowFrame
    fr <- fs[[i]]
    ex <- flowCore::exprs(fr)
    ind_vec <- 1e4 * 1:15 + 1
    ex <- ex[ind_vec, ]
    flowCore::exprs(fr) <- ex
    sample <- fn_vec[i]
    dir_fcs <- file.path(stringr::str_remove(dir_test, "/faustData"), "fcsTest")
    if (!dir.exists(dir_fcs)) dir.create(dir_fcs)
    flowCore::write.FCS(fr, file.path(dir_fcs, sample))

    # get, subset and save annotations
    path_ann <- file.path(
      dir_inst, "sampleData",
      sample, "faustAnnotation.csv"
    )
    faust_ann_tbl <- utils::read.table(
      file = path_ann,
      header = FALSE, sep = "`",
      stringsAsFactors = FALSE
    )[ind_vec, 1, drop = FALSE]
    dir_faust_ann <- file.path(
      dir_test, "sampleData",
      sample
    )
    if (!dir.exists(dir_faust_ann)) dir.create(dir_faust_ann, recursive = TRUE)
    if (!requireNamespace("readr", quietly = TRUE)) {
      install.packages("readr")
    }
    readr::write_csv(
      faust_ann_tbl,
      file.path(dir_faust_ann, "faustAnnotation.csv"),
      append = FALSE, col_names = FALSE
    )
  }

  # save GatingSet
  cs <- flowWorkspace::load_cytoset_from_fcs(list.files(
    file.path(
      stringr::str_remove(dir_test, "/faustData"),
      "fcsTest"
    ),
    pattern = ".fcs", full.names = TRUE
  ))
  gs_test <- flowWorkspace::GatingSet(cs)
  flowWorkspace::save_gs(gs_test, path = here::here("tests/testthat/gs_test"))
}

if (!requireNamespace("here", quietly = TRUE)) {
  utils::install.packages("here")
}
if (!dir.exists(here::here("tests/testthat/gs_test"))) {
  .create_data_for_tests()
}

# ================================
# Test faust_fcs_write function
# ================================

testthat::test_that("faust_fcs_write works correctly", {
  # ==========================================
  # Preparation
  # ==========================================

  # required objects
  dir_proj <- usethis::proj_path(ext = "./tests/testthat")
  gs <- flowWorkspace::load_gs(here::here("tests/testthat/gs_test"))

  # ==========================================
  # Sample1 - CD8-IgD~2~
  # ==========================================

  # run code for one sample
  faust_fcs_write(
    project_path = dir_proj,
    pop = c("CD8-IgD" = 2),
    sample = 1,
    fr_source = gs
  )

  # check that folders are created correctly
  dir_fcs <- file.path(dir_proj, "faustData", "fcsData")
  expect_true(dir.exists(file.path(dir_fcs, "CD8-IgD~2~")))
  expect_identical(list.files(file.path(dir_fcs, "CD8-IgD~2~")), "01-0993 D0 AND 07-1147 DAY0-pid1_ebv-debeaded_2.fcs")

  # check that object was saved correctly
  fr <- flowCore::read.FCS(file.path(dir_fcs, "CD8-IgD~2~", "01-0993 D0 AND 07-1147 DAY0-pid1_ebv-debeaded_2.fcs"))
  ex_2 <- flowCore::exprs(fr)
  expect_identical(nrow(ex_2), 7L)

  # ==========================================
  # Sample1 - CD8-IgD~1~
  # ==========================================

  # run code for same sample but for level 1 of CD8-IgD
  faust_fcs_write(
    project_path = dir_proj,
    pop = c("CD8-IgD" = 1),
    sample = 1,
    fr_source = gs
  )

  # check that folders are created correctly
  dir_fcs <- file.path(dir_proj, "faustData", "fcsData")
  expect_true(dir.exists(file.path(dir_fcs, "CD8-IgD~1~")))
  expect_identical(list.files(file.path(dir_fcs, "CD8-IgD~1~")), "01-0993 D0 AND 07-1147 DAY0-pid1_ebv-debeaded_2.fcs")

  # check that object was saved correctly
  fr <- flowCore::read.FCS(file.path(dir_fcs, "CD8-IgD~1~", "01-0993 D0 AND 07-1147 DAY0-pid1_ebv-debeaded_2.fcs"))
  ex_1 <- flowCore::exprs(fr)
  expect_identical(nrow(ex_1), 8L)
  expect_identical(nrow(ex_2) + nrow(ex_1), 15L)

  unlink(dir_fcs, recursive = TRUE)

  # ==========================================
  # Sample "01-0993 D0 AND 07-1147 DAY0-pid1_ebv-debeaded_2.fcs" - CD3 ~ 2
  # ==========================================

  # run code for same sample but for level 1 of CD8-IgD
  faust_fcs_write(
    project_path = dir_proj,
    pop = c("CD3" = 2),
    sample = "01-0993 D0 AND 07-1147 DAY0-pid1_ebv-debeaded_2.fcs",
    fr_source = gs
  )

  # check that folders are created correctly
  dir_fcs <- file.path(dir_proj, "faustData", "fcsData")
  expect_true(dir.exists(file.path(dir_fcs, "CD3~2~")))
  expect_identical(list.files(file.path(dir_fcs, "CD3~2~")), "01-0993 D0 AND 07-1147 DAY0-pid1_ebv-debeaded_2.fcs")

  # check that object was saved correctly
  fr <- flowCore::read.FCS(file.path(dir_fcs, "CD3~2~", "01-0993 D0 AND 07-1147 DAY0-pid1_ebv-debeaded_2.fcs"))
  ex <- flowCore::exprs(fr)
  expect_identical(nrow(ex), 11L)

  unlink(dir_fcs, recursive = TRUE)

  # ==========================================
  # Both samples
  # ==========================================

  # run code for same sample but for level 1 of CD8-IgD
  faust_fcs_write(
    project_path = dir_proj,
    pop = c("CD4" = 1),
    sample = NULL,
    fr_source = gs
  )

  # check that folders are created correctly
  dir_fcs <- file.path(dir_proj, "faustData", "fcsData")
  expect_true(dir.exists(file.path(dir_fcs, "CD4~1~")))
  expect_identical(list.files(file.path(dir_fcs, "CD4~1~")), c(
    "01-0993 D0 AND 07-1147 DAY0-pid1_ebv-debeaded_2.fcs",
    "01-0993 D0 AND 07-1147 DAY0-pid1_mtbaux-debeaded_2.fcs"
  ))

  # check that object was saved correctly
  fr <- flowCore::read.FCS(file.path(dir_fcs, "CD4~1~", "01-0993 D0 AND 07-1147 DAY0-pid1_mtbaux-debeaded_2.fcs"))
  ex <- flowCore::exprs(fr)
  expect_identical(nrow(ex), 8L)

  fr <- flowCore::read.FCS(file.path(dir_fcs, "CD4~1~", "01-0993 D0 AND 07-1147 DAY0-pid1_ebv-debeaded_2.fcs"))
  ex <- flowCore::exprs(fr)
  expect_identical(nrow(ex), 11L)

  unlink(dir_fcs, recursive = TRUE)

  # ==========================================
  # Sample1 - CD8-IgD~1~ - CD8 level one or two
  # ==========================================

  # run code for same sample but for level 1 of CD8-IgD
  faust_fcs_write(
    project_path = dir_proj,
    pop = list(
      c("CD8-IgD" = 1),
      c("CD8-IgD" = 2)
    ),
    sample = 1,
    fr_source = gs
  )

  # check that folders are created correctly
  dir_fcs <- file.path(dir_proj, "faustData", "fcsData")
  expect_true(dir.exists(file.path(dir_fcs, "CD8-IgD~12~")))
  expect_identical(list.files(file.path(dir_fcs, "CD8-IgD~12~")), "01-0993 D0 AND 07-1147 DAY0-pid1_ebv-debeaded_2.fcs")

  # check that object was saved correctly
  fr <- flowCore::read.FCS(file.path(dir_fcs, "CD8-IgD~12~", "01-0993 D0 AND 07-1147 DAY0-pid1_ebv-debeaded_2.fcs"))
  ex_1 <- flowCore::exprs(fr)
  expect_identical(nrow(ex_1), 15L)

  unlink(dir_fcs, recursive = TRUE)

  # ==========================================
  # Sample1 - CD8-IgD~1~ - CD3 level 2
  # ==========================================

  # preliminary check that first pop is of size 4 and second of size 0
  # run code for same sample but for level 1 of CD8-IgD
  faust_fcs_write(
    project_path = dir_proj,
    pop = c("CD3" = 2, "CD8-IgD" = 1),
    sample = 1,
    fr_source = gs
  )

  # check that folders are created correctly
  dir_fcs <- file.path(dir_proj, "faustData", "fcsData")
  expect_true(dir.exists(file.path(dir_fcs, "CD3~2~CD8-IgD~1~")))
  expect_identical(list.files(file.path(dir_fcs, "CD3~2~CD8-IgD~1~")), "01-0993 D0 AND 07-1147 DAY0-pid1_ebv-debeaded_2.fcs")

  # check that object was saved correctly
  fr <- flowCore::read.FCS(file.path(dir_fcs, "CD3~2~CD8-IgD~1~", "01-0993 D0 AND 07-1147 DAY0-pid1_ebv-debeaded_2.fcs"))
  ex_1 <- flowCore::exprs(fr)
  expect_identical(nrow(ex_1), 4L)

  unlink(dir_fcs, recursive = TRUE)

  # run code for same sample but for level 1 of CD8-IgD
  faust_fcs_write(
    project_path = dir_proj,
    pop = c("CD3" = 1, "CD8-IgD" = 2),
    sample = 1,
    fr_source = gs
  )

  # check that folders are created correctly
  dir_fcs <- file.path(dir_proj, "faustData", "fcsData")
  expect_true(dir.exists(file.path(dir_fcs, "CD3~1~CD8-IgD~2~")))
  expect_identical(list.files(file.path(dir_fcs, "CD3~1~CD8-IgD~2~")), "01-0993 D0 AND 07-1147 DAY0-pid1_ebv-debeaded_2.fcs")

  # check that object was saved correctly
  fr <- flowCore::read.FCS(file.path(dir_fcs, "CD3~1~CD8-IgD~2~", "01-0993 D0 AND 07-1147 DAY0-pid1_ebv-debeaded_2.fcs"))
  ex_1 <- flowCore::exprs(fr)
  expect_identical(nrow(ex_1), 1L)
  expect_identical(ex_1[1, "Ba137Di"][[1]], NaN)

  unlink(dir_fcs, recursive = TRUE)

  # run code for same sample but for level 1 of CD8-IgD
  faust_fcs_write(
    project_path = dir_proj,
    pop = list(
      c("CD3" = 2, "CD8-IgD" = 1),
      c("CD3" = 1, "CD8-IgD" = 2)
    ),
    sample = 1,
    fr_source = gs
  )

  # check that folders are created correctly
  dir_fcs <- file.path(dir_proj, "faustData", "fcsData")
  expect_true(dir.exists(file.path(dir_fcs, "CD3~21~CD8-IgD~12~")))
  expect_identical(list.files(file.path(dir_fcs, "CD3~21~CD8-IgD~12~")), "01-0993 D0 AND 07-1147 DAY0-pid1_ebv-debeaded_2.fcs")

  # check that object was saved correctly
  fr <- flowCore::read.FCS(file.path(dir_fcs, "CD3~21~CD8-IgD~12~", "01-0993 D0 AND 07-1147 DAY0-pid1_ebv-debeaded_2.fcs"))
  ex_1 <- flowCore::exprs(fr)
  expect_identical(nrow(ex_1), 4L)

  unlink(dir_fcs, recursive = TRUE)
})

test_that(".is_faust_ann_a_match_for_marker works when first entry is 0_0_0_0_0", {
  faust_ann <- data.frame(x = c("0_0_0_0_0", "CD3~1~", "CD3~2~"))
  match_vec <- .is_faust_ann_a_match_for_marker(faust_ann = faust_ann, marker = "CD3", level = "1")
  expect_identical(match_vec, c(FALSE, TRUE, FALSE))
})

test_that(".is_faust_ann_a_match_for_marker works when all entries are 0_0_0_0_0", {
  faust_ann <- data.frame(x = c("0_0_0_0_0", "0_0_0_0_0", "0_0_0_0_0"))
  match_vec <- .is_faust_ann_a_match_for_marker(faust_ann = faust_ann, marker = "CD3", level = "1")
  expect_identical(match_vec, rep(FALSE, 3))
})
