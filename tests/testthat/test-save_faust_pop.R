# ================================
# Create data for tests
# ================================

.create_data_for_tests <- function(){
  # Note that initial GatingSet was gs_cytof_acs

  # directories
  # ----------------------------
  dir_test <- DataPackageR::project_path('tests/testthat/faustData')
  if(!dir.exists(dir_test)) dir.create(dir_test)
  dir_inst <- DataPackageR::project_extdata_path('faustData')

  # analysis map
  # -----------------------------
  analysis_map <- readRDS(file.path(dir_inst, "/metaData/analysisMap.rds"))
  analysis_map <- analysis_map[1:2,,drop = FALSE]
  dir_analysis_map <- file.path(dir_test, "metaData")
  if(!dir.exists(dir_analysis_map)) dir.create(dir_analysis_map)
  saveRDS(analysis_map, file.path(dir_analysis_map, "analysisMap.rds"))

  # GatingSet
  # ------------------------------
  dir_gs <- "C:/Users/migue/OneDrive - University of Cape Town/Work/PhD/Data/gs_cytof_acs"
  gs <- flowWorkspace::load_gs(dir_gs)


  # flowFrames and FAUST annotations
  # ------------------------------
  for(i in 1:2){

    # get, subset and save flowFrame
    fr <- flowWorkspace::gh_pop_get_data(gs[[i]])
    ex <- flowCore::exprs(fr)
    ind_vec <- 1e4 * 1:15 + 1
    ex <- ex[ind_vec,]
    flowCore::exprs(fr) <- ex
    sample <- gs[[i]]@name
    dir_fcs <- file.path(stringr::str_remove(dir_test, "/faustData"), "fcsTest")
    if(!dir.exists(dir_fcs)) dir.create(dir_fcs)
    flowCore::write.FCS(fr, file.path(dir_fcs, sample))

    # get, subset and save annotations
    path_ann <- file.path(dir_inst, "sampleData",
                          sample, "faustAnnotation.csv")
    faust_ann_tbl <- utils::read.table(file = path_ann,
                                       header = FALSE, sep = "`",
                                       stringsAsFactors = FALSE)[ind_vec,1,drop = FALSE]
    dir_faust_ann <- file.path(dir_test, "sampleData",
                               sample)
    if(!dir.exists(dir_faust_ann)) dir.create(dir_faust_ann, recursive = TRUE)
    readr::write_csv(faust_ann_tbl, file.path(dir_faust_ann, "faustAnnotation.csv"),
                     append = TRUE, col_names = FALSE)
  }

  # save GatingSet
  fs <- ncdfFlow::read.ncdfFlowSet(list.files(file.path(stringr::str_remove(dir_test, "/faustData"),
                                                        "fcsTest"),
                                              pattern = ".fcs", full.names = TRUE))
  gs_test <- flowWorkspace::GatingSet(fs)
  flowWorkspace::save_gs(gs_test, path = DataPackageR::project_path('tests/testthat/gs_test'))

}

if(FALSE) .create_data_for_tests() # run only if need be

# ================================
# Test save_faust_pop function
# ================================

testthat::test_that("save_faust_pop works correctly",{

  # ==========================================
  # Preparation
  # ==========================================

  # required objects
  dir_proj <- usethis::proj_path(ext = "/tests/testthat")
  gs <- flowWorkspace::load_gs(DataPackageR::project_path('tests/testthat/gs_test'))

  # ==========================================
  # Sample1 - CD8-IgD~2~
  # ==========================================

  # run code for one sample
  save_faust_pop(project_path = dir_proj,
                 pop = list("CD8-IgD" = 2),
                 sample = 1,
                 gs = gs)

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
  save_faust_pop(project_path = dir_proj,
                 pop = list("CD8-IgD" = 1),
                 sample = 1,
                 gs = gs)

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
  save_faust_pop(project_path = dir_proj,
                 pop = list("CD3" = 2),
                 sample = "01-0993 D0 AND 07-1147 DAY0-pid1_ebv-debeaded_2.fcs",
                 gs = gs)

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
  save_faust_pop(project_path = dir_proj,
                 pop = list("CD4" = 1),
                 sample = NULL,
                 gs = gs)

  # check that folders are created correctly
  dir_fcs <- file.path(dir_proj, "faustData", "fcsData")
  expect_true(dir.exists(file.path(dir_fcs, "CD4~1~")))
  expect_identical(list.files(file.path(dir_fcs, "CD4~1~")), c("01-0993 D0 AND 07-1147 DAY0-pid1_ebv-debeaded_2.fcs",
                                                               "01-0993 D0 AND 07-1147 DAY0-pid1_mtbaux-debeaded_2.fcs"))

  # check that object was saved correctly
  fr <- flowCore::read.FCS(file.path(dir_fcs, "CD4~1~", "01-0993 D0 AND 07-1147 DAY0-pid1_mtbaux-debeaded_2.fcs"))
  ex <- flowCore::exprs(fr)
  expect_identical(nrow(ex), 8L)

  fr <- flowCore::read.FCS(file.path(dir_fcs, "CD4~1~", "01-0993 D0 AND 07-1147 DAY0-pid1_ebv-debeaded_2.fcs"))
  ex <- flowCore::exprs(fr)
  expect_identical(nrow(ex), 11L)

  unlink(dir_fcs, recursive = TRUE)

})

