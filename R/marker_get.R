#' @title Get the number of levels used for each marker FAUST used to cluster
#'
#' @description
#' Gives the markers FAUST used and the number of levels used per such marker.
#' Note that it only works after the SCAMP clustering has been performed, i.e.
#' not simply after the annotations have been approved.
#'
#' @inheritParams faust_fcs_write
#'
#' @return \code{Character vector}
#' @export
faust_marker_get_usage <- function(project_path) {
  dir_metadata <- file.path(project_path, "faustData", "metaData")
  cluster_name <- readRDS(file.path(
    dir_metadata, "scampClusterNames.rds"
  ))[1]

  matches <- stringr::str_match_all(cluster_name, "([^~]+)~\\d+~(\\d+)~")[[1]]

  stats::setNames(as.numeric(matches[, 3]), matches[, 2])
}

faust_count_get <- function(project_path, exhaustive = FALSE) {
  path_mat <- file.path(
    project_path,
    switch(as.character(exhaustive),
      "TRUE" = "exhaustiveFaustCountMatrix.rds",
      "FALSE" = "faustCountMatrix.rds",
      stop("exhaustive argument must be TRUE or FALSE")
    )
  )
  readRDS(path_mat)
}
