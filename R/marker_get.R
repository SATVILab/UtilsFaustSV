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

  # all channels considered
  active_chnl_vec <- readRDS(file.path(dir_metadata, "activeChannels.rds"))

  # filter all channels considered to get all channels used
  clustered_chnl_vec <- NULL
  for (x in active_chnl_vec) {
    if (stringr::str_detect(cluster_name, x)) {
      clustered_chnl_vec <- c(clustered_chnl_vec, x)
    }
  }

  # get observed levels for each channel used
  n_level_vec <- stats::setNames(
    rep(NA, length(clustered_chnl_vec)), clustered_chnl_vec
  )
  for (i in seq_along(clustered_chnl_vec)) {
    n_level_loc <- stringr::str_locate(
      cluster_name, clustered_chnl_vec[i]
    )[, "end"][[1]] + 4

    n_level_vec[i] <- as.numeric(
      stringr::str_sub(cluster_name, n_level_loc, n_level_loc)
    )
  }
  n_level_vec
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
