#' @title Save FAUST subset as an FCS file
#'
#' @description
#' Save a specified FAUST-identified population to \code{project_path/faustData/fcsData}
#' as an FCS file for all samples gated or just a specified subset.
#'
#' @param project_path character. Path to directory containing FAUST output.
#' @param gs GatingSet. Original \code{GatingSet} object to which \code{faust} was applied.
#' @param pop \code{list} or \code{named character vector}. If a \code{character vector},
#' then all cells matching the set of marker levels are returned. If a \code{list}, then each
#' element must be a \code{character vector}, and then
#' all cells matching either of these specified sets of marker levels are returned.
#' @param sample integer vector or character vector. If an integer vector, then it
#' specifies the indices of the samples for which to save output. If character,
#' then it specifies the names of the sames to save output for. If \code{NULL},
#' then the output for every sample is saved. Default is \code{NULL}.
#' @param trans function. If supplied, this function is applied to the expression data.
#' Useful for back-transformation, as the original
#' GatingSet . If \code{NULL}, then no transformation is applied. Default is \code{NULL}.
#' Not working at the moment.
#' @param saveFCS Boolean. If \code{TRUE}, then FCS files of the selected cell population(s)
#' and sample(s) are saved to project_path/faustData/fcsData. Default is \code{FALSE}.
#' population of all the selected samplesare saved to project_path/faustData/gsData. Default is \code{FALSE}.
#'
#' @return \code{invisible(TRUE)}. Side effect is the saved FCS file.
#' @examples save_faust_pop(project_path = "", pop = list("CD3" = 2),
#' gs = gs, sample = 1)
#' @export
save_faust_pop <- function(project_path,
                           pop,
                           gs,
                           sample = NULL,
                           trans = NULL){

  # =============================
  # Check
  # =============================

  if(any(c(missing(project_path), missing(gs), missing(pop)))){
    stop("Each of project_path, gs and pop parameters must have arguments.")
  }

  # =============================
  # Preparation
  # =============================

  # get names of samples to save
  # ----------------------------------
  ## analysisMap contains the columns sampleName, experimentalUnit
  # and impH. The sampleName is the FCS file name, the
  # experimental unit specifies the group for gating and
  # the impH column specifies from where missing annotation boundaries
  # are imputed.
  active_sample_vec <- readRDS(paste0(project_path, "/faustData/metaData/analysisMap.rds"))[,"sampleName"]
  sel_sample_vec <- switch(typeof(sample),
                           "integer" = ,
                           "double" = active_sample_vec[sample],
                           "NULL" = active_sample_vec,
                           "character" = active_sample_vec[vapply(active_sample_vec, function(x) x %in% sample,
                                                                  logical(length(sample)))],
                           stop("Incorrect specification of sample parameter."))

  # get vector of sample names in gs
  sample_name_vec <- sapply(seq_along(gs), function(i) gs[[i]]@name)

  # create directory to save to
  # --------------------------------------

  # get concatenated name of population
  pop_name <- ""
  for(i in seq_along(pop)) pop_name <- paste0(pop_name,
                                              names(pop)[i],
                                              "~",
                                              pop[[i]],
                                              "~")

  #  create directory
  dir_save <- file.path(project_path,
                        "faustData",
                        "fcsData",
                        pop_name)
  if(!dir.exists(dir_save)) dir.create(dir_save, recursive = TRUE)

  # =============================
  # Save FAUST pops
  # =============================

  .save_faust_pop(project_path = project_path, gs = gs,
                  sample_name = sample_name_vec,
                  sel_sample = sel_sample_vec,
                  pop = pop,
                  dir_save = dir_save)

}

#' @title Save FAUST subset as an FCS file
#'
#' @inheritParams save_faust_pop
#' @param sel_sample \code{character vector}. Character vector specifying
#' the names of the samples (as saved by FAUST as folder names in the analysis map
#' and equivalent to \code{gs[[i]]@name} where i is an index in \code{gs}.) that
#' are to have their FAUST pops saved.
#' @param dir_save \code{character}. Specifies directory in which files are saved.
#' @param sample_name \code{character vector}. Specifies names of samples in order
#' in which they are found in \code{gs}.
#'
#' @return \code{invisible(TRUE)}.
.save_faust_pop <- function(project_path,
                            gs,
                            sample_name,
                            sel_sample,
                            pop,
                            dir_save){

  # extract data, filter and save
  for(sample in sel_sample){

    # get initial data
    fr <- flowWorkspace::gh_pop_get_data(gs[[which(sample_name == sample)]])
    ex <- flowCore::exprs(fr)

    # get filtered expression matrix
    ex <- .get_faust_pop(sample = sample,
                         project_path = project_path,
                         pop = pop,
                         ex = ex)

    # update flowFrame
    flowCore::exprs(fr) <- ex

    # save flowFrame
    flowCore::write.FCS(fr, file.path(dir_save, sample))
  }

  invisible(TRUE)
}

#' @title Filter an expression matrix to return only cells matching one or more FAUST annotation combinations
#'
#' @param sample \code{character}. Name of sample, as found in directory
#' <project_path>/faustData/sampleData/.
#' @param ex \code{matrix}. Matrix containing marker expression values for \code{sample}.
#' @inheritParams save_faust_pop
#'
#' @return Numeric matrix.
.get_faust_pop <- function(ex, sample, project_path, pop){
  # annotations for each cell for a given sample
  path_ann <- file.path(project_path, "faustData", "sampleData",
                        sample, "faustAnnotation.csv")
  if(!file.exists(path_ann)) stop(paste0("no faust annotation file found for sample", sample))
  faust_ann_tbl <- utils::read.table(file = path_ann,
                                     header = FALSE, sep = "`",
                                     stringsAsFactors = FALSE)[,1,drop = FALSE]

  # get expression matrix for cells in pop
  .get_faust_pop_pop(ex = ex, faust_ann = faust_ann_tbl, pop = pop)
}

#' Return expression matrix for cells that match a single FAUST label
#'
#' @inheritParams save_faust_pop
#' @param ex \code{numeric matrix}. Contains expression data.
#' @param faust_ann \code{character vector}. Specifies FAUST annotation for
#' each cell in \code{ex}.
#'
#' @return \code{numeric matrix} containing the expression levels for the markers specified.
.get_faust_pop_pop <- function(ex, faust_ann, pop){

  # vector to indicate if a match or not
  match <- rep(FALSE, nrow(ex))
  # if not a list, then make pop a list.
  # would not be a list if pop were specified as a character vector,
  # if only one set of markers can be a match.
  # would be a list if either of a set of markers could be a match.
  if(is.character(pop)) pop <- list(pop)
  for(i in seq_along(pop)){
    # set to TRUE if a match for given population
    match <- match | .is_faust_ann_a_match_pop(faust_ann = faust_ann, pop = pop)
  }

  ex[match, , drop = FALSE]
}

#' @title Check if FAUST annotation has a given level for a set of markers
#'
#' @inheritParams .get_faust_pop
#' @inheritParams .get_faust_pop_pop
#' @param pop \code{named character vector}. Names are marker names and values are
#' levels of marker. Values must be of the form or "<num_1>" or "<num_1>~<num_2>", where
#' <num_1> is the level for the marker and <num_2> is the total number of levels for
#' the marker.
#' @importFrom stringr str_locate str_sub
.is_faust_ann_a_match_pop <- function(faust_ann, pop){

  # vector to save if a match or not.
  # initialise to all TRUE, and then set to
  # FALSE if no match for a given marker
  match <- rep(TRUE, length(faust_ann))
  for(i in seq_along(pop)){
    match <- match & .is_faust_ann_a_match_for_marker(faust_ann = faust_ann,
                                                      marker = names(pop)[i],
                                                      level = pop[[i]])
  }
  match
}

#' @title Check if FAUST annotation has as given level for a given marker
#'
#' @inheritParams .get_faust_pop_pop
#' @param marker character. Name of marker.
#' @param level character. Level of marker, e.g. "1", "2" or "3".
.is_faust_ann_a_match_for_marker <- function(faust_ann, marker, level){
  faust_ann_level_loc_start <- stringr::str_locate(faust_ann[[1]], marker)[,"end"][[1]] + 2
  faust_ann_level <- stringr::str_sub(faust_ann[[1]], faust_ann_level_loc_start, faust_ann_level_loc_start + stringr::str_length(level) - 1)
  level == faust_ann_level
}
