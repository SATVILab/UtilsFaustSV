#' @title Save FAUST subset as an FCS file
#'
#' @description
#' Save a specified FAUST-identified population to
#' \code{project_path/faustData/fcsData}
#' as an FCS file for all samples gated or just a specified subset.
#'
#' @param project_path character. FAUST project directory.
#' @param dir_save character.
#' Directory to save to.
#' If \code{NULL} (the default), then FCS files
#' are saved to
#' \code{file.path(project_path, "faustData", "fcsData", <pop_defn>)},
#' where \code{<pop_defn>} is a concatenation of the population
#' definition, e.g list(c("CD3" = 1, "CD4" = 2)) becomes "CD3~1~CD4~2~".
#' If \code{character}, then the FCS files are saved directly to this directory.
#' As stated before, default is \code{NULL}.
#' @param fr_source GatingSet or character vector. If a \code{GatingSet}, then
#' the flowFrames within the GatingSet are used to create the output fcs files.
#' If a character vector,
#' then must specify a directory containing FCS files. These FCS
#' files are used to create the output fcs files.
#' @param pop \code{list} or \code{named character vector}.
#' If a \code{character vector},
#' then all cells matching the set of marker levels are returned.
#' If a \code{list}, then each
#' element must be a \code{character vector}, and then
#' all cells matching either of
#' these specified sets of marker levels are returned.
#' @param sample integer vector or character vector.
#' If an integer vector, then it
#' specifies the indices of the samples for which to save output. If character,
#' then it specifies the names of the sames to save output for. If \code{NULL},
#' then the output for every sample is saved. Default is \code{NULL}.
#' @param trans_fn function. If supplied, this function
#' is applied to the expression data.
#' Useful for back-transformation.
#' If \code{NULL}, then no transformation is applied. Default is \code{NULL}.
#' @param trans_chnl character vector.
#' If specified, \code{trans_fn} is applied to only these channels.
#' If \code{NULL} and if \code{trans_fn} is not \code{NULL},
#' then \code{trans_fn} is applied to entire
#' expression matrix. Default is \code{NULL}.
#'
#' @return \code{invisible(TRUE)}. Side effect is the saved FCS file.
#' @examples faust_fcs_write(
#'   project_path = "", pop = list("CD3" = 2),
#'   gs = gs, sample = 1
#' )
#'
#' @export
faust_fcs_write <- function(project_path,
                            dir_save = NULL,
                            pop,
                            fr_source = NULL,
                            sample = NULL,
                            trans_fn = NULL,
                            trans_chnl = NULL) {
  # =============================
  # Check
  # =============================

  if (any(c(missing(project_path), missing(fr_source), missing(pop)))) {
    stop("Both of project_path and pop parameters must have arguments.")
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
  # get vector of sample names in gs
  # if(!is.null(gs)){
  #  active_sample_vec <- sapply(seq_along(gs), function(i) gs[[i]]@name)# readRDS(paste0(project_path, "/faustData/metaData/analysisMap.rds"))[,"sampleName"]
  #  sel_sample_vec <- switch(typeof(sample),
  #                           "integer" = ,
  #                           "double" = active_sample_vec[sample],
  #                           "NULL" = active_sample_vec,
  #                           "character" = active_sample_vec[vapply(active_sample_vec, function(x) x %in% sample,
  #                                                                  logical(length(sample)))],
  #                           stop("Incorrect specification of sample parameter."))
  #  # get vector of sample names in gs
  #  sample_name_vec <- active_sample_vec
  # } else{
  active_sample_vec <- list.dirs(
    file.path(
      project_path, "faustData",
      "sampleData"
    ),
    recursive = FALSE, full.names = FALSE
  )
  sel_sample_vec <- switch(typeof(sample),
    "integer" = ,
    "double" = active_sample_vec[sample],
    "NULL" = active_sample_vec,
    "character" = active_sample_vec[vapply(
      active_sample_vec, function(x) x %in% sample,
      logical(length(sample))
    )],
    stop("Incorrect specification of sample parameter.")
  )
  # get vector of sample names in GatingSet
  if (class(fr_source) == "GatingSet") {
    sample_name_vec <- flowWorkspace::sampleNames(fr_source)
  } else if (class(fr_source) == "character") {
    sample_name_vec <- lapply(fr_source, function(fr_source_curr) {
      list.files(fr_source_curr, pattern = "fcs$", full.names = FALSE)
    }) %>%
      unlist()
  }




  # create directory to save to and pop list to compare to
  # --------------------------------------

  # get concatenated name of population
  if (!is.list(pop)) {
    if (is.null(names(pop))) {
      pop_name <- pop
    } else {
      pop_name <- ""
      for (i in seq_along(pop)) {
        pop_name <- paste0(
          pop_name,
          names(pop)[i],
          "~",
          pop[[i]],
          "~"
        )
      }
    }
  } else if (
    is.list(pop) &&
      all(purrr::map_lgl(pop, function(x) is.character(x) || is.numeric(x)))
  ) {
    pop_name_vec <- purrr::map(pop, function(pop_curr) names(pop_curr)) %>%
      unlist() %>%
      unique()

    pop_level_list <- purrr::map(pop_name_vec, function(pop_name) {
      pop_level_vec <- purrr::map_chr(pop, function(pop_curr) {
        if (!pop_name %in% names(pop_curr)) {
          return("_")
        }
        stringr::str_sub(as.character(pop_curr[[pop_name]]), end = 1)
      })
      pop_level_vec
    }) %>%
      stats::setNames(pop_name_vec)

    pop_name <- purrr::map_chr(seq_along(pop_level_list), function(i) {
      paste0(
        names(pop_level_list)[i],
        "~",
        paste0(unique(pop_level_list[[i]]),
          collapse = ""
        ),
        "~"
      )
    }) %>%
      paste0(collapse = "")
  }


  #  create directory
  if (is.null(dir_save)) {
    dir_save <- file.path(
      project_path,
      "faustData",
      "fcsData",
      pop_name
    )
  }
  if (!dir.exists(dir_save)) dir.create(dir_save, recursive = TRUE)

  # =============================
  # Save FAUST pops
  # =============================

  .faust_fcs_write(
    project_path = project_path, fr_source = fr_source,
    sample_name = sample_name_vec,
    sel_sample = sel_sample_vec,
    pop = pop,
    dir_save = dir_save,
    trans_fn = trans_fn,
    trans_chnl = trans_chnl
  )
}

#' @title Save FAUST subset as an FCS file
#'
#' @inheritParams faust_fcs_write
#' @param sel_sample \code{character vector}. Character vector specifying
#' the names of the samples
#' (as saved by FAUST as folder names in the analysis map
#' and equivalent to \code{gs[[i]]@name}
#' where i is an index in \code{gs}.) that
#' are to have their FAUST pops saved.
#' @param dir_save \code{character}.
#' Specifies directory in which files are saved.
#' @param sample_name \code{character vector}.
#' Specifies names of samples in order
#' in which they are found in \code{gs}.
#'
#' @return \code{invisible(TRUE)}.
.faust_fcs_write <- function(project_path,
                             fr_source,
                             sample_name,
                             sel_sample,
                             pop,
                             dir_save,
                             trans_fn = NULL,
                             trans_chnl = NULL) {
  # extract data, filter and save
  for (sample in sel_sample) {
    # get initial data
    # if(!is.null(gs)){
    if (class(fr_source) == "GatingSet") {
      fr <- try(flowWorkspace::gh_pop_get_data(
        fr_source[[which(sample_name == sample)]]
      ))
    } else if (class(fr_source) == "character") {
      fr <- purrr::map(fr_source, function(fr_source_curr) {
        fr_path <- file.path(fr_source_curr, sample)
        if (!file.exists(fr_path)) {
          return(NULL)
        }
        flowCore::read.FCS(fr_path)
      }) %>%
        purrr::compact()
      fr <- fr[[1]]
    }

    if (class(fr) == "try-error") {
      print(sample_name)
      print(sample)
      stop(paste0("error in loading sample", sample))
    }


    ex <- flowCore::exprs(fr)


    # annotations for each cell for a given sample
    path_ann <- file.path(
      project_path, "faustData", "sampleData",
      sample, "faustAnnotation.csv"
    )
    if (!file.exists(path_ann)) {
      stop(paste0("no faust annotation file found for sample ", sample))
    }
    faust_ann_tbl <- utils::read.table(
      file = path_ann,
      header = FALSE, sep = "`",
      stringsAsFactors = FALSE
    )[, 1, drop = FALSE]

    if (nrow(ex) != nrow(faust_ann_tbl)) {
      stop(paste0(
        "ex in FlowFrame has different number of rows to the number of faust annotations for sample", # nolint
        sample
      ))
    }

    # get filtered expression matrix
    ex <- .get_faust_pop(
      pop = pop,
      ex = ex,
      faust_ann = faust_ann_tbl
    )

    # transform
    if (!is.null(trans_fn)) {
      if (is.null(trans_chnl)) {
        ex <- trans_fn(ex)
      } else {
        for (nm in trans_chnl) {
          ex[, nm] <- trans_fn(ex[, nm])
        }
      }
    }

    # update flowFrame
    if (class(fr) == "cytoframe") {
      fr <- flowWorkspace::cytoframe_to_flowFrame(fr)
    }
    flowCore::exprs(fr) <- ex

    # save flowFrame
    flowCore::write.FCS(fr, file.path(dir_save, sample))
    TRUE
  }

  invisible(TRUE)
}

#' @title Filter an expression matrix to
#' return only cells matching one or more FAUST annotation combinations
#'
#' @param sample \code{character}. Name of sample, as found in directory
#' <project_path>/faustData/sampleData/.
#' @param ex \code{matrix}.
#' Matrix containing marker expression values for \code{sample}.
#' @inheritParams faust_fcs_write
#'
#' @details
#' This is effectively an inclusive OR statement across the different population
#' annotations.
#'
#' @return Numeric matrix.
.get_faust_pop <- function(ex, pop, faust_ann) {
  # vector to indicate if a match or not
  match <- rep(FALSE, nrow(ex))
  # if not a list, then make pop a list.
  # would not be a list if pop were specified as a character vector,
  # if only one set of markers can be a match.
  # would be a list if either of a set of markers could be a match.
  if (is.character(pop) || is.numeric(pop)) pop <- list(pop)
  for (i in seq_along(pop)) {
    # set to TRUE if a match for given population
    match <- match |
      .is_faust_ann_a_match_pop(faust_ann = faust_ann, pop = pop[[i]])
  }

  if (sum(match) == 0) {
    ex <- ex[1, , drop = FALSE]
    for (j in seq_len(ncol(ex))) ex[1, j] <- NA
    return(ex)
  }

  ex[match, , drop = FALSE]
}

#' @title Check if FAUST annotation has a given level for a set of markers
#'
#' @inheritParams .get_faust_pop
#' @param pop \code{named character vector}.
#' Names are marker names and values are
#' levels of marker.
#' Values must be of the form or "<num_1>" or "<num_1>~<num_2>",
#' where <num_1> is the level for the marker and
#' <num_2> is the total number of levels for
#' the marker.
#'
#' @details
#' This is effectively an AND statement across
#' all markers specified in the
#' single FAUST population annotation.
#'
#' @importFrom stringr str_locate str_sub
.is_faust_ann_a_match_pop <- function(faust_ann, pop) {
  # vector to save if a match or not.
  # initialise to all TRUE, and then set to
  # FALSE if no match for a given marker
  if (is.null(names(pop))) {
    return(faust_ann[[1]] == pop)
  }
  match <- rep(TRUE, length(faust_ann))
  for (i in seq_along(pop)) {
    match <- match & .is_faust_ann_a_match_for_marker(
      faust_ann = faust_ann,
      marker = names(pop)[i],
      level = pop[[i]]
    )
  }
  match
}

#' @title Check if FAUST annotation has as given level for a given marker
#'
#' @inheritParams .get_faust_pop
#' @param marker character. Name of marker.
#' @param level character. Level of marker, e.g. "1", "2" or "3".
.is_faust_ann_a_match_for_marker <- function(faust_ann,
                                             marker,
                                             level) {
  k <- 1
  while (faust_ann[[1]][k] == "0_0_0_0_0") {
    if (k >= nrow(faust_ann)) {
      return(rep(FALSE, nrow(faust_ann)))
    }
    k <- k + 1
  }
  typical_cluster_annotation <- faust_ann[[1]][k]
  faust_ann_level_loc_start <- stringr::str_locate(
    typical_cluster_annotation, marker
  )[, "end"][[1]] + 2
  faust_ann_level <- stringr::str_sub(
    faust_ann[[1]],
    faust_ann_level_loc_start,
    faust_ann_level_loc_start + nchar(level) - 1
  )
  level == faust_ann_level
}
