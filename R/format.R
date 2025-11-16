.faust_pop_format_lvl <- function(pop,
                                  to = "tilde",
                                  project_path = NULL,
                                  marker_to_n_lvl = NULL) {
  # check
  if (!is.character(pop)) {
    stop(paste0(
      "The pop argument must be a character vector, not '", class(pop), "'"
    ))
  }
  if (!to %in% c("tilde", "common")) {
    stop(paste0(
      "The to argument must be one of 'tilde' or 'common', not '", to, "'"
    ))
  }
  switch(to,
    "tilde" = .faust_pop_format_level_to_tilde(
      pop = pop, marker_to_n_lvl = marker_to_n_lvl, project_path = project_path
    ),
    "common" = .faust_pop_format_level_to_common(pop = pop),
    stop(paste0(
      "to argument must be one of 'tilde' or 'common', not '", to, "'"
    ))
  )
}

.faust_pop_format_level_to_tilde <- function(pop,
                                             marker_to_n_lvl = NULL,
                                             project_path = NULL) {
  if (grepl("~\\d~\\d~", pop)) {
    return(pop)
  }
  # get number of levels for the marker
  if (is.null(marker_to_n_lvl)) {
    # if not provided, use FAUST directory
    if (is.null(project_path)) {
      stop(paste0(
        "If marker_to_n_lvl is not provided, then project_path must be provided." # nolint
      ))
    }
    if (!dir.exists(project_path)) {
      stop(paste0(
        "The directory '", project_path, "' does not exist."
      ))
    }
    marker_to_n_lvl <- try(
      faust_marker_get_usage(project_path = project_path),
      silent = TRUE
    )
    if (inherits(marker_to_n_lvl, "try-error")) {
      stop(paste0(
        "Could not get levels from markers in the project directory ", project_path # nolint
      ))
    }
  } else {
    # if provided, check that it's fine
    if (is.null(names(marker_to_n_lvl))) {
      stop(paste0(
        "The marker_to_n_lvl argument must have names"
      ))
    }
  }

  # actually convert
  # ------------------------
  # get marker positions
  marker_to_pos <- lapply(names(marker_to_n_lvl), function(marker) {
    .faust_chr_pos_get(string = pop, pattern = marker)
  }) |>
    stats::setNames(names(marker_to_n_lvl))
  marker_to_pos_start <- sapply(marker_to_pos, function(x) x[1]) |>
    stats::setNames(names(marker_to_n_lvl))
  marker_to_pos_end <- sapply(marker_to_pos, function(x) x[2]) |>
    stats::setNames(names(marker_to_n_lvl))
  marker_vec_ordered <- names(marker_to_pos)[order(unlist(marker_to_pos_start))]
  marker_to_pos_start <- marker_to_pos_start[marker_vec_ordered]
  marker_to_pos_end <- marker_to_pos_end[marker_vec_ordered]
  marker_to_pos_lvl_start <- (marker_to_pos_end + 1) |>
    stats::setNames(marker_vec_ordered)
  marker_to_pos_lvl_end <- c(marker_to_pos_start[-1] - 1, nchar(pop)) |>
    stats::setNames(marker_vec_ordered)
  # order markers by position
  pop_rep <- NULL
  for (i in seq_along(marker_vec_ordered)) {
    marker <- marker_vec_ordered[i]
    lvl <- .faust_marker_format_level_to_tilde(
      marker = marker,
      lvl = substr(
        pop,
        start = marker_to_pos_lvl_start[[marker]],
        stop = marker_to_pos_lvl_end[[marker]]
      ),
      n_lvl = marker_to_n_lvl[[marker]]
    )
    pop_rep <- paste0(pop_rep, marker, lvl)
  }

  # return
  pop_rep
}

.faust_chr_pos_get <- function(string, pattern) {
  # Find the start and end points of the subset in the string
  start_point <- regexpr(pattern, string)
  end_point <- start_point +
    attr(regexpr(pattern, string), "match.length") - 1
  c(start_point, end_point)
}

.faust_pop_format_level_to_common <- function(pop) {
  if (!grepl("~\\d~\\d~", pop)) {
    return(pop)
  }
  for (i in seq_along(tilde_to_common)) {
    pop <- gsub(names(tilde_to_common)[i], tilde_to_common[[i]], pop)
  }
  pop
}

.faust_marker_format_level <- function(marker, lvl,
                                       to = "tilde",
                                       n_lvl = NULL,
                                       project_path = NULL) {
  if (!to %in% c("tilde", "common")) {
    stop(paste0(
      "The to argument must be one of 'tilde' or 'common', not '", to, "'"
    ))
  }
  switch(to,
    "tilde" = .faust_marker_format_level_to_tilde(
      marker = marker, lvl = lvl, n_lvl = n_lvl, project_path = project_path
    ),
    "common" = .faust_marker_format_level_to_common(lvl = lvl),
    stop(paste0(
      "to argument must be one of 'tilde' or 'common', not '", to, "'"
    ))
  )
}

format_vec_valid_common <- c(
  "-", "+", "Dim", "Bright", "MedLow", "MedHigh", "VeryBright"
)

format_vec_valid_tilde <- c(
  "~1~2~", "~2~2~",
  "~1~3~", "~2~3~", "~3~3~",
  "~1~4~", "~2~4~", "~3~4~", "~4~4~"
)
tilde_to_common <- c(
  "~1~2~" = "-", "~2~2~" = "+",
  "~1~3~" = "-", "~2~3~" = "Dim", "~3~3~" = "Bright",
  "~1~4~" = "-", "~2~4~" = "MedLow", "~3~4~" = "MedHigh", "~4~4~" = "VeryBright"
)

.faust_marker_format_level_to_tilde <- function(marker,
                                                lvl,
                                                n_lvl = NULL,
                                                project_path = NULL) {
  # address easy case where it's tilde already
  if (grepl("^~\\d~\\d~$", lvl)) {
    if (!(lvl %in% format_vec_valid_tilde)) {
      stop(paste0(
        "The level '", lvl, "' for marker '", marker, "' is not valid."
      ))
    }
    return(lvl)
  }
  # actually convert
  # ------------------------

  # check
  if (!lvl %in% format_vec_valid_common) {
    stop(paste0(
      "The level '", lvl, "' for marker '",
      marker, "' is not valid."
    ))
  }

  # get number of levels for the marker
  if (is.null(n_lvl)) {
    if (missing(marker)) {
      stop(paste0(
        "If n_lvl is not provided, then marker must be provided."
      ))
    }
    # if not provided, use FAUST directory
    if (is.null(project_path)) {
      stop(paste0(
        "If n_lvl is not provided, then project_path must be provided."
      ))
    }
    if (!dir.exists(project_path)) {
      stop(paste0(
        "The directory '", project_path, "' does not exist."
      ))
    }
    n_lvl <- try(
      faust_marker_get_usage(project_path = project_path)[[marker]],
      silent = TRUE
    )
    if (inherits(n_lvl, "try-error")) {
      stop(paste0(
        "Marker ", marker, " not found in the project directory ", project_path
      ))
    }
    if (is.null(n_lvl)) {
      stop(paste0(
        "Marker ", marker, " not found in the project."
      ))
    }
  } else {
    # if provided, check that it's fine
    n_lvl <- try(as.numeric(n_lvl), silent = TRUE)
    if (!is.numeric(n_lvl)) {
      stop(paste0(
        "The n_lvl argument must be numeric, not '", n_lvl, "'"
      ))
    }
  }

  # actually convert
  common_to_tilde_vec <- switch(as.character(n_lvl),
    "2" = c("-" = "~1~2~", "+" = "~2~2~"),
    "3" = c(
      "-" = "~1~3~", "Dim" = "~2~3~", "Bright" = "~3~3~"
    ),
    "4" = c(
      "-" = "~1~4~", "MedLow" = "~2~4~",
      "MedHigh" = "~3~4~", "VeryBright" = "~4~4~"
    )
  )
  # check it's present
  if (!lvl %in% names(common_to_tilde_vec)) {
    stop(paste0(
      "The level '", lvl, "' for marker '", marker,
      "' is not valid due to the number of levels for marker '", marker,
      "' being '", n_lvl, "'."
    ))
  }
  # return
  common_to_tilde_vec[[lvl]]
}

.faust_marker_format_level_to_common <- function(lvl) {
  # address case where it's common already
  if (!grepl("^~\\d~\\d~$", lvl)) {
    if (!lvl %in% format_vec_valid_common) {
      stop(paste0(
        "The level '", lvl, "' is not valid."
      ))
    }
    return(lvl)
  }
  tilde_to_common[[lvl]]
}
