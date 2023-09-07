#' @title Plot counts of FAUST pops by sample
#'
#' @description Saves plots of counts to \code{project_path/faustData/plotData/pop_stats}
#' for specified FAUST-identified subsets.
#'
#' @inheritParams save_faust_pop
#' @param pop \code{list} or \code{named character vector}. If a \code{character vector},
#' then counts of all subsets matching the set of marker levels are plotted.
#' If a \code{list}, then each list element must be a character vector. In that case,
#' one boxplot of frequencies is plotted for each element of the list of where the frequency is
#' the sum of all cells that match the population specified.
#' If \code{NULL}, then all subsets found are plotted.
#' @param breaks numeric vector. If supplied, then this is the breaks
#' for the x-axis. If not supplied, then the default breaks are used of
#' \code{c(0, 0.001, 0.005, 0.01, 0.02, 0.5, 0.1, 0.5, 1, 2, 5, 10, 20, 50, 100)}.
#' @param p_height \code{numeric}. Height of saved plot in cm. If \code{NULL}, then
#' it is set to \code{max(5, n_pop * 1.25)},
#' where \code{n_pop} is the number of sub-populations in \code{pop}.
#' @param p_width \code{numeric}. Width of saved plot in cm. If \code{NULL}, then set to 40.
#' @param trans_x \code{character} or \code{Transformation} object (from \code{scales} package).
#' the name of a transformation object or the object itself. Built-in
#' transformations include "asinh", "asn", "atanh", "boxcox", "date", "exp",
#' "hms", "identity", "log", "log10", "log1p", "log2", "logit", "modulus",
#' "probability", "probit", "pseudo_log", "reciprocal", "reverse", "sqrt" and
#' "time". Default is \code{asinh}.
#' @param exhaustive logical. If \code{TRUE}, then counts are taken from the exhaustive FAUST count matrix
#' rather than the count matrix after excluding subsets that don't appear in sufficiently many
#' experimental units. Default is \code{FALSE}.
#' @param limitsize logical. If \code{TRUE}, then plots of size 125cm x 125cm or more are
#' not saved but an error is returned when \code{cowplot::ggsave2} is run. Default
#' is \code{FALSE}.
#'
#' @details
#' The key parameter is \code{pop}. If \code{pop} is a list,
#' then each element is treated as a population for which counts
#' are required on aggregate (rather than for individual subsets
#' of that population). If \code{pop} is a character vector, then
#' it is treated as a population for which counts are required for
#' individual subsets.
#'
#' @import ggplot2
#' @import scales
#' @importFrom magrittr %>% %<>%
#' @examples
#' project_path <- usethis::proj_path(ext = "/inst/extdata")
#' # plot all subsets of pop, as it is a character vector
#' pop <- c(
#'   "CD3" = "+", "CD4" = "+", "CD8-IgD" = "-",
#'   "CD20" = "-", "CD33" = "-", "CD14" = "-",
#'   "TCRgd-CD19" = "-"
#' )
#' plot_faust_count(
#'   project_path = project_path,
#'   pop = pop
#' )
#' # plot counts of cells matching annotation of each
#' # element in pop, as pop is a list
#' pop <- list(
#'   c(
#'     "CD3" = "+", "CD4" = "+", "CD8-IgD" = "-",
#'     "CD20" = "-", "CD33" = "-", "CD14" = "-",
#'     "TCRgd-CD19" = "-"
#'   ),
#'   c(
#'     "CD3" = "+", "CD4" = "-", "CD8-IgD" = "+",
#'     "CD20" = "-", "CD33" = "-", "CD14" = "-",
#'     "TCRgd-CD19" = "-"
#'   )
#' )
#' plot_faust_count(
#'   project_path = project_path,
#'   pop = pop
#' )
#' @export
faust_pop_plot_count <- function(project_path,
                                 pop,
                                 font_size = 10,
                                 point_size_max = 2,
                                 breaks = NULL,
                                 trans_x = "asinh",
                                 p_width = NULL,
                                 p_height = NULL,
                                 exhaustive = FALSE,
                                 limitsize = FALSE) {
  # =======================================
  # Checks
  # =======================================

  # check that pop is either a character vector or a list of character vectors
  if (!is.null(pop)) {
    if (is.character(pop)) {
      if (is.null(names(pop))) stop("Character vectors in pop must be named. Names are the markers (e.g. CD4) and values are the FAUST annotations (e.g. '+' or '-').")
    }
    if (!is.character(pop)) {
      if (!is.list(pop)) stop("pop must be either a list of character vectors or a character vector, if not NULL.")
      if (!all(purrr::map_lgl(pop, function(x) is.character(x)))) {
        stop("pop must be either a list of character vectors or a character vector, if not NULL.")
      }
    }
  }


  # =======================================
  # Preparation
  # =======================================

  # break_vec
  if (is.null(breaks)) {
    breaks <- c(
      0, 0.001, 0.005, 0.01, 0.02, 0.5, 0.1, 0.5, 1,
      2, 5, 10, 20, 50, 100
    )
  }

  # base directory
  dir_faust <- file.path(project_path, "faustData")

  # ==============================
  # Subsetting
  # ==============================

  count_tbl_subset <- get_pop_counts(
    project_path = project_path,
    pop = pop,
    dem_col = c("sample", "exp_unit", "tot_count"),
    exhaustive = exhaustive
  )

  # =================================
  # Plot preparation
  # =================================

  # pivot longer and calculate prop and freq
  count_tbl_long <- count_tbl_subset %>%
    tidyr::pivot_longer(-c(sample:tot_count),
      names_to = "pop",
      values_to = "count"
    ) %>%
    dplyr::group_by(sample, tot_count) %>%
    dplyr::mutate(
      prop = count / tot_count,
      perc = prop * 100
    ) %>%
    dplyr::ungroup()

  # remove from sub-population annotation the markers

  # get the name of the overall population
  if (is.null(pop)) {
    title <- "Frequencies of all FAUST subsets"
  } else {
    title <- ifelse(is.list(pop), "Frequencies of subsets",
      paste0(.collapse_pop(pop = pop, search = FALSE),
        collapse = ""
      )
    )
  }

  # get the breaks
  breaks <- c(
    signif(min(count_tbl_long$perc), 2),
    breaks,
    signif(max(count_tbl_long$perc), 2)
  ) %>%
    unique()

  # get the transformation
  if (trans_x == "asinh") trans_x <- trans_asinh

  # =================================
  # Plot
  # =================================

  # plot
  p <- ggplot(count_tbl_long, aes(x = perc, y = pop)) +
    cowplot::theme_cowplot(font_size = font_size) +
    geom_boxplot(outlier.size = -1, size = 1) +
    geom_jitter(aes(col = pop, size = tot_count), alpha = 0.1) +
    guides(colour = "none") +
    labs(x = "Percentage (total classified cells)", y = "Cell population") +
    scale_x_continuous(
      trans = trans_x,
      breaks = function(lims) {
        breaks[breaks >= min(lims) & breaks <= max(lims)]
      }
    ) +
    theme(axis.text.x = element_text(angle = 90)) +
    cowplot::background_grid(major = "x", minor = "x") +
    scale_size_area(
      name = "Total classified cells",
      trans = trans_asinh,
      max_size = point_size_max
    ) +
    labs(title = title)

  # save plot
  dir_save <- file.path(dir_faust, "plotData", "pop_stats")
  if (exhaustive) dir_save <- paste0(dir_save, "-exhaustive")
  if (!dir.exists(dir_save)) dir.create(dir_save, recursive = TRUE)
  fn <- file.path(dir_save, paste0(title, ".png"))
  fn <- suppressWarnings(normalizePath(fn))
  if (file.exists(fn)) file.remove(fn)
  n_pop <- ncol(count_tbl_subset) - 3
  p_height <- ifelse(is.null(p_height), max(5, n_pop * 1.25), p_height)
  p_width <- ifelse(is.null(p_width), 40, p_width)
  cowplot::ggsave2(
    filename = fn,
    plot = p,
    height = p_height, width = p_width, units = "cm",
    limitsize = limitsize
  )

  invisible(TRUE)
}



#' @title Create inverse hyperbolic sin transformation object
trans_asinh <- scales::trans_new(
  "trans_asinh",
  function(x) asinh(x),
  function(x) sinh(x)
)

#' @title Get indices of columns that have a specified annotation
#'
#' @param data dataframe. Columns containg FAUST-population counts.
#' @param pop \code{named character vector}. Names specify marker and values specify level, e.g. c("CD4" = "-", "CD8" = "+").
#'
#' @return \code{Integer vector}.
.get_pop_match_ind <- function(data, pop) {
  pop_search_vec <- .collapse_pop(pop = pop, search = TRUE)
  purrr::map_lgl(colnames(data), function(pop_curr) {
    stringr::str_detect(pop_curr, pop_search_vec) %>% all()
  }) %>%
    which()
}
