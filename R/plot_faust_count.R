#' @title Plot counts of FAUST pops by sample
#'
#' @inheritParams save_faust_pop
#' @param pop \code{list} or \code{named character vector}. If a \code{character vector},
#' then counts of all subsets matching the set of marker levels are plotted.
#' If a \code{list}, then each list element must be a character vector. In that case,
#' one boxplot of frequencies is plotted for each element of the list of where the frequency is
#' the sum of all cells that match the population specified.
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
#' @import ggplot2
#' @import scales
#' @importFrom magrittr %>% %<>%
#' @examples
#' project_path <- usethis::proj_path(ext = "/inst/extdata")
#' pop <- c("CD3" = "+", "CD4" = "+", "CD8-IgD" = "-",
#'          "CD20" = "-", "CD33" = "-", "CD14" = "-",
#'          "TCRgd-CD19" = "-")
#' plot_faust_count(project_path = project_path,
#'                  pop = pop)
#' pop <- list(c("CD3" = "+", "CD4" = "+", "CD8-IgD" = "-",
#' "CD20" = "-", "CD33" = "-", "CD14" = "-",
#' "TCRgd-CD19" = "-"),
#' c("CD3" = "+", "CD4" = "-", "CD8-IgD" = "+",
#'   "CD20" = "-", "CD33" = "-", "CD14" = "-",
#'  "TCRgd-CD19" = "-"))
#' plot_faust_count(project_path = project_path,
#' pop = pop)
plot_faust_count <- function(project_path,
                             pop,
                             font_size = 10,
                             point_size_max = 2,
                             breaks = NULL,
                             trans_x = 'asinh',
                             p_width = NULL,
                             p_height = NULL){
  # =======================================
  # Checks
  # =======================================

  # check that pop is either a character vector or a list of character vectors
  if(!is.character(pop)){
    if(!is.list(pop)) stop('pop must be either a list of character vectors or a character vector.')
    if(!all(purrr::map_lgl(pop, function(x) is.character(x)))){
      stop('pop must be either a list of character vectors or a character vector.')
    }
  }

  # =======================================
  # Preparation
  # =======================================

  # break_vec
  if(is.null(breaks)) breaks <- c(0, 0.001, 0.005, 0.01, 0.02, 0.5, 0.1, 0.5, 1,
                                  2, 5, 10, 20, 50, 100)

  # base directory
  dir_faust <- file.path(project_path, 'faustData')

  # raw data
  count_mat <- readRDS(file.path(dir_faust, 'faustCountMatrix.rds'))
  analysis_map <- readRDS(file.path(dir_faust, 'metaData',
                                    'analysisMap.rds')) %>%
    tibble::as_tibble()

  # bind exp_unit column to count_tbl
  count_tbl <- count_mat %>%
    tibble::as_tibble() %>%
    dplyr::mutate(sampleName = rownames(count_mat)) %>%
    dplyr::select(sampleName, everything()) %>%
    dplyr::left_join(analysis_map %>%
                       dplyr::select(-impH),
                     by = 'sampleName') %>%
    dplyr::select(sampleName, experimentalUnit, everything()) %>%
    dplyr::rename(sample = sampleName, exp_unit = experimentalUnit)

  # calculate total classified cells
  count_tbl_tot <- count_tbl %>%
    tidyr::pivot_longer(-c(sample:exp_unit),
                 names_to = 'pop',
                 values_to = 'count') %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(tot_count = sum(count),
                     .groups = 'drop')

  count_tbl %<>%
    dplyr::left_join(count_tbl_tot,
                     by = 'sample')

  # ==============================
  # Subsetting
  # ==============================

  count_tbl_subset <- .get_pop_counts(data = count_tbl,
                                      pop = pop,
                                      dem_col = c('sample', 'exp_unit', 'tot_count'))

  # =================================
  # Plot preparation
  # =================================

  # pivot longer and calculate prop and freq
  count_tbl_long <- count_tbl_subset %>%
    tidyr::pivot_longer(-c(sample:tot_count),
                        names_to = "pop",
                        values_to = "count") %>%
    dplyr::group_by(sample, tot_count) %>%
    dplyr::mutate(prop = count/tot_count,
                  perc = prop * 100) %>%
    dplyr::ungroup()

  # remove from sub-population annotation the markers

  # get the name of the overall population
  title <- ifelse(is.list(pop), "Total counts by subset",
                  paste0(.collapse_pop(pop = pop, search = FALSE),
                         collapse = ""))
  # get the breaks
  breaks <- c(signif(min(count_tbl_long$perc), 2),
              breaks,
              signif(max(count_tbl_long$perc), 2)) %>%
    unique

  # get the transformation
  if(trans_x == 'asinh') trans_x <- trans_asinh

  # =================================
  # Plot
  # =================================

  # plot
  p <- ggplot(count_tbl_long, aes(x = perc, y = pop)) +
    cowplot::theme_cowplot(font_size = font_size) +
    geom_boxplot(outlier.size = -1, size = 1) +
    geom_jitter(aes(col = pop, size = tot_count), alpha = 0.1) +
    guides(colour = 'none') +
    labs(x = "Percentage (total classified cells)", y = "Cell population") +
    scale_x_continuous(trans = trans_x,
                       breaks = function(lims){
                         breaks[breaks >= min(lims) & breaks <= max(lims)]
                       }) +
    theme(axis.text.x = element_text(angle = 90)) +
    cowplot::background_grid(major = 'x', minor = 'x') +
    scale_size_area(name = 'Total classified cells',
                    trans = trans_asinh,
                    max_size = point_size_max) +
    labs(title = title)

  # save plot
  if(save){
    dir_save <- file.path(dir_faust, 'plotData', 'pop_stats')
    if(!dir.exists(dir_save)) dir.create(dir_save, recursive = TRUE)
    n_pop <- ncol(count_tbl_subset) - 3
    p_height <- ifelse(is.null(p_height), max(5, n_pop * 1.25), p_height)
    p_width <- ifelse(is.null(p_width), 40, p_width)
    cowplot::ggsave2(filename = file.path(dir_save, paste0(title, ".png")),
                     plot = p,
                     height = p_height, width = p_width, units = 'cm')
  }
}

#' @title Get a searchable version of the marker levels
#'
#' @description
#'
#' @inheritParams plot_faust_count
#' @param search logical. If \code{TRUE}, then square brackets are concatenated
#' around each of the elements in the level strings.
#'
#' @return
#'
#' @examples
#' .collapse_pop(c("CD4" = "-", CD3" = "-"), search = FALSE); .collapse_pop(c("CD4" = "-", CD3" = "-"), search = TRUE)
.collapse_pop <- function(pop, search = FALSE){
  purrr::map_chr(seq_along(pop), function(i){
    if(!search) return(paste0(names(pop)[i], pop[[i]]))
    level <- purrr::map_chr(1:stringr::str_length(pop[[i]]), function(j){
      paste0("[[", stringr::str_sub(pop[[i]], j, j), "]]")
    })
    paste0(names(pop)[i], level)
  })
}


#' @title Create inverse hyperbolic sin transformation object
trans_asinh <- scales::trans_new('trans_asinh',
                                 function(x) asinh(x),
                                 function(x) sinh(x))

#' @title Get indices of columns that have a specified annotation
#'
#' @param data dataframe. Columns containg FAUST-population counts.
#' @param pop \code{named character vector}. Names specify marker and values specify level, e.g. c("CD4" = "-", "CD8" = "+").
#'
#' @return \code{Integer vector}.
.get_pop_match_ind <- function(data, pop){
  pop_search_vec <- .collapse_pop(pop = pop, search = TRUE)
  purrr::map_lgl(colnames(data), function(pop_curr){
    stringr::str_detect(pop_curr, pop_search_vec) %>% all
  }) %>%
    which
}

#' @title Get counts of subsets identified by FAUST
#'
#' @param data dataframe. Dataframe where columns specify counts.
#' @param pop \code{named character vector} or \code{list}. If a \code{character vector},
#' then names specify markers and values specify levels, e.g. c("CD4" = "-", "CD8" = "+") means
#' that we are interested in FAUST subsets annotated "CD4-CD8+". Counts of all subsets with these
#' levels for these markers are plotted separately.
#' If a \code{list}, then each element must be a character vector as above. For a given list element,
#' instead of all plotting all subsets individually that have the correct level for the specified markers,
#' we sum over all the subsets matching the specified annotation and plot the final count.
#' @param dem_col \code{character vector}. Specifies names of columns in \code{data} that we wish
#' to keep, regardless of if they match a FAUST annotation.
#'
#' @return A dataframe with columns as specified in \code{dem_col}, \code{pop} (population name) and
#' \code{count} (number of cells in population).
#'
#' @examples
#' .get_pop_counts(data = count_tbl, pop = c("CD4"  = "-", "CD8"  = "+"))
#' .get_pop_counts(data = count_tbl, pop = list(c("CD4"  = "-", "CD8"  = "+"), c("CD8" = "-", 'CD4" = "+")))
.get_pop_counts <- function(data, pop, dem_col = c('sample', 'exp_unit', 'tot_count')){

  # get columns whose indices match the "demographic" info
  dem_col_ind_vec <- which(colnames(data) %in% dem_col)
  data_dem <- data[, dem_col_ind_vec, drop = FALSE]
  # if a single annotation set is all that's specified, return all subsets matching
  # it individually
  if(is.character(pop)){
    pop_col_ind_vec <- .get_pop_match_ind(data = data, pop = pop)
    pop_col_name_vec <- colnames(data[,pop_col_ind_vec])
    # for each column identified, remove each annotation specifying
    # the main subset of which it's a part
    for(i in seq_along(pop_col_name_vec)){
      for(j in seq_along(pop)){
        pop_col_name_vec[i] <- stringr::str_remove(pop_col_name_vec[i],
                                                   paste0(names(pop)[[j]],
                                                          "[[",
                                                          pop[j],
                                                          "]]"))
      }
    }
    # select only columns found to match annotation set
    data_resp <- data[,pop_col_ind_vec, drop = FALSE]
    # rename columns as per above
    colnames(data_resp) <- pop_col_name_vec
    # join dem columns to response columns
    data_out <- dplyr::bind_cols(data_dem, data_resp)
    return(data_out)

  } else if(is.list(pop)){
    # if multiple pops are specified, sum within each
    data_out <- data_dem
    for(pop_curr in pop){
      # get columns matching current annotation set
      pop_col_ind_vec <- .get_pop_match_ind(data = data, pop = pop_curr)

      # add up counts over all population subsets identified
      count_vec <- rep(0, nrow(data_dem))
      for(pop_ind in pop_col_ind_vec){
        count_vec <- count_vec + data[[pop_ind]]
      }

      # append counts to dataframe
      pop_name <- paste0(.collapse_pop(pop = pop_curr, search = FALSE),
                         collapse = "")
      data_out[[pop_name]] <- count_vec
    }
    return(data_out)
  }
}
