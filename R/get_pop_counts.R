#' @title Get counts of subsets identified by FAUST
#'
#' @inheritParams save_faust_pop
#' @param data dataframe. Dataframe where columns specify counts. If \code{NULL},
#' then it is read in from \code{project_path/faustData/faustCountMatrix.rds}.
#' @param pop \code{named character vector} or \code{list}. If a \code{character vector},
#' then names specify markers and values specify levels, e.g. c("CD4" = "-", "CD8" = "+"). Counts of
#' all annotations will be returned that match these annotations separately, i.e. not summed.
#' If a \code{list}, then each element must be a character vector as above. For a given list element,
#' instead of all returning all subsets individually that have the correct level for the specified markers,
#' we sum over all the subsets matching the specified annotation and plot the final count.
#' If \code{NULL}, then all subsets are returned.
#' @param dem_col \code{character vector}. Specifies names of columns in \code{data} that we wish
#' to keep, regardless of if they match a FAUST annotation.
#' @param exhaustive logical. If \code{TRUE}, then counts are taken from the exhaustive FAUST count matrix
#' rather than the count matrix after excluding subsets that don't appear in sufficiently many
#' experimental units. Default is \code{FALSE}.
#'
#' @return A dataframe with columns as specified in \code{dem_col}, \code{tot_count} (total number of cells classified for a
#' sample), \code{pop} (population name) and
#' \code{count} (number of cells in population).
#'
#' @examples
#' # get counts for cells matching one annotation
#' pop <- c("CD4"  = "-", "CD8"  = "+")
#' get_pop_counts(pop = pop)
#' # get counts for cells matching either of the two annotations
#' pop <- list(c("CD4"  = "-", "CD8"  = "+"), c("CD8" = "-", "CD4" = "+"))
#' get_pop_counts(pop = pop)
#' @export
get_pop_counts <- function(project_path = NULL, data = NULL, pop = NULL,
                           dem_col = c('sample', 'exp_unit', 'tot_count', 'tot_count_classified',
                                       'sampleName', 'experimentalUnit'),
                           exhaustive = FALSE){

  # read in data
  if(is.null(data)){
    if(is.null(project_path)) stop("If data is not specified, then project_path must be.")
    # base directory
    dir_faust <- file.path(project_path, 'faustData')

    # raw data
    count_mat_name <- ifelse(exhaustive, 'exhaustiveFaustCountMatrix.rds',
                             'faustCountMatrix.rds')
    count_mat <- readRDS(file.path(dir_faust, count_mat_name))



    # analysis map
    analysis_map <- readRDS(file.path(dir_faust, 'metaData',
                                      'analysisMap.rds')) %>%
      tibble::as_tibble()

    # bind exp_unit column to count_tbl
    data <- count_mat %>%
      tibble::as_tibble() %>%
      dplyr::mutate(sampleName = rownames(count_mat)) %>%
      dplyr::select(sampleName, everything()) %>%
      dplyr::left_join(analysis_map %>%
                         dplyr::select(-impH),
                       by = 'sampleName') %>%
      dplyr::select(sampleName, experimentalUnit, everything()) %>%
      dplyr::rename(sample = sampleName, exp_unit = experimentalUnit)

    # calculate total classified cells
    count_tbl_tot <- data %>%
      tidyr::pivot_longer(-c(sample:exp_unit),
                          names_to = 'pop',
                          values_to = 'count') %>%
      dplyr::group_by(sample) %>%
      dplyr::summarise(tot_count = sum(count),
                       tot_count_classified = sum(count[pop != "0_0_0_0_0"]),
                       .groups = 'drop')

    data %<>%
      dplyr::left_join(count_tbl_tot,
                       by = 'sample')

    if(is.null(pop)) return(data %>% dplyr::select(sample, exp_unit, tot_count, everything()))
  }

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
