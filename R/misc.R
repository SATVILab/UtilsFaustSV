#' @title Get a searchable version of the marker levels
#'
#' @description
#' 
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

    # return the entire string if it's actually letters and not alphanumeric
    if(stringr::str_sub(pop[[i]], 1, 1) %in% c(letters, LETTERS)) return(paste0(names(pop)[i], pop[[i]]))

    level <- purrr::map_chr(1:stringr::str_length(pop[[i]]), function(j){
      paste0("[[", stringr::str_sub(pop[[i]], j, j), "]]")
    })
    paste0(names(pop)[i], level)
  })
}
