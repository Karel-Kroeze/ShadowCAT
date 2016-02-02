#' Maximum Information item selection
#' 
#' Naive item selection based on maximum information only.
#' 
#' Selects the item with the highest information.
#' @param available vector with indeces of available items
#' @param item_information vector with information of each yet available item, with zeros for administered items (returned by get_item_information() with pad = TRUE)
#' @return item index of item with maximum information.
#' @export
get_item_index_max_information <- function(available, item_information) {
  which(item_information == max(item_information[available]))
}


