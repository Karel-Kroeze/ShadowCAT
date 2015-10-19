#' Maximum Information item selection
#' 
#' Naive item selection based on maximum information only.
#' 
#' Selects the item with the highest value of a specified objective function.
#' @param test test object, argument is not used and can be removed
#' @param person person object
#' @param item_information vector with information of each yet available item, with zeros for administered items (returned by objective() with pad = TRUE)
#' @return item index of item with the highest value of the objective function.
#' @export
MI <- function(test, person, item_information) {
  # fetch highest information item
  which(item_information == max(item_information[person$available]))
}


