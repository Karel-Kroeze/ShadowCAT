#' Maximum Information item selection
#' 
#' Naive item selection based on maximum information only.
#' 
#' Selects the item with the highest information.
#' @param available vector with indices of available items
#' @param item_information vector with summarized information of each yet available item, with zeros for administered items (returned by get_summarized_information() with pad = TRUE)
#' @param estimate Vector containing the theta estimate, with its covariance matrix as an attribute
#' @param stop_test rule for when to stop providing new items to patient; should be a list of the form
#' list(target = ..., max_n = ..., min_n = ..., cutoffs = ...), 
#' where max_n = test length at which testing should stop (even if target has not been reached yet in case of variance stopping rule), 
#' target = vector of maximum acceptable variances per dimension; NULL means no variance target,
#' min_n = minimum test length; NULL means no mimimum test length,
#' cutoffs = matrix containing cut off values per dimension (columns) and test iteration (rows). First row contains cut off values for when no items have been
#' administered yet, second row for when one item has been administered, etc. If estimate + 3SE < cutoff for each dimension at certain iteration, test stops; 
#' NULL means no cut off values
#' @param alpha Matrix of alpha parameters, one column per dimension, one row per item. Row names should contain the item keys. Note that so called within-dimensional models still use an alpha matrix, they simply 
#' have only one non-zero loading per item.
#' @param number_answers The number of answers given to the test (length of administered)
#' @return item Index of item with maximum information. For multi dimensional models, items that only load on dimensions for
#' which the variance target has already been reached, will not be returned
#' @export
get_item_index_max_information <- function(available, item_information, estimate, stop_test, alpha, number_answers) {
  result <- function() {
    uncompleted_dimensions <- get_uncompleted_dimensions()
    useful_item_indices <- get_useful_item_indices(uncompleted_dimensions)  
    available_and_useful <- get_available_and_useful_items(useful_item_indices)
    which(item_information == max(item_information[available_and_useful]))
  }
  
  get_uncompleted_dimensions <- function() {
    if (length(estimate) == 1 || is.null(stop_test$target))
      return(NULL)
    which(diag(attr(estimate, "variance")) >= stop_test$target)
  }
  
  get_useful_item_indices <- function(uncompleted_dimensions) {
    if (is.null(uncompleted_dimensions))
      return(1:nrow(alpha))
    loadings_larger_than_zero <- t(apply(alpha, 1, function(x) abs(x) > 1e-5 ))
    which(apply(loadings_larger_than_zero, 1, function(single_item_loadings_larger_than_zero) any(which(single_item_loadings_larger_than_zero) %in% uncompleted_dimensions)))
  }
  
  get_available_and_useful_items <- function(useful_item_indices) {
    available_and_useful <- intersect(useful_item_indices, available)
    if (length(available_and_useful) == 0)
      available
    else
      available_and_useful 
  }
  
  validate <- function() {
    if (is.null(attr(estimate, "variance")))
      add_error("estimate", "should have its variance matrix as an attribute")
  }
  
  invalid_result <- function() {
    list(errors = errors())
  }
  
  validate_and_run()
}


