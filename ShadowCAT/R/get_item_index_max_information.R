#' Item index maximum information
#' 
#' Get the index of the item with maximum information from the available items. 
#' For multi dimensional models, items that only load on dimensions for
#' which the variance target has already been reached, will not be returned.
#' 
#' @param available Vector with indices of available items.
#' @param item_information Vector with summarized information of each yet available item, with zeros for administered items 
#' (as returned by \code{\link{get_summarized_information}} with \code{pad = TRUE}).
#' @param estimate Vector containing the theta estimate, with its covariance matrix as an attribute.
#' @param number_answers The number of answers given thus far (length of \code{administered}).
#' @inheritParams shadowcat
#' @return Index of item with maximum information.
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


