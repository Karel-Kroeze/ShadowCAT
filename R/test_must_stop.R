#' Done
#' 
#' Control function to check if the test is completed.
#' 
#' @param number_answers number of answers given
#' @param estimate vector with current theta estimate, with covariance matrix as its attribute
#' @param min_n minimum test length; NULL means no minimum test length
#' @param max_n test length at which testing should stop (even if target has not been reached yet in case of variance stopping rule) 
#' @param stop_variance_target a vector with the variance target for each dimension; if null, variance is not taken into account
#' @param cutoffs matrix containing cut off values per dimension (columns) and test iteration (rows). First row contains cut off values for when no items have been
#' administered yet, second row for when one item has been administered, etc. If estimate + 3SE < cutoff for each dimension at certain iteration, test stops; NULL means no cut off values
#' @return TRUE if test should stop, FALSE otherwise
#' @export
test_must_stop <- function(number_answers, estimate, min_n, max_n, stop_variance_target = NULL, cutoffs = NULL) {
  min_number_items_reached <- is.null(min_n) || number_answers >= min_n
  max_number_items_reached <- number_answers >= max_n
  variance_target_reached <- !is.null(stop_variance_target) && all(diag(attr(estimate, "variance")) < stop_variance_target)
  estimate_far_enough_below_cutoff <- !is.null(cutoffs) && number_answers < max_n && all(estimate + 3 * sqrt(diag(attr(estimate, "variance"))) < cutoffs[number_answers + 1, ])
  if ((max_number_items_reached || variance_target_reached || estimate_far_enough_below_cutoff) && min_number_items_reached)
    TRUE
  else
    FALSE
}

