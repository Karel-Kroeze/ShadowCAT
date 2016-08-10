#' Terminate test
#' 
#' Control function to check if the test is completed.
#' 
#' @param number_answers Number of answers given.
#' @param estimate Vector with current theta estimate, with covariance matrix as its attribute.
#' @param max_n Test length at which testing should stop (even if target has not been reached yet in case of variance stopping rule).
#' @param min_n Minimum test length; \code{NULL} means no minimum test length.
#' @param variance_target Vector with the variance target for each dimension; \code{NULL} means no variance target.
#' @param cutoffs Matrix containing cut off values per dimension (columns) and test iteration (rows). First row contains cut off values for when no items have been
#' administered yet, second row for when one item has been administered, etc. If estimate + 3SE < cutoff for each dimension at a certain iteration, test stops; \code{NULL} means no cut off values.
#' @return \code{TRUE} if test should stop, \code{FALSE} otherwise.
terminate_test <- function(number_answers, estimate, max_n, min_n = NULL, variance_target = NULL, cutoffs = NULL) {
  min_number_items_reached <- is.null(min_n) || number_answers >= min_n
  max_number_items_reached <- number_answers >= max_n
  variance_target_reached <- !is.null(variance_target) && all(diag(attr(estimate, "variance")) < variance_target)
  estimate_far_enough_below_cutoff <- !is.null(cutoffs) && number_answers < max_n && all(estimate + 3 * sqrt(diag(attr(estimate, "variance"))) < cutoffs[number_answers + 1, ])
  if ((max_number_items_reached || variance_target_reached || estimate_far_enough_below_cutoff) && min_number_items_reached)
    TRUE
  else
    FALSE
}

