#' Done
#' 
#' Control function to check if the test is completed.
#' 
#' @param number_responses number of responses given
#' @param variance_current_estimate variance(s) of the current estimate(s)
#' @param min_n minimum test length; NULL means no minimum test length
#' @param max_n test length at which testing should stop (even if target has not been reached yet in case of variance stopping rule; NULL never stops for max length) 
#' @param stop_variance_target a vector with the variance target for each dimension; if null, variance is not taken into account
#' @return boolean completed
#' @export
test_must_stop <- function(number_responses, variance_current_estimate, min_n, max_n, stop_variance_target = NULL) {
  min_number_items_reached <- is.null(min_n) || number_responses >= min_n
  max_number_items_reached <- !is.null(max_n) && number_responses >= max_n
  variance_target_reached <- !is.null(stop_variance_target) && all(diag(variance_current_estimate) < stop_variance_target)
  if ((max_number_items_reached || variance_target_reached) && min_number_items_reached)
    TRUE
  else
    FALSE
}

