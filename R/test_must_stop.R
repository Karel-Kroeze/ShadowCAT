#' Done
#' 
#' Control function to check if the test is completed.
#' 
#' @param number_responses number of responses given
#' @param variance_current_estimate variance(s) of the current estimate(s)
#' @param max_n test length at which testing should stop (even if target has not been reached yet in case of variance stopping rule; n <= 0 never stops for max length) 
#' @param stop_type one of 'length' or 'variance'
#' @param stop_variance_target if stop_type is 'variance', a vector with the variance target for each dimension
#' @return boolean completed
#' @export
test_must_stop <- function(number_responses, variance_current_estimate, max_n, stop_type, stop_variance_target = NULL) {
  max_number_items_reached <- max_n > 0 && number_responses >= max_n
  variance_target_reached <- stop_type == 'variance' && all(diag(variance_current_estimate) < stop_variance_target)
  if (max_number_items_reached || variance_target_reached)
    TRUE
  else
    FALSE
}

