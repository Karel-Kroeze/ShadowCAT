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
  # throw stopping criteria at the test until we can make something stick
  # max length rule, for use in combination with non-length based stopping rules (<= 0 never stops for max length!)
  if (max_n > 0 && number_responses >= max_n) 
    TRUE
  
  # variance stop rule (target (vector of length Q with 'threshold' variances))
  else if (stop_type == 'variance' && all(diag(variance_current_estimate) < stop_variance_target))
    TRUE

  # if nothing stuck, don't stop!
  else
    FALSE
}

