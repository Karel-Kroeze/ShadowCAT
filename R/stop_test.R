#' Done
#' 
#' Control function to check if the test is completed.
#' 
#' @param test test object
#' @param person person object
#' @return boolean completed
#' @export
stop_test <- function(person, test) {
  # current test duration
  l <- length(person$responses)
  
  # throw stopping criteria at the test until we can make something stick
  # max length rule, for use in combination with non-length based stopping rules (<= 0 never stops for max length!)
  if (test$max_n > 0 && l >= test$max_n) return(TRUE)
  
  # length stop rule ()
  if (test$stop$type == 'length') {
    if (l >= test$stop$n) return(TRUE)
  }
  
  # variance stop rule (target (vector of length Q with 'threshold' variances))
  if (test$stop$type == 'variance') {
    variance <- diag(attr(person$estimate, "variance"))
    if (all(variance < test$stop$target)) return(TRUE)
  }
  
  # if nothing stuck, don't stop!
  return(FALSE)
}

