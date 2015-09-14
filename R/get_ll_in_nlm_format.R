#' Log Likelihood
#' 
#' Internal - fetch appropriate elements from prob for nlm optimizer
#' @param theta
#' @param test
#' @param person
#' @param should values be reversed (useful for minimization, reverses LL as well as derivatives)
#' @param should log of output be returned
#' @return Log-Likelihood, as well as gradient and hessian attributes.
#' @importFrom stats nlm
LL <- function(theta, test, person, minimize = FALSE, log = TRUE) {
  # subset items that have a response
  test$items <- subset(test$items, person$administered)
  
  # get LL and derivatives.
  PROB <- prob(test, person, theta, deriv = TRUE)
  
  # prepare output
  out <- PROB$LL * (-1) ^ minimize
  
  # gradient and hessian for nlm optimizer.
  attr(out, "gradient") <- PROB$d1 * (-1) ^ minimize
  attr(out, "hessian") <- PROB$d2 * (-1) ^ minimize
  
  # return
  if ( ! log) return(invisible(exp(out)))  
  return(invisible(out))
}
