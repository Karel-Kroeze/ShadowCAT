#' initiate person object for use in ShadowCAT
#' 
#' Sets bookkeeping variables and collects person information
#' 
#' 
#' @param items Itembank (required for bookkeeping reasons)
#' @param theta True theta
#' @param prior Prior covariance matrix
#' @return person ShadowCAT.person object
#' @export
initPerson <- function(items, theta = rep(0, items$Q), prior = diag(items$Q)) {
  available <- 1:items$K
  administered <- numeric(0)
  responses <- numeric(0)
  estimate <- rep(0, items$Q)
  variance <- matrix(NA, items$Q, items$Q)
  
  out <- list(theta = theta,
              prior = prior,
              estimate = estimate,
              variance = variance,
              available = available,
              administered = administered,
              responses = responses)
  
  attr(out, "class") <- "ShadowCAT.person"
  
  return(out)
}