#' initiate person object for use in ShadowCAT
#' 
#' Sets bookkeeping variables and collects person information
#' 
#' 
#' @param items Itembank (required for bookkeeping reasons)
#' @param theta True theta; only useful for simulations?
#' @param prior Prior covariance matrix
#' @return person ShadowCAT.person object
#' @export
initPerson <- function(items, theta = rep(0, items$Q), prior = diag(items$Q)) {
  estimate <- rep(0, items$Q)
  attr(estimate, 'variance') <- prior
  
  out <- list(theta = theta,
              prior = prior,
              estimate = estimate,
              variance = prior, # this double variance object (attribute and separate) is confusing
              available = 1:items$K,
              administered = numeric(0),
              responses = numeric(0))
  
  attr(out, "class") <- "ShadowCAT.person"
  
  out
}