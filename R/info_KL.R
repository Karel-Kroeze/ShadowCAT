#' Posterior expected Kullback-Leibler Information
#' 
#' Kullback-Leibler divergence based on the EAP and ML estimates of ability under the posterior distribution of theta.
#' 
#' @param test
#' @param person
#' @param simple
#' @return vector info
#' @export
PEKL <- function(test, person, theta_range = -3:3){
  # collect EAP estimate
  test$estimator <- "EAP"
  theta <- estimate(person, test)$estimate
  
  # we'll perform a very basic integration over the theta range
  # expand the grid for multidimensional models (number of calculations will be length(theta_range)**Q, which can still get quite high for high dimensionalities.)
  grid <- expand.grid(rep(list(theta_range), test$items$Q))

  # sum over values of the grid
  info <- rowSums(apply(grid, 1, KLB, theta0 = theta, test = test, person = person))
  
  return(info)
}

#' Kullback Leibler Divergence for given items and pairs of thetas x posterior density.
#' 
#' @param theta
#' @param theta0
#' @param items
#' @param test
#' @param person
#' @return vector info
KLB <- function(theta, theta0, test, person){
  available_items <- subset(test$items, person$available)
  administered_items <- subset(test$items, person$administered)
  P <- prob(test = test, theta = theta, items = available_items)$P
  P0 <- prob(test = test, theta = theta0, items = available_items)$P
  LL <- prob(test = test, theta = theta, person = person, items = administered_items, deriv = TRUE)$LL
  
  return(rowSums(P0 * (log(P0) - log(P)), na.rm = TRUE) * exp(LL))
}