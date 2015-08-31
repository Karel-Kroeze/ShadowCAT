#' Posterior expected Kullback-Leibler Information
#' 
#' Kullback-Leibler divergence based on the EAP and ML estimates of ability under the posterior distribution of theta.
#' Computes the numerical integral of the expectation of KL under the posterior distribution. The grid is simplified to a
#' range of -3, -2, ..., 3.
#'
#' Note that even with a simplified grid, the number of quadrature points which have to be calculated for each available item, 
#' at each step in the CAT is taken to the power Q. Use of KL information is likely to be slow in 3+ dimensional tests. The theta range to be evaluated can be
#' specified to partially reduce the amount of computations required.
#'   
#' TODO: Add references.
#' 
#' @param test Test object, see \code{\link{initTest}}.
#' @param person Person object, see \code{\link{initPerson}}.
#' @param theta_range Vector of theta values to be evaluated in the numerical integration. Using a sparser range may alleviate stress in higher dimensional tests.
#' @return Vector with PEKL information for each item.
#' @export
PEKL <- function(test, person, theta_range = -3:3){
  # collect EAP estimate
  if (test$estimator == "EAP") {
    theta <- person$estimate
  } else {
    test$estimator <- "EAP"
    theta <- estimate(person, test)$estimate
  }
  
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
  # TODO: wrap this into PEKL, do not recompute P0 for each theta (considering it is constant for the current posterior).
  available_items <- subset(test$items, person$available)
  administered_items <- subset(test$items, person$administered)
  P <- prob(test = test, theta = theta, items = available_items)$P
  P0 <- prob(test = test, theta = theta0, items = available_items)$P
  LL <- prob(test = test, theta = theta, person = person, items = administered_items, deriv = TRUE)$LL
  
  out <- rowSums(P0 * (log(P0) - log(P)), na.rm = TRUE) * exp(LL)
  
  return(out)
}