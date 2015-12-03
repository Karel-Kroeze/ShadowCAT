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
PEKL <- function(test, person, theta_range = -3:3) {
  result <- function() {
    # we'll perform a very basic integration over the theta range
    # expand the grid for multidimensional models (number of calculations will be length(theta_range)**Q, which can still get quite high for high dimensionalities.)
    grid <- expand.grid(rep(list(theta_range), test$items$Q))
    
    # sum over values of the grid
    rowSums(apply(grid, 1, KLB, theta0 = get_theta_estimate(), test = test, person = person))
  }
  
  get_theta_estimate <- function() {
    # collect EAP estimate
    if (test$estimator == "EAP") {
      person$estimate
    } 
    else {
      test$estimator <- "EAP"
      estimate_latent_trait(person$estimate, person$responses, person$prior, test$items$model, person$administered,test$items$Q, test$estimator,test$items$pars$alpha, test$items$pars$beta, test$items$pars$guessing, test$items$pars$m, test$lowerBound, test$upperBound)
    }
  }
  
  #' Kullback Leibler Divergence for given items and pairs of thetas x posterior density.
  #' theta0 is theta estimated with EAP
  #' returns vector containing information for each yet available item
  KLB <- function(theta, theta0, test, person) {
    # TODO: wrap this into PEKL, do not recompute P0 for each theta (considering it is constant for the current posterior).
    P <- probabilities_and_likelihoods(theta, NULL, test$items$model, person$available, test$items$Q, test$estimator, test$items$pars$alpha, test$items$pars$beta, test$items$pars$guessing, output = "probs")
    P0 <- probabilities_and_likelihoods(theta0, NULL, test$items$model, person$available, test$items$Q, test$estimator, test$items$pars$alpha, test$items$pars$beta, test$items$pars$guessing, output = "probs")
    LL <- probabilities_and_likelihoods(theta, person$responses, test$items$model, person$administered, test$items$Q, test$estimator, test$items$pars$alpha, test$items$pars$beta, test$items$pars$guessing, prior = person$prior, output = "likelihoods")
    
    rowSums(P0 * (log(P0) - log(P)), na.rm = TRUE) * exp(LL)
  }
  
  validate <- function() {
    if (is.null(test))
      add_error("test", "is missing")
    if (is.null(person)) 
      add_error("person", "is missing")
  }
  
  invalid_result <- function() {
    list(errors = errors())
  }
  
  validate_and_run()
}

