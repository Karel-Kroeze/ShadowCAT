#' Latent trait estimation
#' 
#' ML, MAP and EAP estimates.
#' 
#' Obtains a latent trait estimate and variance of the estimate.
#' 
#' @section ML and MAP:
#' Maximum Likelihood and Maximum A-Posteriori estimates are based on a Newton-type non-linear minimization algorithm,
#' and handled with package \code{\link{nlm}}.
#'  
#' @section EAP:
#' Expected A-Posteriori estimates require the repeated evaluation of Q nested integrals, where Q is the dimensionality of the test.
#' This is performed with adaptive multidimensional Gauss-Hermite quadrature, and handled by package \code{\link{MultiGHQuad}}, see the documentation there for further details.
#' Note that the number of quadrature points used rises exponentially with the dimensionality of the test - use of EAP estimates with 
#' a 3+ dimensional test may not be a good idea.
#' 
#' @section WML:
#' TODO: UPDATE WITH REFERENCES - MORE PRECISE DETAILS.
#' Note that WML estimation is not included. There is no satisfying solution to multidimensional Weighted Maximum Likelihood Estimation,
#' current WML estimators as used in other sources do not account for the covariance between dimensions. 
#' 
#' @section Variance:
#' Variance of the estimate is added to the estimate as an attribute.
#' 
#' @examples 
#' # create a basic test + person
#' items <- createTestBank("GPCM")
#' test <- initTest(items, estimator = "ML")
#' person <- initPerson(items)
#' 
#' # answer a few items
#' person <- answer(person, test, sample(items$K, 10))
#' 
#' # obtain estimates
#' ML <- estimate(person, test)$estimate
#' test$estimator <- "MAP"
#' MAP <- estimate(person, test)$estimate
#' test$estimator <- "EAP"
#' EAP <- estimate(person, test)$estimate
#' ML; MAP; EAP
#' 
#' # access variance
#' attr(ML, "variance")
#' 
#' # Note that EAP takes considerably more time when dimensionality is higher...
#' items5d <- createTestBank("GPCM", Q=5)
#' test5dEAP <- initTest(items, estimator = "EAP")
#' test5dMAP <- initTest(items, estimator = "MAP")
#' person5d <- answer(initPerson(items), test, sample(items5d$K, 10))
#' 
#' system.time(estimate(person5d, test5dEAP))
#' system.time(estimate(person5d, test5dMAP))
#' 
#' @param person Person object, see \code{\link{initPerson}}.
#' @param test Test object, see \code{\link{initTest}}.
#' @return person object, amended with the new estimate.
#' @importFrom MultiGHQuad init.quad eval.quad
#' @export
estimate <- function(person, test) {
  result <- function() {
    updated_estimate <- get_updated_estimate_and_variance_attribute(test$estimator)     
    person$estimate <- trim_estimate(updated_estimate)    
    invisible(person)
  }
  
  get_updated_estimate_and_variance_ml <- function() {
    # for now, simple nlm (TODO: look at optim, and possible reintroducing pure N-R).
    # We want a maximum, but nlm produces minima -> reverse function call. 
    # LL is the target function, test, person and minimize need to be passed on. We also want the value of the hessian at the final estimate.
    
    # note that prior is applied in LL (incorrectly it seems, but still).
    # suppress warnings and errors and do EAP instead. RM I have removed this option, I don't want users to get something they think
    # is something else. Also, if estimator was ML, the default prior is used which may not make sense.
    person$estimate <- nlm(f = LL, p = person$estimate, test = test, person = person, minimize = TRUE)$estimate # passed on to LL, reverses polarity.
    
    # TODO: We should really store info somewhere so we don't have to redo this (when using FI based selection criteria).
    fisher_information_items <- FI(test, person)
    fisher_information_test_so_far <- apply(fisher_information_items[,,person$administered, drop = FALSE], c(1, 2), sum)

    # inverse
    attr(person$estimate, "variance") <- solve(fisher_information_test_so_far)
    person$estimate
  }
  
  get_updated_estimate_and_variance_map <- function() {
    # for now, simple nlm (TODO: look at optim, and possible reintroducing pure N-R).
    # We want a maximum, but nlm produces minima -> reverse function call. 
    # LL is the target function, test, person and minimize need to be passed on. We also want the value of the hessian at the final estimate.
    
    # note that prior is applied in LL (incorrectly it seems, but still).
    # suppress warnings and errors and do EAP instead. RM I have removed this option, I don't want users to get something they think
    # is something else. Also, if estimator was ML, the default prior is used which may not make sense.
    person$estimate <- nlm(f = LL, p = person$estimate, test = test, person = person, minimize = TRUE)$estimate # passed on to LL, reverses polarity.
    
    # TODO: We should really store info somewhere so we don't have to redo this (when using FI based selection criteria).
    fisher_information_items <- FI(test, person)
    fisher_information_test_so_far <- apply(fisher_information_items[,,person$administered, drop = FALSE], c(1, 2), sum) +
                                      solve(person$prior)
    # inverse
    attr(person$estimate, "variance") <- solve(fisher_information_test_so_far)
    person$estimate
  }
  
  get_updated_estimate_and_variance_eap <- function() {
    # Multidimensional Gauss-Hermite Quadrature
    # TODO: prior mean is currently fixed at zero, update when/if possible.
    # TODO: allow setting ip through internals argument(s)
    adapt <- if (length(person$responses) > 5 & !is.null(attr(person$estimate, 'variance')))
               list(mu = person$estimate, Sigma = as.matrix(attr(person$estimate, "variance")))           
    Q_dim_grid_quad_points <- init.quad(Q = test$items$Q,
                              prior = list(mu = rep(0, test$items$Q), Sigma = person$prior),
                              adapt = adapt,
                              ip = switch(test$items$Q, 50, 15, 6, 4, 3))
    eval.quad(FUN = LL, X = Q_dim_grid_quad_points, test = test, person = person)
  }
  
  get_updated_estimate_and_variance_attribute <- function(estimator) {
    switch(estimator,
           ML = get_updated_estimate_and_variance_ml(),
           MAP = get_updated_estimate_and_variance_map(),
           EAP = get_updated_estimate_and_variance_eap())
  }
  
  trim_estimate <- function(estimate) {
    # enforce boundaries.
    # TODO: make debug output toggleable
    # if (any(person$estimate > test$upperBound | person$estimate < test$lowerBound)) cat("Estimate outside boundaries (k =", length(person$responses), "estimate =", paste0(round(person$estimate, 2), collapse = ", "), ").\n")
    estimate[which(estimate > test$upperBound)] <- test$upperBound[which(estimate > test$upperBound)]
    estimate[which(estimate < test$lowerBound)] <- test$lowerBound[which(estimate < test$lowerBound)]
    estimate
  }
  
  validate <- function() {
    if (is.null(person))
      add_error("person", "is missing")
    if (is.null(test))
      add_error("test", "is missing")
  }
  
  invalid_result <- function() {
    list(errors = errors())
  }
  
  validate_and_run()
}

