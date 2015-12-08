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
#' This is performed with adaptive multidimensional Gauss-Hermite quadrature, and handled by package MultiGHQuad, see the documentation there for further details.
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
#' number_dimensions <- 1
#' estimate <- rep(.3, number_dimensions)
#' model <- "3PLM"
#' number_items <- 50
#' responses <- rep(c(1, 0), 17)
#' administered <- c(6:20, 31:49)
#' alpha <- matrix(runif(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
#' beta <- matrix(rnorm(number_items), nrow = number_items, ncol = 1)
#' guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
#' number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
#' lower_bound <- rep(-3, number_dimensions)
#' upper_bound <- rep(3, number_dimensions)
#' prior <- diag(1)
#' prior_var_safe_ml <- NULL
#'
#' # obtain estimates
#' estimator <- "ML"
#' ML <- estimate_latent_trait(estimate, responses, prior, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, prior_var_safe_ml = NULL)
#' estimator <- "MAP"
#' MAP <- estimate_latent_trait(estimate, responses, prior, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, prior_var_safe_ml = NULL)
#' estimator <- "EAP"
#' # EAP <- estimate_latent_trait(estimate, responses, prior, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, prior_var_safe_ml = NULL)
#' ML; MAP; #EAP
#' 
#' # access variance
#' attr(ML, "variance")
#' 
#' # Note that EAP takes considerably more time when dimensionality is higher...
#' number_dimensions <- 5
#' estimate <- rep(.3, number_dimensions)
#' alpha <- matrix(runif(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
#' lower_bound <- rep(-3, number_dimensions)
#' upper_bound <- rep(3, number_dimensions)
#' prior <- diag(number_dimensions) 
#' 
#' estimator <- "MAP"
#' system.time(estimate_latent_trait(estimate, responses, prior, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, prior_var_safe_ml = NULL))
#' estimator <- "EAP"
#' # system.time(estimate_latent_trait(estimate, responses, prior, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, prior_var_safe_ml = NULL))
#' 
#' @param estimate theta estimate
#' @param responses person responses
#' @param prior prior covariance matrix for theta
#' @param model string, one of '3PLM', 'GPCM', 'SM' or 'GRM', for the three-parameter logistic, generalized partial credit, sequential or graded response model respectively
#' @param administered indeces of administered items
#' @param number_dimensions number of dimensions
#' @param estimator type of estimator to be used, one of "MAP" (Maximum a posteriori estimation), "ML" (maximum likelihood), or "EAP" (Expected A Posteriori Estimation)
#' @param alpha matrix of alpha paramteres
#' @param beta matrix of beta paramteres
#' @param guessing matrix of guessing parameters
#' @param number_itemsteps_per_item vector containing the number of non missing cells per row of the beta matrix
#' @param lower_bound vector with lower bounds for theta per dimension; estimated theta values smaller than the lowerbound values are truncated to the lowerbound values
#' @param upper_bound vector with upper bounds for theta per dimension; estimated theta values larger than the upperbound values are truncated to the upperbound values
#' @param prior_var_safe_ml if not NULL, MAP estimate with prior variance equal to prior_var_safe_ml is computed instead of ML, if ML estimate fails
#' @return person object, amended with the new estimate.
#' @importFrom MultiGHQuad init.quad eval.quad
#' @export
estimate_latent_trait <- function(estimate, responses, prior, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, prior_var_safe_ml = NULL) {
  result <- function() {
    updated_estimate <- get_updated_estimate_and_variance_attribute(estimator)
    trim_estimate(updated_estimate)
  }
  
  get_updated_estimate_and_variance_ml <- function() {
    if (is.null(prior_var_safe_ml))
      get_updated_estimate_and_variance_ml_unsafe()
    else
      get_updated_estimate_and_variance_ml_safe()      
  }
  
  get_updated_estimate_and_variance_ml_unsafe <- function() {
    # for now, simple nlm (TODO: look at optim, and possible reintroducing pure N-R).
    # We want a maximum, but nlm produces minima -> reverse function call.
    estimate <- nlm(f = probabilities_and_likelihoods, p = estimate, responses, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior, inverse_likelihoods = TRUE, output = "likelihoods")$estimate
    
    # TODO: We should really store info somewhere so we don't have to redo this (when using get_fisher_information based selection criteria).
    fisher_information_items <- get_fisher_information(estimate, model, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item)
    fisher_information_test_so_far <- apply(fisher_information_items[,,administered, drop = FALSE], c(1, 2), sum)
    
    attr(estimate, "variance") <- solve(fisher_information_test_so_far)
    estimate
  }
  
  get_updated_estimate_and_variance_ml_safe <- function() { 
    # suppress warnings and errors and do MAP with flat prior instead
    estimate <- tryCatch(nlm(f = probabilities_and_likelihoods, p = estimate, responses, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior, inverse_likelihoods = TRUE, output = "likelihoods")$estimate,
                         error = function(e) { safe_ml() },
                         warning = function(w) { safe_ml() })
    
    fisher_information_items <- get_fisher_information(estimate, model, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item)
    fisher_information_test_so_far <- apply(fisher_information_items[,,administered, drop = FALSE], c(1, 2), sum)
    
    attr(estimate, "variance") <- tryCatch(solve(fisher_information_test_so_far),
                                           error = function(e) { solve(fisher_information_test_so_far + diag(number_dimensions) * prior_var_safe_ml) },
                                           warning = function(w) { solve(fisher_information_test_so_far + diag(number_dimensions) * prior_var_safe_ml) })
    
    estimate
  }
  
  get_updated_estimate_and_variance_map <- function() {
    # note that prior is applied in probabilities_and_likelihoods (incorrectly it seems, but still).
    estimate <- nlm(f = probabilities_and_likelihoods, p = estimate, responses, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior, inverse_likelihoods = TRUE, output = "likelihoods")$estimate
    
    fisher_information_items <- get_fisher_information(estimate, model, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item)
    fisher_information_test_so_far <- apply(fisher_information_items[,,administered, drop = FALSE], c(1, 2), sum) +
      solve(prior)
    # inverse
    attr(estimate, "variance") <- solve(fisher_information_test_so_far)
    estimate
  }
  
  get_updated_estimate_and_variance_eap <- function() {
    # Multidimensional Gauss-Hermite Quadrature
    # TODO: prior mean is currently fixed at zero, update when/if possible.
    # TODO: allow setting ip through internals argument(s)
    adapt <- if (length(responses) > 5 & !is.null(attr(estimate, 'variance'))) list(mu = estimate, Sigma = as.matrix(attr(estimate, "variance")))
    Q_dim_grid_quad_points <- init.quad(Q = number_dimensions,
                                        prior = list(mu = rep(0, number_dimensions), Sigma = prior),
                                        adapt = adapt,
                                        ip = switch(number_dimensions, 50, 15, 6, 4, 3))
    eval.quad(FUN = probabilities_and_likelihoods, X = Q_dim_grid_quad_points, responses, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior, output = "likelihoods")
  }
  
  get_updated_estimate_and_variance_attribute <- function(estimator) {
    switch(estimator,
           ML = get_updated_estimate_and_variance_ml(),
           MAP = get_updated_estimate_and_variance_map(),
           EAP = get_updated_estimate_and_variance_eap())
  }
  
  trim_estimate <- function(estimate) {
    # TODO: make debug output toggleable
    # if (any(person$estimate > test$upperBound | person$estimate < test$lowerBound)) cat("Estimate outside boundaries (k =", length(person$responses), "estimate =", paste0(round(person$estimate, 2), collapse = ", "), ").\n")
    estimate[which(estimate > upper_bound)] <- upper_bound[which(estimate > upper_bound)]
    estimate[which(estimate < lower_bound)] <- lower_bound[which(estimate < lower_bound)]
    estimate
  }
  
  safe_ml <- function() {
    estimator <- "MAP"
    prior <- diag(number_dimensions) * prior_var_safe_ml
    nlm(f = probabilities_and_likelihoods, p = estimate, responses, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior, inverse_likelihoods = TRUE, output = "likelihoods")$estimate
  }
  
  result()
}