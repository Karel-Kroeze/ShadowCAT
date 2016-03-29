#' Latent trait estimation
#' 
#' maximum likelihood, maximum a posteriori and expected a posteriori estimates.
#' 
#' Obtains a latent trait estimate and variance of the estimate.
#' 
#' @section Maximum Likelihood and Maximum A Posteriori:
#' Maximum Likelihood and Maximum A-Posteriori estimates are based on a Newton-type non-linear minimization algorithm,
#' and handled with package \code{\link{nlm}}.
#'  
#' @section Expected A Posteriori:
#' Expected A-Posteriori estimates require the repeated evaluation of Q nested integrals, where Q is the dimensionality of the test.
#' This is performed with adaptive multidimensional Gauss-Hermite quadrature, and handled by package MultiGHQuad, see the documentation there for further details.
#' Note that the number of quadrature points used rises exponentially with the dimensionality of the test - use of EAP estimates with 
#' a 3+ dimensional test may not be a good idea.
#' 
#' @section Weighted Maximum Likelihood:
#' TODO: UPDATE WITH REFERENCES - MORE PRECISE DETAILS.
#' Note that WML estimation is not included. There is no satisfying solution to multidimensional Weighted Maximum Likelihood Estimation,
#' current WML estimators as used in other sources do not account for the covariance between dimensions. 
#' 
#' @section Variance:
#' Covariance matrix of the estimate is added to the estimate as an attribute.
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
#' prior_var_safe_nlm <- diag(number_dimensions)
#'
#' # obtain estimates
#' estimator <- "maximum_likelihood"
#' ML <- estimate_latent_trait(estimate, responses, prior, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, prior_var_safe_nlm)
#' estimator <- "maximum_aposteriori"
#' MAP <- estimate_latent_trait(estimate, responses, prior, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, prior_var_safe_nlm)
#' estimator <- "expected_aposteriori"
#' EAP <- estimate_latent_trait(estimate, responses, prior, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, prior_var_safe_nlm)
#' ML; MAP; EAP
#' 
#' # access variance
#' attr(ML, "variance")
#' 
#' # Note that expected_aposteriori takes considerably more time when dimensionality is higher...
#' number_dimensions <- 5
#' estimate <- rep(.3, number_dimensions)
#' alpha <- matrix(runif(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
#' lower_bound <- rep(-3, number_dimensions)
#' upper_bound <- rep(3, number_dimensions)
#' prior <- diag(number_dimensions) 
#' 
#' estimator <- "maximum_aposteriori"
#' system.time(estimate_latent_trait(estimate, responses, prior, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, prior_var_safe_nlm = diag(number_dimensions)))
#' estimator <- "expected_aposteriori"
#' system.time(estimate_latent_trait(estimate, responses, prior, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound))
#' 
#' @param estimate vector containing theta estimate, with covariance matrix as an attribute
#' @param responses vector with person responses
#' @param prior prior covariance matrix for theta
#' @param model string, one of '3PLM', 'GPCM', 'SM' or 'GRM', for the three-parameter logistic, generalized partial credit, sequential or graded response model respectively
#' @param administered vector with indeces of administered items
#' @param number_dimensions number of dimensions
#' @param estimator type of estimator to be used, one of "maximum_aposteriori", "maximum_likelihood", or "expected_aposteriori"
#' @param alpha matrix of alpha paramteres
#' @param beta matrix of beta paramteres
#' @param guessing matrix of guessing parameters
#' @param number_itemsteps_per_item vector containing the number of non missing cells per row of the beta matrix
#' @param lower_bound vector with lower bounds for theta per dimension; estimated theta values smaller than the lowerbound values are truncated to the lowerbound values
#' @param upper_bound vector with upper bounds for theta per dimension; estimated theta values larger than the upperbound values are truncated to the upperbound values
#' @param prior_var_safe_nlm if not NULL, expected a posteriori estimate with prior variance(s) equal to prior_var_safe_ml is computed instead of maximum_likelihood/maximum_aposteriori, if maximum_likelihood/maximum_aposteriori estimate fails. Can be a scalar 
#' (if variance for each dimension is equal) or vector
#' @return vector containing the updated estimate with the covariance matrix as attribute
#' @importFrom MultiGHQuad init.quad eval.quad
#' @export
estimate_latent_trait <- function(estimate, responses, prior, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, prior_var_safe_nlm = NULL) {
  result <- function() {
    updated_estimate <- get_updated_estimate_and_variance_attribute(estimator)
    trim_estimate(updated_estimate)
  }
  
  get_updated_estimate_and_variance_ml <- function() {
    # for now, simple nlm (TODO: look at optim, and possible reintroducing pure N-R).
    # We want a maximum, but nlm produces minima -> reverse function call.
    estimate <- tryCatch(nlm(f = likelihood_or_post_density, p = estimate, responses, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior, inverse_likelihood_or_post_density = TRUE)$estimate,
                         error = function(e) { switch_to_eap_if_requested() },
                         warning = function(w) { switch_to_eap_if_requested() })
    # TODO: We should really store info somewhere so we don't have to redo this (when using get_fisher_information based selection criteria).
    fisher_information_items <- get_fisher_information(estimate, model, number_dimensions, alpha, beta, guessing, number_itemsteps_per_item)
    fisher_information_test_so_far <- apply(fisher_information_items[,,administered, drop = FALSE], c(1, 2), sum)
    
    attr(estimate, "variance") <- solve(fisher_information_test_so_far)
    estimate
  }
  
  get_updated_estimate_and_variance_map <- function() {
    estimate <- tryCatch(nlm(f = likelihood_or_post_density, p = estimate, responses, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior, inverse_likelihood_or_post_density = TRUE)$estimate,
                         error = function(e) { switch_to_eap_if_requested() },
                         warning = function(w) { switch_to_eap_if_requested() })
    fisher_information_items <- get_fisher_information(estimate, model, number_dimensions, alpha, beta, guessing, number_itemsteps_per_item)
    fisher_information_test_so_far <- apply(fisher_information_items[,,administered, drop = FALSE], c(1, 2), sum) +
      solve(prior)
    attr(estimate, "variance") <- solve(fisher_information_test_so_far)
    estimate
  }
  
  get_updated_estimate_and_variance_eap <- function(prior) {
    # Multidimensional Gauss-Hermite Quadrature
    # TODO: prior mean is currently fixed at zero, update when/if possible.
    adapt <- if (length(responses) > 5 & !is.null(attr(estimate, 'variance'))) list(mu = estimate, Sigma = as.matrix(attr(estimate, "variance")))
    Q_dim_grid_quad_points <- init.quad(Q = number_dimensions,
                                        prior = list(mu = rep(0, number_dimensions), Sigma = prior),
                                        adapt = adapt,
                                        ip = number_gridpoints())
    eval.quad(FUN = likelihood_or_post_density, X = Q_dim_grid_quad_points, responses, model, administered, number_dimensions, estimator = "maximum_likelihood", alpha, beta, guessing)
  }
  
  get_updated_estimate_and_variance_attribute <- function(estimator) {
    switch(estimator,
           maximum_likelihood = get_updated_estimate_and_variance_ml(),
           maximum_aposteriori = get_updated_estimate_and_variance_map(),
           expected_aposteriori = get_updated_estimate_and_variance_eap(prior))
  }
  
  trim_estimate <- function(estimate) {
    estimate[which(estimate > upper_bound)] <- upper_bound[which(estimate > upper_bound)]
    estimate[which(estimate < lower_bound)] <- lower_bound[which(estimate < lower_bound)]
    estimate
  }
  
  switch_to_eap_if_requested <- function() {
    if (is.null(prior_var_safe_nlm))
      stop("something went wrong with nlm maximization and prior_var_safe_ml was set to NULL")
    get_updated_estimate_and_variance_eap(prior = diag(number_dimensions) * prior_var_safe_nlm)
  }
  
  number_gridpoints <- function() {
    if (number_dimensions < 5) 
      switch(number_dimensions, 50, 15, 6, 4)
    else
      3
  }
  
  result()
}
