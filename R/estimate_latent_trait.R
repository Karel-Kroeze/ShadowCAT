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
#' This is performed with an adaptive Riemannsum or multidimensional Gauss-Hermite quadrature, the latter handled by package MultiGHQuad, see the documentation there for further details.
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
#' answers <- rep(c(1, 0), 17)
#' administered <- c(6:20, 31:49)
#' alpha <- matrix(runif(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
#' beta <- matrix(rnorm(number_items), nrow = number_items, ncol = 1)
#' guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
#' number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
#' prior_form <- "normal"
#' prior_parameters <- list(mu = 0, Sigma = diag(1))
#'
#' # obtain estimates
#' estimator <- "maximum_aposteriori" # maximum a posteriori combined with uniform prior is equivalent to maximum likelihood with bounds
#' ML <- estimate_latent_trait(estimate, answers, prior_form = "uniform", prior_parameters = list(lower_bound = -4, upper_bound = 4), model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item)
#' estimator <- "maximum_aposteriori"
#' MAP <- estimate_latent_trait(estimate, answers, prior_form, prior_parameters, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item)
#' estimator <- "expected_aposteriori"
#' EAP <- estimate_latent_trait(estimate, answers, prior_form, prior_parameters, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item)
#' ML; MAP; EAP
#' 
#' # access variance
#' attr(ML, "variance")
#' 
#' # Note that expected_aposteriori takes considerably more time when dimensionality is higher...
#' number_dimensions <- 5
#' estimate <- rep(.3, number_dimensions)
#' alpha <- matrix(runif(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
#' prior_form <- "normal"
#' prior_parameters <- list(mu = rep(0, number_dimensions), Sigma = diag(number_dimensions))
#' 
#' estimator <- "maximum_aposteriori"
#' system.time(estimate_latent_trait(estimate, answers, prior_form, prior_parameters, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item))
#' estimator <- "expected_aposteriori"
#' system.time(estimate_latent_trait(estimate, answers, prior_form, prior_parameters, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item))
#' 
#' @param estimate Vector containing theta estimate, with covariance matrix as an attribute
#' @param answers Vector with person answers
#' @param prior_form String indicating the form of the prior; one of "normal" or "uniform"
#' @param prior_parameters List containing mu and Sigma of the normal prior: list(mu = ..., Sigma = ...), or 
#' the upper and lower bound of the uniform prior: list(lower_bound = ..., upper_bound = ...). Sigma should always
#' be in matrix form.
#' @param model String, one of '3PLM', 'GPCM', 'SM' or 'GRM', for the three-parameter logistic, generalized partial credit, sequential or graded response model respectively
#' @param administered Vector with indices of administered items
#' @param number_dimensions Number of dimensions
#' @param estimator Type of estimator to be used, one of "maximum_aposteriori", "maximum_likelihood", or "expected_aposteriori"
#' @param alpha Matrix of alpha paramteres
#' @param beta Matrix of beta paramteres
#' @param guessing Matrix of guessing parameters
#' @param number_itemsteps_per_item Vector containing the number of non missing cells per row of the beta matrix
#' @param safe_eap Only relevant if estimator is espected_aposteriori. 
#' TRUE if estimator should switch to maximum aposteriori if the integration algorithm results in an error.
#' @param eap_estimation_procedure String indicating the estimation procedure if estimator is expected aposteriori. One of "riemannsum" for integration via Riemannsum or
#' "gauss_hermite_quad" for integration via Gaussian Hermite Quadrature. 
#' @return vector containing the updated estimate with the covariance matrix as attribute
#' @importFrom MultiGHQuad init.quad eval.quad
#' @export
estimate_latent_trait <- function(estimate, answers, prior_form, prior_parameters, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, safe_eap = FALSE, eap_estimation_procedure = "riemannsum") {
  result <- function() {
    get_updated_estimate_and_variance_attribute()
  }
  
  get_updated_estimate_and_variance_attribute <- function() {
    switch(estimator,
           maximum_likelihood = get_updated_estimate_and_variance_ml(),
           maximum_aposteriori = get_updated_estimate_and_variance_map(),
           expected_aposteriori = tryCatch(get_updated_estimate_and_variance_eap(),
                                           error = function(e) { if (safe_eap)
                                                                   get_updated_estimate_and_variance_map()
                                                                 else
                                                                   stop(e) }))
  }
  
  get_updated_estimate_and_variance_ml <- function() {
    # for now, simple nlm (TODO: look at optim, and possible reintroducing pure N-R).
    # We want a maximum, but nlm produces minima -> reverse function call.
    estimate <- nlm(f = likelihood_or_post_density, p = estimate, answers = answers, model = model, items_to_include = administered, number_dimensions = number_dimensions, estimator = estimator, alpha = alpha, beta = beta, guessing = guessing, prior_parameters = prior_parameters, inverse_likelihood_or_post_density = TRUE)$estimate
    fisher_information_items <- get_fisher_information(estimate, model, number_dimensions, alpha, beta, guessing, number_itemsteps_per_item)
    fisher_information_test_so_far <- apply(fisher_information_items[,,administered, drop = FALSE], c(1, 2), sum)
    attr(estimate, "variance") <- solve(fisher_information_test_so_far)
    estimate
  }
  
  get_updated_estimate_and_variance_map <- function() {
    switch(prior_form,
           normal = get_updated_estimate_and_variance_map_normal_prior(),
           uniform = get_updated_estimate_and_variance_map_uniform_prior())
  }
  
  get_updated_estimate_and_variance_eap <- function() {
    switch(eap_estimation_procedure,
           gauss_hermite_quad = get_updated_estimate_and_variance_eap_gauss_hermite_quad(),
           riemannsum = get_updated_estimate_and_variance_eap_riemannsum())
  }
  
  get_updated_estimate_and_variance_map_normal_prior <- function() {
    estimate <- nlm(f = likelihood_or_post_density, p = estimate, 
                    answers = answers, model = model, items_to_include = administered, number_dimensions = number_dimensions, estimator = estimator, alpha = alpha, beta = beta, guessing = guessing, prior_parameters = prior_parameters, inverse_likelihood_or_post_density = TRUE)$estimate
    fisher_information_items <- get_fisher_information(estimate, model, number_dimensions, alpha, beta, guessing, number_itemsteps_per_item)
    fisher_information_test_so_far <- apply(fisher_information_items[,,administered, drop = FALSE], c(1, 2), sum) + solve(prior_parameters$Sigma)
    attr(estimate, "variance") <- solve(fisher_information_test_so_far)
    estimate
  }
  
  get_updated_estimate_and_variance_map_uniform_prior <- function() {
    estimate <- constrOptim(theta = move_values_to_means(values = estimate, means = rowMeans(matrix(c(prior_parameters$lower_bound, prior_parameters$upper_bound), ncol = 2)), amount_change = rep(.001, number_dimensions)), 
                            f = likelihood_or_post_density,
                            grad = function(theta, answers, model, items_to_include, number_dimensions, estimator, alpha, beta, guessing, prior_parameters, return_log_likelihood_or_post_density = TRUE, inverse_likelihood_or_post_density) {
                              likelihood <- likelihood_or_post_density(theta, answers, model, items_to_include, number_dimensions, estimator, alpha, beta, guessing, prior_parameters, return_log_likelihood_or_post_density, inverse_likelihood_or_post_density)
                              attr(likelihood, "gradient")
                            },
                            ui = rbind(diag(number_dimensions), -diag(number_dimensions)), 
                            ci = c(prior_parameters$lower_bound, -prior_parameters$upper_bound),
                            answers = answers, model = model, items_to_include = administered, number_dimensions = number_dimensions, estimator = "maximum_likelihood", alpha = alpha, beta = beta, guessing = guessing, prior_parameters = NULL, inverse_likelihood_or_post_density = TRUE)$par
    fisher_information_items <- get_fisher_information(estimate, model, number_dimensions, alpha, beta, guessing, number_itemsteps_per_item)
    fisher_information_test_so_far <- apply(fisher_information_items[,,administered, drop = FALSE], c(1, 2), sum)
    attr(estimate, "variance") <- solve(fisher_information_test_so_far)
    estimate
  }
  
  get_updated_estimate_and_variance_eap_gauss_hermite_quad <- function() {
    adapt <- if (length(answers) > 5 & !is.null(attr(estimate, 'variance'))) list(mu = estimate, Sigma = as.matrix(attr(estimate, "variance")))
    Q_dim_grid_quad_points <- init.quad(Q = number_dimensions,
                                        prior = prior_parameters,
                                        adapt = adapt,
                                        ip = number_gridpoints(),
                                        forcePD = TRUE)
    eval.quad(FUN = likelihood_or_post_density, X = Q_dim_grid_quad_points, 
              answers = answers, model = model, items_to_include = administered, number_dimensions = number_dimensions, estimator = "maximum_likelihood", alpha = alpha, beta = beta, guessing = guessing)
  }
  
  get_updated_estimate_and_variance_eap_riemannsum <- function() {
    adapt <- if (length(answers) > 5 & !is.null(attr(estimate, 'variance'))) list(mu = estimate, Sigma = as.matrix(attr(estimate, "variance")))
    get_eap_estimate_riemannsum(dimension = number_dimensions, 
               likelihood = likelihood_or_post_density, 
               prior_form = prior_form,
               prior_parameters = prior_parameters,
               adapt = adapt,
               number_gridpoints = number_gridpoints(),
               answers = answers, model = model, items_to_include = administered, number_dimensions = number_dimensions, estimator = "maximum_likelihood", alpha = alpha, beta = beta, guessing = guessing, return_log_likelihood_or_post_density = FALSE)
  }
  
  number_gridpoints <- function() {
    if (number_dimensions < 5) 
      switch(number_dimensions, 50, 15, 6, 4)
    else
      3
  }
  
  result()
}
