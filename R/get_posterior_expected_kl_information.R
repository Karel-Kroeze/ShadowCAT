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
#' @param estimate vector with current theta estimate
#' @param model string, one of '3PLM', 'GPCM', 'SM' or 'GRM', for the three-parameter logistic, generalized partial credit, sequential or graded response model respectively.
#' @param answers vector with person answers
#' @param administered vector with indices of administered items
#' @param available vector with indices of available items
#' @param number_dimensions number of dimensions
#' @param estimator type of estimator to be used, one of "maximum_aposteriori", "maximum_likelihood", or "expected_aposteriori"
#' @param alpha matrix containing the alpha parameters
#' @param beta matrix containing the beta parameters
#' @param guessing matrix containing the quessing parameters
#' @param prior_form String indicating the form of the prior; one of "normal" or "uniform"
#' @param prior_parameters List containing mu and Sigma of the normal prior: list(mu = ..., Sigma = ...), or 
#' the upper and lower bound of the uniform prior: list(lower_bound = ..., upper_bound = ...). Sigma should always
#' be in matrix form.
#' @param number_itemsteps_per_item vector containing the number of non missing cells per row of the beta matrix
#' @param eap_estimation_procedure String indicating the estimation procedure for the expected aposteriori estimate, which is computed
#' here if it is not the requested estimator in shadowcat(). One of "riemannsum" for integration via Riemannsum or
#' "gauss_hermite_quad" for integration via Gaussian Hermite Quadrature. 
#' @return Vector with PEKL information for each yet available item.
#' @export
get_posterior_expected_kl_information <- function(estimate, model, answers, administered, available, number_dimensions, estimator, alpha, beta, guessing, prior_form, prior_parameters, number_itemsteps_per_item, eap_estimation_procedure = "riemannsum") {
  result <- function() {
    # we'll perform a very basic integration over the theta range
    # expand the grid for multidimensional models (number of calculations will be length(theta_values)**Q, which can still get quite high for high dimensionalities.)
    probabilities_given_eap_estimate <- get_probs_and_likelihoods_per_item(get_theta_estimate(), model, get_subset(alpha, available), get_subset(beta, available), get_subset(guessing, available), with_likelihoods = FALSE)$P
    log_probabilities_given_eap_estimate <- log(probabilities_given_eap_estimate)
    theta_grid <- get_theta_grid()
    row_or_vector_sums(apply(theta_grid, 1, kullback_leibler_divergence, probabilities_given_eap_estimate = probabilities_given_eap_estimate, log_probabilities_given_eap_estimate = log_probabilities_given_eap_estimate))
  }
  
  get_theta_estimate <- function() {
    if (estimator == "expected_aposteriori")
      estimate
    else
      estimate_latent_trait(estimate, answers, prior_form, prior_parameters, model, administered, number_dimensions, estimator = "expected_aposteriori", alpha, beta, guessing, number_itemsteps_per_item, eap_estimation_procedure = eap_estimation_procedure)
  }
  
  #' Kullback Leibler Divergence for given items and pairs of thetas x posterior density.
  #' returns vector containing information for each yet available item
  kullback_leibler_divergence <- function(theta, probabilities_given_eap_estimate, log_probabilities_given_eap_estimate) {
    probabilities_given_theta <- get_probs_and_likelihoods_per_item(theta, model, get_subset(alpha, available), get_subset(beta, available), get_subset(guessing, available), with_likelihoods = FALSE)$P
    likelihood_or_post_density_theta <- likelihood_or_post_density(theta, answers, model, administered, number_dimensions, estimator = estimator_likelihood_or_post_density(), alpha, beta, guessing, prior_parameters = prior_parameters, return_log_likelihood_or_post_density = FALSE)
    rowSums(probabilities_given_eap_estimate * (log_probabilities_given_eap_estimate - log(probabilities_given_theta)), na.rm = TRUE) * likelihood_or_post_density_theta
  }
  
  estimator_likelihood_or_post_density <- function() {
    if (prior_form == "uniform")
      "maximum_likelihood"
    else
      "expected_aposteriori"
  }
  
  get_theta_grid <- function() {
    grid_list <- get_grid_list()
    expand.grid(grid_list)
  }
  
  get_grid_list <- function() {
    if (number_dimensions == 1 && prior_form == "uniform")
      list(seq(prior_parameters$lower_bound, prior_parameters$upper_bound, length.out = 21))
    else if (number_dimensions == 1 && prior_form == "normal")
      list(seq(-4, 4, length.out = 21))
    else if (prior_form == "uniform")
      lapply(1:number_dimensions, function(dim) { seq(prior_parameters$lower_bound[dim], prior_parameters$upper_bound[dim], length.out = 9) })
    else
      rep(list(-4:4), number_dimensions)
  }  
  
  result()
}

