#' Posterior expected Kullback-Leibler Information
#' 
#' Kullback-Leibler divergence based on the EAP and ML estimates of ability under the posterior distribution of theta.
#' Computes the numerical integral of the expectation of KL under the posterior distribution.
#' 
#' @details
#' Note that even with a simplified grid, the number of quadrature points which have to be calculated for each available item, 
#' at each step in the CAT is taken to the power Q. Use of KL information is likely to be slow in 3+ dimensional tests.
#'
#' @references
#'  \itemize{
#'  \item Chang, H.-H, & Ying, Z. (1996). A Global Information Approach to Computerized Adaptive 
#'  Testing. Applied Psychological Measurement, 20(3), 213 - 229. doi: 10.1177/014662169602000303.
#'  \item Mulder, J., & van der Linden, W. J. (2010). Multidimentional Adaptive testing with 
#'  Kullback-Leibler information item selection. In W. J. van der Linden & C. A. W. Glas (Eds.), 
#'  Elements of adaptive testing (pp. 77 - 101). New York: Springer.
#'  \item Wang, C., Chang, H.-H., & Boughton, K. A. (2010). Kullback-Leibler Information and Its 
#'  Applications in Multi-Dimensional Adaptive Testing. Psychometrika, 76(1), 13 - 39. doi:10.1007/s11336-010-9186-0.
#' }
#' 
#' @param estimate Vector containing current theta estimate, with covariance matrix as an attribute.
#' @param answers Vector with answers to administered items.
#' @param administered vector with indices of administered items.
#' @param available Vector with indices of available items.
#' @param number_dimensions Number of dimensions of theta.
#' @param prior_form String indicating the form of the prior; one of \code{"normal"} or \code{"uniform"}.
#' @param prior_parameters List containing mu and Sigma of the normal prior: \code{list(mu = ..., Sigma = ...)}, or 
#' the upper and lower bound of the uniform prior: \code{list(lower_bound = ..., upper_bound = ...)}.
#' The list element \code{Sigma} should always be in matrix form. List elements \code{mu}, \code{lower_bound}, and \code{upper_bound} should always be vectors.
#' The length of \code{mu}, \code{lower_bound}, and \code{upper_bound} should be equal to the number of dimensions.
#' For uniform prior in combination with expected aposteriori estimation, true theta should fall within 
#' \code{lower_bound} and \code{upper_bound} and be not too close to one of these bounds, in order to prevent errors. 
#' @param number_itemsteps_per_item Vector containing the number of non missing cells per row of the beta matrix.
#' @inheritParams shadowcat
#' @return Vector with PEKL information for each yet available item.
get_posterior_expected_kl_information <- function(estimate, model, answers, administered, available, number_dimensions, estimator, alpha, beta, guessing, prior_form, prior_parameters, number_itemsteps_per_item, eap_estimation_procedure = "riemannsum") {
  result <- function() {
    # we'll perform a very basic integration over the theta range
    # expand the grid for multidimensional models (number of calculations will be length(theta_values)**Q, which can still get quite high for high dimensionalities.)
    probabilities_given_eap_estimate <- get_probs_and_likelihoods_per_item(get_theta_estimate(), model, get_subset(alpha, available), get_subset(beta, available), get_subset(guessing, available), with_likelihoods = FALSE)
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
    probabilities_given_theta <- get_probs_and_likelihoods_per_item(theta, model, get_subset(alpha, available), get_subset(beta, available), get_subset(guessing, available), with_likelihoods = FALSE)
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

