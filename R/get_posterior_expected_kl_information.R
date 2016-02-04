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
#' @param responses vector with person responses
#' @param administered vector with indeces of administered items
#' @param available vector with indeces of available items
#' @param number_dimensions number of dimensions
#' @param estimator type of estimator to be used, one of "maximum_aposteriori", "maximum_likelihood", or "expected_aposteriori"
#' @param alpha matrix containing the alpha parameters
#' @param beta matrix containing the beta parameters
#' @param guessing matrix containing the quessing parameters
#' @param prior prior covariance matrix for theta
#' @param number_itemsteps_per_item vector containing the number of non missing cells per row of the beta matrix
#' @param lower_bound vector with lower bounds for theta per dimension; estimated theta values smaller than the lowerbound values are truncated to the lowerbound values
#' @param upper_bound vector with upper bounds for theta per dimension; estimated theta values larger than the upperbound values are truncated to the upperbound values
#' @param theta_range Vector of theta values to be evaluated in the numerical integration. Using a sparser range may alleviate stress in higher dimensional tests.
#' @return Vector with PEKL information for each yet available item.
#' @export
get_posterior_expected_kl_information <- function(estimate, model, responses, administered, available, number_dimensions, estimator, alpha, beta, guessing, prior, number_itemsteps_per_item, lower_bound, upper_bound, theta_range = -3:3) {
  result <- function() {
    # we'll perform a very basic integration over the theta range
    # expand the grid for multidimensional models (number of calculations will be length(theta_range)**Q, which can still get quite high for high dimensionalities.)
    grid <- expand.grid(rep(list(theta_range), number_dimensions))

    # sum over values of the grid
    row_or_vector_sums(apply(grid, 1, kullback_leibler_divergence, theta0 = get_theta_estimate()))
  }
  
  get_theta_estimate <- function() {
    # collect expected a posteriori estimate
    if (estimator == "expected_aposteriori")
      estimate
    else
      estimate_latent_trait(estimate, responses, prior, model, administered, number_dimensions, estimator = "expected_aposteriori", alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound)
  }
  
  #' Kullback Leibler Divergence for given items and pairs of thetas x posterior density.
  #' theta0 is theta estimated with expected_aposteriori
  #' returns vector containing information for each yet available item
  kullback_leibler_divergence <- function(theta, theta0) {
    # TODO: wrap this into PEKL, do not recompute P0 for each theta (considering it is constant for the current posterior).
    probabilities_given_theta <- get_probs_and_likelihoods_per_item(theta, model, get_subset(alpha, available), get_subset(beta, available), get_subset(guessing, available), with_likelihoods = FALSE)$P
    probabilities_given_theta0 <- get_probs_and_likelihoods_per_item(theta0, model, get_subset(alpha, available), get_subset(beta, available), get_subset(guessing, available), with_likelihoods = FALSE)$P
    likelihood_or_post_density_theta <- likelihood_or_post_density(theta, responses, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior = prior, return_log_likelihood_or_post_density = FALSE)
    
    rowSums(probabilities_given_theta0 * (log(probabilities_given_theta0) - log(probabilities_given_theta)), na.rm = TRUE) * likelihood_or_post_density_theta
  }
  
  result()
}

