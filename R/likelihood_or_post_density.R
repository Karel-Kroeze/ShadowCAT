#' Get Likelihood or posterior density over all included items, and derivatives, based on a a given set of items
#' 
#' @param theta vector with true or estimated theta
#' @param responses vector with person responses to the administered items
#' @param model string, one of '3PLM', 'GPCM', 'SM' or 'GRM', for the three-parameter logistic, generalized partial credit, sequential or graded response model respectively.
#' @param items_to_include vector with indeces of items to which responses have been given
#' @param number_dimensions number of dimensions
#' @param estimator type of estimator to be used, one of "maximum_aposteriori", "maximum_likelihood", or "expected_aposteriori"
#' @param alpha matrix containing the alpha parameters (for complete test bank)
#' @param beta matrix containing the beta parameters (for complete test bank)
#' @param guessing matrix containing the quessing parameters (for complete test bank)
#' @param prior prior covariance matrix for theta; only required if estimator is "maximum_aposteriori" or "expected_aposteriori" and output is "likelihood" or "both"
#' @param return_log_likelihood_or_post_density if TRUE, log of likelihood or posterior density is returned, else likelihood or posterior density on original scale
#' @param inverse_likelihood_or_post_density should likelihood or posterior density values be reversed (useful for minimization, also reverses derivatives)
#' @return the likelihood (estimator = maximum_likelihood) or posterior density (estimator = maximum_aposteriori or expected_aposteriori) of theta, with first and second derivatives as attributes
#' @export
likelihood_or_post_density <- function(theta, responses = NULL, model, items_to_include, number_dimensions, estimator, alpha, beta, guessing, prior = NULL, return_log_likelihood_or_post_density = TRUE, inverse_likelihood_or_post_density = FALSE) {
  # TODO priors: mean? 
  # priors: Alleen variabele deel van multivariaat normaal verdeling (exp).
  number_items <- length(items_to_include)
  alpha <- get_subset(alpha, items_to_include)
  beta <- get_subset(beta, items_to_include)
  guessing <- get_subset(guessing, items_to_include)
  
  result <- function() {
    probs_and_likelihoods <- get_probs_and_likelihoods_per_item(theta, model, alpha, beta, guessing, responses, TRUE) 
    log_likelihood_or_post_density <- get_log_likelihood_or_post_density(probs_and_likelihoods) * (-1) ^ inverse_likelihood_or_post_density
    attr(log_likelihood_or_post_density, "gradient") <- get_first_derivative(probs_and_likelihoods) * (-1) ^ inverse_likelihood_or_post_density
    attr(log_likelihood_or_post_density, "hessian") <- get_second_derivative(probs_and_likelihoods) * (-1) ^ inverse_likelihood_or_post_density
    if (return_log_likelihood_or_post_density)  
      log_likelihood_or_post_density
    else
      exp(log_likelihood_or_post_density)
  }
  
  validate <- function() {
    if (is.null(prior) && estimator %in% c("maximum_aposteriori", "expected_aposteriori"))
      add_error("prior", "is missing")
  }
  
  invalid_result <- function() {
    list(errors = errors())
  }
  
  get_log_likelihood_or_post_density <- function(probs_and_likelihoods) {
    log_likelihood <- sum(log(probs_and_likelihoods$l))
    if (estimator == "maximum_likelihood")
      log_likelihood
    else
      log_likelihood - (t(theta) %*% solve(prior) %*% theta) / 2    
  }
  
  # TODO: derivatives are correct for a single item, but not for K > 1?
  get_first_derivative <- function(probs_and_likelihoods) {
    derivative1 <- matrix(probs_and_likelihoods$d, nrow = 1) %*% alpha
    if (estimator == "maximum_likelihood")
      derivative1
    else
      derivative1 - t(solve(prior) %*% theta)  
  }
  
  get_second_derivative <- function(probs_and_likelihoods) {
    derivative2 <- sum_loop_outputs(start_object = matrix(0, number_dimensions, number_dimensions), 
                                    loop_vector = 1:number_items, 
                                    FUN = function(item, alpha, D) { alpha[item,] %*% t(alpha[item,]) * D[item] }, 
                                    alpha = alpha, 
                                    D = probs_and_likelihoods$D)
    if (estimator == "maximum_likelihood")
      derivative2
    else
      derivative2 - solve(prior)  
  }
  
  validate_and_run()
}
