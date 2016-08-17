#' Likelihood or posterior density
#'
#' Get Likelihood or posterior density (with normal prior) over all included items, 
#' and derivatives, for a given set of answers to items.
#' 
#' @param theta Vector with true or estimated theta.
#' @param answers Vector with answers to the administered items.
#' @param items_to_include Vector with indices of items to which answers have been given.
#' @param number_dimensions Number of dimensions of theta.
#' @param number_itemsteps_per_item Vector containing the number of non missing cells per row of the beta matrix.
#' @param prior_parameters List containing mu and Sigma of the normal prior: \code{list(mu = ..., Sigma = ...)}.
#' Sigma should always be in matrix form.
#' @param return_log_likelihood_or_post_density If \code{TRUE}, log of likelihood or posterior density is returned, else likelihood or posterior density on original scale.
#' @param inverse_likelihood_or_post_density If \code{TRUE}, likelihood or posterior density value is reversed (useful for minimization, also reverses derivatives).
#' @param with_derivatives If \code{TRUE}, first and second derivatives are added to the return value as attributes.
#' @inheritParams shadowcat
#' @return The likelihood (estimator is maximum_likelihood) or posterior density with normal prior (estimator is not maximum_likelihood) of theta. 
#' If requested, first and second derivatives are added as attributes.
likelihood_or_post_density <- function(theta, answers = NULL, model, items_to_include, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, prior_parameters = NULL, return_log_likelihood_or_post_density = TRUE, inverse_likelihood_or_post_density = FALSE, with_derivatives = TRUE) {
  number_items <- length(items_to_include)
  alpha <- get_subset(alpha, items_to_include)
  beta <- get_subset(beta, items_to_include)
  guessing <- get_subset(guessing, items_to_include)
  number_itemsteps_per_item <- get_subset(number_itemsteps_per_item, items_to_include)
  
  result <- function() {
    probs_and_likelihoods <- get_probs_and_likelihoods_per_item(theta = theta, model = model, alpha = alpha, beta = beta, guessing = guessing, number_dimensions = number_dimensions, number_items = number_items, number_itemsteps_per_item = number_itemsteps_per_item, answers = answers, with_likelihoods = TRUE, with_derivatives = with_derivatives) 
    log_likelihood_or_post_density <- get_log_likelihood_or_post_density(probs_and_likelihoods) * (-1) ^ inverse_likelihood_or_post_density
    if (with_derivatives) {
      attr(log_likelihood_or_post_density, "gradient") <- get_first_derivative(probs_and_likelihoods) * (-1) ^ inverse_likelihood_or_post_density
      attr(log_likelihood_or_post_density, "hessian") <- get_second_derivative(probs_and_likelihoods) * (-1) ^ inverse_likelihood_or_post_density
    }
    if (return_log_likelihood_or_post_density)  
      log_likelihood_or_post_density
    else
      exp(log_likelihood_or_post_density)
  }
  
  validate <- function() {
    if (is.null(prior_parameters) && estimator %in% c("maximum_aposteriori", "expected_aposteriori"))
      add_error("prior_parameters", "is missing")
  }
  
  invalid_result <- function() {
    list(errors = errors())
  }
  
  get_log_likelihood_or_post_density <- function(probs_and_likelihoods) {
    log_likelihood <- sum(log(probs_and_likelihoods$likelihoods))
    if (estimator == "maximum_likelihood")
      log_likelihood
    else
      log_likelihood - (t(theta - prior_parameters$mu) %*% solve(prior_parameters$Sigma) %*% (theta - prior_parameters$mu)) / 2    
  }
  
  get_first_derivative <- function(probs_and_likelihoods) {
    derivative1 <- matrix(probs_and_likelihoods$derivatives1, nrow = 1) %*% alpha
    if (estimator == "maximum_likelihood")
      derivative1
    else
      derivative1 - t(theta) %*% solve(prior_parameters$Sigma) + t(prior_parameters$mu) %*% solve(prior_parameters$Sigma)  
  }
  
  get_second_derivative <- function(probs_and_likelihoods) {
    derivative2 <- sum_loop_outputs(start_object = matrix(0, number_dimensions, number_dimensions), 
                                    loop_vector = 1:number_items, 
                                    FUN = function(item, alpha, D) { alpha[item,] %*% t(alpha[item,]) * D[item] }, 
                                    alpha = alpha, 
                                    D = probs_and_likelihoods$derivatives2)
    if (estimator == "maximum_likelihood")
      derivative2
    else
      derivative2 - solve(prior_parameters$Sigma)  
  }
  
  validate_and_run()
}
