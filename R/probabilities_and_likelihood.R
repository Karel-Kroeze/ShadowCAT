#' Probabilities, likelihood, and derivatives for a given (subset of) items.
#' 
#' @param theta vector with true or estimated theta
#' @param responses vector with person responses to the administered items; only required if output is "likelihood" or "both"
#' @param model string, one of '3PLM', 'GPCM', 'SM' or 'GRM', for the three-parameter logistic, generalized partial credit, sequential or graded response model respectively.
#' @param items_to_include vector with indeces of items to include in computations, usually indeces of either the administered items or all items
#' @param number_dimensions number of dimensions
#' @param estimator type of estimator to be used, one of "maximum_aposteriori", "maximum_likelihood", or "expected_aposteriori"
#' @param alpha matrix containing the alpha parameters
#' @param beta matrix containing the beta parameters
#' @param guessing matrix containing the quessing parameters
#' @param prior prior covariance matrix for theta; only required if estimator is "maximum_aposteriori" or "expected_aposteriori" and output is "likelihood" or "both"
#' @param return_log_likelihood if TRUE, log of likelihood is returned, else likelihood
#' @param inverse_likelihood should likelihood values be reversed (useful for minimization, reverses LL as well as derivatives)
#' @param output string, one of "probs" (return matrix with probabilities for each included item), "likelihood" (return likelihood or posterior density of theta, with first and second derivatives as attributes),
#' or "both" (return list containing both)
#' @return if output = "probs": matrix with for each included item (rows) the probability of scoring in each answer category (columns), given theta
#' if output = "likelihood": the likelihood (estimator = maximum_likelihood) or posterior density (maximum_aposteriori/expected_aposteriori) of theta, with first and second derivatives as attributes
#' if output = "both": al list containing both
#' @importFrom mvtnorm dmvnorm
#' @importFrom Rcpp evalCpp
#' @useDynLib ShadowCAT
#' @export
probabilities_and_likelihood <- function(theta, responses = NULL, model, items_to_include, number_dimensions, estimator, alpha, beta, guessing, prior = NULL, return_log_likelihood = TRUE, inverse_likelihood = FALSE, output = "probs") {
  # TODO: Check input.
  # TODO priors: mean? 
  # priors: Alleen variabele deel van multivariaat normaal verdeling (exp).
  number_items <- length(items_to_include)
  alpha <- get_subset(alpha, items_to_include)
  beta <- get_subset(beta, items_to_include)
  guessing <- get_subset(guessing, items_to_include)

  result <- function() {
    probabilities <- get_probabilities()
    switch(output,
           probs = return_probs(probabilities),
           likelihood = return_likelihood(probabilities),
           both = return_both(probabilities))
  }

  return_probs <- function(probabilities) {
    probabilities$P
  }
  
  return_likelihood <- function(probabilities) {
    log_likelihood <- get_log_likelihood(probabilities) * (-1) ^ inverse_likelihood
    attr(log_likelihood, "gradient") <- get_first_derivative(probabilities) * (-1) ^ inverse_likelihood
    attr(log_likelihood, "hessian") <- get_second_derivative(probabilities) * (-1) ^ inverse_likelihood
    if (return_log_likelihood)  
      log_likelihood
    else
      exp(log_likelihood)
  }

  return_both <- function(probabilities) {
    probs_and_likelihood <- list(probabilities = probabilities$P,
                                 likelihood = if (return_log_likelihood) 
                                                get_log_likelihood(probabilities) * (-1) ^ inverse_likelihood
                                              else
                                                exp(get_log_likelihood(probabilities) * (-1) ^ inverse_likelihood))
    attr(probs_and_likelihood$likelihood, "gradient") <- get_first_derivative(probabilities) * (-1) ^ inverse_likelihood
    attr(probs_and_likelihood$likelihood, "hessian") <- get_second_derivative(probabilities) * (-1) ^ inverse_likelihood
    probs_and_likelihood
  }

  get_responses <- function() {
    if (output != "probs") 
      responses
    else
      numeric(0)
  }
  
  get_probabilities <- function() {
    responses <- get_responses()
    probs <- switch(model,
                    "3PLM" = PROB_3PLM(theta, alpha, beta, guessing, responses, output != "probs"),
                    "GRM" = PROB_GRM(theta, alpha, beta, responses, output != "probs"),
                    "SM" = PROB_SM(theta, alpha, beta, responses, output != "probs"),
                    "GPCM" = PROB_GPCM(theta, alpha, beta, responses, output != "probs"))
    
    # likelihoods can never truly be zero, let alone negative
    probs$P[which(probs$P <= 0)] <- 1e-10
    probs
  }
  
  get_log_likelihood <- function(probabilities) {
    # likelihoods can never truly be zero, let alone negative
    probabilities$l[which(probabilities$l <= 0)] <- 1e-10  
    log_likelihood <- sum(log(probabilities$l))
    if (estimator == "maximum_likelihood")
      log_likelihood
    else
      log_likelihood - (t(theta) %*% solve(prior) %*% theta) / 2    
  }
  
  # TODO: derivatives are correct for a single item, but not for K > 1?
  get_first_derivative <- function(probabilities) {
    derivative1 <- matrix(probabilities$d, nrow = 1) %*% alpha
    if (estimator == "maximum_likelihood")
      derivative1
    else
      derivative1 - t(solve(prior) %*% theta)  
  }
  
  get_second_derivative <- function(probabilities) {
    derivative2 <- sum_loop_outputs(start_object = matrix(0, number_dimensions, number_dimensions), 
                                    loop_vector = 1:number_items, 
                                    FUN = function(item, alpha, D) { alpha[item,] %*% t(alpha[item,]) * D[item] }, 
                                    alpha = alpha, 
                                    D = probabilities$D)
    if (estimator == "maximum_likelihood")
      derivative2
    else
      derivative2 - solve(prior)  
  }
  
  validate <- function() {
    if (is.null(prior) && estimator %in% c("maximum_aposteriori", "expected_aposteriori") && output != "probs")
      add_error("prior", "is missing but required for estimate")
  }
    
  invalid_result <- function() {
    list(P = NA,
         LL = NA,
         d1 = NA,
         d2 = NA,
         errors = errors())
  }
    
  validate_and_run()
}
