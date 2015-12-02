### probability and derivatives

#' Probability and Derivatives for a given (subset of) ShadowCAT.items itembank. It gives the probability of 
#' responses given theta (only useful for simulation I think), and the likelihood of given responses
#' 
#' Generic
#' 
#' Details
#' @param theta true or estimated theta
#' @param responses person responses to the administered items; only required if output is "likelihoods" or "both"
#' @param model string, one of '3PLM', 'GPCM', 'SM' or 'GRM', for the three-parameter logistic, generalized partial credit, sequential or graded response model respectively.
#' @param items_to_include indeces of items to include in computations, usually indeces of either the administered items or all items
#' @param number_dimensions number of dimensions
#' @param estimator type of estimator to be used, one of "MAP" (Maximum a posteriori estimation), "EAP" (Expected A Posteriori Estimation), or "ML" (maximum likelihood)
#' @param alpha matrix containing the alpha parameters
#' @param beta matrix containing the beta parameters
#' @param quessing matrix containing the quessing
#' @param prior prior covariance matrix for theta; only required if estimator is "MAP" or "EAP" and output is "likelihoods" or "both"
#' @param return_log_likelihoods if TRUE, log of likelihoods are returned, else likelihoods
#' @param inverse_likelihoods should likelihood values be reversed (useful for minimization, reverses LL as well as derivatives)
#' @param output string, one of "probs" (return vector of probabilities only), "likelihoods" (return vector of likelihoods with first and second derivatives as attributes),
#' or "both" (return list of both probabilities and likelihoods with first and second derivatives as attributes)
#' @return probabilities and/or likelihoods, depending on output parameter
#' @importFrom mvtnorm dmvnorm
#' @importFrom Rcpp evalCpp
#' @useDynLib ShadowCAT
#' @export
probabilities_and_likelihoods <- function(theta, responses = NULL, model, items_to_include, number_dimensions, estimator, alpha, beta, guessing, prior = NULL, return_log_likelihoods = TRUE, inverse_likelihoods = FALSE, output = "probs") {
  # TODO: Check input.
  # TODO priors: mean? 
  # priors: Alleen variabele deel van multivariaat normaal verdeling (exp).
  number_items <- length(items_to_include)
  alpha <- get_subset(alpha, items_to_include)
  beta <- get_subset(beta, items_to_include)
  quessing <- get_subset(guessing, items_to_include)
  
  result <- function() {
    probabilities <- get_probabilities()
    switch(output,
           probs = return_probs(probabilities),
           likelihoods = return_likelihoods(probabilities),
           both = return_both(probabilities))
  }

  return_probs <- function(probabilities) {
    probabilities$P
  }
  
  return_likelihoods <- function(probabilities) {
    log_likelihoods <- get_log_likelihoods(probabilities) * (-1) ^ inverse_likelihoods
    attr(log_likelihoods, "gradient") <- get_first_derivative(probabilities) * (-1) ^ inverse_likelihoods
    attr(log_likelihoods, "hessian") <- get_second_derivative(probabilities) * (-1) ^ inverse_likelihoods

    if (return_log_likelihoods)  
      log_likelihoods
    else
      exp(log_likelihoods)
  }

  return_both <- function(probabilities) {
    out <- list(probabilities = probabilities$P,
                likelihoods = if (return_log_likelihoods) 
                                get_log_likelihoods(probabilities) * (-1) ^ inverse_likelihoods
                              else
                                exp(get_log_likelihoods(probabilities) * (-1) ^ inverse_likelihoods))
    attr(out$likelihoods, "gradient") <- get_first_derivative(probabilities) * (-1) ^ inverse_likelihoods
    attr(out$likelihoods, "hessian") <- get_second_derivative(probabilities) * (-1) ^ inverse_likelihoods
    out
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
                    "GRM" = PROB_GRM(theta,alpha, beta, responses, output != "probs"),
                    "SM" = PROB_SM(theta, alpha, beta, responses, output != "probs"),
                    "GPCM" = PROB_GPCM(theta, alpha, beta, responses, output != "probs"))
    
    # likelihoods can never truly be zero, let alone negative
    # TODO: make debug output toggleable
    # if (any(out$P <= 0)) cat("\nProbability <= 0 (k =", length(person$responses), ", estimate = ", paste0(round(person$estimate, 2), collapse = ", "), ").")
      
    probs$P[which(probs$P <= 0)] <- 1e-10
    probs
  }
  
  get_log_likelihoods <- function(probabilities) {
    # likelihoods can never truly be zero, let alone negative
    probabilities$l[which(probabilities$l <= 0)] <- 1e-10  
    log_likelihoods <- sum(log(probabilities$l))
    if (estimator == "ML")
      log_likelihoods
    else
      log_likelihoods - (t(theta) %*% solve(prior) %*% theta) / 2    
  }
  
  # TODO: derivatives are correct for a single item, but not for K > 1?
  get_first_derivative <- function(probabilities) {
    derivative1 <- matrix(probabilities$d, nrow = 1) %*% alpha
    if (estimator == "ML")
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
    if (estimator == "ML")
      derivative2
    else
      derivative2 - solve(prior)  
  }
  
  validate <- function() {
    if (!is.null(prior) && estimator == "ML" && output != "probs")
      add_error("prior", "set for ML estimator, this makes no sense")
    if (is.null(prior) && estimator %in% c("MAP", "EAP") && output != "probs")
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
