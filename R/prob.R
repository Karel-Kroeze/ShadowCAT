### probability and derivatives

#' Probability and Derivatives for a given (subset of) ShadowCAT.items itembank. It gives the probability of 
#' responses given theta (only useful for simulation I think), and the likelihood of given responses
#' 
#' Generic
#' 
#' Details
#' @param test 
#' @param person Person object, will use current theta estimate of the person to obtain P values.
#' @param theta Give a specific theta vector to use. Overrides person$estimate if set.
#' @param deriv Should we fetch derivatives and LL?
#' @param prior prior to be applied to derivatives. Overrides person$prior if set
#' @param items overrides test$items if set
#' @return 
#' @importFrom mvtnorm dmvnorm
#' @importFrom Rcpp evalCpp
#' @useDynLib ShadowCAT
#' @export
prob <- function(test, person = NULL, theta = NULL, deriv = FALSE, prior = NULL, items = NULL){
  # TODO: Check input.
  # TODO priors: mean? 
  # priors: Alleen variabele deel van multivariaat normaal verdeling (exp).
  result <- function() {
    items <- get_items()
    theta <- get_theta(items)
    prior <- get_prior(items)
    probabilities <- get_probabilities(items, theta)
    #return(probabilities$d)
    
    if (!deriv)
      list(P = probabilities$P)
    else
      list(P = probabilities$P,
           LL = get_LL(probabilities, prior, theta),
           d1 = get_first_derivative(probabilities, prior, theta, items),
           d2 = get_second_derivative(probabilities, prior, items))  
  }
  
  get_items <- function() {
    if (is.null(items))
      test$items
    else
      items
  }
  
  get_theta <- function(items) {
    if (is.null(person) && is.null(theta)) {
      warning("No person or theta given - defaulting to 0 vector.")
      return(rep(0,items$Q))
    }
    if (is.null(theta)) 
      person$estimate
    else
      theta
  }
  
  get_prior <- function(items) {
    if (is.null(prior) && is.null(person) && test$estimator %in% c("MAP", "EAP") && deriv == TRUE) {
        warning("No person or prior set - defaulting to multivariate standard normal")
        return(diag(items$Q))
    }
    if (is.null(prior) && !is.null(person) && test$estimator %in% c("MAP", "EAP") && deriv == TRUE)
      person$prior
    else
      prior
  }
  
  get_probabilities <- function(items, theta) {
    # This responses definition can be removed once everything has been restructured
    responses <- if (deriv) 
                   person$responses
                 else
                   numeric(0)

    probs <- switch(items$model,
                    "3PLM" = PROB_3PLM(theta, items$pars$alpha, items$pars$beta, items$pars$guessing, responses, deriv),
                    "GRM" = PROB_GRM(theta, items$pars$alpha, items$pars$beta, responses, deriv),
                    "SM" = PROB_SM(theta, items$pars$alpha, items$pars$beta, responses, deriv),
                    "GPCM" = PROB_GPCM(theta, items$pars$alpha, items$pars$beta, responses, deriv))
    
    # likelihoods can never truly be zero, let alone negative
    # TODO: make debug output toggleable
    # if (any(out$P <= 0)) cat("\nProbability <= 0 (k =", length(person$responses), ", estimate = ", paste0(round(person$estimate, 2), collapse = ", "), ").")
    
    # for the 3PLM model P and l are returned in log space, preventing cancellation and overflow. Need to exploit
    # this furher. For now, turn back to non-log space here
    if (items$model == "3PLM") {
     probs$P = as.matrix(exp(probs$P))
     probs$l = exp(probs$l)
    }
      
    probs$P[which(probs$P <= 0)] <- 1e-10
    probs
  }
  
  get_LL <- function(probabilities, prior, theta) {
    # likelihoods can never truly be zero, let alone negative
    probabilities$l[which(probabilities$l <= 0)] <- 1e-10  
    LL <- sum(log(probabilities$l))
    if (test$estimator == "ML")
      LL
    else
      LL - (t(theta) %*% solve(prior) %*% theta) / 2    
  }
  
  # TODO: derivatives are correct for a single item, but not for K > 1?
  get_first_derivative <- function(probabilities, prior, theta, items) {
    derivative1 <- matrix(probabilities$d, nrow = 1) %*% items$pars$alpha
    if (test$estimator == "ML")
      derivative1
    else
      derivative1 - t(solve(prior) %*% theta)  
  }
  
  get_second_derivative <- function(probabilities, prior, items) {
    derivative2 <- sum_loop_outputs(start_object = matrix(0, items$Q, items$Q), 
                                    loop_vector = 1:items$K, 
                                    FUN = function(item, alpha, D) { alpha[item,] %*% t(alpha[item,]) * D[item] }, 
                                    alpha = items$pars$alpha, 
                                    D = probabilities$D)
    if (test$estimator == "ML")
      derivative2
    else
      derivative2 - solve(prior)  
  }
  
  validate <- function() {
    if (!is.null(prior) && test$estimator == "ML" && deriv == TRUE){
      add_error("prior", "set for ML estimator, this makes no sense")
    }
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
