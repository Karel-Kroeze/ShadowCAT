### probability and derivatives

#' Probability and Derivatives for a given (subset of) ShadowCAT.items itembank.
#' 
#' Generic
#' 
#' Details
#' @param items
#' @param person Person object, will use current theta estimate of the person to obtain P values.
#' @param theta Give a specific theta vector to use, overrides person if set.
#' @param deriv Should we fetch derivatives and LL?
#' @param prior If not NULL, prior to be applied to derivatives. 
#' If NULL, prior from person object will be used. 
#' If not set in the person object, a standard normal prior is used in tests with MAP and EAP estimators.
#' Set to FALSE to force not use a prior
#' @importFrom mvtnorm dmvnorm
#' @importFrom Rcpp evalCpp
#' @useDynLib ShadowCAT
#' @export
prob <- function(test, person = NULL, theta = NULL, deriv = FALSE, prior = NULL, items = NULL){
  
  # TODO: Validate input.
  
  # attach item(s) directly
  if (is.null(items)){
    items <- test$items
  } else {
    items <- items
  }
  
  # Make sure a theta exists
  if (is.null(person) && is.null(theta)) {
    warning("No person or theta given - defaulting to 0 vector.")
    theta <- rep(0,items$Q)
  } 
  
  # if theta is not overriden, use the current estimate.
  if (is.null(theta)) { 
    theta <- person$estimate
  }
  
  # if test requires a prior, and it is not set explicitly, fetch from person. 
  # If that isn't set, raise a warning and use a multivariate standard normal
  if (is.null(prior) && test$estimator %in% c("MAP", "EAP") && deriv == TRUE) {
    if (is.null(person)) {
      warning("No person or prior set - defaulting to multivariate standard normal")
      prior <- diag(items$Q)
    } else {
      prior <- person$prior
    }
  }
  
  # if test does not require a prior, but it is set, raise a warning (but do use prior).
  if ( ! is.null(prior) && test$estimator == "ML" && deriv == TRUE){
    warning("Prior set for ML estimator - this is not standard!")
  }
  
  # simplify input
  a <- items$pars$alpha
  b <- items$pars$beta
  c <- items$pars$guessing
  if (deriv) u <- person$responses
  else u <- c(1,2,3,NA,NA,NA,NA,NA) # Batman! (in all seriousness, it just needs to be a numeric vector with one or more NA's).
  
  res <- switch(items$model,
                "3PLM" = PROB_3PLM(theta, a, b, c, u),
                "GRM" = PROB_GRM(theta, a, b, u),
                "SM" = PROB_SM(theta, a, b, u),
                "GPCM" = PROB_GPCM(theta, a, b, u))
  
  out <- list()
  out$P <- res$P
  
  # likelihoods can never truly be zero, let alone negative
  # TODO: make debug output toggleable
  # if (any(out$P <= 0)) cat("\nProbability <= 0 (k =", length(person$responses), ", estimate = ", paste0(round(person$estimate, 2), collapse = ", "), ").")
  out$P[which(out$P <= 0)] <- 1e-10
  
  
  if (deriv){
    # likelihoods can never truly be zero, let alone negative
    res$l[which(res$l <= 0e-10)] <- 1e-10
    
    # create (log)likelihood (L, LL)
    ll <- log(res$l)
    LL <- sum(ll)
    
    # create derivatives
    d1 <- matrix(res$d, nrow = 1) %*% a
    
    # d2
    d2 <- matrix(0,items$Q,items$Q)
    for (i in seq_along(res$D)){
      d2 <- d2 + a[i,] %*% t(a[i,]) * res$D[i]
    }
    
    # prior not null, not false.
    if ( ! is.null(prior) && prior ) {
      # TODO: mean? Prior currently only allows a variable covariance matrix, mean is fixed at 0 (vector)
      # Alleen variabele deel van multivariaat normaal verdeling (exp).
      LL <- LL - (t(theta) %*% solve(prior) %*% theta) / 2
      d1 <- d1 - t(solve(prior) %*% theta)
      d2 <- d2 - solve(prior)
    }
    
    # attach to output
    out$LL <- LL # log likelihood
    out$d1 <- d1 # first derivative of complete set of items at theta (+prior)
    out$d2 <- d2 # second derivative of complete set of items at theta (+prior)
  }
  
  return(invisible(out))
}

#' Log Likelihood
#' 
#' Internal - fetch appropriate elements from prob for nlm optimizer
#' @param theta
#' @param test
#' @param person
#' @param should values be reversed (useful for minimization, reverses LL as well as derivatives)
#' @param posterior If FALSE (default), the result is the Log-Likelihood. 
#' If inclusion of a prior is desired (e.g. for MAP estimates), set to TRUE.
#' @return Log-Likelihood, as well as gradient and hessian attributes.
#' @importFrom stats nlm
LL <- function(theta, test, person, minimize = FALSE, log = TRUE, posterior = FALSE ){
  # subset items that have a response
  test$items <- subset(test$items, person$administered)
  
  # get LL and derivatives.
  if (posterior)
    PROB <- prob(test, person, theta, deriv = TRUE) # prior will be taken from person object if prior = NULL (default).
  else
    PROB <- prob(test, person, theta, deriv = TRUE, prior = FALSE)
  
  # prepare output
  out <- PROB$LL * (-1) ^ minimize
  
  # gradient and hessian for nlm optimizer.
  attr(out, "gradient") <- PROB$d1 * (-1) ^ minimize
  attr(out, "hessian") <- PROB$d2 * (-1) ^ minimize
  
  # return
  if ( ! log) return(invisible(exp(out)))  
  return(invisible(out))
}
