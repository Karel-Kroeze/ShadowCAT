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
#' @importFrom mvtnorm dmvnorm
#' @export
prob <- function(test, person = NULL, theta = NULL, deriv = FALSE, prior = NULL){
  
  # TODO: Check input.
  
  # attach items directly
  items <- test$items
  
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
  
  # if test does not require a theta, but it is set, raise a warning (but do use prior).
  if ( ! is.null(prior) && test$estimator == "ML" && deriv == TRUE){
    warning("Prior set for ML estimator - this is not standard!")
  }
  
  # logistic function
  lf <- function(x){ exp(x)/(1+exp(x)) }
  
  # set up output matrix
  out <- list()
  P <- matrix(NA,items$K,items$M+1) # max categories + 1 for false answer.
  
  # simplify input
  m <- items$pars$m
  a <- items$pars$alpha
  b <- items$pars$beta
  c <- items$pars$guessing
  K <- items$K
  
  # compensatory - inner product of alpha * theta.
  at <- a %*% matrix(theta, ncol = 1)
    
  # simplify input for derivatives.
  if (deriv){
    u <- person$responses
    Q <- items$Q
    l <- d <- D <- numeric(K)
  }
  
  # Three Paramater Logistic (MultiDim) (Segall, 1996)
  if(items$model=="3PLM"){
    aux <- numeric(K)
    for (i in 1:K){
      aux[i] <- -a[i,] %*% (theta - b[i,])
    }
    P[,2] <- c + (1-c)/(1+exp(aux))
    P[,1] <- 1 - P[,2]
    
    if(deriv){
      # Segall (MCAT book, 1996, p.71)
      q <- P[,1]
      p <- P[,2]
      
      l <- p^u * q^(1-u)
      d <- ((p-c) * (u-p)) / ((1-c) * p)
      D <- (q *(p-c)*(c*u-p^2)) / (p^2*(1-c)^2)
    }
  }
  
  # graded response model (Samejima, 1969)
  if(items$model=="GRM"){
    for(i in 1:K){
      Psi <- c(1,lf(at[i]-b[i,1:m[i]]),0)
      P[i,1:(m[i]+1)] <- Psi[1:(m[i]+1)] - Psi[2:(m[i]+2)]
    }
    
    if(deriv){
      # Graded Response Model (Glas & Dagohoy, 2007)
      for(i in 1:K){
        j <- u[i]+1 # no more messing about with 0 based indices.
        l[i] <- P[i,j]
        d[i] <- 1 - Psi[j] - Psi[j+1]
        D[i] <- -(Psi[j] * (1-Psi[j]) + Psi[j+1] * (1-Psi[j+1]))
      }
    }    
  }
  
  # Sequential Model (Tutz, 1990)
  if(items$model=="SM"){
    for(i in 1:K){
      Psi <- c(1,lf(at[i]-b[i,1:m[i]]),0)
      for(j in 1:(m[i]+1)){
        P[i,j] <- prod(Psi[1:j]) * (1 - Psi[j+1])
      }
    }
    
    if (deriv){
      # TODO: Wording in Glas & Dagohoy for dij (SM) is odd. 
      for(i in 1:K){
        j <- u[i]+1 # no more messing about with zero-based indeces
        l[i] <- P[i,j]
        
        # Hacky solution to false answers.
        # NOTE: if u_i == 0, only the -Psi(i(h+1)) term remains.
        if (j == 1) {
          aux = 0
        } else {
          aux = sum(1 - Psi[2:j])
        }
        d[i] <- aux - Psi[j+1]
        D[i] <- -sum(Psi[2:(j+1)] * (1 - Psi[2:(j+1)]))
      }
    }
  }
  
  # Generalised Partial Credit Model (Muraki, 1992)
  if(items$model=="GPCM"){
    for(i in 1:K){
      aux <- exp((1:m[i]) * at[i] - b[i,1:m[i]])
      aux2 <- 1 + sum(aux)
      P[i,2:(m[i]+1)] <- aux / aux2
      P[i,1] <- 1-sum(P[i,],na.rm=TRUE)
    }
    
    if (deriv) {
      # TODO: again, why the extra parentheses? 
      for(i in 1:K){
        mi <- 1:m[i]
        pi <- P[i,mi+1] # remove j = 0, index is now also correct.
        mp <- sum(mi*pi)
        
        l[i] <- P[i,u[i]+1]
        d[i] <- u[i] - mp
        D[i] <- -sum((mi * pi) * (mi - mp))
      }
    }
  }
  
  if (deriv){
    # create (log)likelihood (L, LL)
    LL <- sum(log(l))
    
    # CEES: derivatives are correct for a single item, but not for K > 1?
    # create derivatives
    d1 <- matrix(d, nrow = 1) %*% a
    
    # d2
    d2 <- matrix(0,Q,Q)
    for (i in 1:K){
      d2 <- d2 + a[i,] %*% t(a[i,]) * D[i]
    }
    
    # prior
    if ( ! is.null(prior)) {
      # Alleen variabele deel van multivariaat normaal verdeling (exp).
      LL <- LL - (t(theta) %*% solve(prior) %*% theta) / 2
      d1 <- d1 - t(solve(prior) %*% theta)
      d2 <- d2 - solve(prior)
    }
      
    # attach to output
    out$LL <- LL # log likelihood
    out$d <- d   # individual d terms (befosre summing over alpha)
    out$d1 <- d1 # first derivative of complete set of items at theta (+prior)
    out$D <- D   # individual D terms (before summing over alpha*alpha`)
    out$d2 <- d2 # second derivative of complete set of items at theta (+prior)
  }
  
  out$P <- P
  return(invisible(out))
}

#' Log Likelihood
#' 
#' Internal - fetch appropriate elements from prob for nlm optimizer
#' @param theta
#' @param test
#' @param person
#' @param should values be reversed (useful for minimization, reverses LL as well as derivatives)
#' @return Log-Likelihood, as well as gradient and hessian attributes.
#' @importFrom stats nlm
LL <- function(theta, test, person, minimize = FALSE) {
  # subset items that have a response
  items <- subset(test$items, person$administered)
  
  # get LL and derivatives.
  PROB <- prob(test, person, theta, deriv = TRUE)
  
  # prepare output
  out <- PROB$LL * (-1) ^ minimize
  # used in nlm, buggy!
  # attr(out, "gradient") <- PROB$d1 * (-1) ^ minimize
  # attr(out, "hessian") <- PROB$d2 * (-1) ^ minimize
  
  # return
  return(invisible(out))
}
