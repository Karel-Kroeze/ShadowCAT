#' Fisher Information
#' 
#' Fisher Information (expected information) per item.
#' 
#' The expected information is the sum of second derivatives, weighted by their probability.
#' @param items
#' @param person
#' @return array with an information matrix for each item (QxQxK).
#' @export
FI <- function(test, person) {
  # Fisher Information
  # minus the expectation of the second derivative of the log-likelihood
  # Expectation -> sum of derivatives for each category, 'weighted' by their probability.
  # Notice that the minuses in individual D terms are there - they seem like they shouldn't be, but the results do match up?
  
  # simplify input
  model <- test$items$model
  p <- prob(test$items, person)$P
  a <- test$items$pars$alpha
  b <- test$items$pars$beta
  K <- test$items$K
  m <- test$items$pars$m
  D <- numeric(K)
  Q <- test$items$Q
  
  # define the logistic function
  lf <- function(x) exp(x)/(1+exp(x))
  
  # compensatory -> inner product of alpha * theta
  at <- apply(a * drop(person$estimate),1,sum)
  
  if(model == "3PLM"){
    # exact form Segall 1997, CAT book, p. 72
    c <- test$items$pars$guessing
    q <- p[,1]; p <- p[,2]
    D <- (q/p) * ((p-c)/(1-c))^2
  }
  
  if (model == "GRM"){
    # Graded Response Model (Glas & Dagohoy, 2007)
    for(i in 1:K){
      for(j in 1:(m[i]+1)){
        Psi <- c(1,lf(at[i]-b[i,1:m[i]]),0)
        D[i] <- D[i] + p[i,j] * -(Psi[j] * (1-Psi[j]) + Psi[j+1] * (1-Psi[j+1]))
      }
    }
  }
  
  if (model == "SM") {
    # Sequential Model (Tutz, xxxx)
    # TODO: triple check this.
    for(i in 1:K){
      for(j in 1:(m[i]+1)){
        Psi <- c(1, lf(at[i] - b[i,1:m[i]]), 0) # basically: 1, lf(at - b), 0.
        D[i] <- D[i] + p[i,j] * -sum(Psi[2:(j+1)] * (1 - Psi[2:(j+1)]))
      }
    }
  }
  
  if (model == "GPCM"){
    # Generalized Partial Credit Model (Muraki, 1992)
    # TODO: again, why the extra parentheses? 
    # article calls for j=0 ~ m, but j=0 term falls off?
    for(i in 1:K){
      mi <- 1:m[i]
      for(j in (1:m[i])+1){ # sequence is offset -> 2:m+1, since false answer is unnecesary. mi accounts for weighting, not j.
        mp <- sum(mi*p[i,j])
        D[i] <- D[i] + p[i,j] * -sum((mi * p[i,j]) * (mi - mp))
      }
    }
  }
  
  # just dump everything back, how to deal with information is up to caller function.
  out <- array(0,c(Q, Q, K))
  for (i in 1:K) {
    out[,,i] <- (a[i,] %*% t(a[i,])) * D[i]
  }
  
  # return
  return(invisible(out))
}