### probability and derivatives

#' Probability and Derivatives for a given (subset of) ShadowCAT.items itembank.
#' 
#' Generic
#' 
#' Details
#' @param test
#' @param person Person object, will use current theta estimate of the person to obtain P values.
#' @param theta Give a specific theta vector to use, overrides person if set.
#' @export
prob <- function(test, person = NULL,theta = NULL){
  if (is.null(person) && is.null(theta)) {
    warning("No person or theta given - defaulting to 0 vector.")
    theta <- rep(0,test$items$Q)
  } 
  if (is.null(theta)) {
    theta <- person$estimate
  }
    # TODO: Check input.
  
  # logistic function
  lf <- function(x){ exp(x)/(1+exp(x)) }
  
  # set up output matrix
  out <- list()
  P <- matrix(NA,test$items$K,test$items$M+1) # max categories + 1 for false answer.
  
  # simplify input
  m <- test$items$pars$m
  a <- test$items$pars$alpha
  b <- test$items$pars$beta
  c <- test$items$pars$guessing
  K <- test$items$K
  
  # only do derivatives when we have info - so after we have administered some items. Post items.
  post = FALSE
  # simplify input for derivatives.
  if ( ! is.null(person) && ! is.null(person$resp)){
    post = TRUE
    u <- person$resp
    Q <- test$items$Q
    l <- d <- D <- numeric(K)
  }

  # compensatory - inner product of alpha * theta.
  at <- apply(a * drop(theta),1,sum)
  
  # Three Paramater Logistic (MultiDim) (Segall, 1996)
  if(test$items$model=="3PLM"){
    aux <- numeric(K)
    for (i in 1:K){
      aux[i] <- -a[i,] %*% (theta - b[i,])
    }
    P[,2] <- c + (1-c)/(1+exp(aux))
    P[,1] <- 1 - P[,2]
    
    if(post){
      # Segall (MCAT book, 1996, p.71)
      q <- P[,1]
      p <- P[,2]
      
      l <- p^u * q^(1-u)
      d <- ((p-c) * (u-p)) / ((1-c) * p)
      D <- (q *(p-c)*(c*u-p^2)) / (p^2*(1-c)^2)
    }
  }
  
  # graded response model (Samejima, 1969)
  if(test$items$model=="GRM"){
    for(i in 1:K){
      Psi <- c(1,lf(at[i]-b[i,1:m[i]]),0)
      P[i,1:(m[i]+1)] <- Psi[1:(m[i]+1)] - Psi[2:(m[i]+2)]
    }
    
    if(post){
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
  if(test$items$model=="SM"){
    for(i in 1:K){
      Psi <- c(1,lf(at[i]-b[i,1:m[i]]),0)
      for(j in 1:(m[i]+1)){
        P[i,j] <- prod(Psi[1:j]) * (1 - Psi[j+1])
      }
    }
    
    if (post){
      # TODO: Wording in Glas & Dagohoy for dij (SM) is odd. So Psi_0 = 1, right? Why 1-Psi_(m+1)=1 and not Psi_(m+1)=0?
      # TODO: also, why is there an extra set of parenthesis in the formula for dij?
      # TODO: Finally, it doesn't work. See comments in estiamte.R.
      for(i in 1:K){
        j <- u[i]+1 # no more messing about with 0 based indices.
        l[i] <- P[i,j]
        d[i] <- sum((1-Psi[2:j]) - Psi[(2:j)+1])
        D[i] <- -sum(Psi[2:(j+1)] * (1 - Psi[2:(j+1)]))
      }
    }
  }
  
  # Generalised Partial Credit Model (Muraki, 1992)
  if(test$items$model=="GPCM"){
    for(i in 1:K){
      aux <- exp((1:m[i]) * at[i] - b[i,1:m[i]])
      aux2 <- 1 + sum(aux)
      P[i,2:(m[i]+1)] <- aux / aux2
      P[i,1] <- 1-sum(P[i,],na.rm=TRUE)
    }
    
    if (post) {
      # TODO: again, why the extra parentheses? 
      for(i in 1:K){
        mi <- 1:m[i]
        pi <- P[i,mi+1] # remove j = 0, index is now also correct.
        mp <- sum(mi*pi)
        
        l <- P[i,u[i]+1]
        d[i] <- u[i] - mp
        D[i] <- -sum((mi * pi) * (mi - mp))
      }
    }
  }
  
  if (post){
    # create (log)likelihood (L, LL)
    LL <- sum(log(l))
    L <- exp(LL)
    
    # create derivatives
    d1 <- apply(a * d,2,sum)
    
    # d2
    d2 <- matrix(0,Q,Q)
    for (i in 1:K){
      d2 <- d2 + a[i,] %*% t(a[i,]) * D[i]
    }
    
    # prior
    if ( ! is.null(prior)) d1 <- d1 - solve(prior)%*%theta
    if ( ! is.null(prior)) d2 <- d2 - solve(prior)
    
    # attach to output
    out$LL <- LL
    out$L <- L
    out$d1 <- d1
    out$d2 <- d2
  }
  
  out$P <- P
  return(invisible(out))
}

#### Everything below here is testing code. TODO: remove!
testit <- function(model="GRM",alpha=1,beta=0){
  theta <- seq(-3,3,length.out = 100)
  item <- initItembank(model, alpha, beta,silent=T)
  p <- matrix(0,100,item$M+1)
  
  for (i in seq_along(theta)){
    p[i,] <- prob(item,theta=theta[i])$P
  }
  
  matplot(theta,p,type='l',main=model)
}

# par(mfrow=c(2,2))
# testit('GPCM',1,matrix(c(1),1))
# testit('GRM',1,matrix(c(1),1))
# testit('SM',1,matrix(c(1),1))
# testit('3PLM',1,1)
