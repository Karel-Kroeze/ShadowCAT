#' Latent trait estimation
#' 
#' ML, MAP, EAP
#' 
#' 
#' @param person
#' @param test
#' @param ...
#' @return person
#' @importFrom MultiGHQuad init.quad eval.quad
#' @export
estimate <- function(person, test, ...) {
  # catch additional arguments
  args <- list(...)
  
  if (test$estimator %in% c("ML", "MAP")){
    # for now, simple nlm (TODO: look at optim, and possible reintroducing pure N-R).
    # We want a maximum, but nlm produces minima -> reverse function call. 
    # LL is the target function, test, person and minimize need to be passed on. We also want the value of the hessian at the final estimate.
    
    # note that prior is applied in LL (incorrectly it seems, but still).
    minimum <- nlm(f = LL, p = person$estimate, test = test, person = person,
                   minimize = TRUE)
    
    person$estimate <- minimum$estimate
  }
  
  if (test$estimator == "EAP"){
    # Multidimensional Gauss-Hermite Quadrature
    Q <- test$items$Q
    QP <- init.quad(Q = Q, Sigma = person$prior, ip = max(ceiling(20/(Q*2)),3))
    person$estimate <- eval.quad(FUN = LL, X = QP, test = test, person = person)
  }
  
  return(invisible(person))
}
