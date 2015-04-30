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
    # suppress warnings and errors and do EAP instead.
    person$estimate <- tryCatch( nlm(f = LL, p = person$estimate, test = test, person = person, minimize = TRUE)$estimate, # passed on to LL, reverses polarity.
                         error = function(e) {
                           #message(paste0(test$estimator, " failed, trying EAP estimate."))
                           test$estimator <- "EAP"
                           return(estimate(person, test)$estimate)
                         },
                         warning = function(w) {
                           #message(paste0(test$estimator, " failed, trying EAP estimate."))
                           test$estimator <- "EAP"
                           return(estimate(person, test)$estimate)
                         })
    
    # variance
    # get FI
    # TODO: We should really store info somewhere so we don't have to redo this (when using FI based selection criteria).
    info <- FI(test, person)
    so_far <- apply(info[,,person$administered, drop = FALSE], c(1,2), sum)
  
    #add prior
    if (test$estimator == "MAP") so_far <- so_far + solve(person$prior)
    
    # inverse
    attr(person$estimate, "variance") <- solve(so_far)
  }
  
  if (test$estimator == "EAP"){
    # Multidimensional Gauss-Hermite Quadrature
    Q <- test$items$Q
    # TODO: prior mean is currently fixed at zero, update when/if possible.
    # TODO: allow setting ip through internals argument(s)
    adapt <- NULL
    if (length(person$responses) > 5 & !is.null(attr(person$estimate, 'variance'))) adapt <- list(mu = person$estimate, Sigma = as.matrix(attr(person$estimate, "variance")))
    
    QP <- init.quad(Q = Q,
                    prior = list(mu = rep(0, test$items$Q), Sigma = person$prior),
                    adapt = adapt,
                    ip = switch(Q, 50, 15, 6, 4, 3))
    person$estimate <- eval.quad(FUN = LL, X = QP, test = test, person = person, ...)
  }
  
  # enforce boundaries.
  # TODO: make debug output toggleable
  # if (any(person$estimate > test$upperBound | person$estimate < test$lowerBound)) cat("Estimate outside boundaries (k =", length(person$responses), "estimate =", paste0(round(person$estimate, 2), collapse = ", "), ").\n")
  person$estimate[which(person$estimate > test$upperBound)] <- test$upperBound[which(person$estimate > test$upperBound)]
  person$estimate[which(person$estimate < test$lowerBound)] <- test$lowerBound[which(person$estimate < test$lowerBound)]
  
  return(invisible(person))
}

