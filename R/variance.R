#' Variance
#' 
#' Computes the (asymptotic) variance-covariance matrix for the current estimate
#' 
#' @param person
#' @param test
#' @param ... additional internal arguments
#' @return person
#' @export
variance <- function(person, test, ...) {
  
  if ( test$estimator %in% c("MAP", "ML")){
    info <- FI(test, person)
    total <- apply(info[,,person$administered, drop = FALSE], c(1,2), sum)
    
    if (test$estimator == "MAP") total <- total + solve(person$prior)
    
    person$variance <- try(solve(total))
  } else if (test$estimator == "EAP") {
    
  }
}