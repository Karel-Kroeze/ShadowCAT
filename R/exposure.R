#' exposure control
#' 
#' Sympson-Hetter style exposure control. Van der Linden & Veldkamp, 2004.
#' 
#' This is half-implemented. The routines are here, but there is no mechanism for tracking exposure over time. 
#' Administration of exposure will depend on the test setting, and is left to the user. 
#' 
#' ONLY WORKS WITH SHADOW TESTING!
#' 
#' @param test
#' @param exposure vector of actual exposure counts.    
#' @param eligible vector of eligibility counts (including lifting of restrictions.)
#' @param feasible count of respondents for whom a feasible solution was found.
#' @param total count of total respondents so far.
#' @param target (vector) target exposure rate(s). If scalar, will be used for all.
#' @return test
#' @export
#' @importFrom lpSolve lp
createExposureConstraint <- function(test, exposure, eligible, feasible, total, target) {
  
  # probability of being eligible
  PE <- pmin(1 - total / feasible + ( total * eligible * target ) / ( exposure * feasible), 1)
  # early on will be 0/0 -> NaN 
  PE[is.nan(PE)] <- 1
  
  # toss a coin
  ineligible <- as.integer(runif(length(PE)) > PE)
  
  # we need to check if this solution is feasible.
  if(lp(direction = 'max',
        const.mat = as.matrix(cbind(test$constraints$lp_chars, ineligible)),
        const.dir = c(test$constraints$constraints$op, '='),
        const.rhs = c(test$constraints$constraints$target, 0),
        objective.in = rep(1, test$items$K),
        all.bin = TRUE,
        transpose.constraints = FALSE)$status)
  {
    # if status != 0, there was no feasible solution, therefore we shall remove the ineligibility constraints
    ineligible <- rep(0, test$items$K)
    feasible <- FALSE
  } else {
    feasible <- TRUE
  }
  
  # add constraints
  test$constraints$lp_chars <- cbind(test$constraints$lp_chars, ineligible = ineligible)
  test$constraints$constraints <- rbind(test$constraints$constraints, list("ineligible", "=", 0))
    
  # administrative fluff
  test$exposure <- list(feasible = feasible, eligible = !ineligible)
  
  return(test)
}