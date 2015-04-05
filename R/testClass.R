#' Initiate ShadowCAT test object
#' 
#' Creates an object to be used with ShadowCAT
#' 
#' @param items
#' @param start
#' @param stop
#' @param max_n
#' @param estimator
#' @param objective
#' @param selection
#' @param constraints
#' @param lowerBound
#' @param upperBound
#' @param ...
#' @return ShadowCAT.test
#' @export
initTest <- function(items, 
                     start = list(type = 'random', n = 5), 
                     stop = list(type = 'length', n = 30),
                     max_n = 50, # utter maximum
                     estimator = 'MAP',
                     objective = 'PD',
                     selection = 'MI',
                     constraints = NULL,
                     lowerBound = rep(-6, items$Q),
                     upperBound = rep(6, items$Q),
                     ...)
  {

  # attach everything
  out <- list(items = items,
              start = start,
              stop = stop,
              max_n = max_n,
              lowerBound = lowerBound,
              upperBound = upperBound,
              estimator = estimator,
              objective = objective,
              selection = selection,
              constraints = constraints,
              internal = list(...))
  
  attr(out, 'class') <- c("ShadowCAT.test")  
  
  # set up default constraints
  out$constraints <- createConstraints(out, constraints)
  
  return(invisible(out))
}