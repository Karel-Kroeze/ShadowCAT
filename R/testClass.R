#' Initiate ShadowCAT test object
#' 
#' Creates an object to be used with ShadowCAT
#' 
#' @param items
#' @param start
#' @param stop
#' @param estimator
#' @param objective
#' @param selection
#' @param constraints
#' @param ...
#' @return ShadowCAT.test
#' @export
initTest <- function(items, start = NULL, stop = NULL,
                     estimator = 'MAP', objective = 'PFI', selection = 'MI', constraints = NULL,
                     lowerBound = rep(-6, items$Q), upperBound = rep(6, items$Q), ...){
  # TODO: input validation
  
  # default start is randomly selecting 5 items.
  if (is.null(start)) start <- list(type = 'random', n = 5) 
  
  # default stop rule is length of 30.
  if (is.null(stop)) stop <- list(type = 'length', n = 30)
    
  # attach everything
  out <- list(items = items,
              start = start,
              stop = stop,
              lowerBound = lowerBound,
              upperBound = upperBound,
              estimator = estimator,
              objective = objective,
              selection = selection,
              constraints = constraints,
              internal = list(...))
  
  attr(out, 'class') <- c("ShadowCAT.test")
  
  return(invisible(out))
}