#' Shadow Test item selection
#' 
#' Based on the Shadow Testing adaptive test assembly procedure by van der Linden (XXX)
#' 
#' At each step in a CAT, a Shadow Test is created consisting of an optimal complete test meeting all constraints.
#' The best item that was not already administrated is then selected from this test, and presented to the respondent.
#' 
#' @param test
#' @param person
#' @param objective
#' @return integer item index
#' @importFrom lpSolve lp
Shadow <- function(test, person, objective) {
  # produce a small warning about Shadow Testing without constraints.
  if (is.null(test$constraints)) stop("No constraints matrix given. 
                                      Shadow Testing without constraints is equivalent to maximum information selection, but with more overhead.")
  
  # proceed to the meat of the thing.
  # user created constraints are combined with constraint to select all administered items.
  # CRITICAL ERROR: TODO: update administered constraint
  possible <- solution <- with(test$constraints, lp(direction = 'max',
                            const.mat = cbind(lp_chars, person$administered),
                            const.dir = c(constraints$operator, '=='),
                            const.rhs = c(constraints$target, length(person$responses)),
                            objective.in = objective,
                            all.bin = TRUE,
                            transpose.constraints = FALSE))$solution
  
  # get the item with the highest value of the objective function
  possible[person$administered] <- FALSE
  
  out <- which(objective == max(objective[possible]))
  
  # return index
  return(out)
}
