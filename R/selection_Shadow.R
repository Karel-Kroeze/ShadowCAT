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
  # We need to include all previously administered items in the test.
  administered_binair <- rep(0L, test$items$K)
  administered_binair[person$administered] <- 1L
    
  # user created constraints are combined with constraint to select all administered items.
  solution <- lp(direction = 'max',
                            const.mat = as.matrix(cbind(test$constraints$lp_chars, administered_binair)),
                            const.dir = c(test$constraints$constraints$op, '='),
                            const.rhs = c(test$constraints$constraints$target, length(person$responses)),
                            objective.in = objective,
                            all.bin = TRUE,
                            transpose.constraints = FALSE)$solution
  
  # get items from the pool that have the max value from the solution set.
  out <- which(objective == max(objective[as.logical(solution)]))
  
  # since it is theoretically possible that items not in the set share this value, cull them.
  if (length(out) > 1) out <- intersect(out, person$available)
  
  # return index
  return(out)
}

# ### TEST
# items <- createTestBank("GPCM")
# test <- initTest(items, stop = list(type = 'length', n = 20), estimator = 'EAP', objective = 'PEKL', selection = 'Shadow')
# test <- createConstraints(test)
# person <- initPerson(items)
# 
# Shadow(test, person, abs(rnorm(test$items$K)))

