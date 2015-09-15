#' Shadow Test item selection
#' 
#' Based on the Shadow Testing adaptive test assembly procedure by van der Linden (XXX)
#' 
#' At each step in a CAT, a Shadow Test is created consisting of an optimal complete test meeting all constraints.
#' The best item that was not already administrated is then selected from this test, and presented to the respondent.
#' 
#' @param test test object
#' @param person person object
#' @param item_information item_information vector with information of each yet available item, with zeros for administered items (returned by objective() with pad = TRUE)
#' @return integer item index of item with maximum information within constraints
#' @importFrom lpSolve lp
Shadow <- function(test, person, item_information) {
  # produce a small warning about Shadow Testing without constraints.
  if (is.null(test$constraints)) stop("No constraints matrix given. 
                                      Shadow Testing without constraints is equivalent to maximum information selection, but with more overhead.")
  
  # user created constraints are combined with constraint to select all administered items.
  # RM: I have removed the administered_binair = length(responses) and length = nr_items restrictions; 
  # I think the first one is not useful (it is already accounted for in the objective function) and may cause problems, the latter one always causes problems since it
  # does not allow variables (Xi) to be equal to 0; i.e., the whole solution vector can only be one or zero (= no solution)
  solution <- lp(direction = 'max',
                 objective.in = item_information,
                 const.mat = as.matrix(test$constraints$lp_chars)[, -1],
                 const.dir = test$constraints$constraints$op[-1],
                 const.rhs = test$constraints$constraints$target[-1],
                 all.bin = TRUE,
                 transpose.constraints = FALSE)$solution
  
  # get items from the pool that have the max value from the solution set.
  out <- which(item_information == max(item_information[as.logical(solution)]))
  
  # since it is theoretically possible that items not in the set share this value, cull them.
  if (length(out) > 1) out <- intersect(out, person$available)
  
  # return index
  return(out)
}