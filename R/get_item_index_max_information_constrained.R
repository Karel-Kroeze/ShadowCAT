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
  result <- function() {
    administered_binary <- sapply(1:test$items$K, FUN = function(x) { if (x %in% person$administered) 1 else 0 } )
    
    # user created constraints are combined with constraint to select all administered items.
    # RM: I have removed the length = nr_items restrictions restriction; 
    # the restriction causes problems since it does not allow variables (Xi) to be equal to 0; i.e., the whole solution vector can only be one or zero (= no solution)
    solution <- lp(direction = 'max',
                   objective.in = item_information,
                   const.mat = as.matrix(cbind(test$constraints$lp_chars[, -1], administered_binary)),
                   const.dir = c(test$constraints$constraints$op[-1], "="),
                   const.rhs = c(test$constraints$constraints$target[-1], length(person$responses)),
                   all.bin = TRUE,
                   transpose.constraints = FALSE)$solution
    
    # get items from the pool that have the max value from the solution set.
    which(item_information == max(item_information[as.logical(solution)]))
    
    # since it is theoretically possible that items not in the set share this value, cull them.
    # RM administered items have a value of 0, so this won't happen
    # if (length(out) > 1) out <- intersect(out, person$available)
  }
  
  validate <- function() {
    # RM will improve this if statement later
    if (ncol(as.matrix(test$constraints$lp_chars)[, -1]) == 0)
      add_error("constraints", "matrix not given. Shadow Testing without constraints is equivalent to maximum information selection, but with more overhead.")
    if (is.null(test))
      add_error("test", "is missing")
    if (is.null(person))
      add_error("person", "is missing")
    if (is.null(item_information))
      add_error("item_information", "is missing")
  }
  
  invalid_result <- function() {
    list(errors = errors())
  }
  
  validate_and_run()
}