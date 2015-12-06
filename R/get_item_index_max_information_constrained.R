#' Shadow Test item selection
#'
#' Based on the Shadow Testing adaptive test assembly procedure by van der Linden (XXX)
#'
#' At each step in a CAT, a Shadow Test is created consisting of an optimal complete test meeting all constraints.
#' The best item that was not already administrated is then selected from this test, and presented to the respondent.
#'
#' @param number_items number of items
#' @param administered vector with indeces of administered items
#' @param available vector with indeces of yet available items
#' @param responses vector with person responses
#' @param constraints data frame with constraints in lp format: the constraints list returned by constraints_correct_format()
#' @param lp_characters data frame with constraint characters in lp format: the lp_chars list returned by constraints_correct_format()
#' @param item_information item_information vector with information of each yet available item, with zeros for administered items (returned by objective() with pad = TRUE)
#' @return integer item index of item with maximum information within constraints
#' @importFrom lpSolve lp
get_item_index_max_information_constrained <- function(number_items, administered, available, responses, constraints, lp_characters, item_information) {
  administered_binary <- sapply(1:number_items, FUN = function(x) { if (x %in% administered) 1 else 0 } )
  
  result <- function() {
    solution <- get_lp_solution()
    
    # get items from the pool that have the max value from the solution set.
    item_with_max_information <- which(item_information == max(item_information[as.logical(solution)]))

    # since it is theoretically possible that items not in the set share this value, cull them.
    if (length(item_with_max_information) > 1)
      intersect(item_with_max_information, available)
    else
      item_with_max_information
  }
  
  get_lp_solution <- function() {
    # user created constraints are combined with constraint to select all administered items.
    lp(direction = 'max',
       objective.in = item_information,
       const.mat = as.matrix(cbind(lp_characters, administered_binary)),
       const.dir = c(constraints$op, "="),
       const.rhs = c(constraints$target, length(responses)),
       all.bin = TRUE,
       transpose.constraints = FALSE)$solution
  }

  validate <- function() {
    if (is.null(constraints))
      add_error("constraints", "matrix not given. Shadow Testing without constraints is equivalent to maximum information selection, but with more overhead.")
  }

  invalid_result <- function() {
    list(errors = errors())
  }

  validate_and_run()
}
