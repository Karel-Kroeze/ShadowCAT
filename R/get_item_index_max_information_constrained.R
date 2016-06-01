#' Shadow Test item selection
#'
#' Based on the Shadow Testing adaptive test assembly procedure by van der Linden (XXX)
#'
#' At each step in a CAT, a Shadow Test is created consisting of an optimal complete test meeting all constraints.
#' The best item that was not already administrated is then selected from this test, and presented to the respondent.
#'
#' @param number_items number of items in the test bank
#' @param administered vector with indices of administered items
#' @param available vector with indices of yet available items
#' @param answers vector with person answers
#' @param lp_constraints data frame with constraints in lp format: the lp_constraints from the list returned by constraints_lp_format()
#' @param lp_characters data frame with characteristics in lp format: the lp_chars from the list returned by constraints_lp_format()
#' @param item_information vector with summarized information of each yet available item, with zeros for administered items (returned by get_summarized_information() with pad = TRUE)
#' @return integer item index of item with maximum information within constraints
#' @importFrom lpSolve lp
get_item_index_max_information_constrained <- function(number_items, administered, available, answers, lp_constraints, lp_characters, item_information) {
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
       const.dir = c(lp_constraints$op, "="),
       const.rhs = c(lp_constraints$target, length(answers)),
       all.bin = TRUE,
       transpose.constraints = FALSE)$solution
  }

  validate <- function() {
    if (is.null(lp_constraints))
      add_error("lp_constraints", "matrix not given. Shadow Testing without constraints is equivalent to maximum information selection, but with more overhead.")
  }

  invalid_result <- function() {
    list(errors = errors())
  }

  validate_and_run()
}
