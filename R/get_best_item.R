#' Next item in an adaptive test. Takes a person and test object, and returns the index of the next item based 
#' 
#' This function is a wrapper that sends the actual work to the correct subroutines.
#' 
#' @param test test object
#' @param person person object
#' @return integer item index
#' @export
best_item <- function(person, test) {
  # TODO: make selection with 0 responses work as expected
  result <- function() {
    item_with_max_information <- get_item_with_max_information()
    
    # if there's no answer for some obfuscated reason, return a random available item rather than crashing and burning.
    # RM try to find another way of dealing with this situation
    if (is.na(item_with_max_information) || length(item_with_max_information) == 0) 
      sample(person$available, 1)
    else
      item_with_max_information
  }
  
  get_item_with_max_information <- function() {
    # get the values of the objective function for this test/person combo
    item_information <- get_item_information(test$objective, person$estimate, test$items$model, person$responses, person$prior, person$available, person$administered, test$items$K, test$items$Q, test$estimator, test$items$pars$alpha, test$items$pars$beta, test$items$pars$guessing, test$items$pars$m, test$lowerBound, test$upperBound, pad = TRUE)
    
    # find item with largest information; 'MI' is a simple maximum, 'Shadow' does ShadowTesting.
    item_with_max_information <- switch(test$selection,
                                        "MI" = get_item_index_max_information(person$available, item_information),
                                        "Shadow" = get_item_index_max_information_constrained(test$items$K, person$administered, person$available, person$responses, test$constraints$constraints, test$constraints$lp_chars, item_information))
    
    # if there's more than one, select one at random (edge case)
    if (length(item_with_max_information) > 1) 
      sample(item_with_max_information, 1)
    else
      item_with_max_information
  }
  
  validate <- function() {
    if (is.null(test))
      add_error("test", "is missing")
    if (is.null(person))
      add_error("person", "is missing")
  }
  
  invalid_result <- function() {
    list(errors = errors())
  }
  
  validate_and_run()
}
