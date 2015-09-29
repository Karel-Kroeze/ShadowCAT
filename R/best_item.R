#' Next item in an adaptive test. Takes a person and test object, and returns the index of the next item based 
#' 
#' This function is a wrapper that sends the actual work to the correct subroutines.
#' 
#' @param test
#' @param person
#' @return integer item index
#' @export
best_item <- function(person, test) {
  # TODO: make selection with 0 responses work as expected
  # get the values of the objective function for this test/person combo
  objective <- objective(test, person, TRUE)
  
  # pawn out the actual work, 'MI' is a simple maximum, 'Shadow' does ShadowTesting.
  out <- switch(test$selection,
                "MI" = MI(test, person, objective),
                "Shadow" = Shadow(test, person, objective))
  
  # if there's more than one, select one at random (edge case)
  if (length(out) > 1) out <- sample(out, 1)
  
  # if there's no answer for some obfuscated reason, return a random available item rather than crashing and burning.
  if (length(out) == 0 || is.na(out)) out <- sample(person$available, 1)
  
  return(out)
}
