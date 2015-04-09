#' Next item in an adaptive test. Takes a person and test object, and returns the index of the next item based 
#' 
#' This function is a wrapper that sends the actual work to the correct subroutines.
#' 
#' @param test
#' @param person
#' @return integer item index
#' @export
best_item <- function(person, test) {
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


#' Obtain a vector of values of the objective function. 
#' 
#' @param test
#' @param person
#' @param pad Should the return vector be padded with zeros for items that have already been administered?
#' @return vector
#' @export
objective <- function(test, person, pad = TRUE) {
  if (test$objective %in% c('A','D','PA','PD')) {
    # Fisher information based criteria
    info <- FI(test, person)
    
    so_far <- apply(info[,,person$administered, drop = FALSE], c(1,2), sum)
    
    # add prior
    if (test$objective %in% c('PA', 'PD')) so_far <- so_far + solve(person$prior)
    
    # determinant
    if (test$objective %in% c('PD', 'D')) {
      out <- apply(info[,,person$available, drop = FALSE], 3, function(x) det(so_far + x))
    }
    # trace
    if (test$objective %in% c('PA','A')){
      out <- apply(info[,,person$available, drop = FALSE], 3, function(x) sum(diag(so_far + x)))
    }
  }
  else if (test$objective == "PEKL") {
    out <- PEKL(test, person)
  } 
  else stop("unknown objective function.")
  
  # If all objective values are 0, something went horribly wrong.
  # This is made worse by lpSolve -> it will give back a full vector, not respecting constraints.
  # TODO: check if this is ok.
  if (all(out == 0)) {
    out <- rep(1, length(out)) # so replace them by 1's.
    warning("Objective is (computationally) zero for all items. Are you sure appropriate starting items are selected?")
  } 
  
  # pad to full length K objective vector.
  if (pad) {
    full <- rep(0, test$items$K)
    full[person$available] <- out
    out <- full
  }
  
  # set missings to 0. I'm hoping this is underflow.
  # TODO: investigate, (mostly occurs in 3PLM weirdly enough.)
  out[is.na(out)] <- 0  
  
  return(out)
}
