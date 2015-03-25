#' Next item in an adaptive test. Takes a person and test object, and returns the index of the next item based 
#' 
#' This function is a wrapper that sends the actual work to the (hopefully) correct subroutines.
#' 
#' @param test
#' @param person
#' @return integer item index
#' @export
next_item <- function(test, person) {
  # get the values of the objective function for this test/person combo
  objective <- objective(test, person, TRUE)
  
  # pawn out the actual work, 'MI' is a simple maximum, 'Shadow' does ShadowTesting.
  out <- switch(test$selection,
                "MI" = MI(test, person, objective),
                "Shadow" = Shadow(test, person, objective))
  
  #cat("\n Next item: ", out)
  if (is.na(out)) {
    out <- sample(person$available, 1)
  #  cat("\n Backup item: ", out)
  }
  
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
  if (test$objective %in% c('TFI','DFI','TPFI','DPFI')) {
    # Fisher information based criteria
    info <- FI(test, person)
    so_far <- apply(info[,,person$administered, drop = FALSE], c(1,2), sum)
    
    # add prior
    if (test$objective %in% c('TPFI', 'DPFI')) so_far <- so_far + solve(person$prior)
    
    # determinant
    if (test$objective %in% c('DFI', 'DPFI')) {
      out <- apply(info[,,person$available, drop = FALSE], 3, function(x) det(so_far + x))
    }
    # trace
    if (test$objective %in% c('TFI','TPFI')){
      out <- apply(info[,,person$available, drop = FALSE], 3, function(x) sum(diag(so_far + x)))
    }
  }
  else if (test$objective == "PEKL") {
    out <- PEKL(test, person)
  } 
  else stop("unknown objective function.")
  
  if (pad) {
    full <- rep(0, test$items$K)
    full[person$available] <- out
    out <- full
  }
  
  out[is.na(out)] <- 0 # set missings to 0. I'm hoping this is underflow. TODO: investigate, (mostly occurs in 3PL weirdly enough.)
  
  return(out)
}
