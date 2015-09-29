#' Obtain a vector with information for each available item (values of the objective function). 
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
    cat("\nObjective is (computationally) zero for all items.")
  } 
  
  # pad to full length K objective vector.
  if (pad) {
    full <- rep(0, test$items$K)
    full[person$available] <- out
    out <- full
  }
  
  # set missings to 0. I'm hoping this is underflow.
  # TODO: investigate / remove, (mostly occurs in 3PLM weirdly enough.)
  if (any(is.na(out))) cat("\nMissing values in objective function.")
  out[is.na(out)] <- 0  
  
  return(out)
}
