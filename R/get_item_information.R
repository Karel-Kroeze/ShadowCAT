#' Obtain a vector with information for each available item (values of the objective function). 
#' 
#' @param test
#' @param person
#' @param pad Should the return vector be padded with zeros for items that have already been administered?
#' @return vector
#' @export
objective <- function(test, person, pad = TRUE) {
  result <- function() {
    item_information <- get_item_information_switch()
    item_information_imputed_missings <- impute_zero_for_na(item_information)   
    if (pad) 
      pad_zeros(item_information_imputed_missings)
    else
      item_information_imputed_missings  
  }
  
  pad_zeros <- function(item_information) {
    item_information_padded <- rep(0, test$items$K)
    item_information_padded[person$available] <- item_information
    item_information_padded
  }
  
  # set missings to 0. I'm hoping this is underflow.
  # TODO: investigate / remove, (mostly occurs in 3PLM weirdly enough.)
  impute_zero_for_na <- function(item_information) {
    if (any(is.na(item_information))) {
      cat("\nMissing values in objective function.\n")
      item_information[is.na(item_information)] <- 0
    }
    item_information
  }
  
  get_item_information_switch <- function() {
    item_information <- switch(test$objective,
                               "A" = item_information_trace(),
                               "PA" = item_information_post_trace(),
                               "D" = item_information_determinant(),
                               "PD" = item_information_post_determinant(),
                               "PEKL" = item_information_pekl())
    
    # If all objective values are 0, something went horribly wrong.
    # This is made worse by lpSolve -> it will give back a full vector, not respecting constraints.
    # TODO: check if this is ok.
    if (all(item_information == 0)) {
      cat("\nObjective is (computationally) zero for all items.\n")
      rep(1, length(item_information))   
    }
    else {
      item_information
    }        
  }
  
  # A
  item_information_trace <- function() {
    fisher_information <- FI(test, person)
    information_administered <- apply(fisher_information[,,person$administered, drop = FALSE], c(1, 2), sum)
    apply(fisher_information[,,person$available, drop = FALSE], 3, function(x) sum(diag(information_administered + x)))
  }
  
  # PA
  item_information_post_trace <- function() {
    fisher_information <- FI(test, person)
    information_administered <- apply(fisher_information[,,person$administered, drop = FALSE], c(1, 2), sum) + solve(person$prior)
    apply(fisher_information[,,person$available, drop = FALSE], 3, function(x) sum(diag(information_administered + x)))
  }
  
  # D
  item_information_determinant <- function() {
    fisher_information <- FI(test, person)
    information_administered <- apply(fisher_information[,,person$administered, drop = FALSE], c(1, 2), sum)
    apply(fisher_information[,,person$available, drop = FALSE], 3, function(x) det(information_administered + x))
  }
  
  # PD
  item_information_post_determinant <- function() {
    fisher_information <- FI(test, person)
    information_administered <- apply(fisher_information[,,person$administered, drop = FALSE], c(1, 2), sum) + solve(person$prior)
    apply(fisher_information[,,person$available, drop = FALSE], 3, function(x) det(information_administered + x))
  }
  
  # PEKL
  item_information_pekl <- function() {
    PEKL(test, person)
  }
  
  validate <- function() {
    if (test$objective %not_in% c("A", "PA", "D", "PD", "PEKL"))
      add_error("objective", "of unknown type")
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
