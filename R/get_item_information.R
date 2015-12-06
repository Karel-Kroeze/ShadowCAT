#' Obtain a vector with information for each available item (values of the objective function). 
#'
#' @param information_summary called "objective" by Kroeze; how to summarize information; one of
#' "D" = determinant: compute determinant(info_sofar_QxQ + info_QxQ_k) for each yet available item k
#' "PD" = posterior determinant: compute determinant(info_sofar_QxQ_plus_prior + info_QxQ_k) for each yet available item k
#' "A" = trace: compute trace((info_sofar_QxQ + info_QxQ_k) for each yet available item k
#' "PA" = posterior trace: compute trace(info_sofar_QxQ_plus_prior + info_QxQ_k) for each yet available item k
#' "PEKL" = compute Posterior expected Kullback-Leibler Information
#' @param estimate current theta estimate
#' @param model string, one of '3PLM', 'GPCM', 'SM' or 'GRM', for the three-parameter logistic, generalized partial credit, sequential or graded response model respectively.
#' @param responses vector with person responses
#' @param prior prior covariance matrix for theta; only required if estimator is "MAP" or "EAP" and output is "likelihoods" or "both"
#' @param available vector with indeces of yet available items
#' @param administered vector with indeces of administered items
#' @param number_items number of items
#' @param number_dimensions number of dimensions
#' @param estimator type of estimator to be used, one of "MAP" (Maximum a posteriori estimation), "EAP" (Expected A Posteriori Estimation), or "ML" (maximum likelihood)
#' @param alpha matrix containing the alpha parameters
#' @param beta matrix containing the beta parameters
#' @param guessing matrix containing the quessing
#' @param number_itemsteps_per_item number_itemsteps_per_item vector containing the number of non missing cells per row of the beta matrix
#' @param lower_bound vector with lower bounds for theta per dimension; estimated theta values smaller than the lowerbound values are truncated to the lowerbound values
#' @param upper_bound vector with upper bounds for theta per dimension; estimated theta values larger than the upperbound values are truncated to the upperbound values
#' @param pad Should the return vector be padded with zeros for items that have already been administered?
#' @return vector
#' @export
get_item_information <- function(information_summary, estimate, model, responses, prior, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, pad = TRUE) {
  fisher_information <- get_fisher_information(estimate, model, number_dimensions, estimator, alpha, beta, guessing,  number_itemsteps_per_item)
  result <- function() {
    item_information <- get_item_information_switch()
    item_information_imputed_missings <- impute_zero_for_na(item_information)   
    if (pad) 
      pad_zeros(item_information_imputed_missings)
    else
      item_information_imputed_missings  
  }
  
  pad_zeros <- function(item_information) {
    item_information_padded <- rep(0, number_items)
    item_information_padded[available] <- item_information
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
    item_information <- switch(information_summary,
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
    information_administered <- apply(fisher_information[,,administered, drop = FALSE], c(1, 2), sum)
    apply(fisher_information[,,available, drop = FALSE], 3, function(x) sum(diag(information_administered + x)))
  }
  
  # PA
  item_information_post_trace <- function() {
    information_administered <- apply(fisher_information[,,administered, drop = FALSE], c(1, 2), sum) + solve(prior)
    apply(fisher_information[,,available, drop = FALSE], 3, function(x) sum(diag(information_administered + x)))
  }
  
  # D
  item_information_determinant <- function() {
    information_administered <- apply(fisher_information[,,administered, drop = FALSE], c(1, 2), sum)
    apply(fisher_information[,,available, drop = FALSE], 3, function(x) det(information_administered + x))
  }
  
  # PD
  item_information_post_determinant <- function() {
    information_administered <- apply(fisher_information[,,administered, drop = FALSE], c(1, 2), sum) + solve(prior)
    apply(fisher_information[,,available, drop = FALSE], 3, function(x) det(information_administered + x))
  }
  
  # PEKL
  item_information_pekl <- function() {
    get_posterior_expected_kl_information(estimate, model, responses, administered, available, number_dimensions, estimator, alpha, beta, guessing, prior, number_itemsteps_per_item, lower_bound, upper_bound)
  }
  
  validate <- function() {
    if (information_summary %not_in% c("A", "PA", "D", "PD", "PEKL"))
      add_error("objective", "of unknown type")
  }
  
  invalid_result <- function() {
    list(errors = errors())
  }
  
  validate_and_run()
}
