#' Next item in an adaptive test. Takes a person and test object, and returns the index of the next item based 
#' 
#' This function is a wrapper that sends the actual work to the correct subroutines.
#' 
#' @param information_summary called "objective" by Kroeze; how to summarize information; one of
#' "determinant": compute determinant(info_sofar_QxQ + info_QxQ_k) for each yet available item k
#' "posterior_determinant": compute determinant(info_sofar_QxQ_plus_prior + info_QxQ_k) for each yet available item k
#' "trace": compute trace((info_sofar_QxQ + info_QxQ_k) for each yet available item k
#' "posterior_trace": compute trace(info_sofar_QxQ_plus_prior + info_QxQ_k) for each yet available item k
#' "posterior_expected_kullback_leibler" = compute Posterior expected Kullback-Leibler Information
#' @param lp_constraints data frame with constraints in lp format: the lp_constraints from the list returned by constraints_lp_format(); NULL means no constraints
#' @param lp_characters data frame with constraint characters in lp format: the lp_chars from the list returned by constraints_lp_format(); NULL means no constraints
#' @param estimate vector containing current theta estimate
#' @param model string, one of '3PLM', 'GPCM', 'SM' or 'GRM', for the three-parameter logistic, generalized partial credit, sequential or graded response model respectively.
#' @param responses vector with person responses
#' @param prior prior covariance matrix for theta
#' @param available vector with indeces of yet available items
#' @param administered vector with indeces of administered items
#' @param number_items number of items in the test bank
#' @param number_dimensions number of dimensions
#' @param estimator type of estimator to be used, one of "maximum_aposteriori", "maximum_likelihood", or "expected_aposteriori"
#' @param alpha matrix containing the alpha parameters
#' @param beta matrix containing the beta parameters
#' @param guessing matrix containing the quessing parameters
#' @param number_itemsteps_per_item vector containing the number of non missing cells per row of the beta matrix
#' @param lower_bound vector with lower bounds for theta per dimension; estimated theta values smaller than the lowerbound values are truncated to the lowerbound values
#' @param upper_bound vector with upper bounds for theta per dimension; estimated theta values larger than the upperbound values are truncated to the upperbound values
#' @return integer item index
#' @export
get_best_item <- function(information_summary, lp_constraints, lp_characters, estimate, model, responses, prior, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound) {
  # TODO: make selection with 0 responses work as expected
  result <- function() {
    item_with_max_information <- get_item_with_max_information()
    
    # if there's no answer for some obfuscated reason, return a random available item rather than crashing and burning.
    # RM try to find another way of dealing with this situation
    if (is.na(item_with_max_information) || length(item_with_max_information) == 0) 
      sample(available, 1)
    else
      item_with_max_information
  }
  
  get_item_with_max_information <- function() {
    # get the values of the objective function for this test/person combo
    item_information <- get_item_information(information_summary, estimate, model, responses, prior, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, pad = TRUE)
    
    # find item with largest information; 'MI' is a simple maximum, 'Shadow' does ShadowTesting.
    item_with_max_information <- switch(item_selection_type(),
                                        "maximum_information_only" = get_item_index_max_information(available, item_information),
                                        "with_constraints" = get_item_index_max_information_constrained(number_items, administered, available, responses, lp_constraints, lp_characters, item_information))
    
    # if there's more than one, select one at random (edge case)
    if (length(item_with_max_information) > 1) 
      sample(item_with_max_information, 1)
    else
      item_with_max_information
  }
  
  item_selection_type <- function() {
    if (is.null(lp_constraints))
      "maximum_information_only"
    else
      "with_constraints"
  }
  
  result()
}
