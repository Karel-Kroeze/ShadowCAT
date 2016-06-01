#' Best item in an adaptive test.
#' 
#' This function is a wrapper that sends the actual work to the correct subroutines.
#' 
#' @param information_summary called "objective" by Kroeze; how to summarize information; one of
#' "determinant": compute determinant(info_sofar_QxQ + info_QxQ_k) for each yet available item k
#' "posterior_determinant": compute determinant(info_sofar_QxQ_plus_prior_information + info_QxQ_k) for each yet available item k
#' "trace": compute trace((info_sofar_QxQ + info_QxQ_k) for each yet available item k
#' "posterior_trace": compute trace(info_sofar_QxQ_plus_prior_information + info_QxQ_k) for each yet available item k
#' "posterior_expected_kullback_leibler" = compute Posterior expected Kullback-Leibler Information
#' @param lp_constraints data frame with constraints in lp format: the lp_constraints from the list returned by constraints_lp_format(); NULL means no constraints
#' @param lp_characters data frame with constraint characters in lp format: the lp_chars from the list returned by constraints_lp_format(); NULL means no constraints
#' @param estimate vector containing current theta estimate
#' @param model string, one of '3PLM', 'GPCM', 'SM' or 'GRM', for the three-parameter logistic, generalized partial credit, sequential or graded response model respectively.
#' @param answers vector with person answers
#' @param prior_form String indicating the form of the prior; one of "normal" or "uniform"
#' @param prior_parameters List containing mu and Sigma of the normal prior: list(mu = ..., Sigma = ...), or 
#' the upper and lower bound of the uniform prior: list(lower_bound = ..., upper_bound = ...). Sigma should always
#' be in matrix form.
#' @param available vector with indices of yet available items
#' @param administered vector with indices of administered items
#' @param number_items number of items in the test bank
#' @param number_dimensions number of dimensions
#' @param estimator type of estimator to be used, one of "maximum_aposteriori", "maximum_likelihood", or "expected_aposteriori"
#' @param alpha matrix containing the alpha parameters
#' @param beta matrix containing the beta parameters
#' @param guessing matrix containing the quessing parameters
#' @param number_itemsteps_per_item vector containing the number of non missing cells per row of the beta matrix
#' @param stop_test rule for when to stop providing new items to patient; should be a list of the form
#' list(target = ..., max_n = ..., min_n = ..., cutoffs = ...), 
#' where max_n = test length at which testing should stop (even if target has not been reached yet in case of variance stopping rule), 
#' target = vector of maximum acceptable variances per dimension; NULL means no variance target,
#' min_n = minimum test length; NULL means no mimimum test length,
#' cutoffs = matrix containing cut off values per dimension (columns) and test iteration (rows). First row contains cut off values for when no items have been
#' administered yet, second row for when one item has been administered, etc. If estimate + 3SE < cutoff for each dimension at certain iteration, test stops; 
#' NULL means no cut off values
#' @param eap_estimation_procedure String indicating the estimation procedure for the expected aposteriori estimate, which is computed
#' in get_posterior_expected_kl_information() if it is not the requested estimator in shadowcat(). One of "riemannsum" for integration via Riemannsum or
#' "gauss_hermite_quad" for integration via Gaussian Hermite Quadrature. Only important here if information_summary is posterior_expected_kl_information.
#' @return integer item index of best item
#' @export
get_best_item <- function(information_summary, lp_constraints, lp_characters, estimate, model, answers, prior_form, prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, stop_test, eap_estimation_procedure = "riemannsum") {
  # TODO: make selection with 0 answers work as expected
  result <- function() {
    item_with_max_information <- get_item_with_max_information()
    if (is.na(item_with_max_information) || length(item_with_max_information) == 0) 
      stop("item with maximum information could not be found")
    item_with_max_information
  }
  
  get_item_with_max_information <- function() {
    item_information <- get_summarized_information(information_summary, estimate, model, answers, prior_form, prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, pad = TRUE, eap_estimation_procedure = eap_estimation_procedure)
    item_with_max_information <- switch(item_selection_type(),
                                        "maximum_information_only" = get_item_index_max_information(available, item_information, estimate, stop_test, alpha, number_answers = length(administered)),
                                        "with_constraints" = get_item_index_max_information_constrained(number_items, administered, available, answers, lp_constraints, lp_characters, item_information))
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
