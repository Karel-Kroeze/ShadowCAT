#' Next item
#' 
#' Control function to select the next item, takes account of starting conditions.
#' 
#' @param start_items items that are shown to the patient before adaptive proces starts; one of
#' list(type = 'random', n)
#' list(type = 'fixed', indices, n)
#' list(type = 'random_by_dimension', n_by_dimension, n)
#' where n = total number of initial items, indices = vector of initial item indices, 
#' n_by_dimension = scalar of number of initial items per dimension, or vector with number of initial items for each dimension
#' @param information_summary called "objective" by Kroeze; how to summarize information; one of
#' "determinant": compute determinant(info_sofar_QxQ + info_QxQ_k) for each yet available item k
#' "posterior_determinant": compute determinant(info_sofar_QxQ_plus_prior_information + info_QxQ_k) for each yet available item k
#' "trace": compute trace((info_sofar_QxQ + info_QxQ_k) for each yet available item k
#' "posterior_trace": compute trace(info_sofar_QxQ_plus_prior_information + info_QxQ_k) for each yet available item k
#' "posterior_expected_kullback_leibler" = compute Posterior expected Kullback-Leibler Information
#' @param lp_constraints data frame with constraints in lp format: the lp_constraints from the list returned by constraints_lp_format(); NULL means no constraints
#' @param lp_characters data frame with constraint characters in lp format: the lp_chars from the list returned by constraints_lp_format(); NULL means no constraints
#' @param estimate vector with current theta estimate
#' @param model string, one of '3PLM', 'GPCM', 'SM' or 'GRM', for the three-parameter logistic, generalized partial credit, sequential or graded response model respectively.
#' @param answers vector with person answers
#' @param prior_form String indicating the form of the prior; one of "normal" or "uniform"
#' @param prior_parameters List containing mu and Sigma of the normal prior: list(mu = ..., Sigma = ...), or 
#' the upper and lower bound of the uniform prior: list(lower_bound = ..., upper_bound = ...). Sigma should always
#' be in matrix form.
#' @param available vector with indices of yet available items
#' @param administered vector with indices of administered items
#' @param number_items number of items in test bank
#' @param number_dimensions number of dimensions
#' @param estimator type of estimator to be used, one of "maximum_aposteriori", "expected_aposteriori", or "maximum_likelihood"
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
#' @return integer item index next item
#' @export
get_next_item <- function(start_items, information_summary, lp_constraints, lp_characters, estimate, model, answers, prior_form, prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, stop_test, eap_estimation_procedure = "riemannsum") {
  result <- function() {
    if (length(answers) < start_items$n)
      get_start_item(start_items$type)
    else
      get_best_item(information_summary, lp_constraints, lp_characters, estimate, model, answers, prior_form, prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, stop_test, eap_estimation_procedure = eap_estimation_procedure)
  }
  
  get_start_item <- function(start_type) {
    switch(start_type,
           "random" = get_start_item_random(),
           "fixed" = get_start_item_fixed(),
           "random_by_dimension" = get_start_item_random_by_dimension())
  }
  
  get_start_item_random <- function() {
    sample(available, 1)
  }
  
  get_start_item_fixed <- function() {
    start_items$indices[length(answers) + 1]
  }
  
  # picks n_by_dimension starting items per dimension (or n_i if n_by_dimension is a length Q vector), assumes that items load on a single dimension, 
  # if any item has a non-zero loading on a dimension, it is considered to be part of that dimension.
  get_start_item_random_by_dimension <- function() {
    n_by_dimension_vector <- get_n_by_dimension_vector()
    design_matrix_item_loadings <- alpha > 0
    dimension <- find_dimension_to_draw_from(n_by_dimension_vector)
     
    # get items for this dimension, cull administered items, and select at random.
    items_with_positive_loading_on_dimension <- (1:number_items)[design_matrix_item_loadings[,dimension]]
    positive_loading_and_available <- intersect(items_with_positive_loading_on_dimension, available)
    sample(positive_loading_and_available, 1)
  }
  
  get_n_by_dimension_vector <- function() {
    if (length(start_items$n_by_dimension) == 1) 
      rep(start_items$n_by_dimension, number_dimensions)
    else
      start_items$n_by_dimension
  }
  
  # if enough items from first dimension are drawn, items from next dimension are drawn, etc.
  find_dimension_to_draw_from <- function(n_by_dimension_vector) {
    sum(length(answers) >= cumsum(n_by_dimension_vector)) + 1
  }
  
  validate <- function() {
    if (start_items$type == "random_by_dimension" && length(start_items$n_by_dimension) > 1 && start_items$n != sum(start_items$n_by_dimension))
      return(add_error("start_items", "contains inconsistent information. Total length of start phase and sum of length per dimension do not match (n != sum(n_by_dimension)"))
    if (start_items$type == "random_by_dimension" && length(start_items$n_by_dimension) == 1 && start_items$n != sum(rep(start_items$n_by_dimension, number_dimensions)))
      return(add_error("start_items", "contains inconsistent information. Total length of start phase and sum of length per dimension do not match"))
  }
  
  invalid_result <- function() {
    list(errors = errors())
  }
  
  validate_and_run()
}
