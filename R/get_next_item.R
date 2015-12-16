#' Next item
#' 
#' Control function to select the next item, takes account of starting conditions.
#' 
#' @param start_items items that are shown to the patient before adaptive proces starts; one of
#' list(type = 'random', n)
#' list(type = 'fixed', indices, n)
#' list(type = 'randomByDimension', nByDimension, n)
#' where n = total number of initial items, indices = vector of initial item indeces, 
#' nByDimension = scalar of number of initial items per dimension, or vector with number of initial items for each dimension
#' @param item_selection selection criterion; one of "MI" (maximum information) or "Shadow" (maximum information and take constraints into account)
#' @param information_summary called "objective" by Kroeze; how to summarize information; one of
#' "D" = determinant: compute determinant(info_sofar_QxQ + info_QxQ_k) for each yet available item k
#' "PD" = posterior determinant: compute determinant(info_sofar_QxQ_plus_prior + info_QxQ_k) for each yet available item k
#' "A" = trace: compute trace((info_sofar_QxQ + info_QxQ_k) for each yet available item k
#' "PA" = posterior trace: compute trace(info_sofar_QxQ_plus_prior + info_QxQ_k) for each yet available item k
#' "PEKL" = compute Posterior expected Kullback-Leibler Information
#' @param lp_constraints data frame with constraints in lp format: the lp_constraints from the list returned by constraints_lp_format()
#' @param lp_characters data frame with constraint characters in lp format: the lp_chars from the list returned by constraints_lp_format()
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
#' @return integer item index next item
#' @export
get_next_item <- function(start_items, item_selection, information_summary, lp_constraints, lp_characters, estimate, model, responses, prior, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound) {
  result <- function() {
    if (length(responses) < start_items$n)
      get_start_item(start_items$type)
    else
      get_best_item(item_selection, information_summary, lp_constraints, lp_characters, estimate, model, responses, prior, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound)
  }
  
  get_start_item <- function(start_type) {
    switch(start_type,
           "random" = get_start_item_random(),
           "fixed" = get_start_item_fixed(),
           "randomByDimension" = get_start_item_randomByDimension())
  }
  
  get_start_item_random <- function() {
    sample(available, 1)
  }
  
  get_start_item_fixed <- function() {
    start_items$indices[length(responses) + 1]
  }
  
  # picks nByDimension starting items per dimension (or n_i if nByDimension is a length Q vector), assumes between models, 
  # if any item has a non-zero loading on a dimension, it is considered to be part of that dimension. 
  # they CAN overlap, which may cause unwanted side effects, and in within models the result is identical to 'normal' random starting.
  get_start_item_randomByDimension <- function() {
    n_by_dimension_vector <- get_n_by_dimension_vector()
    design_matrix_item_loadings <- alpha > 0
    
    # if enough items from first dimension are drawn, items from next dimension are drawn, etc.
    dimension <- find_dimension_to_draw_from(n_by_dimension_vector)
     
    # get items for this dimension, cull non-applicable, and select at random.
    items_with_positive_loading_on_dimension <- (1:number_items)[design_matrix_item_loadings[,dimension]]
    positive_loading_and_available <- intersect(items_with_positive_loading_on_dimension, available)
    sample(positive_loading_and_available, 1)
  }
  
  get_n_by_dimension_vector <- function() {
    if (length(start_items$nByDimension) == 1) 
      rep(start_items$nByDimension, number_dimensions)
    else
      start_items$nByDimension
  }
  
  find_dimension_to_draw_from <- function(n_by_dimension_vector) {
    sum(length(responses) >= cumsum(n_by_dimension_vector)) + 1
  }
  
  validate <- function() {
    if (start_items$type == "randomByDimension" && length(start_items$nByDimension) > 1 && start_items$n != sum(start_items$nByDimension))
      return(add_error("start_items", "contains inconsistent information. Total length of start phase and sum of length per dimension do not match (n != sum(nByDimension)"))
    if (start_items$type == "randomByDimension" && length(start_items$nByDimension) == 1 && start_items$n != sum(rep(start_items$nByDimension, number_dimensions)))
      return(add_error("start_items", "contains inconsistent information. Total length of start phase and sum of length per dimension do not match"))
  }
  
  invalid_result <- function() {
    list(errors = errors())
  }
  
  validate_and_run()
}
