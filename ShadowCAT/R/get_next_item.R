#' Index next item
#' 
#' Get the index of the next item to adminster, taking starting (burn in) conditions into account.
#' 
#' @param lp_constraints Data frame with constraints in \code{\link{lp}} format: the \code{lp_constraints} from the list returned by \code{\link{constraints_lp_format}}. \code{NULL} means no constraints.
#' @param lp_characters Data frame with constraint characteristics in \code{\link{lp}} format: the \code{lp_chars} from the list returned by \code{\link{constraints_lp_format}}. \code{NULL} means no constraints.
#' @param estimate Vector containing current theta estimate, with covariance matrix as an attribute.
#' @param answers Vector with answers to administered items.
#' @param available vector with indices of yet available items.
#' @param administered Vector with indices of administered items.
#' @param number_items Number of items in test bank.
#' @param number_dimensions Number of dimensions in theta.
#' @param number_itemsteps_per_item Vector containing the number of non missing cells per row of the beta matrix
#' @inheritParams shadowcat
#' @return integer item index next item
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
    indices <- match(start_items$item_keys, rownames(alpha))
    indices[length(answers) + 1]
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
