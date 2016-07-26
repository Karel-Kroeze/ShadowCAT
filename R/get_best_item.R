#' Select best item
#' 
#' Get the index of the item with maximum information from the available items, possibly taking constraints into account.
#' 
#' @param lp_constraints Data frame with constraints in \code{\link{lp}} format: the \code{lp_constraints} from the list returned by \code{\link{constraints_lp_format}}. \code{NULL} means no constraints.
#' @param lp_characters Data frame with constraint characteristics in \code{\link{lp}} format: the \code{lp_chars} from the list returned by \code{\link{constraints_lp_format}}. \code{NULL} means no constraints.
#' @param estimate Vector containing current theta estimate, with covariance matrix as an attribute.
#' @param answers Vector with answers to administered items.
#' @param available Vector with indices of yet available items.
#' @param administered Vector with indices of administered items.
#' @param number_items Number of items in the test bank.
#' @param number_dimensions Number of dimensions of theta.
#' @param number_itemsteps_per_item Vector containing the number of non missing cells per row of the beta matrix.
#' @inheritParams shadowcat
#' @return Item index of best item.
get_best_item <- function(information_summary, lp_constraints, lp_characters, estimate, model, answers, prior_form, prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, stop_test, eap_estimation_procedure = "riemannsum") {
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
