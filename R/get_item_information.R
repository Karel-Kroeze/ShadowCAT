#' Obtain a vector with information for each available item (values of the objective function). 
#'
#' @param information_summary called "objective" by Kroeze; how to summarize information; one of
#' "determinant": compute determinant(info_sofar_QxQ + info_QxQ_k) for each yet available item k
#' "posterior_determinant": compute determinant(info_sofar_QxQ_plus_prior + info_QxQ_k) for each yet available item k
#' "trace": compute trace((info_sofar_QxQ + info_QxQ_k) for each yet available item k
#' "posterior_trace": compute trace(info_sofar_QxQ_plus_prior + info_QxQ_k) for each yet available item k
#' "posterior_expected_kullback_leibler" = compute Posterior expected Kullback-Leibler Information
#' @param estimate vector with current theta estimate
#' @param model string, one of '3PLM', 'GPCM', 'SM' or 'GRM', for the three-parameter logistic, generalized partial credit, sequential or graded response model respectively.
#' @param responses vector with person responses
#' @param prior prior covariance matrix for theta
#' @param available vector with indeces of yet available items
#' @param administered vector with indeces of administered items
#' @param number_items number of items in test bank
#' @param number_dimensions number of dimensions
#' @param estimator type of estimator to be used, one of "maximum_aposteriori", "maximum_likelihood", or "expected_aposteriori"
#' @param alpha matrix containing the alpha parameters
#' @param beta matrix containing the beta parameters
#' @param guessing matrix containing the quessing parameters
#' @param number_itemsteps_per_item vector containing the number of non missing cells per row of the beta matrix
#' @param lower_bound vector with lower bounds for theta per dimension; estimated theta values smaller than the lowerbound values are truncated to the lowerbound values
#' @param upper_bound vector with upper bounds for theta per dimension; estimated theta values larger than the upperbound values are truncated to the upperbound values
#' @param pad Should the return vector be padded with zeros for items that have already been administered?
#' @return vector with information for each available item
#' @export
get_item_information <- function(information_summary, estimate, model, responses, prior, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, pad = TRUE) {
  fisher_information <- get_fisher_information(estimate, model, number_dimensions, alpha, beta, guessing, number_itemsteps_per_item)
  result <- function() {
    item_information <- get_item_information_switch()
    if (pad) 
      pad_zeros(item_information)
    else
      item_information  
  }
  
  pad_zeros <- function(item_information) {
    item_information_padded <- rep(0, number_items)
    item_information_padded[available] <- item_information
    item_information_padded
  }
  
  get_item_information_switch <- function() {
    item_information <- switch(information_summary,
                               "trace" = item_information_trace(),
                               "posterior_trace" = item_information_post_trace(),
                               "determinant" = item_information_determinant(),
                               "posterior_determinant" = item_information_post_determinant(),
                               "posterior_expected_kullback_leibler" = item_information_pekl())
    if (any(is.na(item_information)))
      stop("Information is NA from some items ")
    if (all(item_information == 0))
      stop("Information is zero for all items.\n")
    item_information      
  }
  
  item_information_trace <- function() {
    information_administered <- apply(fisher_information[,,administered, drop = FALSE], c(1, 2), sum)
    apply(fisher_information[,,available, drop = FALSE], 3, function(x) sum(diag(information_administered + x)))
  }
  
  item_information_post_trace <- function() {
    information_administered <- apply(fisher_information[,,administered, drop = FALSE], c(1, 2), sum) + solve(prior)
    apply(fisher_information[,,available, drop = FALSE], 3, function(x) sum(diag(information_administered + x)))
  }
  
  item_information_determinant <- function() {
    information_administered <- apply(fisher_information[,,administered, drop = FALSE], c(1, 2), sum)
    apply(fisher_information[,,available, drop = FALSE], 3, function(x) det(information_administered + x))
  }

  item_information_post_determinant <- function() {
    information_administered <- apply(fisher_information[,,administered, drop = FALSE], c(1, 2), sum) + solve(prior)
    apply(fisher_information[,,available, drop = FALSE], 3, function(x) det(information_administered + x))
  }

  item_information_pekl <- function() {
    get_posterior_expected_kl_information(estimate, model, responses, administered, available, number_dimensions, estimator, alpha, beta, guessing, prior, number_itemsteps_per_item, lower_bound, upper_bound)
  }
  
  validate <- function() {
    if (information_summary %not_in% c("trace", "posterior_trace", "determinant", "posterior_determinant", "posterior_expected_kullback_leibler"))
      add_error("information_summary", "of unknown type")
  }
  
  invalid_result <- function() {
    list(errors = errors())
  }
  
  validate_and_run()
}
