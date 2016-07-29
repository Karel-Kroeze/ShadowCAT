#' Summarize Fisher Information
#' 
#' Obtain a vector with information, summarized into one value, for each available item.
#' 
#' @references
#' \itemize{
#' \item Segall, D. O. (1996). Multidimensional adaptive testing. Psychometrika, 61(2), 331 - 354. doi:10.1007/BF02294343.
#' \item Segall, D. O. (2000). Principles of multidimensional adaptive testing. In W. J. van der 
#'  Linden & en C. A. W. Glas (Eds.), Computerized adaptive testing: Theory and practice (pp. 53 - 74). Dordrecht: Kluwer Academic Publishers.
#' \item Van der Linden, W. J. (1999). Multidimensional Adaptive Testing with a Minimum Error-
#' Variance Criterion. Journal of Educational and Behavioral Statistics, 24(4), 398 - 412. doi:10.3102/10769986024004398.
#' }
#' 
#' @param estimate Vector with current theta estimate.
#' @param answers Vector with answers to administered items.
#' @param available Vector with indices of yet available items.
#' @param administered Vector with indices of administered items.
#' @param number_items Number of items in test bank.
#' @param number_dimensions Number of dimensions of theta.
#' @param number_itemsteps_per_item Vector containing the number of non missing cells per row of the beta matrix.
#' @param pad If \code{TRUE}, the return vector is padded with zeros for items that have already been administered.
#' @inheritParams shadowcat
#' @return Vector with summarized information for each available item.
get_summarized_information <- function(information_summary, estimate, model, answers, prior_form, prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, pad = TRUE, eap_estimation_procedure = "riemannsum") {
  fisher_information <- get_fisher_information(estimate, model, number_dimensions, alpha, beta, guessing, number_itemsteps_per_item)
  result <- function() {
    item_information <- get_item_information_switch()
    if (pad) 
      pad_zeros(item_information)
    else
      item_information  
  }
  
  redefine_information_summary <- function() {
    if (is.null(prior_form) && information_summary %in% c("posterior_determinant", "posterior_trace", "posterior_expected_kl_information"))
      stop("Estimator is maximum likelihood, information summary should be trace or determinant")
    if (prior_form == "uniform" && information_summary == "posterior_determinant")
      "determinant"
    else if (prior_form == "uniform" && information_summary == "posterior_trace")
      "trace"
    else
      information_summary   
  }
  
  pad_zeros <- function(item_information) {
    item_information_padded <- rep(0, number_items)
    item_information_padded[available] <- item_information
    item_information_padded
  }
  
  get_item_information_switch <- function() {
    item_information <-  switch(redefine_information_summary(),
                                "trace" = item_information_trace(),
                                "posterior_trace" = item_plus_prior_information_trace(),
                                "determinant" = item_information_determinant(),
                                "posterior_determinant" = item_plus_prior_information_determinant(),
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
  
  item_plus_prior_information_trace <- function() {
    information_administered <- apply(fisher_information[,,administered, drop = FALSE], c(1, 2), sum) + solve(prior_parameters$Sigma)
    apply(fisher_information[,,available, drop = FALSE], 3, function(x) sum(diag(information_administered + x)))
  }
  
  item_information_determinant <- function() {
    information_administered <- apply(fisher_information[,,administered, drop = FALSE], c(1, 2), sum)
    apply(fisher_information[,,available, drop = FALSE], 3, function(x) det(information_administered + x))
  }

  item_plus_prior_information_determinant <- function() {
    information_administered <- apply(fisher_information[,,administered, drop = FALSE], c(1, 2), sum) + solve(prior_parameters$Sigma)
    apply(fisher_information[,,available, drop = FALSE], 3, function(x) det(information_administered + x))
  }

  item_information_pekl <- function() {
    get_posterior_expected_kl_information(estimate, model, answers, administered, available, number_dimensions, estimator, alpha, beta, guessing, prior_form, prior_parameters, number_itemsteps_per_item, eap_estimation_procedure = eap_estimation_procedure)
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
