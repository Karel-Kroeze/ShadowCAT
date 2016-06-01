#' Obtain a vector with information, summarized into one value, for each available item (values of the objective function).
#' 
#' Information of the test so far (including all administered items) is added to the information of the available
#' item before the summary is computed 
#'
#' @param information_summary called "objective" by Kroeze; how to summarize information; one of
#' "determinant": compute determinant(info_sofar_QxQ + info_QxQ_k) for each yet available item k
#' "posterior_determinant": compute determinant(info_sofar_QxQ_plus_prior_information + info_QxQ_k) for each yet available item k
#' "trace": compute trace((info_sofar_QxQ + info_QxQ_k) for each yet available item k
#' "posterior_trace": compute trace(info_sofar_QxQ_plus_prior_information + info_QxQ_k) for each yet available item k
#' "posterior_expected_kullback_leibler" = compute Posterior expected Kullback-Leibler Information
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
#' @param estimator type of estimator to be used, one of "maximum_aposteriori", "maximum_likelihood", or "expected_aposteriori"
#' @param alpha matrix containing the alpha parameters
#' @param beta matrix containing the beta parameters
#' @param guessing matrix containing the quessing parameters
#' @param number_itemsteps_per_item vector containing the number of non missing cells per row of the beta matrix
#' @param pad Should the return vector be padded with zeros for items that have already been administered?
#' @param eap_estimation_procedure String indicating the estimation procedure for the expected aposteriori estimate, which is computed
#' in get_posterior_expected_kl_information() if it is not the requested estimator in shadowcat(). One of "riemannsum" for integration via Riemannsum or
#' "gauss_hermite_quad" for integration via Gaussian Hermite Quadrature. Only important here if information_summary is posterior_expected_kl_information.
#' @return vector with information for each available item
#' @export
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
