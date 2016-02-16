#' Returns a list with the key of the next item to be administered given a new response,
#' an updated estimate of theta, and the responses to the administered items
#'
#' @param responses named list of previous responses and new response, with names being the item keys; should be initialized with NULL
#' @param estimate estimate of latent trait theta, with covariance matrix as its attribute
#' @param model String, one of '3PLM', 'GPCM', 'SM' or 'GRM', for the three-parameter logistic, generalized partial credit, sequential or graded response model respectively.
#' @param alpha Matrix of alpha parameters, one column per dimension, one row per item. Row names should contain the item keys. Note that so called within-dimensional models still use an alpha matrix, they simply 
#' have only one non-zero loading per item.
#' @param beta Matrix of beta parameters, one column per item step, one row per item. Row names should contain the item keys. Note that ShadowCAT expects response categories to be sequential,
#' and without gaps. That is, the weight parameter in the GPCM model is assumed to be sequential, and equal to the position of the 'location' of the beta parameter in the Beta matrix.
#' The matrix will have a number of columns equal to the largest number of response categories, items with fewer response categories should be 
#' right-padded with \code{NA}. \code{NA} values between response categories are not allowed, and will lead to errors.
#' Beta matrix can be set to NULL if model is GPCM and eta is defined
#' More flexibility in Beta parameters might be added in future versions.
#' @param start_items items that are shown to the patient before adaptive proces starts; one of
#' list(type = 'random', n)
#' list(type = 'fixed', indeces, n)
#' list(type = 'random_by_dimension', n_by_dimension, n)
#' where n = total number of initial items, indeces = vector of initial item indeces, 
#' n_by_dimension = scalar of number of initial items per dimension, or vector with number of initial items for each dimension
#' If n is 0, only n needs to be defined
#' 'random_by_dimension' assumes that items load on a single dimension, if any item has a non-zero loading on a dimension, it is considered to be part of that dimension. 
#' @param stop_test rule for when to stop providing new items to patient; should be a list of the form
#' list(target = ..., max_n = ..., min_n = ..., cutoffs = ...), 
#' where max_n = test length at which testing should stop (even if target has not been reached yet in case of variance stopping rule), 
#' target = vector of maximum acceptable variances per dimension; if target = NULL, only max_n is taken into account,
#' min_n = minimum test length; NULL means no mimimum test length,
#' cutoffs = matrix containing cut off values per dimension (columns) and test iteration (rows). First row contains cut off values for when no items have been
#' administered yet, second row for when one item has been administered, etc. If estimate + 3SE < cutoff for each dimension at certain iteration, test stops; 
#' NULL means no cut off values
#' @param estimator type of estimator to be used, one of "maximum_aposteriori", "maximum_likelihood", or "expected_aposteriori"
#' @param information_summary called "objective" by Kroeze; how to summarize information; one of
#' "determinant": compute determinant(info_sofar_QxQ + info_QxQ_k) for each yet available item k
#' "posterior_determinant": compute determinant(info_sofar_QxQ_plus_prior + info_QxQ_k) for each yet available item k
#' "trace": compute trace((info_sofar_QxQ + info_QxQ_k) for each yet available item k
#' "posterior_trace": compute trace(info_sofar_QxQ_plus_prior + info_QxQ_k) for each yet available item k
#' "posterior_expected_kullback_leibler" = compute Posterior expected Kullback-Leibler Information
#' @param prior covariance matrix of the (multi variate) normal prior for theta; mean vector is fixed at zero; not used for maximum likelihood estimator
#' #' note that this prior should be a square matrix with number of rows and columns equal to the number of dimensions; values on the diagonal should be larger than 0
#' @param guessing matrix with one column of guessing parameters per item. Row names should contain the item keys. Optionally used in 3PLM model, ignored for all others.
#' @param eta Matrix of location parameters, optionally used in GPCM model, ignored for all others. Row names should contain the item keys.
#' @param constraints_and_characts list with constraints and characteristics; NULL means no constraints
#' constraints should be specified as a list of constraints, each constraint is a list with three named values;
#' name: the column name of the characteristic this constraint applies to. For categorical characteristics the level should be specified as name/value.
#' op: the logical operator to be used. Valid options are "<", "=", ">" and "><".
#' target: the target value, numeric. If the operator is "><", this should be a length two vector in between which the target should fall.
#' characteristics should be a data.frame with characteristics, one row per item (rows in the same order as they are in alpha, beta, etc.), one column per characteristic.
#' See constraints_lp_format() for details
#' @param lower_bound vector with lower bounds for theta per dimension; estimated theta values smaller than the lowerbound values are truncated to the lowerbound values 
#' @param upper_bound vector with upper bounds for theta per dimension; estimated theta values larger than the upperbound values are truncated to the upperbound values
#' @param prior_var_safe_ml if not NULL, expected a posteriori estimate with prior variance equal to prior_var_safe_ml (scalar or vector) is computed instead of maximum likelihood/maximum a posteriori, if maximum likelihood/maximum a posteriori estimate fails
#' @return a list containing the key of the next item to be administered given a new response (or "stop_test"), 
#' updated estimate of theta, and the responses to the administered items (named list)
#' @export
shadowcat <- function(responses, estimate, model, alpha, beta, start_items, stop_test, estimator, information_summary, prior = NULL, guessing = NULL, eta = NULL, constraints_and_characts = NULL, lower_bound = rep(-3, ncol(alpha)), upper_bound = rep(3, ncol(alpha)), prior_var_safe_ml = NULL) {      
  result <- function() {
    beta <- get_beta(model, beta, eta)
    guessing <- get_guessing(guessing, beta) 
    number_items <- nrow(beta)
    number_dimensions <- ncol(alpha)
    number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
    lp_constraints_and_characts <- constraints_lp_format(stop_test$max_n, number_items, constraints_and_characts$characteristics, constraints_and_characts$constraints)
    item_keys <- rownames(alpha)
    item_keys_administered <- names(responses)
    item_keys_available <- get_item_keys_available(item_keys_administered, item_keys)
    
    estimate <- update_person_estimate(estimate, unlist(responses), match(item_keys_administered, item_keys), number_dimensions, alpha, beta, guessing, number_itemsteps_per_item)
    continue_test <- !test_must_stop(length(responses), estimate, stop_test$min_n, stop_test$max_n, stop_test$target, stop_test$cutoffs)
    if (continue_test) {
      index_new_item <- get_next_item(start_items, information_summary, lp_constraints_and_characts$lp_constraints, lp_constraints_and_characts$lp_chars, estimate, model, unlist(responses), prior, match(item_keys_available, item_keys), match(item_keys_administered, item_keys), number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound)
      key_new_item <- item_keys[index_new_item]
    }
    else {
      key_new_item <- "stop_test"
    }
    
    list(key_new_item = key_new_item,
         estimate = estimate,
         responses = responses)
  }
  
  # if inititial items have been administered (so we are in the CAT phase), update person estimate after each newly answered item
  update_person_estimate <- function(estimate, responses_vector, item_indeces_administered, number_dimensions, alpha, beta, guessing, number_itemsteps_per_item) { 
    if (length(responses) > start_items$n)
      estimate_latent_trait(estimate, responses_vector, prior, model, item_indeces_administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, prior_var_safe_ml)
    else
      estimate
  }
  
  get_item_keys_available <- function(item_keys_administered, item_keys) {
    if (is.null(item_keys_administered))
      item_keys
    else
      item_keys[-which(item_keys %in% item_keys_administered)]
  }
  
  validate <- function() {
    if (is.null(estimate))
      return(add_error("estimate", "is missing"))
    if (is.null(attr(estimate, "variance")))
      return(add_error("variance", "is missing as an attribute of estimate"))
    if (is.null(model))
      return(add_error("model", "is missing"))
    if (is.null(alpha))
      return(add_error("alpha", "is missing"))
    if (is.null(start_items))
      return(add_error("start_items", "is missing"))
    if (is.null(stop_test))
      return(add_error("stop_test", "is missing"))
    if (is.null(estimator))
      return(add_error("estimator", "is missing"))
    if (is.null(information_summary))
      return(add_error("information_summary", "is missing"))
    if (!is.matrix(alpha) || is.null(rownames(alpha)))
      return(add_error("alpha", "should be a matrix with item keys as row names"))
    if (!is.null(beta) && (!is.matrix(beta) || is.null(rownames(beta))))
      return(add_error("beta", "should be a matrix with item keys as row names"))
    if (!is.null(eta) && (!is.matrix(eta) || is.null(rownames(eta))))
      return(add_error("eta", "should be a matrix with item keys as row names"))
    if (!is.null(guessing) && (!is.matrix(guessing) || ncol(guessing) != 1 || is.null(rownames(guessing))))
      return(add_error("guessing", "should be a single column matrix with item keys as row names"))
    if (!row_names_are_equal(rownames(alpha), list(alpha, beta, eta, guessing)))
      add_error("alpha_beta_eta_guessing", "should have equal row names, in same order")
    if (model != "GPCM" && is.null(beta))
      add_error("beta", "is missing")
    if (model == "GPCM" && is.null(beta) && is.null(eta))
      add_error("beta_and_eta", "are both missing; define at least one of them")
    if (model == "GPCM" && !is.null(beta) && !is.null(eta) && !all(row_cumsum(eta) == beta))
      add_error("beta_and_eta", "objects do not match")
    if ((estimator != "maximum_likelihood" || information_summary %in% c("posterior_determinant", "posterior_trace")) && is.null(prior))
      add_error("prior", "is missing")
    if (is.null(stop_test$max_n))
      add_error("stop_test", "contains no max_n")
    if (start_items$n == 0 && information_summary == "posterior_expected_kullback_leibler")
      add_error("start_items", "requires n > 0 for posterior expected kullback leibler information summary")
    if (!is.null(stop_test$cutoffs) && !is.matrix(stop_test$cutoffs))
      add_error("stop_test", "contains cutoff values in non-matrix format")
  }
  
  invalid_result <- function() {
    list(errors = errors())
  }
  
  validate_and_run()
}

#' get beta matrix from beta or eta
#' 
#' @param model String, one of '3PLM', 'GPCM', 'SM' or 'GRM', for the three-parameter logistic, generalized partial credit, sequential or graded response model respectively.
#' @param beta matrix containing beta parameters, may be NULL if model is GPCM and eta is defined
#' @param eta matrix containing eta parameters, may be NULL if beta is defined
#' @return beta matrix obtained from beta or eta
#' @export
get_beta <- function(model, beta, eta) {
  # allow calculating beta from eta.
  if (model == "GPCM" && is.null(beta) && !is.null(eta))
    row_cumsum(eta)
  else
    beta
}

#' get guessing matrix
#' 
#' @param guessing vector containing guessing parameters; may be NULL in case of zero guessing parameters
#' @param beta matrix containing beta parameters
#' @return matrix containing guessing parameters
#' @export
get_guessing <- function(guessing, beta) {
  if (is.null(guessing))
    matrix(0, nrow(as.matrix(beta)), 1, dimnames = list(rownames(beta), NULL))
  else
    guessing
}
