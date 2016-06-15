#' Returns a list with the key of the next item to be administered given a new answer,
#' an updated estimate of theta, updated covariance matrix of theta converted to a vector, and the answers to the administered items
#'
#' @param answers named list of previous answers and new answer, with names being the item keys; should be initialized with NULL
#' @param estimate vector with estimate of latent trait theta
#' @param variance (co)variance matrix of the estimate, as vector
#' @param model One of '3PLM', 'GPCM', 'SM' or 'GRM', for the three-parameter logistic, generalized partial credit, sequential or graded response model respectively.
#' @param alpha Matrix of alpha parameters, one column per dimension, one row per item. Row names should contain the item keys. Note that so called within-dimensional models still use an alpha matrix, they simply 
#' have only one non-zero loading per item.
#' @param beta Matrix of beta parameters, one column per item step, one row per item. Row names should contain the item keys. Note that ShadowCAT expects answer categories to be sequential,
#' and without gaps. That is, the weight parameter in the GPCM model is assumed to be sequential, and equal to the position of the 'location' of the beta parameter in the Beta matrix.
#' The matrix will have a number of columns equal to the largest number of item steps over items, items with fewer answer categories should be 
#' right-padded with \code{NA}. \code{NA} values between answer categories are not allowed, and will lead to errors.
#' Beta matrix can be set to NULL if model is GPCM and eta is defined
#' More flexibility in Beta parameters might be added in future versions.
#' @param start_items List indicating the items that are shown to the respondent before adaptive proces starts; one of
#' list(type = 'random', n)
#' list(type = 'fixed', item_keys, n)
#' list(type = 'random_by_dimension', n_by_dimension, n)
#' where n = total number of initial items, item_keys = character vector with keys of initial items, 
#' n_by_dimension = scalar of number of initial items per dimension, or vector with number of initial items for each dimension
#' If n is 0, only n needs to be defined
#' 'random_by_dimension' assumes that items load on a single dimension, if any item has a non-zero loading on a dimension, it is considered to be part of that dimension. 
#' @param stop_test List indicating rules for when to stop providing new items to respondent; should be a list of the form
#' list(target = ..., max_n = ..., min_n = ..., cutoffs = ...), 
#' where max_n = test length at which testing should stop (even if target has not been reached yet in case of variance stopping rule), 
#' target = vector of maximum acceptable variances per dimension; NULL means no variance target,
#' min_n = minimum test length; NULL means no mimimum test length,
#' cutoffs = matrix containing cut off values per dimension (columns) and test iteration (rows). First row contains cut off values for when no items have been
#' administered yet, second row for when one item has been administered, etc. If estimate + 3SE < cutoff for each dimension at certain iteration, test stops; 
#' NULL means no cut off values
#' @param estimator Type of estimator to be used, one of "maximum_aposteriori", "maximum_likelihood", or "expected_aposteriori"
#' @param information_summary How to summarize Fisher information, used for selection of item with maximum information; one of
#' "determinant": compute determinant(info_sofar_QxQ + info_QxQ_k) for each yet available item k
#' "posterior_determinant": compute determinant(info_sofar_QxQ_plus_prior_information + info_QxQ_k) for each yet available item k
#' "trace": compute trace((info_sofar_QxQ + info_QxQ_k) for each yet available item k
#' "posterior_trace": compute trace(info_sofar_QxQ_plus_prior_information + info_QxQ_k) for each yet available item k
#' "posterior_expected_kullback_leibler" = compute Posterior expected Kullback-Leibler Information
#' @param prior_form String indicating the form of the prior; one of "normal" or "uniform". Not required if estimator is maximum_likelihood
#' @param prior_parameters List containing mu and Sigma of the normal prior: list(mu = ..., Sigma = ...), or 
#' the upper and lower bound of the uniform prior: list(lower_bound = ..., upper_bound = ...). . Not required if estimator is maximum_likelihood.
#' Sigma should always be in matrix form. The length of lower_bound and upper_bound should be equal to the number of dimensions.
#' For uniform prior in combination with expected aposteriori estimation, true theta should fall within 
#' lower_bound and upper_bound and be not too close to one of these bounds. Setting the function argument safe_eap to TRUE
#' ensures that the estimation switches to maximum aposteriori if the expected aposteriori estimate fails. 
#' @param guessing Matrix with one column of guessing parameters per item. Row names should contain the item keys. Optionally used in 3PLM model, ignored for all others.
#' @param eta Matrix of location parameters, optionally used in GPCM model, ignored for all others. Row names should contain the item keys.
#' @param constraints_and_characts list with constraints and characteristics; NULL means no constraints
#' constraints should be specified as a list of constraints, each constraint is a list with three named values;
#' name: the column name of the characteristic this constraint applies to. For categorical characteristics the level should be specified as name/value.
#' op: the logical operator to be used. Valid options are "<", "=", ">" and "><".
#' target: the target value, numeric. If the operator is "><", this should be a length two vector in between which the target should fall.
#' characteristics should be a data.frame with characteristics, one row per item (rows in the same order as they are in alpha, beta, etc.), one column per characteristic.
#' See constraints_lp_format() for details
#' @param lower_bound Vector with lower bounds for theta per dimension; estimated theta values smaller than the lowerbound values are truncated to the lowerbound values.
#' Can only be defined when estimator is maximum_likelihood. Setting bounds with maximum likelihood estimation is equivalent to
#' using maximum aposteriori estimation with a uniform prior. 
#' @param upper_bound Vector with upper bounds for theta per dimension; estimated theta values larger than the upperbound values are truncated to the upperbound values
#' Can only be defined when estimator is maximum_likelihood. Setting bounds with maximum likelihood estimation is equivalent to
#' using maximum aposteriori estimation with a uniform prior.
#' @param safe_eap Only relevant if estimator is expected_aposteriori. 
#' TRUE if estimator should switch to maximum aposteriori if the integration algorithm results in an error.
#' An error may occur if the prior is uniform, estimator is expected aposteriori, and the bounds do not exceed the true theta value, or are too close to it.
#' @param eap_estimation_procedure String indicating the estimation procedure if estimator is expected aposteriori and prior form is normal. One of "riemannsum" for integration via Riemannsum or
#' "gauss_hermite_quad" for integration via Gaussian Hermite Quadrature. If prior form is uniform, estimation procedure should always be "riemannsum".
#' @return List containing the key of the next item to be administered given a new answer, an indicator for whether the test should be continued,
#' updated estimate of theta, updated covariance matrix of theta converted to a vector, and the answers to the administered items (named list)
#' @examples
#' alpha_beta <- simulate_testbank(model = "GPCM", number_items = 100, number_dimensions = 3, number_itemsteps = 3)
#' model <- "GPCM"
#' start_items <- list(type = 'fixed', item_keys = c("item33", "item5", "item23"), n = 3)
#' stop_test <- list(min_n = 4, max_n = 30, target = c(.1, .1, .1))
#' estimator <- "maximum_aposteriori"
#' information_summary <- "posterior_determinant"
#' prior_form <- "normal"
#' prior_parameters <- list(mu = c(0, 0, 0), Sigma = diag(3))
#' 
#' # Initial call: get key of first item to adminster
#' call1 <- shadowcat(answers = NULL, estimate = c(0, 0, 0), variance = as.vector(diag(3) * 25), model = model, alpha = alpha_beta$alpha, beta = alpha_beta$beta, start_items = start_items, stop_test = stop_test, estimator = estimator, information_summary = information_summary, prior_form = prior_form, prior_parameters = prior_parameters)
#' # Second to fourth call: number of start items is set to 3, so no update in theta estimate yet
#' call2 <- shadowcat(answers = list(item33 = 2), estimate = call1$estimate, variance = call1$variance, model = model, alpha = alpha_beta$alpha, beta = alpha_beta$beta, start_items = start_items, stop_test = stop_test, estimator = estimator, information_summary = information_summary, prior_form = prior_form, prior_parameters = prior_parameters)
#' call3 <- shadowcat(answers = list(item33 = 2, item5 = 3), estimate = call2$estimate, variance = call2$variance, model = model, alpha = alpha_beta$alpha, beta = alpha_beta$beta, start_items = start_items, stop_test = stop_test, estimator = estimator, information_summary = information_summary, prior_form = prior_form, prior_parameters = prior_parameters)
#' call4 <- shadowcat(answers = list(item33 = 2, item5 = 3, item23 = 3), estimate = call3$estimate, variance = call3$variance, model = model, alpha = alpha_beta$alpha, beta = alpha_beta$beta, start_items = start_items, stop_test = stop_test, estimator = estimator, information_summary = information_summary, prior_form = prior_form, prior_parameters = prior_parameters)
#' # Fifth call: first time theta estimate is updated
#' call5 <- shadowcat(answers = list(item33 = 2, item5 = 3, item23 = 3, item84 = 1), estimate = call4$estimate, variance = call4$variance, model = model, alpha = alpha_beta$alpha, beta = alpha_beta$beta, start_items = start_items, stop_test = stop_test, estimator = estimator, information_summary = information_summary, prior_form = prior_form, prior_parameters = prior_parameters)
#' # Sixth call: use the updated estimate and variance as the current values for estimate and variance
#' call6 <- shadowcat(answers = list(item33 = 2, item5 = 3, item23 = 3, item84 = 1, item36 = 2), estimate = call5$estimate, variance = call5$variance, model = model, alpha = alpha_beta$alpha, beta = alpha_beta$beta, start_items = start_items, stop_test = stop_test, estimator = estimator, information_summary = information_summary, prior_form = prior_form, prior_parameters = prior_parameters)
#' 
#' # With constraints
#' 
#' constraints_and_characteristics <- list(characteristics = data.frame(content = sample(c('algebra','physics','calculus'), size = 100, replace = TRUE),
#'                                                                      time = runif(100),
#'                                                                      exclusive = sapply(1:100, FUN = function (x) { if (x %in% sample(1:100, size = 4)) 1 else 0 })),
#'                                         constraints = list(list(name = 'content/algebra',
#'                                                                 op = '><',
#'                                                                 target = c(5, 10)), # ensure number of algebra items is between 5 and 10
#'                                                            list(name = 'content/physics',
#'                                                                 op = '><',
#'                                                                 target = c(2, 5)), # ensure number of physics items is between 2 and 5
#'                                                            list(name = 'time',
#'                                                                 op = '<',
#'                                                                 target = 20), # Ensure total tests takes no longer than 20 minutes
#'                                                            list(name = 'exclusive',
#'                                                                 op = '<',
#'                                                                 target = 2))) # Ensure number of exclusive items equals 2
#'                                                              
#' shadowcat(answers = NULL, estimate = c(0, 0, 0), variance = as.vector(diag(3) * 25), model = model, alpha = alpha_beta$alpha, beta = alpha_beta$beta, start_items = start_items, stop_test = stop_test, estimator = estimator, information_summary = information_summary, prior_form = prior_form, prior_parameters = prior_parameters, constraints_and_characts = constraints_and_characteristics)
#' @importFrom matrixcalc is.positive.definite
#' @export
shadowcat <- function(answers, estimate, variance, model, alpha, beta, start_items, stop_test, estimator, information_summary, prior_form = NULL, prior_parameters = NULL, guessing = NULL, eta = NULL, constraints_and_characts = NULL, lower_bound = NULL, upper_bound = NULL, safe_eap = FALSE, eap_estimation_procedure = "riemannsum") {      
  result <- function() {
    switch_to_maximum_aposteriori <- estimator == "maximum_likelihood" && !is.null(lower_bound) && !is.null(upper_bound)
    estimator <- get_estimator(switch_to_maximum_aposteriori = switch_to_maximum_aposteriori)
    prior_form <- get_prior_form(switch_to_maximum_aposteriori = switch_to_maximum_aposteriori)
    prior_parameters <- get_prior_parameters(switch_to_maximum_aposteriori = switch_to_maximum_aposteriori)
    beta <- get_beta()
    guessing <- get_guessing() 
    number_items <- nrow(alpha)
    number_dimensions <- ncol(alpha)
    number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
    lp_constraints_and_characts <- get_lp_constraints_and_characts(number_items = number_items)
    item_keys <- rownames(alpha)
    item_keys_administered <- names(answers)
    item_keys_available <- get_item_keys_available(item_keys_administered = item_keys_administered, item_keys = item_keys)
    attr(estimate, "variance") <- matrix(variance, ncol = number_dimensions)    
    
    estimate <- update_person_estimate(estimate = estimate, answers_vector = unlist(answers), item_indices_administered = match(item_keys_administered, item_keys), number_dimensions = number_dimensions, alpha = alpha, beta = beta, guessing = guessing, number_itemsteps_per_item = number_itemsteps_per_item, estimator = estimator, prior_form = prior_form, prior_parameters = prior_parameters)
    continue_test <- !test_must_stop(number_answers = length(answers), estimate = estimate, min_n = stop_test$min_n, max_n = stop_test$max_n, stop_variance_target = stop_test$target, cutoffs = stop_test$cutoffs)
    if (continue_test) {
      index_new_item <- get_next_item(start_items = start_items, information_summary = information_summary, lp_constraints = lp_constraints_and_characts$lp_constraints, lp_characters = lp_constraints_and_characts$lp_chars, estimate = estimate, model = model, answers = unlist(answers), prior_form = prior_form, prior_parameters = prior_parameters, 
                                      available = match(item_keys_available, item_keys), administered = match(item_keys_administered, item_keys), number_items = number_items, number_dimensions = number_dimensions, estimator = estimator, alpha = alpha, beta = beta, guessing = guessing, number_itemsteps_per_item = number_itemsteps_per_item, stop_test = stop_test, eap_estimation_procedure = eap_estimation_procedure)
      key_new_item <- item_keys[index_new_item]
    }
    else {
      key_new_item <- NULL
    }
    
    list(key_new_item = as.scalar2(key_new_item),
         continue_test = as.scalar2(continue_test),
         estimate = as.vector(estimate),
         variance = as.vector(attr(estimate, "variance")),
         answers = answers)
  }
  
  # if inititial items have been administered (so we are in the CAT phase), update person estimate after each newly answered item
  update_person_estimate <- function(estimate, answers_vector, item_indices_administered, number_dimensions, alpha, beta, guessing, number_itemsteps_per_item, estimator, prior_form, prior_parameters) { 
    if (length(answers) > start_items$n)
      estimate_latent_trait(estimate = estimate, answers = answers_vector, prior_form = prior_form, prior_parameters = prior_parameters, model = model, administered = item_indices_administered, number_dimensions = number_dimensions, estimator = estimator, alpha = alpha, beta = beta, guessing = guessing, number_itemsteps_per_item = number_itemsteps_per_item, safe_eap = safe_eap, eap_estimation_procedure = eap_estimation_procedure)
    else
      estimate
  }
  
  get_item_keys_available <- function(item_keys_administered, item_keys) {
    if (is.null(item_keys_administered))
      item_keys
    else
      item_keys[-which(item_keys %in% item_keys_administered)]
  }
  
  get_beta <- function() {
    # allow calculating beta from eta.
    if (model == "GPCM" && is.null(beta) && !is.null(eta))
      row_cumsum(eta)
    else
      beta
  }
  
  get_guessing <- function() {
    if (is.null(guessing))
      matrix(0, nrow = nrow(as.matrix(alpha)), ncol = 1, dimnames = list(rownames(alpha), NULL))
    else
      guessing
  }
  
  get_estimator <- function(switch_to_maximum_aposteriori) {
    if (switch_to_maximum_aposteriori)
      "maximum_aposteriori"
    else
      estimator
  }
  
  get_prior_form <- function(switch_to_maximum_aposteriori) {
    if (switch_to_maximum_aposteriori)
      "uniform"
    else
      prior_form
  }
  
  get_prior_parameters <- function(switch_to_maximum_aposteriori) {
    if (switch_to_maximum_aposteriori)
      list(lower_bound = lower_bound, upper_bound = upper_bound)
    else
      prior_parameters
  }
  
  get_lp_constraints_and_characts <- function(number_items) {
    if (is.null(constraints_and_characts))
      NULL
    else
      constraints_lp_format(max_n = stop_test$max_n, number_items = number_items, characteristics = constraints_and_characts$characteristics, constraints = constraints_and_characts$constraints) 
  }
    
  validate <- function() {
    if (is.null(estimate))
      return(add_error("estimate", "is missing"))
    if (is.null(variance))
      return(add_error("variance", "is missing"))
    if (!is.vector(variance))
      return(add_error("variance", "should be entered as vector"))
    if (sqrt(length(variance)) != round(sqrt(length(variance))))
      return(add_error("variance", "should be a covariance matrix turned into a vector"))
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
    if (!is.null(start_items$type) && start_items$type == "random_by_dimension" && length(start_items$n_by_dimension) %not_in% c(1, length(estimate)))
      return(add_error("start_items", "length of n_by_dimension should be a scalar or vector of the length of estimate"))
    if (!row_names_are_equal(rownames(alpha), list(alpha, beta, eta, guessing)))
      add_error("alpha_beta_eta_guessing", "should have equal row names, in same order")
    if (!is.null(beta) && !na_only_end_rows(beta))
      add_error("beta", "can only contain NA at the end of rows, no values allowed after an NA in a row")
    if (!is.null(eta) && !na_only_end_rows(eta))
      add_error("eta", "can only contain NA at the end of rows, no values allowed after an NA in a row")
    if (length(estimate) != ncol(alpha))
      add_error("estimate", "length should be equal to the number of columns of the alpha matrix")
    if (length(estimate)^2 != length(variance))
      add_error("variance", "should have a length equal to the length of estimate squared")
    if (is.null(answers) && !is.positive.definite(matrix(variance, ncol = sqrt(length(variance)))))
      add_error("variance", "matrix is not positive definite")
    if (model %not_in% c("3PLM", "GPCM", "SM", "GRM"))
      add_error("model", "of unknown type")
    if (model != "GPCM" && is.null(beta))
      add_error("beta", "is missing")
    if (model == "GPCM" && is.null(beta) && is.null(eta))
      add_error("beta_and_eta", "are both missing; define at least one of them")
    if (model == "GPCM" && !is.null(beta) && !is.null(eta) && !all(row_cumsum(eta) == beta))
      add_error("beta_and_eta", "objects do not match")
    if (estimator != "maximum_likelihood" && is.null(prior_form))
      add_error("prior_form", "is missing")
    if (estimator != "maximum_likelihood" && is.null(prior_parameters))
      add_error("prior_parameters", "is missing")
    if (!is.null(prior_form) && prior_form %not_in% c("normal", "uniform"))
      add_error("prior_form", "of unknown type")
    if (!is.null(prior_form) && !is.null(prior_parameters) && prior_form == "uniform" && (is.null(prior_parameters$lower_bound) || is.null(prior_parameters$upper_bound)))
      add_error("prior_form_is_uniform", "so prior_parameters should contain lower_bound and upper_bound")
    if (!is.null(prior_form) && !is.null(prior_parameters) && prior_form == "normal" && (is.null(prior_parameters$mu) || is.null(prior_parameters$Sigma)))
      add_error("prior_form_is_normal", "so prior_parameters should contain mu and Sigma")
    if (!is.null(prior_parameters$mu) && length(prior_parameters$mu) != length(estimate))
      add_error("prior_parameters_mu", "should have same length as estimate")    
    if (!is.null(prior_parameters$Sigma) && (!is.matrix(prior_parameters$Sigma) || !all(dim(prior_parameters$Sigma) == c(length(estimate), length(estimate))) || !is.positive.definite(prior_parameters$Sigma)))
        add_error("prior_parameters_sigma", "should be a square positive definite matrix, with dimensions equal to the length of estimate")
    if (!is.null(prior_parameters$lower_bound) && !is.null(prior_parameters$upper_bound) && (length(prior_parameters$lower_bound) != length(estimate) || length(prior_parameters$upper_bound) != length(estimate)))
      add_error("prior_parameters_bounds", "should contain lower and upper bound of the same length as estimate")
    if (is.null(stop_test$max_n))
      add_error("stop_test", "contains no max_n")
    if (!is.null(stop_test$max_n) && stop_test$max_n > nrow(alpha))
      add_error("stop_test_max_n", "is larger than the number of items in the item bank")
    if (!is.null(stop_test$max_n) && !is.null(stop_test$cutoffs) && (!is.matrix(stop_test$cutoffs) || nrow(stop_test$cutoffs) < stop_test$max_n || ncol(stop_test$cutoffs) != length(estimate) || any(is.na(stop_test$cutoffs))))
      add_error("stop_test_cutoffs", "should be a matrix without missing values, and number of rows equal to max_n and number of columns equal to the number of dimensions")
    if (start_items$n == 0 && information_summary == "posterior_expected_kullback_leibler")
      add_error("start_items", "requires n > 0 for posterior expected kullback leibler information summary")
    if (!is.null(start_items$type) && start_items$type == "random_by_dimension" && length(start_items$n_by_dimension) == length(estimate) && start_items$n != sum(start_items$n_by_dimension))
      add_error("start_items_n", "contains inconsistent information. Total length of start phase and sum of length per dimension do not match (n != sum(n_by_dimension)")
    if (!is.null(start_items$type) && start_items$type == "random_by_dimension" && length(start_items$n_by_dimension) == 1 && start_items$n != sum(rep(start_items$n_by_dimension, length(estimate))))
      add_error("start_items_n", "contains inconsistent information. Total length of start phase and sum of length per dimension do not match")
    if (!is.null(stop_test$cutoffs) && !is.matrix(stop_test$cutoffs))
      add_error("stop_test", "contains cutoff values in non-matrix format")
    if (!all(names(answers) %in% rownames(alpha)))
      add_error("answers", "contains non-existing key")
    if (estimator %not_in% c("maximum_likelihood", "maximum_aposteriori", "expected_aposteriori"))
      add_error("estimator", "of unknown type")
    if (information_summary %not_in% c("determinant", "posterior_determinant", "trace", "posterior_trace", "posterior_expected_kullback_leibler"))
      add_error("information_summary", "of unknown type")
    if (estimator == "maximum_likelihood" && information_summary %in% c("posterior_determinant", "posterior_trace", "posterior_expected_kullback_leibler"))
      add_error("estimator_is_maximum_likelihood", "so using a posterior information summary makes no sense")
    if (estimator != "maximum_likelihood" && (!is.null(lower_bound) || !is.null(upper_bound)))
      add_error("bounds", "can only be defined if estimator is maximum likelihood")
    if (!is.null(lower_bound) && length(lower_bound) %not_in% c(1, length(estimate)))
      add_error("lower_bound", "length of lower bound should be a scalar or vector of the length of estimate")
    if (!is.null(upper_bound) && length(upper_bound) %not_in% c(1, length(estimate)))
      add_error("upper_bound", "length of upper bound should be a scalar or vector of the length of estimate")
    if (!no_missing_information(constraints_and_characts$characteristics, constraints_and_characts$constraints))
      add_error("constraints_and_characts", "constraints and characteristics should either be defined both or not at all")
    if (!characteristics_correct_format(constraints_and_characts$characteristics, number_items = nrow(alpha)))
      add_error("characteristics", "should be a data frame with number of rows equal to the number of items in the item bank")
    if (!constraints_correct_structure(constraints_and_characts$constraints))
      add_error("constraints_structure", "should be a list of length three lists, with elements named 'name', 'op', 'target'")
    if (!constraints_correct_names(constraints_and_characts$constraints, constraints_and_characts$characteristics))
      add_error("constraints_name_elements", "should be defined as described in the details section of constraints_lp_format()")
    if (!constraints_correct_operators(constraints_and_characts$constraints))
      add_error("constraints_operator_elements", "should be defined as described in the details section of constraints_lp_format()")
    if (!constraints_correct_targets(constraints_and_characts$constraints))
      add_error("constraints_target_elements", "should be defined as described in the details section of constraints_lp_format()")
  }
  
  invalid_result <- function() {
    list(errors = errors())
  }
  
  validate_and_run()
}

