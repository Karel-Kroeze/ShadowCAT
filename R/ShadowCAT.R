#' Returns a list with the index of the next item to be administered given a new response, and an updated person object
#'
#' @param new_response new response from respondent, should be initialized with NULL
#' @param estimate estimate of latent trait theta
#' @param responses vector of given responses; numeric(0) at first iteration
#' @param administered vector containing indeces of administered items; numeric(0) at first iteration
#' @param available vector containing indeces of yet available items
#' @param prior covariance matrix of the multi variate normal prior for theta; mean vector is fixed at zero; only used when estimator type is MAP or EAP, but at this point should always be defined
#' #' note that this prior should be a square matrix with number of rows and columns equal to the number of dimensions; values on the diagonal should be larger than 1
#' @param model String, one of '3PLM', 'GPCM', 'SM' or 'GRM', for the three-parameter logistic, generalized partial credit, sequential or graded response model respectively.
#' @param alpha Matrix of alpha parameters, one column per dimension, one row per item. Note that so called within-dimensional models still use an alpha matrix, they simply 
#' have only one non-zero loading per item.
#' @param beta Matrix of beta parameters, one column per item step, one row per item. Note that ShadowCAT expects response categories to be sequential,
#' and without gaps. That is, the weight parameter in the GPCM model is assumed to be sequential, and equal to the position of the 'location' of the beta parameter in the Beta matrix.
#' The matrix will have a number of columns equal to the largest number of response categories, items with fewer response categories should be 
#' right-padded with \code{NA}. \code{NA} values between response categories are not allowed, and will lead to errors.
#' More flexibility in Beta parameters might be added in future versions.
#' @param guessing vector of guessing parameters per item. Optionally used in 3PLM model, ignored for all others.
#' @param eta Matrix of location parameters, optionally used in GPCM model, ignored for all others.
#' @param start_items items that are shown to the patient before adaptive proces starts; one of
#' list(type = 'random', n)
#' list(type = 'fixed', indices, n)
#' list(type = 'randomByDimension', nByDimension, n)
#' where n = total number of initial items, indices = vector of initial item indeces, 
#' nByDimension = scalar of number of initial items per dimension, or vector with number of initial items for each dimension
#' @param stop_test rule for when to stop providing new items to patient; should be a list of the form
#' list(target = ..., n = ...), 
#' where n = test length at which testing should stop (even if target has not been reached yet in case of variance stopping rule; n <= 0 never stops for max length), 
#' target = vector of maximum acceptable variances per dimension; if target = NULL, only n is taken into account
#' @param estimator type of estimator to be used, one of "MAP" (Maximum a posteriori estimation) or "ML" (maximum likelihood); 
#' "EAP" (Expected A Posteriori Estimation) is currently not working due to problems with the MultiGHQuad package
#' @param information_summary called "objective" by Kroeze; how to summarize information; one of
#' "D" = determinant: compute determinant(info_sofar_QxQ + info_QxQ_k) for each yet available item k
#' "PD" = posterior determinant: compute determinant(info_sofar_QxQ_plus_prior + info_QxQ_k) for each yet available item k
#' "A" = trace: compute trace((info_sofar_QxQ + info_QxQ_k) for each yet available item k
#' "PA" = posterior trace: compute trace(info_sofar_QxQ_plus_prior + info_QxQ_k) for each yet available item k
#' "PEKL" = compute Posterior expected Kullback-Leibler Information
#' @param item_selection selection criterion; one of "MI" (maximum information) or "Shadow" (maximum information and take constraints into account)
#' @param constraints_and_characts list with constraints and characteristics
#' constraints should be specified as a list of constraints, each constraint is a list with three named values;
#' name: the column name of the characteristic this constraint applies to. For categorical characteristics the level should be specified as name/value.
#' op: the logical operator to be used. Valid options are "<", "=", ">" and "><".
#' target: the target value, numeric. If the operator is "><", this should be a length two vector in between which the target should fall.
#' characteristics should be a data.frame with characteristics, one row per item, one column per characteristic.
#' See constraints_correct_format() for details
#' @param lower_bound vector with lower bounds for theta per dimension; estimated theta values smaller than the lowerbound values are truncated to the lowerbound values 
#' @param upper_bound vector with upper bounds for theta per dimension; estimated theta values larger than the upperbound values are truncated to the upperbound values
#' @param prior_var_safe_ml if not NULL, MAP estimate with prior variance equal to prior_var_safe_ml is computed instead of ML, if ML estimate fails
#' @return a list containing the index of the next item to be administered given a new response (or "stop_test"), 
#' updated estimate of theta, responses, indeces of administered items, and indeces of available items
#' @export
shadowcat_roqua <- function(new_response, estimate, responses, administered, available, prior, model, alpha, beta, guessing = NULL, eta = NULL, start_items, stop_test, estimator, information_summary, item_selection = "MI", constraints_and_characts = NULL, lower_bound = rep(-3, ncol(alpha)), upper_bound = rep(3, ncol(alpha)), prior_var_safe_ml = NULL) {    
  alpha <- as.matrix(alpha)
  beta <- get_beta(model, beta, eta)
  guessing <- get_guessing(guessing, beta) 
  number_items <- nrow(beta)
  number_dimensions <- ncol(alpha)
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  lp_constraints_and_characts <- constraints_correct_format(stop_test$n, number_items, constraints_and_characts$characteristics, constraints_and_characts$constraints)
  
  result <- function() {
    if (is.null(new_response)) { # first iteration: no responses given yet
      index_new_item <- get_next_item(start_items, item_selection, information_summary, lp_constraints_and_characts$constraints, lp_constraints_and_characts$lp_chars, estimate, model, responses, prior, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound)
      return(list(index_new_item = index_new_item,
                  estimate = estimate,
                  responses = responses,
                  administered = index_new_item,
                  available = available[-which(available %in% index_new_item)]))
    }     
    
    responses <- c(responses, new_response)
    estimate <- update_person_estimate(estimate, responses, administered, available)
    continue_test <- !test_must_stop(length(responses), attr(estimate, 'variance'), stop_test$n, stop_test$target)

    if (continue_test) {
      index_new_item <- get_next_item(start_items, item_selection, information_summary, lp_constraints_and_characts$constraints, lp_constraints_and_characts$lp_chars, estimate, model, responses, prior, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound)
      list(index_new_item = index_new_item,
           estimate = estimate,
           responses = responses,
           administered = c(administered, index_new_item),
           available = available[-which(available %in% index_new_item)])
    }
    else {
      list(index_new_item = "stop_test",
           estimate = estimate,
           responses = responses,
           administered = administered,
           available = available)
    } 
  }
  
  # if inititial items have been administered (so we are in the CAT phase), update person estimate after each newly answered item
  update_person_estimate <- function(estimate, responses, administered, available) { 
    if (length(responses) > start_items$n)
      estimate_latent_trait(estimate, responses, prior, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, prior_var_safe_ml)
    else
      estimate
  }
  
  validate <- function() {
    if (is.null(model))
      return(add_error("model", "is missing"))
    if (is.null(alpha))
      return(add_error("alpha", "is missing"))
    if (is.null(beta))
      return(add_error("beta", "is missing"))
    if (is.null(start_items))
      return(add_error("start_items", "is missing"))
    if (is.null(stop_test))
      return(add_error("stop_test", "is missing"))
    if (is.null(estimator))
      return(add_error("estimator", "is missing"))
    if (is.null(information_summary))
      return(add_error("information_summary", "is missing"))
    if (model != "GPCM" && is.null(beta))
      add_error("beta", "is missing")
    if (model == "GPCM" && is.null(beta) && is.null(eta))
      add_error("beta_and_eta", "are both missing; define at least one of them")
    if (model == "GPCM" && !is.null(beta) && !is.null(eta) && !all.equal(row_cumsum(eta), as.matrix(beta)))
      add_error("beta_and_eta", "objects do not match.")
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
    as.matrix(beta)
}

#' get guessing matrix
#' 
#' @param guessing vector containing guessing parameters; may be NULL in case of zero guessing parameters
#' @param beta matrix containing beta parameters
#' @return matrix containing guessing parameters
#' @export
get_guessing <- function(guessing, beta) {
  if (is.null(guessing))
    matrix(0, nrow(as.matrix(beta)), 1)
  else
    as.matrix(guessing)
}

