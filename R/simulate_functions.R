#' Simulate alpha and beta matrices
#' 
#' Simulate quick and simple itembanks. Only for testing purposess
#' 
#' @param model String, one of '3PLM', 'GPCM', 'SM' or 'GRM', for the three-parameter logistic, generalized partial credit,
#'  sequential or graded response model respectively.
#' @param number_items Number of items. May be increased if number_items / number_dimensions is not an integer.
#' @param number_dimensions Number of dimensions.
#' @param number_itemsteps Number of item steps (number of categories minus 1); forced to 1 if model is 3PLM.
#' @param items_load_one_dimension If TRUE, force items to load on one dimension each.
#' @param varying_number_item_steps If TRUE, some item steps are set to NA; in this case number_itemsteps 
#' is the maximum number of itemsteps.
#' @return List containing simulated alpha and beta matrix
#' @importFrom stringr str_c
#' @importFrom stats runif rnorm
#' @export
simulate_testbank <- function(model, number_items = 50, number_dimensions = 1, number_itemsteps = 4, items_load_one_dimension = FALSE, varying_number_item_steps = FALSE){
  result <- function() {
    number_items <- get_number_items()
    number_itemsteps <- get_number_itemsteps()
    alpha <- get_alpha(number_items, number_itemsteps)
    rownames(alpha) <- str_c("item", 1:number_items)
    beta <- get_beta(number_items, number_itemsteps)
    rownames(beta) <- str_c("item", 1:number_items)
    list(alpha = alpha, beta = beta)
  }
  
  get_number_itemsteps <- function() {
    if (model == "3PLM") 
      1
    else
      number_itemsteps
  }
  
  get_number_items <- function() {
    # make sure the number of items is divisible by the number of dimensions
    ceiling(number_items / number_dimensions) * number_dimensions
  }
  
  get_alpha <- function(number_items, number_itemsteps) {
    alpha <- matrix(runif(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
    if (number_dimensions > 1 && items_load_one_dimension) {
      set <- number_items / number_dimensions
      for (dimension in 1:number_dimensions)
        alpha[((dimension - 1) * set + 1):(dimension * set), (1:number_dimensions)[-dimension]] <- 0
      alpha
    }
    else {
      alpha
    }
  }
  
  get_beta <- function(number_items, number_itemsteps) {
    beta_one_itemstep <- matrix(rnorm(number_items), nrow = number_items, ncol = 1)
    if (number_itemsteps == 1) {
      beta_one_itemstep
    }
    else {
      # spread polytomous items cats -2 to +2.
      spread <- seq(-2, 2, length.out = number_itemsteps)
      beta_multiple_itemsteps <- if (model == "GPCM") 
                                   row_cumsum(t(apply(beta_one_itemstep, 1, function(x) x + spread))) 
                                  else
                                    beta_multiple_itemsteps <- t(apply(beta_one_itemstep, 1, function(x) x + spread)) 
      if (varying_number_item_steps)
        induce_varying_number_item_steps(beta_multiple_itemsteps = beta_multiple_itemsteps, number_items = number_items, number_itemsteps = number_itemsteps)
      else
        beta_multiple_itemsteps
    }  
  }
  
  induce_varying_number_item_steps <- function(beta_multiple_itemsteps, number_items, number_itemsteps) {
    beta_multiple_itemsteps[sample(1:number_items, size = ceiling(number_items / 10)), number_itemsteps] <- NA
    if (number_itemsteps > 2)
      beta_multiple_itemsteps[sample(1:number_items, size = ceiling(number_items / 10)), (number_itemsteps - 1):number_itemsteps] <- NA
    beta_multiple_itemsteps
  }
  
  result()
}


#' Simulates answers
#' 
#' Simulates answers on specified items, given true theta. Only for testing purposes
#' 
#' @param theta Vector with true theta.
#' @param model One of '3PLM', 'GPCM', 'SM' or 'GRM', for the three-parameter logistic, generalized partial credit, sequential or graded response model respectively.
#' @param alpha Matrix of alpha parameters. See \code{shadowcat} for details.
#' @param beta Matrix of beta parameters, one column per item step, one row per item. See \code{shadowcat} for details.
#' @param guessing One column matrix of guessing parameters per item. Row names should contain the item keys. Optionally used in 3PLM model, ignored for all others.
#' @param item_keys Character vector of item keys for which answers should be simulated.
#' @return Vector with responses
#' @examples 
#' item_bank <- simulate_testbank("GPCM")
#' simulate_answer(theta = .3, model = "GPCM", alpha = item_bank$alpha, beta = item_bank$beta, guessing = NULL, item_keys = "item3")
#' @importFrom stats runif
#' @export
simulate_answer <- function(theta, model, alpha, beta, guessing, item_keys) {
  indices <- match(item_keys, rownames(alpha))
  result <- function() {
    guessing <- get_guessing()
    # probabilities, generated with true theta.
    probabilities <- get_probs_and_likelihoods_per_item(theta = theta, model = model, alpha = get_subset(alpha, indices), beta = get_subset(beta, indices), guessing = get_subset(guessing, indices), with_likelihoods = FALSE)$P
    cumulative_probabilities <- row_cumsum(probabilities) 
    random_numbers <- runif(length(indices))
    
    # answer is the number of categories that have a cumulative probability smaller than random_numbers
    apply(random_numbers > cumulative_probabilities, 1, sum, na.rm=TRUE)
  }
  
  get_guessing <- function() {
    if (is.null(guessing))
      matrix(0, nrow = nrow(as.matrix(alpha)), ncol = 1, dimnames = list(rownames(alpha), NULL))
    else
      guessing
  }
  
  result()
}


#' Simulate a testing routine with shadowcat
#' 
#' @param true_theta True theta value or vector
#' @param prior_form String indicating the form of the prior; one of "normal" or "uniform"
#' @param prior_parameters List containing mu and Sigma of the normal prior: list(mu = ..., Sigma = ...), or 
#' the upper and lower bound of the uniform prior: list(lower_bound = ..., upper_bound = ...). Sigma should always
#' be in matrix form. The length of lower_bound and upper_bound should be equal to the number of dimensions
#' @param model String, one of '3PLM', 'GPCM', 'SM' or 'GRM', for the three-parameter logistic, generalized partial credit, sequential or graded response model respectively.
#' @param alpha Matrix of alpha parameters. See \code{shadowcat} for details.
#' @param beta Matrix of beta parameters. See \code{shadowcat} for details.
#' @param guessing Matrix with one column of guessing parameters per item. See \code{shadowcat} for details.
#' @param eta Matrix of location parameters, optionally used in GPCM model, ignored for all others. See \code{shadowcat} for details.
#' @param start_items List indicating the items that are shown to the respondent before adaptive proces starts. One of
#' list(type = 'random', n)
#' list(type = 'fixed', indices, n)
#' list(type = 'random_by_dimension', n_by_dimension, n)
#' See \code{shadowcat} for details.
#' @param stop_test List indicating rules for when to stop providing new items to respondent. Should be of the form
#' list(target = ..., max_n = ..., min_n = ..., cutoffs = ...). See \code{shadowcat} for details.
#' @param estimator Type of estimator to be used, one of "maximum_aposteriori", "maximum_likelihood", or "expected_aposteriori".
#' @param information_summary How to summarize Fisher information, used for selection of item with maximum information. 
#' One of "determinant", "posterior_determinant", "trace", "posterior_trace", or "posterior_expected_kullback_leibler". See \code{shadowcat} for details.
#' @param constraints_and_characts list with constraints and characteristics. See \code{shadowcat} for details.
#' @param lower_bound Vector with lower bounds for theta per dimension; estimated theta values smaller than the lowerbound values are truncated to the lowerbound values.
#' Can only be defined when estimator is maximum_likelihood. Setting bounds with maximum likelihood estimation is equivalent to
#' using maximum aposteriori estimation with a uniform prior. 
#' @param upper_bound Vector with upper bounds for theta per dimension; estimated theta values larger than the upperbound values are truncated to the upperbound values
#' Can only be defined when estimator is maximum_likelihood. Setting bounds with maximum likelihood estimation is equivalent to
#' using maximum aposteriori estimation with a uniform prior.
#' @param safe_eap Only relevant if estimator is expected_aposteriori. 
#' TRUE if estimator should switch to maximum aposteriori if the integration algorithm results in an error.
#' An error may occur if the prior is uniform, estimator is expected aposteriori, and the bounds do not exceed the true theta value, or are too close to it.
#' @param initital_estimate Vector containing the initial theta estimates (starting values)
#' @param initial_variance Matrix containing the initial covariance matrix (staring values)
#' @param eap_estimation_procedure String indicating the estimation procedure if estimator is expected aposteriori. One of "riemannsum" for integration via Riemannsum or
#' "gauss_hermite_quad" for integration via Gaussian Hermite Quadrature. 
#' @return List as returned by shadowcat after the last answer, with variance element turned into matrix
#' @export
test_shadowcat <- function(true_theta, prior_form, prior_parameters, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, constraints_and_characts = NULL, lower_bound = NULL, upper_bound = NULL, safe_eap = FALSE, initital_estimate = rep(0, ncol(alpha)), initial_variance = diag(ncol(alpha)) * 25, eap_estimation_procedure = "riemannsum") {
  answers <- NULL
  next_item_and_test_outcome <- shadowcat(answers = answers, estimate = initital_estimate, variance = as.vector(initial_variance), model = model, alpha = alpha, beta = beta, start_items = start_items, stop_test = stop_test, estimator = estimator, information_summary = information_summary, prior_form = prior_form, prior_parameters = prior_parameters, guessing = guessing, eta = eta, constraints_and_characts = constraints_and_characts, lower_bound = lower_bound, upper_bound = upper_bound, safe_eap = safe_eap, eap_estimation_procedure = eap_estimation_procedure)
  
  while (next_item_and_test_outcome$continue_test) {
    new_answer <- simulate_answer(theta = true_theta, model = model, alpha = alpha, beta = beta, guessing = guessing, item_keys = next_item_and_test_outcome$key_new_item)
    next_item_and_test_outcome$answers[[next_item_and_test_outcome$key_new_item]] <- new_answer
    next_item_and_test_outcome <- shadowcat(answers = as.list(next_item_and_test_outcome$answers), estimate = next_item_and_test_outcome$estimate, variance = next_item_and_test_outcome$variance, model = model, alpha = alpha, beta = beta, start_items = start_items, stop_test = stop_test, estimator = estimator, information_summary = information_summary, prior_form = prior_form, prior_parameters = prior_parameters, guessing = guessing, eta = eta, constraints_and_characts = constraints_and_characts, lower_bound = lower_bound, upper_bound = upper_bound, safe_eap = safe_eap, eap_estimation_procedure = eap_estimation_procedure)  
  }
  
  attr(next_item_and_test_outcome$estimate, "variance") <- matrix(next_item_and_test_outcome$variance, ncol = ncol(alpha))
  next_item_and_test_outcome
}


