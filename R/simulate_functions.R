#' Simulate alpha and beta matrices
#' 
#' Simulate quick and simple itembanks.
#' 
#' @param model Model for which the item bank should be simulated. One of \code{"3PLM"}, \code{"GPCM"}, \code{"SM"} or \code{"GRM"}, 
#' for the three-parameter logistic, generalized partial credit, sequential or graded response model respectively. 
#' If model is 3PLM, number of item steps is forced to 1. 
#' @param number_items Number of items.
#' @param number_dimensions Number of dimensions.
#' @param number_itemsteps Number of item steps (number of categories minus 1); forced to 1 if model is 3PLM.
#' @param items_load_one_dimension If \code{TRUE}, force items to load on one dimension each.
#' @param varying_number_item_steps If \code{TRUE}, some item steps are set to NA; in this case \code{number_itemsteps}
#' is the maximum number of itemsteps.
#' @param alpha_bounds Vector containing lower and upper bound, respectively, of the uniform distribution from which the alpha values are drawn.
#' @return Named list containing simulated alpha and beta matrix.
#' @examples simulate_testbank(model = "GPCM", number_items = 50, number_dimensions = 2, number_itemsteps = 3)
#' simulate_testbank(model = "GPCM", number_items = 50, number_dimensions = 3, number_itemsteps = 4, items_load_one_dimension = TRUE, varying_number_item_steps = TRUE)
#' @importFrom stringr str_c
#' @importFrom stats runif rnorm
#' @export
simulate_testbank <- function(model, number_items = 50, number_dimensions = 1, number_itemsteps = 4, items_load_one_dimension = FALSE, varying_number_item_steps = FALSE, alpha_bounds = c(.3, 1.5)) {
  result <- function() {
    number_itemsteps <- get_number_itemsteps()
    alpha <- get_alpha(number_itemsteps)
    rownames(alpha) <- str_c("item", 1:number_items)
    beta <- get_beta(number_itemsteps)
    rownames(beta) <- str_c("item", 1:number_items)
    list(alpha = alpha, beta = beta)
  }
  
  get_number_itemsteps <- function() {
    if (model == "3PLM") 
      1
    else
      number_itemsteps
  }
  
  get_alpha <- function(number_itemsteps) {
    alpha <- matrix(runif(number_items * number_dimensions, alpha_bounds[1], alpha_bounds[2]), nrow = number_items, ncol = number_dimensions)
    if (number_dimensions > 1 && items_load_one_dimension)
      make_alpha_load_one_dimension(alpha = alpha)
    else
      alpha
  }
  
  get_beta <- function(number_itemsteps) {
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
        induce_varying_number_item_steps(beta_multiple_itemsteps = beta_multiple_itemsteps, number_itemsteps = number_itemsteps)
      else
        beta_multiple_itemsteps
    }  
  }
  
  make_alpha_load_one_dimension <- function(alpha) {
    set_size <- ceiling(number_items / number_dimensions)
    last_set_size <- number_items - (number_dimensions - 1) * set_size
    for (dimension in 1:(number_dimensions - 1))
      alpha[((dimension - 1) * set_size + 1):(dimension * set_size), (1:number_dimensions)[-dimension]] <- 0
    alpha[((number_dimensions - 1) * set_size + 1):nrow(alpha), (1:number_dimensions)[-number_dimensions]] <- 0
    alpha
  }
  
  induce_varying_number_item_steps <- function(beta_multiple_itemsteps, number_itemsteps) {
    beta_multiple_itemsteps[sample(1:number_items, size = ceiling(number_items / 10)), number_itemsteps] <- NA
    if (number_itemsteps > 2)
      beta_multiple_itemsteps[sample(1:number_items, size = ceiling(number_items / 10)), (number_itemsteps - 1):number_itemsteps] <- NA
    beta_multiple_itemsteps
  }
  
  result()
}


#' Simulate answer
#' 
#' Simulate answer on specified items, given true theta.
#' 
#' @param theta Vector with true theta.
#' @param model String, one of \code{"3PLM"}, \code{"GPCM"}, \code{"SM"} or \code{"GRM"}, for the three-parameter logistic, generalized partial credit, sequential or graded response model respectively.
#' @param alpha Matrix of alpha parameters. See \code{\link{shadowcat}} for details.
#' @param beta Matrix of beta parameters, one column per item step, one row per item. See \code{\link{shadowcat}} for details.
#' @param guessing One column matrix of guessing parameters per item. Row names should contain the item keys. Optionally used in 3PLM model, ignored for all others.
#' @param item_keys Character vector of item keys for which answers should be simulated.
#' @return Vector with responses.
#' @examples 
#' alpha_beta <- simulate_testbank(model = "3PLM", number_items = 50, number_dimensions = 1, number_itemsteps = 1)
#' guessing <- matrix(rep(.5, 50), dimnames = list(rownames(alpha_beta$alpha), NULL))
#' 
#' # Without guessing parameter
#' simulate_answer(theta = .3, model = "3PLM", alpha = alpha_beta$alpha, beta = alpha_beta$beta, guessing = NULL, item_keys = "item3")
#' 
#' # With guessing parameter
#' simulate_answer(theta = .3, model = "3PLM", alpha = alpha_beta$alpha, beta = alpha_beta$beta, guessing = guessing, item_keys = "item3")
#' 
#' # Simulate answers for more than one item
#' simulate_answer(theta = .3, model = "3PLM", alpha = alpha_beta$alpha, beta = alpha_beta$beta, guessing = NULL, item_keys = c("item5", "item2", "item8", "item1", "item18"))
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
    apply(random_numbers > cumulative_probabilities, 1, sum, na.rm = TRUE)
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
#' Simulate several iterations of getting the key of the new item to administer and the updated theta estimate, and getting an answer to the new item.
#' 
#' @param true_theta True theta value or vector.
#' @param prior_form String indicating the form of the prior; one of \code{"normal"} or \code{"uniform"}. Not required if estimator is maximum likelihood.
#' @param prior_parameters List containing mu and Sigma of the normal prior: \code{list(mu = ..., Sigma = ...)}, or 
#' the upper and lower bound of the uniform prior: \code{list(lower_bound = ..., upper_bound = ...)}. Not required if estimator is maximum likelihood.
#' See \code{\link{shadowcat}} for details.
#' @param model One of \code{"3PLM"}, \code{"GPCM"}, \code{"SM"} or \code{"GRM"}, for the three-parameter logistic, generalized partial credit, sequential or graded response model, respectively.
#' @param alpha Matrix of alpha parameters. See \code{\link{shadowcat}} for details.
#' @param beta Matrix of beta parameters. See \code{\link{shadowcat}} for details.
#' @param guessing Matrix with one column of guessing parameters per item. See \code{\link{shadowcat}} for details.
#' @param eta Matrix of location parameters, optionally used in GPCM model, ignored for all others. See \code{\link{shadowcat}} for details.
#' @param start_items List indicating the items that should be shown to the respondent before the theta estimate will be updated
#' for the first time. See \code{\link{shadowcat}} for details.
#' @param stop_test List indicating rules for when to terminate the test. See \code{\link{shadowcat}} for details.
#' @param estimator Type of estimator to be used, one of \code{"maximum_likelihood"}, \code{"maximum_aposteriori"}, or \code{"expected_aposteriori"}.
#' @param information_summary How to summarize Fisher information, used for selection of item with maximum information. One of
#' \code{"determinant"}, \code{"posterior_determinant"}, \code{"trace"}, \code{"posterior_trace"}, or \code{"posterior_expected_kullback_leibler"}.
#' @param constraints_and_characts List with constraints and characteristics. See \code{\link{shadowcat}} for details.
#' @param lower_bound Vector with lower bounds for theta per dimension. See \code{\link{shadowcat}} for details.
#' @param upper_bound Vector with upper bounds for theta per dimension. See \code{\link{shadowcat}} for details.
#' @param safe_eap Only relevant if estimator is expected aposteriori. 
#' Set to \code{TRUE} if estimator should switch to maximum aposteriori if the integration algorithm results in an error.
#' An error may occur if the prior is uniform, estimator is expected aposteriori, and the bounds of the prior do not exceed the true theta value, or are too close to it.
#' @param initital_estimate Vector containing the initial theta estimate, before any items have been administered.
#' @param initial_variance Matrix containing the initial covariance matrix, before any items have been administered.
#' @param eap_estimation_procedure String indicating the estimation procedure if estimator is expected aposteriori and prior form is normal. One of \code{"riemannsum"} 
#' for integration via Riemannsum or \code{"gauss_hermite_quad"} for integration via Gaussian Hermite Quadrature. 
#' If prior form is uniform, estimation procedure should always be \code{"riemannsum"}.
#' @return List as returned by \code{\link{shadowcat}} after test is terminated, with \code{variance} element turned into matrix.
#' @examples 
#' # One dimension
#' alpha_beta_one_dim <- simulate_testbank(model = "GPCM", number_items = 50, number_dimensions = 1, number_itemsteps = 3)
#' test_shadowcat(true_theta = 2, prior_form = "normal", prior_parameters = list(mu = 0, Sigma = diag(1)), model = "SM", alpha = alpha_beta_one_dim$alpha, beta = alpha_beta_one_dim$beta, guessing = NULL, start_items = list(type = 'random', n = 3), stop_test = list(max_n = 20, target = 0.1), estimator = "maximum_aposteriori", information_summary = "posterior_determinant")
#' 
#' # Three dimensions
#' alpha_beta_three_dim <- simulate_testbank(model = "GPCM", number_items = 100, number_dimensions = 3, number_itemsteps = 3)
#' test_shadowcat(true_theta = c(0, 1, -.5), prior_form = "normal", prior_parameters = list(mu = c(0, 0, 0), Sigma = diag(3)), model = "SM", alpha = alpha_beta_three_dim$alpha, beta = alpha_beta_three_dim$beta, guessing = NULL, start_items = list(type = 'random', n = 3), stop_test = list(max_n = 60, target = c(.1, .1, .1)), estimator = "maximum_aposteriori", information_summary = "posterior_determinant")
#' @export
test_shadowcat <- function(true_theta, prior_form, prior_parameters, model, alpha, beta, guessing, eta = NULL, start_items, stop_test, estimator, information_summary, constraints_and_characts = NULL, lower_bound = NULL, upper_bound = NULL, safe_eap = FALSE, initital_estimate = rep(0, ncol(alpha)), initial_variance = diag(ncol(alpha)) * 25, eap_estimation_procedure = "riemannsum") {
  answers <- NULL
  next_item_and_theta_estimate <- shadowcat(answers = answers, estimate = initital_estimate, variance = as.vector(initial_variance), model = model, alpha = alpha, beta = beta, start_items = start_items, stop_test = stop_test, estimator = estimator, information_summary = information_summary, prior_form = prior_form, prior_parameters = prior_parameters, guessing = guessing, eta = eta, constraints_and_characts = constraints_and_characts, lower_bound = lower_bound, upper_bound = upper_bound, safe_eap = safe_eap, eap_estimation_procedure = eap_estimation_procedure)
  
  while (next_item_and_theta_estimate$continue_test) {
    new_answer <- simulate_answer(theta = true_theta, model = model, alpha = alpha, beta = beta, guessing = guessing, item_keys = next_item_and_theta_estimate$key_new_item)
    next_item_and_theta_estimate$answers[[next_item_and_theta_estimate$key_new_item]] <- new_answer
    next_item_and_theta_estimate <- shadowcat(answers = as.list(next_item_and_theta_estimate$answers), estimate = next_item_and_theta_estimate$estimate, variance = next_item_and_theta_estimate$variance, model = model, alpha = alpha, beta = beta, start_items = start_items, stop_test = stop_test, estimator = estimator, information_summary = information_summary, prior_form = prior_form, prior_parameters = prior_parameters, guessing = guessing, eta = eta, constraints_and_characts = constraints_and_characts, lower_bound = lower_bound, upper_bound = upper_bound, safe_eap = safe_eap, eap_estimation_procedure = eap_estimation_procedure)  
  }
  
  next_item_and_theta_estimate$variance <- matrix(next_item_and_theta_estimate$variance, ncol = ncol(alpha))
  next_item_and_theta_estimate
}


