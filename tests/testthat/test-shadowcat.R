# only for whithin R:
'
library(testthat)
library(pbapply)
library(stringr)
library(matrixcalc)
'

make_random_seed_exist <- rnorm(1)

#' simulate a testing routine with shadowcat
#' 
#' @param true_theta true theta value or vector
#' @param prior_form String indicating the form of the prior; one of "normal" or "uniform"
#' @param prior_parameters List containing mu and Sigma of the normal prior: list(mu = ..., Sigma = ...), or 
#' the upper and lower bound of the uniform prior: list(lower_bound = ..., upper_bound = ...). Sigma should always
#' be in matrix form. The length of lower_bound and upper_bound should be equal to the number of dimensions
#' @param model String, one of '3PLM', 'GPCM', 'SM' or 'GRM', for the three-parameter logistic, generalized partial credit, sequential or graded response model respectively.
#' @param alpha Matrix of alpha parameters, one column per dimension, one row per item. Row names should contain the item keys. Note that so called within-dimensional models still use an alpha matrix, they simply 
#' have only one non-zero loading per item.
#' @param beta Matrix of beta parameters, one column per item step, one row per item. Row names should contain the item keys. Note that ShadowCAT expects answer categories to be sequential,
#' and without gaps. That is, the weight parameter in the GPCM model is assumed to be sequential, and equal to the position of the 'location' of the beta parameter in the Beta matrix.
#' The matrix will have a number of columns equal to the largest number of answer categories, items with fewer answer categories should be 
#' right-padded with \code{NA}. \code{NA} values between answer categories are not allowed, and will lead to errors.
#' Beta matrix can be set to NULL if model is GPCM and eta is defined
#' @param guessing matrix with one column of guessing parameters per item. Row names should contain the item keys. Optionally used in 3PLM model, ignored for all others.
#' @param eta Matrix of location parameters, optionally used in GPCM model, ignored for all others.
#' @param start_items items that are shown to the patient before adaptive proces starts; one of
#' list(type = 'random', n)
#' list(type = 'fixed', indices, n)
#' list(type = 'random_by_dimension', n_by_dimension, n)
#' where n = total number of initial items, indices = vector of initial item indices, 
#' n_by_dimension = scalar of number of initial items per dimension, or vector with number of initial items for each dimension
#' If n is 0, only n needs to be defined
#' @param stop_test rule for when to stop providing new items to patient; should be a list of the form
#' list(target = ..., max_n = ..., min_n = ..., cutoffs = ...), 
#' where max_n = test length at which testing should stop (even if target has not been reached yet in case of variance stopping rule), 
#' target = vector of maximum acceptable variances per dimension; if target = NULL, only max_n is taken into account,
#' min_n = minimum test length; NULL means no mimimum test length,
#' cutoffs = matrix containing cut off values per dimension (columns) and test replication (rows). First row contains cut off values for when no items have been
#' administered yet, second row for when one item has been administered, etc. If estimate + 3SE < cutoff for each dimension at certain replication, test stops; 
#' NULL means no cut off values
#' @param estimator type of estimator to be used, one of "maximum_aposteriori", "maximum_likelihood", or "expected_aposteriori"
#' @param information_summary called "objective" by Kroeze; how to summarize information; one of
#' "determinant": compute determinant(info_sofar_QxQ + info_QxQ_k) for each yet available item k
#' "posterior_determinant": compute determinant(info_sofar_QxQ_plus_prior + info_QxQ_k) for each yet available item k
#' "trace": compute trace((info_sofar_QxQ + info_QxQ_k) for each yet available item k
#' "posterior_trace": compute trace(info_sofar_QxQ_plus_prior + info_QxQ_k) for each yet available item k
#' "posterior_expected_kullback_leibler" = compute Posterior expected Kullback-Leibler Information
#' @param constraints_and_characts list with constraints and characteristics: constraints_and_characts = list(constraints = ..., characteristics = ...)
#' constraints should be specified as a list of constraints, each constraint is a list with three named values;
#' name: the column name of the characteristic this constraint applies to. For categorical characteristics the level should be specified as name/value.
#' op: the logical operator to be used. Valid options are "<", "=", ">" and "><".
#' target: the target value, numeric. If the operator is "><", this should be a length two vector in between which the target should fall.
#' characteristics should be a data.frame with characteristics, one row per item, one column per characteristic.
#' See constraints_lp_format() for details
#' @param lower_bound Vector with lower bounds for theta per dimension; estimated theta values smaller than the lowerbound values are truncated to the lowerbound values.
#' Can only be defined when estimator is maximum_likelihood. Setting bounds with maximum likelihood estimation is equivalent to
#' using maximum aposteriori estimation with a uniform prior. 
#' @param upper_bound Vector with upper bounds for theta per dimension; estimated theta values larger than the upperbound values are truncated to the upperbound values
#' Can only be defined when estimator is maximum_likelihood. Setting bounds with maximum likelihood estimation is equivalent to
#' using maximum aposteriori estimation with a uniform prior.
#' @param safe_eap Only relevant if estimator is espected_aposteriori. 
#' TRUE if estimator should switch to maximum aposteriori if the integration algorithm results in an error.
#' @param initital_estimate vector containing the initial theta estimates (starting values)
#' @param initial_variance matrix containing the initial covariance matrix (staring values)
#' @param eap_estimation_procedure String indicating the estimation procedure if estimator is expected aposteriori. One of "riemannsum" for integration via Riemannsum or
#' "gauss_hermite_quad" for integration via Gaussian Hermite Quadrature. 
#' @return
test_shadowcat <- function(true_theta, prior_form, prior_parameters, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, constraints_and_characts = NULL, lower_bound = NULL, upper_bound = NULL, safe_eap = FALSE, initital_estimate = rep(0, ncol(alpha)), initial_variance = diag(ncol(alpha)) * 25, eap_estimation_procedure = "riemannsum") {
  item_keys <- rownames(alpha)
  answers <- NULL
  next_item_and_test_outcome <- shadowcat(answers, estimate = initital_estimate, variance = as.vector(initial_variance), model, alpha, beta, start_items, stop_test, estimator, information_summary, prior_form, prior_parameters, guessing, eta, constraints_and_characts, lower_bound, upper_bound, safe_eap, eap_estimation_procedure)

  while (next_item_and_test_outcome$continue_test) {
    new_answer <- simulate_answer(true_theta, model, ncol(alpha), estimator, alpha, beta, guessing, ncol(beta), match(next_item_and_test_outcome$key_new_item, item_keys))
    next_item_and_test_outcome$answers[[next_item_and_test_outcome$key_new_item]] <- new_answer
    next_item_and_test_outcome$answers <- as.list(next_item_and_test_outcome$answers)
    next_item_and_test_outcome <- shadowcat(next_item_and_test_outcome$answers, next_item_and_test_outcome$estimate, next_item_and_test_outcome$variance, model, alpha, beta, start_items, stop_test, estimator, information_summary, prior_form, prior_parameters, guessing, eta, constraints_and_characts, lower_bound, upper_bound, safe_eap, eap_estimation_procedure)  
  }
  
  attr(next_item_and_test_outcome$estimate, "variance") <- matrix(next_item_and_test_outcome$variance, ncol = ncol(alpha))
  next_item_and_test_outcome
}

#' get conditions for simulation
#' 
#' @param true_theta_vec vector containing true theta values. If number dimensions is 1, simulations are performed for each true theta value in true_theta_vec;
#' if number dimensions is larger than 1, true_theta_vec is interpreted as containing the true thetas for each dimension
#' @param number_items_vec vector containing conditions for number of test bank items to be simulated. Simulations are performed for each value in number_items_vec.
#' If constraints are added, number_items_vec can only have length 1
#' @param number_answer_categories_vec vector containing conditions for the number of answer categories to be simulated. 
#' Simulations are performed for each value in number_answer_categories_vec
#' @param model_vec vector containing the conditions for the model to be used. Simulations are performed for each model in model_vec. Model options are
#' '3PLM', 'GPCM', 'SM' and 'GRM'
#' @param estimator_vec vector containing the conditions for the estimator to be used. Simulations are performed for each estimator in estimator_vec.
#' Options are "maximum_likelihood", "maximum_aposteriori", and "expected_aposteriori"
#' @param information_summary_vec vector containing the conditions for the information_summary to be used. Simulations are performed for each model summary in estimator_vec,
#' Options are "determinant", "posterior_determinant", "trace", "posterior_trace", and "posterior_expected_kullback_leibler" 
#' @param replications_per_unique_condition number of replications to be performed within each unique condition
#' @param number_dimensions the number of dimensions of the model (either 1 or the length of true_theta_vec)
get_conditions <- function(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, replications_per_unique_condition, number_dimensions) {
  if (number_dimensions == 1) {
    conditions <- expand.grid(1:replications_per_unique_condition, true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec)
    colnames(conditions) <- c("replication", "true_theta", "number_items", "number_answer_categories", "model", "estimator", "information_summary")
  }
  else {
    conditions <- expand.grid(1:replications_per_unique_condition, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec)
    colnames(conditions) <- c("replication", "number_items", "number_answer_categories", "model", "estimator", "information_summary") 
  }
  conditions
}

#' simulate testing routines with shadowcat, for several conditions
#' 
#' @param true_theta_vec numeric vector containing true theta values. If number dimensions is 1, simulations are performed for each true theta value in true_theta_vec;
#' if number dimensions is larger than 1, true_theta_vec is interpreted as containing the true thetas for each dimension
#' @param number_items_vec numeric vector containing conditions for number of test bank items to be simulated. Simulations are performed for each value in number_items_vec.
#' If constraints_and_characts is defined, number_items_vec can only have length 1
#' @param number_answer_categories_vec numeric vector containing conditions for the number of answer categories to be simulated. 
#' Simulations are performed for each value in number_answer_categories_vec
#' @param model_vec character vector containing the conditions for the model to be used. Simulations are performed for each model in model_vec. Model options are
#' '3PLM', 'GPCM', 'SM' and 'GRM'
#' @param estimator_vec character vector containing the conditions for the estimator to be used. Simulations are performed for each estimator in estimator_vec.
#' Options are "maximum_likelihood", "maximum_aposteriori", and "expected_aposteriori"
#' @param information_summary_vec character vector containing the conditions for the information_summary to be used. Simulations are performed for each model summary in estimator_vec,
#' Options are "determinant", "posterior_determinant", "trace", "posterior_trace", and "posterior_expected_kullback_leibler"
#' @param start_items list indicating the items that are shown to the patient before adaptive proces starts; one of
#' list(type = 'random', n)
#' list(type = 'fixed', indices, n)
#' list(type = 'random_by_dimension', n_by_dimension, n)
#' where n = total number of initial items, indices = vector of initial item indices, 
#' n_by_dimension = scalar of number of initial items per dimension, or vector with number of initial items for each dimension
#' If n is 0, only n needs to be defined
#' @param variance_target value equal to the variance of theta at which testing should stop
#' @param replications_per_unique_condition value equal to the number of replications to be performed within each unique condition
#' @param number_dimensions value equal to the the number of dimensions of the model (either 1 or the length of true_theta_vec)
#' @param constraints_and_characts list with constraints and characteristics: constraints_and_characts = list(constraints = ..., characteristics = ...)
#' constraints should be specified as a list of constraints, each constraint is a list with three named values;
#' name: the column name of the characteristic this constraint applies to. For categorical characteristics the level should be specified as name/value.
#' op: the logical operator to be used. Valid options are "<", "=", ">" and "><".
#' target: the target value, numeric. If the operator is "><", this should be a length two vector in between which the target should fall.
#' characteristics should be a data.frame with characteristics, one row per item, one column per characteristic.
#' See constraints_lp_format() for details
#' @param guessing matrix with one column of guessing parameters per item. Row names should contain the item keys. Optionally used in 3PLM model, ignored for all others.
#' @param items_load_one_dimension if TRUE, items are simulated which load on one dimension. If FALSE, items are simulated which load on all dimensions
#' @param prior_form String indicating the form of the prior; one of "normal" or "uniform"
#' @param prior_parameters List containing mu and Sigma of the normal prior: list(mu = ..., Sigma = ...), or 
#' the upper and lower bound of the uniform prior: list(lower_bound = ..., upper_bound = ...). Sigma should always
#' be in matrix form. The length of lower_bound and upper_bound should be equal to the number of dimensions
#' @param lower_bound Vector with lower bounds for theta per dimension; estimated theta values smaller than the lowerbound values are truncated to the lowerbound values.
#' Is only used when estimator is maximum_likelihood. Setting bounds with maximum likelihood estimation is equivalent to
#' using maximum aposteriori estimation with a uniform prior. 
#' @param upper_bound Vector with upper bounds for theta per dimension; estimated theta values larger than the upperbound values are truncated to the upperbound values
#' Is only used when estimator is maximum_likelihood. Setting bounds with maximum likelihood estimation is equivalent to
#' using maximum aposteriori estimation with a uniform prior.
#' @param safe_eap Only relevant if estimator is espected_aposteriori. 
#' TRUE if estimator should switch to maximum aposteriori if the integration algorithm results in an error.
#' @param return_administered_item_indices if TRUE, indices of administered items are added to the output
#' @param min_n value equal to the minimum number of items to administer
#' @param max_n value equal to the maximum number of items to administer (test stops at this number, even if variance target has not been reached). NULL means max_n is equal to number of items in test bank
#' @param varying_number_item_steps if TRUE, the simulated number of item steps differs over items. In that case, number_answer_categories_vec (number of itemsteps + 1)
#' is considered the maxixmum number of categories
#' @param eap_estimation_procedure String indicating the estimation procedure if estimator is expected aposteriori. One of "riemannsum" for integration via Riemannsum or
#' "gauss_hermite_quad" for integration via Gaussian Hermite Quadrature. 
#' @return matrix with in each row (= one condition): named vector containing estimated theta, variance of the estimate, and if return_administered_item_indices is TRUE, the indices of the administered items
run_simulation <- function(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, start_items, variance_target, replications_per_unique_condition, number_dimensions, constraints_and_characts = NULL, guessing = NULL, items_load_one_dimension = TRUE, prior_form = "normal", prior_parameters = list(mu = rep(0, length(number_dimensions)), Sigma = diag(number_dimensions) * 20), lower_bound = NULL, upper_bound = NULL, safe_eap = FALSE, return_administered_item_indices = FALSE, min_n = NULL, max_n = NULL, varying_number_item_steps = FALSE, eap_estimation_procedure = "riemannsum") {                   
  if (number_dimensions > 1 && number_dimensions != length(true_theta_vec))
    stop("number_dimensions is larger than 1 but not equal to the length of true_theta_vec")
  conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, replications_per_unique_condition, number_dimensions)
  
  pbapply::pbsapply(1:nrow(conditions), 
                    FUN = function(condition) {
                      if (is.null(max_n))
                        max_n <- conditions[condition, "number_items"] 
                      stop_test <- list(target = variance_target, max_n = max_n, min_n = min_n)
                      true_theta <- ( if (number_dimensions == 1) 
                                        conditions[condition, "true_theta"]
                                      else
                                        true_theta_vec )
                      alpha_beta <- simulate_testbank(model = as.character(conditions[condition, "model"]), number_items = conditions[condition, "number_items"], number_dimensions = number_dimensions, number_itemsteps = conditions[condition, "number_answer_categories"] - 1, items_load_one_dimension = items_load_one_dimension, varying_number_item_steps = varying_number_item_steps)
                      estimate_theta <- tryCatch(test_shadowcat(true_theta, prior_form, prior_parameters, as.character(conditions[condition, "model"]), alpha_beta$alpha, alpha_beta$beta, guessing, eta = NULL, start_items, stop_test, as.character(conditions[condition, "estimator"]), as.character(conditions[condition, "information_summary"]), constraints_and_characts, lower_bound, upper_bound, safe_eap, eap_estimation_procedure = eap_estimation_procedure),
                                                 error = function(e) e)
                      
                      if (return_administered_item_indices)
                        c("estimated_theta" = estimate_theta$estimate,
                          "variance_estimate" = attr(estimate_theta$estimate, "variance"),
                          "items_administered" = as.numeric(sapply(names(estimate_theta$answers), substring, 5)))
                      else
                        c("estimated_theta" = estimate_theta$estimate,
                          "variance_estimate" =  attr(estimate_theta$estimate, "variance"))               
                    })
}

context("validate shadowcat single conditions")

test_that("true theta is 2, estimator is maximum_aposteriori, N(0, 5) prior", {
  # define true theta for simulation of answers
  true_theta <- 2
  
  # define item characteristics
  number_items <- 100
  number_dimensions <- 1
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  item_keys <- str_c("item", 1:number_items)
  
  guessing <- matrix(c(rep(.1, number_items / 2), rep(.2, number_items / 2)), ncol = 1, dimnames = list(item_keys, NULL))
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions, dimnames = list(item_keys, NULL))
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1, dimnames = list(item_keys, NULL))
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(type = 'random', n = 3)
  stop_test <- list(max_n = 100)
  estimator <- 'maximum_aposteriori'
  information_summary <- 'posterior_determinant'
  
  prior_form <- "normal"  
  prior_parameters <- list(mu = 0, Sigma = matrix(5))
  
  test_outcome <- with_random_seed(2, test_shadowcat)(true_theta, prior_form, prior_parameters, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary)
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), 1.938)
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3)), .113)
  expect_equal(length(test_outcome$answers), 100)
})

test_that("true theta is 2, estimator is maximum_aposteriori, N(-1, 1) prior", {
  # define true theta for simulation of answers
  true_theta <- 2
  
  # define item characteristics
  number_items <- 100
  number_dimensions <- 1
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  item_keys <- str_c("item", 1:number_items)
  
  guessing <- matrix(c(rep(.1, number_items / 2), rep(.2, number_items / 2)), ncol = 1, dimnames = list(item_keys, NULL))
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions, dimnames = list(item_keys, NULL))
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1, dimnames = list(item_keys, NULL))
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(type = 'random', n = 3)
  stop_test <- list(max_n = 100)
  estimator <- 'maximum_aposteriori'
  information_summary <- 'posterior_determinant'
  
  prior_form <- "normal"  
  prior_parameters <- list(mu = -1, Sigma = matrix(1))
  
  test_outcome <- with_random_seed(2, test_shadowcat)(true_theta, prior_form, prior_parameters, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary)
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), 1.494)
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3)), .089)
  expect_equal(length(test_outcome$answers), 100)
})

test_that("true theta is 2, estimator is maximum_aposteriori, U(-1, 1) prior", {
  # define true theta for simulation of answers
  true_theta <- 2
  
  # define item characteristics
  number_items <- 100
  number_dimensions <- 1
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  item_keys <- str_c("item", 1:number_items)
  
  guessing <- matrix(c(rep(.1, number_items / 2), rep(.2, number_items / 2)), ncol = 1, dimnames = list(item_keys, NULL))
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions, dimnames = list(item_keys, NULL))
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1, dimnames = list(item_keys, NULL))
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(type = 'random', n = 3)
  stop_test <- list(max_n = 100)
  estimator <- 'maximum_aposteriori'
  information_summary <- 'posterior_determinant'
  
  prior_form <- "uniform"  
  prior_parameters <- list(lower_bound = -1, upper_bound = 1)
  
  test_outcome <- with_random_seed(2, test_shadowcat)(true_theta, prior_form, prior_parameters, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary)
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), 1)
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3)), .087)
  expect_equal(length(test_outcome$answers), 100)
})

test_that("true theta is 2, estimator is maximum_likelihood", {
  # define true theta for simulation of answers
  true_theta <- 2
  
  # define item characteristics
  number_items <- 100
  number_dimensions <- 1
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  item_keys <- str_c("item", 1:number_items)
  
  guessing <- matrix(c(rep(.1, number_items / 2), rep(.2, number_items / 2)), dimnames = list(item_keys, NULL))
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions, dimnames = list(item_keys, NULL))
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1, dimnames = list(item_keys, NULL))
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(type = 'random', n = 3)
  stop_test <- list(max_n = 100)
  estimator <- 'maximum_likelihood'
  information_summary <- 'determinant'
  lower_bound <- -3
  upper_bound <- 3
  
  test_outcome <- with_random_seed(2, test_shadowcat)(true_theta, prior_form = NULL, prior_parameters = NULL, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, lower_bound = lower_bound, upper_bound = upper_bound)
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), 2.169)
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3)), .129)
  expect_equal(length(test_outcome$answers), 100) 
})


test_that("true theta is 2, estimator is expected_aposteriori, N(0, 5) prior", {
  # define true theta for simulation of answers
  true_theta <- 2
  
  # define item characteristics
  number_items <- 100
  number_dimensions <- 1
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  item_keys <- str_c("item", 1:number_items)
  
  guessing <- NULL
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions, dimnames = list(item_keys, NULL))
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1, dimnames = list(item_keys, NULL))
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(type = 'random', n = 3)
  stop_test <- list(max_n = 100)
  estimator <- 'expected_aposteriori'
  information_summary <- 'posterior_determinant'
  
  prior_form <- "normal"  
  prior_parameters <- list(mu = 0, Sigma = matrix(5))
  
  test_outcome_gauss_hermite <- with_random_seed(2, test_shadowcat)(true_theta, prior_form = prior_form, prior_parameters = prior_parameters, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, eap_estimation_procedure = "gauss_hermite_quad")
  test_outcome_riemann <- with_random_seed(2, test_shadowcat)(true_theta, prior_form = prior_form, prior_parameters = prior_parameters, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary)
  
  # gauss hermite
  expect_equal(as.vector(round(test_outcome_gauss_hermite$estimate, 3)), 1.833)
  expect_equal(as.vector(round(attr(test_outcome_gauss_hermite$estimate, "variance"), 3)), .087)
  expect_equal(length(test_outcome_gauss_hermite$answers), 100)
  
  # riemann
  expect_equal(as.vector(round(test_outcome_riemann$estimate, 3)), 1.833)
  expect_equal(as.vector(round(attr(test_outcome_riemann$estimate, "variance"), 3)), .087)
  expect_equal(length(test_outcome_riemann$answers), 100)
})

test_that("true theta is 2, estimator is expected_aposteriori, N(-1, 1) prior", {
  # define true theta for simulation of answers
  true_theta <- 2
  
  # define item characteristics
  number_items <- 100
  number_dimensions <- 1
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  item_keys <- str_c("item", 1:number_items)
  
  guessing <- NULL
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions, dimnames = list(item_keys, NULL))
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1, dimnames = list(item_keys, NULL))
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(type = 'random', n = 3)
  stop_test <- list(max_n = 100)
  estimator <- 'expected_aposteriori'
  information_summary <- 'posterior_determinant'
  
  prior_form <- "normal"  
  prior_parameters <- list(mu = -1, Sigma = matrix(1))
  
  test_outcome_gauss_hermite <- with_random_seed(2, test_shadowcat)(true_theta, prior_form = prior_form, prior_parameters = prior_parameters, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, eap_estimation_procedure = "gauss_hermite_quad")
  test_outcome_riemann <- with_random_seed(2, test_shadowcat)(true_theta, prior_form = prior_form, prior_parameters = prior_parameters, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary)
  
  # gauss hermite
  expect_equal(as.vector(round(test_outcome_gauss_hermite$estimate, 3)), 1.522)
  expect_equal(as.vector(round(attr(test_outcome_gauss_hermite$estimate, "variance"), 3)), .072)
  expect_equal(length(test_outcome_gauss_hermite$answers), 100)
  
  # riemann
  expect_equal(as.vector(round(test_outcome_riemann$estimate, 3)), 1.522)
  expect_equal(as.vector(round(attr(test_outcome_riemann$estimate, "variance"), 3)), .072)
  expect_equal(length(test_outcome_riemann$answers), 100)
})

test_that("true theta is 2, estimator is expected_aposteriori, U(-4, 4) prior", {
  # define true theta for simulation of answers
  true_theta <- 2
  
  # define item characteristics
  number_items <- 100
  number_dimensions <- 1
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  item_keys <- str_c("item", 1:number_items)
  
  guessing <- NULL
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions, dimnames = list(item_keys, NULL))
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1, dimnames = list(item_keys, NULL))
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(type = 'random', n = 3)
  stop_test <- list(max_n = 100)
  estimator <- 'expected_aposteriori'
  information_summary <- 'posterior_determinant'
  
  prior_form <- "uniform"  
  prior_parameters <- list(lower_bound = -4, upper_bound = 4)
  
  test_outcome_riemann <- with_random_seed(2, test_shadowcat)(true_theta, prior_form = prior_form, prior_parameters = prior_parameters, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary)
  
  # riemann
  expect_equal(as.vector(round(test_outcome_riemann$estimate, 3)), 2.283)
  expect_equal(as.vector(round(attr(test_outcome_riemann$estimate, "variance"), 3)), .112)
  expect_equal(length(test_outcome_riemann$answers), 100)
})

test_that("true theta is 2, estimator is expected_aposteriori, U(-1, 1) prior, with safe_eap", {
  # define true theta for simulation of answers
  true_theta <- 2
  
  # define item characteristics
  number_items <- 100
  number_dimensions <- 1
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  item_keys <- str_c("item", 1:number_items)
  
  guessing <- NULL
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions, dimnames = list(item_keys, NULL))
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1, dimnames = list(item_keys, NULL))
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(type = 'random', n = 3)
  stop_test <- list(max_n = 100)
  estimator <- 'expected_aposteriori'
  information_summary <- 'posterior_determinant'
  
  prior_form <- "uniform"  
  prior_parameters <- list(lower_bound = -1, upper_bound = 1)
  
  test_outcome_riemann <- with_random_seed(2, test_shadowcat)(true_theta, prior_form = prior_form, prior_parameters = prior_parameters, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, safe_eap = TRUE)
  
  # riemann
  expect_equal(as.vector(round(test_outcome_riemann$estimate, 3)), .963)
  expect_equal(as.vector(round(attr(test_outcome_riemann$estimate, "variance"), 3)), 0)
  expect_equal(length(test_outcome_riemann$answers), 100)
})

test_that("true theta is 1, 0, 2, estimator is maximum_aposteriori, N(c(0,0,0), diag(3) * 20) prior", {  
  # define true theta for simulation of answers
  true_theta <- c(1, 0, 2)
  
  # define item characteristics
  number_items <- 300
  number_dimensions <- 3
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  item_keys <- str_c("item", 1:number_items)
  
  guessing <- NULL
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions, dimnames = list(item_keys, NULL))
  alpha[1:100,2:3] <- 0
  alpha[101:200,c(1,3)] <- 0
  alpha[201:300,1:2] <- 0
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1, dimnames = list(item_keys, NULL))
  
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(type = 'random', n = 3)
  stop_test <- list(max_n = 300)
  estimator <- 'maximum_aposteriori'
  information_summary <- 'posterior_determinant'
  
  # define prior covariance matrix
  prior_form <- "normal"  
  prior_parameters <- list(mu = rep(0, number_dimensions), Sigma = diag(number_dimensions) * 20)
  
  test_outcome <- with_random_seed(3, test_shadowcat)(true_theta, prior_form = prior_form, prior_parameters = prior_parameters, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary)
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), c(.841, -.123, 1.947))
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3))[1:3],c(.064, .000, .000))
  expect_equal(length(test_outcome$answers), 300)
})

test_that("true theta is 1, 0, 2, estimator is maximum_aposteriori, N(c(-1,2,-2), diag(3)) prior", {  
  # define true theta for simulation of answers
  true_theta <- c(1, 0, 2)
  
  # define item characteristics
  number_items <- 300
  number_dimensions <- 3
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  item_keys <- str_c("item", 1:number_items)
  
  guessing <- NULL
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions, dimnames = list(item_keys, NULL))
  alpha[1:100,2:3] <- 0
  alpha[101:200,c(1,3)] <- 0
  alpha[201:300,1:2] <- 0
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1, dimnames = list(item_keys, NULL))
  
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(type = 'random', n = 3)
  stop_test <- list(max_n = 300)
  estimator <- 'maximum_aposteriori'
  information_summary <- 'posterior_determinant'
  
  # define prior covariance matrix
  prior_form <- "normal"  
  prior_parameters <- list(mu = c(-1, 2, -2), Sigma = diag(number_dimensions))
  
  test_outcome <- with_random_seed(3, test_shadowcat)(true_theta, prior_form = prior_form, prior_parameters = prior_parameters, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary)
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), c(.672, .219, 1.682))
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3))[1:3],c(.059, .000, .000))
  expect_equal(length(test_outcome$answers), 300)
})

test_that("true theta is 1, 0, 2, estimator is maximum_aposteriori, U(c(-1, -1, -1), c(1, 1, 1)) prior", {  
  # define true theta for simulation of answers
  true_theta <- c(1, 0, 2)
  
  # define item characteristics
  number_items <- 300
  number_dimensions <- 3
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  item_keys <- str_c("item", 1:number_items)
  
  guessing <- NULL
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions, dimnames = list(item_keys, NULL))
  alpha[1:100,2:3] <- 0
  alpha[101:200,c(1,3)] <- 0
  alpha[201:300,1:2] <- 0
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1, dimnames = list(item_keys, NULL))
  
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(type = 'random', n = 3)
  stop_test <- list(max_n = 300)
  estimator <- 'maximum_aposteriori'
  information_summary <- 'posterior_determinant'
  
  # define prior covariance matrix
  prior_form <- "uniform"  
  prior_parameters <- list(lower_bound = c(-1, -1, -1), upper_bound = c(1, 1, 1))
  
  test_outcome <- with_random_seed(3, test_shadowcat)(true_theta, prior_form = prior_form, prior_parameters = prior_parameters, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary)
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), c(.677, .059, 1.000))
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3))[1:3],c(.063, .000, .000))
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3))[7:9],c(.000, .000, .068))
  expect_equal(length(test_outcome$answers), 300)
})


test_that("true theta is 1, 0, 2, estimator is maximum_likelihood", {  
  # define true theta for simulation of answers
  true_theta <- c(1, 0, 2)
  
  # define item characteristics
  number_items <- 300
  number_dimensions <- 3
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  item_keys <- str_c("item", 1:number_items)
  
  guessing <- NULL
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions, dimnames = list(item_keys, NULL))
  alpha[1:100,2:3] <- 0
  alpha[101:200,c(1,3)] <- 0
  alpha[201:300,1:2] <- 0
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1, dimnames = list(item_keys, NULL))
  
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(type = 'random', n = 3)
  stop_test <- list(max_n = 300)
  estimator <- 'maximum_likelihood'
  information_summary <- 'determinant'
  
  test_outcome <- with_random_seed(3, test_shadowcat)(true_theta, prior_form = NULL, prior_parameters = NULL, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary)
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), c(.985, .019, 2.024))
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3))[1:3],c(.066, .000, .000))
  expect_equal(length(test_outcome$answers), 300)
  
  # defining prior has no effect on outcome
  test_outcome <- with_random_seed(3, test_shadowcat)(true_theta, prior_form = "normal", prior_parameters = list(mu = c(0, 0, 0), Sigma = diag(number_dimensions) * 2), model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary)
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), c(.985, .019, 2.024))
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3))[1:3],c(.066, .000, .000))
  expect_equal(length(test_outcome$answers), 300)  
})

test_that("true theta is 1, 0, 2, estimator is expected_aposteriori, N(c(0, 0, 0), diag(3) * 20) prior", {  
  # define true theta for simulation of answers
  true_theta <- c(1, 0, 2)
  
  # define item characteristics
  number_items <- 300
  number_dimensions <- 3
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  item_keys <- str_c("item", 1:number_items)
  
  guessing <- NULL
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions, dimnames = list(item_keys, NULL))
  alpha[1:100,2:3] <- 0
  alpha[101:200,c(1,3)] <- 0
  alpha[201:300,1:2] <- 0
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1, dimnames = list(item_keys, NULL))
  
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(type = 'random', n = 3)
  stop_test <- list(max_n = 300)
  estimator <- 'expected_aposteriori'
  information_summary <- 'posterior_determinant'
  
  prior_form <- "normal"  
  prior_parameters <- list(mu = rep(0, number_dimensions), Sigma = diag(number_dimensions) * 20)
  
  test_outcome_gauss_hermite <- with_random_seed(3, test_shadowcat)(true_theta, prior_form = prior_form, prior_parameters = prior_parameters, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, eap_estimation_procedure = "gauss_hermite_quad")
  test_outcome_riemann <- with_random_seed(3, test_shadowcat)(true_theta, prior_form = prior_form, prior_parameters = prior_parameters, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, eap_estimation_procedure = "riemannsum")
  
  # gauss hermite
  expect_equal(as.vector(round(test_outcome_gauss_hermite$estimate, 3)), c(.878, .305, 1.455))
  expect_equal(as.vector(round(attr(test_outcome_gauss_hermite$estimate, "variance"), 3))[1:3],c(.065, .000, .000))
  expect_equal(length(test_outcome_gauss_hermite$answers), 300)
  
  # riemann
  expect_equal(as.vector(round(test_outcome_riemann$estimate, 3)), c(.649, .073, 1.660))
  expect_equal(as.vector(round(attr(test_outcome_riemann$estimate, "variance"), 3))[1:3],c(.064, .000, .000))
  expect_equal(length(test_outcome_riemann$answers), 300)
})

test_that("true theta is 1, 0, 2, estimator is expected_aposteriori, N(c(-1, 1, -2), diag(3)) prior", {  
  # define true theta for simulation of answers
  true_theta <- c(1, 0, 2)
  
  # define item characteristics
  number_items <- 300
  number_dimensions <- 3
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  item_keys <- str_c("item", 1:number_items)
  
  guessing <- NULL
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions, dimnames = list(item_keys, NULL))
  alpha[1:100,2:3] <- 0
  alpha[101:200,c(1,3)] <- 0
  alpha[201:300,1:2] <- 0
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1, dimnames = list(item_keys, NULL))
  
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(type = 'random', n = 3)
  stop_test <- list(max_n = 300)
  estimator <- 'expected_aposteriori'
  information_summary <- 'posterior_determinant'
  
  prior_form <- "normal"  
  prior_parameters <- list(mu = c(-1, 1, -2), Sigma = diag(number_dimensions))
  
  test_outcome_gauss_hermite <- with_random_seed(3, test_shadowcat)(true_theta, prior_form = prior_form, prior_parameters = prior_parameters, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, eap_estimation_procedure = "gauss_hermite_quad")
  test_outcome_riemann <- with_random_seed(3, test_shadowcat)(true_theta, prior_form = prior_form, prior_parameters = prior_parameters, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, eap_estimation_procedure = "riemannsum")
  
  # gauss hermite
  expect_equal(as.vector(round(test_outcome_gauss_hermite$estimate, 3)), c(.583, .131, 1.609))
  expect_equal(as.vector(round(attr(test_outcome_gauss_hermite$estimate, "variance"), 3))[1:3],c(.059, .000, .000))
  expect_equal(length(test_outcome_gauss_hermite$answers), 300)
  
  # riemann
  expect_equal(as.vector(round(test_outcome_riemann$estimate, 3)), c(.627, .206, 1.697))
  expect_equal(as.vector(round(attr(test_outcome_riemann$estimate, "variance"), 3))[1:3],c(.061, .000, .000))
  expect_equal(length(test_outcome_riemann$answers), 300)
})

test_that("true theta is 1, 0, 2, estimator is expected_aposteriori, U(c(-4, -4, -4), c(4, 4, 4)) prior", {  
  # define true theta for simulation of answers
  true_theta <- c(1, 0, 2)
  
  # define item characteristics
  number_items <- 300
  number_dimensions <- 3
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  item_keys <- str_c("item", 1:number_items)
  
  guessing <- NULL
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions, dimnames = list(item_keys, NULL))
  alpha[1:100,2:3] <- 0
  alpha[101:200,c(1,3)] <- 0
  alpha[201:300,1:2] <- 0
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1, dimnames = list(item_keys, NULL))
  
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(type = 'random', n = 3)
  stop_test <- list(max_n = 300)
  estimator <- 'expected_aposteriori'
  information_summary <- 'posterior_determinant'
  
  prior_form <- "uniform"  
  prior_parameters <- list(lower_bound = c(-4, -4, -4), upper_bound = c(4, 4, 4))
  
  test_outcome_riemann <- with_random_seed(3, test_shadowcat)(true_theta, prior_form = prior_form, prior_parameters = prior_parameters, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, eap_estimation_procedure = "riemannsum")
  
  # riemann
  expect_equal(as.vector(round(test_outcome_riemann$estimate, 3)), c(1.067, .196, 1.987))
  expect_equal(as.vector(round(attr(test_outcome_riemann$estimate, "variance"), 3))[1:3],c(.068, .000, .000))
  expect_equal(length(test_outcome_riemann$answers), 300)
})

test_that("true theta is 1, 0, 2, estimator is expected_aposteriori, U(c(-1, -1, -1), c(1, 1, 1)) prior, with safe_eap", {  
  # define true theta for simulation of answers
  true_theta <- c(1, 0, 2)
  
  # define item characteristics
  number_items <- 300
  number_dimensions <- 3
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  item_keys <- str_c("item", 1:number_items)
  
  guessing <- NULL
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions, dimnames = list(item_keys, NULL))
  alpha[1:100,2:3] <- 0
  alpha[101:200,c(1,3)] <- 0
  alpha[201:300,1:2] <- 0
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1, dimnames = list(item_keys, NULL))
  
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(type = 'random', n = 3)
  stop_test <- list(max_n = 300)
  estimator <- 'expected_aposteriori'
  information_summary <- 'posterior_determinant'
  
  prior_form <- "uniform"  
  prior_parameters <- list(lower_bound = c(-1, -1, -1), upper_bound = c(1, 1, 1))
  
  test_outcome_riemann <- with_random_seed(3, test_shadowcat)(true_theta, prior_form = prior_form, prior_parameters = prior_parameters, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, eap_estimation_procedure = "riemannsum", safe_eap = TRUE)
  
  # riemann
  expect_equal(as.vector(round(test_outcome_riemann$estimate, 3)), c(.876, -.349, .939))
  expect_equal(as.vector(round(attr(test_outcome_riemann$estimate, "variance"), 3))[1:3],c(.000, .000, .000))
  expect_equal(as.vector(round(attr(test_outcome_riemann$estimate, "variance"), 3))[7:9],c(.000, .000, .000))
  expect_equal(length(test_outcome_riemann$answers), 300)
})


test_that("items load on three dimensions", {  
  # define true theta for simulation of answers
  true_theta <- c(1, 0, 2)
  
  # define item characteristics
  number_items <- 300
  number_dimensions <- 3
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  item_keys <- str_c("item", 1:number_items)
  
  guessing <- NULL
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions, dimnames = list(item_keys, NULL))
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1, dimnames = list(item_keys, NULL))
  
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(type = 'random', n = 3)
  stop_test <- list(max_n = 300)
  estimator <- 'maximum_aposteriori'
  information_summary <- 'posterior_determinant'
  
  # define prior covariance matrix
  prior_form = "normal"
  prior_parameters = list(mu = rep(0, number_dimensions), Sigma = diag(number_dimensions) * 20)
  
  test_outcome <- with_random_seed(3, test_shadowcat)(true_theta, prior_form = prior_form, prior_parameters = prior_parameters, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary)
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), c(.454, .719, 1.745))
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3))[1:3],c(.178, -.088, -.081))
  expect_equal(length(test_outcome$answers), 300)
})

test_that("true theta is 2, 2, 2", {  
  # define true theta for simulation of answers
  true_theta <- c(2, 2, 2)
  
  # define item characteristics
  number_items <- 300
  number_dimensions <- 3
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  item_keys <- str_c("item", 1:number_items)
  
  guessing <- NULL
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions, dimnames = list(item_keys, NULL))
  alpha[1:100,2:3] <- 0
  alpha[101:200,c(1,3)] <- 0
  alpha[201:300,1:2] <- 0
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1, dimnames = list(item_keys, NULL))
  
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(type = 'random', n = 3)
  stop_test <- list(max_n = 300)
  estimator <- 'maximum_aposteriori'
  information_summary <- 'posterior_determinant'
  
  prior_form = "normal"
  prior_parameters = list(mu = rep(0, number_dimensions), Sigma = diag(number_dimensions) * 20)
  
  test_outcome <- with_random_seed(3, test_shadowcat)(true_theta, prior_form = prior_form, prior_parameters = prior_parameters, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary)
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), c(2.287, 1.825, 1.672))
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3))[1:3],c(.11, .000, .000))
  expect_equal(length(test_outcome$answers), 300)
})

test_that("with constraints max_n 260", {  
  # define true theta for simulation of answers
  true_theta <- c(-2, 1, 2)
  
  # define item characteristics
  number_items <- 300
  number_dimensions <- 3
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  item_keys <- str_c("item", 1:number_items)
  
  guessing <- NULL
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions, dimnames = list(item_keys, NULL))
  alpha[1:100,2:3] <- 0
  alpha[101:200,c(1,3)] <- 0
  alpha[201:300,1:2] <- 0
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1, dimnames = list(item_keys, NULL))
  
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(type = 'random_by_dimension', n_by_dimension = 3, n = 9)
  stop_test <- list(max_n = 260)
  estimator <- 'maximum_aposteriori'
  information_summary <- 'posterior_determinant'
  
  prior_form = "normal"
  prior_parameters = list(mu = rep(0, number_dimensions), Sigma = diag(number_dimensions) * 20)
  
  #create item characteristics and constraints
  characteristics <- data.frame(content = c(rep('depression', number_items / 3), rep('anxiety', number_items / 3), rep('somatic', number_items / 3)))
  constraints <- list(list(name = 'content/depression',
                           op = '><',
                           target = c(50, 75)),
                      list(name = 'content/somatic',
                           op = '><',
                           target = c(75, 90)))
  
  test_outcome <- with_random_seed(3, test_shadowcat)(true_theta, prior_form = prior_form, prior_parameters = prior_parameters, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, constraints_and_characts = list(characteristics = characteristics, constraints = constraints))
  indices_administered <- as.numeric(sapply(names(test_outcome$answers), substring, 5))
  
  number_depression_items <- sum(indices_administered <= 100)
  number_anxiety_items <- sum(indices_administered > 100 & indices_administered <= 200)
  number_somatic_items <- sum(indices_administered > 200)
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), c(-1.900, 0.949, 1.724))
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3))[1:3], c(.104, .000, .000))
  expect_equal(length(test_outcome$answers), 260)
  expect_equal(number_depression_items, 75)
  expect_equal(number_anxiety_items, 95)
  expect_equal(number_somatic_items, 90) 
})

test_that("with constraints max_n 130", {  
  # define true theta for simulation of answers
  true_theta <- c(-2, 1, 2)
  
  # define item characteristics
  number_items <- 300
  number_dimensions <- 3
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  item_keys <- str_c("item", 1:number_items)
  
  guessing <- NULL
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions, dimnames = list(item_keys, NULL))
  alpha[1:100,2:3] <- 0
  alpha[101:200,c(1,3)] <- 0
  alpha[201:300,1:2] <- 0
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1, dimnames = list(item_keys, NULL))
  
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(type = 'random_by_dimension', n_by_dimension = 3, n = 9)
  stop_test <- list(max_n = 130)
  estimator <- 'maximum_aposteriori'
  information_summary <- 'posterior_determinant'
  
  prior_form = "normal"
  prior_parameters = list(mu = rep(0, number_dimensions), Sigma = diag(number_dimensions) * 20)
  
  #create item characteristics and constraints
  characteristics <- data.frame(content = c(rep('depression', number_items / 3), rep('anxiety', number_items / 3), rep('somatic', number_items / 3)))
  constraints <- list(list(name = 'content/depression',
                           op = '><',
                           target = c(50, 75)),
                      list(name = 'content/somatic',
                           op = '><',
                           target = c(75, 90)))
  
  test_outcome <- with_random_seed(3, test_shadowcat)(true_theta, prior_form = prior_form, prior_parameters = prior_parameters, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, constraints_and_characts = list(characteristics = characteristics, constraints = constraints))
  indices_administered <- as.numeric(sapply(names(test_outcome$answers), substring, 5))
  
  number_depression_items <- sum(indices_administered <= 100)
  number_anxiety_items <- sum(indices_administered > 100 & indices_administered <= 200)
  number_somatic_items <- sum(indices_administered > 200)
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), c(-2.303, .774, 1.822))
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3))[1:3], c(.146, .000, .000))
  expect_equal(length(test_outcome$answers), 130)
  expect_equal(number_depression_items, 50)
  expect_equal(number_anxiety_items, 5)
  expect_equal(number_somatic_items, 75) 
})

test_that("start n is zero, no constraints", {  
  # define true theta for simulation of answers
  true_theta <- c(1, 0, 2)
  
  # define item characteristics
  number_items <- 300
  number_dimensions <- 3
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  item_keys <- str_c("item", 1:number_items)
  
  guessing <- NULL
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions, dimnames = list(item_keys, NULL))
  alpha[1:100,2:3] <- 0
  alpha[101:200,c(1,3)] <- 0
  alpha[201:300,1:2] <- 0
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1, dimnames = list(item_keys, NULL))
  
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(n = 0)
  stop_test <- list(max_n = 300)
  estimator <- 'maximum_aposteriori'
  information_summary <- 'posterior_determinant'
  
  prior_form = "normal"
  prior_parameters = list(mu = rep(0, number_dimensions), Sigma = diag(number_dimensions) * 20)
  
  test_outcome <- with_random_seed(3, test_shadowcat)(true_theta, prior_form = prior_form, prior_parameters = prior_parameters, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, initital_estimate = rep(.2, number_dimensions), initial_variance = diag(number_dimensions) * 20)
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), c(1.114, -0.018, 1.725))
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3))[1:3],c(.068, .000, .000))
  expect_equal(length(test_outcome$answers), 300)
})

test_that("start n is zero, with constraints", {  
  # define true theta for simulation of answers
  true_theta <- c(-2, 1, 2)
  
  # define item characteristics
  number_items <- 300
  number_dimensions <- 3
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  item_keys <- str_c("item", 1:number_items)
  
  guessing <- NULL
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions, dimnames = list(item_keys, NULL))
  alpha[1:100,2:3] <- 0
  alpha[101:200,c(1,3)] <- 0
  alpha[201:300,1:2] <- 0
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1, dimnames = list(item_keys, NULL))
  
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(n = 0)
  stop_test <- list(max_n = 130)
  estimator <- 'maximum_aposteriori'
  information_summary <- 'posterior_determinant'
  
  prior_form = "normal"
  prior_parameters = list(mu = rep(0, number_dimensions), Sigma = diag(number_dimensions) * 20)
  
  #create item characteristics and constraints
  characteristics <- data.frame(content = c(rep('depression', number_items / 3), rep('anxiety', number_items / 3), rep('somatic', number_items / 3)))
  constraints <- list(list(name = 'content/depression',
                           op = '><',
                           target = c(50, 75)),
                      list(name = 'content/somatic',
                           op = '><',
                           target = c(75, 90)))
  
  test_outcome <- with_random_seed(3, test_shadowcat)(true_theta, prior_form = prior_form, prior_parameters = prior_parameters, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, constraints_and_characts = list(characteristics = characteristics, constraints = constraints, initital_estimate = rep(.2, number_dimensions), initial_variance = diag(number_dimensions) * 20))
  indices_administered <- as.numeric(sapply(names(test_outcome$answers), substring, 5))
  
  number_depression_items <- sum(indices_administered <= 100)
  number_anxiety_items <- sum(indices_administered > 100 & indices_administered <= 200)
  number_somatic_items <- sum(indices_administered > 200)
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), c(-1.975, .897, 2.496))
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3))[1:3],c(.122, .000, .000))
  expect_equal(length(test_outcome$answers), 130)
  expect_equal(number_depression_items, 50)
  expect_equal(number_anxiety_items, 5)
  expect_equal(number_somatic_items, 75) 
})

context("check stop rule")

test_that("stop rule is number of items", {
  # define true theta for simulation of answers
  true_theta <- 0
  
  # define item characteristics
  number_items <- 50
  number_dimensions <- 1
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  item_keys <- str_c("item", 1:number_items)
  
  guessing <- matrix(c(rep(.1, number_items / 2), rep(.2, number_items / 2)), ncol = 1, dimnames = list(item_keys, NULL))
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions, dimnames = list(item_keys, NULL))
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1, dimnames = list(item_keys, NULL))
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(type = 'random', n = 5)
  stop_test <- list(max_n = 10)
  estimator <- 'maximum_aposteriori'
  information_summary <- 'posterior_determinant'
  
  prior_form = "normal"
  prior_parameters = list(mu = rep(0, number_dimensions), Sigma = diag(number_dimensions))
  
  test_outcome <- with_random_seed(2, test_shadowcat)(true_theta, prior_form = prior_form, prior_parameters = prior_parameters, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary)
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), .405)
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3)), .315)
  expect_equal(names(test_outcome$answers), str_c("item", c(10, 30, 48, 7, 24, 5, 17, 6, 45, 47)))
  expect_equal(unname(unlist(test_outcome$answers)), c(1, 0, 1, 1, 0, 1, 0, 1, 0, 1))
})

test_that("stop rule is variance", {
  # define true theta for simulation of answers
  true_theta <- 0
  
  # define item characteristics
  number_items <- 50
  number_dimensions <- 1
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  item_keys <- str_c("item", 1:number_items)
  
  guessing <- matrix(c(rep(.1, number_items / 2), rep(.2, number_items / 2)), ncol = 1, dimnames = list(item_keys, NULL))
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions, dimnames = list(item_keys, NULL))
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1, dimnames = list(item_keys, NULL))
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(type = 'random', n = 5)
  stop_test <- list(target = .5, max_n = 50)
  estimator <- 'maximum_aposteriori'
  information_summary <- 'posterior_determinant'
  
  prior_form = "normal"
  prior_parameters = list(mu = rep(0, number_dimensions), Sigma = diag(number_dimensions))
  
  test_outcome <- with_random_seed(2, test_shadowcat)(true_theta, prior_form = prior_form, prior_parameters = prior_parameters, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary)
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), 0.329)
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3)), 0.47)
  expect_equal(names(test_outcome$answers), str_c("item", c(10, 30, 48, 7, 24, 5, 17)))
  expect_equal(unname(unlist(test_outcome$answers)), c(1, 0, 1, 1, 0, 1, 0))
})

test_that("items that load on finished dimensions are not administered", {
  # define true theta for simulation of answers
  true_theta <- c(0, 0, 0)
  
  # define item characteristics
  number_items <- 300
  number_dimensions <- 3
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  item_keys <- str_c("item", 1:number_items)
  
  guessing <- matrix(c(rep(.1, number_items / 2), rep(.2, number_items / 2)), ncol = 1, dimnames = list(item_keys, NULL))
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions, dimnames = list(item_keys, NULL))
  alpha[1:100, 1:2] <- 0
  alpha[101:300, 3] <- 0
  
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1, dimnames = list(item_keys, NULL))
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(type = 'random_by_dimension', n_by_dimension = 2, n = 6)
  stop_test <- list(target = c(.1, .1, 2), max_n = 100)
  estimator <- 'maximum_aposteriori'
  information_summary <- 'posterior_determinant'
  
  prior_form = "normal"
  prior_parameters = list(mu = rep(0, number_dimensions), Sigma = diag(number_dimensions) * 20)
  
  test_outcome <- with_random_seed(2, test_shadowcat)(true_theta, prior_form = prior_form, prior_parameters = prior_parameters, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary)
  item_numbers <- sapply(names(test_outcome$answers), function(item) { as.numeric(substr(x = item, start = 5, stop = nchar(item))) })
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), c(-.635, .487, .482))
  expect_equal(diag(round(attr(test_outcome$estimate, "variance"), 3)), c(.250, .221, 1.754))
  expect_equal(all(item_numbers[9:100] > 100), TRUE)
})


test_that("stop rule is variance and minimum number of items is taken into account", {
  # define true theta for simulation of answers
  true_theta <- 0
  
  # define item characteristics
  number_items <- 50
  number_dimensions <- 1
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  item_keys <- str_c("item", 1:number_items)
  
  guessing <- matrix(c(rep(.1, number_items / 2), rep(.2, number_items / 2)), ncol = 1, dimnames = list(item_keys, NULL))
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions, dimnames = list(item_keys, NULL))
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1, dimnames = list(item_keys, NULL))
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(type = 'random', n = 5)
  stop_test <- list(target = .5, max_n = 50, min_n = 10)
  estimator <- 'maximum_aposteriori'
  information_summary <- 'posterior_determinant'
  
  prior_form = "normal"
  prior_parameters = list(mu = rep(0, number_dimensions), Sigma = diag(number_dimensions))
  
  test_outcome <- with_random_seed(2, test_shadowcat)(true_theta, prior_form = prior_form, prior_parameters = prior_parameters, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary)
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), 0.405)
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3)), .315)
  expect_equal(names(test_outcome$answers), str_c("item", c(10, 30, 48, 7, 24, 5, 17, 6, 45, 47)))
  expect_equal(unname(unlist(test_outcome$answers)), c(1, 0, 1, 1, 0, 1, 0, 1, 0, 1))
})

test_that("stop rule is cutoff", {  
  # define true theta for simulation of answers
  true_theta <- c(.2, 0, -2)
  
  # define item characteristics
  number_items <- 300
  number_dimensions <- 3
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  item_keys <- str_c("item", 1:number_items)
  
  guessing <- NULL
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions, dimnames = list(item_keys, NULL))
  alpha[1:100,2:3] <- 0
  alpha[101:200,c(1,3)] <- 0
  alpha[201:300,1:2] <- 0
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1, dimnames = list(item_keys, NULL))
  
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(type = 'random_by_dimension', n_by_dimension = 3, n = 9)
  stop_test <- list(max_n = 300, cutoffs = with_random_seed(2, matrix)(runif(900, 1, 2), ncol = 3))
  estimator <- 'maximum_aposteriori'
  information_summary <- 'posterior_determinant'
  
  prior_form = "normal"
  prior_parameters = list(mu = rep(0, number_dimensions), Sigma = diag(number_dimensions) * 20)
  
  test_outcome <- with_random_seed(3, test_shadowcat)(true_theta, prior_form = prior_form, prior_parameters = prior_parameters, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary)
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), c(.088, .047, -6.075))
  expect_equal(diag(round(attr(test_outcome$estimate, "variance"), 3)), c(.214, .223, 5.241))
  expect_equal(length(test_outcome$answers), 32)
})

test_that("invalid input", {  
  # define item characteristics
  number_items <- 300
  number_dimensions <- 3
  number_answer_categories <- 2
  item_keys <- str_c("item", 1:number_items)
  guessing <- NULL
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions, dimnames = list(item_keys, NULL))
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1, dimnames = list(item_keys, NULL))
  beta_wrong_na_pattern <- matrix(with_random_seed(2, rnorm)(number_items * 2), nrow = number_items, ncol = 2, dimnames = list(item_keys, NULL))
  beta_wrong_na_pattern[1,1] <- NA
  
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(type = 'random_by_dimension', n_by_dimension = 3, n = 9)
  stop_test <- list(max_n = 300, cutoffs = with_random_seed(2, matrix)(runif(903, 1, 2), ncol = 3))
  estimator <- 'maximum_aposteriori'
  information_summary <- 'posterior_determinant'
  
  # define prior covariance matrix
  prior_form = "normal"
  prior_parameters = list(mu = rep(0, number_dimensions), Sigma = diag(number_dimensions) * 20)
  
  estimate <- rep(.5, number_dimensions)
  variance <- as.vector(diag(number_dimensions) * 2)
  
  characteristics <- data.frame(time = with_random_seed(2, rnorm)(200),
                                type = with_random_seed(2, sample)(c('depression', 'insomnia', 'anxiety'), 200, TRUE),
                                stressful = rep(c(0, 0, 1, 0, 0), 200/5))
  
  constraints1 <- list(list(name = 'time',
                            op = c('>', '<')),
                       list(name = 'depression',
                            op = c('>', '<'),
                            target = c(2, 5)),
                       list(name = 'insomnia',
                            op = c('>', '<'),
                            target = c(2, 5)),
                       list(name = 'stressful',
                            op = '<',
                            target = "3"))
  
  constraints2 <- list(list(name = 'time',
                            op = c('>', '<'),
                            target = c(6, 10)),
                       list(name = 'depression',
                            op = c('>', '<'),
                            target = c(2, 5)),
                       list(name = 'insomnia',
                            op = c('>', '<'),
                            target = c(2, 5)),
                       list(name = 'stressful',
                            op = '<',
                            target = "3"))
   
  error_message_estimate <- shadowcat(answers = NULL, estimate = NULL, variance, model, alpha, beta, start_items, stop_test, estimator, information_summary, prior_form, prior_parameters, guessing, eta)
  error_message_variance <- shadowcat(answers = NULL, estimate, variance = NULL, model, alpha, beta, start_items, stop_test, estimator, information_summary, prior_form, prior_parameters, guessing, eta)
  error_message_variance_class <- shadowcat(answers = NULL, estimate, variance = diag(3), model, alpha, beta, start_items, stop_test, estimator, information_summary, prior_form, prior_parameters, guessing, eta)
  error_message_variance_length <- shadowcat(answers = NULL, estimate, variance = c(2, 2, 2), model, alpha, beta, start_items, stop_test, estimator, information_summary, prior_form, prior_parameters, guessing, eta) 
  error_message_model <- shadowcat(answers = NULL, estimate, variance, model = NULL, alpha, beta, start_items, stop_test, estimator, information_summary, prior_form, prior_parameters, guessing, eta)
  error_message_alpha <- shadowcat(answers = NULL, estimate, variance, model, alpha = NULL, beta, start_items, stop_test, estimator, information_summary, prior_form, prior_parameters, guessing, eta)
  error_message_start_items <- shadowcat(answers = NULL, estimate, variance, model, alpha, beta, start_items = NULL, stop_test, estimator, information_summary, prior_form, prior_parameters, guessing, eta)
  error_message_stop_test <- shadowcat(answers = NULL, estimate, variance, model, alpha, beta, start_items, stop_test = NULL, estimator, information_summary, prior_form, prior_parameters, guessing, eta)
  error_message_estimator <- shadowcat(answers = NULL, estimate, variance, model, alpha, beta, start_items, stop_test, estimator = NULL, information_summary, prior_form, prior_parameters, guessing, eta)
  error_message_information_summary <- shadowcat(answers = NULL, estimate, variance, model, alpha, beta, start_items, stop_test, estimator, information_summary = NULL, prior_form, prior_parameters, guessing, eta)
  error_message_n_by_dimension <- shadowcat(answers = NULL, estimate, variance, model, alpha, beta, start_items = list(type = "random_by_dimension", n_by_dimension = c(2, 3)), stop_test, estimator, information_summary, prior_form, prior_parameters, guessing, eta)
  error_message_alpha_matrix <- shadowcat(answers = NULL, estimate, variance, model, alpha = 1:10, beta, start_items, stop_test, estimator, information_summary, prior_form, prior_parameters, guessing, eta)
  error_message_alpha_rownames <- shadowcat(answers = NULL, estimate, variance, model, alpha = matrix(1:10), beta, start_items, stop_test, estimator, information_summary, prior_form, prior_parameters, guessing, eta)
  error_message_beta_matrix <- shadowcat(answers = NULL, estimate, variance, model, alpha, beta = 1:10, start_items, stop_test, estimator, information_summary, prior_form, prior_parameters, guessing, eta)
  error_message_beta_rownames <- shadowcat(answers = NULL, estimate, variance, model, alpha, beta = matrix(1:10), start_items, stop_test, estimator, information_summary, prior_form, prior_parameters, guessing, eta)
  error_message_eta_matrix <- shadowcat(answers = NULL, estimate, variance, model, alpha, beta, start_items, stop_test, estimator, information_summary, prior_form, prior_parameters, guessing, eta = 1:10)
  error_message_eta_rownames <- shadowcat(answers = NULL, estimate, variance, model, alpha, beta, start_items, stop_test, estimator, information_summary, prior_form, prior_parameters, guessing, eta = matrix(1:10))
  error_message_guessing_matrix <- shadowcat(answers = NULL, estimate, variance, model, alpha, beta, start_items, stop_test, estimator, information_summary, prior_form, prior_parameters, guessing = 1:10, eta)
  error_message_guessing_rownames <- shadowcat(answers = NULL, estimate, variance, model, alpha, beta, start_items, stop_test, estimator, information_summary, prior_form, prior_parameters, guessing = matrix(1:10), eta)
  error_message_guessing_ncol <- shadowcat(answers = NULL, estimate, variance, model, alpha, beta, start_items, stop_test, estimator, information_summary, prior_form, prior_parameters, guessing = matrix(1:10, ncol = 2, dimnames = list(c("a", "b", "c", "d", "e"), NULL)), eta)
  error_message_unequal_rownames <- shadowcat(answers = NULL, estimate, variance, model, alpha, beta, start_items, stop_test, estimator, information_summary, prior_form, prior_parameters, guessing = matrix(1:300, ncol = 1, dimnames = list(str_c("it", 1:300), NULL)), eta)
  error_message_beta_na_pattern <- shadowcat(answers = NULL, estimate, variance, model, alpha, beta = beta_wrong_na_pattern, start_items, stop_test, estimator, information_summary, prior_form, prior_parameters, guessing, eta)
  error_message_eta_na_pattern <- shadowcat(answers = NULL, estimate, variance, model, alpha, beta, start_items, stop_test, estimator, information_summary, prior_form, prior_parameters, guessing, eta = beta_wrong_na_pattern)
  error_message_incorrect_dim_estimate <- shadowcat(answers = NULL, estimate = c(.5, .5), variance = as.vector(diag(2)), model, alpha, beta, start_items, stop_test, estimator, information_summary, prior_form, prior_parameters, guessing, eta)
  error_message_incorrect_dim_variance <- shadowcat(answers = NULL, estimate, variance = as.vector(diag(2)), model, alpha, beta, start_items, stop_test, estimator, information_summary, prior_form, prior_parameters, guessing, eta)
  error_message_variance_not_positive_definite <- shadowcat(answers = NULL, estimate, variance = as.vector(diag(3)) * -1, model, alpha, beta, start_items, stop_test, estimator, information_summary, prior_form, prior_parameters, guessing, eta)
  error_message_model_unknown <- shadowcat(answers = NULL, estimate, variance, model = "PLM", alpha, beta, start_items, stop_test, estimator, information_summary, prior_form, prior_parameters, guessing, eta)
  error_message_beta <- shadowcat(answers = NULL, estimate, variance, model, alpha, beta = NULL, start_items, stop_test, estimator, information_summary, prior_form, prior_parameters, guessing, eta)
  error_message_beta_theta <- shadowcat(answers = NULL, estimate, variance, model = "GPCM", alpha, beta = NULL, start_items, stop_test, estimator, information_summary, prior_form, prior_parameters, guessing, eta) 
  error_message_beta_theta_mismatch <- shadowcat(answers = NULL, estimate, variance, model = "GPCM", alpha, beta = cbind(beta, beta + .1), start_items, stop_test, estimator, information_summary, prior_form, prior_parameters, guessing, eta = cbind(beta, beta + .1))
  error_message_prior_form <- shadowcat(answers = NULL, estimate, variance, model, alpha, beta, start_items, stop_test, estimator = "maximum_aposteriori", information_summary, prior_form = NULL, prior_parameters, guessing, eta)
  error_message_prior_parameters <- shadowcat(answers = NULL, estimate, variance, model, alpha, beta, start_items, stop_test, estimator = "maximum_aposteriori", information_summary, prior_form, prior_parameters = NULL, guessing, eta)
  error_message_prior_form_type <- shadowcat(answers = NULL, estimate, variance, model, alpha, beta, start_items, stop_test, estimator = "maximum_aposteriori", information_summary, prior_form = "gamma", prior_parameters, guessing, eta)
  error_message_prior_form_par_mismatch_uniform <- shadowcat(answers = NULL, estimate, variance, model, alpha, beta, start_items, stop_test, estimator = "maximum_aposteriori", information_summary, prior_form = "uniform", prior_parameters, guessing, eta)
  error_message_prior_form_par_mismatch_normal <- shadowcat(answers = NULL, estimate, variance, model, alpha, beta, start_items, stop_test, estimator = "maximum_aposteriori", information_summary, prior_form = "normal", prior_parameters = list(lower_bound = rep(-3, 3), upper_bound = rep(3, 3)), guessing, eta)
  error_message_prior_parameters_mu <- shadowcat(answers = NULL, estimate, variance, model, alpha, beta, start_items, stop_test, estimator = "maximum_aposteriori", information_summary, prior_form = "normal", prior_parameters = list(mu = 0, Sigma = diag(3)), guessing, eta)
  error_message_prior_parameters_sigma1 <- shadowcat(answers = NULL, estimate, variance, model, alpha, beta, start_items, stop_test, estimator = "maximum_aposteriori", information_summary, prior_form = "normal", prior_parameters = list(mu = rep(0, 3), Sigma = rep(.1, 3)), guessing, eta)
  error_message_prior_parameters_sigma2 <- shadowcat(answers = NULL, estimate, variance, model, alpha, beta, start_items, stop_test, estimator = "maximum_aposteriori", information_summary, prior_form = "normal", prior_parameters = list(mu = rep(0, 3), Sigma = as.matrix(rep(.1, 3))), guessing, eta)
  error_message_prior_parameters_sigma3 <- shadowcat(answers = NULL, estimate, variance, model, alpha, beta, start_items, stop_test, estimator = "maximum_aposteriori", information_summary, prior_form = "normal", prior_parameters = list(mu = rep(0, 3), Sigma = diag(c(-.1, .1, .5))), guessing, eta)
  error_message_prior_parameters_bounds <- shadowcat(answers = NULL, estimate, variance, model, alpha, beta, start_items, stop_test, estimator = "maximum_aposteriori", information_summary, prior_form = "uniform", prior_parameters = list(lower_bound = -3, upper_bound = 3), guessing, eta)
  error_message_max_n <- shadowcat(answers = NULL, estimate, variance, model = "GPCM", alpha, beta, start_items, stop_test = list("variance" = .6), estimator, information_summary, prior_form, prior_parameters, guessing, eta)
  error_message_max_n_too_large <- shadowcat(answers = NULL, estimate, variance, model = "GPCM", alpha, beta, start_items, stop_test = list(max_n = 500), estimator, information_summary, prior_form, prior_parameters, guessing, eta)
  error_message_invalid_cutoffs <- shadowcat(answers = NULL, estimate, variance, model = "GPCM", alpha, beta, start_items, stop_test = list(max_n = 300, cutoffs = matrix(c(NA, .3, .5, .9), ncol = 1)), estimator, information_summary, prior_form, prior_parameters, guessing, eta)
  error_message_start_items_0 <- shadowcat(answers = NULL, estimate, variance, model = "GPCM", alpha, beta, start_items = list(n = 0), stop_test, estimator, information_summary = "posterior_expected_kullback_leibler", prior_form, prior_parameters, guessing, eta)
  error_message_start_items_random_by_dim_scalar <- shadowcat(answers = NULL, estimate, variance, model = "GPCM", alpha, beta, start_items = list(type = "random_by_dimension", n_by_dimension = 2, n = 2), stop_test, estimator, information_summary, prior_form, prior_parameters, guessing, eta)
  error_message_start_items_random_by_dim_vector <- shadowcat(answers = NULL, estimate, variance, model = "GPCM", alpha, beta, start_items = list(type = "random_by_dimension", n_by_dimension = c(2, 3, 2), n = 9), stop_test, estimator, information_summary, prior_form, prior_parameters, guessing, eta)
  error_message_cutoffs <- shadowcat(answers = NULL, estimate, variance, model = "GPCM", alpha, beta, start_items, stop_test = list(max_n = 70, cutoffs = 1:70), estimator, information_summary, prior_form, prior_parameters, guessing, eta)
  error_message_answers <- shadowcat(answers = list(item1 = 0, item22 = 1, item301 = 0), estimate, variance, model, alpha, beta, start_items, stop_test, estimator, information_summary, prior_form, prior_parameters, guessing, eta)
  error_message_estimator_unknown <- shadowcat(answers = NULL, estimate, variance, model, alpha, beta, start_items, stop_test, estimator = "ML", information_summary, prior_form, prior_parameters, guessing, eta)
  error_message_information_summary_unknown <- shadowcat(answers = NULL, estimate, variance, model, alpha, beta, start_items, stop_test, estimator, information_summary = "PT", prior_form, prior_parameters, guessing, eta)
  error_message_estimator_ml_posterior_mixed <- shadowcat(answers = NULL, estimate, variance, model, alpha, beta, start_items, stop_test, estimator = "maximum_likelihood", information_summary = "posterior_determinant", prior_form, prior_parameters, guessing, eta)
  error_message_bounds_should_be_null <- shadowcat(answers = NULL, estimate, variance, model, alpha, beta, start_items, stop_test, estimator, information_summary, prior_form, prior_parameters, guessing, eta, lower_bound = rep(-3, 3), upper_bound = rep(3, 3))
  error_message_lower_bound <- shadowcat(answers = NULL, estimate, variance, model, alpha, beta, start_items, stop_test, estimator, information_summary, prior_form, prior_parameters, guessing, eta, lower_bound = c(-2,-3))
  error_message_upper_bound <- shadowcat(answers = NULL, estimate, variance, model, alpha, beta, start_items, stop_test, estimator, information_summary, prior_form, prior_parameters, guessing, eta, upper_bound = c(-2,-3))
  error_message_only_characteristics <- shadowcat(answers = NULL, estimate, variance, model, alpha, beta, start_items, stop_test, estimator, information_summary, prior_form, prior_parameters, guessing, eta, constraints_and_characts = list(characteristics = characteristics))
  error_message_constraints_structure <- shadowcat(answers = NULL, estimate, variance, model, alpha, beta, start_items, stop_test, estimator, information_summary, prior_form, prior_parameters, guessing, eta, constraints_and_characts = list(characteristics = characteristics, constraints = constraints1))
  error_message_wrong_characts_and_constraints <- shadowcat(answers = NULL, estimate, variance, model, alpha, beta, start_items, stop_test, estimator, information_summary, prior_form, prior_parameters, guessing, eta, constraints_and_characts = list(characteristics = characteristics, constraints = constraints2))
  
  
  expect_equal(error_message_estimate$errors$estimate, "is missing")
  expect_equal(error_message_variance$errors$variance, "is missing")
  expect_equal(error_message_variance_class$errors$variance, "should be entered as vector")
  expect_equal(error_message_variance_length$errors$variance, "should be a covariance matrix turned into a vector")
  expect_equal(error_message_model$errors$model, "is missing")
  expect_equal(error_message_alpha$errors$alpha, "is missing")
  expect_equal(error_message_start_items$errors$start_items, "is missing")
  expect_equal(error_message_stop_test$errors$stop_test, "is missing")
  expect_equal(error_message_estimator$errors$estimator, "is missing")
  expect_equal(error_message_information_summary$errors$information_summary, "is missing")
  expect_equal(error_message_n_by_dimension$errors$start_items, "length of n_by_dimension should be a scalar or vector of the length of estimate")
  expect_equal(error_message_alpha_matrix$errors$alpha, "should be a matrix with item keys as row names")
  expect_equal(error_message_alpha_rownames$errors$alpha, "should be a matrix with item keys as row names")
  expect_equal(error_message_beta_matrix$errors$beta, "should be a matrix with item keys as row names")
  expect_equal(error_message_beta_rownames$errors$beta, "should be a matrix with item keys as row names")
  expect_equal(error_message_eta_matrix$errors$eta, "should be a matrix with item keys as row names")
  expect_equal(error_message_eta_rownames$errors$eta, "should be a matrix with item keys as row names")
  expect_equal(error_message_guessing_matrix$errors$guessing, "should be a single column matrix with item keys as row names")
  expect_equal(error_message_guessing_rownames$errors$guessing, "should be a single column matrix with item keys as row names")
  expect_equal(error_message_guessing_ncol$errors$guessing, "should be a single column matrix with item keys as row names")
  expect_equal(error_message_unequal_rownames$errors$alpha_beta_eta_guessing, "should have equal row names, in same order")
  expect_equal(error_message_beta_na_pattern$errors$beta, "can only contain NA at the end of rows, no values allowed after an NA in a row")
  expect_equal(error_message_eta_na_pattern$errors$eta, "can only contain NA at the end of rows, no values allowed after an NA in a row")
  expect_equal(error_message_incorrect_dim_estimate$errors$estimate, "length should be equal to the number of columns of the alpha matrix")
  expect_equal(error_message_incorrect_dim_variance$errors$variance, "should have a length equal to the length of estimate squared")
  expect_equal(error_message_variance_not_positive_definite$errors$variance, "matrix is not positive definite") 
  expect_equal(error_message_model_unknown$errors$model, "of unknown type")
  expect_equal(error_message_beta$errors$beta, "is missing")
  expect_equal(error_message_beta_theta$errors$beta_and_eta, "are both missing; define at least one of them")
  expect_equal(error_message_beta_theta_mismatch$errors$beta_and_eta, "objects do not match")
  expect_equal(error_message_prior_form$errors$prior_form, "is missing")
  expect_equal(error_message_prior_parameters$errors$prior_parameters, "is missing")
  expect_equal(error_message_prior_form_type$errors$prior_form, "of unknown type")
  expect_equal(error_message_prior_form_par_mismatch_uniform$errors$prior_form_is_uniform, "so prior_parameters should contain lower_bound and upper_bound")
  expect_equal(error_message_prior_form_par_mismatch_normal$errors$prior_form_is_normal, "so prior_parameters should contain mu and Sigma")
  expect_equal(error_message_prior_parameters_mu$errors$prior_parameters_mu, "should have same length as estimate")
  expect_equal(error_message_prior_parameters_sigma1$errors$prior_parameters_sigma, "should be a square positive definite matrix, with dimensions equal to the length of estimate")
  expect_equal(error_message_prior_parameters_sigma2$errors$prior_parameters_sigma, "should be a square positive definite matrix, with dimensions equal to the length of estimate")
  expect_equal(error_message_prior_parameters_sigma3$errors$prior_parameters_sigma, "should be a square positive definite matrix, with dimensions equal to the length of estimate")
  expect_equal(error_message_prior_parameters_bounds$errors$prior_parameters_bounds, "should contain lower and upper bound of the same length as estimate")
  expect_equal(error_message_max_n$errors$stop_test, "contains no max_n")
  expect_equal(error_message_max_n_too_large$errors$stop_test_max_n, "is larger than the number of items in the item bank")
  expect_equal(error_message_invalid_cutoffs$errors$stop_test_cutoffs, "should be a matrix without missing values, and number of rows equal to max_n and number of columns equal to the number of dimensions")
  expect_equal(error_message_start_items_0$errors$start_items, "requires n > 0 for posterior expected kullback leibler information summary")
  expect_equal(error_message_start_items_random_by_dim_scalar$errors$start_items_n, "contains inconsistent information. Total length of start phase and sum of length per dimension do not match")
  expect_equal(error_message_start_items_random_by_dim_vector$errors$start_items_n, "contains inconsistent information. Total length of start phase and sum of length per dimension do not match (n != sum(n_by_dimension)")
  expect_equal(error_message_cutoffs$errors$stop_test, "contains cutoff values in non-matrix format")
  expect_equal(error_message_answers$errors$answers, "contains non-existing key")
  expect_equal(error_message_estimator_unknown$errors$estimator, "of unknown type")
  expect_equal(error_message_information_summary_unknown$errors$information_summary, "of unknown type")
  expect_equal(error_message_estimator_ml_posterior_mixed$errors$estimator_is_maximum_likelihood, "so using a posterior information summary makes no sense")
  expect_equal(error_message_bounds_should_be_null$errors$bounds, "can only be defined if estimator is maximum likelihood")
  expect_equal(error_message_only_characteristics$errors$constraints_and_characts, "constraints and characteristics should either be defined both or not at all")  
  expect_equal(error_message_constraints_structure$errors$constraints_structure, "should be a list of length three lists, with elements named 'name', 'op', 'target'")  
  expect_equal(error_message_wrong_characts_and_constraints$errors$characteristics, "should be a data frame with number of rows equal to the number of items in the item bank")
  expect_equal(error_message_wrong_characts_and_constraints$errors$constraints_name_elements, "should be defined as described in the details section of constraints_lp_format()")
  expect_equal(error_message_wrong_characts_and_constraints$errors$constraints_operator_elements, "should be defined as described in the details section of constraints_lp_format()")
  expect_equal(error_message_wrong_characts_and_constraints$errors$constraints_target_elements, "should be defined as described in the details section of constraints_lp_format()")  
 })



# these simulations take a long time to run, if (FALSE) ensures that they are not run each time the tests are run
if (FALSE) {
  context("simulations")
  
  test_that("Simulations for context Jan Bebber", {
    replications_per_unique_condition <- 100
    true_theta_vec <- c(-2, 1, 3)
    number_items_vec <- 100
    number_answer_categories_vec <- 2
    number_dimensions <- 1
    
    start_items <- list(type = 'random', n = 0)
    variance_target <- .45^2
    model_vec <- "GPCM"
    estimator_vec <- "expected_aposteriori"
    information_summary_vec <- c("determinant", "posterior_determinant", "trace", "posterior_trace")
    prior_form <- "normal"
    prior_parameters <- list(mu = 0, Sigma = diag(1))
    min_n <- 4
    max_n <- 12
    
    estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, start_items, variance_target, replications_per_unique_condition, number_dimensions, prior_form = prior_form, prior_parameters = prior_parameters, min_n = min_n, max_n = max_n)
    
    conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, replications_per_unique_condition, number_dimensions)
    condition_vector <- sort(rep(1:(ncol(estimates_and_variance)/replications_per_unique_condition), replications_per_unique_condition))
    
    estimates_and_conditions <- cbind(t(estimates_and_variance), conditions[, c("true_theta", "replication", "number_items", "number_answer_categories", "model", "estimator", "information_summary")], condition_vector)
    
    average_per_condition_true_minus2 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == -2), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == -2)]), "mean")
    sd_per_condition_true_minus2 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == -2), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == -2)]), "sd")
    
    average_per_condition_true_1 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 1), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 1)]), "mean")
    sd_per_condition_true_1 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 1), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 1)]), "sd") 
    
    average_per_condition_true_3 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 3), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 3)]), "mean")
    sd_per_condition_true_3 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 3), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 3)]), "sd") 
      
    expect_equal(round(average_per_condition_true_minus2$x, 3), c(-1.533, -1.574, -1.616, -1.573)) 
    # prior sd = 3 and 1000 replications: -2.144 -2.176 -2.132 -2.143; small bias is due to left skewed posterior
    # prior sd = 3 and item bank 500 and 1000 replications: -2.059 -2.077 -2.098 -2.079
    # prior sd = 3 and map estimator: -1.957 -2.147 -2.090 -2.043
    # prior sd = 3 and larger number of items: -2.067 -2.091 -2.085 -2.048
    expect_equal(round(average_per_condition_true_1$x, 3), c(.803, .780, .892, .831))
    # prior sd = 3: 0.941 1.082 1.010 1.043
    # prior sd = 3 and item bank 500: 0.994 1.013 0.959 0.961
    # prior sd = 3 and map estimator: 0.995 0.932 1.034 1.042
    # prior sd = 3 and larger number of items: 0.996 0.973 1.024 0.985
    expect_equal(round(average_per_condition_true_3$x, 3), c(2.109, 2.057, 2.133, 2.067))
    # prior sd = 3: 3.105 3.188 3.270 3.061
    # prior sd = 3 and item bank 500: 2.828 2.813 2.788 2.793
    # prior sd = 3 and map estimator: 2.954 2.962 2.941 2.897
    
    # Observed standard deviation of the estimates per condition
    expect_equal(round(sd_per_condition_true_minus2$x, 3), c(.357, .361, .348, .385))
    # prior sd = 3: 0.668 0.602 0.615 0.754
    expect_equal(round(sd_per_condition_true_1$x, 3), c(.401, .406, .394, .423))
    # prior sd = 3: 0.463 0.561 0.539 0.521
    expect_equal(round(sd_per_condition_true_3$x, 3), c(.326, .358, .289, .319))
    # prior sd = 3: 0.985 0.852 0.815 0.780
    
    # Five number summary of reported posterior sd
    expect_equal(round(sqrt(fivenum(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == -2), "variance_estimate"])), 3), c(.424, .450, .470, .496, .578))
    # prior sd = 3: 0.460 0.555 0.616 0.696 1.172
    expect_equal(round(sqrt(fivenum(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 1), "variance_estimate"])), 3), c(.426, .438, .443, .448, .513))
    # prior sd = 3: 0.432 0.464 0.484 0.514 0.856
    expect_equal(round(sqrt(fivenum(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 3), "variance_estimate"])), 3), c(.434, .497, .524, .553, .607))
    # prior sd = 3: 0.472 0.700 0.807 0.977 1.660
  })
    
  test_that("one dimension, one replication per condition, expected_aposteriori, estimation procedure is gauss hermite quadrature", {
    replications_per_unique_condition <- 1
    true_theta_vec <- c(-2, 1, 3)
    number_items_vec <- c(50, 200)
    number_answer_categories_vec <- c(2, 4)
    number_dimensions <- 1
    
    start_items <- list(type = 'random', n = 3)
    variance_target <- .1^2
    model_vec <- c("3PLM", "GRM", "GPCM", "SM")
    estimator_vec <- "expected_aposteriori"
    information_summary_vec <- c("determinant", "posterior_determinant", "trace", "posterior_trace", "posterior_expected_kullback_leibler") 
    prior_form <- "normal"
    prior_parameters <- list(mu = 0, Sigma = matrix(9))
    
    estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, start_items, variance_target, replications_per_unique_condition, number_dimensions, prior_form = prior_form, prior_parameters = prior_parameters, eap_estimation_procedure = "gauss_hermite_quad")
    
    conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, replications_per_unique_condition, number_dimensions)
    
    estimates_and_conditions <- cbind(t(estimates_and_variance), conditions[, c("true_theta", "number_items", "number_answer_categories", "model", "estimator", "information_summary")])
    
    average_per_true_theta_and_number_items <- aggregate(estimates_and_conditions[,"estimated_theta"], list(estimates_and_conditions[,"true_theta"], estimates_and_conditions[,"number_items"]), "mean")
    sd_per_true_theta_and_number_items <- aggregate(estimates_and_conditions[,"estimated_theta"], list(estimates_and_conditions[,"true_theta"], estimates_and_conditions[,"number_items"]), "sd") 
    
    expect_equal(round(average_per_true_theta_and_number_items[,"x"], 3), c(-2.013, .998, 3.167, -1.985, 1.039, 3.079))
    expect_equal(round(sd_per_true_theta_and_number_items[,"x"], 3), c(.536, .388, .624, .205, .189, .382))
    expect_equal(round(sqrt(fivenum(estimates_and_conditions[,"variance_estimate"])), 3), c(.109, .184, .285, .379, 1.220))    
  })
  
  test_that("one dimension, one replication per condition, expected_aposteriori, estimation procedure is riemansumm", {
    replications_per_unique_condition <- 1
    true_theta_vec <- c(-2, 1, 3)
    number_items_vec <- c(50, 200)
    number_answer_categories_vec <- c(2, 4)
    number_dimensions <- 1
    
    start_items <- list(type = 'random', n = 3)
    variance_target <- .1^2
    model_vec <- c("3PLM", "GRM", "GPCM", "SM")
    estimator_vec <- "expected_aposteriori"
    information_summary_vec <- c("determinant", "posterior_determinant", "trace", "posterior_trace", "posterior_expected_kullback_leibler") 
    prior_form <- "normal"
    prior_parameters <- list(mu = 0, Sigma = matrix(9))
    
    estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, start_items, variance_target, replications_per_unique_condition, number_dimensions, prior_form = prior_form, prior_parameters = prior_parameters, eap_estimation_procedure = "riemannsum")
    
    conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, replications_per_unique_condition, number_dimensions)
    
    estimates_and_conditions <- cbind(t(estimates_and_variance), conditions[, c("true_theta", "number_items", "number_answer_categories", "model", "estimator", "information_summary")])
    
    average_per_true_theta_and_number_items <- aggregate(estimates_and_conditions[,"estimated_theta"], list(estimates_and_conditions[,"true_theta"], estimates_and_conditions[,"number_items"]), "mean")
    sd_per_true_theta_and_number_items <- aggregate(estimates_and_conditions[,"estimated_theta"], list(estimates_and_conditions[,"true_theta"], estimates_and_conditions[,"number_items"]), "sd") 
    
    expect_equal(round(average_per_true_theta_and_number_items[,"x"], 3), c(-2.013, .998, 3.151, -1.985, 1.039, 3.079))
    expect_equal(round(sd_per_true_theta_and_number_items[,"x"], 3), c(.536, .388, .635, .205, .189, .382))
    expect_equal(round(sqrt(fivenum(estimates_and_conditions[,"variance_estimate"])), 3), c(.109, .184, .285, .379, 1.219))    
  })
  
  test_that("one dimension, one replication per condition, maximum_aposteriori", {
    replications_per_unique_condition <- 1
    true_theta_vec <- c(-2, 1, 3)
    number_items_vec <- c(50, 200)
    number_answer_categories_vec <- c(2, 4)
    number_dimensions <- 1
    
    start_items <- list(type = 'random', n = 3)
    variance_target <- .1^2
    model_vec <- c("3PLM","GRM","GPCM","SM")
    estimator_vec <- "maximum_aposteriori"
    information_summary_vec <- c("determinant", "posterior_determinant", "trace", "posterior_trace", "posterior_expected_kullback_leibler")  
    prior_form <- "normal"
    prior_parameters <- list(mu = 0, Sigma = matrix(9))
    
    estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, start_items, variance_target, replications_per_unique_condition, number_dimensions, prior_form = prior_form, prior_parameters = prior_parameters, eap_estimation_procedure = "gauss_hermite_quad")
    
    conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, replications_per_unique_condition, number_dimensions)
    
    estimates_and_conditions <- cbind(t(estimates_and_variance), conditions[, c("true_theta", "number_items", "number_answer_categories", "model", "estimator", "information_summary")])
    
    average_per_true_theta_and_number_items <- aggregate(estimates_and_conditions[,"estimated_theta"], list(estimates_and_conditions[,"true_theta"], estimates_and_conditions[,"number_items"]), "mean")
    sd_per_true_theta_and_number_items <- aggregate(estimates_and_conditions[,"estimated_theta"], list(estimates_and_conditions[,"true_theta"], estimates_and_conditions[,"number_items"]), "sd") 
    
    expect_equal(round(average_per_true_theta_and_number_items[,"x"], 3), c(-1.988, .974, 3.035, -1.972, 1.026, 3.006))
    expect_equal(round(sd_per_true_theta_and_number_items[,"x"], 3), c(.530, .352, .525, .213, .166, .317))
    expect_equal(round(sqrt(fivenum(estimates_and_conditions[,"variance_estimate"])), 3), c(.110, .185, .281, .372, 1.004))
  })
  
  test_that("one dimension, one replication per condition, maximum_likelihood", {
    replications_per_unique_condition <- 1
    true_theta_vec <- c(-2, 1, 3)
    number_items_vec <- c(50, 200)
    number_answer_categories_vec <- c(2, 4)
    number_dimensions <- 1
    lowerbound <- -6
    upperbound <- 6
    
    start_items <- list(type = 'random', n = 3)
    variance_target <- .1^2
    model_vec <- c("3PLM","GRM","GPCM","SM")
    estimator_vec <- "maximum_likelihood"
    information_summary_vec <- c("determinant", "trace")
    
    estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, start_items, variance_target, replications_per_unique_condition, number_dimensions, lower_bound = lowerbound, upper_bound = upperbound)
    
    conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, replications_per_unique_condition, number_dimensions)
    
    estimates_and_conditions <- cbind(t(estimates_and_variance), conditions[, c("true_theta", "number_items", "number_answer_categories", "model", "estimator", "information_summary")])
    
    average_per_true_theta_and_number_items <- aggregate(estimates_and_conditions[,"estimated_theta"], list(estimates_and_conditions[,"true_theta"], estimates_and_conditions[,"number_items"]), "mean")
    sd_per_true_theta_and_number_items <- aggregate(estimates_and_conditions[,"estimated_theta"], list(estimates_and_conditions[,"true_theta"], estimates_and_conditions[,"number_items"]), "sd") 
    
    expect_equal(round(average_per_true_theta_and_number_items[,"x"], 3), c(-2.030, 1.079, 3.252, -1.960, 1.048, 3.161))
    expect_equal(round(sd_per_true_theta_and_number_items[,"x"], 3), c(.419, .451, .843, .219, .137, .516))
    expect_equal(round(sqrt(fivenum(estimates_and_conditions[,"variance_estimate"])), 3), c(.111, .184, .282, .406, 1.742))
  })
  
  
  test_that("one dimension, estimator maximum_likelihood, 100 replications per condition", {
    #run again
    replications_per_unique_condition <- 100
    true_theta_vec <- c(-2, 1, 3)
    number_items_vec <- c(50, 200)
    number_answer_categories_vec <- c(2, 4)
    number_dimensions <- 1
    lowerbound <- -6
    upperbound <- 6
    
    start_items <- list(type = 'random', n = 3)
    variance_target <- .1^2
    model_vec <- c("3PLM","GRM","GPCM","SM")
    estimator_vec <- "maximum_likelihood"
    information_summary_vec <- c("determinant", "posterior_determinant", "trace", "posterior_trace")  # "posterior_expected_kullback_leibler" 
    
    estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, start_items, variance_target, replications_per_unique_condition, number_dimensions, lowerbound = lowerbound, upperbound = upperbound)    
    
    conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, replications_per_unique_condition, number_dimensions)
    
    condition_vector <- sort(rep(1:(ncol(estimates_and_variance)/replications_per_unique_condition), replications_per_unique_condition))
    
    estimates_and_conditions <- cbind(t(estimates_and_variance)[, c("estimated_theta", "variance_estimate")], 
                                      conditions[, c("true_theta", "replication", "number_items", "number_answer_categories", "model", "estimator", "information_summary")], 
                                      condition_vector)
    
    average_thetaminus2_nitems50 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 50), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 50)]), "mean")
    average_theta1_nitems50 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 50), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 50)]), "mean")
    average_theta3_nitems50 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 50), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 50)]), "mean")
    
    sd_thetaminus2_nitems50 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 50), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 50)]), "sd")
    sd_theta1_nitems50 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 50), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 50)]), "sd")
    sd_theta3_nitems50 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 50), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 50)]), "sd")
    
    average_thetaminus2_nitems200 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 200), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 200)]), "mean")
    average_theta1_nitems200 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 200), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 200)]), "mean")
    average_theta3_nitems200 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 200), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 200)]), "mean")
    
    sd_thetaminus2_nitems200 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 200), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 200)]), "sd")
    sd_theta1_nitems200 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 200), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 200)]), "sd")
    sd_theta3_nitems200 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 200), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 200)]), "sd")
       
    # five number summary of average theta estimate per condition, with true theta is -2 and 50 items
    expect_equal(round(fivenum(average_thetaminus2_nitems50[,"x"]), 3), c(-2.159, -2.095, -2.048, -2.020, -1.979))
    # five number summary of average theta estimate per condition, with true theta is 1 and 50 items
    expect_equal(round(fivenum(average_theta1_nitems50[,"x"]), 3), c(.969, .995, 1.028, 1.063, 1.108))
    # five number summary of average theta estimate per condition, with true theta is 3 and 50 items
    expect_equal(round(fivenum(average_theta3_nitems50[,"x"]), 3), c(2.984, 3.032, 3.102, 3.203, 3.340))
    # five number summary of average theta estimate per condition, with true theta is -2 and 200 items
    expect_equal(round(fivenum(average_thetaminus2_nitems200[,"x"]), 3), c(-2.059, -2.031, -2.010, -2.001, -1.975))
    # five number summary of average theta estimate per condition, with true theta is 1 and 200 items
    expect_equal(round(fivenum(average_theta1_nitems200[,"x"]), 3), c(.958, .991, 1.003, 1.015, 1.037))
    # five number summary of average theta estimate per condition, with true theta is 3 and 200 items
    expect_equal(round(fivenum(average_theta3_nitems200[,"x"]), 3), c(2.980, 3.004, 3.019, 3.031, 3.094))
    
    # five number summary of observed sd of the theta estimates within each condition, with true theta is -2 en 50 items
    expect_equal(round(fivenum(sd_thetaminus2_nitems50[,"x"]), 3), c(.238, .291, .454, .512, .617))
    # five number summary of observed sd of the theta estimates within each condition, with true theta is 1 en 50 items
    expect_equal(round(fivenum(sd_theta1_nitems50[,"x"]), 3), c(.228, .258, .350, .381, .452))
    # five number summary of observed sd of the theta estimates within each condition, with true theta is 3 en 50 items
    expect_equal(round(fivenum(sd_theta3_nitems50[,"x"]), 3), c(.304, .367, .694, .777, .901))
    # five number summary of observed sd of the theta estimates within each condition, with true theta is -2 en 200 items
    expect_equal(round(fivenum(sd_thetaminus2_nitems200[,"x"]), 3), c(.118, .145, .221, .239, .271))
    # five number summary of observed sd of the theta estimates within each condition, with true theta is 1 en 200 items
    expect_equal(round(fivenum(sd_theta1_nitems200[,"x"]), 3), c(.106, .121, .169, .183, .201))
    # five number summary of observed sd of the theta estimates within each condition, with true theta is 3 en 200 items
    expect_equal(round(fivenum(sd_theta3_nitems200[,"x"]), 3), c(.136, .173, .311, .339, .369))
        
    # five number summary of reported sd of the theta estimate within each condition where number of items is 50 and 200, respectively
    expect_equal(round(sqrt(fivenum(estimates_and_conditions[which(estimates_and_conditions[,"number_items"] == 50), "variance_estimate"])), 3), c(.198, .312, .376, .506, 18991.320))
    expect_equal(round(sqrt(fivenum(estimates_and_conditions[which(estimates_and_conditions[,"number_items"] == 200), "variance_estimate"])), 3), c(.105, .156, .183, .246, .571))
  })
  
test_that("one dimension, estimator maximum_aposteriori, 100 replications per condition, prior is normal(0, 9)", {
  #run again
  replications_per_unique_condition <- 100
  true_theta_vec <- c(-2, 1, 3)
  number_items_vec <- c(50, 200)
  number_answer_categories_vec <- c(2, 4)
  number_dimensions <- 1
  
  start_items <- list(type = 'random', n = 3)
  variance_target <- .1^2
  model_vec <- c("3PLM","GRM","GPCM","SM")
  estimator_vec <- "maximum_aposteriori"
  information_summary_vec <- c("determinant", "posterior_determinant", "trace", "posterior_trace", "posterior_expected_kullback_leibler")
  lowerbound <- -20
  upperbound <- 20
  prior <- diag(1) * 9
  
  estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, start_items, variance_target, replications_per_unique_condition, number_dimensions, lowerbound = lowerbound, upperbound = upperbound, prior = prior)
  
  conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, replications_per_unique_condition, number_dimensions)
  
  condition_vector <- sort(rep(1:(ncol(estimates_and_variance)/replications_per_unique_condition), replications_per_unique_condition))
  
  estimates_and_conditions <- cbind(t(estimates_and_variance)[, c("estimated_theta", "variance_estimate")], 
                                    conditions[, c("true_theta", "replication", "number_items", "number_answer_categories", "model", "estimator", "information_summary")], 
                                    condition_vector)
  
  average_thetaminus2_nitems50 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 50), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 50)]), "mean")
  average_theta1_nitems50 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 50), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 50)]), "mean")
  average_theta3_nitems50 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 50), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 50)]), "mean")
  
  sd_thetaminus2_nitems50 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 50), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 50)]), "sd")
  sd_theta1_nitems50 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 50), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 50)]), "sd")
  sd_theta3_nitems50 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 50), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 50)]), "sd")
  
  average_thetaminus2_nitems200 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 200), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 200)]), "mean")
  average_theta1_nitems200 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 200), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 200)]), "mean")
  average_theta3_nitems200 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 200), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 200)]), "mean")
  
  sd_thetaminus2_nitems200 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 200), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 200)]), "sd")
  sd_theta1_nitems200 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 200), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 200)]), "sd")
  sd_theta3_nitems200 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 200), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 200)]), "sd")
  
  
  # five number summary of average theta estimate per condition, with true theta is -2 and 50 items
  expect_equal(round(fivenum(average_thetaminus2_nitems50[,"x"]), 3), c(-2.116, -2.042, -2.007, -1.988, -1.916))
  # five number summary of average theta estimate per condition, with true theta is 1 and 50 items
  expect_equal(round(fivenum(average_theta1_nitems50[,"x"]), 3), c(.952, .984, 1.008, 1.051, 1.086))
  # five number summary of average theta estimate per condition, with true theta is 3 and 50 items
  expect_equal(round(fivenum(average_theta3_nitems50[,"x"]), 3), c(2.881, 2.954, 2.982, 3.017, 3.098))
  # five number summary of average theta estimate per condition, with true theta is -2 and 200 items
  expect_equal(round(fivenum(average_thetaminus2_nitems200[,"x"]), 3), c(-2.040, -2.016, -1.995, -1.988, -1.962))
  # five number summary of average theta estimate per condition, with true theta is 1 and 200 items
  expect_equal(round(fivenum(average_theta1_nitems200[,"x"]), 3), c(.954, .986, .998, 1.006, 1.032))
  # five number summary of average theta estimate per condition, with true theta is 3 and 200 items
  expect_equal(round(fivenum(average_theta3_nitems200[,"x"]), 3), c(2.947, 2.981, 2.991, 3.006, 3.061))
  
  # five number summary of observed sd of the theta estimates within each condition, with true theta is -2 en 50 items
  expect_equal(round(fivenum(sd_thetaminus2_nitems50[,"x"]), 3), c(.239, .286, .439, .477, .603))
  # five number summary of observed sd of the theta estimates within each condition, with true theta is 1 en 50 items
  expect_equal(round(fivenum(sd_theta1_nitems50[,"x"]), 3), c(.220, .252, .344, .365, .441))
  # five number summary of observed sd of the theta estimates within each condition, with true theta is 3 en 50 items
  expect_equal(round(fivenum(sd_theta3_nitems50[,"x"]), 3), c(.279, .345, .576, .636, .756))
  # five number summary of observed sd of the theta estimates within each condition, with true theta is -2 en 200 items
  expect_equal(round(fivenum(sd_thetaminus2_nitems200[,"x"]), 3), c(.120, .142, .213, .227, .265))
  # five number summary of observed sd of the theta estimates within each condition, with true theta is 1 en 200 items
  expect_equal(round(fivenum(sd_theta1_nitems200[,"x"]), 3), c(.104, .129, .171, .182, .202))
  # five number summary of observed sd of the theta estimates within each condition, with true theta is 3 en 200 items
  expect_equal(round(fivenum(sd_theta3_nitems200[,"x"]), 3), c(.148, .169, .301, .324, .370))
  
  # five number summary of reported sd of the theta estimate within each condition where number of items is 50 and 200, respectively
  expect_equal(round(sqrt(fivenum(estimates_and_conditions[which(estimates_and_conditions[,"number_items"] == 50), "variance_estimate"])), 3), c(.197, .308, .370, .487, 1.497))
  expect_equal(round(sqrt(fivenum(estimates_and_conditions[which(estimates_and_conditions[,"number_items"] == 200), "variance_estimate"])), 3), c(.105, .156, .183, .243, .531))
})

test_that("one dimension, estimator maximum_aposteriori, 100 replications per condition, prior is uniform(-4, 4)", {
  replications_per_unique_condition <- 100
  true_theta_vec <- c(-2, 1, 3)
  number_items_vec <- c(50, 200)
  number_answer_categories_vec <- c(2, 4)
  number_dimensions <- 1
  
  start_items <- list(type = 'random', n = 3)
  variance_target <- .1^2
  model_vec <- c("3PLM","GRM","GPCM","SM")
  estimator_vec <- "maximum_aposteriori"
  information_summary_vec <- c("determinant", "posterior_determinant", "trace", "posterior_trace", "posterior_expected_kullback_leibler")
  prior_form <- "uniform"
  prior_parameters <- list(lower_bound = -4, upper_bound = 4)
  
  estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, start_items, variance_target, replications_per_unique_condition, number_dimensions, prior_form = prior_form, prior_parameters = prior_parameters)
  
  conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, replications_per_unique_condition, number_dimensions)
  
  condition_vector <- sort(rep(1:(ncol(estimates_and_variance)/replications_per_unique_condition), replications_per_unique_condition))
  
  estimates_and_conditions <- cbind(t(estimates_and_variance)[, c("estimated_theta", "variance_estimate")], 
                                    conditions[, c("true_theta", "replication", "number_items", "number_answer_categories", "model", "estimator", "information_summary")], 
                                    condition_vector)
  
  average_thetaminus2_nitems50 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 50), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 50)]), "mean")
  average_theta1_nitems50 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 50), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 50)]), "mean")
  average_theta3_nitems50 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 50), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 50)]), "mean")
  
  sd_thetaminus2_nitems50 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 50), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 50)]), "sd")
  sd_theta1_nitems50 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 50), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 50)]), "sd")
  sd_theta3_nitems50 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 50), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 50)]), "sd")
  
  average_thetaminus2_nitems200 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 200), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 200)]), "mean")
  average_theta1_nitems200 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 200), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 200)]), "mean")
  average_theta3_nitems200 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 200), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 200)]), "mean")
  
  sd_thetaminus2_nitems200 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 200), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 200)]), "sd")
  sd_theta1_nitems200 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 200), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 200)]), "sd")
  sd_theta3_nitems200 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 200), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 200)]), "sd")
  
  
  # five number summary of average theta estimate per condition, with true theta is -2 and 50 items
  expect_equal(round(fivenum(average_thetaminus2_nitems50[,"x"]), 3), c(-2.201, -2.083, -2.045, -2.018, -1.956))
  # five number summary of average theta estimate per condition, with true theta is 1 and 50 items
  expect_equal(round(fivenum(average_theta1_nitems50[,"x"]), 3), c(.972, .996, 1.030, 1.060, 1.108))
  # five number summary of average theta estimate per condition, with true theta is 3 and 50 items
  expect_equal(round(fivenum(average_theta3_nitems50[,"x"]), 3), c(2.942 ,3.015, 3.047, 3.105, 3.212))
  # five number summary of average theta estimate per condition, with true theta is -2 and 200 items
  expect_equal(round(fivenum(average_thetaminus2_nitems200[,"x"]), 3), c(-2.055, -2.030, -2.006, -1.992, -1.969))
  # five number summary of average theta estimate per condition, with true theta is 1 and 200 items
  expect_equal(round(fivenum(average_theta1_nitems200[,"x"]), 3), c(.964, .990, 1.002, 1.013, 1.035))
  # five number summary of average theta estimate per condition, with true theta is 3 and 200 items
  expect_equal(round(fivenum(average_theta3_nitems200[,"x"]), 3), c(2.964, 2.992, 3.012, 3.038, 3.101))
  
  # five number summary of observed sd of the theta estimates within each condition, with true theta is -2 en 50 items
  expect_equal(round(fivenum(sd_thetaminus2_nitems50[,"x"]), 3), c(.234, .288, .460, .501, .604))
  # five number summary of observed sd of the theta estimates within each condition, with true theta is 1 en 50 items
  expect_equal(round(fivenum(sd_theta1_nitems50[,"x"]), 3), c(.228, .259, .353, .383, .454))
  # five number summary of observed sd of the theta estimates within each condition, with true theta is 3 en 50 items
  expect_equal(round(fivenum(sd_theta3_nitems50[,"x"]), 3), c(.301, .360, .575, .601, .659))
  # five number summary of observed sd of the theta estimates within each condition, with true theta is -2 en 200 items
  expect_equal(round(fivenum(sd_thetaminus2_nitems200[,"x"]), 3), c(.121, .142, .215, .240, .261))
  # five number summary of observed sd of the theta estimates within each condition, with true theta is 1 en 200 items
  expect_equal(round(fivenum(sd_theta1_nitems200[,"x"]), 3), c(.106, .130, .166, .183, .199))
  # five number summary of observed sd of the theta estimates within each condition, with true theta is 3 en 200 items
  expect_equal(round(fivenum(sd_theta3_nitems200[,"x"]), 3), c(.126, .171, .314, .334, .365))
  
  # five number summary of reported sd of the theta estimate within each condition where number of items is 50 and 200, respectively
  expect_equal(round(sqrt(fivenum(estimates_and_conditions[which(estimates_and_conditions[,"number_items"] == 50), "variance_estimate"])), 3), c(.198, .312, .376, .505, 1.138))
  expect_equal(round(sqrt(fivenum(estimates_and_conditions[which(estimates_and_conditions[,"number_items"] == 200), "variance_estimate"])), 3), c(.105, .156, .183, .245, .485))
})


test_that("one dimension, estimator expected_aposteriori, 100 replications per condition, prior is normal(0, 9), estimation procedure is gauss hermite quadrature", {
  #run again
  replications_per_unique_condition <- 100
  true_theta_vec <- c(-2, 1, 3)
  number_items_vec <- c(50, 200)
  number_answer_categories_vec <- c(2, 4)
  number_dimensions <- 1
  
  start_items <- list(type = 'random', n = 3)
  variance_target <- .1^2
  model_vec <- c("3PLM","GRM","GPCM","SM")
  estimator_vec <- "expected_aposteriori"
  information_summary_vec <- c("determinant", "posterior_determinant", "trace", "posterior_trace", "posterior_expected_kullback_leibler")
  lowerbound <- -20
  upperbound <- 20
  prior <- diag(1) * 9  
  estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, start_items, variance_target, replications_per_unique_condition, number_dimensions, lowerbound = lowerbound, upperbound = upperbound, prior = prior, eap_estimation_procedure = "gauss_hermite_quad")
  
  conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, replications_per_unique_condition, number_dimensions)
  
  condition_vector <- sort(rep(1:(ncol(estimates_and_variance)/replications_per_unique_condition), replications_per_unique_condition))
  
  estimates_and_conditions <- cbind(t(estimates_and_variance)[, c("estimated_theta", "variance_estimate")], 
                                    conditions[, c("true_theta", "replication", "number_items", "number_answer_categories", "model", "estimator", "information_summary")], 
                                    condition_vector)
  
  average_thetaminus2_nitems50 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 50), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 50)]), "mean")
  average_theta1_nitems50 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 50), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 50)]), "mean")
  average_theta3_nitems50 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 50), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 50)]), "mean")
  
  sd_thetaminus2_nitems50 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 50), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 50)]), "sd")
  sd_theta1_nitems50 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 50), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 50)]), "sd")
  sd_theta3_nitems50 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 50), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 50)]), "sd")
  
  average_thetaminus2_nitems200 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 200), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 200)]), "mean")
  average_theta1_nitems200 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 200), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 200)]), "mean")
  average_theta3_nitems200 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 200), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 200)]), "mean")
  
  sd_thetaminus2_nitems200 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 200), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 200)]), "sd")
  sd_theta1_nitems200 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 200), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 200)]), "sd")
  sd_theta3_nitems200 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 200), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 200)]), "sd")
  
  # five number summary of average theta estimate per condition, with true theta is -2 and 50 items
  expect_equal(round(fivenum(average_thetaminus2_nitems50[,"x"]), 3), c(-2.168, -2.096, -2.049, -2.008, -1.947))
  # five number summary of average theta estimate per condition, with true theta is 1 and 50 items
  expect_equal(round(fivenum(average_theta1_nitems50[,"x"]), 3), c(.960, .994, 1.029, 1.059, 1.127))
  # five number summary of average theta estimate per condition, with true theta is 3 and 50 items
  expect_equal(round(fivenum(average_theta3_nitems50[,"x"]), 3), c(2.959, 3.013, 3.086, 3.148, 3.247))
  # five number summary of average theta estimate per condition, with true theta is -2 and 200 items
  expect_equal(round(fivenum(average_thetaminus2_nitems200[,"x"]), 3), c(-2.059, -2.030, -2.004, -1.992, -1.976))
  # five number summary of average theta estimate per condition, with true theta is 1 and 200 items
  expect_equal(round(fivenum(average_theta1_nitems200[,"x"]), 3), c(.972, .990, 1.002, 1.016, 1.039))
  # five number summary of average theta estimate per condition, with true theta is 3 and 200 items
  expect_equal(round(fivenum(average_theta3_nitems200[,"x"]), 3), c(2.976, 2.998, 3.019, 3.038, 3.122))
  
  # five number summary of observed sd of the theta estimates within each condition, with true theta is -2 en 50 items
  expect_equal(round(fivenum(sd_thetaminus2_nitems50[,"x"]), 3), c(.242, .304, .427, .507, .569))
  # five number summary of observed sd of the theta estimates within each condition, with true theta is 1 en 50 items
  expect_equal(round(fivenum(sd_theta1_nitems50[,"x"]), 3), c(.228, .256, .357, .371, .460))
  # five number summary of observed sd of the theta estimates within each condition, with true theta is 3 en 50 items
  expect_equal(round(fivenum(sd_theta3_nitems50[,"x"]), 3), c(.301, .358, .629, .683, .800))
  # five number summary of observed sd of the theta estimates within each condition, with true theta is -2 en 200 items
  expect_equal(round(fivenum(sd_thetaminus2_nitems200[,"x"]), 3), c(.121, .142, .220, .239, .274))
  # five number summary of observed sd of the theta estimates within each condition, with true theta is 1 en 200 items
  expect_equal(round(fivenum(sd_theta1_nitems200[,"x"]), 3), c(.105, .132, .170, .178, .199))
  # five number summary of observed sd of the theta estimates within each condition, with true theta is 3 en 200 items
  expect_equal(round(fivenum(sd_theta3_nitems200[,"x"]), 3), c(.138, .171, .301, .330, .359))
  
  # five number summary of reported sd of the theta estimate within each condition where number of items is 50 and 200, respectively
  expect_equal(round(sqrt(fivenum(estimates_and_conditions[which(estimates_and_conditions[,"number_items"] == 50), "variance_estimate"])), 3), c(.198, .311, .376, .500, 1.557))
  expect_equal(round(sqrt(fivenum(estimates_and_conditions[which(estimates_and_conditions[,"number_items"] == 200), "variance_estimate"])), 3), c(.105, .156, .183, .246, .543))
})

test_that("one dimension, estimator expected_aposteriori, 100 replications per condition, prior is normal(0, 9), estimation procedure is riemannsum", {
  #run again
  replications_per_unique_condition <- 100
  true_theta_vec <- c(-2, 1, 3)
  number_items_vec <- c(50, 200)
  number_answer_categories_vec <- c(2, 4)
  number_dimensions <- 1
  
  start_items <- list(type = 'random', n = 3)
  variance_target <- .1^2
  model_vec <- c("3PLM","GRM","GPCM","SM")
  estimator_vec <- "expected_aposteriori"
  information_summary_vec <- c("determinant", "posterior_determinant", "trace", "posterior_trace", "posterior_expected_kullback_leibler")
  lowerbound <- -20
  upperbound <- 20
  prior <- diag(1) * 9  
  estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, start_items, variance_target, replications_per_unique_condition, number_dimensions, lowerbound = lowerbound, upperbound = upperbound, prior = prior, eap_estimation_procedure = "riemannsum")
  
  conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, replications_per_unique_condition, number_dimensions)
  
  condition_vector <- sort(rep(1:(ncol(estimates_and_variance)/replications_per_unique_condition), replications_per_unique_condition))
  
  estimates_and_conditions <- cbind(t(estimates_and_variance)[, c("estimated_theta", "variance_estimate")], 
                                    conditions[, c("true_theta", "replication", "number_items", "number_answer_categories", "model", "estimator", "information_summary")], 
                                    condition_vector)
  
  average_thetaminus2_nitems50 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 50), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 50)]), "mean")
  average_theta1_nitems50 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 50), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 50)]), "mean")
  average_theta3_nitems50 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 50), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 50)]), "mean")
  
  sd_thetaminus2_nitems50 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 50), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 50)]), "sd")
  sd_theta1_nitems50 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 50), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 50)]), "sd")
  sd_theta3_nitems50 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 50), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 50)]), "sd")
  
  average_thetaminus2_nitems200 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 200), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 200)]), "mean")
  average_theta1_nitems200 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 200), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 200)]), "mean")
  average_theta3_nitems200 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 200), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 200)]), "mean")
  
  sd_thetaminus2_nitems200 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 200), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 200)]), "sd")
  sd_theta1_nitems200 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 200), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 200)]), "sd")
  sd_theta3_nitems200 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 200), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 200)]), "sd")
  
  # five number summary of average theta estimate per condition, with true theta is -2 and 50 items
  expect_equal(round(fivenum(average_thetaminus2_nitems50[,"x"]), 3), c(-2.170, -2.095, -2.048, -2.008, -1.947))
  # five number summary of average theta estimate per condition, with true theta is 1 and 50 items
  expect_equal(round(fivenum(average_theta1_nitems50[,"x"]), 3), c(.958, .993, 1.029, 1.059, 1.127))
  # five number summary of average theta estimate per condition, with true theta is 3 and 50 items
  expect_equal(round(fivenum(average_theta3_nitems50[,"x"]), 3), c(2.959, 3.013, 3.085, 3.148, 3.241))
  # five number summary of average theta estimate per condition, with true theta is -2 and 200 items
  expect_equal(round(fivenum(average_thetaminus2_nitems200[,"x"]), 3), c(-2.057, -2.030, -2.004, -1.992, -1.975))
  # five number summary of average theta estimate per condition, with true theta is 1 and 200 items
  expect_equal(round(fivenum(average_theta1_nitems200[,"x"]), 3), c(.972, .990, 1.003, 1.015, 1.039))
  # five number summary of average theta estimate per condition, with true theta is 3 and 200 items
  expect_equal(round(fivenum(average_theta3_nitems200[,"x"]), 3), c(2.969, 2.998, 3.019, 3.036, 3.119))
  
  # five number summary of observed sd of the theta estimates within each condition, with true theta is -2 en 50 items
  expect_equal(round(fivenum(sd_thetaminus2_nitems50[,"x"]), 3), c(.242, .305, .426, .506, .567))
  # five number summary of observed sd of the theta estimates within each condition, with true theta is 1 en 50 items
  expect_equal(round(fivenum(sd_theta1_nitems50[,"x"]), 3), c(.227, .256, .357, .371, .464))
  # five number summary of observed sd of the theta estimates within each condition, with true theta is 3 en 50 items
  expect_equal(round(fivenum(sd_theta3_nitems50[,"x"]), 3), c(.301, .358, .629, .687, .800))
  # five number summary of observed sd of the theta estimates within each condition, with true theta is -2 en 200 items
  expect_equal(round(fivenum(sd_thetaminus2_nitems200[,"x"]), 3), c(.121, .143, .219, .239, .272))
  # five number summary of observed sd of the theta estimates within each condition, with true theta is 1 en 200 items
  expect_equal(round(fivenum(sd_theta1_nitems200[,"x"]), 3), c(.105, .132, .170, .178, .199))
  # five number summary of observed sd of the theta estimates within each condition, with true theta is 3 en 200 items
  expect_equal(round(fivenum(sd_theta3_nitems200[,"x"]), 3), c(.138, .171, .298, .331, .359))
  
  # five number summary of reported sd of the theta estimate within each condition where number of items is 50 and 200, respectively
  expect_equal(round(sqrt(fivenum(estimates_and_conditions[which(estimates_and_conditions[,"number_items"] == 50), "variance_estimate"])), 3), c(.198, .311, .376, .500, 1.555))
  expect_equal(round(sqrt(fivenum(estimates_and_conditions[which(estimates_and_conditions[,"number_items"] == 200), "variance_estimate"])), 3), c(.105, .156, .183, .246, .543))
})

test_that("one dimension, estimator expected_aposteriori, 100 replications per condition, prior is uniform(-4, 4)", {
  replications_per_unique_condition <- 100
  true_theta_vec <- c(-2, 1, 3)
  number_items_vec <- c(50, 200)
  number_answer_categories_vec <- c(2, 4)
  number_dimensions <- 1
  
  start_items <- list(type = 'random', n = 3)
  variance_target <- .1^2
  model_vec <- c("3PLM","GRM","GPCM","SM")
  estimator_vec <- "expected_aposteriori"
  information_summary_vec <- c("determinant", "posterior_determinant", "trace", "posterior_trace", "posterior_expected_kullback_leibler")
  prior_form <- "uniform"
  prior_parameters <- list(lower_bound = -4, upper_bound = 4)
  
  estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, start_items, variance_target, replications_per_unique_condition, number_dimensions, prior_form = prior_form, prior_parameters = prior_parameters, eap_estimation_procedure = "riemannsum")
  
  conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, replications_per_unique_condition, number_dimensions)
  
  condition_vector <- sort(rep(1:(ncol(estimates_and_variance)/replications_per_unique_condition), replications_per_unique_condition))
  
  estimates_and_conditions <- cbind(t(estimates_and_variance)[, c("estimated_theta", "variance_estimate")], 
                                    conditions[, c("true_theta", "replication", "number_items", "number_answer_categories", "model", "estimator", "information_summary")], 
                                    condition_vector)
  
  average_thetaminus2_nitems50 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 50), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 50)]), "mean")
  average_theta1_nitems50 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 50), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 50)]), "mean")
  average_theta3_nitems50 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 50), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 50)]), "mean")
  
  sd_thetaminus2_nitems50 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 50), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 50)]), "sd")
  sd_theta1_nitems50 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 50), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 50)]), "sd")
  sd_theta3_nitems50 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 50), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 50)]), "sd")
  
  average_thetaminus2_nitems200 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 200), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 200)]), "mean")
  average_theta1_nitems200 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 200), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 200)]), "mean")
  average_theta3_nitems200 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 200), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 200)]), "mean")
  
  sd_thetaminus2_nitems200 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 200), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == -2 & estimates_and_conditions[,"number_items"] == 200)]), "sd")
  sd_theta1_nitems200 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 200), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 1 & estimates_and_conditions[,"number_items"] == 200)]), "sd")
  sd_theta3_nitems200 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 200), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 3 & estimates_and_conditions[,"number_items"] == 200)]), "sd")
  
  # five number summary of average theta estimate per condition, with true theta is -2 and 50 items
  expect_equal(round(fivenum(average_thetaminus2_nitems50[,"x"]), 3), c(-2.222, -2.135, -2.073, -2.038, -1.946))
  # five number summary of average theta estimate per condition, with true theta is 1 and 50 items
  expect_equal(round(fivenum(average_theta1_nitems50[,"x"]), 3), c(.974, 1.007, 1.041, 1.077, 1.136))
  # five number summary of average theta estimate per condition, with true theta is 3 and 50 items
  expect_equal(round(fivenum(average_theta3_nitems50[,"x"]), 3), c(2.898, 2.963, 2.987, 3.014, 3.052))
  # five number summary of average theta estimate per condition, with true theta is -2 and 200 items
  expect_equal(round(fivenum(average_thetaminus2_nitems200[,"x"]), 3), c(-2.067, -2.047, -2.017, -2.004, -1.981))
  # five number summary of average theta estimate per condition, with true theta is 1 and 200 items
  expect_equal(round(fivenum(average_theta1_nitems200[,"x"]), 3), c(.966, .993, 1.003, 1.021, 1.058))
  # five number summary of average theta estimate per condition, with true theta is 3 and 200 items
  expect_equal(round(fivenum(average_theta3_nitems200[,"x"]), 3), c(2.987, 3.007, 3.027, 3.044, 3.108))
  
  # five number summary of observed sd of the theta estimates within each condition, with true theta is -2 en 50 items
  expect_equal(round(fivenum(sd_thetaminus2_nitems50[,"x"]), 3), c(.237, .308, .451, .482, .567))
  # five number summary of observed sd of the theta estimates within each condition, with true theta is 1 en 50 items
  expect_equal(round(fivenum(sd_theta1_nitems50[,"x"]), 3), c(.228, .260, .354, .386, .449))
  # five number summary of observed sd of the theta estimates within each condition, with true theta is 3 en 50 items
  expect_equal(round(fivenum(sd_theta3_nitems50[,"x"]), 3), c(.272, .316, .378, .395, .451))
  # five number summary of observed sd of the theta estimates within each condition, with true theta is -2 en 200 items
  expect_equal(round(fivenum(sd_thetaminus2_nitems200[,"x"]), 3), c(.117, .142, .218, .242, .269))
  # five number summary of observed sd of the theta estimates within each condition, with true theta is 1 en 200 items
  expect_equal(round(fivenum(sd_theta1_nitems200[,"x"]), 3), c(.103, .130, .174, .182, .206))
  # five number summary of observed sd of the theta estimates within each condition, with true theta is 3 en 200 items
  expect_equal(round(fivenum(sd_theta3_nitems200[,"x"]), 3), c(.141, .176, .276, .293, .319))
  
  # five number summary of reported sd of the theta estimate within each condition where number of items is 50 and 200, respectively
  expect_equal(round(sqrt(fivenum(estimates_and_conditions[which(estimates_and_conditions[,"number_items"] == 50), "variance_estimate"])), 3), c(.198, .309, .370, .473, 0.618))
  expect_equal(round(sqrt(fivenum(estimates_and_conditions[which(estimates_and_conditions[,"number_items"] == 200), "variance_estimate"])), 3), c(.105, .157, .184, .247, .358))
})


test_that("three dimensions, maximum_likelihood, 100 replications per condition", {
  # To be run yet
  replications_per_unique_condition <- 100 
  true_theta_vec <- c(-2, 1, 3)
  number_items_vec <- c(300)
  number_answer_categories_vec <- c(2, 4)
  number_dimensions <- 3
  
  start_items <- list(type = 'random_by_dimension', n_by_dimension = 3, n = 9)
  variance_target <- .1^2
  model_vec <- c("3PLM","GRM","GPCM","SM")
  estimator_vec <- "maximum_likelihood"
  information_summary_vec <- c("determinant", "trace") 
  
  estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, start_items, variance_target, replications_per_unique_condition, number_dimensions)
  estimates_and_variance_without_errors <- sapply(estimates_and_variance, 
                                                  FUN = function(x) { 
                                                    if (is.null(x)) 
                                                      matrix(rep(NA, number_dimensions + number_dimensions^2), nrow = 1, dimnames = list(NULL, names(estimates_and_variance[[1]])))
                                                    else 
                                                      x } ) 
  
  conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, replications_per_unique_condition, number_dimensions)
  
  condition_vector <- sort(rep(1:(ncol(estimates_and_variance_without_errors)/replications_per_unique_condition), replications_per_unique_condition))
  
  estimates_and_conditions <- cbind(t(estimates_and_variance_without_errors)[, c("estimated_theta1", "estimated_theta2", "estimated_theta3", "variance_estimate1", "variance_estimate5", "variance_estimate9")], 
                                    conditions[, c("replication", "number_items", "number_answer_categories", "model", "estimator", "information_summary")], 
                                    condition_vector)
  
  average_per_condition_dim1 <- aggregate(estimates_and_conditions[, "estimated_theta1"], list(condition_vector), "mean", na.rm = TRUE)
  average_per_condition_dim2 <- aggregate(estimates_and_conditions[, "estimated_theta2"], list(condition_vector), "mean", na.rm = TRUE)
  average_per_condition_dim3 <- aggregate(estimates_and_conditions[, "estimated_theta3"], list(condition_vector), "mean", na.rm = TRUE)
  
  sd_per_condition_dim1 <- aggregate(estimates_and_conditions[, "estimated_theta1"], list(condition_vector), "sd", na.rm = TRUE)
  sd_per_condition_dim2 <- aggregate(estimates_and_conditions[, "estimated_theta2"], list(condition_vector), "sd", na.rm = TRUE)
  sd_per_condition_dim3 <- aggregate(estimates_and_conditions[, "estimated_theta3"], list(condition_vector), "sd", na.rm = TRUE)
  
  number_na_per_condition <- aggregate(estimates_and_conditions[, "estimated_theta1"], list(condition_vector), FUN = function(x) { sum(is.na(x)) })
  
  # five number summary of average theta estimate per condition, dimension 1 (true theta is -2)
  expect_equal(round(fivenum(average_per_condition_dim1[,"x"]), 3), c(-2.102, -2.043, -2.026, -2.013, -1.971))
  # five number summary of average theta estimate per condition, dimension 2 (true theta is 1)
  expect_equal(round(fivenum(average_per_condition_dim2[,"x"]), 3), c(.965, .993, 1.004, 1.016, 1.044))
  # five number summary of average theta estimate per condition, dimension 3 (true theta is 2)
  expect_equal(round(fivenum(average_per_condition_dim3[,"x"]), 3), c(1.954, 1.992, 2.013, 2.048, 2.090))
  
  # five number summary of observed sd of the theta estimates within each condition, for dimension 1, 2, and 3, respectively
  expect_equal(round(fivenum(sd_per_condition_dim1[,"x"]), 3), c(.182, .213, .326, .351, .387))
  expect_equal(round(fivenum(sd_per_condition_dim2[,"x"]), 3), c(.151, .190, .251, .265, .291))
  expect_equal(round(fivenum(sd_per_condition_dim3[,"x"]), 3), c(.166, .210, .297, .329, .370))
  
  # five number summary of reported sd of the theta estimate within each condition, for dimension 1, 2, and 3, respectively
  expect_equal(round(sqrt(fivenum(estimates_and_conditions[, "variance_estimate1"])), 3), c(.163, .262, .291, .317, .605))
  expect_equal(round(sqrt(fivenum(estimates_and_conditions[, "variance_estimate5"])), 3), c(.145, .208, .247, .277, .361))
  expect_equal(round(sqrt(fivenum(estimates_and_conditions[, "variance_estimate9"])), 3), c(.159, .207, .289, .316, .507))
  
  # some errors/missings for maximum_likelihood estimator because prior_var_safe_ml = NULL
  expect_equal(number_na_per_condition[, "x"], 
               c(rep(0, 18), 1, 0, 0, 0, 2, rep(0, 7), 1, 0))
})

test_that("three dimensions, maximum_aposteriori, 100 replications per condition, prior is normal()", {
  # To be run yet
  replications_per_unique_condition <- 100 
  true_theta_vec <- c(-2, 1, 3)
  number_items_vec <- c(300)
  number_answer_categories_vec <- c(2, 4)
  number_dimensions <- 3
  
  start_items <- list(type = 'random_by_dimension', n_by_dimension = 3, n = 9)
  variance_target <- .1^2
  model_vec <- c("3PLM","GRM","GPCM","SM")
  estimator_vec <- "maximum_aposteriori"
  information_summary_vec <- c("determinant", "posterior_determinant", "trace", "posterior_trace", "posterior_expected_kullback_leibler")
  
  estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, start_items, variance_target, replications_per_unique_condition, number_dimensions)
 
  conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, replications_per_unique_condition, number_dimensions)
  
  condition_vector <- sort(rep(1:(ncol(estimates_and_variance)/replications_per_unique_condition), replications_per_unique_condition))
  
  estimates_and_conditions <- cbind(t(estimates_and_variance)[, c("estimated_theta1", "estimated_theta2", "estimated_theta3", "variance_estimate1", "variance_estimate5", "variance_estimate9")], 
                                    conditions[, c("replication", "number_items", "number_answer_categories", "model", "estimator", "information_summary")], 
                                    condition_vector)
  
  average_per_condition_dim1 <- aggregate(estimates_and_conditions[, "estimated_theta1"], list(condition_vector), "mean", na.rm = TRUE)
  average_per_condition_dim2 <- aggregate(estimates_and_conditions[, "estimated_theta2"], list(condition_vector), "mean", na.rm = TRUE)
  average_per_condition_dim3 <- aggregate(estimates_and_conditions[, "estimated_theta3"], list(condition_vector), "mean", na.rm = TRUE)
  
  sd_per_condition_dim1 <- aggregate(estimates_and_conditions[, "estimated_theta1"], list(condition_vector), "sd", na.rm = TRUE)
  sd_per_condition_dim2 <- aggregate(estimates_and_conditions[, "estimated_theta2"], list(condition_vector), "sd", na.rm = TRUE)
  sd_per_condition_dim3 <- aggregate(estimates_and_conditions[, "estimated_theta3"], list(condition_vector), "sd", na.rm = TRUE)
  
  number_na_per_condition <- aggregate(estimates_and_conditions[, "estimated_theta1"], list(condition_vector), FUN = function(x) { sum(is.na(x)) })
  
  # five number summary of average theta estimate per condition, dimension 1 (true theta is -2)
  expect_equal(round(fivenum(average_per_condition_dim1[,"x"]), 3), c(-2.062, -2.034, -2.016, -1.996, -1.935))
  # five number summary of average theta estimate per condition, dimension 2 (true theta is 1)
  expect_equal(round(fivenum(average_per_condition_dim2[,"x"]), 3), c(.966, .990, 1.002, 1.024, 1.084))
  # five number summary of average theta estimate per condition, dimension 3 (true theta is 2)
  expect_equal(round(fivenum(average_per_condition_dim3[,"x"]), 3), c(1.956, 1.992, 2.009, 2.034, 2.104))
  
  # five number summary of observed sd of the theta estimates within each condition, for dimension 1, 2, and 3, respectively
  expect_equal(round(fivenum(sd_per_condition_dim1[,"x"]), 3), c(.167, .211, .312, .334, .385))
  expect_equal(round(fivenum(sd_per_condition_dim2[,"x"]), 3), c(.148, .177, .231, .248, .277))
  expect_equal(round(fivenum(sd_per_condition_dim3[,"x"]), 3), c(.170, .213, .300, .320, .384))
  
  # five number summary of reported sd of the theta estimate within each condition, for dimension 1, 2, and 3, respectively
  expect_equal(round(sqrt(fivenum(estimates_and_conditions[, "variance_estimate1"])), 3), c(.160, .261, .289, .315, .496))
  expect_equal(round(sqrt(fivenum(estimates_and_conditions[, "variance_estimate5"])), 3), c(.143, .208, .247, .277, .392))
  expect_equal(round(sqrt(fivenum(estimates_and_conditions[, "variance_estimate9"])), 3), c(.161, .213, .288, .314, .760)) 
})

test_that("three dimensions, maximum_aposteriori, 100 replications per condition, prior is uniform((-4, -4, -4), (4, 4, 4))", {
  replications_per_unique_condition <- 100 
  true_theta_vec <- c(-2, 1, 3)
  number_items_vec <- c(300)
  number_answer_categories_vec <- c(2, 4)
  number_dimensions <- 3
  
  start_items <- list(type = 'random_by_dimension', n_by_dimension = 3, n = 9)
  variance_target <- .1^2
  model_vec <- c("3PLM","GRM","GPCM","SM")
  estimator_vec <- "maximum_aposteriori"
  information_summary_vec <- c("determinant", "posterior_determinant", "trace", "posterior_trace", "posterior_expected_kullback_leibler")
  
  prior_form <- "uniform"
  prior_parameters <- list(lower_bound = rep(-4, 3), upper_bound = rep(4, 3))
  
  estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, start_items, variance_target, replications_per_unique_condition, number_dimensions, prior_form = prior_form, prior_parameters = prior_parameters, safe_eap = FALSE)
  
  conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, replications_per_unique_condition, number_dimensions)
  
  condition_vector <- sort(rep(1:(ncol(estimates_and_variance)/replications_per_unique_condition), replications_per_unique_condition))
  
  estimates_and_conditions <- cbind(t(estimates_and_variance)[, c("estimated_theta1", "estimated_theta2", "estimated_theta3", "variance_estimate1", "variance_estimate5", "variance_estimate9")], 
                                    conditions[, c("replication", "number_items", "number_answer_categories", "model", "estimator", "information_summary")], 
                                    condition_vector)
  
  average_per_condition_dim1 <- aggregate(estimates_and_conditions[, "estimated_theta1"], list(condition_vector), "mean", na.rm = TRUE)
  average_per_condition_dim2 <- aggregate(estimates_and_conditions[, "estimated_theta2"], list(condition_vector), "mean", na.rm = TRUE)
  average_per_condition_dim3 <- aggregate(estimates_and_conditions[, "estimated_theta3"], list(condition_vector), "mean", na.rm = TRUE)
  
  sd_per_condition_dim1 <- aggregate(estimates_and_conditions[, "estimated_theta1"], list(condition_vector), "sd", na.rm = TRUE)
  sd_per_condition_dim2 <- aggregate(estimates_and_conditions[, "estimated_theta2"], list(condition_vector), "sd", na.rm = TRUE)
  sd_per_condition_dim3 <- aggregate(estimates_and_conditions[, "estimated_theta3"], list(condition_vector), "sd", na.rm = TRUE)
    
  # five number summary of average theta estimate per condition, dimension 1 (true theta is -2)
  expect_equal(round(fivenum(average_per_condition_dim1[,"x"]), 3), c(-2.111, -2.044, -2.022, -2.003, -1.951))
  # five number summary of average theta estimate per condition, dimension 2 (true theta is 1)
  expect_equal(round(fivenum(average_per_condition_dim2[,"x"]), 3), c(.963, .986, 1.003, 1.021, 1.046))
  # five number summary of average theta estimate per condition, dimension 3 (true theta is 2)
  expect_equal(round(fivenum(average_per_condition_dim3[,"x"]), 3), c(2.929, 3.008, 3.037, 3.079, 3.173))
  
  # five number summary of observed sd of the theta estimates within each condition, for dimension 1, 2, and 3, respectively
  expect_equal(round(fivenum(sd_per_condition_dim1[,"x"]), 3), c(.180, .226, .316, .354, .427))
  expect_equal(round(fivenum(sd_per_condition_dim2[,"x"]), 3), c(.141, .180, .251, .275, .294))
  expect_equal(round(fivenum(sd_per_condition_dim3[,"x"]), 3), c(.216, .248, .435, .458, .526))
  
  # five number summary of reported sd of the theta estimate within each condition, for dimension 1, 2, and 3, respectively
  expect_equal(round(sqrt(fivenum(estimates_and_conditions[, "variance_estimate1"])), 3), c(.161, .209, .298, .334, .563))
  expect_equal(round(sqrt(fivenum(estimates_and_conditions[, "variance_estimate5"])), 3), c(.143, .186, .241, .257, .331))
  expect_equal(round(sqrt(fivenum(estimates_and_conditions[, "variance_estimate9"])), 3), c(.184, .252, .401, .479, .785)) 
})


test_that("three dimensions, expected_aposteriori, 100 replications per condition, prior is normal(), estimation procedure is gauss hermite quadrature", {
  #run again
  replications_per_unique_condition <- 100 
  true_theta_vec <- c(-2, 1, 3)
  number_items_vec <- 300
  number_answer_categories_vec <- c(2, 4)
  number_dimensions <- 3
  
  start_items <- list(type = 'random_by_dimension', n_by_dimension = 3, n = 9)
  variance_target <- .1^2
  model_vec <- c("3PLM","GRM","GPCM","SM")
  estimator_vec <- "expected_aposteriori"
  information_summary_vec <- c("determinant", "posterior_determinant", "trace", "posterior_trace", "posterior_expected_kullback_leibler")
  lowerbound <- -20
  upperbound <- 20
  prior <- diag(number_dimensions) * 9 
  
  estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, start_items, variance_target, replications_per_unique_condition, number_dimensions, lowerbound = lowerbound, upperbound = upperbound, prior = prior, eap_estimation_procedure = "gauss_hermite_quad")
  conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, replications_per_unique_condition, number_dimensions)
  
  condition_vector <- sort(rep(1:(ncol(estimates_and_variance)/replications_per_unique_condition), replications_per_unique_condition))
  
  estimates_and_conditions <- cbind(t(estimates_and_variance)[, c("estimated_theta1", "estimated_theta2", "estimated_theta3", "variance_estimate1", "variance_estimate5", "variance_estimate9")], 
                                    conditions[, c("replication", "number_items", "number_answer_categories", "model", "estimator", "information_summary")], 
                                    condition_vector)
  
  average_per_condition_dim1 <- aggregate(estimates_and_conditions[, "estimated_theta1"], list(condition_vector), "mean", na.rm = TRUE)
  average_per_condition_dim2 <- aggregate(estimates_and_conditions[, "estimated_theta2"], list(condition_vector), "mean", na.rm = TRUE)
  average_per_condition_dim3 <- aggregate(estimates_and_conditions[, "estimated_theta3"], list(condition_vector), "mean", na.rm = TRUE)
  
  sd_per_condition_dim1 <- aggregate(estimates_and_conditions[, "estimated_theta1"], list(condition_vector), "sd", na.rm = TRUE)
  sd_per_condition_dim2 <- aggregate(estimates_and_conditions[, "estimated_theta2"], list(condition_vector), "sd", na.rm = TRUE)
  sd_per_condition_dim3 <- aggregate(estimates_and_conditions[, "estimated_theta3"], list(condition_vector), "sd", na.rm = TRUE)
  
  # five number summary of average theta estimate per condition, dimension 1 (true theta is -2)
  expect_equal(round(fivenum(average_per_condition_dim1[,"x"]), 3), c(-2.113, -2.043, -2.023, -2.002, -1.961))
  # five number summary of average theta estimate per condition, dimension 2 (true theta is 1)
  expect_equal(round(fivenum(average_per_condition_dim2[,"x"]), 3), c(.963, .995, 1.005, 1.026, 1.057))
  # five number summary of average theta estimate per condition, dimension 3 (true theta is 2)
  expect_equal(round(fivenum(average_per_condition_dim3[,"x"]), 3), c(2.957, 2.996, 3.027, 3.068, 3.197))
  
  # five number summary of observed sd of the theta estimates within each condition, for dimension 1, 2, and 3, respectively
  expect_equal(round(fivenum(sd_per_condition_dim1[,"x"]), 3), c(.179, .207, .325, .339, .444))
  expect_equal(round(fivenum(sd_per_condition_dim2[,"x"]), 3), c(.145, .184, .239, .272, .339))
  expect_equal(round(fivenum(sd_per_condition_dim3[,"x"]), 3), c(.228, .268, .432, .481, .661))
  
  # five number summary of reported sd of the theta estimate within each condition, for dimension 1, 2, and 3, respectively
  expect_equal(round(sqrt(fivenum(estimates_and_conditions[, "variance_estimate1"])), 3), c(.007, .209, .298, .334, .677))
  expect_equal(round(sqrt(fivenum(estimates_and_conditions[, "variance_estimate5"])), 3), c(.012, .186, .242, .258, .343))
  expect_equal(round(sqrt(fivenum(estimates_and_conditions[, "variance_estimate9"])), 3), c(.055, .251, .399, .467, 1.025))
})

test_that("three dimensions, expected_aposteriori, 100 replications per condition, prior is normal(), estimation procedure is riemannsum", {
  #run again
  replications_per_unique_condition <- 100 
  true_theta_vec <- c(-2, 1, 3)
  number_items_vec <- 300
  number_answer_categories_vec <- c(2, 4)
  number_dimensions <- 3
  
  start_items <- list(type = 'random_by_dimension', n_by_dimension = 3, n = 9)
  variance_target <- .1^2
  model_vec <- c("3PLM","GRM","GPCM","SM")
  estimator_vec <- "expected_aposteriori"
  information_summary_vec <- c("determinant", "posterior_determinant", "trace", "posterior_trace", "posterior_expected_kullback_leibler")
  lowerbound <- -20
  upperbound <- 20
  prior <- diag(number_dimensions) * 9 
  
  estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, start_items, variance_target, replications_per_unique_condition, number_dimensions, lowerbound = lowerbound, upperbound = upperbound, prior = prior, eap_estimation_procedure = "riemannsum")
  conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, replications_per_unique_condition, number_dimensions)
  
  condition_vector <- sort(rep(1:(ncol(estimates_and_variance)/replications_per_unique_condition), replications_per_unique_condition))
  
  estimates_and_conditions <- cbind(t(estimates_and_variance)[, c("estimated_theta1", "estimated_theta2", "estimated_theta3", "variance_estimate1", "variance_estimate5", "variance_estimate9")], 
                                    conditions[, c("replication", "number_items", "number_answer_categories", "model", "estimator", "information_summary")], 
                                    condition_vector)
  
  average_per_condition_dim1 <- aggregate(estimates_and_conditions[, "estimated_theta1"], list(condition_vector), "mean", na.rm = TRUE)
  average_per_condition_dim2 <- aggregate(estimates_and_conditions[, "estimated_theta2"], list(condition_vector), "mean", na.rm = TRUE)
  average_per_condition_dim3 <- aggregate(estimates_and_conditions[, "estimated_theta3"], list(condition_vector), "mean", na.rm = TRUE)
  
  sd_per_condition_dim1 <- aggregate(estimates_and_conditions[, "estimated_theta1"], list(condition_vector), "sd", na.rm = TRUE)
  sd_per_condition_dim2 <- aggregate(estimates_and_conditions[, "estimated_theta2"], list(condition_vector), "sd", na.rm = TRUE)
  sd_per_condition_dim3 <- aggregate(estimates_and_conditions[, "estimated_theta3"], list(condition_vector), "sd", na.rm = TRUE)
  
  # five number summary of average theta estimate per condition, dimension 1 (true theta is -2)
  expect_equal(round(fivenum(average_per_condition_dim1[,"x"]), 3), c(-2.129, -2.061, -2.036, -2.019, -1.954))
  # five number summary of average theta estimate per condition, dimension 2 (true theta is 1)
  expect_equal(round(fivenum(average_per_condition_dim2[,"x"]), 3), c(.976, .999, 1.027, 1.041, 1.112))
  # five number summary of average theta estimate per condition, dimension 3 (true theta is 2)
  expect_equal(round(fivenum(average_per_condition_dim3[,"x"]), 3), c(2.942, 2.998, 3.037, 3.068, 3.150))
  
  # five number summary of observed sd of the theta estimates within each condition, for dimension 1, 2, and 3, respectively
  expect_equal(round(fivenum(sd_per_condition_dim1[,"x"]), 3), c(.185, .214, .315, .341, .358))
  expect_equal(round(fivenum(sd_per_condition_dim2[,"x"]), 3), c(.167, .244, .260, .284, .410))
  expect_equal(round(fivenum(sd_per_condition_dim3[,"x"]), 3), c(.218, .262, .440, .480, .584))
  
  # five number summary of reported sd of the theta estimate within each condition, for dimension 1, 2, and 3, respectively
  expect_equal(round(sqrt(fivenum(estimates_and_conditions[, "variance_estimate1"])), 3), c(.000, .212, .304, .340, .531))
  expect_equal(round(sqrt(fivenum(estimates_and_conditions[, "variance_estimate5"])), 3), c(.001, .189, .246, .262, .356))
  expect_equal(round(sqrt(fivenum(estimates_and_conditions[, "variance_estimate9"])), 3), c(.001, .254, .405, .477, 1.060))
})

test_that("three dimensions, expected_aposteriori, 100 replications per condition, prior is uniform((-4, -4, -4), (4, 4, 4))", {
  replications_per_unique_condition <- 100 
  true_theta_vec <- c(-2, 1, 3)
  number_items_vec <- c(300)
  number_answer_categories_vec <- c(2, 4)
  number_dimensions <- 3
  
  start_items <- list(type = 'random_by_dimension', n_by_dimension = 3, n = 9)
  variance_target <- .1^2
  model_vec <- c("3PLM","GRM","GPCM","SM")
  estimator_vec <- "expected_aposteriori"
  information_summary_vec <- c("determinant", "posterior_determinant", "trace", "posterior_trace", "posterior_expected_kullback_leibler")
  
  prior_form <- "uniform"
  prior_parameters <- list(lower_bound = rep(-4, 3), upper_bound = rep(4, 3))
  
  estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, start_items, variance_target, replications_per_unique_condition, number_dimensions, prior_form = prior_form, prior_parameters = prior_parameters, safe_eap = FALSE)
  
  conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, replications_per_unique_condition, number_dimensions)
  
  condition_vector <- sort(rep(1:(ncol(estimates_and_variance)/replications_per_unique_condition), replications_per_unique_condition))
  
  estimates_and_conditions <- cbind(t(estimates_and_variance)[, c("estimated_theta1", "estimated_theta2", "estimated_theta3", "variance_estimate1", "variance_estimate5", "variance_estimate9")], 
                                    conditions[, c("replication", "number_items", "number_answer_categories", "model", "estimator", "information_summary")], 
                                    condition_vector)
  
  average_per_condition_dim1 <- aggregate(estimates_and_conditions[, "estimated_theta1"], list(condition_vector), "mean", na.rm = TRUE)
  average_per_condition_dim2 <- aggregate(estimates_and_conditions[, "estimated_theta2"], list(condition_vector), "mean", na.rm = TRUE)
  average_per_condition_dim3 <- aggregate(estimates_and_conditions[, "estimated_theta3"], list(condition_vector), "mean", na.rm = TRUE)
  
  sd_per_condition_dim1 <- aggregate(estimates_and_conditions[, "estimated_theta1"], list(condition_vector), "sd", na.rm = TRUE)
  sd_per_condition_dim2 <- aggregate(estimates_and_conditions[, "estimated_theta2"], list(condition_vector), "sd", na.rm = TRUE)
  sd_per_condition_dim3 <- aggregate(estimates_and_conditions[, "estimated_theta3"], list(condition_vector), "sd", na.rm = TRUE)
  
  # five number summary of average theta estimate per condition, dimension 1 (true theta is -2)
  expect_equal(round(fivenum(average_per_condition_dim1[,"x"]), 3), c(-2.124, -2.071, -2.048, -2.020, -1.992))
  # five number summary of average theta estimate per condition, dimension 2 (true theta is 1)
  expect_equal(round(fivenum(average_per_condition_dim2[,"x"]), 3), c(.968, .995, 1.018, 1.032, 1.100))
  # five number summary of average theta estimate per condition, dimension 3 (true theta is 2)
  expect_equal(round(fivenum(average_per_condition_dim3[,"x"]), 3), c(2.971, 3.005, 3.027, 3.045, 3.109))
  
  # five number summary of observed sd of the theta estimates within each condition, for dimension 1, 2, and 3, respectively
  expect_equal(round(fivenum(sd_per_condition_dim1[,"x"]), 3), c(.168, .210, .321, .347, .373))
  expect_equal(round(fivenum(sd_per_condition_dim2[,"x"]), 3), c(.147, .196, .244, .264, .312))
  expect_equal(round(fivenum(sd_per_condition_dim3[,"x"]), 3), c(.222, .247, .311, .348, .388))
  
  # five number summary of reported sd of the theta estimate within each condition, for dimension 1, 2, and 3, respectively
  expect_equal(round(sqrt(fivenum(estimates_and_conditions[, "variance_estimate1"])), 3), c(.162, .210, .301, .341, .474))
  expect_equal(round(sqrt(fivenum(estimates_and_conditions[, "variance_estimate5"])), 3), c(.142, .187, .243, .259, .338))
  expect_equal(round(sqrt(fivenum(estimates_and_conditions[, "variance_estimate9"])), 3), c(.138, .249, .356, .406, .487)) 
})

  
  test_that("items load on all dimensions prior_var_safe_ml is 100", {
    #run again
    replications_per_unique_condition <- 100 
    true_theta_vec <- c(2, -1, -2)
    number_items_vec <- c(300)
    number_answer_categories_vec <- c(2, 4)
    number_dimensions <- 3
    
    start_items <- list(type = 'random_by_dimension', n_by_dimension = 3, n = 9)
    variance_target <- .1^2
    model_vec <- c("3PLM","GRM","GPCM","SM")
    estimator_vec <- c("maximum_likelihood", "maximum_aposteriori") # AEP
    information_summary_vec <- c("determinant", "posterior_determinant", "trace", "posterior_trace") # "posterior_expected_kullback_leibler"
    
    estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, start_items, variance_target, replications_per_unique_condition, number_dimensions, items_load_one_dimension = FALSE, prior_var_safe_ml = 100)
    #save(estimates_and_variance, file = "/Users/rivkadevries/Desktop/simulationsCAT/estimates_and_variance_within.R")
    estimates_and_variance_without_errors <- sapply(estimates_and_variance, 
                                                    FUN = function(x) { 
                                                      if (is.list(x)) 
                                                        matrix(rep(NA, number_dimensions + number_dimensions^2), nrow = 1, dimnames = list(NULL, names(estimates_and_variance[[1]])))
                                                      else 
                                                        x } ) 
    
    conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, replications_per_unique_condition, number_dimensions)
    
    condition_vector <- sort(rep(1:(ncol(estimates_and_variance_without_errors)/replications_per_unique_condition), replications_per_unique_condition))
    
    estimates_and_conditions <- cbind(t(estimates_and_variance_without_errors)[, c("estimated_theta1", "estimated_theta2", "estimated_theta3", "variance_estimate1", "variance_estimate5", "variance_estimate9")], 
                                      conditions[, c("replication", "number_items", "number_answer_categories", "model", "estimator", "information_summary")], 
                                      condition_vector)
    
    average_per_condition_dim1 <- aggregate(estimates_and_conditions[, "estimated_theta1"], list(condition_vector), "mean", na.rm = TRUE)
    average_per_condition_dim2 <- aggregate(estimates_and_conditions[, "estimated_theta2"], list(condition_vector), "mean", na.rm = TRUE)
    average_per_condition_dim3 <- aggregate(estimates_and_conditions[, "estimated_theta3"], list(condition_vector), "mean", na.rm = TRUE)
    
    sd_per_condition_dim1 <- aggregate(estimates_and_conditions[, "estimated_theta1"], list(condition_vector), "sd", na.rm = TRUE)
    sd_per_condition_dim2 <- aggregate(estimates_and_conditions[, "estimated_theta2"], list(condition_vector), "sd", na.rm = TRUE)
    sd_per_condition_dim3 <- aggregate(estimates_and_conditions[, "estimated_theta3"], list(condition_vector), "sd", na.rm = TRUE)
    
    number_na_per_condition <- aggregate(estimates_and_conditions[, "estimated_theta1"], list(condition_vector), FUN = function(x) { sum(is.na(x)) })
    
    # five number summary of average theta estimate per condition, dimension 1 (true theta is -2)
    expect_equal(round(fivenum(average_per_condition_dim1[,"x"]), 3), c(1.933, 1.989, 2.006, 2.042, 2.154))
    # five number summary of average theta estimate per condition, dimension 2 (true theta is 1)
    expect_equal(round(fivenum(average_per_condition_dim2[,"x"]), 3), c(-1.083, -1.027, -.995, -.976, -.714))
    # five number summary of average theta estimate per condition, dimension 3 (true theta is 2)
    expect_equal(round(fivenum(average_per_condition_dim3[,"x"]), 3), c(-2.103, -2.034, -2.012, -1.985, -1.907))
    
    # five number summary of observed sd of the theta estimates within each condition, for dimension 1, 2, and 3, respectively
    expect_equal(round(fivenum(sd_per_condition_dim1[,"x"]), 3), c(.189, .276, .375, .406, .660))
    expect_equal(round(fivenum(sd_per_condition_dim2[,"x"]), 3), c(.203, .261, .346, .420, 1.190))
    expect_equal(round(fivenum(sd_per_condition_dim3[,"x"]), 3), c(.205, .252, .375, .419, .836))
    
    # five number summary of reported sd of the theta estimate within each condition, for dimension 1, 2, and 3, respectively
    expect_equal(round(sqrt(fivenum(estimates_and_conditions[, "variance_estimate1"])), 3), c(.100, .377, .423, .568, 2.124))
    expect_equal(round(sqrt(fivenum(estimates_and_conditions[, "variance_estimate5"])), 3), c(.100, .365, .419, .549, 2.178))
    expect_equal(round(sqrt(fivenum(estimates_and_conditions[, "variance_estimate9"])), 3), c(.100, .382, .427, .544, 1.992))
    
    # no errors/missings for maximum_aposteriori estimator
    expect_equal(number_na_per_condition[which(conditions[seq(1, 6400, 100), "estimator"] == "maximum_aposteriori"), "x"], rep(0, 32))
    # some errors/missings for maximum_likelihood estimator
    expect_equal(number_na_per_condition[which(conditions[seq(1, 6400, 100), "estimator"] == "maximum_likelihood"), "x"], 
                 c(0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 1, 0))
    
  })
  
  test_that("items load on all dimensions prior_var_safe_ml is 1", {
    #run again
    replications_per_unique_condition <- 100 
    true_theta_vec <- c(2, -1, -2)
    number_items_vec <- c(300)
    number_answer_categories_vec <- c(2, 4)
    number_dimensions <- 3
    
    start_items <- list(type = 'random_by_dimension', n_by_dimension = 3, n = 9)
    variance_target <- .1^2
    model_vec <- c("3PLM","GRM","GPCM","SM")
    estimator_vec <- c("maximum_likelihood", "maximum_aposteriori") # AEP
    information_summary_vec <- c("determinant", "posterior_determinant", "trace", "posterior_trace") # PEKL 
    
    estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, start_items, variance_target, replications_per_unique_condition, number_dimensions, items_load_one_dimension = FALSE, prior_var_safe_ml = 1)
    #save(estimates_and_variance, file = "/Users/rivkadevries/Desktop/simulationsCAT/estimates_and_variance_within_safe_var1.R")
    estimates_and_variance_without_errors <- sapply(estimates_and_variance, 
                                                    FUN = function(x) { 
                                                      if (is.list(x)) 
                                                        matrix(rep(NA, number_dimensions + number_dimensions^2), nrow = 1, dimnames = list(NULL, names(estimates_and_variance[[1]])))
                                                      else 
                                                        x } ) 
    
    conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, replications_per_unique_condition, number_dimensions)
    
    condition_vector <- sort(rep(1:(ncol(estimates_and_variance_without_errors)/replications_per_unique_condition), replications_per_unique_condition))
    
    estimates_and_conditions <- cbind(t(estimates_and_variance_without_errors)[, c("estimated_theta1", "estimated_theta2", "estimated_theta3", "variance_estimate1", "variance_estimate5", "variance_estimate9")], 
                                      conditions[, c("replication", "number_items", "number_answer_categories", "model", "estimator", "information_summary")], 
                                      condition_vector)
    
    average_per_condition_dim1 <- aggregate(estimates_and_conditions[, "estimated_theta1"], list(condition_vector), "mean", na.rm = TRUE)
    average_per_condition_dim2 <- aggregate(estimates_and_conditions[, "estimated_theta2"], list(condition_vector), "mean", na.rm = TRUE)
    average_per_condition_dim3 <- aggregate(estimates_and_conditions[, "estimated_theta3"], list(condition_vector), "mean", na.rm = TRUE)
    
    sd_per_condition_dim1 <- aggregate(estimates_and_conditions[, "estimated_theta1"], list(condition_vector), "sd", na.rm = TRUE)
    sd_per_condition_dim2 <- aggregate(estimates_and_conditions[, "estimated_theta2"], list(condition_vector), "sd", na.rm = TRUE)
    sd_per_condition_dim3 <- aggregate(estimates_and_conditions[, "estimated_theta3"], list(condition_vector), "sd", na.rm = TRUE)
    
    number_na_per_condition <- aggregate(estimates_and_conditions[, "estimated_theta1"], list(condition_vector), FUN = function(x) { sum(is.na(x)) })
    
    # five number summary of average theta estimate per condition, dimension 1 (true theta is -2)
    expect_equal(round(fivenum(average_per_condition_dim1[,"x"]), 3), c(1.944, 1.983, 2.008, 2.023, 2.155))
    # five number summary of average theta estimate per condition, dimension 2 (true theta is 1)
    expect_equal(round(fivenum(average_per_condition_dim2[,"x"]), 3), c(-1.088, -1.028, -1.009, -.984, -.930))
    # five number summary of average theta estimate per condition, dimension 3 (true theta is 2)
    expect_equal(round(fivenum(average_per_condition_dim3[,"x"]), 3), c(-2.140, -2.026, -2.005, -1.981, -1.907))
    
    # five number summary of observed sd of the theta estimates within each condition, for dimension 1, 2, and 3, respectively
    expect_equal(round(fivenum(sd_per_condition_dim1[,"x"]), 3), c(.206, .266, .358, .391, .440))
    expect_equal(round(fivenum(sd_per_condition_dim2[,"x"]), 3), c(.189, .260, .347, .378, .447))
    expect_equal(round(fivenum(sd_per_condition_dim3[,"x"]), 3), c(.196, .264, .355, .390, .444))
    
    # five number summary of reported sd of the theta estimate within each condition, for dimension 1, 2, and 3, respectively
    expect_equal(round(sqrt(fivenum(estimates_and_conditions[, "variance_estimate1"])), 3), c(.197, .379, .425, .576, 2.190))
    expect_equal(round(sqrt(fivenum(estimates_and_conditions[, "variance_estimate5"])), 3), c(.195, .367, .420, .557, 2.749))
    expect_equal(round(sqrt(fivenum(estimates_and_conditions[, "variance_estimate9"])), 3), c(.195, .383, .429, .553, 2.143))
    
    # no errors/missings for maximum_aposteriori estimator
    expect_equal(number_na_per_condition[which(conditions[seq(1, 6400, 100), "estimator"] == "maximum_aposteriori"), "x"], rep(0, 32))
    # some errors/missings for maximum_likelihood estimator
    expect_equal(number_na_per_condition[which(conditions[seq(1, 6400, 100), "estimator"] == "maximum_likelihood"), "x"], 
                 c(0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0))
    
  })
  
  
  test_that("simulate with constraints, max_n 130", {
    #run again
    replications_per_unique_condition <- 100
    true_theta_vec <- c(-2, 1, 2)
    number_items_vec <- 300 # can only have length one here (with item characteristics) and should be divisible by 3, to keep things simple
    number_answer_categories_vec <- 4
    number_dimensions <- 3
    
    start_items <- list(type = 'random_by_dimension', n_by_dimension = 3, n = 9)
    variance_target <- .001^2
    model_vec <- "SM"
    estimator_vec <- c("maximum_likelihood", "maximum_aposteriori") # AEP
    information_summary_vec <- "determinant" 
    max_n = 130
    
    #create item characteristics and constraints
    characteristics <- data.frame(content = c(rep('depression', number_items_vec / 3), rep('anxiety', number_items_vec / 3), rep('somatic', number_items_vec / 3)))
    constraints <- list(list(name = 'content/depression',
                             op = '><',
                             target = c(50, 75)),
                        list(name = 'content/somatic',
                             op = '><',
                             target = c(75, 100)))
    
    estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, start_items, variance_target, replications_per_unique_condition, number_dimensions, constraints_and_characts = list(characteristics = characteristics, constraints = constraints), prior_var_safe_ml = 100, return_administered_item_indices = TRUE, max_n = max_n)
    #save(estimates_and_variance, file = "/Users/rivkadevries/Desktop/simulationsCAT/estimates_and_variance_3dim_constraints_maxn130.R")
    
    conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, replications_per_unique_condition, number_dimensions)
    
    condition_vector <- sort(rep(1:(ncol(estimates_and_variance)/replications_per_unique_condition), replications_per_unique_condition))
    
    estimates_and_conditions <- cbind(t(estimates_and_variance)[, c("estimated_theta1", "estimated_theta2", "estimated_theta3", "variance_estimate1", "variance_estimate5", "variance_estimate9",
                                                                    str_c("items_administered", 1:max_n))], 
                                      conditions[, c("replication", "number_items", "number_answer_categories", "model", "estimator", "information_summary")], 
                                      condition_vector)
    
    average_per_condition_dim1 <- aggregate(estimates_and_conditions[, "estimated_theta1"], list(condition_vector), "mean", na.rm = TRUE)
    average_per_condition_dim2 <- aggregate(estimates_and_conditions[, "estimated_theta2"], list(condition_vector), "mean", na.rm = TRUE)
    average_per_condition_dim3 <- aggregate(estimates_and_conditions[, "estimated_theta3"], list(condition_vector), "mean", na.rm = TRUE)
    
    sd_per_condition_dim1 <- aggregate(estimates_and_conditions[, "estimated_theta1"], list(condition_vector), "sd", na.rm = TRUE)
    sd_per_condition_dim2 <- aggregate(estimates_and_conditions[, "estimated_theta2"], list(condition_vector), "sd", na.rm = TRUE)
    sd_per_condition_dim3 <- aggregate(estimates_and_conditions[, "estimated_theta3"], list(condition_vector), "sd", na.rm = TRUE)
    
    number_na_per_condition <- aggregate(estimates_and_conditions[, "estimated_theta1"], list(condition_vector), FUN = function(x) { sum(is.na(x)) })
    
    number_depression_items <- apply(estimates_and_conditions[,str_c("items_administered", 1:max_n)], 
                                     1, 
                                     FUN = function(x) {sum(x <= 100)} )
    number_anxiety_items <- apply(estimates_and_conditions[,str_c("items_administered", 1:max_n)], 
                                  1, 
                                  FUN = function(x) {sum((x > 100 & x <= 200))} )
    number_somatic_items <- apply(estimates_and_conditions[,str_c("items_administered", 1:max_n)], 
                                  1, 
                                  FUN = function(x) {sum(x > 200)} )
    
    # average theta estimate per condition, dimension 1 (true theta is -2)
    expect_equal(round(average_per_condition_dim1[,"x"], 3), c(-2.006, -1.992))
    # average theta estimate per condition, dimension 2 (true theta is 1)
    expect_equal(round(average_per_condition_dim2[,"x"], 3), c(.996, 1.056))
    # average theta estimate per condition, dimension 3 (true theta is 2)
    expect_equal(round(average_per_condition_dim3[,"x"], 3), c(2.000, 2.027))
    
    # five number summary of observed sd of the theta estimates within each condition, for dimension 1, 2, and 3, respectively
    expect_equal(round(sd_per_condition_dim1[,"x"], 3), c(.286, .272))
    expect_equal(round(sd_per_condition_dim2[,"x"], 3), c(.756, .626))
    expect_equal(round(sd_per_condition_dim3[,"x"], 3), c(.188, .216))
    
    # five number summary of reported sd of the theta estimate within each condition, for dimension 1, 2, and 3, respectively
    expect_equal(round(sqrt(fivenum(estimates_and_conditions[, "variance_estimate1"])), 3), c(.240, .277, .295, .311, .397))
    expect_equal(round(sqrt(fivenum(estimates_and_conditions[, "variance_estimate5"])), 3), c(.546, .618, .674, .711, .841))
    expect_equal(round(sqrt(fivenum(estimates_and_conditions[, "variance_estimate9"])), 3), c(.172, .182, .188, .193, .214))
    
    # no errors/missings for maximum_aposteriori estimator
    expect_equal(number_na_per_condition[which(conditions[seq(1, 6400, 100), "estimator"] == "maximum_aposteriori"), "x"], 0)
    # no errors/missings for maximum_likelihood estimator
    expect_equal(number_na_per_condition[which(conditions[seq(1, 6400, 100), "estimator"] == "maximum_likelihood"), "x"], 0)
    
    expect_equal(number_depression_items, rep(50, 200))
    expect_equal(number_anxiety_items, rep(5, 200))
    expect_equal(number_somatic_items, rep(75, 200)) 
  })
  
  test_that("simulate with constraints, max_n 260", {
    #run again
    replications_per_unique_condition <- 100
    true_theta_vec <- c(-2, 1, 2)
    number_items_vec <- 300 # can only have length one here (with item characteristics) and should be divisible by 3, to keep things simple
    number_answer_categories_vec <- 4
    number_dimensions <- 3
    
    start_items <- list(type = 'random_by_dimension', n_by_dimension = 3, n = 9)
    variance_target <- .001^2
    model_vec <- "SM"
    estimator_vec <- c("maximum_likelihood", "maximum_aposteriori") # AEP
    information_summary_vec <- "determinant"
    max_n = 260
    
    #create item characteristics and constraints
    characteristics <- data.frame(content = c(rep('depression', number_items_vec / 3), rep('anxiety', number_items_vec / 3), rep('somatic', number_items_vec / 3)))
    constraints <- list(list(name = 'content/depression',
                             op = '><',
                             target = c(50, 75)),
                        list(name = 'content/somatic',
                             op = '><',
                             target = c(75, 90)))
    
    estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, start_items, variance_target, replications_per_unique_condition, number_dimensions, constraints_and_characts = list(characteristics = characteristics, constraints = constraints), prior_var_safe_ml = 100, return_administered_item_indices = TRUE, max_n = max_n)
    #save(estimates_and_variance, file = "/Users/rivkadevries/Desktop/simulationsCAT/estimates_and_variance_3dim_constraints_maxn260.R")
    
    conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, replications_per_unique_condition, number_dimensions)
    
    condition_vector <- sort(rep(1:(ncol(estimates_and_variance)/replications_per_unique_condition), replications_per_unique_condition))
    
    estimates_and_conditions <- cbind(t(estimates_and_variance)[, c("estimated_theta1", "estimated_theta2", "estimated_theta3", "variance_estimate1", "variance_estimate5", "variance_estimate9",
                                                                    str_c("items_administered", 1:max_n))], 
                                      conditions[, c("replication", "number_items", "number_answer_categories", "model", "estimator", "information_summary")], 
                                      condition_vector)
    
    number_depression_items <- apply(estimates_and_conditions[,str_c("items_administered", 1:max_n)], 
                                     1, 
                                     FUN = function(x) {sum(x <= 100)} )
    number_anxiety_items <- apply(estimates_and_conditions[,str_c("items_administered", 1:max_n)], 
                                  1, 
                                  FUN = function(x) {sum((x > 100 & x <= 200))} )
    number_somatic_items <- apply(estimates_and_conditions[,str_c("items_administered", 1:max_n)], 
                                  1, 
                                  FUN = function(x) {sum(x > 200)} )
    
    # with max_n small, the minimum of 50 depression and 75 somatic items is not reached
    expect_equal(number_depression_items, rep(75, 200))
    expect_equal(all(number_anxiety_items < 98 & number_anxiety_items > 94), TRUE)
    expect_equal(all(number_somatic_items < 91 & number_somatic_items > 87), TRUE)
  })
  
  test_that("maximum_aposteriori with informative prior", {
    # run again with PEKL
    replications_per_unique_condition <- 100
    true_theta_vec <- c(-2, 1, 2)
    number_items_vec <- c(300)
    number_answer_categories_vec <- c(2, 4)
    number_dimensions <- 3
    
    start_items <- list(type = 'random_by_dimension', n_by_dimension = 3, n = 9)
    variance_target <- .1^2
    model_vec <- c("3PLM","GRM","GPCM","SM")
    estimator_vec <- "maximum_aposteriori"
    information_summary_vec <- c("determinant", "posterior_determinant", "trace", "posterior_trace") # PEKL 
    
    estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, start_items, variance_target, replications_per_unique_condition, number_dimensions, prior = diag(3) * .5, prior_var_safe_ml = 100)
    
    conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, replications_per_unique_condition, number_dimensions)
    
    condition_vector <- sort(rep(1:(ncol(estimates_and_variance)/replications_per_unique_condition), replications_per_unique_condition))
    
    estimates_and_conditions <- cbind(t(estimates_and_variance)[, c("estimated_theta1", "estimated_theta2", "estimated_theta3", "variance_estimate1", "variance_estimate5", "variance_estimate9")], 
                                      conditions[, c("replication", "number_items", "number_answer_categories", "model", "estimator", "information_summary")], 
                                      condition_vector)
    
    average_per_condition_dim1 <- aggregate(estimates_and_conditions[, "estimated_theta1"], list(condition_vector), "mean", na.rm = TRUE)
    average_per_condition_dim2 <- aggregate(estimates_and_conditions[, "estimated_theta2"], list(condition_vector), "mean", na.rm = TRUE)
    average_per_condition_dim3 <- aggregate(estimates_and_conditions[, "estimated_theta3"], list(condition_vector), "mean", na.rm = TRUE)
    
    sd_per_condition_dim1 <- aggregate(estimates_and_conditions[, "estimated_theta1"], list(condition_vector), "sd", na.rm = TRUE)
    sd_per_condition_dim2 <- aggregate(estimates_and_conditions[, "estimated_theta2"], list(condition_vector), "sd", na.rm = TRUE)
    sd_per_condition_dim3 <- aggregate(estimates_and_conditions[, "estimated_theta3"], list(condition_vector), "sd", na.rm = TRUE)
    
    number_na_per_condition <- aggregate(estimates_and_conditions[, "estimated_theta1"], list(condition_vector), FUN = function(x) { sum(is.na(x)) })
    
    # five number summary of average theta estimate per condition, dimension 1 (true theta is -2)
    expect_equal(round(fivenum(average_per_condition_dim1[,"x"]), 3), c(-1.903, -1.856, -1.714, -1.690, -1.659))
    # five number summary of average theta estimate per condition, dimension 2 (true theta is 1)
    expect_equal(round(fivenum(average_per_condition_dim2[,"x"]), 3), c(.862, .890, .913, .936, .989))
    # five number summary of average theta estimate per condition, dimension 3 (true theta is 2)
    expect_equal(round(fivenum(average_per_condition_dim3[,"x"]), 3), c(1.644, 1.693, 1.713, 1.840, 1.899))
    
    # five number summary of observed sd of the theta estimates within each condition, for dimension 1, 2, and 3, respectively
    expect_equal(round(fivenum(sd_per_condition_dim1[,"x"]), 3), c(.150, .185, .219, .230, .255))
    expect_equal(round(fivenum(sd_per_condition_dim2[,"x"]), 3), c(.150, .178, .201, .215, .247))
    expect_equal(round(fivenum(sd_per_condition_dim3[,"x"]), 3), c(.150, .192, .218, .235, .261))
    
    # five number summary of reported sd of the theta estimate within each condition, for dimension 1, 2, and 3, respectively
    expect_equal(round(sqrt(fivenum(estimates_and_conditions[, "variance_estimate1"])), 3), c(.157, .237, .255, .270, .336))
    expect_equal(round(sqrt(fivenum(estimates_and_conditions[, "variance_estimate5"])), 3), c(.141, .197, .231, .249, .294))
    expect_equal(round(sqrt(fivenum(estimates_and_conditions[, "variance_estimate9"])), 3), c(.158, .199, .254, .269, .335))
  })
  
}

