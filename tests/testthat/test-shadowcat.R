# only for whithin R:
'
library(testthat)
library(pbapply)
library(stringr)
'

make_random_seed_exist <- rnorm(1)

#' simulate a testing routine with shadowcat
#' 
#' @param true_theta true theta value or vector
#' @param prior covariance matrix of the (multi variate) normal prior for theta; mean vector is fixed at zero; not used for maximum_likelihood estimator
#' @param model String, one of '3PLM', 'GPCM', 'SM' or 'GRM', for the three-parameter logistic, generalized partial credit, sequential or graded response model respectively.
#' @param alpha Matrix of alpha parameters, one column per dimension, one row per item. Row names should contain the item keys. Note that so called within-dimensional models still use an alpha matrix, they simply 
#' have only one non-zero loading per item.
#' @param beta Matrix of beta parameters, one column per item step, one row per item. Row names should contain the item keys. Note that ShadowCAT expects response categories to be sequential,
#' and without gaps. That is, the weight parameter in the GPCM model is assumed to be sequential, and equal to the position of the 'location' of the beta parameter in the Beta matrix.
#' The matrix will have a number of columns equal to the largest number of response categories, items with fewer response categories should be 
#' right-padded with \code{NA}. \code{NA} values between response categories are not allowed, and will lead to errors.
#' Beta matrix can be set to NULL if model is GPCM and eta is defined
#' @param guessing matrix with one column of guessing parameters per item. Row names should contain the item keys. Optionally used in 3PLM model, ignored for all others.
#' @param eta Matrix of location parameters, optionally used in GPCM model, ignored for all others.
#' @param start_items items that are shown to the patient before adaptive proces starts; one of
#' list(type = 'random', n)
#' list(type = 'fixed', indeces, n)
#' list(type = 'random_by_dimension', n_by_dimension, n)
#' where n = total number of initial items, indeces = vector of initial item indeces, 
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
#' @param lower_bound vector with lower bounds for theta per dimension; estimated theta values smaller than the lowerbound values are truncated to the lowerbound values 
#' @param upper_bound vector with upper bounds for theta per dimension; estimated theta values larger than the upperbound values are truncated to the upperbound values
#' @param prior_var_safe_ml if not NULL, expected_aposteriori estimate with prior variance equal to prior_var_safe_ml (scalar or vector) is computed instead of maximum_likelihood/maximum_aposteriori, if maximum_likelihood/maximum_aposteriori estimate fails.
#' @param initital_estimate vector containing the initial theta estimates (starting values)
#' @param initial_variance matrix containing the initial covariance matrix (staring values)
#' @return
test_shadowcat <- function(true_theta, prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, constraints_and_characts = NULL, lowerbound = rep(-3, ncol(alpha)), upperbound = rep(3, ncol(alpha)), prior_var_safe_ml = NULL, initital_estimate = rep(0, ncol(alpha)), initial_variance = prior) {
  item_keys <- rownames(alpha)
  responses <- NULL
  attr(initital_estimate, 'variance') <- initial_variance
  next_item_and_test_outcome <- shadowcat(responses, estimate = initital_estimate, model, alpha, beta, start_items, stop_test, estimator, information_summary, prior, guessing, eta, constraints_and_characts, lowerbound, upperbound, prior_var_safe_ml)

  while (next_item_and_test_outcome$key_new_item != "stop_test") {
    new_response <- simulate_answer(true_theta, model, ncol(alpha), estimator, alpha, beta, guessing, ncol(beta), match(next_item_and_test_outcome$key_new_item, item_keys))
    next_item_and_test_outcome$responses[[next_item_and_test_outcome$key_new_item]] <- new_response
    next_item_and_test_outcome$responses <- as.list(next_item_and_test_outcome$responses)
    next_item_and_test_outcome <- shadowcat(next_item_and_test_outcome$responses, next_item_and_test_outcome$estimate, model, alpha, beta, start_items, stop_test, estimator, information_summary, prior,  guessing, eta, constraints_and_characts, lowerbound, upperbound, prior_var_safe_ml)  
  }
  
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
#' list(type = 'fixed', indeces, n)
#' list(type = 'random_by_dimension', n_by_dimension, n)
#' where n = total number of initial items, indeces = vector of initial item indeces, 
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
#' @param lower_bound vector with lower bounds for theta per dimension; estimated theta values smaller than the lowerbound values are truncated to the lowerbound values 
#' @param upper_bound vector with upper bounds for theta per dimension; estimated theta values larger than the upperbound values are truncated to the upperbound values
#' @param prior covariance matrix of the (multi variate) normal prior for theta; mean vector is fixed at zero; not used for maximum_likelihood estimator
#' @param prior_var_safe_ml if not NULL, expected_aposteriori estimate with prior variance equal to prior_var_safe_ml (scalar or vector) is computed instead of maximum_likelihood/maximum_aposteriori, if maximum_likelihood/maximum_aposteriori estimate fails.
#' @param return_administered_item_indeces if TRUE, indeces of administered items are added to the output
#' @param min_n value equal to the minimum number of items to administer
#' @param max_n value equal to the maximum number of items to administer (test stops at this number, even if variance target has not been reached). NULL means max_n is equal to number of items in test bank
#' @param varying_number_item_steps if TRUE, the simulated number of item steps differs over items. In that case, number_answer_categories_vec (number of itemsteps + 1)
#' is considered the maxixmum number of categories
#' @return matrix with in each row (= one condition): named vector containing estimated theta, variance of the estimate, and if return_administered_item_indeces is TRUE, the indeces of the administered items
run_simulation <- function(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, start_items, variance_target, replications_per_unique_condition, number_dimensions, constraints_and_characts = NULL, guessing = NULL, items_load_one_dimension = TRUE, lowerbound = rep(-3, number_dimensions), upperbound = rep(3, number_dimensions), prior = diag(number_dimensions) * 20, prior_var_safe_ml = NULL, return_administered_item_indeces = FALSE, min_n = NULL, max_n = NULL, varying_number_item_steps = FALSE) {                   
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
                      estimate_theta <- tryCatch(test_shadowcat(true_theta, prior, as.character(conditions[condition, "model"]), alpha_beta$alpha, alpha_beta$beta, guessing, eta = NULL, start_items, stop_test, as.character(conditions[condition, "estimator"]), as.character(conditions[condition, "information_summary"]), constraints_and_characts, lowerbound, upperbound, prior_var_safe_ml),
                                                 error = function(e) e)

                      if (return_administered_item_indeces)
                        c("estimated_theta" = estimate_theta$estimate,
                          "variance_estimate" = attr(estimate_theta$estimate, "variance"),
                          "items_administered" = as.numeric(sapply(names(estimate_theta$responses), substring, 5)))
                      else
                        c("estimated_theta" = estimate_theta$estimate,
                          "variance_estimate" =  attr(estimate_theta$estimate, "variance"))               
                    })
}

context("validate shadowcat single conditions")

test_that("true theta is 2, estimator is maximum_aposteriori", {
  # define true theta for simulation of responses
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
  
  # define prior covariance matrix
  prior <- diag(number_dimensions) * 5
  
  test_outcome <- with_random_seed(2, test_shadowcat)(true_theta, prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary)
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), 1.938)
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3)), .113)
  expect_equal(length(test_outcome$responses), 100)
})

test_that("true theta is 2, estimator is maximum_likelihood", {
  # define true theta for simulation of responses
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
  
  test_outcome <- with_random_seed(2, test_shadowcat)(true_theta, prior = NULL, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, initial_variance = diag(1))
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), 2.169)
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3)), .129)
  expect_equal(length(test_outcome$responses), 100)
})

test_that("true theta is 2, estimator is expected_aposteriori", {
  # define true theta for simulation of responses
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
  
  prior <- diag(number_dimensions) * 5
  
  test_outcome <- with_random_seed(2, test_shadowcat)(true_theta, prior = prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, initial_variance = diag(1))
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), 1.833)
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3)), .087)
  expect_equal(length(test_outcome$responses), 100)
})


test_that("true theta is 1, 0, 2, estimator is maximum_aposteriori", {  
  # define true theta for simulation of responses
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
  prior <- diag(number_dimensions) * 20
  
  test_outcome <- with_random_seed(3, test_shadowcat)(true_theta, prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary)
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), c(.841, -.123, 1.947))
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3))[1:3],c(.064, .000, .000))
  expect_equal(length(test_outcome$responses), 300)
})

test_that("true theta is 1, 0, 2, estimator is maximum_likelihood", {  
  # define true theta for simulation of responses
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
  
  test_outcome <- with_random_seed(3, test_shadowcat)(true_theta, prior = NULL, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, initial_variance = diag(1))
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), c(.755, -.070, 2.221))
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3))[1:3],c(.063, .000, .000))
  expect_equal(length(test_outcome$responses), 300)
  
  # defining prior has no effect on outcome
  test_outcome <- with_random_seed(3, test_shadowcat)(true_theta, prior = diag(number_dimensions) * 2, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, initial_variance = diag(1))
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), c(.755, -.070, 2.221))
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3))[1:3],c(.063, .000, .000))
  expect_equal(length(test_outcome$responses), 300)  
})

test_that("true theta is 1, 0, 2, estimator is expected_aposteriori", {  
  # define true theta for simulation of responses
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
  
  # define prior covariance matrix
  prior <- diag(number_dimensions) * 20
  
  test_outcome <- with_random_seed(3, test_shadowcat)(true_theta, prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary)
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), c(1.423, -.087, 1.849))
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3))[1:3],c(.075, .000, .000))
  expect_equal(length(test_outcome$responses), 300)
})


test_that("items load on three dimensions", {  
  # define true theta for simulation of responses
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
  prior <- diag(number_dimensions) * 20
  
  test_outcome <- with_random_seed(3, test_shadowcat)(true_theta, prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary)
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), c(.454, .719, 1.745))
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3))[1:3],c(.178, -.088, -.081))
  expect_equal(length(test_outcome$responses), 300)
})

test_that("true theta is 2, 2, 2", {  
  # define true theta for simulation of responses
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
  
  # define prior covariance matrix
  prior <- diag(number_dimensions) * 20
  
  test_outcome <- with_random_seed(3, test_shadowcat)(true_theta, prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary)
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), c(2.287, 1.825, 1.672))
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3))[1:3],c(.11, .000, .000))
  expect_equal(length(test_outcome$responses), 300)
})

test_that("with constraints max_n 260", {  
  # define true theta for simulation of responses
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
  
  # define prior covariance matrix
  prior <- diag(number_dimensions) * 20
  
  #create item characteristics and constraints
  characteristics <- data.frame(content = c(rep('depression', number_items / 3), rep('anxiety', number_items / 3), rep('somatic', number_items / 3)))
  constraints <- list(list(name = 'content/depression',
                           op = '><',
                           target = c(50, 75)),
                      list(name = 'content/somatic',
                           op = '><',
                           target = c(75, 90)))
  
  test_outcome <- with_random_seed(3, test_shadowcat)(true_theta, prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, constraints_and_characts = list(characteristics = characteristics, constraints = constraints))
  indeces_administered <- as.numeric(sapply(names(test_outcome$responses), substring, 5))
  
  number_depression_items <- sum(indeces_administered <= 100)
  number_anxiety_items <- sum(indeces_administered > 100 & indeces_administered <= 200)
  number_somatic_items <- sum(indeces_administered > 200)
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), c(-2.689, .697, 2.260))
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3))[1:3],c(.155, .000, .000))
  expect_equal(length(test_outcome$responses), 260)
  expect_equal(number_depression_items, 75)
  expect_equal(number_anxiety_items, 95)
  expect_equal(number_somatic_items, 90) 
})

test_that("with constraints max_n 130", {  
  # define true theta for simulation of responses
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
  
  # define prior covariance matrix
  prior <- diag(number_dimensions) * 20
  
  #create item characteristics and constraints
  characteristics <- data.frame(content = c(rep('depression', number_items / 3), rep('anxiety', number_items / 3), rep('somatic', number_items / 3)))
  constraints <- list(list(name = 'content/depression',
                           op = '><',
                           target = c(50, 75)),
                      list(name = 'content/somatic',
                           op = '><',
                           target = c(75, 90)))
  
  test_outcome <- with_random_seed(3, test_shadowcat)(true_theta, prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, constraints_and_characts = list(characteristics = characteristics, constraints = constraints))
  indeces_administered <- as.numeric(sapply(names(test_outcome$responses), substring, 5))
  
  number_depression_items <- sum(indeces_administered <= 100)
  number_anxiety_items <- sum(indeces_administered > 100 & indeces_administered <= 200)
  number_somatic_items <- sum(indeces_administered > 200)
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), c(-1.577, .071, 1.906))
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3))[1:3],c(.106, .000, .000))
  expect_equal(length(test_outcome$responses), 130)
  expect_equal(number_depression_items, 50)
  expect_equal(number_anxiety_items, 5)
  expect_equal(number_somatic_items, 75) 
})

test_that("start n is zero, no constraints", {  
  # define true theta for simulation of responses
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
  
  # define prior covariance matrix
  prior <- diag(number_dimensions) * 20
  
  test_outcome <- with_random_seed(3, test_shadowcat)(true_theta, prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, initital_estimate = rep(.2, number_dimensions), initial_variance = diag(number_dimensions) * 20)
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), c(1.114, -0.018, 1.725))
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3))[1:3],c(.068, .000, .000))
  expect_equal(length(test_outcome$responses), 300)
})

test_that("start n is zero, with constraints", {  
  # define true theta for simulation of responses
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
  
  # define prior covariance matrix
  prior <- diag(number_dimensions) * 20
  
  #create item characteristics and constraints
  characteristics <- data.frame(content = c(rep('depression', number_items / 3), rep('anxiety', number_items / 3), rep('somatic', number_items / 3)))
  constraints <- list(list(name = 'content/depression',
                           op = '><',
                           target = c(50, 75)),
                      list(name = 'content/somatic',
                           op = '><',
                           target = c(75, 90)))
  
  test_outcome <- with_random_seed(3, test_shadowcat)(true_theta, prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, constraints_and_characts = list(characteristics = characteristics, constraints = constraints, initital_estimate = rep(.2, number_dimensions), initial_variance = diag(number_dimensions) * 20))
  indeces_administered <- as.numeric(sapply(names(test_outcome$responses), substring, 5))
  
  number_depression_items <- sum(indeces_administered <= 100)
  number_anxiety_items <- sum(indeces_administered > 100 & indeces_administered <= 200)
  number_somatic_items <- sum(indeces_administered > 200)
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), c(-1.975, .897, 2.496))
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3))[1:3],c(.122, .000, .000))
  expect_equal(length(test_outcome$responses), 130)
  expect_equal(number_depression_items, 50)
  expect_equal(number_anxiety_items, 5)
  expect_equal(number_somatic_items, 75) 
})

context("check stop rule")

test_that("stop rule is number of items", {
  # define true theta for simulation of responses
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
  
  # define prior covariance matrix
  prior <- diag(number_dimensions)
  
  test_outcome <- with_random_seed(2, test_shadowcat)(true_theta, prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary)
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), .405)
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3)), .315)
  expect_equal(names(test_outcome$responses), str_c("item", c(10, 30, 48, 7, 24, 5, 17, 6, 45, 47)))
  expect_equal(unname(unlist(test_outcome$responses)), c(1, 0, 1, 1, 0, 1, 0, 1, 0, 1))
})

test_that("stop rule is variance", {
  # define true theta for simulation of responses
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
  
  # define prior covariance matrix
  prior <- diag(number_dimensions)
  
  test_outcome <- with_random_seed(2, test_shadowcat)(true_theta, prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary)
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), 0.329)
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3)), 0.47)
  expect_equal(names(test_outcome$responses), str_c("item", c(10, 30, 48, 7, 24, 5, 17)))
  expect_equal(unname(unlist(test_outcome$responses)), c(1, 0, 1, 1, 0, 1, 0))
})

test_that("stop rule is variance and minimum number of items is taken into account", {
  # define true theta for simulation of responses
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
  
  # define prior covariance matrix
  prior <- diag(number_dimensions)
  
  test_outcome <- with_random_seed(2, test_shadowcat)(true_theta, prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary)
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), 0.405)
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3)), .315)
  expect_equal(names(test_outcome$responses), str_c("item", c(10, 30, 48, 7, 24, 5, 17, 6, 45, 47)))
  expect_equal(unname(unlist(test_outcome$responses)), c(1, 0, 1, 1, 0, 1, 0, 1, 0, 1))
})

test_that("stop rule is cutoff", {  
  # define true theta for simulation of responses
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
  stop_test <- list(max_n = 300, cutoffs = with_random_seed(2, matrix)(runif(903, 1, 2), ncol = 3))
  estimator <- 'maximum_aposteriori'
  information_summary <- 'posterior_determinant'
  
  # define prior covariance matrix
  prior <- diag(number_dimensions) * 20
  
  test_outcome <- with_random_seed(3, test_shadowcat)(true_theta, prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary)
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), c(.105, -.110, -1.794))
  expect_equal(diag(round(attr(test_outcome$estimate, "variance"), 3)),c(.234, .233, .351))
  expect_equal(length(test_outcome$responses), 32)
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
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(type = 'random_by_dimension', n_by_dimension = 3, n = 9)
  stop_test <- list(max_n = 300, cutoffs = with_random_seed(2, matrix)(runif(903, 1, 2), ncol = 3))
  estimator <- 'maximum_aposteriori'
  information_summary <- 'posterior_determinant'
  
  # define prior covariance matrix
  prior <- diag(number_dimensions) * 20
  estimate <- .3
  attr(estimate, "variance") <- .5
  estimate_without_attr <- .3
  
  error_message_estimate <- shadowcat(responses = NULL, estimate = NULL, model, alpha, beta, start_items, stop_test, estimator, information_summary, prior, guessing, eta)
  error_message_variance <- shadowcat(responses = NULL, estimate_without_attr, model, alpha, beta, start_items, stop_test, estimator, information_summary, prior, guessing, eta)
  error_message_model <- shadowcat(responses = NULL, estimate, model = NULL, alpha, beta, start_items, stop_test, estimator, information_summary, prior, guessing, eta)
  error_message_alpha <- shadowcat(responses = NULL, estimate, model, alpha = NULL, beta, start_items, stop_test, estimator, information_summary, prior, guessing, eta)
  error_message_start_items <- shadowcat(responses = NULL, estimate, model, alpha, beta, start_items = NULL, stop_test, estimator, information_summary, prior, guessing, eta)
  error_message_stop_test <- shadowcat(responses = NULL, estimate, model, alpha, beta, start_items, stop_test = NULL, estimator, information_summary, prior, guessing, eta)
  error_message_estimator <- shadowcat(responses = NULL, estimate, model, alpha, beta, start_items, stop_test, estimator = NULL, information_summary, prior, guessing, eta)
  error_message_information_summary <- shadowcat(responses = NULL, estimate, model, alpha, beta, start_items, stop_test, estimator, information_summary = NULL, prior, guessing, eta)
  error_message_alpha_matrix <- shadowcat(responses = NULL, estimate, model, alpha = 1:10, beta, start_items, stop_test, estimator, information_summary, prior, guessing, eta)
  error_message_alpha_rownames <- shadowcat(responses = NULL, estimate, model, alpha = matrix(1:10), beta, start_items, stop_test, estimator, information_summary, prior, guessing, eta)
  error_message_beta_matrix <- shadowcat(responses = NULL, estimate, model, alpha, beta = 1:10, start_items, stop_test, estimator, information_summary, prior, guessing, eta)
  error_message_beta_rownames <- shadowcat(responses = NULL, estimate, model, alpha, beta = matrix(1:10), start_items, stop_test, estimator, information_summary, prior, guessing, eta)
  error_message_eta_matrix <- shadowcat(responses = NULL, estimate, model, alpha, beta, start_items, stop_test, estimator, information_summary, prior, guessing, eta = 1:10)
  error_message_eta_rownames <- shadowcat(responses = NULL, estimate, model, alpha, beta, start_items, stop_test, estimator, information_summary, prior, guessing, eta = matrix(1:10))
  error_message_guessing_matrix <- shadowcat(responses = NULL, estimate, model, alpha, beta, start_items, stop_test, estimator, information_summary, prior, guessing = 1:10, eta)
  error_message_guessing_rownames <- shadowcat(responses = NULL, estimate, model, alpha, beta, start_items, stop_test, estimator, information_summary, prior, guessing = matrix(1:10), eta)
  error_message_guessing_ncol <- shadowcat(responses = NULL, estimate, model, alpha, beta, start_items, stop_test, estimator, information_summary, prior, guessing = matrix(1:10, ncol = 2, dimnames = list(c("a", "b", "c", "d", "e"), NULL)), eta)
  error_message_unequal_rownames <- shadowcat(responses = NULL, estimate, model, alpha, beta, start_items, stop_test, estimator, information_summary, prior, guessing = matrix(1:300, ncol = 1, dimnames = list(str_c("it", 1:300), NULL)), eta)
  error_message_beta <- shadowcat(responses = NULL, estimate, model, alpha, beta = NULL, start_items, stop_test, estimator, information_summary, prior, guessing, eta)
  error_message_beta_theta <- shadowcat(responses = NULL, estimate, model = "GPCM", alpha, beta = NULL, start_items, stop_test, estimator, information_summary, prior, guessing, eta)
  error_message_beta_theta_mismatch <- shadowcat(responses = NULL, estimate, model = "GPCM", alpha, beta = cbind(beta, beta + .1), start_items, stop_test, estimator, information_summary, prior, guessing, eta = cbind(beta, beta + .1))
  error_message_prior1 <- shadowcat(responses = NULL, estimate, model = "GPCM", alpha, beta, start_items, stop_test, estimator, information_summary = "determinant", prior = NULL, guessing, eta)
  error_message_prior2 <- shadowcat(responses = NULL, estimate, model = "GPCM", alpha, beta, start_items, stop_test, estimator = "maximum_likelihood", information_summary = "posterior_trace", prior = NULL, guessing, eta)
  error_message_max_n <- shadowcat(responses = NULL, estimate, model = "GPCM", alpha, beta, start_items, stop_test = list("variance" = .6), estimator, information_summary, prior, guessing, eta)
  error_message_start_items_0 <- shadowcat(responses = NULL, estimate, model = "GPCM", alpha, beta, start_items = list(n = 0), stop_test, estimator, information_summary = "posterior_expected_kullback_leibler", prior, guessing, eta)
  error_message_cutoffs <- shadowcat(responses = NULL, estimate, model = "GPCM", alpha, beta, start_items, stop_test = list(max_n = 70, cutoffs = 1:70), estimator, information_summary, prior, guessing, eta)
  
  expect_equal(error_message_estimate$errors$estimate, "is missing")
  expect_equal(error_message_variance$errors$variance, "is missing as an attribute of estimate")
  expect_equal(error_message_model$errors$model, "is missing")
  expect_equal(error_message_alpha$errors$alpha, "is missing")
  expect_equal(error_message_start_items$errors$start_items, "is missing")
  expect_equal(error_message_stop_test$errors$stop_test, "is missing")
  expect_equal(error_message_estimator$errors$estimator, "is missing")
  expect_equal(error_message_information_summary$errors$information_summary, "is missing")
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
  expect_equal(error_message_beta$errors$beta, "is missing")
  expect_equal(error_message_beta_theta$errors$beta_and_eta, "are both missing; define at least one of them")
  expect_equal(error_message_beta_theta_mismatch$errors$beta_and_eta, "objects do not match")
  expect_equal(error_message_prior1$errors$prior, "is missing")
  expect_equal(error_message_prior2$errors$prior, "is missing")
  expect_equal(error_message_max_n$errors$stop_test, "contains no max_n")
  expect_equal(error_message_start_items_0$errors$start_items, "requires n > 0 for posterior expected kullback leibler information summary")
  expect_equal(error_message_cutoffs$errors$stop_test, "contains cutoff values in non-matrix format")
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
    lowerbound <- -20
    upperbound <- 20
    
    start_items <- list(type = 'random', n = 0)
    variance_target <- .45^2
    model_vec <- "GPCM"
    estimator_vec <- "expected_aposteriori"
    information_summary_vec <- c("determinant", "posterior_determinant", "trace", "posterior_trace")
    prior <- diag(number_dimensions)
    min_n <- 4
    max_n <- 12
    
    estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, start_items, variance_target, replications_per_unique_condition, number_dimensions, lowerbound = lowerbound, upperbound = upperbound, prior = prior, min_n = min_n, max_n = max_n)
    
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
  
  test_that("simulations with empirical itembank", {
    # itembank.RData is private, these simulations cannot be replicated without itembank.Rdata
    prior <- diag(1)
    model <- "GPCM"
    alpha <- itembank$alpha
    beta <- itembank$beta
    guessing <- NULL
    eta <- NULL
    start_items <- list(type = 'random', n = 0)
    stop_test <- list(target = .45^2, min_n = 4, max_n = 12)
    estimator <- "expected_aposteriori"
    information_summary <- "posterior_determinant"
    lowerbound <- -20
    upperbound <- 20
    
    number_repliations <- 100
    estimated_theta <- numeric(number_repliations)
    sd_estimate <- numeric(number_repliations)
    
    # True theta -2
    true_theta <- -2
    for (s in 1:number_repliations) {
      test_outcome <- with_random_seed(s, test_shadowcat)(true_theta, prior = prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, lowerbound = lowerbound, upperbound = upperbound)$estimate
      estimated_theta[s] <- test_outcome
      sd_estimate[s] <- sqrt(attr(test_outcome, "variance"))
    }
    expect_equal(round(mean(estimated_theta), 3), -1.393)  
    # MAP: -1.314
    expect_equal(round(fivenum(estimated_theta), 3), c(-1.835, -1.835, -1.411, -1.090, -.476))
    expect_equal(round(fivenum(sd_estimate), 3), c(.448, .503, .545, .609, .609))
    
    # True theta 1
    true_theta <- 1
    for (s in 1:number_repliations) {
      test_outcome <- with_random_seed(s, test_shadowcat)(true_theta, prior = prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, lowerbound = lowerbound, upperbound = upperbound)$estimate
      estimated_theta[s] <- test_outcome
      sd_estimate[s] <- sqrt(attr(test_outcome, "variance"))
    }
    expect_equal(round(mean(estimated_theta), 3), .828)  
    # MAP: .836
    expect_equal(round(fivenum(estimated_theta), 3), c(-.191, .569, .866, 1.160, 1.820))
    expect_equal(round(fivenum(sd_estimate), 3), c(.416, .432, .436, .443, .449))
    
    # True theta 3
    true_theta <- 3
    for (s in 1:number_repliations) {
      test_outcome <- with_random_seed(s, test_shadowcat)(true_theta, prior = prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, lowerbound = lowerbound, upperbound = upperbound)$estimate
      estimated_theta[s] <- test_outcome
      sd_estimate[s] <- sqrt(attr(test_outcome, "variance"))
    }
    expect_equal(round(mean(estimated_theta), 3), 2.343)  
    expect_equal(round(fivenum(estimated_theta), 3), c(1.321, 2.045, 2.493, 2.614, 2.929))
    expect_equal(round(fivenum(sd_estimate), 3), c(.418, .443, .448, .449, .491))
  })
  
  test_that("one dimension, one replication per condition, expected_aposteriori", {
    replications_per_unique_condition <- 1
    true_theta_vec <- c(-2, 1, 3)
    number_items_vec <- c(50, 200)
    number_answer_categories_vec <- c(2, 4)
    number_dimensions <- 1
    lowerbound <- -20
    upperbound <- 20
    
    start_items <- list(type = 'random', n = 3)
    variance_target <- .1^2
    model_vec <- c("3PLM", "GRM", "GPCM", "SM")
    estimator_vec <- "expected_aposteriori"
    information_summary_vec <- c("determinant", "posterior_determinant", "trace", "posterior_trace", "posterior_expected_kullback_leibler") 
    prior <- diag(number_dimensions) * 9
    
    estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, start_items, variance_target, replications_per_unique_condition, number_dimensions, lowerbound = lowerbound, upperbound = upperbound, prior = prior)
    
    conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, replications_per_unique_condition, number_dimensions)
    
    estimates_and_conditions <- cbind(t(estimates_and_variance), conditions[, c("true_theta", "number_items", "number_answer_categories", "model", "estimator", "information_summary")])
    
    average_per_true_theta_and_number_items <- aggregate(estimates_and_conditions[,"estimated_theta"], list(estimates_and_conditions[,"true_theta"], estimates_and_conditions[,"number_items"]), "mean")
    sd_per_true_theta_and_number_items <- aggregate(estimates_and_conditions[,"estimated_theta"], list(estimates_and_conditions[,"true_theta"], estimates_and_conditions[,"number_items"]), "sd") 
    
    expect_equal(round(average_per_true_theta_and_number_items[,"x"], 3), c(-2.004, .990, 3.171, -1.964, 1.034, 3.078))
    expect_equal(round(sd_per_true_theta_and_number_items[,"x"], 3), c(.533, .386, .633, .200, .179, .375))
    expect_equal(round(sqrt(fivenum(estimates_and_conditions[,"variance_estimate"])), 3), c(.109, .184, .285, .379, 1.220))    
  })
  
  test_that("one dimension, one replication per condition, maximum_aposteriori", {
    replications_per_unique_condition <- 1
    true_theta_vec <- c(-2, 1, 3)
    number_items_vec <- c(50, 200)
    number_answer_categories_vec <- c(2, 4)
    number_dimensions <- 1
    lowerbound <- -20
    upperbound <- 20
    
    start_items <- list(type = 'random', n = 3)
    variance_target <- .1^2
    model_vec <- c("3PLM","GRM","GPCM","SM")
    estimator_vec <- "maximum_aposteriori"
    information_summary_vec <- c("determinant", "posterior_determinant", "trace", "posterior_trace", "posterior_expected_kullback_leibler")  
    prior <- diag(1) * 9
    
    estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, start_items, variance_target, replications_per_unique_condition, number_dimensions, lowerbound = lowerbound, upperbound = upperbound)
    
    conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, replications_per_unique_condition, number_dimensions)
    
    estimates_and_conditions <- cbind(t(estimates_and_variance), conditions[, c("true_theta", "number_items", "number_answer_categories", "model", "estimator", "information_summary")])
    
    average_per_true_theta_and_number_items <- aggregate(estimates_and_conditions[,"estimated_theta"], list(estimates_and_conditions[,"true_theta"], estimates_and_conditions[,"number_items"]), "mean")
    sd_per_true_theta_and_number_items <- aggregate(estimates_and_conditions[,"estimated_theta"], list(estimates_and_conditions[,"true_theta"], estimates_and_conditions[,"number_items"]), "sd") 
    
    expect_equal(round(average_per_true_theta_and_number_items[,"x"], 3), c(-2.047, .985, 3.175, -1.957, 1.023, 3.094))
    expect_equal(round(sd_per_true_theta_and_number_items[,"x"], 3), c(.580, .404, .647, .200, .167, .403))
    expect_equal(round(sqrt(fivenum(estimates_and_conditions[,"variance_estimate"])), 3), c(.109, .184, .288, .392, 1.232))
  })
  
  test_that("one dimension, one replication per condition, maximum_likelihood", {
    # maximum_likelihood and PEKL do not go well together; makes sense to me, I think maximum_likelihood should be combined with D or A information summary
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
    information_summary_vec <- c("determinant", "posterior_determinant", "trace", "posterior_trace")  # "posterior_expected_kullback_leibler"
    
    estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, start_items, variance_target, replications_per_unique_condition, number_dimensions, lowerbound = lowerbound, upperbound = upperbound)
    
    conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, replications_per_unique_condition, number_dimensions)
    
    estimates_and_conditions <- cbind(t(estimates_and_variance), conditions[, c("true_theta", "number_items", "number_answer_categories", "model", "estimator", "information_summary")])
    
    average_per_true_theta_and_number_items <- aggregate(estimates_and_conditions[,"estimated_theta"], list(estimates_and_conditions[,"true_theta"], estimates_and_conditions[,"number_items"]), "mean")
    sd_per_true_theta_and_number_items <- aggregate(estimates_and_conditions[,"estimated_theta"], list(estimates_and_conditions[,"true_theta"], estimates_and_conditions[,"number_items"]), "sd") 
    
    expect_equal(round(average_per_true_theta_and_number_items[,"x"], 3), c(-1.994, .989, 3.172, -1.984, 1.019, 3.124))
    expect_equal(round(sd_per_true_theta_and_number_items[,"x"], 3), c(.456, .400, .737, .193, .165, .450))
    expect_equal(round(sqrt(fivenum(estimates_and_conditions[,"variance_estimate"])), 3), c(.110, .183, .285, .394, 1.764))
  })
  
  
  test_that("one dimension, estimator maximum_likelihood, 100 replications per condition", {
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
  
test_that("one dimension, estimator maximum_aposteriori, no constraints on item selection, 100 replications per condition", {
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

test_that("one dimension, estimator expected_aposteriori, no constraints on item selection, 100 replications per condition", {
  replications_per_unique_condition <- 100
  true_theta_vec <- c(-2, 1)
  number_items_vec <- c(100, 300)
  number_answer_categories_vec <- c(2, 4)
  number_dimensions <- 1
  
  start_items <- list(type = 'random', n = 3)
  variance_target <- .1^2
  model_vec <- c("3PLM","GRM","GPCM","SM")
  estimator_vec <- "expected_aposteriori"
  information_summary_vec <- c("determinant", "posterior_determinant", "trace", "posterior_trace", "posterior_expected_kullback_leibler")
  
  estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, start_items, variance_target, replications_per_unique_condition, number_dimensions)
  
  conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, replications_per_unique_condition, number_dimensions)
  condition_vector <- sort(rep(1:(ncol(estimates_and_variance)/replications_per_unique_condition), replications_per_unique_condition))
  
  estimates_and_conditions <- cbind(t(estimates_and_variance), conditions[, c("true_theta", "replication", "number_items", "number_answer_categories", "model", "estimator", "information_summary")], condition_vector)
  
  average_per_condition_true_minus2 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == -2), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == -2)]), "mean")
  sd_per_condition_true_minus2 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == -2), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == -2)]), "sd")
  
  average_per_condition_true_1 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 1), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 1)]), "mean")
  sd_per_condition_true_1 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 1), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 1)]), "sd") 
  
  # five number summary of average theta estimate per condition, with true theta is -2
  expect_equal(round(fivenum(average_per_condition_true_minus2[,"x"]), 3), c(-1.669, -1.537, -1.363, -1.265, -1.170))
  # five number summary of average theta estimate per condition, with true theta is 1
  expect_equal(round(fivenum(average_per_condition_true_1[,"x"]), 3), c(1.286, 1.388, 1.537, 1.666, 1.757))
  
  # five number summary of observed sd of the theta estimates within each condition, with true theta is -2
  expect_equal(round(fivenum(sd_per_condition_true_minus2[,"x"]), 3), c(.305, .355, .385, .407, .452))  
  # five number summary of observed sd of the theta estimates within each condition, with true theta is 1
  expect_equal(round(fivenum(sd_per_condition_true_1[,"x"]), 3), c(.262, .349, .425, .486, .551))  
  
  # five number summary of reported sd of the theta estimate within each condition where max number of items is 50 and 100, respectively
  expect_equal(round(sqrt(fivenum(estimates_and_conditions[which(estimates_and_conditions[,"number_items"] == 50), "variance_estimate"])), 3), c(.044, .098, .098, .099, .129))
  expect_equal(round(sqrt(fivenum(estimates_and_conditions[which(estimates_and_conditions[,"number_items"] == 100), "variance_estimate"])), 3), c(.054, .097, .098, .099, .100))  
})

test_that("three dimensions, maximum_likelihood, information summary D, PD, A, and PA, no constraints on item selection, 100 replications per condition", {
  replications_per_unique_condition <- 100 
  true_theta_vec <- c(-2, 1, 2)
  number_items_vec <- c(300)
  number_answer_categories_vec <- c(2, 4)
  number_dimensions <- 3
  
  start_items <- list(type = 'random_by_dimension', n_by_dimension = 3, n = 9)
  variance_target <- .1^2
  model_vec <- c("3PLM","GRM","GPCM","SM")
  estimator_vec <- "maximum_likelihood"
  information_summary_vec <- c("determinant", "posterior_determinant", "trace", "posterior_trace") 
  
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

test_that("three dimensions, maximum_aposteriori, no constraints on item selection, 100 replications per condition", {
  # To be run yet
  replications_per_unique_condition <- 100 
  true_theta_vec <- c(-2, 1, 2)
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

test_that("three dimensions, expected_aposteriori, no constraints on item selection, 100 replications per condition", {
  # To be run yet
  replications_per_unique_condition <- 100 
  true_theta_vec <- c(-2, 1, 2)
  number_items_vec <- c(300)
  number_answer_categories_vec <- c(2, 4)
  number_dimensions <- 3
  
  start_items <- list(type = 'random_by_dimension', n_by_dimension = 3, n = 9)
  variance_target <- .1^2
  model_vec <- c("3PLM","GRM","GPCM","SM")
  estimator_vec <- "expected_aposteriori"
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
  expect_equal(round(fivenum(average_per_condition_dim1[,"x"]), 3), )
  # five number summary of average theta estimate per condition, dimension 2 (true theta is 1)
  expect_equal(round(fivenum(average_per_condition_dim2[,"x"]), 3), )
  # five number summary of average theta estimate per condition, dimension 3 (true theta is 2)
  expect_equal(round(fivenum(average_per_condition_dim3[,"x"]), 3), )
  
  # five number summary of observed sd of the theta estimates within each condition, for dimension 1, 2, and 3, respectively
  expect_equal(round(fivenum(sd_per_condition_dim1[,"x"]), 3), )
  expect_equal(round(fivenum(sd_per_condition_dim2[,"x"]), 3), )
  expect_equal(round(fivenum(sd_per_condition_dim3[,"x"]), 3), )
  
  # five number summary of reported sd of the theta estimate within each condition, for dimension 1, 2, and 3, respectively
  expect_equal(round(sqrt(fivenum(estimates_and_conditions[, "variance_estimate1"])), 3), )
  expect_equal(round(sqrt(fivenum(estimates_and_conditions[, "variance_estimate5"])), 3), )
  expect_equal(round(sqrt(fivenum(estimates_and_conditions[, "variance_estimate9"])), 3), )
})

  test_that("test prior_var_safe_ml is 100", {
    # run again with new safe_ml code in estimate_latent_trait()
    # run maximum_likelihood only
    replications_per_unique_condition <- 100
    true_theta_vec <- c(-2, 1, 2)
    number_items_vec <- c(300)
    number_answer_categories_vec <- c(2, 4)
    number_dimensions <- 3
    
    start_items <- list(type = 'random_by_dimension', n_by_dimension = 3, n = 9)
    variance_target <- .1^2
    model_vec <- c("3PLM","GRM","GPCM","SM")
    estimator_vec <- c("maximum_likelihood", "maximum_aposteriori")
    information_summary_vec <- c("determinant", "posterior_determinant", "trace", "posterior_trace")
    
    estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, start_items, variance_target, replications_per_unique_condition, number_dimensions, prior_var_safe_ml = 100)
    
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
    expect_equal(round(fivenum(average_per_condition_dim1[,"x"]), 3), c(-2.108, -2.045, -2.017, -1.992, -1.950))
    # five number summary of average theta estimate per condition, dimension 2 (true theta is 1)
    expect_equal(round(fivenum(average_per_condition_dim2[,"x"]), 3), c(.954, .990, 1.011, 1.020, 1.077))
    # five number summary of average theta estimate per condition, dimension 3 (true theta is 2)
    expect_equal(round(fivenum(average_per_condition_dim3[,"x"]), 3), c(1.967, 1.999, 2.016, 2.044, 2.123))
    
    # five number summary of observed sd of the theta estimates within each condition, for dimension 1, 2, and 3, respectively
    expect_equal(round(fivenum(sd_per_condition_dim1[,"x"]), 3), c(.176, .213, .304, .340, .392))
    expect_equal(round(fivenum(sd_per_condition_dim2[,"x"]), 3), c(.157, .195, .244, .264, .332))
    expect_equal(round(fivenum(sd_per_condition_dim3[,"x"]), 3), c(.171, .206, .314, .340, .395))
    
    # five number summary of reported sd of the theta estimate within each condition, for dimension 1, 2, and 3, respectively
    expect_equal(round(sqrt(fivenum(estimates_and_conditions[, "variance_estimate1"])), 3), c(.100, .262, .289, .316, .525))
    expect_equal(round(sqrt(fivenum(estimates_and_conditions[, "variance_estimate5"])), 3), c(.100, .206, .248, .277, .382))
    expect_equal(round(sqrt(fivenum(estimates_and_conditions[, "variance_estimate9"])), 3), c(.100, .210, .289, .316, .522))
    
    # no errors/missings for maximum_aposteriori estimator
    expect_equal(number_na_per_condition[which(conditions[seq(1, 6400, 100), "estimator"] == "maximum_aposteriori"), "x"], rep(0, 32))
    # no errors/missings for maximum_likelihood estimator
    expect_equal(number_na_per_condition[which(conditions[seq(1, 6400, 100), "estimator"] == "maximum_likelihood"), "x"], rep(0, 32))
  })
  
  test_that("items load on all dimensions prior_var_safe_ml is 100", {
    # run again with new safe_ml code in estimate_latent_trait()
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
    # run again with new safe_ml code in estimate_latent_trait()
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
    # run again with new safe_ml code in estimate_latent_trait()
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
    
    estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, start_items, variance_target, replications_per_unique_condition, number_dimensions, constraints_and_characts = list(characteristics = characteristics, constraints = constraints), prior_var_safe_ml = 100, return_administered_item_indeces = TRUE, max_n = max_n)
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
    # run again with new safe_ml code in estimate_latent_trait
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
    
    estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, start_items, variance_target, replications_per_unique_condition, number_dimensions, constraints_and_characts = list(characteristics = characteristics, constraints = constraints), prior_var_safe_ml = 100, return_administered_item_indeces = TRUE, max_n = max_n)
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

