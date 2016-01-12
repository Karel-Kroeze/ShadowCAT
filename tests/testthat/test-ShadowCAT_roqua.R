# only for whithin R:
'
library(testthat)
library(pbapply)
library(stringr) # only for test constraints simulations
'

make_random_seed_exist <- rnorm(1)

#' simulate a testing routine with shadowcat
#' 
#' @param true_theta true theta value or vector
#' @param prior covariance matrix of the (multi variate) normal prior for theta; mean vector is fixed at zero; not used for ML estimator
#' @param model String, one of '3PLM', 'GPCM', 'SM' or 'GRM', for the three-parameter logistic, generalized partial credit, sequential or graded response model respectively.
#' @param alpha Matrix of alpha parameters, one column per dimension, one row per item. Note that so called within-dimensional models still use an alpha matrix, they simply 
#' have only one non-zero loading per item.
#' @param beta Matrix of beta parameters, one column per item step, one row per item. Note that ShadowCAT expects response categories to be sequential,
#' and without gaps. That is, the weight parameter in the GPCM model is assumed to be sequential, and equal to the position of the 'location' of the beta parameter in the Beta matrix.
#' The matrix will have a number of columns equal to the largest number of response categories, items with fewer response categories should be 
#' right-padded with \code{NA}. \code{NA} values between response categories are not allowed, and will lead to errors.
#' @param guessing vector of guessing parameters per item. Optionally used in 3PLM model, ignored for all others.
#' @param eta Matrix of location parameters, optionally used in GPCM model, ignored for all others.
#' @param start_items items that are shown to the patient before adaptive proces starts; one of
#' list(type = 'random', n)
#' list(type = 'fixed', indices, n)
#' list(type = 'random_by_dimension', n_by_dimension, n)
#' where n = total number of initial items, indices = vector of initial item indeces, 
#' n_by_dimension = scalar of number of initial items per dimension, or vector with number of initial items for each dimension
#' If n is 0, only n needs to be defined
#' @param stop_test rule for when to stop providing new items to patient; should be a list of the form
#' list(target = ..., max_n = ..., min_n = ..., cutoffs = ...), 
#' where max_n = test length at which testing should stop (even if target has not been reached yet in case of variance stopping rule), 
#' target = vector of maximum acceptable variances per dimension; if target = NULL, only max_n is taken into account,
#' min_n = minimum test length; NULL means no mimimum test length,
#' cutoffs = matrix containing cut off values per dimension (columns) and test iteration (rows). First row contains cut off values for when no items have been
#' administered yet, second row for when one item has been administered, etc. If estimate + 3SE < cutoff for each dimension at certain iteration, test stops; 
#' NULL means no cut off values
#' @param estimator type of estimator to be used, one of "MAP" (Maximum a posteriori estimation) or "ML" (maximum likelihood); 
#' "EAP" (Expected A Posteriori Estimation) is currently not working due to problems with the MultiGHQuad package
#' @param information_summary called "objective" by Kroeze; how to summarize information; one of
#' "D" = determinant: compute determinant(info_sofar_QxQ + info_QxQ_k) for each yet available item k
#' "PD" = posterior determinant: compute determinant(info_sofar_QxQ_plus_prior + info_QxQ_k) for each yet available item k
#' "A" = trace: compute trace((info_sofar_QxQ + info_QxQ_k) for each yet available item k
#' "PA" = posterior trace: compute trace(info_sofar_QxQ_plus_prior + info_QxQ_k) for each yet available item k
#' "PEKL" = compute Posterior expected Kullback-Leibler Information
#' @param item_selection selection criterion; one of "MI" (maximum information) or "Shadow" (maximum information and take constraints into account)
#' @param constraints_and_characts list with constraints and characteristics: constraints_and_characts = list(constraints = ..., characteristics = ...)
#' constraints should be specified as a list of constraints, each constraint is a list with three named values;
#' name: the column name of the characteristic this constraint applies to. For categorical characteristics the level should be specified as name/value.
#' op: the logical operator to be used. Valid options are "<", "=", ">" and "><".
#' target: the target value, numeric. If the operator is "><", this should be a length two vector in between which the target should fall.
#' characteristics should be a data.frame with characteristics, one row per item, one column per characteristic.
#' See constraints_lp_format() for details
#' @param lower_bound vector with lower bounds for theta per dimension; estimated theta values smaller than the lowerbound values are truncated to the lowerbound values 
#' @param upper_bound vector with upper bounds for theta per dimension; estimated theta values larger than the upperbound values are truncated to the upperbound values
#' @param prior_var_safe_ml if not NULL, EAP estimate with prior variance equal to prior_var_safe_ml (scalar or vector) is computed instead of ML/MAP, if ML/MAP estimate fails.
#' @param initital_estimate vector containing the initial theta estimates (starting values)
#' @param initial_variance matrix containing the initial covariance matrix (staring values)
#' @return
test_shadowcat_roqua <- function(true_theta, prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, item_selection = "MI", constraints_and_characts = NULL, lowerbound = rep(-3, ncol(alpha)), upperbound = rep(3, ncol(alpha)), prior_var_safe_nlm = NULL, initital_estimate = rep(0, ncol(alpha)), initial_variance = prior) {
  new_response <- NULL
  attr(initital_estimate, 'variance') <- initial_variance
  next_item_and_test_outcome <- shadowcat_roqua(new_response, estimate = initital_estimate , responses = numeric(0), administered = numeric(0), available = 1:nrow(beta), model, alpha, beta, start_items, stop_test, estimator, information_summary, prior, guessing, eta, item_selection, constraints_and_characts, lowerbound, upperbound, prior_var_safe_nlm)

  while (next_item_and_test_outcome$index_new_item != "stop_test") {
    new_response <- simulate_answer(true_theta, model, ncol(alpha), estimator, alpha, beta, guessing, ncol(beta), next_item_and_test_outcome$index_new_item)
    next_item_and_test_outcome <- shadowcat_roqua(new_response, next_item_and_test_outcome$estimate, next_item_and_test_outcome$responses, next_item_and_test_outcome$administered, next_item_and_test_outcome$available, model, alpha, beta, start_items, stop_test, estimator, information_summary, prior,  guessing, eta, item_selection, constraints_and_characts, lowerbound, upperbound, prior_var_safe_nlm)  
  }
  
  next_item_and_test_outcome
}

#' get conditions for simulation
#' 
#' @param true_theta_vec vector containing true theta values. If number dimensions is 1, simulations are performed for each true theta value in true_theta_vec;
#' if number dimensions is larger than 1, true_theta_vec is interpreted as containing the true thetas for each dimension
#' @param number_items_vec vector containing conditions for number of test bank items to be simulated. Simulations are performed for each value in number_items_vec.
#' If item_selection is "Shadow", number_items_vec can only have length 1
#' @param number_answer_categories_vec vector containing conditions for the number of answer categories to be simulated. 
#' Simulations are performed for each value in number_answer_categories_vec
#' @param model_vec vector containing the conditions for the model to be used. Simulations are performed for each model in model_vec. Model options are
#' '3PLM', 'GPCM', 'SM' and 'GRM'
#' @param estimator_vec vector containing the conditions for the estimator to be used. Simulations are performed for each estimator in estimator_vec.
#' Options are "ML", "MAP", and "EAP"
#' @param information_summary_vec vector containing the conditions for the information_summary to be used. Simulations are performed for each model summary in estimator_vec,
#' Options are "D", "PD", "A", "PA", and "PEKL"
#' @param item_selection selection criterion; one of "MI" (maximum information) or "Shadow" (maximum information and take constraints into account).
#' In case of "Shadow", number_items_vec should have length 1 
#' @param iterations_per_unique_condition number of iterations to be performed within each unique condition
#' @param number_dimensions the number of dimensions of the model (either 1 or the length of true_theta_vec)
get_conditions <- function(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection, iterations_per_unique_condition, number_dimensions) {
  if (number_dimensions == 1) {
    conditions <- expand.grid(1:iterations_per_unique_condition, true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection)
    colnames(conditions) <- c("iteration", "true_theta", "number_items", "number_answer_categories", "model", "estimator", "information_summary", "item_selection")
  }
  else {
    conditions <- expand.grid(1:iterations_per_unique_condition, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection)
    colnames(conditions) <- c("iteration", "number_items", "number_answer_categories", "model", "estimator", "information_summary", "item_selection") 
  }
  conditions
}

# item_selection can be "MI" or "Shadow". In case of "Shadow", constraints_and_characts should be defined, and number_items_vec can only have length 1
# max_n can only have length 1; if null, max_n is set to the number of items (which may differ accross conditions)

#' simulate testing routines with shadowcat, for several conditions
#' 
#' @param true_theta_vec vector containing true theta values. If number dimensions is 1, simulations are performed for each true theta value in true_theta_vec;
#' if number dimensions is larger than 1, true_theta_vec is interpreted as containing the true thetas for each dimension
#' @param number_items_vec vector containing conditions for number of test bank items to be simulated. Simulations are performed for each value in number_items_vec.
#' If item_selection is "Shadow", number_items_vec can only have length 1
#' @param number_answer_categories_vec vector containing conditions for the number of answer categories to be simulated. 
#' Simulations are performed for each value in number_answer_categories_vec
#' @param model_vec vector containing the conditions for the model to be used. Simulations are performed for each model in model_vec. Model options are
#' '3PLM', 'GPCM', 'SM' and 'GRM'
#' @param estimator_vec vector containing the conditions for the estimator to be used. Simulations are performed for each estimator in estimator_vec.
#' Options are "ML", "MAP", and "EAP"
#' @param information_summary_vec vector containing the conditions for the information_summary to be used. Simulations are performed for each model summary in estimator_vec,
#' Options are "D", "PD", "A", "PA", and "PEKL"
#' @param item_selection selection criterion; one of "MI" (maximum information) or "Shadow" (maximum information and take constraints_and_characts into account).
#' In case of "Shadow", number_items_vec should have length 1
#' @param start_items items that are shown to the patient before adaptive proces starts; one of
#' list(type = 'random', n)
#' list(type = 'fixed', indices, n)
#' list(type = 'random_by_dimension', n_by_dimension, n)
#' where n = total number of initial items, indices = vector of initial item indeces, 
#' n_by_dimension = scalar of number of initial items per dimension, or vector with number of initial items for each dimension
#' If n is 0, only n needs to be defined
#' @param variance_target variance of theta at which testing should stop
#' @param iterations_per_unique_condition number of iterations to be performed within each unique condition
#' @param number_dimensions the number of dimensions of the model (either 1 or the length of true_theta_vec)
#' @param constraints_and_characts list with constraints and characteristics: constraints_and_characts = list(constraints = ..., characteristics = ...)
#' constraints should be specified as a list of constraints, each constraint is a list with three named values;
#' name: the column name of the characteristic this constraint applies to. For categorical characteristics the level should be specified as name/value.
#' op: the logical operator to be used. Valid options are "<", "=", ">" and "><".
#' target: the target value, numeric. If the operator is "><", this should be a length two vector in between which the target should fall.
#' characteristics should be a data.frame with characteristics, one row per item, one column per characteristic.
#' See constraints_lp_format() for details
#' @param guessing vector of guessing parameters per item. Optionally used in 3PLM model, ignored for all others.
#' @param items_load_one_dimension if TRUE, items are simulated which load on one dimension. If FALSE, items are simulated which load on all dimensions
#' @param lower_bound vector with lower bounds for theta per dimension; estimated theta values smaller than the lowerbound values are truncated to the lowerbound values 
#' @param upper_bound vector with upper bounds for theta per dimension; estimated theta values larger than the upperbound values are truncated to the upperbound values
#' @param prior covariance matrix of the (multi variate) normal prior for theta; mean vector is fixed at zero; not used for ML estimator
#' @param prior_var_safe_ml if not NULL, EAP estimate with prior variance equal to prior_var_safe_ml (scalar or vector) is computed instead of ML/MAP, if ML/MAP estimate fails.
#' @param return_administered_item_indeces if TRUE, indeces of administered items are added to the output
#' @param max_n the maxixmum number of items to administer (test stops at this number, even if variance target has not been reached)
#' @param varying_number_item_steps if TRUE, the simulated number of item steps differs over items. In that case, number_answer_categories_vec (number of itemsteps + 1)
#' is considered the maxixmum number of categories
#' @return matrix with in each row (= one condition): named vector containing estimated theta, variance of the estimate, and if return_administered_item_indeces is TRUE, the indeces of the administered items
run_simulation <- function(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection, start_items, variance_target, iterations_per_unique_condition, number_dimensions, constraints_and_characts = NULL, guessing = NULL, items_load_one_dimension = TRUE, lowerbound = rep(-3, number_dimensions), upperbound = rep(3, number_dimensions), prior = diag(number_dimensions) * 20, prior_var_safe_nlm = NULL, return_administered_item_indeces = FALSE, max_n = NULL, varying_number_item_steps = FALSE) {                   
  conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection, iterations_per_unique_condition, number_dimensions)
  
  pbapply::pbsapply(1:nrow(conditions), 
                    FUN = function(condition) {
                      if (is.null(max_n))
                        max_n <- conditions[condition, "number_items"] 
                      stop_test <- list(target = variance_target, max_n = max_n)
                      true_theta <- ( if (number_dimensions == 1) 
                                        conditions[condition, "true_theta"]
                                      else
                                        true_theta_vec )
                      alpha_beta <- simulate_testbank(model = as.character(conditions[condition, "model"]), number_items = conditions[condition, "number_items"], number_dimensions = number_dimensions, number_itemsteps = conditions[condition, "number_answer_categories"] - 1, items_load_one_dimension = items_load_one_dimension, return_testbank_properties = FALSE, varying_number_item_steps = varying_number_item_steps)
                      estimate_theta <- tryCatch(test_shadowcat_roqua(true_theta, prior, as.character(conditions[condition, "model"]), alpha_beta$alpha, alpha_beta$beta, guessing, eta = NULL, start_items, stop_test, as.character(conditions[condition, "estimator"]), as.character(conditions[condition, "information_summary"]), item_selection, constraints_and_characts, lowerbound, upperbound, prior_var_safe_nlm),
                                                 error = function(e) e)

                      if (return_administered_item_indeces)
                        c("estimated_theta" = estimate_theta$estimate,
                          "variance_estimate" = attr(estimate_theta$estimate, "variance"),
                          "items_administered" = estimate_theta$administered)
                      else
                        c("estimated_theta" = estimate_theta$estimate,
                          "variance_estimate" =  attr(estimate_theta$estimate, "variance"))               
                    })
}

context("validate shadowcat_roqua single conditions")

test_that("true theta is 2, estimator is MAP", {
  # define true theta for simulation of responses
  true_theta <- 2
  
  # define item characteristics
  number_items <- 100
  number_dimensions <- 1
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(type = 'random', n = 3)
  stop_test <- list(max_n = 100)
  estimator <- 'MAP'
  information_summary <- 'PD'
  item_selection <- 'MI'
  
  # define prior covariance matrix
  prior <- diag(number_dimensions) * 5
  
  test_outcome <- with_random_seed(2, test_shadowcat_roqua)(true_theta, prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary)
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), 1.938)
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3)), .113)
  expect_equal(test_outcome$available, numeric(0))
  expect_equal(length(test_outcome$administered), 100)
})

test_that("true theta is 2, estimator is ML", {
  # define true theta for simulation of responses
  true_theta <- 2
  
  # define item characteristics
  number_items <- 100
  number_dimensions <- 1
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(type = 'random', n = 3)
  stop_test <- list(max_n = 100)
  estimator <- 'ML'
  information_summary <- 'D'
  item_selection <- 'MI'
  
  test_outcome <- with_random_seed(2, test_shadowcat_roqua)(true_theta, prior = NULL, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, initial_variance = diag(1))
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), 2.169)
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3)), .129)
  expect_equal(test_outcome$available, numeric(0))
  expect_equal(length(test_outcome$administered), 100)
})

test_that("true theta is 2, estimator is EAP", {
  # define true theta for simulation of responses
  true_theta <- 2
  
  # define item characteristics
  number_items <- 100
  number_dimensions <- 1
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  
  guessing <- NULL
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(type = 'random', n = 3)
  stop_test <- list(max_n = 100)
  estimator <- 'EAP'
  information_summary <- 'PD'
  item_selection <- 'MI'
  
  prior <- diag(number_dimensions) * 5
  
  test_outcome <- with_random_seed(2, test_shadowcat_roqua)(true_theta, prior = prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, initial_variance = diag(1))
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), 2.913)
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3)), .002)
  expect_equal(test_outcome$available, numeric(0))
  expect_equal(length(test_outcome$administered), 100)
})


test_that("true theta is 1, 0, 2, estimator is MAP", {  
  # define true theta for simulation of responses
  true_theta <- c(1, 0, 2)
  
  # define item characteristics
  number_items <- 300
  number_dimensions <- 3
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  
  guessing <- NULL
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  alpha[1:100,2:3] <- 0
  alpha[101:200,c(1,3)] <- 0
  alpha[201:300,1:2] <- 0
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(type = 'random', n = 3)
  stop_test <- list(max_n = 300)
  estimator <- 'MAP'
  information_summary <- 'PD'
  item_selection <- 'MI'
  
  # define prior covariance matrix
  prior <- diag(number_dimensions) * 20
  
  test_outcome <- with_random_seed(3, test_shadowcat_roqua)(true_theta, prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary)
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), c(.841, -.123, 1.947))
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3))[1:3],c(.064, .000, .000))
  expect_equal(length(test_outcome$administered), 300)
})

test_that("true theta is 1, 0, 2, estimator is ML", {  
  # define true theta for simulation of responses
  true_theta <- c(1, 0, 2)
  
  # define item characteristics
  number_items <- 300
  number_dimensions <- 3
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  
  guessing <- NULL
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  alpha[1:100,2:3] <- 0
  alpha[101:200,c(1,3)] <- 0
  alpha[201:300,1:2] <- 0
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(type = 'random', n = 3)
  stop_test <- list(max_n = 300)
  estimator <- 'ML'
  information_summary <- 'D'
  item_selection <- 'MI'
  
  test_outcome <- with_random_seed(3, test_shadowcat_roqua)(true_theta, prior = NULL, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, initial_variance = diag(1))
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), c(.755, -.070, 2.221))
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3))[1:3],c(.063, .000, .000))
  expect_equal(length(test_outcome$administered), 300)
  
  # defining prior has no effect on outcome
  test_outcome <- with_random_seed(3, test_shadowcat_roqua)(true_theta, prior = diag(number_dimensions) * 2, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, initial_variance = diag(1))
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), c(.755, -.070, 2.221))
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3))[1:3],c(.063, .000, .000))
  expect_equal(length(test_outcome$administered), 300)  
})

test_that("true theta is 1, 0, 2, estimator is EAP", {  
  # define true theta for simulation of responses
  true_theta <- c(1, 0, 2)
  
  # define item characteristics
  number_items <- 300
  number_dimensions <- 3
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  
  guessing <- NULL
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  alpha[1:100,2:3] <- 0
  alpha[101:200,c(1,3)] <- 0
  alpha[201:300,1:2] <- 0
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(type = 'random', n = 3)
  stop_test <- list(max_n = 300)
  estimator <- 'MAP'
  information_summary <- 'PD'
  item_selection <- 'MI'
  
  # define prior covariance matrix
  prior <- diag(number_dimensions) * 20
  
  test_outcome <- with_random_seed(3, test_shadowcat_roqua)(true_theta, prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary)
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), c(.841, -.123, 1.947))
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3))[1:3],c(.064, .000, .000))
  expect_equal(length(test_outcome$administered), 300)
})


test_that("items load on three dimensions", {  
  # define true theta for simulation of responses
  true_theta <- c(1, 0, 2)
  
  # define item characteristics
  number_items <- 300
  number_dimensions <- 3
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  
  guessing <- NULL
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(type = 'random', n = 3)
  stop_test <- list(max_n = 300)
  estimator <- 'MAP'
  information_summary <- 'PD'
  item_selection <- 'MI'
  
  # define prior covariance matrix
  prior <- diag(number_dimensions) * 20
  
  test_outcome <- with_random_seed(3, test_shadowcat_roqua)(true_theta, prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary)
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), c(.454, .719, 1.745))
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3))[1:3],c(.178, -.088, -.081))
  expect_equal(length(test_outcome$administered), 300)
})

test_that("true theta is 2, 2, 2", {  
  # define true theta for simulation of responses
  true_theta <- c(2, 2, 2)
  
  # define item characteristics
  number_items <- 300
  number_dimensions <- 3
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  
  guessing <- NULL
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  alpha[1:100,2:3] <- 0
  alpha[101:200,c(1,3)] <- 0
  alpha[201:300,1:2] <- 0
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(type = 'random', n = 3)
  stop_test <- list(max_n = 300)
  estimator <- 'MAP'
  information_summary <- 'PD'
  item_selection <- 'MI'
  
  # define prior covariance matrix
  prior <- diag(number_dimensions) * 20
  
  test_outcome <- with_random_seed(3, test_shadowcat_roqua)(true_theta, prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary)
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), c(2.287, 1.825, 1.672))
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3))[1:3],c(.11, .000, .000))
  expect_equal(length(test_outcome$administered), 300)
})

test_that("with constraints max_n 260", {  
  # define true theta for simulation of responses
  true_theta <- c(-2, 1, 2)
  
  # define item characteristics
  number_items <- 300
  number_dimensions <- 3
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  
  guessing <- NULL
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  alpha[1:100,2:3] <- 0
  alpha[101:200,c(1,3)] <- 0
  alpha[201:300,1:2] <- 0
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(type = 'random_by_dimension', n_by_dimension = 3, n = 9)
  stop_test <- list(max_n = 260)
  estimator <- 'MAP'
  information_summary <- 'PD'
  item_selection <- 'Shadow'  
  
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
  
  test_outcome <- with_random_seed(3, test_shadowcat_roqua)(true_theta, prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, item_selection =  item_selection, constraints_and_characts = list(characteristics = characteristics, constraints = constraints))
  
  number_depression_items <- sum(test_outcome$administered <= 100)
  number_anxiety_items <- sum(test_outcome$administered > 100 & test_outcome$administered <= 200)
  number_somatic_items <- sum(test_outcome$administered > 200)
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), c(-2.689, .697, 2.260))
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3))[1:3],c(.155, .000, .000))
  expect_equal(length(test_outcome$administered), 260)
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
  
  guessing <- NULL
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  alpha[1:100,2:3] <- 0
  alpha[101:200,c(1,3)] <- 0
  alpha[201:300,1:2] <- 0
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(type = 'random_by_dimension', n_by_dimension = 3, n = 9)
  stop_test <- list(max_n = 130)
  estimator <- 'MAP'
  information_summary <- 'PD'
  item_selection <- 'Shadow'  
  
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
  
  test_outcome <- with_random_seed(3, test_shadowcat_roqua)(true_theta, prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, item_selection =  item_selection, constraints_and_characts = list(characteristics = characteristics, constraints = constraints))
  
  number_depression_items <- sum(test_outcome$administered <= 100)
  number_anxiety_items <- sum(test_outcome$administered > 100 & test_outcome$administered <= 200)
  number_somatic_items <- sum(test_outcome$administered > 200)
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), c(-1.577, .071, 1.906))
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3))[1:3],c(.106, .000, .000))
  expect_equal(length(test_outcome$administered), 130)
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
  
  guessing <- NULL
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  alpha[1:100,2:3] <- 0
  alpha[101:200,c(1,3)] <- 0
  alpha[201:300,1:2] <- 0
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(n = 0)
  stop_test <- list(max_n = 300)
  estimator <- 'MAP'
  information_summary <- 'PD'
  item_selection <- 'MI'
  
  # define prior covariance matrix
  prior <- diag(number_dimensions) * 20
  
  test_outcome <- with_random_seed(3, test_shadowcat_roqua)(true_theta, prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, initital_estimate = rep(.2, number_dimensions), initial_variance = diag(number_dimensions) * 20)
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), c(1.114, -0.018, 1.725))
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3))[1:3],c(.068, .000, .000))
  expect_equal(length(test_outcome$administered), 300)
})

test_that("start n is zero, with constraints", {  
  # define true theta for simulation of responses
  true_theta <- c(-2, 1, 2)
  
  # define item characteristics
  number_items <- 300
  number_dimensions <- 3
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  
  guessing <- NULL
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  alpha[1:100,2:3] <- 0
  alpha[101:200,c(1,3)] <- 0
  alpha[201:300,1:2] <- 0
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(n = 0)
  stop_test <- list(max_n = 130)
  estimator <- 'MAP'
  information_summary <- 'PD'
  item_selection <- 'Shadow'  
  
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
  
  test_outcome <- with_random_seed(3, test_shadowcat_roqua)(true_theta, prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, item_selection =  item_selection, constraints_and_characts = list(characteristics = characteristics, constraints = constraints, initital_estimate = rep(.2, number_dimensions), initial_variance = diag(number_dimensions) * 20))
  
  number_depression_items <- sum(test_outcome$administered <= 100)
  number_anxiety_items <- sum(test_outcome$administered > 100 & test_outcome$administered <= 200)
  number_somatic_items <- sum(test_outcome$administered > 200)
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), c(-1.975, .897, 2.496))
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3))[1:3],c(.122, .000, .000))
  expect_equal(length(test_outcome$administered), 130)
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
  
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(type = 'random', n = 5)
  stop_test <- list(max_n = 10)
  estimator <- 'MAP'
  information_summary <- 'PD'
  item_selection <- 'MI'
  
  # define prior covariance matrix
  prior <- diag(number_dimensions)
  
  test_outcome <- with_random_seed(2, test_shadowcat_roqua)(true_theta, prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary)
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), .405)
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3)), .315)
  expect_equal(test_outcome$available, c(1:4, 8:9, 11:16, 18:23, 25:29, 31:44, 46, 49:50))
  expect_equal(test_outcome$administered, c(10, 30, 48, 7, 24, 5, 17, 6, 45, 47))
  expect_equal(test_outcome$responses, c(1, 0, 1, 1, 0, 1, 0, 1, 0, 1))
})

test_that("stop rule is variance", {
  # define true theta for simulation of responses
  true_theta <- 0
  
  # define item characteristics
  number_items <- 50
  number_dimensions <- 1
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(type = 'random', n = 5)
  stop_test <- list(target = .5, max_n = 50)
  estimator <- 'MAP'
  information_summary <- 'PD'
  item_selection <- 'MI'
  
  # define prior covariance matrix
  prior <- diag(number_dimensions)
  
  test_outcome <- with_random_seed(2, test_shadowcat_roqua)(true_theta, prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary)
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), 0.329)
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3)), 0.47)
  expect_equal(test_outcome$available, c(1:4, 6, 8:9, 11:16, 18:23, 25:29, 31:47, 49:50))
  expect_equal(test_outcome$administered, c(10, 30, 48, 7, 24, 5, 17))
  expect_equal(test_outcome$responses, c(1, 0, 1, 1, 0, 1, 0))
})

test_that("stop rule is variance and minimum number of items is taken into account", {
  # define true theta for simulation of responses
  true_theta <- 0
  
  # define item characteristics
  number_items <- 50
  number_dimensions <- 1
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(type = 'random', n = 5)
  stop_test <- list(target = .5, max_n = 50, min_n = 10)
  estimator <- 'MAP'
  information_summary <- 'PD'
  item_selection <- 'MI'
  
  # define prior covariance matrix
  prior <- diag(number_dimensions)
  
  test_outcome <- with_random_seed(2, test_shadowcat_roqua)(true_theta, prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary)
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), 0.405)
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3)), .315)
  expect_equal(test_outcome$administered, c(10, 30, 48, 7, 24, 5, 17, 6, 45, 47))
  expect_equal(test_outcome$responses, c(1, 0, 1, 1, 0, 1, 0, 1, 0, 1))
})

test_that("stop rule is cutoff", {  
  # define true theta for simulation of responses
  true_theta <- c(.2, 0, -2)
  
  # define item characteristics
  number_items <- 300
  number_dimensions <- 3
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  
  guessing <- NULL
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  alpha[1:100,2:3] <- 0
  alpha[101:200,c(1,3)] <- 0
  alpha[201:300,1:2] <- 0
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(type = 'random_by_dimension', n_by_dimension = 3, n = 9)
  stop_test <- list(max_n = 300, cutoffs = with_random_seed(2, matrix)(runif(903, 1, 2), ncol = 3))
  estimator <- 'MAP'
  information_summary <- 'PD'
  item_selection <- 'MI'
  
  # define prior covariance matrix
  prior <- diag(number_dimensions) * 20
  
  test_outcome <- with_random_seed(3, test_shadowcat_roqua)(true_theta, prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary)
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), c(.105, -.110, -1.794))
  expect_equal(diag(round(attr(test_outcome$estimate, "variance"), 3)),c(.234, .233, .351))
  expect_equal(length(test_outcome$administered), 32)
})


# these simulations take a long time to run, if (FALSE) ensures that they are not each time the tests are run
if (FALSE) {
  context("simulations")
  
  test_that("one dimension, no constraints on item selection, one iteration per condition, EAP", {
    iterations_per_unique_condition <- 1
    true_theta_vec <- c(-2, 1)
    number_items_vec <- c(15, 50)
    number_answer_categories_vec <- c(2, 4)
    number_dimensions <- 1
    lowerbound <- -3
    upperbound <- 3
    
    start_items <- list(type = 'random', n = 3)
    variance_target <- .1^2
    model_vec <- c("3PLM", "GRM", "GPCM", "SM")
    estimator_vec <- "EAP"
    information_summary_vec <- c("D", "PD", "A", "PA", "PEKL") 
    item_selection <- 'MI'
    prior <- diag(number_dimensions) * 100
    
    estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection, start_items, variance_target, iterations_per_unique_condition, number_dimensions, lowerbound = lowerbound, upperbound = upperbound, prior = prior)
    
    conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection, iterations_per_unique_condition, number_dimensions)
    
    estimates_and_conditions <- cbind(t(estimates_and_variance), conditions[, c("true_theta", "number_items", "number_answer_categories", "model", "estimator", "information_summary")])
    
    average_per_true_theta_and_number_items <- aggregate(estimates_and_conditions[,"estimated_theta"], list(estimates_and_conditions[,"true_theta"], estimates_and_conditions[,"number_items"]), "mean")
    sd_per_true_theta_and_number_items <- aggregate(estimates_and_conditions[,"estimated_theta"], list(estimates_and_conditions[,"true_theta"], estimates_and_conditions[,"number_items"]), "sd") 
    
    expect_equal(round(average_per_true_theta_and_number_items[,"x"], 3), c(-1.668, 1.910, -1.535, 1.917))
    expect_equal(round(sd_per_true_theta_and_number_items[,"x"], 3), c(.808, .679, .586, .426))
    expect_equal(round(fivenum(estimates_and_conditions[,"variance_estimate"]), 3), c(.000, .008, .010, .044, .228))    
  })
  
  test_that("one dimension, no constraints on item selection, one iteration per condition, MAP", {
    iterations_per_unique_condition <- 1
    true_theta_vec <- c(-2, 1)
    number_items_vec <- c(15, 50)
    number_answer_categories_vec <- c(2, 4)
    number_dimensions <- 1
    lowerbound <- -3
    upperbound <- 3
    
    start_items <- list(type = 'random', n = 3)
    variance_target <- .1^2
    model_vec <- c("3PLM","GRM","GPCM","SM")
    estimator_vec <- "MAP"
    information_summary_vec <- c("D", "PD", "A", "PA", "PEKL") 
    item_selection <- 'MI'
    
    estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection, start_items, variance_target, iterations_per_unique_condition, number_dimensions, lowerbound = lowerbound, upperbound = upperbound)
    
    conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection, iterations_per_unique_condition, number_dimensions)
    
    estimates_and_conditions <- cbind(t(estimates_and_variance), conditions[, c("true_theta", "number_items", "number_answer_categories", "model", "estimator", "information_summary")])
    
    average_per_true_theta_and_number_items <- aggregate(estimates_and_conditions[,"estimated_theta"], list(estimates_and_conditions[,"true_theta"], estimates_and_conditions[,"number_items"]), "mean")
    sd_per_true_theta_and_number_items <- aggregate(estimates_and_conditions[,"estimated_theta"], list(estimates_and_conditions[,"true_theta"], estimates_and_conditions[,"number_items"]), "sd") 
    
    expect_equal(round(average_per_true_theta_and_number_items[,"x"], 3), c(-1.884, 1.145, -2.046, .996))
    expect_equal(round(sd_per_true_theta_and_number_items[,"x"], 3), c(.606, .696, .378, .386))
    expect_equal(fivenum(round(estimates_and_conditions[,"variance_estimate"], 3)), c(.049, .122, .212, .383, 4.547))
  })
  
  test_that("one dimension, no constraints on item selection, one iteration per condition, ML", {
    # ML and PEKL do not go well together; makes sense to me, I think ML should be combined with D or A information summary
    iterations_per_unique_condition <- 1
    true_theta_vec <- c(-2, 1)
    number_items_vec <- c(15, 50)
    number_answer_categories_vec <- c(2, 4)
    number_dimensions <- 1
    lowerbound <- -3
    upperbound <- 3
    
    start_items <- list(type = 'random', n = 3)
    variance_target <- .1^2
    model_vec <- c("3PLM","GRM","GPCM","SM")
    estimator_vec <- "ML"
    information_summary_vec <- c("D", "PD", "A", "PA") # PEKL
    item_selection <- 'MI'
    
    estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection, start_items, variance_target, iterations_per_unique_condition, number_dimensions, lowerbound = lowerbound, upperbound = upperbound)
    
    conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection, iterations_per_unique_condition, number_dimensions)
    
    estimates_and_conditions <- cbind(t(estimates_and_variance), conditions[, c("true_theta", "number_items", "number_answer_categories", "model", "estimator", "information_summary")])
    
    average_per_true_theta_and_number_items <- aggregate(estimates_and_conditions[,"estimated_theta"], list(estimates_and_conditions[,"true_theta"], estimates_and_conditions[,"number_items"]), "mean")
    sd_per_true_theta_and_number_items <- aggregate(estimates_and_conditions[,"estimated_theta"], list(estimates_and_conditions[,"true_theta"], estimates_and_conditions[,"number_items"]), "sd") 
    
    expect_equal(round(average_per_true_theta_and_number_items[,"x"], 3), c(-1.920, 1.231, -2.074, 1.012))
                 expect_equal(round(sd_per_true_theta_and_number_items[,"x"], 3), c(.606, .705, .404, .401))
                 expect_equal(fivenum(round(estimates_and_conditions[,"variance_estimate"], 3)), c(4.9e-02, 1.24e-01, 2.085e-01, 3.915e-01, 111839145))
                 # This extremely large variance corresponds to an estimated theta at the lowerbound of -3; setting the lowerbound to -5 does not help
  })
  
  
  test_that("one dimension, estimator ML, information summary D, PD, A, and PA, no constraints on item selection, 100 iterations per condition", {
    # run one more time to check
    iterations_per_unique_condition <- 100
    true_theta_vec <- c(-2, 1)
    number_items_vec <- c(50, 100)
    number_answer_categories_vec <- c(2, 4)
    number_dimensions <- 1
    
    start_items <- list(type = 'random', n = 3)
    variance_target <- .1^2
    model_vec <- c("3PLM","GRM","GPCM","SM")
    estimator_vec <- "ML"
    information_summary_vec <- c("D", "PD", "A", "PA") # PEKL 
    item_selection <- 'MI'
    
    estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection, start_items, variance_target, iterations_per_unique_condition, number_dimensions)
    estimates_and_variance_without_errors <- sapply(estimates_and_variance, 
                                                    FUN = function(x) { 
                                                      if (is.null(x)) 
                                                        matrix(rep(NA, 2), nrow = 1, dimnames = list(NULL, names(estimates_and_variance[[1]])))
                                                      else 
                                                        x } ) 
    
    conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection, iterations_per_unique_condition, number_dimensions)
    
    condition_vector <- sort(rep(1:(ncol(estimates_and_variance_without_errors)/iterations_per_unique_condition), iterations_per_unique_condition))
    
    estimates_and_conditions <- cbind(t(estimates_and_variance_without_errors)[, c("estimated_theta", "variance_estimate")], 
                                      conditions[, c("true_theta", "iteration", "number_items", "number_answer_categories", "model", "estimator", "information_summary")], 
                                      condition_vector)
    
    average_per_condition_true_minus2 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == -2), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == -2)]), "mean")
    sd_per_condition_true_minus2 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == -2), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == -2)]), "sd")
    
    average_per_condition_true_1 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 1), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 1)]), "mean")
    sd_per_condition_true_1 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 1), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 1)]), "sd") 
    
    number_na_per_condition <- aggregate(estimates_and_conditions[, "estimated_theta"], list(condition_vector), FUN = function(x) { sum(is.na(x)) })

    
    # five number summary of average theta estimate per condition, with true theta is -2
    expect_equal(round(fivenum(average_per_condition_true_minus2[,"x"]), 3), c(-2.140, -2.056, -2.025, -2.001, -1.937))
    # five number summary of average theta estimate per condition, with true theta is 1
    expect_equal(round(fivenum(average_per_condition_true_1[,"x"]), 3), c(.939, .989, 1.004, 1.022, 1.078))
    
    # five number summary of observed sd of the theta estimates within each condition, with true theta is -2
    expect_equal(round(fivenum(sd_per_condition_true_minus2[,"x"]), 3), c(.173, .289, .321, .440, .548))  
    # five number summary of observed sd of the theta estimates within each condition, with true theta is 1
    expect_equal(round(fivenum(sd_per_condition_true_1[,"x"]), 3), c(.162, .228, .259, .352, .455))  
    
    # five number summary of reported sd of the theta estimate within each condition where max number of items is 50 and 100, respectively
    expect_equal(round(sqrt(fivenum(estimates_and_conditions[which(estimates_and_conditions[,"number_items"] == 50), "variance_estimate"])), 3), c(.188, .284, .352, .420, 1.585))
    expect_equal(round(sqrt(fivenum(estimates_and_conditions[which(estimates_and_conditions[,"number_items"] == 100), "variance_estimate"])), 3), c(.142, .200, .249, .299, .569))

    expect_equal(number_na_per_condition[ ,"x"], c(rep(0, 126), 1, 0))
  })
  
test_that("one dimension, estimator MAP, no constraints on item selection, 100 iterations per condition", {
  iterations_per_unique_condition <- 100
  true_theta_vec <- c(-2, 1)
  number_items_vec <- c(50, 100)
  number_answer_categories_vec <- c(2, 4)
  number_dimensions <- 1
  
  start_items <- list(type = 'random', n = 3)
  variance_target <- .1^2
  model_vec <- c("3PLM","GRM","GPCM","SM")
  estimator_vec <- "MAP"
  information_summary_vec <- c("D", "PD", "A", "PA", "PEKL")
  item_selection <- 'MI'
  
  estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection, start_items, variance_target, iterations_per_unique_condition, number_dimensions)
  
  conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection, iterations_per_unique_condition, number_dimensions)
  condition_vector <- sort(rep(1:(ncol(estimates_and_variance)/iterations_per_unique_condition), iterations_per_unique_condition))
  
  estimates_and_conditions <- cbind(t(estimates_and_variance), conditions[, c("true_theta", "iteration", "number_items", "number_answer_categories", "model", "estimator", "information_summary")], condition_vector)
  
  average_per_condition_true_minus2 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == -2), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == -2)]), "mean")
  sd_per_condition_true_minus2 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == -2), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == -2)]), "sd")
  
  average_per_condition_true_1 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 1), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 1)]), "mean")
  sd_per_condition_true_1 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 1), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 1)]), "sd") 
  
  # five number summary of average theta estimate per condition, with true theta is -2
  expect_equal(round(fivenum(average_per_condition_true_minus2[,"x"]), 3), c(-2.126, -2.036, -2.016, -1.997, -1.934))
  # five number summary of average theta estimate per condition, with true theta is 1
  expect_equal(round(fivenum(average_per_condition_true_1[,"x"]), 3), c(.934, .986, 1.002, 1.018, 1.084))
  
  # five number summary of observed sd of the theta estimates within each condition, with true theta is -2
  expect_equal(round(fivenum(sd_per_condition_true_minus2[,"x"]), 3), c(.179, .278, .320, .439, .542))  
  # five number summary of observed sd of the theta estimates within each condition, with true theta is 1
  expect_equal(round(fivenum(sd_per_condition_true_1[,"x"]), 3), c(.154, .226, .256, .346, .434))  
  
  # five number summary of reported sd of the theta estimate within each condition where max number of items is 50 and 100, respectively
  expect_equal(round(sqrt(fivenum(estimates_and_conditions[which(estimates_and_conditions[,"number_items"] == 50), "variance_estimate"])), 3), c(.188, .283, .350, .418, 1.307))
  expect_equal(round(sqrt(fivenum(estimates_and_conditions[which(estimates_and_conditions[,"number_items"] == 100), "variance_estimate"])), 3), c(.142, .200, .248, .298, .543))  
})

test_that("one dimension, estimator EAP, no constraints on item selection, 100 iterations per condition", {
  iterations_per_unique_condition <- 100
  true_theta_vec <- c(-2, 1)
  number_items_vec <- c(50, 100)
  number_answer_categories_vec <- c(2, 4)
  number_dimensions <- 1
  
  start_items <- list(type = 'random', n = 3)
  variance_target <- .1^2
  model_vec <- c("3PLM","GRM","GPCM","SM")
  estimator_vec <- "EAP"
  information_summary_vec <- c("D", "PD", "A", "PA", "PEKL")
  item_selection <- 'MI'
  
  estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection, start_items, variance_target, iterations_per_unique_condition, number_dimensions)
  
  conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection, iterations_per_unique_condition, number_dimensions)
  condition_vector <- sort(rep(1:(ncol(estimates_and_variance)/iterations_per_unique_condition), iterations_per_unique_condition))
  
  estimates_and_conditions <- cbind(t(estimates_and_variance), conditions[, c("true_theta", "iteration", "number_items", "number_answer_categories", "model", "estimator", "information_summary")], condition_vector)
  
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

test_that("three dimensions, ML, information summary D, PD, A, and PA, no constraints on item selection, 100 iterations per condition", {
  iterations_per_unique_condition <- 100 
  true_theta_vec <- c(-2, 1, 2)
  number_items_vec <- c(300)
  number_answer_categories_vec <- c(2, 4)
  number_dimensions <- 3
  
  start_items <- list(type = 'random_by_dimension', n_by_dimension = 3, n = 9)
  variance_target <- .1^2
  model_vec <- c("3PLM","GRM","GPCM","SM")
  estimator_vec <- "ML"
  information_summary_vec <- c("D", "PD", "A", "PA") # PEKL 
  item_selection <- 'MI'
  
  estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection, start_items, variance_target, iterations_per_unique_condition, number_dimensions)
  estimates_and_variance_without_errors <- sapply(estimates_and_variance, 
                                                  FUN = function(x) { 
                                                    if (is.null(x)) 
                                                      matrix(rep(NA, number_dimensions + number_dimensions^2), nrow = 1, dimnames = list(NULL, names(estimates_and_variance[[1]])))
                                                    else 
                                                      x } ) 
  
  conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection, iterations_per_unique_condition, number_dimensions)
  
  condition_vector <- sort(rep(1:(ncol(estimates_and_variance_without_errors)/iterations_per_unique_condition), iterations_per_unique_condition))
  
  estimates_and_conditions <- cbind(t(estimates_and_variance_without_errors)[, c("estimated_theta1", "estimated_theta2", "estimated_theta3", "variance_estimate1", "variance_estimate5", "variance_estimate9")], 
                                    conditions[, c("iteration", "number_items", "number_answer_categories", "model", "estimator", "information_summary")], 
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
  
  # some errors/missings for ML estimator because prior_var_safe_nlm = NULL
  expect_equal(number_na_per_condition[, "x"], 
               c(rep(0, 18), 1, 0, 0, 0, 2, rep(0, 7), 1, 0))
})

test_that("three dimensions, MAP, no constraints on item selection, 100 iterations per condition", {
  # To be run yet
  iterations_per_unique_condition <- 100 
  true_theta_vec <- c(-2, 1, 2)
  number_items_vec <- c(300)
  number_answer_categories_vec <- c(2, 4)
  number_dimensions <- 3
  
  start_items <- list(type = 'random_by_dimension', n_by_dimension = 3, n = 9)
  variance_target <- .1^2
  model_vec <- c("3PLM","GRM","GPCM","SM")
  estimator_vec <- "MAP"
  information_summary_vec <- c("D", "PD", "A", "PA", "PEKL")
  item_selection <- 'MI'
  
  estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection, start_items, variance_target, iterations_per_unique_condition, number_dimensions)
 
  conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection, iterations_per_unique_condition, number_dimensions)
  
  condition_vector <- sort(rep(1:(ncol(estimates_and_variance)/iterations_per_unique_condition), iterations_per_unique_condition))
  
  estimates_and_conditions <- cbind(t(estimates_and_variance)[, c("estimated_theta1", "estimated_theta2", "estimated_theta3", "variance_estimate1", "variance_estimate5", "variance_estimate9")], 
                                    conditions[, c("iteration", "number_items", "number_answer_categories", "model", "estimator", "information_summary")], 
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

test_that("three dimensions, EAP, no constraints on item selection, 100 iterations per condition", {
  # To be run yet
  iterations_per_unique_condition <- 100 
  true_theta_vec <- c(-2, 1, 2)
  number_items_vec <- c(300)
  number_answer_categories_vec <- c(2, 4)
  number_dimensions <- 3
  
  start_items <- list(type = 'random_by_dimension', n_by_dimension = 3, n = 9)
  variance_target <- .1^2
  model_vec <- c("3PLM","GRM","GPCM","SM")
  estimator_vec <- "AEP"
  information_summary_vec <- c("D", "PD", "A", "PA", "PEKL") 
  item_selection <- 'MI'
  
  estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection, start_items, variance_target, iterations_per_unique_condition, number_dimensions)
  estimates_and_variance_without_errors <- sapply(estimates_and_variance, 
                                                  FUN = function(x) { 
                                                    if (is.list(x)) 
                                                      matrix(rep(NA, number_dimensions + number_dimensions^2), nrow = 1, dimnames = list(NULL, names(estimates_and_variance[[1]])))
                                                    else 
                                                      x } ) 
  
  conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection, iterations_per_unique_condition, number_dimensions)
  
  condition_vector <- sort(rep(1:(ncol(estimates_and_variance_without_errors)/iterations_per_unique_condition), iterations_per_unique_condition))
  
  estimates_and_conditions <- cbind(t(estimates_and_variance_without_errors)[, c("estimated_theta1", "estimated_theta2", "estimated_theta3", "variance_estimate1", "variance_estimate5", "variance_estimate9")], 
                                    conditions[, c("iteration", "number_items", "number_answer_categories", "model", "estimator", "information_summary")], 
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

  test_that("test prior_var_safe_nlm is 100", {
    # run again with new safe_ml code in estimate_latent_trait()
    # run ML only
    iterations_per_unique_condition <- 100
    true_theta_vec <- c(-2, 1, 2)
    number_items_vec <- c(300)
    number_answer_categories_vec <- c(2, 4)
    number_dimensions <- 3
    
    start_items <- list(type = 'random_by_dimension', n_by_dimension = 3, n = 9)
    variance_target <- .1^2
    model_vec <- c("3PLM","GRM","GPCM","SM")
    estimator_vec <- c("ML", "MAP")
    information_summary_vec <- c("D", "PD", "A", "PA")
    item_selection <- 'MI'
    
    estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection, start_items, variance_target, iterations_per_unique_condition, number_dimensions, prior_var_safe_nlm = 100)
    
    conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection_vec, iterations_per_unique_condition, number_dimensions)
    
    condition_vector <- sort(rep(1:(ncol(estimates_and_variance)/iterations_per_unique_condition), iterations_per_unique_condition))
    
    estimates_and_conditions <- cbind(t(estimates_and_variance)[, c("estimated_theta1", "estimated_theta2", "estimated_theta3", "variance_estimate1", "variance_estimate5", "variance_estimate9")], 
                                      conditions[, c("iteration", "number_items", "number_answer_categories", "model", "estimator", "information_summary")], 
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
    
    # no errors/missings for MAP estimator
    expect_equal(number_na_per_condition[which(conditions[seq(1, 6400, 100), "estimator"] == "MAP"), "x"], rep(0, 32))
    # no errors/missings for ML estimator
    expect_equal(number_na_per_condition[which(conditions[seq(1, 6400, 100), "estimator"] == "ML"), "x"], rep(0, 32))
  })
  
  test_that("items load on all dimensions prior_var_safe_nlm is 100", {
    # run again with new safe_ml code in estimate_latent_trait()
    iterations_per_unique_condition <- 100 
    true_theta_vec <- c(2, -1, -2)
    number_items_vec <- c(300)
    number_answer_categories_vec <- c(2, 4)
    number_dimensions <- 3
    
    start_items <- list(type = 'random_by_dimension', n_by_dimension = 3, n = 9)
    variance_target <- .1^2
    model_vec <- c("3PLM","GRM","GPCM","SM")
    estimator_vec <- c("ML", "MAP") # AEP
    information_summary_vec <- c("D", "PD", "A", "PA") # PEKL 
    item_selection <- 'MI'
    
    estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection, start_items, variance_target, iterations_per_unique_condition, number_dimensions, items_load_one_dimension = FALSE, prior_var_safe_nlm = 100)
    #save(estimates_and_variance, file = "/Users/rivkadevries/Desktop/simulationsCAT/estimates_and_variance_within.R")
    estimates_and_variance_without_errors <- sapply(estimates_and_variance, 
                                                    FUN = function(x) { 
                                                      if (is.list(x)) 
                                                        matrix(rep(NA, number_dimensions + number_dimensions^2), nrow = 1, dimnames = list(NULL, names(estimates_and_variance[[1]])))
                                                      else 
                                                        x } ) 
    
    conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection, iterations_per_unique_condition, number_dimensions)
    
    condition_vector <- sort(rep(1:(ncol(estimates_and_variance_without_errors)/iterations_per_unique_condition), iterations_per_unique_condition))
    
    estimates_and_conditions <- cbind(t(estimates_and_variance_without_errors)[, c("estimated_theta1", "estimated_theta2", "estimated_theta3", "variance_estimate1", "variance_estimate5", "variance_estimate9")], 
                                      conditions[, c("iteration", "number_items", "number_answer_categories", "model", "estimator", "information_summary")], 
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
    
    # no errors/missings for MAP estimator
    expect_equal(number_na_per_condition[which(conditions[seq(1, 6400, 100), "estimator"] == "MAP"), "x"], rep(0, 32))
    # some errors/missings for ML estimator
    expect_equal(number_na_per_condition[which(conditions[seq(1, 6400, 100), "estimator"] == "ML"), "x"], 
                 c(0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 1, 0))
    
  })
  
  test_that("items load on all dimensions prior_var_safe_nlm is 1", {
    # run again with new safe_ml code in estimate_latent_trait()
    iterations_per_unique_condition <- 100 
    true_theta_vec <- c(2, -1, -2)
    number_items_vec <- c(300)
    number_answer_categories_vec <- c(2, 4)
    number_dimensions <- 3
    
    start_items <- list(type = 'random_by_dimension', n_by_dimension = 3, n = 9)
    variance_target <- .1^2
    model_vec <- c("3PLM","GRM","GPCM","SM")
    estimator_vec <- c("ML", "MAP") # AEP
    information_summary_vec <- c("D", "PD", "A", "PA") # PEKL 
    item_selection <- 'MI'
    
    estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection, start_items, variance_target, iterations_per_unique_condition, number_dimensions, items_load_one_dimension = FALSE, prior_var_safe_nlm = 1)
    #save(estimates_and_variance, file = "/Users/rivkadevries/Desktop/simulationsCAT/estimates_and_variance_within_safe_var1.R")
    estimates_and_variance_without_errors <- sapply(estimates_and_variance, 
                                                    FUN = function(x) { 
                                                      if (is.list(x)) 
                                                        matrix(rep(NA, number_dimensions + number_dimensions^2), nrow = 1, dimnames = list(NULL, names(estimates_and_variance[[1]])))
                                                      else 
                                                        x } ) 
    
    conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection, iterations_per_unique_condition, number_dimensions)
    
    condition_vector <- sort(rep(1:(ncol(estimates_and_variance_without_errors)/iterations_per_unique_condition), iterations_per_unique_condition))
    
    estimates_and_conditions <- cbind(t(estimates_and_variance_without_errors)[, c("estimated_theta1", "estimated_theta2", "estimated_theta3", "variance_estimate1", "variance_estimate5", "variance_estimate9")], 
                                      conditions[, c("iteration", "number_items", "number_answer_categories", "model", "estimator", "information_summary")], 
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
    
    # no errors/missings for MAP estimator
    expect_equal(number_na_per_condition[which(conditions[seq(1, 6400, 100), "estimator"] == "MAP"), "x"], rep(0, 32))
    # some errors/missings for ML estimator
    expect_equal(number_na_per_condition[which(conditions[seq(1, 6400, 100), "estimator"] == "ML"), "x"], 
                 c(0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0))
    
  })
  
  
  test_that("simulate with constraints, max_n 130", {
    # run again with new safe_ml code in estimate_latent_trait()
    iterations_per_unique_condition <- 100
    true_theta_vec <- c(-2, 1, 2)
    number_items_vec <- 300 # can only have length one here (with item characteristics) and should be divisible by 3, to keep things simple
    number_answer_categories_vec <- 4
    number_dimensions <- 3
    
    start_items <- list(type = 'random_by_dimension', n_by_dimension = 3, n = 9)
    variance_target <- .001^2
    model_vec <- "SM"
    estimator_vec <- c("ML", "MAP") # AEP
    information_summary_vec <- "D" # PEKL 
    item_selection <- 'Shadow' 
    max_n = 130
    
    #create item characteristics and constraints
    characteristics <- data.frame(content = c(rep('depression', number_items_vec / 3), rep('anxiety', number_items_vec / 3), rep('somatic', number_items_vec / 3)))
    constraints <- list(list(name = 'content/depression',
                             op = '><',
                             target = c(50, 75)),
                        list(name = 'content/somatic',
                             op = '><',
                             target = c(75, 100)))
    
    estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection, start_items, variance_target, iterations_per_unique_condition, number_dimensions, constraints_and_characts = list(characteristics = characteristics, constraints = constraints), prior_var_safe_nlm = 100, return_administered_item_indeces = TRUE, max_n = max_n)
    #save(estimates_and_variance, file = "/Users/rivkadevries/Desktop/simulationsCAT/estimates_and_variance_3dim_constraints_maxn130.R")
    
    conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection, iterations_per_unique_condition, number_dimensions)
    
    condition_vector <- sort(rep(1:(ncol(estimates_and_variance)/iterations_per_unique_condition), iterations_per_unique_condition))
    
    estimates_and_conditions <- cbind(t(estimates_and_variance)[, c("estimated_theta1", "estimated_theta2", "estimated_theta3", "variance_estimate1", "variance_estimate5", "variance_estimate9",
                                                                    str_c("items_administered", 1:max_n))], 
                                      conditions[, c("iteration", "number_items", "number_answer_categories", "model", "estimator", "information_summary")], 
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
    
    # no errors/missings for MAP estimator
    expect_equal(number_na_per_condition[which(conditions[seq(1, 6400, 100), "estimator"] == "MAP"), "x"], 0)
    # no errors/missings for ML estimator
    expect_equal(number_na_per_condition[which(conditions[seq(1, 6400, 100), "estimator"] == "ML"), "x"], 0)
    
    expect_equal(number_depression_items, rep(50, 200))
    expect_equal(number_anxiety_items, rep(5, 200))
    expect_equal(number_somatic_items, rep(75, 200)) 
  })
  
  test_that("simulate with constraints, max_n 260", {
    # run again with new safe_ml code in estimate_latent_trait
    iterations_per_unique_condition <- 100
    true_theta_vec <- c(-2, 1, 2)
    number_items_vec <- 300 # can only have length one here (with item characteristics) and should be divisible by 3, to keep things simple
    number_answer_categories_vec <- 4
    number_dimensions <- 3
    
    start_items <- list(type = 'random_by_dimension', n_by_dimension = 3, n = 9)
    variance_target <- .001^2
    model_vec <- "SM"
    estimator_vec <- c("ML", "MAP") # AEP
    information_summary_vec <- "D" # PEKL 
    item_selection <- 'Shadow' 
    max_n = 260
    
    #create item characteristics and constraints
    characteristics <- data.frame(content = c(rep('depression', number_items_vec / 3), rep('anxiety', number_items_vec / 3), rep('somatic', number_items_vec / 3)))
    constraints <- list(list(name = 'content/depression',
                             op = '><',
                             target = c(50, 75)),
                        list(name = 'content/somatic',
                             op = '><',
                             target = c(75, 90)))
    
    estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection, start_items, variance_target, iterations_per_unique_condition, number_dimensions, constraints_and_characts = list(characteristics = characteristics, constraints = constraints), prior_var_safe_nlm = 100, return_administered_item_indeces = TRUE, max_n = max_n)
    #save(estimates_and_variance, file = "/Users/rivkadevries/Desktop/simulationsCAT/estimates_and_variance_3dim_constraints_maxn260.R")
    
    conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection, iterations_per_unique_condition, number_dimensions)
    
    condition_vector <- sort(rep(1:(ncol(estimates_and_variance)/iterations_per_unique_condition), iterations_per_unique_condition))
    
    estimates_and_conditions <- cbind(t(estimates_and_variance)[, c("estimated_theta1", "estimated_theta2", "estimated_theta3", "variance_estimate1", "variance_estimate5", "variance_estimate9",
                                                                    str_c("items_administered", 1:max_n))], 
                                      conditions[, c("iteration", "number_items", "number_answer_categories", "model", "estimator", "information_summary")], 
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
  
  test_that("MAP with informative prior", {
    # run again with PEKL
    iterations_per_unique_condition <- 100
    true_theta_vec <- c(-2, 1, 2)
    number_items_vec <- c(300)
    number_answer_categories_vec <- c(2, 4)
    number_dimensions <- 3
    
    start_items <- list(type = 'random_by_dimension', n_by_dimension = 3, n = 9)
    variance_target <- .1^2
    model_vec <- c("3PLM","GRM","GPCM","SM")
    estimator_vec <- "MAP"
    information_summary_vec <- c("D", "PD", "A", "PA") # PEKL 
    item_selection <- 'MI'
    
    estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection, start_items, variance_target, iterations_per_unique_condition, number_dimensions, prior = diag(3) * .5, prior_var_safe_nlm = 100)
    
    conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection_vec, iterations_per_unique_condition, number_dimensions)
    
    condition_vector <- sort(rep(1:(ncol(estimates_and_variance)/iterations_per_unique_condition), iterations_per_unique_condition))
    
    estimates_and_conditions <- cbind(t(estimates_and_variance)[, c("estimated_theta1", "estimated_theta2", "estimated_theta3", "variance_estimate1", "variance_estimate5", "variance_estimate9")], 
                                      conditions[, c("iteration", "number_items", "number_answer_categories", "model", "estimator", "information_summary")], 
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


