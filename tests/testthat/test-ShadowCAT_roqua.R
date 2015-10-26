# only for whithin R:
'
library(testthat)
library(pbapply)
'

# Tests take a while, select shorter ones later
if (FALSE) {
  
make_random_seed_exist <- rnorm(1)

test_shadowcat_roqua <- function(true_theta, prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, item_selection = "MI", constraints = NULL, lowerbound = rep(-3, ncol(alpha)), upperbound = rep(3, ncol(alpha))) {
  new_response <- NULL
  next_item_and_test_outcome <- shadowcat_roqua(new_response, prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, item_selection, constraints, lowerbound, upperbound)
  
  while (next_item_and_test_outcome$index_new_item != "stop_test") {
    person_updated_after_new_response$theta <- true_theta
    new_response <- tail(answer(person = person_updated_after_new_response, 
                                test = initTest(items = initItembank(model, alpha, beta, guessing, eta, silent = TRUE), 
                                                start = start_items, 
                                                stop = stop_test,
                                                max_n = stop_test$n,
                                                estimator = estimator,
                                                objective = information_summary,
                                                selection = item_selection), 
                                indeces = next_item_and_test_outcome$index_new_item)$responses, 1) # simulate person response on new item
    
    next_item_and_test_outcome <- shadowcat_roqua(new_response, prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, item_selection, constraints, lowerbound, upperbound)   
  }
  
  next_item_and_test_outcome$person_updated_after_new_response
}

get_conditions <- function(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection_vec, iterations_per_unique_condition, number_dimensions) {
  if (number_dimensions == 1) {
    conditions <- expand.grid(1:iterations_per_unique_condition, true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection_vec)
    colnames(conditions) <- c("iteration", "true_theta", "number_items", "number_answer_categories", "model", "estimator", "information_summary", "item_selection")
  }
  else {
    conditions <- expand.grid(1:iterations_per_unique_condition, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection_vec)
    colnames(conditions) <- c("iteration", "number_items", "number_answer_categories", "model", "estimator", "information_summary", "item_selection") 
  }
  conditions
}

run_simulation <- function(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection_vec, prior, start_items, variance_target, iterations_per_unique_condition, number_dimensions, item_selection = "MI", constraints = NULL, guessing = NULL, items_load_one_dimension = TRUE, lowerbound = rep(-3, number_dimensions), upperbound = rep(3, number_dimensions)) {                   
  conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection_vec, iterations_per_unique_condition, number_dimensions)
  
  pbapply::pbsapply(1:nrow(conditions), 
                    FUN = function(condition) {
                    prior <- diag(number_dimensions) * 20
                    stop_test <- list(type = 'variance', target = variance_target, n = conditions[condition, "number_items"])
                    true_theta <- ( if (number_dimensions == 1) 
                                      conditions[condition, "true_theta"]
                                    else
                                      true_theta_vec )
                    alpha_beta <- createTestBank(model = as.character(conditions[condition, "model"]), K = conditions[condition, "number_items"], Q = number_dimensions, M = conditions[condition, "number_answer_categories"] - 1, between = items_load_one_dimension, run_initItembank = FALSE)
                    estimate_theta <- tryCatch(test_shadowcat_roqua(true_theta, prior, as.character(conditions[condition, "model"]), alpha_beta$alpha, alpha_beta$beta, guessing, eta = NULL, start_items, stop_test, as.character(conditions[condition, "estimator"]), as.character(conditions[condition, "information_summary"]), item_selection, constraints, lowerbound, upperbound)$estimate,
                                               error = function(e) e)
           
                    c("estimated_theta" = estimate_theta,
                      "variance_estimate" = attr(estimate_theta, "variance"))
                    })
}

context("simulations")

test_that("one dimension, no constraints on item selection, one iteration per condition", {
  # EAP estimation does not work
  # PEKL information summaries give errors because it also makes use of EAP estimation
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
  estimator_vec <- c("ML", "MAP") # AEP
  information_summary_vec <- c("D", "PD", "A", "PA") # PEKL 
  item_selection_vec <- 'MI'
  
  estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection_vec, prior, start_items, variance_target, iterations_per_unique_condition, number_dimensions, lowerbound = lowerbound, upperbound = upperbound)
  
  conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection_vec, iterations_per_unique_condition, number_dimensions)
  
  estimates_and_conditions <- cbind(t(estimates_and_variance), conditions[, c("true_theta", "number_items", "number_answer_categories", "model", "estimator", "information_summary")])
  
  average_per_true_theta_and_number_items <- aggregate(estimates_and_conditions[,"estimated_theta"], list(estimates_and_conditions[,"true_theta"], estimates_and_conditions[,"number_items"]), "mean")
  sd_per_true_theta_and_number_items <- aggregate(estimates_and_conditions[,"estimated_theta"], list(estimates_and_conditions[,"true_theta"], estimates_and_conditions[,"number_items"]), "sd") 
  
  expect_equal(round(average_per_true_theta_and_number_items[,"x"], 3), c(-1.943, 1.157, -2.084, 1.053))
  expect_equal(round(sd_per_true_theta_and_number_items[,"x"], 3), c(.661, .593, .426, .376))
  expect_equal(fivenum(round(estimates_and_conditions[,"variance_estimate"], 3)), c(4.4e-02, 1.2e-01, 2.06e-01, 3.76e-01, 138953441))
  # This extremely large variance corresponds to an estimated theta at the lowerbound of -3; setting the lowerbound to -5 does not help
  
})

test_that("one dimension, no constraints on item selection, two iterations per condition", {
  # EAP estimation does not work
  # PEKL information summaries give errors because it also makes use of EAP estimation
  iterations_per_unique_condition <- 100 # extend to 100 when time permits and when EAP issue is fixed
  true_theta_vec <- c(-2, 1)
  number_items_vec <- c(50, 100)
  number_answer_categories_vec <- c(2, 4)
  number_dimensions <- 1
  
  start_items <- list(type = 'random', n = 3)
  variance_target <- .1^2
  model_vec <- c("3PLM","GRM","GPCM","SM")
  estimator_vec <- c("ML", "MAP") # AEP
  information_summary_vec <- c("D", "PD", "A", "PA") # PEKL 
  item_selection_vec <- 'MI'
  
  estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection_vec, prior, start_items, variance_target, iterations_per_unique_condition, number_dimensions)
  
  conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection_vec, iterations_per_unique_condition, number_dimensions)
  condition_vector <- sort(rep(1:(ncol(estimates_and_variance)/iterations_per_unique_condition), iterations_per_unique_condition))
  
  estimates_and_conditions <- cbind(t(estimates_and_variance), conditions[, c("true_theta", "iteration", "number_items", "number_answer_categories", "model", "estimator", "information_summary")], condition_vector)
  
  average_per_condition_true_minus2 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == -2), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == -2)]), "mean")
  sd_per_condition_true_minus2 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == -2), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == -2)]), "sd")
  
  average_per_condition_true_1 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 1), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 1)]), "mean")
  sd_per_condition_true_1 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 1), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 1)]), "sd") 
  
  expect_equal(round(average_per_condition_true_minus2[,"x"], 3)[1:4], c(-2.268, -1.733, -1.171, -2.122))
  expect_equal(round(average_per_condition_true_1[,"x"], 3)[1:4], c(.768, .779, 1.435, 1.016))
  expect_equal(round(sd_per_condition_true_minus2[,"x"], 3)[1:4], c(1.035, .407, .100, .102))
  expect_equal(round(sd_per_condition_true_1[,"x"], 3)[1:4], c(.278, .178, .403, .081))
  
  expect_equal(round(fivenum(average_per_condition_true_minus2[,"x"]), 3), c(-2.717, -2.174, -1.999, -1.863, -1.171))
  expect_equal(round(fivenum(average_per_condition_true_1[,"x"]), 3), c(.442, .877, 1.005, 1.147, 1.735))
  expect_equal(round(range(sd_per_condition_true_minus2[,"x"]), 3), c(.009, 1.035))
  expect_equal(round(range(sd_per_condition_true_1[,"x"]), 3), c(.001, .827))
  
  expect_equal(round(sqrt(fivenum(estimates_and_conditions[which(estimates_and_conditions[,"number_items"] == 50), "variance_estimate"])), 3), c(.213, .282, .355, .419, .901))
  expect_equal(round(sqrt(fivenum(estimates_and_conditions[which(estimates_and_conditions[,"number_items"] == 100), "variance_estimate"])), 3), c(.152, .198, .247, .292, .499))
 
})

test_that("three dimensions, no constraints on item selection, ten iterations per condition", {
  # EAP estimation does not work
  # PEKL information summaries give errors because it also makes use of EAP estimation
  # ML estimate gives error here
  iterations_per_unique_condition <- 100 # extend to 100 when time permits and when EAP issue is fixed
  true_theta_vec <- c(-2, 1, 2)
  number_items_vec <- c(300)
  number_answer_categories_vec <- c(2, 4)
  number_dimensions <- 3
  
  start_items <- list(type = 'random', n = 3)
  variance_target <- .1^2
  model_vec <- c("3PLM","GRM","GPCM","SM")
  estimator_vec <- c("ML", "MAP") # AEP
  information_summary_vec <- c("D", "PD", "A", "PA") # PEKL 
  item_selection_vec <- 'MI'
  
  estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection_vec, prior, start_items, variance_target, iterations_per_unique_condition, number_dimensions)
  
  conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection_vec, iterations_per_unique_condition, number_dimensions)
  condition_vector <- sort(rep(1:(ncol(estimates_and_variance)/iterations_per_unique_condition), iterations_per_unique_condition))
  
  estimates_and_conditions <- cbind(t(estimates_and_variance)[, c("estimated_theta1", "estimated_theta2", "estimated_theta3", "variance_estimate1", "variance_estimate5", "variance_estimate9")], 
                                    conditions[, c("iteration", "number_items", "number_answer_categories", "model", "estimator", "information_summary")], 
                                    condition_vector)
  
  average_per_condition_dim1 <- aggregate(estimates_and_conditions[, "estimated_theta1"], list(condition_vector), "mean")
  average_per_condition_dim2 <- aggregate(estimates_and_conditions[, "estimated_theta2"], list(condition_vector), "mean")
  average_per_condition_dim3 <- aggregate(estimates_and_conditions[, "estimated_theta3"], list(condition_vector), "mean")
  
  sd_per_condition_dim1 <- aggregate(estimates_and_conditions[, "estimated_theta1"], list(condition_vector), "sd")
  sd_per_condition_dim2 <- aggregate(estimates_and_conditions[, "estimated_theta2"], list(condition_vector), "sd")
  sd_per_condition_dim3 <- aggregate(estimates_and_conditions[, "estimated_theta3"], list(condition_vector), "sd")
  
  expect_equal(round(fivenum(average_per_condition_dim1[,"x"]), 3), c())
  expect_equal(round(fivenum(average_per_condition_dim2[,"x"]), 3), c())
  expect_equal(round(fivenum(average_per_condition_dim3[,"x"]), 3), c())
  
  expect_equal(round(range(sd_per_condition_dim1[,"x"]), 3), c())
  expect_equal(round(range(sd_per_condition_dim2[,"x"]), 3), c())
  expect_equal(round(range(sd_per_condition_dim3[,"x"]), 3), c())
  
  expect_equal(round(sqrt(fivenum(estimates_and_conditions[, "variance_estimate1"])), 3), c())
  expect_equal(round(sqrt(fivenum(estimates_and_conditions[, "variance_estimate5"])), 3), c())
  expect_equal(round(sqrt(fivenum(estimates_and_conditions[, "variance_estimate9"])), 3), c())
  
})

context("validate shadowcat_roqua single conditions")

test_that("true theta is 2", {
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
  stop_test <- list(type = 'length', n = 100)
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

test_that("true theta is 1, 0, 2", {  
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
  stop_test <- list(type = 'length', n = 300)
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
  stop_test <- list(type = 'length', n = 300)
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


test_that("validate with ShadowCAT Karel Kroeze", {
  # define true theta for simulation of responses
  true_theta <- 2
  
  # define item characteristics
  number_items <- 50
  number_dimensions <- 1
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  eta <- NULL # only relevant for GPCM model
  
  model <- '3PLM'
  start_items <- list(type = 'random', n = 1)
  stop_test <- list(type = 'length', n = 50)
  estimator <- 'MAP'
  information_summary <- 'PD'
  item_selection <- 'MI'
  
  # define prior covariance matrix
  prior <- diag(number_dimensions) 

  # make person and test object for ShadowCAT from Kroeze
  item_characteristics_shadowcat_format <- initItembank(model = model, alpha = alpha, beta = beta, guessing = guessing, eta = eta, silent = TRUE)
  
  person <- initPerson(items = item_characteristics_shadowcat_format,
                       theta = true_theta,
                       prior = prior)
  
  test <- initTest(items = item_characteristics_shadowcat_format, 
                   start = start_items, 
                   stop = stop_test,
                   max_n = stop_test$n,
                   estimator = estimator,
                   objective = information_summary,
                   selection = item_selection,
                   constraints = NULL,
                   exposure = NULL,
                   lowerBound = rep(-3, item_characteristics_shadowcat_format$Q), 
                   upperBound = rep(3, item_characteristics_shadowcat_format$Q))
  
  test_outcome_roqua <- with_random_seed(2, test_shadowcat_roqua)(true_theta, prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary)
  person_kroeze <- with_random_seed(2, ShadowCAT)(person, test)
  
  expect_equal(person_kroeze$estimate, test_outcome_roqua$estimate)
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
  stop_test <- list(type = 'length', n = 10)
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
  stop_test <- list(type = 'variance', target = .5, n = 50)
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

# gives error at this point due to bug in updated MultiGHQuad package
if (FALSE) {  
  context("problematic")

  # This one gets stuck after the 248st administered item. At that point, the estimated variance matrix becomes
  # non-invertible (has an eigenvalue of -5.530436e-18). See the problems in the trans function within init.quad().
  # Also, when running the 248 items, the resulting estimates are poor. 
  # After update of MultiGHQuad the EAP estimate gives an error
  test_that("true theta is 2, -2, 3", {
    # define true theta for simulation of responses
    true_theta <- c(2, -2, 3)
    
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
    stop_test <- list(type = 'length', n = 248)
    estimator <- 'EAP'
    information_summary <- 'PD'
    item_selection <- 'MI'
    
    # define prior covariance matrix
    prior <- diag(number_dimensions) * 5
    
    test_outcome <- with_random_seed(2, test_shadowcat_roqua)(true_theta, prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary)
     
    expect_equal(as.vector(round(test_outcome$estimate, 3)), c(.984, .171, -.055))
    expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3))[1:3],c(.238, .259, .122))
    expect_equal(length(test_outcome$administered), 248)
  })
}

}
