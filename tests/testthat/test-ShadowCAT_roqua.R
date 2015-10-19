# only for whithin R:
'
library(testthat)
'
make_random_seed_exist <- rnorm(1)

test_shadowcat_roqua <- function(true_theta, prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, item_selection = "MI", constraints = NULL) {
  new_response <- NULL
  next_item_and_test_outcome <- shadowcat_roqua(new_response, prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, item_selection, constraints)
  
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
    
    next_item_and_test_outcome <- shadowcat_roqua(new_response, prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, item_selection, constraints)   
  }
  
  next_item_and_test_outcome$person_updated_after_new_response
}

run_simulation <- function(true_theta_vec, number_items_vec, number_dimensions_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection_vec, prior, start_items, variance_target, item_selection = "MI", constraints = NULL, guessing = NULL) {
  conditions <- expand.grid(true_theta_vec, number_items_vec, number_dimensions_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection_vec)
  colnames(conditions) <- c("true_theta", "number_items", "number_dimensions", "number_answer_categories", "model", "estimator", "information_summary", "item_selection")
    
  sapply(1:nrow(conditions), 
         FUN = function(condition) {
           stop_test <- list(type = 'variance', target = variance_target, n = conditions[condition, "number_items"])
           true_theta <- conditions[condition, "true_theta"]
           alpha_beta <- createTestBank(model = as.character(conditions[condition, "model"]), K = conditions[condition, "number_items"], Q = conditions[condition, "number_dimensions"], M = conditions[condition, "number_answer_categories"] - 1, run_initItembank = FALSE)
           estimate_theta <- tryCatch(test_shadowcat_roqua(true_theta, prior, as.character(conditions[condition, "model"]), alpha_beta$alpha, alpha_beta$beta, guessing, eta = NULL, start_items, stop_test, as.character(conditions[condition, "estimator"]), as.character(conditions[condition, "information_summary"]), item_selection, constraints)$estimate,
                                      error = function(e) e)
           
           c("estimated_theta" = estimate_theta,
             "variance_estimate" = attr(estimate_theta, "variance"))
         })
}

context("validate shadowcat_roqua")

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
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), 1.912)
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3)), .112)
  expect_equal(test_outcome$available, numeric(0))
  expect_equal(length(test_outcome$administered), 100)
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
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), -0.022)
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3)), .315)
  expect_equal(test_outcome$available, c(1:4, 8:9, 11:16, 18:23, 25:29, 31:40, 42:46, 49:50))
  expect_equal(test_outcome$administered, c(10, 30, 48, 7, 24, 5, 17, 6, 47, 41))
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
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), .608)
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3)), 0.468)
  expect_equal(test_outcome$available, c(1:4, 6, 8:9, 11:16, 18:23, 25:29, 31:47, 49:50))
  expect_equal(test_outcome$administered, c(10, 30, 48, 7, 24, 5, 17))
  expect_equal(test_outcome$responses, c(1, 0, 1, 1, 0, 1, 0))
})

# not run because it takes too much time
if (FALSE) {
  
  context("simulations")
  
  test_that("one dimension, no constraints on item selection", {
    # not all information summaries work well, EAP estimation seems to give either 3 or -3 as theta estimate
    true_theta_vec <- c(-2, 1)
    number_items_vec <- c(15, 50)
    number_answer_categories_vec <- c(2, 4)
    number_dimensions_vec <- 1
    
    start_items <- list(type = 'random', n = 3)
    variance_target <- .1
    model_vec <- c("3PLM","GRM","GPCM","SM")
    estimator_vec <- c("ML", "MAP", "EAP")
    information_summary_vec <- c("D", "PD", "A", "AD")
    item_selection_vec <- 'MI'
    
    prior <- diag(number_dimensions) * 5
    
    estimates_and_variance <- run_simulation(true_theta_vec, number_items_vec, number_dimensions_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection_vec, prior, start_items, variance_target)
    
    conditions <- expand.grid(true_theta_vec, number_items_vec, number_dimensions_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection_vec)
    colnames(conditions) <- c("true_theta", "number_items", "number_dimensions", "number_answer_categories", "model", "estimator", "information_summary", "item_selection")
    
    estimates_and_conditions <- cbind(t(estimates_and_variance), conditions[, c("true_theta", "number_items", "number_dimensions", "number_answer_categories", "model", "estimator", "information_summary")])
    
  })
  
  
  context("problematic")

  # This one gets stuck after the 248st administered item. At that point, the estimated variance matrix becomes
  # non-invertible (has an eigenvalue of -5.530436e-18). See the problems in the trans function within init.quad().
  # Also, when running the 248 items, the resulting estimates are poor. 
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


