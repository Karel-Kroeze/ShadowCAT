# only for whithin R:
'
library(testthat)
library(pbapply)
library(stringr) # only for test constraints simulations
'

make_random_seed_exist <- rnorm(1)

test_shadowcat_roqua <- function(true_theta, prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, item_selection = "MI", constraints = NULL, lowerbound = rep(-3, ncol(alpha)), upperbound = rep(3, ncol(alpha)), prior_var_safe_ml = NULL, initital_estimate = rep(0, ncol(alpha)), initial_variance = prior) {
  new_response <- NULL
  attr(initital_estimate, 'variance') <- initial_variance
  next_item_and_test_outcome <- shadowcat_roqua(new_response, estimate = initital_estimate , responses = numeric(0), administered = numeric(0), available = 1:nrow(beta), model, alpha, beta, start_items, stop_test, estimator, information_summary, prior, guessing, eta, item_selection, constraints, lowerbound, upperbound, prior_var_safe_ml)

  while (next_item_and_test_outcome$index_new_item != "stop_test") {
    new_response <- answer(true_theta, model, ncol(alpha), estimator, alpha, beta, guessing, ncol(beta), next_item_and_test_outcome$index_new_item)
    next_item_and_test_outcome <- shadowcat_roqua(new_response, next_item_and_test_outcome$estimate, next_item_and_test_outcome$responses, next_item_and_test_outcome$administered, next_item_and_test_outcome$available, model, alpha, beta, start_items, stop_test, estimator, information_summary, prior,  guessing, eta, item_selection, constraints, lowerbound, upperbound, prior_var_safe_ml)  
  }
  
  next_item_and_test_outcome
}

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

# if number dimensions is 1, simulations are performed for each true theta value in true_theta_vec
# if number dimensions is larger than 1, true_theta_vec is interpreted as containing the true thetas for each dimension
# item_selection can be "MI" or "Shadow". In case of "Shadow", constraints should be defined, and number_items_vec can only have length 1
# max_n can only have length 1; if null, max_n is set to the number of items (which may differ accross conditions)
run_simulation <- function(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection, start_items, variance_target, iterations_per_unique_condition, number_dimensions, constraints = NULL, guessing = NULL, items_load_one_dimension = TRUE, lowerbound = rep(-3, number_dimensions), upperbound = rep(3, number_dimensions), prior = diag(number_dimensions) * 20, prior_var_safe_ml = NULL, return_administered_item_indeces = FALSE, max_n = NULL) {                   
  conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection, iterations_per_unique_condition, number_dimensions)
  
  pbapply::pbsapply(1:nrow(conditions), 
                    FUN = function(condition) {
                    prior <- prior
                    if (is.null(max_n))
                      max_n <- conditions[condition, "number_items"] 
                    stop_test <- list(target = variance_target, n = max_n)
                    true_theta <- ( if (number_dimensions == 1) 
                                      conditions[condition, "true_theta"]
                                    else
                                      true_theta_vec )
                    alpha_beta <- createTestBank(model = as.character(conditions[condition, "model"]), K = conditions[condition, "number_items"], Q = number_dimensions, M = conditions[condition, "number_answer_categories"] - 1, between = items_load_one_dimension, run_initItembank = FALSE)
                    estimate_theta <- tryCatch(test_shadowcat_roqua(true_theta, prior, as.character(conditions[condition, "model"]), alpha_beta$alpha, alpha_beta$beta, guessing, eta = NULL, start_items, stop_test, as.character(conditions[condition, "estimator"]), as.character(conditions[condition, "information_summary"]), item_selection, constraints, lowerbound, upperbound, prior_var_safe_ml),
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
  stop_test <- list(n = 100)
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
  stop_test <- list(n = 100)
  estimator <- 'ML'
  information_summary <- 'D'
  item_selection <- 'MI'
  
  test_outcome <- with_random_seed(2, test_shadowcat_roqua)(true_theta, prior = NULL, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, initial_variance = diag(1))
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), 2.169)
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3)), .129)
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
  stop_test <- list(n = 300)
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
  stop_test <- list(n = 300)
  estimator <- 'ML'
  information_summary <- 'D'
  item_selection <- 'MI'
  
  test_outcome <- with_random_seed(3, test_shadowcat_roqua)(true_theta, prior = NULL, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, initial_variance = diag(1))
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), c(.755, -.070, 2.221))
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3))[1:3],c(.063, .000, .000))
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
  stop_test <- list(n = 300)
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
  stop_test <- list(n = 300)
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
  start_items <- list(type = 'randomByDimension', nByDimension = 3, n = 9)
  stop_test <- list(n = 260)
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
  
  test_outcome <- with_random_seed(3, test_shadowcat_roqua)(true_theta, prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, item_selection =  item_selection, constraints = list(characteristics = characteristics, constraints = constraints))
  
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
  start_items <- list(type = 'randomByDimension', nByDimension = 3, n = 9)
  stop_test <- list(n = 130)
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
  
  test_outcome <- with_random_seed(3, test_shadowcat_roqua)(true_theta, prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary, item_selection =  item_selection, constraints = list(characteristics = characteristics, constraints = constraints))
  
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
  stop_test <- list(n = 10)
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
  stop_test <- list(target = .5, n = 50)
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
    stop_test <- list(n = 248)
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

# these simulations take a long time to run, if (FALSE) ensures that they are not each time the tests are run
if (FALSE) {
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
    item_selection <- 'MI'
    
    estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection, start_items, variance_target, iterations_per_unique_condition, number_dimensions, lowerbound = lowerbound, upperbound = upperbound)
    
    conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection_vec, iterations_per_unique_condition, number_dimensions)
    
    estimates_and_conditions <- cbind(t(estimates_and_variance), conditions[, c("true_theta", "number_items", "number_answer_categories", "model", "estimator", "information_summary")])
    
    average_per_true_theta_and_number_items <- aggregate(estimates_and_conditions[,"estimated_theta"], list(estimates_and_conditions[,"true_theta"], estimates_and_conditions[,"number_items"]), "mean")
    sd_per_true_theta_and_number_items <- aggregate(estimates_and_conditions[,"estimated_theta"], list(estimates_and_conditions[,"true_theta"], estimates_and_conditions[,"number_items"]), "sd") 
    
    expect_equal(round(average_per_true_theta_and_number_items[,"x"], 3), c(-1.943, 1.157, -2.084, 1.053))
    expect_equal(round(sd_per_true_theta_and_number_items[,"x"], 3), c(.661, .593, .426, .376))
    expect_equal(fivenum(round(estimates_and_conditions[,"variance_estimate"], 3)), c(4.4e-02, 1.2e-01, 2.06e-01, 3.76e-01, 138953441))
    # This extremely large variance corresponds to an estimated theta at the lowerbound of -3; setting the lowerbound to -5 does not help
    
  })
  
  test_that("one dimension, estimates ML and MAP, information summary D, PD, A, and PA, no constraints on item selection, 100 iterations per condition", {
    # EAP estimation does not work
    # PEKL information summaries give errors because it also makes use of EAP estimation
    iterations_per_unique_condition <- 100
    true_theta_vec <- c(-2, 1)
    number_items_vec <- c(50, 100)
    number_answer_categories_vec <- c(2, 4)
    number_dimensions <- 1
    
    start_items <- list(type = 'random', n = 3)
    variance_target <- .1^2
    model_vec <- c("3PLM","GRM","GPCM","SM")
    estimator_vec <- c("ML", "MAP") # AEP
    information_summary_vec <- c("D", "PD", "A", "PA") # PEKL 
    item_selection <- 'MI'
    
    estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection, start_items, variance_target, iterations_per_unique_condition, number_dimensions)
    
    conditions <- get_conditions(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection_vec, iterations_per_unique_condition, number_dimensions)
    condition_vector <- sort(rep(1:(ncol(estimates_and_variance)/iterations_per_unique_condition), iterations_per_unique_condition))
    
    estimates_and_conditions <- cbind(t(estimates_and_variance), conditions[, c("true_theta", "iteration", "number_items", "number_answer_categories", "model", "estimator", "information_summary")], condition_vector)
    
    average_per_condition_true_minus2 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == -2), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == -2)]), "mean")
    sd_per_condition_true_minus2 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == -2), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == -2)]), "sd")
    
    average_per_condition_true_1 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 1), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 1)]), "mean")
    sd_per_condition_true_1 <- aggregate(estimates_and_conditions[which(estimates_and_conditions[,"true_theta"] == 1), "estimated_theta"], list(condition_vector[which(estimates_and_conditions[,"true_theta"] == 1)]), "sd") 
    
    # five number summary of average theta estimate per condition, with true theta is -2
    expect_equal(round(fivenum(average_per_condition_true_minus2[,"x"]), 3), c(-2.182, -2.049, -2.024, -1.999, -1.938))
    # five number summary of average theta estimate per condition, with true theta is 1
    expect_equal(round(fivenum(average_per_condition_true_1[,"x"]), 3), c(.934, .986, 1.003, 1.025, 1.084))
    
    # five number summary of observed sd of the theta estimates within each condition, with true theta is -2
    expect_equal(round(fivenum(sd_per_condition_true_minus2[,"x"]), 3), c(.156, .278, .322, .435, .548))  
    # five number summary of observed sd of the theta estimates within each condition, with true theta is 1
    expect_equal(round(fivenum(sd_per_condition_true_1[,"x"]), 3), c(.156, .227, .257, .346, .434))  
    
    # five number summary of reported sd of the theta estimate within each condition where max number of items is 50 and 100, respectively
    expect_equal(round(sqrt(fivenum(estimates_and_conditions[which(estimates_and_conditions[,"number_items"] == 50), "variance_estimate"])), 3), c(.188, .283, .351, .420, 1.300))
    expect_equal(round(sqrt(fivenum(estimates_and_conditions[which(estimates_and_conditions[,"number_items"] == 100), "variance_estimate"])), 3), c(.142, .200, .249, .298, .569))
    
  })
  
  test_that("three dimensions, estimates ML and MAP, information summary D, PD, A, and PA, no constraints on item selection, 100 iterations per condition", {
    # EAP estimation does not work
    # PEKL information summaries give errors because it also makes use of EAP estimation
    iterations_per_unique_condition <- 100 
    true_theta_vec <- c(-2, 1, 2)
    number_items_vec <- c(300)
    number_answer_categories_vec <- c(2, 4)
    number_dimensions <- 3
    
    start_items <- list(type = 'randomByDimension', nByDimension = 3, n = 9)
    variance_target <- .1^2
    model_vec <- c("3PLM","GRM","GPCM","SM")
    estimator_vec <- c("ML", "MAP") # AEP
    information_summary_vec <- c("D", "PD", "A", "PA") # PEKL 
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
    expect_equal(round(fivenum(average_per_condition_dim1[,"x"]), 3), c(-2.108, -2.045, -2.015, -1.992, -1.950))
    # five number summary of average theta estimate per condition, dimension 2 (true theta is 1)
    expect_equal(round(fivenum(average_per_condition_dim2[,"x"]), 3), c(.954, .990, 1.011, 1.019, 1.058))
    # five number summary of average theta estimate per condition, dimension 3 (true theta is 2)
    expect_equal(round(fivenum(average_per_condition_dim3[,"x"]), 3), c(1.967, 1.999, 2.016, 2.040, 2.123))
    
    # five number summary of observed sd of the theta estimates within each condition, for dimension 1, 2, and 3, respectively
    expect_equal(round(fivenum(sd_per_condition_dim1[,"x"]), 3), c(.176, .211, .304, .340, .387))
    expect_equal(round(fivenum(sd_per_condition_dim2[,"x"]), 3), c(.157, .195, .243, .260, .301))
    expect_equal(round(fivenum(sd_per_condition_dim3[,"x"]), 3), c(.171, .206, .312, .337, .385))
    
    # five number summary of reported sd of the theta estimate within each condition, for dimension 1, 2, and 3, respectively
    expect_equal(round(sqrt(fivenum(estimates_and_conditions[, "variance_estimate1"])), 3), c(.160, .262, .289, .316, .525))
    expect_equal(round(sqrt(fivenum(estimates_and_conditions[, "variance_estimate5"])), 3), c(.146, .206, .248, .277, .382))
    expect_equal(round(sqrt(fivenum(estimates_and_conditions[, "variance_estimate9"])), 3), c(.158, .212, .289, .316, .522))
    
    # no errors/missings for MAP estimator
    expect_equal(number_na_per_condition[which(conditions[seq(1, 6400, 100), "estimator"] == "MAP"), "x"], rep(0, 32))
    # some errors/missings for ML estimator because prior_var_safe_ml = NULL
    expect_equal(number_na_per_condition[which(conditions[seq(1, 6400, 100), "estimator"] == "ML"), "x"], 
                 c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0))
    
  })
  
  test_that("test prior_var_safe_ml is 100", {
    # EAP estimation does not work
    # PEKL information summaries give errors because it also makes use of EAP estimation
    # ML estimate gives error here
    iterations_per_unique_condition <- 100
    true_theta_vec <- c(-2, 1, 2)
    number_items_vec <- c(300)
    number_answer_categories_vec <- c(2, 4)
    number_dimensions <- 3
    
    start_items <- list(type = 'randomByDimension', nByDimension = 3, n = 9)
    variance_target <- .1^2
    model_vec <- c("3PLM","GRM","GPCM","SM")
    estimator_vec <- c("ML", "MAP") # AEP
    information_summary_vec <- c("D", "PD", "A", "PA") # PEKL 
    item_selection <- 'MI'
    
    estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection, start_items, variance_target, iterations_per_unique_condition, number_dimensions, prior_var_safe_ml = 100)
    
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
  
  
  test_that("items load on all dimensions prior_var_safe_ml is 100", {
    # find out where the occasional errors come from when estimator is ML
    # EAP estimation does not work
    # PEKL information summaries give errors because it also makes use of EAP estimation
    iterations_per_unique_condition <- 100 
    true_theta_vec <- c(2, -1, -2)
    number_items_vec <- c(300)
    number_answer_categories_vec <- c(2, 4)
    number_dimensions <- 3
    
    start_items <- list(type = 'randomByDimension', nByDimension = 3, n = 9)
    variance_target <- .1^2
    model_vec <- c("3PLM","GRM","GPCM","SM")
    estimator_vec <- c("ML", "MAP") # AEP
    information_summary_vec <- c("D", "PD", "A", "PA") # PEKL 
    item_selection <- 'MI'
    
    estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection, start_items, variance_target, iterations_per_unique_condition, number_dimensions, items_load_one_dimension = FALSE, prior_var_safe_ml = 100)
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
  
  test_that("items load on all dimensions prior_var_safe_ml is 1", {
    # find out where the occasional errors come from when estimator is ML
    # EAP estimation does not work
    # PEKL information summaries give errors because it also makes use of EAP estimation
    iterations_per_unique_condition <- 100 
    true_theta_vec <- c(2, -1, -2)
    number_items_vec <- c(300)
    number_answer_categories_vec <- c(2, 4)
    number_dimensions <- 3
    
    start_items <- list(type = 'randomByDimension', nByDimension = 3, n = 9)
    variance_target <- .1^2
    model_vec <- c("3PLM","GRM","GPCM","SM")
    estimator_vec <- c("ML", "MAP") # AEP
    information_summary_vec <- c("D", "PD", "A", "PA") # PEKL 
    item_selection <- 'MI'
    
    estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection, start_items, variance_target, iterations_per_unique_condition, number_dimensions, items_load_one_dimension = FALSE, prior_var_safe_ml = 1)
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
    # EAP estimation does not work
    # PEKL information summaries give errors because it also makes use of EAP estimation
    # ML estimate gives error here
    iterations_per_unique_condition <- 100
    true_theta_vec <- c(-2, 1, 2)
    number_items_vec <- 300 # can only have length one here (with item characteristics) and should be divisible by 3, to keep things simple
    number_answer_categories_vec <- 4
    number_dimensions <- 3
    
    start_items <- list(type = 'randomByDimension', nByDimension = 3, n = 9)
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
    
    estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection, start_items, variance_target, iterations_per_unique_condition, number_dimensions, constraints = list(characteristics = characteristics, constraints = constraints), prior_var_safe_ml = 100, return_administered_item_indeces = TRUE, max_n = max_n)
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
    # EAP estimation does not work
    # PEKL information summaries give errors because it also makes use of EAP estimation
    # ML estimate gives error here
    iterations_per_unique_condition <- 100
    true_theta_vec <- c(-2, 1, 2)
    number_items_vec <- 300 # can only have length one here (with item characteristics) and should be divisible by 3, to keep things simple
    number_answer_categories_vec <- 4
    number_dimensions <- 3
    
    start_items <- list(type = 'randomByDimension', nByDimension = 3, n = 9)
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
    
    estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection, start_items, variance_target, iterations_per_unique_condition, number_dimensions, constraints = list(characteristics = characteristics, constraints = constraints), prior_var_safe_ml = 100, return_administered_item_indeces = TRUE, max_n = max_n)
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
    # EAP estimation does not work
    # PEKL information summaries give errors because it also makes use of EAP estimation
    # ML estimate gives error here
    iterations_per_unique_condition <- 100
    true_theta_vec <- c(-2, 1, 2)
    number_items_vec <- c(300)
    number_answer_categories_vec <- c(2, 4)
    number_dimensions <- 3
    
    start_items <- list(type = 'randomByDimension', nByDimension = 3, n = 9)
    variance_target <- .1^2
    model_vec <- c("3PLM","GRM","GPCM","SM")
    estimator_vec <- "MAP"
    information_summary_vec <- c("D", "PD", "A", "PA") # PEKL 
    item_selection <- 'MI'
    
    estimates_and_variance <- with_random_seed(2, run_simulation)(true_theta_vec, number_items_vec, number_answer_categories_vec, model_vec, estimator_vec, information_summary_vec, item_selection, start_items, variance_target, iterations_per_unique_condition, number_dimensions, prior = diag(3) * .5, prior_var_safe_ml = 100)
    
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


