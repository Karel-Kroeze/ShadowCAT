# only for whithin R:
'
library(testthat)
'
make_random_seed_exist <- rnorm(1)

test_that("validate", {
  # create item characteristics
  model <- '3PLM'
  number_items <- 50
  number_dimensions <- 1
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  eta <- NULL # only relevant for GPCM model
  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  
  item_characteristics_shadowcat_format <- initItembank(model = model, alpha = alpha, beta = beta, guessing = guessing, eta = eta, silent = TRUE)
  
  # get initiated test
  initiated_test <- initTest(item_characteristics_shadowcat_format, 
                             start = list(type = 'random', n = 1), 
                             stop = list(type = 'length', n = 50),
                             max_n = 50, # utter maximum
                             estimator = 'MAP',
                             objective = 'PD',
                             selection = 'MI',
                             constraints = NULL,
                             exposure = NULL,
                             lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                             upperBound = rep(3, item_characteristics_shadowcat_format$Q))
  
  # get initiated person
  initiated_person <- initPerson(item_characteristics_shadowcat_format, prior = diag(item_characteristics_shadowcat_format$Q))
  initiated_person$administered <- 1
  initiated_person$responses <- 1
  initiated_person$available <- 2:50
  
  new_response = NULL
  
  test_shadowcat_roqua <- function() {
    person_next_shadow_item <- 0 # just need a start value for the while loop
    while (!is.list(person_next_shadow_item)) {
      person_next_shadow_item <- with_random_seed(2, ShadowCAT_roqua)(new_response, initiated_person, initiated_test) 
      if (!is.list(person_next_shadow_item))
        new_response <- tail(answer(person_updated_after_new_response, initiated_test, indeces = person_next_shadow_item)$responses, 1) # simulate person response on new item  
    }
    person_next_shadow_item
  }
  
  person_roqua <- with_random_seed(2, test_shadowcat_roqua)()
  person_kroeze <- with_random_seed(2, ShadowCAT)(initiated_person, initiated_test)
  
  expect_equal(person_kroeze$estimate, person_roqua$estimate)
})

context("check stop rule")

test_that("stop rule is number of items", {
  # create item characteristics
  model <- '3PLM'
  number_items <- 50
  number_dimensions <- 1
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  eta <- NULL # only relevant for GPCM model
  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  
  item_characteristics_shadowcat_format <- initItembank(model = model, alpha = alpha, beta = beta, guessing = guessing, eta = eta, silent = TRUE)
  
  # get initiated test
  initiated_test <- initTest(item_characteristics_shadowcat_format, 
                             start = list(type = 'random', n = 5), 
                             stop = list(type = 'length', n = 10),
                             max_n = 50, # utter maximum
                             estimator = 'MAP',
                             objective = 'PD',
                             selection = 'MI',
                             constraints = NULL,
                             exposure = NULL,
                             lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                             upperBound = rep(3, item_characteristics_shadowcat_format$Q))
  
  # get initiated person
  initiated_person <- initPerson(item_characteristics_shadowcat_format, prior = diag(item_characteristics_shadowcat_format$Q))
  new_response = NULL
  
  test_shadowcat_roqua <- function() {
    person_next_shadow_item <- 0 # just need a start value for the while loop
    while (!is.list(person_next_shadow_item)) {
      person_next_shadow_item <- with_random_seed(2, ShadowCAT_roqua)(new_response, initiated_person, initiated_test) 
      if (!is.list(person_next_shadow_item))
        new_response <- tail(answer(person_updated_after_new_response, initiated_test, indeces = person_next_shadow_item)$responses, 1) # simulate person response on new item  
    }
    person_next_shadow_item
  }
  
  person_next_shadow_item <- with_random_seed(2, test_shadowcat_roqua)()
  
  expect_equal(as.vector(round(person_next_shadow_item$estimate, 3)), .219)
  expect_equal(as.vector(round(attr(person_next_shadow_item$estimate, "variance"), 3)), .276)
  expect_equal(person_next_shadow_item$available, c(1:4, 7, 14:40, 42:46, 48:50))
  expect_equal(person_next_shadow_item$administered, c(10, 11, 9, 12, 13, 5, 41, 47, 8, 6))
  expect_equal(person_next_shadow_item$responses, c(0, 1, 0, 0, 1, 1, 0, 1, 1, 1))
})

test_that("stop rule is variance", {
  # create item characteristics
  model <- '3PLM'
  number_items <- 50
  number_dimensions <- 1
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  eta <- NULL # only relevant for GPCM model
  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  
  item_characteristics_shadowcat_format <- initItembank(model = model, alpha = alpha, beta = beta, guessing = guessing, eta = eta, silent = TRUE)
  
  # get initiated test
  initiated_test <- initTest(item_characteristics_shadowcat_format, 
                             start = list(type = 'random', n = 5), 
                             stop = list(type = 'variance', target = .5),
                             max_n = 50, # utter maximum
                             estimator = 'MAP',
                             objective = 'PD',
                             selection = 'MI',
                             constraints = NULL,
                             exposure = NULL,
                             lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                             upperBound = rep(3, item_characteristics_shadowcat_format$Q))
  
  # get initiated person
  initiated_person <- initPerson(item_characteristics_shadowcat_format, prior = diag(item_characteristics_shadowcat_format$Q))
  new_response = NULL
  
  test_shadowcat_roqua <- function() {
    person_next_shadow_item <- 0 # just need a start value for the while loop
    while (!is.list(person_next_shadow_item)) {
      person_next_shadow_item <- with_random_seed(2, ShadowCAT_roqua)(new_response, initiated_person, initiated_test) 
      if (!is.list(person_next_shadow_item))
        new_response <- tail(answer(person_updated_after_new_response, initiated_test, indeces = person_next_shadow_item)$responses, 1) # simulate person response on new item  
    }
    person_next_shadow_item
  }
  
  person_next_shadow_item <- with_random_seed(2, test_shadowcat_roqua)()
  
  expect_equal(as.vector(round(person_next_shadow_item$estimate, 3)), -.255)
  expect_equal(as.vector(round(attr(person_next_shadow_item$estimate, "variance"), 3)), 0.465)
  expect_equal(person_next_shadow_item$available, c(1:4, 6:8, 14:50))
  expect_equal(person_next_shadow_item$administered, c(10, 11, 9, 12, 13, 5))
  expect_equal(person_next_shadow_item$responses, c(0, 1, 0, 0, 1, 1))
})



