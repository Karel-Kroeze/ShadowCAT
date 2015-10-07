# only for whithin R:
'
library(testthat)
'
make_random_seed_exist <- rnorm(1)

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
  person <- initPerson(item_characteristics_shadowcat_format, prior = diag(item_characteristics_shadowcat_format$Q))
  new_response = NULL
  
  person_next_shadow_item <- 0 # just need a start value for the while loop
  while (!is.list(person_next_shadow_item)) {
    person_next_shadow_item <- with_random_seed(2, ShadowCAT_roqua)(new_response, person, initiated_test) 
    if (!is.list(person_next_shadow_item))
      new_response <- with_random_seed(2, answer)(person, initiated_test, indeces = person_next_shadow_item)$responses[tail(1)] # simulate person response on new item  
  }
  
  expect_equal(as.vector(round(person_next_shadow_item$estimate, 3)), -1.188)
  expect_equal(as.vector(round(attr(person_next_shadow_item$estimate, "variance"), 3)), .387)
  expect_equal(person_next_shadow_item$available, c(1:4, 6:7, 14:36, 38:40, 42:49))
  expect_equal(person_next_shadow_item$administered, c(10, 11, 9, 12, 13, 5, 37, 50, 41, 8))
  expect_equal(person_next_shadow_item$responses, c(rep(0, 6), 1, 1, 0, 0))
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
  person <- initPerson(item_characteristics_shadowcat_format, prior = diag(item_characteristics_shadowcat_format$Q))
  new_response = NULL
  
  person_next_shadow_item <- 0 # just need a start value for the while loop
  while (!is.list(person_next_shadow_item)) {
    person_next_shadow_item <- with_random_seed(2, ShadowCAT_roqua)(new_response, person, initiated_test) 
    if (!is.list(person_next_shadow_item))
      new_response <- with_random_seed(2, answer)(person, initiated_test, indeces = person_next_shadow_item)$responses[tail(1)] # simulate person response on new item  
  }
  
  expect_equal(as.vector(round(person_next_shadow_item$estimate, 3)), -.895)
  expect_equal(as.vector(round(attr(person_next_shadow_item$estimate, "variance"), 3)), .426)
  expect_equal(person_next_shadow_item$available, c(1:4, 6:8, 14:36, 38:49))
  expect_equal(person_next_shadow_item$administered, c(10, 11, 9, 12, 13, 5, 37, 50))
  expect_equal(person_next_shadow_item$responses, c(rep(0, 6), 1, 1))
})
