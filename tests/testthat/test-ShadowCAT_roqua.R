# only for whithin R:
'
library(testthat)
'
make_random_seed_exist <- rnorm(1)

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
  
  # initialize new response
  new_response = NULL
  
  test_shadowcat_roqua <- function() {
    next_item_and_test_outcome <- shadowcat_roqua(new_response, prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary)
   
    while (next_item_and_test_outcome$index_new_item != "stop_test") {
      person_updated_after_new_response$theta <- true_theta
      new_response <- tail(answer(person_updated_after_new_response, 
                                  initTest(items = initItembank(model, alpha, beta, guessing, eta, silent = TRUE), 
                                           start = start_items, 
                                           stop = stop_test,
                                           max_n = stop_test$n,
                                           estimator = estimator,
                                           objective = information_summary,
                                           selection = item_selection), 
                                  indeces = next_item_and_test_outcome$index_new_item)$responses, 1) # simulate person response on new item
        
      next_item_and_test_outcome <- shadowcat_roqua(new_response, prior, model, alpha, beta, guessing, eta, start_items, stop_test, estimator, information_summary)   
    }
    
    next_item_and_test_outcome$person_updated_after_new_response
  }
  
  # running 
  test_outcome <- with_random_seed(2, test_shadowcat_roqua)()
  
  expect_equal(as.vector(round(test_outcome$estimate, 3)), 1.912)
  expect_equal(as.vector(round(attr(test_outcome$estimate, "variance"), 3)), .112)
  expect_equal(test_outcome$available, numeric(0))
  expect_equal(length(test_outcome$administered), 100)
})


# do not run yet: should first be adapted to new structure in shadowcat_roqua
if (FALSE) {
test_that("validate with ShadowCAT Karel Kroeze", {
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
      person_next_shadow_item <- shadowcat_roqua(new_response, initiated_person, initiated_test) 
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
      person_next_shadow_item <- with_random_seed(2, shadowcat_roqua)(new_response, initiated_person, initiated_test) 
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
      person_next_shadow_item <- with_random_seed(2, shadowcat_roqua)(new_response, initiated_person, initiated_test) 
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

# not run because it takes too much time
if (FALSE) {
  context("problematic")

  # This one gets stuck after the 248st administered item. At that point, the estimated variance matrix becomes
  # non-invertible (has an eigenvalue of -5.530436e-18). See the problems in the trans function within init.quad().
  # Also, when running the 248 items, the resulting estimates are poor. 
  test_that("true theta is 2, -2, 3", {
    # create item characteristics
    model <- '3PLM'
    number_items <- 300
    number_dimensions <- 3
    number_answer_categories <- 2 # can only be 2 for 3PLM model
    guessing <- NULL
    eta <- NULL # only relevant for GPCM model
    
    alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
    alpha[1:100,2:3] <- 0
    alpha[101:200,c(1,3)] <- 0
    alpha[201:300,1:2] <- 0
    beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
    
    item_characteristics_shadowcat_format <- initItembank(model = model, alpha = alpha, beta = beta, guessing = guessing, eta = eta, silent = TRUE)
    
    # get initiated test
    initiated_test <- initTest(item_characteristics_shadowcat_format, 
                               start = list(type = 'random', n = 3), 
                               stop = list(type = 'length', n = 248),
                               max_n = 300, # utter maximum
                               estimator = 'EAP',
                               objective = 'PD',
                               selection = 'MI',
                               constraints = NULL,
                               exposure = NULL,
                               lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                               upperBound = rep(3, item_characteristics_shadowcat_format$Q))
    
    # get initiated person
    initiated_person <- initPerson(item_characteristics_shadowcat_format, theta = c(2, -2, 3), prior = diag(item_characteristics_shadowcat_format$Q) * 5)
    
    new_response = NULL
    
    test_shadowcat_roqua <- function() {
      person_next_shadow_item <- 0 # just need a start value for the while loop
      while (!is.list(person_next_shadow_item)) {
        person_next_shadow_item <- shadowcat_roqua(new_response, initiated_person, initiated_test)
        if (!is.list(person_next_shadow_item))
          new_response <- tail(answer(person_updated_after_new_response, initiated_test, indeces = person_next_shadow_item)$responses, 1) # simulate person response on new item  
      }
      person_next_shadow_item
    }
    
    person_next_shadow_item <- with_random_seed(2, test_shadowcat_roqua)()
    #person_kroeze <- with_random_seed(2, ShadowCAT)(initiated_person, initiated_test)
     
    expect_equal(as.vector(round(person_next_shadow_item$estimate, 3)), c(.984, .171, -.055))
    expect_equal(as.vector(round(attr(person_next_shadow_item$estimate, "variance"), 3))[1:3],c(.238, .259, .122))
    expect_equal(length(person_next_shadow_item$administered), 248)
  })
}
}
