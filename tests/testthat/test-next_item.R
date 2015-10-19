# only for whithin R:
'
library(testthat)
'

make_random_seed_exist <- rnorm(1)

context("number of responses is at least as large as number of required starting items")
  
test_that("number of responses is 5, number of required starting items 5", {
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
                             stop = list(type = 'length', n = 30),
                             max_n = 50, # utter maximum
                             estimator = 'MAP',
                             objective = 'PD',
                             selection = 'MI',
                             constraints = NULL,
                             exposure = NULL,
                             lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                             upperBound = rep(3, item_characteristics_shadowcat_format$Q))
  
  initiated_person <- initPerson(item_characteristics_shadowcat_format, theta = rep(.2, item_characteristics_shadowcat_format$Q), prior = .4)
  initiated_person$responses <- rep(1, 5)
  initiated_person$available <- c(3:10, 14:50)
  initiated_person$administered <- c(1, 3, 11, 12, 13)
  
  item_next <- next_item(initiated_person, initiated_test)
  expect_equal(item_next, 5)
})

test_that("number of responses is 7, number of required starting items 5", {
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
                             stop = list(type = 'length', n = 30),
                             max_n = 50, # utter maximum
                             estimator = 'MAP',
                             objective = 'PD',
                             selection = 'MI',
                             constraints = NULL,
                             exposure = NULL,
                             lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                             upperBound = rep(3, item_characteristics_shadowcat_format$Q))
  
  initiated_person <- initPerson(item_characteristics_shadowcat_format, theta = rep(.2, item_characteristics_shadowcat_format$Q), prior = .4)
  initiated_person$responses <- rep(1, 7)
  initiated_person$available <- c(3:4, 7:10, 14:50)
  initiated_person$administered <- c(1, 3, 5:6, 11, 12, 13)
  
  item_next <- next_item(initiated_person, initiated_test)
  expect_equal(item_next, 47)
})


context("start type is random")

test_that("number of responses is 0, number of required starting items 5, start type random", {
  initiated_test <- list(start = list(type = 'random', n = 5))
  initiated_person <- list(available = 1:50,
                           responses = numeric(0))
  item_next <- with_random_seed(2, next_item)(initiated_person, initiated_test)
  expect_equal(item_next, 10)
})

test_that("number of responses is 2, number of required starting items 5, start type random", {
  initiated_test <- list(start = list(type = 'random', n = 5))
  initiated_person <- list(available = 3:50,
                           responses = 1:2)
  item_next <- with_random_seed(2, next_item)(initiated_person, initiated_test)
  expect_equal(item_next, 11)
})


context("start type is fixed")

test_that("number of responses is 0, number of required starting items 5, start type fixed", {
  initiated_test <- list(start = list(type = 'fixed', indices = c(2, 4, 7:9), n = 5))
  initiated_person <- list(responses = numeric(0))
  item_next <- next_item(initiated_person, initiated_test)
  expect_equal(item_next, 2)
})

test_that("number of responses is 4, number of required starting items 5, start type fixed", {
  initiated_test <- list(start = list(type = 'fixed', indices = c(2, 4, 7:9), n = 5))
  initiated_person <- list(responses = c(2, 4, 7, 8))
  item_next <- next_item(initiated_person, initiated_test)
  expect_equal(item_next, 9)
})

context("start type is random by dimension")

test_that("number of responses is 0, number of required starting items 3 per dimension, random by dimension", {
  # create item characteristics
  model <- '3PLM'
  number_items <- 50
  number_dimensions <- 3
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  eta <- NULL # only relevant for GPCM model
  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  alpha[1:15,2:3] <- 0
  alpha[16:30,c(1,3)] <- 0
  alpha[31:50,1:2] <- 0
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  
  item_characteristics_shadowcat_format <- initItembank(model = model, alpha = alpha, beta = beta, guessing = guessing, eta = eta, silent = TRUE)
  
  # get initiated test
  initiated_test <- initTest(item_characteristics_shadowcat_format, 
                             start = list(type = 'randomByDimension', nByDimension = 3, n = 9), 
                             stop = list(type = 'length', n = 30),
                             max_n = 50, # utter maximum
                             estimator = 'MAP',
                             objective = 'PD',
                             selection = 'MI',
                             constraints = NULL,
                             exposure = NULL,
                             lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                             upperBound = rep(3, item_characteristics_shadowcat_format$Q))
  
  initiated_person <- list(available = 1:50,
                           responses = numeric(0))
  
  item_next <- with_random_seed(2, next_item)(initiated_person, initiated_test)
  expect_equal(item_next, 3)
})

test_that("number of responses is 6, number of required starting items 3 per dimension, random by dimension", {
  # create item characteristics
  model <- '3PLM'
  number_items <- 50
  number_dimensions <- 3
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  eta <- NULL # only relevant for GPCM model
  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  alpha[1:15,2:3] <- 0
  alpha[16:30,c(1,3)] <- 0
  alpha[31:50,1:2] <- 0
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  
  item_characteristics_shadowcat_format <- initItembank(model = model, alpha = alpha, beta = beta, guessing = guessing, eta = eta, silent = TRUE)
  
  # get initiated test
  initiated_test <- initTest(item_characteristics_shadowcat_format, 
                             start = list(type = 'randomByDimension', nByDimension = 3, n = 9), 
                             stop = list(type = 'length', n = 30),
                             max_n = 50, # utter maximum
                             estimator = 'MAP',
                             objective = 'PD',
                             selection = 'MI',
                             constraints = NULL,
                             exposure = NULL,
                             lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                             upperBound = rep(3, item_characteristics_shadowcat_format$Q))
  
  initiated_person <- list(available = c(2, 4:5, 7:15, 17:18, 20:26, 28:50),
                           responses = c(1, 3, 6, 16, 19, 27))
  
  item_next <- with_random_seed(2, next_item)(initiated_person, initiated_test)
  expect_equal(item_next, 34)
})

test_that("number of responses is 0, number of required starting items for each dimension is 2, 4, 2, random by dimension", {
  # create item characteristics
  model <- '3PLM'
  number_items <- 50
  number_dimensions <- 3
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  eta <- NULL # only relevant for GPCM model
  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  alpha[1:15,2:3] <- 0
  alpha[16:30,c(1,3)] <- 0
  alpha[31:50,1:2] <- 0
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  
  item_characteristics_shadowcat_format <- initItembank(model = model, alpha = alpha, beta = beta, guessing = guessing, eta = eta, silent = TRUE)
  
  # get initiated test
  initiated_test <- initTest(item_characteristics_shadowcat_format, 
                             start = list(type = 'randomByDimension', nByDimension = c(2, 4, 2), n = 8), 
                             stop = list(type = 'length', n = 30),
                             max_n = 50, # utter maximum
                             estimator = 'MAP',
                             objective = 'PD',
                             selection = 'MI',
                             constraints = NULL,
                             exposure = NULL,
                             lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                             upperBound = rep(3, item_characteristics_shadowcat_format$Q))
  
  initiated_person <- list(available = 1:50,
                           responses = numeric(0))
  
  item_next <- with_random_seed(2, next_item)(initiated_person, initiated_test)
  expect_equal(item_next, 3)
})

test_that("number of responses is 6, number of required starting items for each dimension is 2, 5, 1, random by dimension", {
  # create item characteristics
  model <- '3PLM'
  number_items <- 50
  number_dimensions <- 3
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  eta <- NULL # only relevant for GPCM model
  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  alpha[1:15,2:3] <- 0
  alpha[16:30,c(1,3)] <- 0
  alpha[31:50,1:2] <- 0
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  
  item_characteristics_shadowcat_format <- initItembank(model = model, alpha = alpha, beta = beta, guessing = guessing, eta = eta, silent = TRUE)
  
  # get initiated test
  initiated_test <- initTest(item_characteristics_shadowcat_format, 
                             start = list(type = 'randomByDimension', nByDimension = c(2, 5, 1), n = 8), 
                             stop = list(type = 'length', n = 30),
                             max_n = 50, # utter maximum
                             estimator = 'MAP',
                             objective = 'PD',
                             selection = 'MI',
                             constraints = NULL,
                             exposure = NULL,
                             lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                             upperBound = rep(3, item_characteristics_shadowcat_format$Q))
  
  initiated_person <- list(available = c(2, 4:5, 7:15, 17:18, 20:26, 28:50),
                           responses = c(1, 3, 6, 16, 19, 27))
  
  item_next <- with_random_seed(2, next_item)(initiated_person, initiated_test)
  expect_equal(item_next, 20)
})


test_that("items load on all three dimension while start type is nByDimension", {
  # create item characteristics
  model <- '3PLM'
  number_items <- 50
  number_dimensions <- 3
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  eta <- NULL # only relevant for GPCM model
  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  
  item_characteristics_shadowcat_format <- initItembank(model = model, alpha = alpha, beta = beta, guessing = guessing, eta = eta, silent = TRUE)
  
  # get initiated test
  initiated_test <- initTest(item_characteristics_shadowcat_format, 
                             start = list(type = 'randomByDimension', nByDimension = c(2, 5, 1), n = 8), 
                             stop = list(type = 'length', n = 30),
                             max_n = 50, # utter maximum
                             estimator = 'MAP',
                             objective = 'PD',
                             selection = 'MI',
                             constraints = NULL,
                             exposure = NULL,
                             lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                             upperBound = rep(3, item_characteristics_shadowcat_format$Q))
  
  initiated_person <- list(available = c(2, 4:5, 7:15, 17:18, 20:26, 28:50),
                           responses = c(1, 3, 6, 16, 19, 27))
  
  item_next <- with_random_seed(2, next_item)(initiated_person, initiated_test)
  expect_equal(item_next, 12)
})

context("invalid input")

test_that("sum nByDimension is not equal to n, same n for each dimension", {
  initiated_test <- list(start = list(type = 'randomByDimension', nByDimension = 2, n = 9),
                         items = list(Q = 3))
  initiated_person <- list(available = c(2, 4:5, 7:15, 17:18, 20:26, 28:50),
                           responses = c(1, 3, 6, 16, 19, 27))
  item_next <- next_item(initiated_person, initiated_test)
  expect_equal(item_next$errors$start, "contains inconsistent information. Total length of start phase and sum of length per dimension do not match")
})

test_that("sum nByDimension is not equal to n, different n for each dimension", {
  initiated_test <- list(start = list(type = 'randomByDimension', nByDimension = c(2, 5, 1), n = 9),
                         items = list(Q = 3))
  initiated_person <- list(available = c(2, 4:5, 7:15, 17:18, 20:26, 28:50),
                           responses = c(1, 3, 6, 16, 19, 27))
  item_next <- next_item(initiated_person, initiated_test)
  expect_equal(item_next$errors$start, "contains inconsistent information. Total length of start phase and sum of length per dimension do not match (n != sum(nByDimension)")
})


