# only for whithin R:
'
library(testthat)
'

context("stop rule is variance")

test_that("stop rule is variance, targets not reached, one dimension", {
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
  
  initiated_test <- initTest(item_characteristics_shadowcat_format,
                             start = list(type = 'fixed', indices = c(2, 4, 5), n = 3),
                             stop = list(type = 'variance', target = .2),
                             max_n = 20)
  
  # get initiated person
  initiated_person <- initPerson(item_characteristics_shadowcat_format, theta = rep(.2, item_characteristics_shadowcat_format$Q), prior = matrix(.4))
  initiated_person$responses <- rep(1, 15)
  
  should_stop <- stop_test(initiated_person, initiated_test)
  expect_equal(should_stop, FALSE)
})

test_that("stop rule is variance, targets not reached, three dimensions", {
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
  
  initiated_test <- initTest(item_characteristics_shadowcat_format,
                             start = list(type = 'fixed', indices = c(2, 4, 5), n = 3),
                             stop = list(type = 'variance', target = .2),
                             max_n = 20)
  
  # get initiated person
  initiated_person <- initPerson(item_characteristics_shadowcat_format, theta = rep(.2, item_characteristics_shadowcat_format$Q), prior = diag(c(.4, .5, .1)))
  initiated_person$responses <- rep(1, 15)
  
  should_stop <- stop_test(initiated_person, initiated_test)
  expect_equal(should_stop, FALSE)
})

test_that("stop rule is variance, variance reached, three dimensions", {
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
  
  initiated_test <- initTest(item_characteristics_shadowcat_format,
                             start = list(type = 'fixed', indices = c(2, 4, 5), n = 3),
                             stop = list(type = 'variance', target = .2),
                             max_n = 20)
  
  # get initiated person
  initiated_person <- initPerson(item_characteristics_shadowcat_format, theta = rep(.2, item_characteristics_shadowcat_format$Q), prior = diag(c(.19, .1, .1)))
  initiated_person$responses <- rep(1, 15)
  
  should_stop <- stop_test(initiated_person, initiated_test)
  expect_equal(should_stop, TRUE)
})

test_that("stop rule is variance, maximum number of items reached, three dimensions", {
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
  
  initiated_test <- initTest(item_characteristics_shadowcat_format,
                             start = list(type = 'fixed', indices = c(2, 4, 5), n = 3),
                             stop = list(type = 'variance', target = .2),
                             max_n = 20)
  
  # get initiated person
  initiated_person <- initPerson(item_characteristics_shadowcat_format, theta = rep(.2, item_characteristics_shadowcat_format$Q), prior = diag(c(.6, .1, .1)))
  initiated_person$responses <- rep(1, 20)
  
  should_stop <- stop_test(initiated_person, initiated_test)
  expect_equal(should_stop, TRUE)
})

context("stop rule is number of items")

test_that("stop rule is number of items, target not reached", {
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
  
  initiated_test <- initTest(item_characteristics_shadowcat_format,
                             start = list(type = 'fixed', indices = c(2, 4, 5), n = 3),
                             stop = list(type = 'length', n = 15),
                             max_n = 20)
  
  # get initiated person
  initiated_person <- initPerson(item_characteristics_shadowcat_format, theta = rep(.2, item_characteristics_shadowcat_format$Q), prior = diag(c(.1, .1, .1)))
  initiated_person$responses <- rep(1, 14)
  
  should_stop <- stop_test(initiated_person, initiated_test)
  expect_equal(should_stop, FALSE)
})

test_that("stop rule is number of items, target reached", {
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
  
  initiated_test <- initTest(item_characteristics_shadowcat_format,
                             start = list(type = 'fixed', indices = c(2, 4, 5), n = 3),
                             stop = list(type = 'length', n = 15),
                             max_n = 20)
  
  # get initiated person
  initiated_person <- initPerson(item_characteristics_shadowcat_format, theta = rep(.2, item_characteristics_shadowcat_format$Q), prior = diag(c(.1, .1, .1)))
  initiated_person$responses <- rep(1, 15)
  
  should_stop <- stop_test(initiated_person, initiated_test)
  expect_equal(should_stop, TRUE)
})
