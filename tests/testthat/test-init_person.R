# only for whithin R:
'
library(testthat)
'

test_that("1 dimension, 2 categories", {
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
  
  # initiate person
  initiated_person <- initPerson(item_characteristics_shadowcat_format, theta = rep(0, item_characteristics_shadowcat_format$Q), prior = diag(item_characteristics_shadowcat_format$Q))
  
  expect_equal(initiated_person$theta, 0)
  expect_equal(initiated_person$prior, matrix(1))
  expect_equal(initiated_person$estimate[1], 0)
  expect_equal(attr(initiated_person$estimate, "variance"), matrix(1))
  expect_equal(initiated_person$available, 1:50)
  expect_equal(initiated_person$administered, numeric(0))
  expect_equal(initiated_person$responses, numeric(0)) 
})

test_that("3 dimensions, 4 categories", {
  # create item characteristics
  model <- 'GPCM'
  number_items <- 50
  number_dimensions <- 3
  number_answer_categories <- 4
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  eta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (number_answer_categories - 1))))
  beta <- t(apply(eta, 1, cumsum))
  
  item_characteristics_shadowcat_format <- initItembank(model = model, alpha = alpha, beta = beta, guessing = guessing, eta = eta, silent = TRUE)
  
  # initiate person
  initiated_person <- initPerson(item_characteristics_shadowcat_format, theta = rep(.2, item_characteristics_shadowcat_format$Q), prior = diag(item_characteristics_shadowcat_format$Q))
  
  expect_equal(initiated_person$theta, c(.2, .2, .2))
  expect_equal(initiated_person$prior, diag(3))
  expect_equal(initiated_person$estimate[1:3], c(0, 0, 0))
  expect_equal(attr(initiated_person$estimate, "variance"), diag(3))
  expect_equal(initiated_person$available, 1:50)
  expect_equal(initiated_person$administered, numeric(0))
  expect_equal(initiated_person$responses, numeric(0)) 
})