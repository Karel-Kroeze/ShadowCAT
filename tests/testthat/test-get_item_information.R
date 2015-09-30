# only for whithin R:
'
library(testthat)
library(lpSolve)
'

make_random_seed_exist <- rnorm(1)

context("objective is D (determinant)")

test_that("First 1 item administered, objective is D (determinant), one dimension, with padding", {
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
                             max_n = 20, # utter maximum
                             estimator = 'MAP',
                             objective = 'D',
                             selection = 'MI',
                             constraints = NULL,
                             exposure = NULL,
                             lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                             upperBound = rep(3, item_characteristics_shadowcat_format$Q))
  
  # get initiated person
  initiated_person <- initPerson(item_characteristics_shadowcat_format, theta = rep(.2, item_characteristics_shadowcat_format$Q), prior = .4)
  initiated_person$responses <- 1
  initiated_person$available <- 2:50
  initiated_person$administered <- 1

  item_information <- objective(initiated_test, initiated_person, pad = TRUE)

  expect_equal(length(item_information), 50)
  expect_equal(round(item_information[c(1:5, 35:38,47:50)], 3), c(.000, .313, .139, .104, .477, .190, .236, .170, .123, .419, .142, .193, .237))
})

test_that("First 1 item administered, objective is D (determinant), three dimensions, with padding", {
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
                             max_n = 20, # utter maximum
                             estimator = 'MAP',
                             objective = 'D',
                             selection = 'MI',
                             constraints = NULL,
                             exposure = NULL,
                             lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                             upperBound = rep(3, item_characteristics_shadowcat_format$Q))
  
  # get initiated person
  initiated_person <- initPerson(item_characteristics_shadowcat_format, theta = rep(.5, item_characteristics_shadowcat_format$Q), prior = diag(c(.4, .8, 1.5)))
  initiated_person$responses <- 1
  initiated_person$available <- 2:50
  initiated_person$administered <- 1
  
  item_information <- objective(initiated_test, initiated_person, pad = TRUE)
  
  expect_equal(length(item_information), 50)
  expect_equal(round(item_information[c(1:5, 35:38,47:50)], 23), c(0, -4.050e-21, 0, -2.687e-20, -6.800e-22,  3.100e-22, -3.858e-20, -2.400e-22, -5.120e-21,  4.880e-21,  1.100e-22,  1.980e-21, -6.410e-21))
})

test_that("First 10 item administered, objective is D (determinant), three dimensions, with padding", {
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
                             max_n = 20, # utter maximum
                             estimator = 'MAP',
                             objective = 'D',
                             selection = 'MI',
                             constraints = NULL,
                             exposure = NULL,
                             lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                             upperBound = rep(3, item_characteristics_shadowcat_format$Q))
  
  # get initiated person
  initiated_person <- initPerson(item_characteristics_shadowcat_format, theta = rep(.5, item_characteristics_shadowcat_format$Q), prior = diag(c(.4, .8, 1.5)))
  initiated_person$responses <- rep(c(1, 0), 5)
  initiated_person$available <- 11:50
  initiated_person$administered <- 1:10
  
  item_information <- objective(initiated_test, initiated_person, pad = TRUE)
  
  expect_equal(length(item_information), 50)
  expect_equal(round(item_information[c(1:5, 35:38,47:50)], 3), c(.000, .000, .000, .000, .000, .036, .040, .033, .050, .045, .039, .047, .034))
})

context("objective is PD (posterior determinant)")

test_that("First 1 item administered, objective is PD (posterior determinant), one dimension, no padding", {
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
                             max_n = 20, # utter maximum
                             estimator = 'MAP',
                             objective = 'PD',
                             selection = 'MI',
                             constraints = NULL,
                             exposure = NULL,
                             lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                             upperBound = rep(3, item_characteristics_shadowcat_format$Q))
  
  # get initiated person
  initiated_person <- initPerson(item_characteristics_shadowcat_format, theta = rep(.2, item_characteristics_shadowcat_format$Q), prior = .4)
  initiated_person$responses <- 1
  initiated_person$available <- 2:50
  initiated_person$administered <- 1
  
  item_information <- objective(initiated_test, initiated_person, pad = FALSE)
  
  expect_equal(length(item_information), 49)
  expect_equal(round(item_information[c(1:5, 35:38,47:49)], 3), c(2.813, 2.639, 2.604, 2.977, 2.963, 2.736, 2.670, 2.623, 2.754, 2.642, 2.693, 2.737))
})

test_that("First 1 item administered, objective is PD (posterior determinant), three dimensions, no padding", {
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
                             max_n = 20, # utter maximum
                             estimator = 'MAP',
                             objective = 'PD',
                             selection = 'MI',
                             constraints = NULL,
                             exposure = NULL,
                             lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                             upperBound = rep(3, item_characteristics_shadowcat_format$Q))
  
  # get initiated person
  initiated_person <- initPerson(item_characteristics_shadowcat_format, theta = rep(.5, item_characteristics_shadowcat_format$Q), prior = diag(c(.4, .8, 1.5)))
  initiated_person$responses <- 1
  initiated_person$available <- 2:50
  initiated_person$administered <- 1
  
  item_information <- objective(initiated_test, initiated_person, pad = FALSE)
  
  expect_equal(length(item_information), 49)
  expect_equal(round(item_information[c(1:5, 35:38,47:49)], 3), c(2.906, 2.283, 2.483, 3.047, 3.685, 2.879, 2.293, 2.672, 2.754, 2.645, 2.856, 2.388))
})

test_that("First 10 item administered, objective is PD (posterior determinant), three dimensions, no padding", {
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
                             max_n = 20, # utter maximum
                             estimator = 'MAP',
                             objective = 'PD',
                             selection = 'MI',
                             constraints = NULL,
                             exposure = NULL,
                             lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                             upperBound = rep(3, item_characteristics_shadowcat_format$Q))
  
  # get initiated person
  initiated_person <- initPerson(item_characteristics_shadowcat_format, theta = rep(.5, item_characteristics_shadowcat_format$Q), prior = diag(c(.4, .8, 1.5)))
  initiated_person$responses <- rep(c(1, 0), 5)
  initiated_person$available <- 11:50
  initiated_person$administered <- 1:10
  
  item_information <- objective(initiated_test, initiated_person, pad = FALSE)
  
  expect_equal(length(item_information), 40)
  expect_equal(round(item_information[c(1:5, 35:40)], 3), c(9.370, 8.728, 9.485, 9.083, 8.596, 8.824, 8.596, 9.845, 9.095, 9.439, 8.723))
})

context("objective is A (trace)")

test_that("First 1 item administered, objective is A (trace), one dimension, with padding", {
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
                             max_n = 20, # utter maximum
                             estimator = 'MAP',
                             objective = 'A',
                             selection = 'MI',
                             constraints = NULL,
                             exposure = NULL,
                             lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                             upperBound = rep(3, item_characteristics_shadowcat_format$Q))
  
  # get initiated person
  initiated_person <- initPerson(item_characteristics_shadowcat_format, theta = rep(.2, item_characteristics_shadowcat_format$Q), prior = .4)
  initiated_person$responses <- 1
  initiated_person$available <- 2:50
  initiated_person$administered <- 1
  
  item_information <- objective(initiated_test, initiated_person, pad = TRUE)
  
  expect_equal(length(item_information), 50)
  expect_equal(round(item_information[c(1:5, 35:38,47:50)], 3), c(.000, .313, .139, .104, .477, .190, .236, .170, .123, .419, .142, .193, .237))
})

test_that("First 1 item administered, objective is A (trace), three dimensions, with padding", {
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
                             max_n = 20, # utter maximum
                             estimator = 'MAP',
                             objective = 'A',
                             selection = 'MI',
                             constraints = NULL,
                             exposure = NULL,
                             lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                             upperBound = rep(3, item_characteristics_shadowcat_format$Q))
  
  # get initiated person
  initiated_person <- initPerson(item_characteristics_shadowcat_format, theta = rep(.5, item_characteristics_shadowcat_format$Q), prior = diag(c(.4, .8, 1.5)))
  initiated_person$responses <- 1
  initiated_person$available <- 2:50
  initiated_person$administered <- 1
  
  item_information <- objective(initiated_test, initiated_person, pad = TRUE)
  
  expect_equal(length(item_information), 50)
  expect_equal(round(item_information[c(1:5, 35:38,47:50)], 3), c(.000, .488, .101, .190, .691, .266, .417, .107, .258, .684, .295, .401, .176))
})

test_that("First 10 item administered, objective is A (trace), three dimensions, with padding", {
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
                             max_n = 20, # utter maximum
                             estimator = 'MAP',
                             objective = 'A',
                             selection = 'MI',
                             constraints = NULL,
                             exposure = NULL,
                             lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                             upperBound = rep(3, item_characteristics_shadowcat_format$Q))
  
  # get initiated person
  initiated_person <- initPerson(item_characteristics_shadowcat_format, theta = rep(.5, item_characteristics_shadowcat_format$Q), prior = diag(c(.4, .8, 1.5)))
  initiated_person$responses <- rep(c(1, 0), 5)
  initiated_person$available <- 11:50
  initiated_person$administered <- 1:10
  
  item_information <- objective(initiated_test, initiated_person, pad = TRUE)
  
  expect_equal(length(item_information), 50)
  expect_equal(round(item_information[c(1:5, 35:38,47:50)], 3), c(.000, .000, .000, .000, .000, 3.678, 3.829, 3.520, 3.670, 4.096, 3.707, 3.813, 3.588))
})

context("objective is PA (posterior trace)")

test_that("First 1 item administered, objective is PA (posterior trace), one dimension, with padding", {
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
                             max_n = 20, # utter maximum
                             estimator = 'MAP',
                             objective = 'PA',
                             selection = 'MI',
                             constraints = NULL,
                             exposure = NULL,
                             lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                             upperBound = rep(3, item_characteristics_shadowcat_format$Q))
  
  # get initiated person
  initiated_person <- initPerson(item_characteristics_shadowcat_format, theta = rep(.2, item_characteristics_shadowcat_format$Q), prior = .4)
  initiated_person$responses <- 1
  initiated_person$available <- 2:50
  initiated_person$administered <- 1
  
  item_information <- objective(initiated_test, initiated_person, pad = TRUE)
  
  expect_equal(length(item_information), 50)
  expect_equal(round(item_information[c(1:5, 35:38,47:50)], 3), c(.000, 2.813, 2.639, 2.604, 2.977, 2.690, 2.736, 2.670, 2.623, 2.919, 2.642, 2.693, 2.737))
})

test_that("First 1 item administered, objective is PA (posterior trace), three dimensions, with padding", {
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
                             max_n = 20, # utter maximum
                             estimator = 'MAP',
                             objective = 'PA',
                             selection = 'MI',
                             constraints = NULL,
                             exposure = NULL,
                             lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                             upperBound = rep(3, item_characteristics_shadowcat_format$Q))
  
  # get initiated person
  initiated_person <- initPerson(item_characteristics_shadowcat_format, theta = rep(.5, item_characteristics_shadowcat_format$Q), prior = diag(c(.4, .8, 1.5)))
  initiated_person$responses <- 1
  initiated_person$available <- 2:50
  initiated_person$administered <- 1
  
  item_information <- objective(initiated_test, initiated_person, pad = TRUE)
  
  expect_equal(length(item_information), 50)
  expect_equal(round(item_information[c(1:5, 35:38,47:50)], 3), c(.000, 4.904, 4.518, 4.606, 5.107, 4.683, 4.834, 4.524, 4.675, 5.101, 4.711, 4.818, 4.592))
})

test_that("First 10 item administered, objective is PA (posterior trace), three dimensions, with padding", {
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
                             max_n = 20, # utter maximum
                             estimator = 'MAP',
                             objective = 'PA',
                             selection = 'MI',
                             constraints = NULL,
                             exposure = NULL,
                             lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                             upperBound = rep(3, item_characteristics_shadowcat_format$Q))
  
  # get initiated person
  initiated_person <- initPerson(item_characteristics_shadowcat_format, theta = rep(.5, item_characteristics_shadowcat_format$Q), prior = diag(c(.4, .8, 1.5)))
  initiated_person$responses <- rep(c(1, 0), 5)
  initiated_person$available <- 11:50
  initiated_person$administered <- 1:10
  
  item_information <- objective(initiated_test, initiated_person, pad = TRUE)
  
  expect_equal(length(item_information), 50)
  expect_equal(round(item_information[c(1:5, 35:38,47:50)], 3), c(.000, .000, .000, .000, .000, 8.095, 8.246, 7.936, 8.087, 8.513, 8.124, 8.230, 8.005))
})

context("objective is PEKL")

test_that("First 1 item administered, objective is PEKL, one dimension, with padding", {
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
                             max_n = 20, # utter maximum
                             estimator = 'MAP',
                             objective = 'PEKL',
                             selection = 'MI',
                             constraints = NULL,
                             exposure = NULL,
                             lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                             upperBound = rep(3, item_characteristics_shadowcat_format$Q))
  
  # get initiated person
  initiated_person <- initPerson(item_characteristics_shadowcat_format, theta = rep(.2, item_characteristics_shadowcat_format$Q), prior = .4)
  initiated_person$responses <- 1
  initiated_person$available <- 2:50
  initiated_person$administered <- 1
  
  item_information <- objective(initiated_test, initiated_person, pad = TRUE)
  
  expect_equal(length(item_information), 50)
  expect_equal(round(item_information[c(1:5, 35:38,47:50)], 3), c(.000, .049, .019, .010, .075, .025, .034, .023, .013, .063, .017, .026, .034))
})


context("many item characteristics equal")

test_that("objective is PD", {
  # create item characteristics
  model <- '3PLM'
  number_items <- 50
  number_dimensions <- 1
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  eta <- NULL # only relevant for GPCM model
  
  alpha <- matrix(rep(1, number_items * number_dimensions), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(rep(1, number_items), nrow = number_items, ncol = 1)
  
  item_characteristics_shadowcat_format <- initItembank(model = model, alpha = alpha, beta = beta, guessing = guessing, eta = eta, silent = TRUE)
  
  initiated_test <- initTest(item_characteristics_shadowcat_format,
                             start = list(type = 'fixed', indices = c(2, 4, 5), n = 3),
                             stop = list(type = 'variance', target = .2),
                             max_n = 20, # utter maximum
                             estimator = 'MAP',
                             objective = 'PD',
                             selection = 'MI',
                             constraints = NULL,
                             exposure = NULL,
                             lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                             upperBound = rep(3, item_characteristics_shadowcat_format$Q))
  
  # get initiated person
  initiated_person <- initPerson(item_characteristics_shadowcat_format, theta = rep(.2, item_characteristics_shadowcat_format$Q), prior = .4)
  initiated_person$responses <- c(1, 0)
  initiated_person$available <- 3:50
  initiated_person$administered <- 1:2
  
  item_information <- objective(initiated_test, initiated_person, pad = TRUE)
  
  expect_equal(length(item_information), 50)
  expect_equal(round(item_information[c(1:5, 35:38,47:50)], 3), c(.000, .000, 2.917, 2.917, 2.917, 2.880, 2.880, 2.880, 2.880, 2.880, 2.880, 2.880, 2.880))
})

test_that("objective is PEKL", {
  # create item characteristics
  model <- '3PLM'
  number_items <- 50
  number_dimensions <- 1
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  eta <- NULL # only relevant for GPCM model
  
  alpha <- matrix(rep(1, number_items * number_dimensions), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(rep(1, number_items), nrow = number_items, ncol = 1)
  
  item_characteristics_shadowcat_format <- initItembank(model = model, alpha = alpha, beta = beta, guessing = guessing, eta = eta, silent = TRUE)
  
  initiated_test <- initTest(item_characteristics_shadowcat_format,
                             start = list(type = 'fixed', indices = c(2, 4, 5), n = 3),
                             stop = list(type = 'variance', target = .2),
                             max_n = 20, # utter maximum
                             estimator = 'MAP',
                             objective = 'PEKL',
                             selection = 'MI',
                             constraints = NULL,
                             exposure = NULL,
                             lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                             upperBound = rep(3, item_characteristics_shadowcat_format$Q))
  
  # get initiated person
  initiated_person <- initPerson(item_characteristics_shadowcat_format, theta = rep(.2, item_characteristics_shadowcat_format$Q), prior = .4)
  initiated_person$responses <- c(1, 0)
  initiated_person$available <- 3:50
  initiated_person$administered <- 1:2
  
  item_information <- objective(initiated_test, initiated_person, pad = TRUE)
  
  expect_equal(length(item_information), 50)
  expect_equal(round(item_information[c(1:5, 35:38,47:50)], 3), c(.000, .000, .009, .009, .009, .007, .007, .007, .007, .007, .007, .007, .007))
})

context("invalid input")

test_that("objective is of unknown type", {
  # create item characteristics
  model <- '3PLM'
  number_items <- 50
  number_dimensions <- 1
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  eta <- NULL # only relevant for GPCM model
  
  alpha <- matrix(rep(1, number_items * number_dimensions), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(rep(1, number_items), nrow = number_items, ncol = 1)
  
  item_characteristics_shadowcat_format <- initItembank(model = model, alpha = alpha, beta = beta, guessing = guessing, eta = eta, silent = TRUE)
  
  initiated_test <- initTest(item_characteristics_shadowcat_format,
                             start = list(type = 'fixed', indices = c(2, 4, 5), n = 3),
                             stop = list(type = 'variance', target = .2),
                             max_n = 20, # utter maximum
                             estimator = 'MAP',
                             objective = 'PPP',
                             selection = 'MI',
                             constraints = NULL,
                             exposure = NULL,
                             lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                             upperBound = rep(3, item_characteristics_shadowcat_format$Q))
  
  # get initiated person
  initiated_person <- initPerson(item_characteristics_shadowcat_format, theta = rep(.2, item_characteristics_shadowcat_format$Q), prior = .4)
  initiated_person$responses <- c(1, 0)
  initiated_person$available <- 3:50
  initiated_person$administered <- 1:2
  
  item_information <- objective(initiated_test, initiated_person, pad = TRUE)
  
  expect_equal(item_information$errors$objective, "of unknown type")
})
