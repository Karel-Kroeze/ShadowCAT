# only for whithin R:
'
library(testthat)
'

context("3PLM model")

test_that("model is 3PLM, 1 dimensions, 2 categories", {
  # create item characteristics
  model <- "3PLM"
  number_items <- 50
  number_dimensions <- 1
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  eta <- NULL # only relevant for GPCM model
  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  
  item_characteristics_shadowcat_format <- initItembank(model = model, alpha = alpha, beta = beta, guessing = guessing, eta = eta, silent = TRUE)
  
  # get initiated test: adaptive test rules
  initiated_test <- initTest(item_characteristics_shadowcat_format, 
                             start = list(type = 'random', n = 5), 
                             stop = list(type = 'length', n = 30),
                             max_n = 50, # utter maximum
                             estimator = 'ML',
                             objective = 'PD',
                             selection = 'MI',
                             constraints = NULL,
                             exposure = NULL,
                             lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                             upperBound = rep(3, item_characteristics_shadowcat_format$Q))
  
  # get initiated person
  initiated_person <- initPerson(item_characteristics_shadowcat_format, theta = rep(.2, item_characteristics_shadowcat_format$Q), prior = .4, responses = rep(c(1, 0), 25))
  
  fisher_information <- FI(initiated_test, initiated_person)
  
  expect_equal(dim(fisher_information), c(1, 1, 50))
  expect_equal(round(fisher_information[,,c(1,3, 40, 50)], 3), c(.055, .085, .039, .182))
})

test_that("model is 3PLM, 3 dimensions, 2 categories", {
  # create item characteristics
  model <- "3PLM"
  number_items <- 50
  number_dimensions <- 3
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  eta <- NULL # only relevant for GPCM model
  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  
  item_characteristics_shadowcat_format <- initItembank(model = model, alpha = alpha, beta = beta, guessing = guessing, eta = eta, silent = TRUE)
  
  # get initiated test: adaptive test rules
  initiated_test <- initTest(item_characteristics_shadowcat_format, 
                             start = list(type = 'random', n = 5), 
                             stop = list(type = 'length', n = 30),
                             max_n = 50, # utter maximum
                             estimator = 'ML',
                             objective = 'PD',
                             selection = 'MI',
                             constraints = NULL,
                             exposure = NULL,
                             lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                             upperBound = rep(3, item_characteristics_shadowcat_format$Q))
  
  # get initiated person
  initiated_person <- initPerson(item_characteristics_shadowcat_format, theta = rep(.2, item_characteristics_shadowcat_format$Q), prior = matrix(c(1.2, 1.5, 1.7, 1.5, .9, 1.5, 1.7, 1.5, 1.1), ncol = 3), responses = rep(c(1, 0), 25))
  
  fisher_information <- FI(initiated_test, initiated_person)
  
  expect_equal(dim(fisher_information), c(3, 3, 50))
  expect_equal(dim(fisher_information[,,2]), c(3, 3))
  expect_equal(round(fisher_information[,1,2], 3), c(.244, .068, .174))
  expect_equal(round(fisher_information[,2,3], 3), c(0, 0, 0))
  expect_equal(round(fisher_information[1,,24], 3), c(0.000, 0.000, 0.001))
  expect_equal(round(fisher_information[3,,40], 3), c(.063, .093, .102))
})

test_that("model is GPCM, 1 dimensions, 2 categories", {
  # create item characteristics
  model <- 'GPCM'
  number_items <- 50
  number_dimensions <- 1
  number_answer_categories <- 2
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  
  item_characteristics_shadowcat_format <- initItembank(model = model, alpha = alpha, beta = beta, guessing = guessing, eta = NULL, silent = TRUE)
  
  # get initiated test: adaptive test rules
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
  
  # get initiated person
  initiated_person <- initPerson(item_characteristics_shadowcat_format, theta = .5, prior = .6, responses = rep(c(0,1), 50))
  
  fisher_information <- FI(initiated_test, initiated_person)
  
  expect_equal(dim(fisher_information), c(1, 1, 50))
  expect_equal(round(fisher_information[,,c(1,3, 40, 50)], 3), c(.056, .138, .057, .288))
}) 


test_that("model is GPCM, 3 dimensions, varying numbers of categories", {
  # create item characteristics
  model <- 'GPCM'
  number_items <- 50
  number_dimensions <- 3
  max_number_answer_categories <- 5
  guessing <- NULL
  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  eta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (max_number_answer_categories - 1))))
  eta[c(1, 5:10), 3:4] <- NA
  eta[c(40:45), 4] <- NA
  beta <- NULL
  
  item_characteristics_shadowcat_format <- initItembank(model = model, alpha = alpha, beta = beta, guessing = guessing, eta = eta, silent = TRUE)
  
  # get initiated test: adaptive test rules
  initiated_test <- initTest(item_characteristics_shadowcat_format, 
                             start = list(type = 'random', n = 5), 
                             stop = list(type = 'length', n = 30),
                             max_n = 50, # utter maximum
                             estimator = 'ML',
                             objective = 'PD',
                             selection = 'MI',
                             constraints = NULL,
                             exposure = NULL,
                             lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                             upperBound = rep(3, item_characteristics_shadowcat_format$Q))
  
  # get initiated person
  initiated_person <- initPerson(item_characteristics_shadowcat_format, theta = rep(.2, item_characteristics_shadowcat_format$Q), prior = matrix(c(1.2, 1.5, 1.7, 1.5, .9, 1.5, 1.7, 1.5, 1.1), ncol = 3), responses = rep(c(1, 0), 25))
  
  fisher_information <- FI(initiated_test, initiated_person)
  
  expect_equal(dim(fisher_information), c(3, 3, 50))
  expect_equal(dim(fisher_information[,,2]), c(3, 3))
  expect_equal(round(fisher_information[,1,2], 3), c(.952, .265, .678))
  expect_equal(round(fisher_information[,2,3], 3), c(.643, .729, .962))
  expect_equal(round(fisher_information[1,,24], 3), c(.116, .126, .330))
  expect_equal(round(fisher_information[3,,40], 3), c(.216, .319, .351))
})
context("SM model")

test_that("model is SM 1 dimensions, 2 categories", {
  # create item characteristics
  model <- "SM"
  number_items <- 50
  number_dimensions <- 1
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  eta <- NULL # only relevant for GPCM model
  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  
  item_characteristics_shadowcat_format <- initItembank(model = model, alpha = alpha, beta = beta, guessing = guessing, eta = eta, silent = TRUE)
  
  # get initiated test: adaptive test rules
  initiated_test <- initTest(item_characteristics_shadowcat_format, 
                             start = list(type = 'random', n = 5), 
                             stop = list(type = 'length', n = 30),
                             max_n = 50, # utter maximum
                             estimator = 'ML',
                             objective = 'PD',
                             selection = 'MI',
                             constraints = NULL,
                             exposure = NULL,
                             lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                             upperBound = rep(3, item_characteristics_shadowcat_format$Q))
  
  # get initiated person
  initiated_person <- initPerson(item_characteristics_shadowcat_format, theta = rep(.2, item_characteristics_shadowcat_format$Q), prior = .4, responses = rep(c(1, 0), 25))
  
  fisher_information <- FI(initiated_test, initiated_person)
  
  expect_equal(dim(fisher_information), c(1, 1, 50))
  expect_equal(round(fisher_information[,,c(1,3, 40, 50)], 3), c(.056, .138, .057, .288))
})

test_that("model is SM, 3 dimensions, 3 dimensions, 4 categories", {
  # create item characteristics
  model <- 'SM'
  number_items <- 50
  number_dimensions <- 3
  number_answer_categories <- 4
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  beta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (number_answer_categories - 1))))
  
  item_characteristics_shadowcat_format <- initItembank(model = model, alpha = alpha, beta = beta, guessing = guessing, eta = NULL, silent = TRUE)
  
  # get initiated test: adaptive test rules
  initiated_test <- initTest(item_characteristics_shadowcat_format, 
                             start = list(type = 'random', n = 5), 
                             stop = list(type = 'length', n = 30),
                             max_n = 50, # utter maximum
                             estimator = 'ML',
                             objective = 'PD',
                             selection = 'MI',
                             constraints = NULL,
                             exposure = NULL,
                             lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                             upperBound = rep(3, item_characteristics_shadowcat_format$Q))
  
  # get initiated person
  initiated_person <- initPerson(item_characteristics_shadowcat_format, theta = rep(.2, item_characteristics_shadowcat_format$Q), prior = matrix(c(1.2, 1.5, 1.7, 1.5, .9, 1.5, 1.7, 1.5, 1.1), ncol = 3), responses = rep(c(1, 0), 25))
  
  fisher_information <- FI(initiated_test, initiated_person)
  
  expect_equal(dim(fisher_information), c(3, 3, 50))
  expect_equal(dim(fisher_information[,,2]), c(3, 3))
  expect_equal(round(fisher_information[,1,2], 3), c(.482, .134, .343))
  expect_equal(round(fisher_information[,2,3], 3), c(.362, .410, .541))
  expect_equal(round(fisher_information[1,,24], 3), c(.071, .077, .201))
  expect_equal(round(fisher_information[3,,40], 3), c(.140, .206, .227))
})

test_that("model is SM, 3 dimensions, varying number of categories", {
  # create item characteristics
  model <- 'SM'
  number_items <- 50
  number_dimensions <- 3
  max_number_answer_categories <- 5
  guessing <- NULL
  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  eta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (max_number_answer_categories - 1))))
  eta[c(1, 5:10), 3:4] <- NA
  eta[c(40:45), 4] <- NA
  beta <- row_cumsum(eta)
  
  item_characteristics_shadowcat_format <- initItembank(model = model, alpha = alpha, beta = beta, guessing = guessing, silent = TRUE)
  
  # get initiated test: adaptive test rules
  initiated_test <- initTest(item_characteristics_shadowcat_format, 
                             start = list(type = 'random', n = 5), 
                             stop = list(type = 'length', n = 30),
                             max_n = 50, # utter maximum
                             estimator = 'ML',
                             objective = 'PD',
                             selection = 'MI',
                             constraints = NULL,
                             exposure = NULL,
                             lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                             upperBound = rep(3, item_characteristics_shadowcat_format$Q))
  # get initiated person
  initiated_person <- initPerson(item_characteristics_shadowcat_format, theta = rep(.2, item_characteristics_shadowcat_format$Q), prior = matrix(c(1.2, 1.5, 1.7, 1.5, .9, 1.5, 1.7, 1.5, 1.1), ncol = 3), responses = rep(c(1, 0), 25))
  
  fisher_information <- FI(initiated_test, initiated_person)
  
  expect_equal(dim(fisher_information), c(3, 3, 50))
  expect_equal(dim(fisher_information[,,2]), c(3, 3))
  expect_equal(round(fisher_information[,1,2], 3), c(.589, .164, .419))
  expect_equal(round(fisher_information[,2,3], 3), c(.435, .494, .651))
  expect_equal(round(fisher_information[1,,24], 3), c(.079, .085, .223))
  expect_equal(round(fisher_information[3,,40], 3), c(.064, .095, 0.104))
})

context("GRM model")

test_that("model is GRM 1 dimensions, 2 categories", {
  # create item characteristics
  model <- "GRM"
  number_items <- 50
  number_dimensions <- 1
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  eta <- NULL # only relevant for GPCM model
  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  
  item_characteristics_shadowcat_format <- initItembank(model = model, alpha = alpha, beta = beta, guessing = guessing, eta = eta, silent = TRUE)
  
  # get initiated test: adaptive test rules
  initiated_test <- initTest(item_characteristics_shadowcat_format, 
                             start = list(type = 'random', n = 5), 
                             stop = list(type = 'length', n = 30),
                             max_n = 50, # utter maximum
                             estimator = 'ML',
                             objective = 'PD',
                             selection = 'MI',
                             constraints = NULL,
                             exposure = NULL,
                             lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                             upperBound = rep(3, item_characteristics_shadowcat_format$Q))
  
  # get initiated person
  initiated_person <- initPerson(item_characteristics_shadowcat_format, theta = rep(.2, item_characteristics_shadowcat_format$Q), prior = .4, responses = rep(c(1, 0), 25))
  
  fisher_information <- FI(initiated_test, initiated_person)
  
  expect_equal(dim(fisher_information), c(1, 1, 50))
  expect_equal(round(fisher_information[,,c(1,3, 40, 50)], 3), c(.056, .138, .057, .288))
})

test_that("model is GRM, 3 dimensions, varying numbers of categories", {
  # create item characteristics
  model <- 'GRM'
  number_items <- 50
  number_dimensions <- 3
  number_answer_categories <- 4
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  beta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (number_answer_categories - 1))))
  
  item_characteristics_shadowcat_format <- initItembank(model = model, alpha = alpha, beta = beta, guessing = guessing, eta = NULL, silent = TRUE)
  
  # get initiated test: adaptive test rules
  initiated_test <- initTest(item_characteristics_shadowcat_format, 
                             start = list(type = 'random', n = 5), 
                             stop = list(type = 'length', n = 30),
                             max_n = 50, # utter maximum
                             estimator = 'ML',
                             objective = 'PD',
                             selection = 'MI',
                             constraints = NULL,
                             exposure = NULL,
                             lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                             upperBound = rep(3, item_characteristics_shadowcat_format$Q))
  
  # get initiated person
  initiated_person <- initPerson(item_characteristics_shadowcat_format, theta = rep(.2, item_characteristics_shadowcat_format$Q), prior = matrix(c(1.2, 1.5, 1.7, 1.5, .9, 1.5, 1.7, 1.5, 1.1), ncol = 3), responses = rep(c(1, 0), 25))
  
  fisher_information <- FI(initiated_test, initiated_person)
  
  expect_equal(dim(fisher_information), c(3, 3, 50))
  expect_equal(dim(fisher_information[,,2]), c(3, 3))
  expect_equal(round(fisher_information[,1,2], 3), c(.385, .107, .274))
  expect_equal(round(fisher_information[,2,3], 3), c(.315, .357, .470))
  expect_equal(round(fisher_information[1,,24], 3), c(.063, .069, .180))
  expect_equal(round(fisher_information[3,,40], 3), c(.110, .163, .179))
})

test_that("model is GRM, 3 dimensions, varying number of categories", {
  # create item characteristics
  model <- 'GRM'
  number_items <- 50
  number_dimensions <- 3
  max_number_answer_categories <- 5
  guessing <- NULL
  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  eta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (max_number_answer_categories - 1))))
  eta[c(1, 5:10), 3:4] <- NA
  eta[c(40:45), 4] <- NA
  beta <- row_cumsum(eta)
  
  item_characteristics_shadowcat_format <- initItembank(model = model, alpha = alpha, beta = beta, guessing = guessing, silent = TRUE)
  
  # get initiated test: adaptive test rules
  initiated_test <- initTest(item_characteristics_shadowcat_format, 
                             start = list(type = 'random', n = 5), 
                             stop = list(type = 'length', n = 30),
                             max_n = 50, # utter maximum
                             estimator = 'ML',
                             objective = 'PD',
                             selection = 'MI',
                             constraints = NULL,
                             exposure = NULL,
                             lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                             upperBound = rep(3, item_characteristics_shadowcat_format$Q))
  # get initiated person
  initiated_person <- initPerson(item_characteristics_shadowcat_format, theta = rep(.2, item_characteristics_shadowcat_format$Q), prior = matrix(c(1.2, 1.5, 1.7, 1.5, .9, 1.5, 1.7, 1.5, 1.1), ncol = 3), responses = rep(c(1, 0), 25))
  
  fisher_information <- FI(initiated_test, initiated_person)
  
  expect_equal(dim(fisher_information), c(3, 3, 50))
  expect_equal(dim(fisher_information[,,2]), c(3, 3))
  expect_equal(round(fisher_information[,1,2], 3), c(.382, .106, .272))
  expect_equal(round(fisher_information[,2,3], 3), c(.330, .374, .493))
  expect_equal(round(fisher_information[1,,24], 3), c(.066, .071, .186))
  expect_equal(round(fisher_information[3,,40], 3), c(.024, .035, .039))
})

