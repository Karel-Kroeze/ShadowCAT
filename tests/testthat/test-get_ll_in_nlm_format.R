# only for whithin R:
'
library(testthat)
'

make_random_seed_exist <- rnorm(1)

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
initiated_person <- initPerson(item_characteristics_shadowcat_format, theta = rep(.2, item_characteristics_shadowcat_format$Q), prior = .4, responses = rep(c(1, 0), 17))
initiated_person$available <- c(1:5, 21:30, 50)
initiated_person$administered <- c(6:20, 31:49)

test_that("minimize is false, log is true", {    
  ll_in_nlm_format <- LL(theta = .3, initiated_test, initiated_person, minimize = FALSE, return_log = TRUE)
  
  expect_equal(round(as.vector(ll_in_nlm_format), 3), -25.917)
  expect_equal(round(attr(ll_in_nlm_format, "gradient"), 3), matrix(-5.396))
  expect_equal(round(attr(ll_in_nlm_format, "hessian"), 3), matrix(-5.354))
})

test_that("minimize is true, log is true", {
  ll_in_nlm_format <- LL(theta = .3, initiated_test, initiated_person, minimize = TRUE, return_log = TRUE)
  
  expect_equal(round(as.vector(ll_in_nlm_format), 3), 25.917)
  expect_equal(round(attr(ll_in_nlm_format, "gradient"), 3), matrix(5.396))
  expect_equal(round(attr(ll_in_nlm_format, "hessian"), 3), matrix(5.354))
})

test_that("minimize is true, log is false", {
  ll_in_nlm_format <- LL(theta = .3, initiated_test, initiated_person, minimize = TRUE, return_log = FALSE)
  
  expect_equal(round(as.vector(ll_in_nlm_format), 3), 180095565136)
  expect_equal(round(attr(ll_in_nlm_format, "gradient"), 3), matrix(5.396))
  expect_equal(round(attr(ll_in_nlm_format, "hessian"), 3), matrix(5.354))
})

