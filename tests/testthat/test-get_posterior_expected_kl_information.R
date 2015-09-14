# only for whithin R:
'
library(testthat)
'

make_random_seed_exist <- rnorm(1)

test_that("1 dimensions, 2 categories, estimator is ML", {
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
  initiated_person <- initPerson(item_characteristics_shadowcat_format, theta = rep(.2, item_characteristics_shadowcat_format$Q), prior = .4, responses = rep(c(1, 0), 17))
  initiated_person$available <- c(1:5, 21:30, 50)
  initiated_person$administered <- c(6:20, 31:49)
  
  posterior_expected_kl_information <- PEKL(initiated_test, initiated_person)
  
  expect_equal(length(posterior_expected_kl_information), 16)
  expect_equal(round(posterior_expected_kl_information[c(1, 3, 16)], 19), c(4.112e-13, 3.457e-13, 1.753e-12))
})

test_that("1 dimensions, 2 categories, estimator is EAP", {
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
                             estimator = 'EAP',
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
  
  posterior_expected_kl_information <- PEKL(initiated_test, initiated_person)
  
  expect_equal(length(posterior_expected_kl_information), 16)
  expect_equal(round(posterior_expected_kl_information[c(1, 3, 16)], 19), c(4.112e-13, 3.457e-13, 1.753e-12))
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
  initiated_person <- initPerson(item_characteristics_shadowcat_format, theta = rep(.2, item_characteristics_shadowcat_format$Q), responses = rep(c(1, 0), 17))
  initiated_person$available <- c(1:5, 21:30, 50)
  initiated_person$administered <- c(6:20, 31:49)
  
  posterior_expected_kl_information <- PEKL(initiated_test, initiated_person)
  
  expect_equal(length(posterior_expected_kl_information), 16)
  expect_equal(round(posterior_expected_kl_information[c(1, 3, 16)], 19), c(8.690e-17, 5.240e-17, 5.737e-16))
})
