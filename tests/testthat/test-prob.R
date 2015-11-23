# only for whithin R:
'
library(testthat)
'

make_random_seed_exist <- rnorm(1)

context("3PLM model")

test_that("model is 3PLM, 1 dimension, 2 categories, estimator is MAP, deriv is false", {
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
  
  initiated_person <- initPerson(item_characteristics_shadowcat_format, prior = NULL)
  
  probabilities <- prob(initiated_test, person = NULL, theta = .5, deriv = FALSE, prior = NULL, items = NULL)
  
  expect_equal(round(probabilities$P[1,], 3), c(.293, .707))
  expect_equal(round(probabilities$P[40,], 3), c(.329, .671))
  expect_equal(dim(probabilities$P), c(50, 2))
  expect_equal(length(probabilities), 1)
})

test_that("model is 3PLM, 2 dimensions, 2 categories, estimator is ML, deriv is true", {
  # create item characteristics
  model <- "3PLM"
  number_items <- 50
  number_dimensions <- 2
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
  initiated_person <- initPerson(item_characteristics_shadowcat_format, prior = NULL)
  initiated_person$responses <- rep(c(1, 0), 25)
  
  probabilities <- prob(initiated_test, person = initiated_person, theta = NULL, deriv = TRUE, prior = NULL, items = NULL)
  
  expect_equal(round(probabilities$P[1,], 3), c(.290, .710))
  expect_equal(round(probabilities$P[40,], 3), c(.342, .658))
  expect_equal(dim(probabilities$P), c(50, 2))
  expect_equal(length(probabilities$LL), 1)
  expect_equal(dim(probabilities$d1), c(1, 2))
  expect_equal(dim(probabilities$d2), c(2, 2))
  expect_equal(length(probabilities), 4)
  
  expect_equal(round(probabilities$LL, 3), -61.592)
  expect_equal(round(probabilities$d1, 3), matrix(c(-5.115, -7.757), ncol = 2))
  expect_equal(round(probabilities$d2, 3), matrix(c(-2.986, -3.166, -3.166, -3.833), ncol = 2))

})

test_that("model is 3PLM, 1 dimension, 2 categories, estimator is MAP, deriv is true", {
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
  initiated_person <- initPerson(item_characteristics_shadowcat_format, prior = diag(item_characteristics_shadowcat_format$Q))
  initiated_person$responses <- rep(c(1, 0), 25)
  
  probabilities <- prob(test = initiated_test, person = initiated_person, theta = NULL, deriv = TRUE, prior = NULL, items = NULL)
  
  expect_equal(round(probabilities$P[1,], 3), c(.347, .653))
  expect_equal(round(probabilities$P[40,], 3), c(.376, .624))
  expect_equal(dim(probabilities$P), c(50, 2))
  expect_equal(dim(probabilities$LL), c(1, 1))
  expect_equal(dim(probabilities$d1), c(1, 1))
  expect_equal(dim(probabilities$d2), c(1, 1))
  expect_equal(length(probabilities), 4)
  
  expect_equal(round(probabilities$LL, 3), matrix(-47.383))
  expect_equal(round(probabilities$d1, 3), matrix(-2.446))
  expect_equal(round(probabilities$d2, 3), matrix(-5.916))

})

# test_that("model is 3PLM, 3 dimensions, 2 categories, estimator is EAP, deriv is true", {
#   # create item characteristics
#   model <- '3PLM'
#   number_items <- 50
#   number_dimensions <- 3
#   number_answer_categories <- 2 # can only be 2 for 3PLM model
#   guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
#   eta <- NULL # only relevant for GPCM model
#   
#   alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
#   beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
#   
#   item_characteristics_shadowcat_format <- initItembank(model = model, alpha = alpha, beta = beta, guessing = guessing, eta = eta, silent = TRUE)
#   
#   # get initiated test: adaptive test rules
#   initiated_test <- initTest(item_characteristics_shadowcat_format, 
#                              start = list(type = 'random', n = 5), 
#                              stop = list(type = 'length', n = 30),
#                              max_n = 50, # utter maximum
#                              estimator = 'EAP',
#                              objective = 'PD',
#                              selection = 'MI',
#                              constraints = NULL,
#                              exposure = NULL,
#                              lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
#                              upperBound = rep(3, item_characteristics_shadowcat_format$Q))
#   
#   # get initiated person
#   initiated_person <- initPerson(item_characteristics_shadowcat_format, prior = matrix(c(1.2, 1.5, 1.7, 1.5, .9, 1.5, 1.7, 1.5, 1.1), ncol = 3))
#   initiated_person$responses <- rep(c(1, 0), 25)
#   
#   probabilities <- prob(test = initiated_test, person = initiated_person, theta = NULL, deriv = TRUE, prior = NULL, items = NULL)
#   
#   expect_equal(round(probabilities$P[1,], 3), c(.204, .796))
#   expect_equal(round(probabilities$P[40,], 3), c(.305, .695))
#   expect_equal(dim(probabilities$P), c(50, 2))
#   expect_equal(dim(probabilities$LL), c(1, 1))
#   expect_equal(dim(probabilities$d1), c(1, 3))
#   expect_equal(dim(probabilities$d2), c(3, 3))
#   expect_equal(length(probabilities), 4)
# })

context("model is GPCM")

test_that("model is GPCM, 1 dimensions, 2 categories, estimator is MAP, deriv is true", {
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
  initiated_person <- initPerson(item_characteristics_shadowcat_format, prior = .6)
  initiated_person$responses <- rep(c(0, 1), 25)
  
  probabilities <- prob(initiated_test, person = initiated_person, theta = NULL, deriv = TRUE, prior = NULL, items = NULL)
  
  expect_equal(round(probabilities$P[1,], 3), c(.290, .710))
  expect_equal(round(probabilities$P[40,], 3), c(.439, .561))
  expect_equal(dim(probabilities$P), c(50, 2))
  expect_equal(dim(probabilities$LL), c(1, 1))
  expect_equal(dim(probabilities$d1), c(1, 1))
  expect_equal(dim(probabilities$d2), c(1, 1))
  expect_equal(length(probabilities), 4)
  
  expect_equal(round(probabilities$LL, 3), matrix(-33.772))
  expect_equal(round(probabilities$d1, 3), matrix(-1.824))
  expect_equal(round(probabilities$d2, 3), matrix(-11.325))
}) 

test_that("model is GPCM, 3 dimensions, 4 categories, estimator is EAP, deriv is true", {
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
  
  item_characteristics_shadowcat_format <- initItembank(model = model, alpha = alpha, beta = beta, guessing = guessing, eta = NULL, silent = TRUE)
  
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
  initiated_person <- initPerson(item_characteristics_shadowcat_format, prior = matrix(c(1.2, 1.5, 1.7, 1.5, .9, 1.5, 1.7, 1.5, 1.1), ncol = 3))
  initiated_person$responses = c(rep(c(0, 1, 2), 16), 1, 1)
  
  probabilities <- prob(initiated_test, person = initiated_person, theta = NULL, deriv = TRUE, prior = NULL, items = NULL)
  
  expect_equal(round(probabilities$P[1,], 3), c(.013, .231, .567, .188))
  expect_equal(round(probabilities$P[40,], 3), c(.041, .384, .491, .085))
  expect_equal(dim(probabilities$P), c(50, 4))
  expect_equal(dim(probabilities$LL), c(1, 1))
  expect_equal(dim(probabilities$d1), c(1, 3))
  expect_equal(dim(probabilities$d2), c(3, 3))
  expect_equal(length(probabilities), 4)
  
  expect_equal(round(probabilities$LL, 3), matrix(-101.251))
  expect_equal(round(probabilities$d1, 3), matrix(c(-20.747, -21.071, -21.938), ncol = 3))
  expect_equal(round(probabilities$d2[1,], 3), c(-20.739, -19.060, -19.274))
}) 

test_that("model is GPCM, 3 dimensions, varying number of categories, estimator is EAP, deriv is true", {
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
  beta <- row_cumsum(eta)
  
  item_characteristics_shadowcat_format <- initItembank(model = model, alpha = alpha, beta = beta, guessing = guessing, silent = TRUE)
  
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
  initiated_person <- initPerson(item_characteristics_shadowcat_format, theta = rep(.3, item_characteristics_shadowcat_format$Q), prior = matrix(c(1.2, 1.5, 1.7, 1.5, .9, 1.5, 1.7, 1.5, 1.1), ncol = 3))
  initiated_person$responses = c(rep(c(0, 1, 2), 16), 1, 1)
  
  probabilities <- prob(initiated_test, person = initiated_person, theta = NULL, deriv = TRUE, prior = NULL, items = NULL)
  
  expect_equal(round(probabilities$P[1,], 3), c(.009, .171, .819, NA,  NA))
  expect_equal(round(probabilities$P[40,], 3), c(.020, .191, .476, .313,  NA))
  expect_equal(round(probabilities$P[50,], 3), c(.002, .043, .276, .469, .211))
  expect_equal(dim(probabilities$P), c(50, 5))
  expect_equal(dim(probabilities$LL), c(1, 1))
  expect_equal(dim(probabilities$d1), c(1, 3))
  expect_equal(dim(probabilities$d2), c(3, 3))
  expect_equal(length(probabilities), 4)
  
  expect_equal(round(probabilities$LL, 3), matrix(-123.724))
  expect_equal(round(probabilities$d1, 3), matrix(c(-39.473, -39.336, -40.724), ncol = 3))
  expect_equal(round(probabilities$d2[1,], 3), c(-27.036, -23.968, -25.218))
}) 


test_that("model is GPCM, 3 dimensions, 4 categories, estimator is ML, deriv is true", {
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
  initiated_person <- initPerson(item_characteristics_shadowcat_format, prior = NULL)
  initiated_person$responses <- c(rep(c(0, 1, 2), 16), 1, 1)
  
  probabilities <- prob(initiated_test, person = initiated_person, theta = NULL, deriv = TRUE, prior = NULL, items = NULL)
  
  expect_equal(round(probabilities$P[1,], 3), c(.013, .231, .567, .188))
  expect_equal(round(probabilities$P[40,], 3), c(.041, .384, .491, .085))
  expect_equal(dim(probabilities$P), c(50, 4))
  expect_equal(length(probabilities$LL), 1)
  expect_equal(dim(probabilities$d1), c(1, 3))
  expect_equal(dim(probabilities$d2), c(3, 3))
  expect_equal(length(probabilities), 4)
}) 

test_that("model is GPCM, 3 dimensions, 4 categories, estimator is EAP, deriv is false", {
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
  
  item_characteristics_shadowcat_format <- initItembank(model = model, alpha = alpha, beta = beta, guessing = guessing, eta = NULL, silent = TRUE)
  
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
  initiated_person <- initPerson(item_characteristics_shadowcat_format, prior = matrix(c(1.2, 1.5, 1.7, 1.5, .9, 1.5, 1.7, 1.5, 1.1), ncol = 3))
  initiated_person$responses <- c(rep(c(0, 1, 2), 16), 1, 1)
  
  probabilities <- prob(initiated_test, person = initiated_person, theta = NULL, deriv = FALSE, prior = NULL, items = NULL)
  
  expect_equal(round(probabilities$P[1,], 3), c(.013, .231, .567, .188))
  expect_equal(round(probabilities$P[40,], 3), c(.041, .384, .491, .085))
  expect_equal(dim(probabilities$P), c(50, 4))
  expect_equal(length(probabilities), 1)
}) 

context("GRM model")

test_that("model is GRM, 3 dimensions, 4 categories, estimator is EAP, deriv is true", {
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
                             estimator = 'EAP',
                             objective = 'PD',
                             selection = 'MI',
                             constraints = NULL,
                             exposure = NULL,
                             lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                             upperBound = rep(3, item_characteristics_shadowcat_format$Q))
  
  # get initiated person
  initiated_person <- initPerson(item_characteristics_shadowcat_format, prior = matrix(c(1.2, 1.5, 1.7, 1.5, .9, 1.5, 1.7, 1.5, 1.1), ncol = 3))
  initiated_person$responses = c(rep(c(0, 1, 2), 50), 1, 1)
  
  probabilities <- prob(initiated_test, person = initiated_person, theta = NULL, deriv = TRUE, prior = NULL, items = NULL)
  initiated_person$responses <- c(rep(c(0, 1, 2), 16), 1, 1)
  
  expect_equal(round(probabilities$P[1,], 3), c(.052, .237, .461, .249))
  expect_equal(round(probabilities$P[40,], 3), c(.096, .343, .414, .148))
  expect_equal(dim(probabilities$P), c(50, 4))
  expect_equal(dim(probabilities$LL), c(1, 1))
  expect_equal(dim(probabilities$d1), c(1, 3))
  expect_equal(dim(probabilities$d2), c(3, 3))
  expect_equal(length(probabilities), 4)
  
  expect_equal(round(probabilities$LL, 3), matrix(-84.705))
  expect_equal(round(probabilities$d1, 3), matrix(c(-11.185, -10.826, -11.789), ncol = 3))
  expect_equal(round(probabilities$d2[1,], 3), c(-10.704, -10.179, -10.393))
}) 

test_that("model is GRM, 3 dimensions, 4 categories, estimator is ML, deriv is true", {
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
  initiated_person <- initPerson(item_characteristics_shadowcat_format, prior = NULL)
  initiated_person$responses <- c(rep(c(0, 1, 2), 16), 1, 1)
  
  probabilities <- prob(initiated_test, person = initiated_person, theta = NULL, deriv = TRUE, prior = NULL, items = NULL)
  
  expect_equal(round(probabilities$P[1,], 3), c(.052, .237, .461, .249))
  expect_equal(round(probabilities$P[40,], 3), c(.096, .343, .414, .148))
  expect_equal(dim(probabilities$P), c(50, 4))
  expect_equal(length(probabilities$LL), 1)
  expect_equal(dim(probabilities$d1), c(1, 3))
  expect_equal(dim(probabilities$d2), c(3, 3))
  expect_equal(length(probabilities), 4)
}) 

test_that("model is GRM, 3 dimensions, varying number of categories, estimator is EAP, deriv is true", {
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
                             estimator = 'EAP',
                             objective = 'PD',
                             selection = 'MI',
                             constraints = NULL,
                             exposure = NULL,
                             lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                             upperBound = rep(3, item_characteristics_shadowcat_format$Q))
  
  # get initiated person
  initiated_person <- initPerson(item_characteristics_shadowcat_format, theta = rep(.3, item_characteristics_shadowcat_format$Q), prior = matrix(c(1.2, 1.5, 1.7, 1.5, .9, 1.5, 1.7, 1.5, 1.1), ncol = 3))
  initiated_person$responses = c(rep(c(0, 1, 2), 16), 1, 1)
  
  probabilities <- prob(initiated_test, person = initiated_person, theta = NULL, deriv = TRUE, prior = NULL, items = NULL)
  
  expect_equal(round(probabilities$P[1,], 3), c(.052, .000, .989, NA, NA))
  expect_equal(round(probabilities$P[40,], 3), c(.096, .000, .020, .939, NA))
  expect_equal(round(probabilities$P[50,], 3), c(.039, .000, .000, .005, .992))
  expect_equal(dim(probabilities$P), c(50, 5))
  expect_equal(dim(probabilities$LL), c(1, 1))
  expect_equal(dim(probabilities$d1), c(1, 3))
  expect_equal(dim(probabilities$d2), c(3, 3))
  expect_equal(length(probabilities), 4)
  
  expect_equal(round(probabilities$LL, 3), matrix(-427.562))
  expect_equal(round(probabilities$d1, 3), matrix(c(-26.099, -25.069, -25.222), ncol = 3))
  expect_equal(round(probabilities$d2[1,], 3), c(-7.377, -7.087, -8.154))
}) 

test_that("model is GRM, 3 dimensions, 4 categories, estimator is EAP, deriv is false", {
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
                             estimator = 'EAP',
                             objective = 'PD',
                             selection = 'MI',
                             constraints = NULL,
                             exposure = NULL,
                             lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                             upperBound = rep(3, item_characteristics_shadowcat_format$Q))
  
  # get initiated person
  initiated_person <- initPerson(item_characteristics_shadowcat_format, prior = matrix(c(1.2, 1.5, 1.7, 1.5, .9, 1.5, 1.7, 1.5, 1.1), ncol = 3))
  
  probabilities <- prob(initiated_test, person = initiated_person, theta = NULL, deriv = FALSE, prior = NULL, items = NULL)
  initiated_person$responses <- c(rep(c(0, 1, 2), 16), 1, 1)
  
  expect_equal(round(probabilities$P[1,], 3), c(.052, .237, .461, .249))
  expect_equal(round(probabilities$P[40,], 3), c(.096, .343, .414, .148))
  expect_equal(dim(probabilities$P), c(50, 4))
  expect_equal(length(probabilities), 1)
}) 

context("SM model")

test_that("model is SM, 3 dimensions, 4 categories, estimator is EAP, deriv is true", {
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
                             estimator = 'EAP',
                             objective = 'PD',
                             selection = 'MI',
                             constraints = NULL,
                             exposure = NULL,
                             lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                             upperBound = rep(3, item_characteristics_shadowcat_format$Q))
  
  # get initiated person
  initiated_person <- initPerson(item_characteristics_shadowcat_format, prior = matrix(c(1.2, 1.5, 1.7, 1.5, .9, 1.5, 1.7, 1.5, 1.1), ncol = 3))
  initiated_person$responses = c(rep(c(0, 1, 2), 50), 1, 1)
  
  probabilities <- prob(initiated_test, person = initiated_person, theta = NULL, deriv = TRUE, prior = NULL, items = NULL)
  initiated_person$responses <- c(rep(c(0, 1, 2), 16), 1, 1)
  
  expect_equal(round(probabilities$P[1,], 3), c(.052, .275, .505, .168))
  expect_equal(round(probabilities$P[40,], 3), c(.096, .397, .433, .075))
  expect_equal(dim(probabilities$P), c(50, 4))
  expect_equal(dim(probabilities$LL), c(1, 1))
  expect_equal(dim(probabilities$d1), c(1, 3))
  expect_equal(dim(probabilities$d2), c(3, 3))
  expect_equal(length(probabilities), 4)
  
  expect_equal(round(probabilities$LL, 3), matrix(-84.888))
  expect_equal(round(probabilities$d1, 3), matrix(c(-7.74, -7.28, -8.297), ncol = 3))
  expect_equal(round(probabilities$d2[1,], 3), c(-12.905, -12.184, -12.239))
}) 

test_that("model is SM, 3 dimensions, varying number of categories, estimator is EAP, deriv is true", {
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
                             estimator = 'EAP',
                             objective = 'PD',
                             selection = 'MI',
                             constraints = NULL,
                             exposure = NULL,
                             lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                             upperBound = rep(3, item_characteristics_shadowcat_format$Q))
  
  # get initiated person
  initiated_person <- initPerson(item_characteristics_shadowcat_format, theta = rep(.3, item_characteristics_shadowcat_format$Q), prior = matrix(c(1.2, 1.5, 1.7, 1.5, .9, 1.5, 1.7, 1.5, 1.1), ncol = 3))
  initiated_person$responses <- c(rep(c(0, 1, 2), 16), 1, 1)
  
  probabilities <- prob(initiated_test, person = initiated_person, theta = NULL, deriv = TRUE, prior = NULL, items = NULL)
  
  expect_equal(round(probabilities$P[1,], 3), c(.052, .011, .937, NA, NA))
  expect_equal(round(probabilities$P[40,], 3), c(.096, .037, .053, .815, NA))
  expect_equal(round(probabilities$P[50,], 3), c(.039, 0.006, 0.004, 0.008, 0.943))
  expect_equal(dim(probabilities$P), c(50, 5))
  expect_equal(dim(probabilities$LL), c(1, 1))
  expect_equal(dim(probabilities$d1), c(1, 3))
  expect_equal(dim(probabilities$d2), c(3, 3))
  expect_equal(length(probabilities), 4)
  
  expect_equal(round(probabilities$LL, 3), matrix(-133.928))
  expect_equal(round(probabilities$d1, 3), matrix(c(-22.653, -21.523, -21.73), ncol = 3))
  expect_equal(round(probabilities$d2[1,], 3), c(-9.577, -9.092, -10.001))
}) 

test_that("model is SM, 3 dimensions, 4 categories, estimator is ML, deriv is true", {
  # create item characteristics
  model <- "SM"
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
  initiated_person <- initPerson(item_characteristics_shadowcat_format, prior = NULL)
  initiated_person$responses <- c(rep(c(0, 1, 2), 16), 1, 1)
  
  probabilities <- prob(initiated_test, person = initiated_person, theta = NULL, deriv = TRUE, prior = NULL, items = NULL)
  
  expect_equal(round(probabilities$P[1,], 3), c(.052, .275, .505, .168))
  expect_equal(round(probabilities$P[40,], 3), c(.096, .397, .433, .075))
  expect_equal(dim(probabilities$P), c(50, 4))
  expect_equal(length(probabilities$LL), 1)
  expect_equal(dim(probabilities$d1), c(1, 3))
  expect_equal(dim(probabilities$d2), c(3, 3))
  expect_equal(length(probabilities), 4)
}) 

test_that("model is SM, 3 dimensions, 4 categories, estimator is ML, deriv is false", {
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
  initiated_person <- initPerson(item_characteristics_shadowcat_format, prior = NULL)
  
  probabilities <- prob(initiated_test, person = initiated_person, theta = NULL, deriv = FALSE, prior = NULL, items = NULL)
  
  expect_equal(round(probabilities$P[1,], 3), c(.052, .275, .505, .168))
  expect_equal(round(probabilities$P[40,], 3), c(.096, .397, .433, .075))
  expect_equal(dim(probabilities$P), c(50, 4))
  expect_equal(length(probabilities), 1)
}) 
