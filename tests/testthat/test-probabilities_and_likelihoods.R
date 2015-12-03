# only for whithin R:
'
library(testthat)
'

make_random_seed_exist <- rnorm(1)

context("3PLM model")

test_that("model is 3PLM, 1 dimension, 2 categories, estimator is MAP, deriv is false", {
  theta <- .5
  model <- "3PLM"
  administered <- 1:50
  number_dimensions <- 1
  estimator <- "MAP"
 
  number_items <- 50
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  
  probabilities <- probabilities_and_likelihoods(theta, responses = NULL, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior = NULL, output = "probs")
  
  expect_equal(round(probabilities[1,], 3), c(.293, .707))
  expect_equal(round(probabilities[40,], 3), c(.329, .671))
  expect_equal(dim(probabilities), c(50, 2))
})

test_that("model is 3PLM, 2 dimensions, 2 categories, estimator is ML, deriv is true", {
  theta <- c(0, 0)
  model <- "3PLM"
  administered <- 1:50
  number_dimensions <- 2
  estimator <- "ML"
  responses <- rep(c(1, 0), 25)
  
  number_items <- 50
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  
  probabilities <- probabilities_and_likelihoods(theta, responses = responses, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior = NULL, output = "both")
    
  expect_equal(round(probabilities$probabilities[1,], 3), c(.290, .710))
  expect_equal(round(probabilities$probabilities[40,], 3), c(.342, .658))
  expect_equal(dim(probabilities$probabilities), c(50, 2))
  expect_equal(length(probabilities$likelihoods), 1)
  expect_equal(dim(attr(probabilities$likelihoods, "gradient")), c(1, 2))
  expect_equal(dim(attr(probabilities$likelihoods, "hessian")), c(2, 2))
  expect_equal(length(probabilities), 2)
  
  expect_equal(round(probabilities$likelihoods[1], 3), -61.592)
  expect_equal(round(attr(probabilities$likelihoods, "gradient"), 3), matrix(c(-5.115, -7.757), ncol = 2))
  expect_equal(round(attr(probabilities$likelihoods, "hessian"), 3), matrix(c(-2.986, -3.166, -3.166, -3.833), ncol = 2))
})

test_that("model is 3PLM, 1 dimension, 2 categories, estimator is MAP, deriv is true", {
  theta <- 0
  model <- "3PLM"
  administered <- 1:50
  number_dimensions <- 1
  estimator <- "MAP"
  responses <- rep(c(1, 0), 25)
  prior <- diag(1)
  
  number_items <- 50
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  
  probabilities <- probabilities_and_likelihoods(theta, responses, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior = prior, output = "both")
    
  expect_equal(round(probabilities$probabilities[1,], 3), c(.347, .653))
  expect_equal(round(probabilities$probabilities[40,], 3), c(.376, .624))
  expect_equal(dim(probabilities$probabilities), c(50, 2))
  expect_equal(dim(probabilities$likelihoods), c(1, 1))
  expect_equal(dim(attr(probabilities$likelihoods, "gradient")), c(1, 1))
  expect_equal(dim(attr(probabilities$likelihoods, "hessian")), c(1, 1))
  expect_equal(length(probabilities), 2)
  
  expect_equal(round(probabilities$likelihoods[1], 3), -47.383)
  expect_equal(round(attr(probabilities$likelihoods, "gradient"), 3), matrix(-2.446))
  expect_equal(round(attr(probabilities$likelihoods, "hessian"), 3), matrix(-5.916))
})

test_that("model is 3PLM, 3 dimensions, 2 categories, estimator is EAP, deriv is true", {
  theta <- c(0, 0, 0)
  model <- "3PLM"
  administered <- 1:50
  number_dimensions <- 3
  estimator <- "EAP"
  responses <- rep(c(1, 0), 25)
  prior <- matrix(c(1.2, 1.5, 1.7, 1.5, .9, 1.5, 1.7, 1.5, 1.1), ncol = 3)
  
  number_items <- 50
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  
  probabilities <- probabilities_and_likelihoods(theta, responses, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior = prior, output = "both")
  
  expect_equal(round(probabilities$probabilities[1,], 3), c(.204, .796))
  expect_equal(round(probabilities$probabilities[40,], 3), c(.305, .695))
  expect_equal(dim(probabilities$probabilities), c(50, 2))
  expect_equal(dim(probabilities$likelihoods), c(1, 1))
  expect_equal(dim(attr(probabilities$likelihoods, "gradient")), c(1, 3))
  expect_equal(dim(attr(probabilities$likelihoods, "hessian")), c(3, 3))
  expect_equal(length(probabilities), 2)
})

context("model is GPCM")

test_that("model is GPCM, 1 dimensions, 2 categories, estimator is MAP, deriv is true", {
  theta <- 0
  model <- "GPCM"
  administered <- 1:50
  number_dimensions <- 1
  estimator <- "MAP"
  responses <- rep(c(0, 1), 25)
  prior <- diag(1) * .6
  
  number_items <- 50
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  
  probabilities <- probabilities_and_likelihoods(theta, responses, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior = prior, output = "both")
    
  expect_equal(round(probabilities$probabilities[1,], 3), c(.290, .710))
  expect_equal(round(probabilities$probabilities[40,], 3), c(.439, .561))
  expect_equal(dim(probabilities$probabilities), c(50, 2))
  expect_equal(dim(probabilities$likelihoods), c(1, 1))
  expect_equal(dim(attr(probabilities$likelihoods, "gradient")), c(1, 1))
  expect_equal(dim(attr(probabilities$likelihoods, "hessian")), c(1, 1))
  expect_equal(length(probabilities), 2)
  
  expect_equal(round(probabilities$likelihoods[1], 3), -33.772)
  expect_equal(round(attr(probabilities$likelihoods, "gradient"), 3), matrix(-1.824))
  expect_equal(round(attr(probabilities$likelihoods, "hessian"), 3), matrix(-11.325))
}) 

test_that("model is GPCM, 3 dimensions, 4 categories, estimator is EAP, deriv is true", {
  theta <- c(0, 0, 0)
  model <- "GPCM"
  administered <- 1:50
  number_dimensions <- 3
  estimator <- "EAP"
  responses <- c(rep(c(0, 1, 2), 16), 1, 1)
  prior <-  matrix(c(1.2, 1.5, 1.7, 1.5, .9, 1.5, 1.7, 1.5, 1.1), ncol = 3)
  
  number_items <- 50
  number_answer_categories <- 4
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  eta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (number_answer_categories - 1))))
  beta <- t(apply(eta, 1, cumsum))
  guessing <- NULL
  
  probabilities <- probabilities_and_likelihoods(theta, responses, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior = prior, output = "both")
  
  expect_equal(round(probabilities$probabilities[1,], 3), c(.013, .231, .567, .188))
  expect_equal(round(probabilities$probabilities[40,], 3), c(.041, .384, .491, .085))
  expect_equal(dim(probabilities$probabilities), c(50, 4))
  expect_equal(dim(probabilities$likelihoods), c(1, 1))
  expect_equal(dim(attr(probabilities$likelihoods, "gradient")), c(1, 3))
  expect_equal(dim(attr(probabilities$likelihoods, "hessian")), c(3, 3))
  expect_equal(length(probabilities), 2)
  
  expect_equal(round(probabilities$likelihoods[1], 3), -101.251)
  expect_equal(round(attr(probabilities$likelihoods, "gradient"), 3), matrix(c(-20.747, -21.071, -21.938), ncol = 3))
  expect_equal(round(attr(probabilities$likelihoods, "hessian")[1,], 3), c(-20.739, -19.060, -19.274))
}) 

test_that("model is GPCM, 3 dimensions, varying number of categories, estimator is EAP, deriv is true", {
  theta <- c(0, 0, 0)
  model <- "GPCM"
  administered <- 1:50
  number_dimensions <- 3
  estimator <- "EAP"
  responses <- c(rep(c(0, 1, 2), 16), 1, 1)
  prior <-  matrix(c(1.2, 1.5, 1.7, 1.5, .9, 1.5, 1.7, 1.5, 1.1), ncol = 3)
  
  number_items <- 50
  max_number_answer_categories <- 5
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  eta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (max_number_answer_categories - 1))))
  eta[c(1, 5:10), 3:4] <- NA
  eta[c(40:45), 4] <- NA
  beta <- row_cumsum(eta)
  guessing <- NULL
  
  probabilities <- probabilities_and_likelihoods(theta, responses, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior = prior, output = "both")
  
  expect_equal(round(probabilities$probabilities[1,], 3), c(.009, .171, .819, NA,  NA))
  expect_equal(round(probabilities$probabilities[40,], 3), c(.020, .191, .476, .313,  NA))
  expect_equal(round(probabilities$probabilities[50,], 3), c(.002, .043, .276, .469, .211))
  expect_equal(dim(probabilities$probabilities), c(50, 5))
  expect_equal(dim(probabilities$likelihoods), c(1, 1))
  expect_equal(dim(attr(probabilities$likelihoods, "gradient")), c(1, 3))
  expect_equal(dim(attr(probabilities$likelihoods, "hessian")), c(3, 3))
  expect_equal(length(probabilities), 2)
  
  expect_equal(round(probabilities$likelihoods[1], 3), -123.724)
  expect_equal(round(attr(probabilities$likelihoods, "gradient"), 3), matrix(c(-39.473, -39.336, -40.724), ncol = 3))
  expect_equal(round(attr(probabilities$likelihoods, "hessian")[1,], 3), c(-27.036, -23.968, -25.218))
}) 


test_that("model is GPCM, 3 dimensions, 4 categories, estimator is ML, deriv is true", {
  theta <- c(0, 0, 0)
  model <- "GPCM"
  administered <- 1:50
  number_dimensions <- 3
  estimator <- "ML"
  responses <- c(rep(c(0, 1, 2), 16), 1, 1)
  prior <- NULL
  
  number_items <- 50
  number_answer_categories <- 4
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  eta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (number_answer_categories - 1))))
  beta <- t(apply(eta, 1, cumsum))
  guessing <- NULL
  
  probabilities <- probabilities_and_likelihoods(theta, responses, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior = prior, output = "both")
  
  expect_equal(round(probabilities$probabilities[1,], 3), c(.013, .231, .567, .188))
  expect_equal(round(probabilities$probabilities[40,], 3), c(.041, .384, .491, .085))
  expect_equal(dim(probabilities$probabilities), c(50, 4))
  expect_equal(length(probabilities$likelihoods), 1)
  expect_equal(dim(attr(probabilities$likelihoods, "gradient")), c(1, 3))
  expect_equal(dim(attr(probabilities$likelihoods, "hessian")), c(3, 3))
  expect_equal(length(probabilities), 2)
}) 

test_that("model is GPCM, 3 dimensions, 4 categories, estimator is EAP, deriv is false", {
  theta <- c(0, 0, 0)
  model <- "GPCM"
  administered <- 1:50
  number_dimensions <- 3
  estimator <- "EAP"
  responses <- c(rep(c(0, 1, 2), 16), 1, 1)
  prior <- NULL
  
  number_items <- 50
  number_answer_categories <- 4
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  eta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (number_answer_categories - 1))))
  beta <- t(apply(eta, 1, cumsum))
  guessing <- NULL
  
  probabilities <- probabilities_and_likelihoods(theta, responses, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior = prior, output = "probs")
  
  expect_equal(round(probabilities[1,], 3), c(.013, .231, .567, .188))
  expect_equal(round(probabilities[40,], 3), c(.041, .384, .491, .085))
  expect_equal(dim(probabilities), c(50, 4))
}) 

context("GRM model")

test_that("model is GRM, 3 dimensions, 4 categories, estimator is EAP, deriv is true", {
  theta <- c(0, 0, 0)
  model <- "GRM"
  administered <- 1:50
  number_dimensions <- 3
  estimator <- "EAP"
  responses <- c(rep(c(0, 1, 2), 16), 0, 1)
  prior <- matrix(c(1.2, 1.5, 1.7, 1.5, .9, 1.5, 1.7, 1.5, 1.1), ncol = 3)
  
  number_items <- 50
  number_answer_categories <- 4
  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  beta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (number_answer_categories - 1))))
  guessing <- NULL
  
  probabilities <- probabilities_and_likelihoods(theta, responses, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior = prior, output = "both")
  
  expect_equal(round(probabilities$probabilities[1,], 3), c(.052, .237, .461, .249))
  expect_equal(round(probabilities$probabilities[40,], 3), c(.096, .343, .414, .148))
  expect_equal(dim(probabilities$probabilities), c(50, 4))
  expect_equal(dim(probabilities$likelihoods), c(1, 1))
  expect_equal(dim(attr(probabilities$likelihoods, "gradient")), c(1, 3))
  expect_equal(dim(attr(probabilities$likelihoods, "hessian")), c(3, 3))
  expect_equal(length(probabilities), 2)
  
  expect_equal(round(probabilities$likelihoods[1], 3), -84.705)
  expect_equal(round(attr(probabilities$likelihoods, "gradient"), 3), matrix(c(-11.185, -10.826, -11.789), ncol = 3))
  expect_equal(round(attr(probabilities$likelihoods, "hessian")[1,], 3), c(-10.704, -10.179, -10.393))
}) 

test_that("model is GRM, 3 dimensions, 4 categories, estimator is ML, deriv is true", {
  theta <- c(0, 0, 0)
  model <- "GRM"
  administered <- 1:50
  number_dimensions <- 3
  estimator <- "ML"
  responses <- c(rep(c(0, 1, 2), 16), 1, 1)
  prior <- NULL
  
  number_items <- 50
  number_answer_categories <- 4
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  beta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (number_answer_categories - 1))))
  guessing <- NULL
  
  probabilities <- probabilities_and_likelihoods(theta, responses, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior = prior, output = "both")
  
  expect_equal(round(probabilities$probabilities[1,], 3), c(.052, .237, .461, .249))
  expect_equal(round(probabilities$probabilities[40,], 3), c(.096, .343, .414, .148))
  expect_equal(dim(probabilities$probabilities), c(50, 4))
  expect_equal(length(probabilities$likelihoods), 1)
  expect_equal(dim(attr(probabilities$likelihoods, "gradient")), c(1, 3))
  expect_equal(dim(attr(probabilities$likelihoods, "hessian")), c(3, 3))
  expect_equal(length(probabilities), 2)
}) 

test_that("model is GRM, 3 dimensions, varying number of categories, estimator is EAP, deriv is true", {
  theta <- c(0, 0, 0)
  model <- "GRM"
  administered <- 1:50
  number_dimensions <- 3
  estimator <- "EAP"
  responses <- c(rep(c(0, 1, 2), 16), 1, 1)
  prior <- matrix(c(1.2, 1.5, 1.7, 1.5, .9, 1.5, 1.7, 1.5, 1.1), ncol = 3)
  
  number_items <- 50
  max_number_answer_categories <- 5
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  eta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (max_number_answer_categories - 1))))
  eta[c(1, 5:10), 3:4] <- NA
  eta[c(40:45), 4] <- NA
  beta <- row_cumsum(eta)
  guessing <- NULL
  
  probabilities <- probabilities_and_likelihoods(theta, responses, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior = prior, output = "both")
  
  expect_equal(round(probabilities$probabilities[1,], 3), c(.052, .000, .989, NA, NA))
  expect_equal(round(probabilities$probabilities[40,], 3), c(.096, .000, .020, .939, NA))
  expect_equal(round(probabilities$probabilities[50,], 3), c(.039, .000, .000, .005, .992))
  expect_equal(dim(probabilities$probabilities), c(50, 5))
  expect_equal(dim(probabilities$likelihoods), c(1, 1))
  expect_equal(dim(attr(probabilities$likelihoods, "gradient")), c(1, 3))
  expect_equal(dim(attr(probabilities$likelihoods, "hessian")), c(3, 3))
  expect_equal(length(probabilities), 2)
  
  expect_equal(round(probabilities$likelihoods[1], 3), -427.562)
  expect_equal(round(attr(probabilities$likelihoods, "gradient"), 3), matrix(c(-26.099, -25.069, -25.222), ncol = 3))
  expect_equal(round(attr(probabilities$likelihoods, "hessian")[1,], 3), c(-7.377, -7.087, -8.154))
}) 

test_that("model is GRM, 3 dimensions, 4 categories, estimator is EAP, deriv is false", {
  theta <- c(0, 0, 0)
  model <- "GRM"
  administered <- 1:50
  number_dimensions <- 3
  estimator <- "EAP"
  responses <- c(rep(c(0, 1, 2), 16), 1, 1)
  prior <- NULL
  
  number_items <- 50
  number_answer_categories <- 4
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  beta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (number_answer_categories - 1))))
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  
  probabilities <- probabilities_and_likelihoods(theta, responses, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior = prior, output = "probs")
  
  expect_equal(round(probabilities[1,], 3), c(.052, .237, .461, .249))
  expect_equal(round(probabilities[40,], 3), c(.096, .343, .414, .148))
  expect_equal(dim(probabilities), c(50, 4))
}) 

context("SM model")

test_that("model is SM, 3 dimensions, 4 categories, estimator is EAP, deriv is true", {
  theta <- c(0, 0, 0)
  model <- "SM"
  administered <- 1:50
  number_dimensions <- 3
  estimator <- "EAP"
  responses <- c(rep(c(0, 1, 2), 16), 0, 1)
  prior <- matrix(c(1.2, 1.5, 1.7, 1.5, .9, 1.5, 1.7, 1.5, 1.1), ncol = 3)
  
  number_items <- 50
  number_answer_categories <- 4
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  beta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (number_answer_categories - 1))))
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  
  probabilities <- probabilities_and_likelihoods(theta, responses, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior = prior, output = "both")
  
  expect_equal(round(probabilities$probabilities[1,], 3), c(.052, .275, .505, .168))
  expect_equal(round(probabilities$probabilities[40,], 3), c(.096, .397, .433, .075))
  expect_equal(dim(probabilities$probabilities), c(50, 4))
  expect_equal(dim(probabilities$likelihoods), c(1, 1))
  expect_equal(dim(attr(probabilities$likelihoods, "gradient")), c(1, 3))
  expect_equal(dim(attr(probabilities$likelihoods, "hessian")), c(3, 3))
  expect_equal(length(probabilities), 2)
  
  expect_equal(round(probabilities$likelihoods[1], 3), -84.888)
  expect_equal(round(attr(probabilities$likelihoods, "gradient"), 3), matrix(c(-7.74, -7.28, -8.297), ncol = 3))
  expect_equal(round(attr(probabilities$likelihoods, "hessian")[1,], 3), c(-12.905, -12.184, -12.239))
}) 

test_that("model is SM, 3 dimensions, varying number of categories, estimator is EAP, deriv is true", {
  theta <- c(0, 0, 0)
  model <- "SM"
  administered <- 1:50
  number_dimensions <- 3
  estimator <- "EAP"
  responses <- c(rep(c(0, 1, 2), 16), 1, 1)
  prior <- matrix(c(1.2, 1.5, 1.7, 1.5, .9, 1.5, 1.7, 1.5, 1.1), ncol = 3)
  
  number_items <- 50
  max_number_answer_categories <- 5
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  eta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (max_number_answer_categories - 1))))
  eta[c(1, 5:10), 3:4] <- NA
  eta[c(40:45), 4] <- NA
  beta <- row_cumsum(eta)
  guessing <- NULL
  
  probabilities <- probabilities_and_likelihoods(theta, responses, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior = prior, output = "both")
  
  expect_equal(round(probabilities$probabilities[1,], 3), c(.052, .011, .937, NA, NA))
  expect_equal(round(probabilities$probabilities[40,], 3), c(.096, .037, .053, .815, NA))
  expect_equal(round(probabilities$probabilities[50,], 3), c(.039, 0.006, 0.004, 0.008, 0.943))
  expect_equal(dim(probabilities$probabilities), c(50, 5))
  expect_equal(dim(probabilities$likelihoods), c(1, 1))
  expect_equal(dim(attr(probabilities$likelihoods, "gradient")), c(1, 3))
  expect_equal(dim(attr(probabilities$likelihoods, "hessian")), c(3, 3))
  expect_equal(length(probabilities), 2)
  
  expect_equal(round(probabilities$likelihoods[1], 3), -133.928)
  expect_equal(round(attr(probabilities$likelihoods, "gradient"), 3), matrix(c(-22.653, -21.523, -21.73), ncol = 3))
  expect_equal(round(attr(probabilities$likelihoods, "hessian")[1,], 3), c(-9.577, -9.092, -10.001))
}) 

test_that("model is SM, 3 dimensions, 4 categories, estimator is ML, deriv is true", {
  theta <- c(0, 0, 0)
  model <- "SM"
  administered <- 1:50
  number_dimensions <- 3
  estimator <- "ML"
  responses <- c(rep(c(0, 1, 2), 16), 1, 1)
  prior <- NULL
  
  number_items <- 50
  number_answer_categories <- 4
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  beta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (number_answer_categories - 1))))
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  
  probabilities <- probabilities_and_likelihoods(theta, responses, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior = prior, output = "both")
  
  expect_equal(round(probabilities$probabilities[1,], 3), c(.052, .275, .505, .168))
  expect_equal(round(probabilities$probabilities[40,], 3), c(.096, .397, .433, .075))
  expect_equal(dim(probabilities$probabilities), c(50, 4))
  expect_equal(length(probabilities$likelihoods), 1)
  expect_equal(dim(attr(probabilities$likelihoods, "gradient")), c(1, 3))
  expect_equal(dim(attr(probabilities$likelihoods, "hessian")), c(3, 3))
  expect_equal(length(probabilities), 2)
}) 

test_that("model is SM, 3 dimensions, 4 categories, estimator is ML, deriv is false", {
  theta <- c(0, 0, 0)
  model <- "SM"
  administered <- 1:50
  number_dimensions <- 3
  estimator <- "ML"
  responses <- c(rep(c(0, 1, 2), 16), 1, 1)
  prior <- NULL
  
  number_items <- 50
  number_answer_categories <- 4
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  beta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (number_answer_categories - 1))))
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  
  probabilities <- probabilities_and_likelihoods(theta, responses, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior = prior, output = "probs")
    
  expect_equal(round(probabilities[1,], 3), c(.052, .275, .505, .168))
  expect_equal(round(probabilities[40,], 3), c(.096, .397, .433, .075))
  expect_equal(dim(probabilities), c(50, 4))
}) 
