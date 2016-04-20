# only for whithin R:
'
library(testthat)
'

make_random_seed_exist <- rnorm(1)

context("3PLM model")

test_that("model is 3PLM, 2 dimensions, 2 categories, estimator is maximum_likelihood", {
  theta <- c(.3, .5)
  model <- "3PLM"
  administered <- 1:50
  number_dimensions <- 2
  estimator <- "maximum_likelihood"
  answers <- rep(c(1, 0), 25)
  
  number_items <- 50
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  
  likelihood <- likelihood_or_post_density(theta, answers = answers, model, administered, number_dimensions, estimator, alpha, beta, guessing)
    
  expect_equal(length(likelihood), 1)
  expect_equal(dim(attr(likelihood, "gradient")), c(1, 2))
  expect_equal(dim(attr(likelihood, "hessian")), c(2, 2))
  
  expect_equal(round(likelihood[1], 3), -68.063)
  expect_equal(round(attr(likelihood, "gradient"), 3), matrix(c(-7.519, -10.47), ncol = 2))
  expect_equal(round(attr(likelihood, "hessian"), 3), matrix(c(-3.005, -2.816, -2.816, -3.373), ncol = 2))
})

test_that("model is 3PLM, 1 dimension, 2 categories, estimator is maximum_aposteriori", {
  theta <- .2
  model <- "3PLM"
  administered <- 1:50
  number_dimensions <- 1
  estimator <- "maximum_aposteriori"
  answers <- rep(c(1, 0), 25)
  prior_parameters <- list(mu = .5, Sigma = diag(1))
  
  number_items <- 50
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  
  post_density <- likelihood_or_post_density(theta, answers, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior_parameters = prior_parameters)
    
  expect_equal(dim(post_density), c(1, 1))
  expect_equal(dim(attr(post_density, "gradient")), c(1, 1))
  expect_equal(dim(attr(post_density, "hessian")), c(1, 1))
  
  expect_equal(round(post_density[1], 3), -48.018)
  expect_equal(round(attr(post_density, "gradient"), 3), matrix(-3.175))
  expect_equal(round(attr(post_density, "hessian"), 3), matrix(-6.35))
})

test_that("model is 3PLM, 3 dimensions, 2 categories, estimator is expected_aposteriori", {
  theta <- c(.1, -.4, 2)
  model <- "3PLM"
  administered <- 1:50
  number_dimensions <- 3
  estimator <- "expected_aposteriori"
  answers <- rep(c(1, 0), 25)
  prior_parameters <- list(mu = c(-1, .5, 0), Sigma = matrix(c(1.2, 1.5, 1.7, 1.5, .9, 1.5, 1.7, 1.5, 1.1), ncol = 3)) 
  
  number_items <- 50
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  
  post_density <- likelihood_or_post_density(theta, answers, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior_parameters = prior_parameters)
  
  expect_equal(dim(post_density), c(1, 1))
  expect_equal(dim(attr(post_density, "gradient")), c(1, 3))
  expect_equal(dim(attr(post_density, "hessian")), c(3, 3))
  
  expect_equal(round(post_density[1], 3), -90.823)
  expect_equal(round(attr(post_density, "gradient"), 3), matrix(c(-9.467, -16.565, -10.403), ncol = 3))
  expect_equal(round(attr(post_density, "hessian")[2,], 3), c(-3.018, -1.180, -2.903))
})

context("model is GPCM")

test_that("model is GPCM, 1 dimensions, 2 categories, estimator is maximum_aposteriori", {
  theta <- .6
  model <- "GPCM"
  administered <- 1:50
  number_dimensions <- 1
  estimator <- "maximum_aposteriori"
  answers <- rep(c(0, 1), 25)
  prior_parameters <- list(mu = -1.5, Sigma  = diag(1) * .6)
  
  number_items <- 50
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  
  post_density <- likelihood_or_post_density(theta, answers, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior_parameters = prior_parameters)

  expect_equal(dim(post_density), c(1, 1))
  expect_equal(dim(attr(post_density, "gradient")), c(1, 1))
  expect_equal(dim(attr(post_density, "hessian")), c(1, 1))
  
  expect_equal(round(post_density[1], 3), -40.266)
  expect_equal(round(attr(post_density, "gradient"), 3), matrix(-11.011))
  expect_equal(round(attr(post_density, "hessian"), 3), matrix(-10.722))
}) 

test_that("model is GPCM, 3 dimensions, 4 categories, estimator is expected_aposteriori", {
  theta <- c(-2, 1, .4)
  model <- "GPCM"
  administered <- 1:50
  number_dimensions <- 3
  estimator <- "expected_aposteriori"
  answers <- c(rep(c(0, 1, 2), 16), 1, 1)
  prior_parameters <-  list(mu = c(.1, -.3, 2), Sigma = matrix(c(1.2, 1.5, 1.7, 1.5, .9, 1.5, 1.7, 1.5, 1.1), ncol = 3))
  
  number_items <- 50
  number_answer_categories <- 4
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  eta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (number_answer_categories - 1))))
  beta <- t(apply(eta, 1, cumsum))
  guessing <- NULL
  
  post_density <- likelihood_or_post_density(theta, answers, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior_parameters = prior_parameters)
  
  expect_equal(dim(post_density), c(1, 1))
  expect_equal(dim(attr(post_density, "gradient")), c(1, 3))
  expect_equal(dim(attr(post_density, "hessian")), c(3, 3))
  
  expect_equal(round(post_density[1], 3), -93.505)
  expect_equal(round(attr(post_density, "gradient"), 3), matrix(c(-6.379, -7.14, -13.055), ncol = 3))
  expect_equal(round(attr(post_density, "hessian")[1,], 3), c(-18.092, -17.893, -17.446))
}) 

test_that("model is GPCM, 3 dimensions, varying number of categories, estimator is expected_aposteriori", {
  theta <- c(.7, -1.2, 2)
  model <- "GPCM"
  administered <- 1:50
  number_dimensions <- 3
  estimator <- "expected_aposteriori"
  answers <- c(rep(0,5), rep(1,5), rep(c(1:3),9), rep(c(0:2), 4), 3)
  prior_parameters <-  list(mu = c(1, 1, 1), Sigma = matrix(c(1.2, 1.5, 1.7, 1.5, .9, 1.5, 1.7, 1.5, 1.1), ncol = 3))
  
  number_items <- 50
  max_number_answer_categories <- 5
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  eta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (max_number_answer_categories - 1))))
  eta[c(1, 5:10), 3:4] <- NA
  eta[c(40:45), 4] <- NA
  beta <- row_cumsum(eta)
  guessing <- NULL
  
  post_density <- likelihood_or_post_density(theta, answers, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior_parameters = prior_parameters)
  
  expect_equal(dim(post_density), c(1, 1))
  expect_equal(dim(attr(post_density, "gradient")), c(1, 3))
  expect_equal(dim(attr(post_density, "hessian")), c(3, 3))
  
  expect_equal(round(post_density[1], 3), -153.806)
  expect_equal(round(attr(post_density, "gradient"), 3), matrix(c(-56.46, -53.098, -57.615), ncol = 3))
  expect_equal(round(attr(post_density, "hessian")[1,], 3), c(-20.055, -18.670, -19.012))
}) 


context("GRM model")

test_that("model is GRM, 1 dimension, 2 categories, estimator is expected_aposteriori", {
  theta <- .6
  model <- "GRM"
  administered <- 1:50
  number_dimensions <- 1
  estimator <- "expected_aposteriori"
  answers <- rep(c(1, 0), 25)
  prior_parameters <- list(mu = 0, Sigma = matrix(1.2))
  
  number_items <- 50
  number_answer_categories <- 2
  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  
  post_density <- likelihood_or_post_density(theta, answers, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior_parameters = prior_parameters)
  
  expect_equal(dim(post_density), c(1, 1))
  expect_equal(dim(attr(post_density, "gradient")), c(1, 1))
  expect_equal(dim(attr(post_density, "hessian")), c(1, 1))
  
  expect_equal(round(post_density[1], 3), -49.666)
  expect_equal(round(attr(post_density, "gradient"), 3), matrix(-3.256))
  expect_equal(round(attr(post_density, "hessian"), 3), matrix(-9.889))
}) 

test_that("model is GRM, 3 dimensions, 4 categories, estimator is expected_aposteriori", {
  theta <- c(.2, 2, 2.5)
  model <- "GRM"
  administered <- 1:50
  number_dimensions <- 3
  estimator <- "expected_aposteriori"
  answers <- c(rep(c(0, 1, 2), 16), 0, 1)
  prior_parameters <-  list(mu = c(0, 0, 0), Sigma = matrix(c(1.2, 1.5, 1.7, 1.5, .9, 1.5, 1.7, 1.5, 1.1), ncol = 3))
  
  number_items <- 50
  number_answer_categories <- 4
  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  beta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (number_answer_categories - 1))))
  guessing <- NULL
  
  post_density <- likelihood_or_post_density(theta, answers, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior_parameters = prior_parameters)

  expect_equal(dim(post_density), c(1, 1))
  expect_equal(dim(attr(post_density, "gradient")), c(1, 3))
  expect_equal(dim(attr(post_density, "hessian")), c(3, 3))
  
  expect_equal(round(post_density[1], 3), -218.875)
  expect_equal(round(attr(post_density, "gradient"), 3), matrix(c(-43.099, -38.046, -40.282), ncol = 3))
  expect_equal(round(attr(post_density, "hessian")[1,], 3), c(-1.864, -3.085, -2.745))
}) 

 
test_that("model is GRM, 3 dimensions, varying number of categories, estimator is expected_aposteriori", {
  theta <- c(-2, -1.1, 2)
  model <- "GRM"
  administered <- 1:50
  number_dimensions <- 3
  estimator <- "expected_aposteriori"
  answers <- c(rep(0,5), rep(1,5), rep(c(1:3),9), rep(c(0:2), 4), 3)
  prior_parameters <-  list(mu = c(0, 0, 0), Sigma = matrix(c(1.2, 1.5, 1.7, 1.5, .9, 1.5, 1.7, 1.5, 1.1), ncol = 3))
  
  number_items <- 50
  max_number_answer_categories <- 5
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  eta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (max_number_answer_categories - 1))))
  eta[c(1, 5:10), 3:4] <- NA
  eta[c(40:45), 4] <- NA
  beta <- row_cumsum(eta)
  guessing <- NULL
  
  post_density <- likelihood_or_post_density(theta, answers, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior_parameters = prior_parameters)
  
  expect_equal(dim(post_density), c(1, 1))
  expect_equal(dim(attr(post_density, "gradient")), c(1, 3))
  expect_equal(dim(attr(post_density, "hessian")), c(3, 3))
  
  expect_equal(round(post_density[1], 3), -416.437)
  expect_equal(round(attr(post_density, "gradient"), 3), matrix(c(-14.884, -15.49, -11.676), ncol = 3))
  expect_equal(round(attr(post_density, "hessian")[1,], 3), c(-11.650, -11.245, -11.237))
}) 

context("SM model")

test_that("model is SM, 1 dimension, 2 categories, estimator is expected_aposteriori", {
  theta <- .6
  model <- "SM"
  administered <- c(4, 19, 38)
  number_dimensions <- 1
  estimator <- "expected_aposteriori"
  answers <- rep(c(1, 0), 25)
  prior_parameters <- list(mu = 0, Sigma = matrix(1.2))
  
  number_items <- 50
  number_answer_categories <- 2
  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  
  post_density <- likelihood_or_post_density(theta, answers, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior_parameters = prior_parameters)
  
  expect_equal(dim(post_density), c(1, 1))
  expect_equal(dim(attr(post_density, "gradient")), c(1, 1))
  expect_equal(dim(attr(post_density, "hessian")), c(1, 1))
  
  expect_equal(round(post_density[1], 3), -1.077)
  expect_equal(round(attr(post_density, "gradient"), 3), matrix(-0.577))
  expect_equal(round(attr(post_density, "hessian"), 3), matrix(-1.105))
}) 

test_that("model is SM, 3 dimensions, 4 categories, estimator is expected_aposteriori", {
  theta <- c(2, 1.2, -2.3)
  model <- "SM"
  administered <- 1:50
  number_dimensions <- 3
  estimator <- "expected_aposteriori"
  answers <- c(rep(c(0, 1, 2), 16), 0, 1)
  prior_parameters <-  list(mu = c(0, 0, 0), Sigma = matrix(c(1.2, 1.5, 1.7, 1.5, .9, 1.5, 1.7, 1.5, 1.1), ncol = 3))
  
  number_items <- 50
  number_answer_categories <- 4
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  beta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (number_answer_categories - 1))))
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  
  post_density <- likelihood_or_post_density(theta, answers, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior_parameters = prior_parameters)

  expect_equal(dim(post_density), c(1, 1))
  expect_equal(dim(attr(post_density, "gradient")), c(1, 3))
  expect_equal(dim(attr(post_density, "hessian")), c(3, 3))
  
  expect_equal(round(post_density[1], 3), -93.544)
  expect_equal(round(attr(post_density, "gradient"), 3), matrix(c(-16.913, -14.88, -18.074), ncol = 3))
  expect_equal(round(attr(post_density, "hessian")[1,], 3), c(-9.229, -9.184, -10.147))
}) 

test_that("model is SM, 3 dimensions, varying number of categories, estimator is expected_aposteriori, deriv is true", {
  theta <- c(1.8, -1.4, 2.5)
  model <- "SM"
  administered <- 1:50
  number_dimensions <- 3
  estimator <- "expected_aposteriori"
  answers <- c(rep(0,5), rep(1,5), rep(c(1:3),9), rep(c(0:2), 4), 3)
  prior_parameters <-  list(mu = c(0, 0, 0), Sigma = matrix(c(1.2, 1.5, 1.7, 1.5, .9, 1.5, 1.7, 1.5, 1.1), ncol = 3))
  
  number_items <- 50
  max_number_answer_categories <- 5
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  eta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (max_number_answer_categories - 1))))
  eta[c(1, 5:10), 3:4] <- NA
  eta[c(40:45), 4] <- NA
  beta <- row_cumsum(eta)
  guessing <- NULL
  
  post_density <- likelihood_or_post_density(theta, answers, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior_parameters = prior_parameters)
  
  expect_equal(dim(post_density), c(1, 1))
  expect_equal(dim(attr(post_density, "gradient")), c(1, 3))
  expect_equal(dim(attr(post_density, "hessian")), c(3, 3))
  
  expect_equal(round(post_density[1], 3), -226.666)
  expect_equal(round(attr(post_density, "gradient"), 3), matrix(c(-37.975, -42.21, -36.426), ncol = 3))
  expect_equal(round(attr(post_density, "hessian")[1,], 3), c(-1.872, -3.647, -3.411))
}) 


