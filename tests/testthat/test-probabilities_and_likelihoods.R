# only for whithin R:
'
library(testthat)
'

make_random_seed_exist <- rnorm(1)

context("3PLM model")

test_that("model is 3PLM, 1 dimension, 2 categories, estimator is maximum_aposteriori, deriv is false", {
  theta <- .5
  model <- "3PLM"
  administered <- 1:50
  number_dimensions <- 1
  estimator <- "maximum_aposteriori"
 
  number_items <- 50
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  
  probabilities <- probabilities_and_likelihood(theta, responses = NULL, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior = NULL, output = "probs")
  
  expect_equal(round(probabilities[1,], 3), c(.293, .707))
  expect_equal(round(probabilities[40,], 3), c(.329, .671))
  expect_equal(dim(probabilities), c(50, 2))
})

test_that("model is 3PLM, 2 dimensions, 2 categories, estimator is maximum_likelihood, deriv is true", {
  theta <- c(.3, .5)
  model <- "3PLM"
  administered <- 1:50
  number_dimensions <- 2
  estimator <- "maximum_likelihood"
  responses <- rep(c(1, 0), 25)
  
  number_items <- 50
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  
  probabilities <- probabilities_and_likelihood(theta, responses = responses, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior = NULL, output = "both")
    
  expect_equal(round(probabilities$probabilities[1,], 3), c(.232, .768))
  expect_equal(round(probabilities$probabilities[40,], 3), c(.25, .75))
  expect_equal(dim(probabilities$probabilities), c(50, 2))
  expect_equal(length(probabilities$likelihood), 1)
  expect_equal(dim(attr(probabilities$likelihood, "gradient")), c(1, 2))
  expect_equal(dim(attr(probabilities$likelihood, "hessian")), c(2, 2))
  expect_equal(length(probabilities), 2)
  
  expect_equal(round(probabilities$likelihood[1], 3), -68.063)
  expect_equal(round(attr(probabilities$likelihood, "gradient"), 3), matrix(c(-7.519, -10.47), ncol = 2))
  expect_equal(round(attr(probabilities$likelihood, "hessian"), 3), matrix(c(-3.005, -2.816, -2.816, -3.373), ncol = 2))
})

test_that("model is 3PLM, 1 dimension, 2 categories, estimator is maximum_aposteriori, deriv is true", {
  theta <- .2
  model <- "3PLM"
  administered <- 1:50
  number_dimensions <- 1
  estimator <- "maximum_aposteriori"
  responses <- rep(c(1, 0), 25)
  prior <- diag(1)
  
  number_items <- 50
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  
  probabilities <- probabilities_and_likelihood(theta, responses, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior = prior, output = "both")
    
  expect_equal(round(probabilities$probabilities[1,], 3), c(.325, .675))
  expect_equal(round(probabilities$probabilities[40,], 3), c(.357, .643))
  expect_equal(dim(probabilities$probabilities), c(50, 2))
  expect_equal(dim(probabilities$likelihood), c(1, 1))
  expect_equal(dim(attr(probabilities$likelihood, "gradient")), c(1, 1))
  expect_equal(dim(attr(probabilities$likelihood, "hessian")), c(1, 1))
  expect_equal(length(probabilities), 2)
  
  expect_equal(round(probabilities$likelihood[1], 3), -47.993)
  expect_equal(round(attr(probabilities$likelihood, "gradient"), 3), matrix(-3.675))
  expect_equal(round(attr(probabilities$likelihood, "hessian"), 3), matrix(-6.35))
})

test_that("model is 3PLM, 3 dimensions, 2 categories, estimator is expected_aposteriori, deriv is true", {
  theta <- c(.1, -.4, 2)
  model <- "3PLM"
  administered <- 1:50
  number_dimensions <- 3
  estimator <- "expected_aposteriori"
  responses <- rep(c(1, 0), 25)
  prior <- matrix(c(1.2, 1.5, 1.7, 1.5, .9, 1.5, 1.7, 1.5, 1.1), ncol = 3)
  
  number_items <- 50
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  
  probabilities <- probabilities_and_likelihood(theta, responses, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior = prior, output = "both")
  
  expect_equal(round(probabilities$probabilities[1,], 3), c(.087, .913))
  expect_equal(round(probabilities$probabilities[40,], 3), c(.113, .887))
  expect_equal(dim(probabilities$probabilities), c(50, 2))
  expect_equal(dim(probabilities$likelihood), c(1, 1))
  expect_equal(dim(attr(probabilities$likelihood, "gradient")), c(1, 3))
  expect_equal(dim(attr(probabilities$likelihood, "hessian")), c(3, 3))
  expect_equal(length(probabilities), 2)
})

context("model is GPCM")

test_that("model is GPCM, 1 dimensions, 2 categories, estimator is maximum_aposteriori, deriv is true", {
  theta <- .6
  model <- "GPCM"
  administered <- 1:50
  number_dimensions <- 1
  estimator <- "maximum_aposteriori"
  responses <- rep(c(0, 1), 25)
  prior <- diag(1) * .6
  
  number_items <- 50
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  
  probabilities <- probabilities_and_likelihood(theta, responses, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior = prior, output = "both")
    
  expect_equal(round(probabilities$probabilities[1,], 3), c(.23, .77))
  expect_equal(round(probabilities$probabilities[40,], 3), c(.369, .631))
  expect_equal(dim(probabilities$probabilities), c(50, 2))
  expect_equal(dim(probabilities$likelihood), c(1, 1))
  expect_equal(dim(attr(probabilities$likelihood, "gradient")), c(1, 1))
  expect_equal(dim(attr(probabilities$likelihood, "hessian")), c(1, 1))
  expect_equal(length(probabilities), 2)
  
  expect_equal(round(probabilities$likelihood[1], 3), -36.891)
  expect_equal(round(attr(probabilities$likelihood, "gradient"), 3), matrix(-8.511))
  expect_equal(round(attr(probabilities$likelihood, "hessian"), 3), matrix(-10.722))
}) 

test_that("model is GPCM, 3 dimensions, 4 categories, estimator is expected_aposteriori, deriv is true", {
  theta <- c(-2, 1, .4)
  model <- "GPCM"
  administered <- 1:50
  number_dimensions <- 3
  estimator <- "expected_aposteriori"
  responses <- c(rep(c(0, 1, 2), 16), 1, 1)
  prior <-  matrix(c(1.2, 1.5, 1.7, 1.5, .9, 1.5, 1.7, 1.5, 1.1), ncol = 3)
  
  number_items <- 50
  number_answer_categories <- 4
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  eta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (number_answer_categories - 1))))
  beta <- t(apply(eta, 1, cumsum))
  guessing <- NULL
  
  probabilities <- probabilities_and_likelihood(theta, responses, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior = prior, output = "both")
  
  expect_equal(round(probabilities$probabilities[1,], 3), c(.033, .352, .514, .102))
  expect_equal(round(probabilities$probabilities[40,], 3), c(.037, .370, .501, .092))
  expect_equal(dim(probabilities$probabilities), c(50, 4))
  expect_equal(dim(probabilities$likelihood), c(1, 1))
  expect_equal(dim(attr(probabilities$likelihood, "gradient")), c(1, 3))
  expect_equal(dim(attr(probabilities$likelihood, "hessian")), c(3, 3))
  expect_equal(length(probabilities), 2)
  
  expect_equal(round(probabilities$likelihood[1], 3), -95.13)
  expect_equal(round(attr(probabilities$likelihood, "gradient"), 3), matrix(c(-7.362, -9.081, -10.707), ncol = 3))
  expect_equal(round(attr(probabilities$likelihood, "hessian")[1,], 3), c(-18.092, -17.893, -17.446))
}) 

test_that("model is GPCM, 3 dimensions, varying number of categories, estimator is expected_aposteriori, deriv is true", {
  theta <- c(.7, -1.2, 2)
  model <- "GPCM"
  administered <- 1:50
  number_dimensions <- 3
  estimator <- "expected_aposteriori"
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
  
  probabilities <- probabilities_and_likelihood(theta, responses, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior = prior, output = "both")
  
  expect_equal(round(probabilities$probabilities[1,], 3), c(.001, .066, .932, NA,  NA))
  expect_equal(round(probabilities$probabilities[40,], 3), c(.002, .047, .332, .620,  NA))
  expect_equal(round(probabilities$probabilities[50,], 3), c(.000, .003, .066, .374, .556))
  expect_equal(dim(probabilities$probabilities), c(50, 5))
  expect_equal(dim(probabilities$likelihood), c(1, 1))
  expect_equal(dim(attr(probabilities$likelihood, "gradient")), c(1, 3))
  expect_equal(dim(attr(probabilities$likelihood, "hessian")), c(3, 3))
  expect_equal(length(probabilities), 2)
  
  expect_equal(round(probabilities$likelihood[1], 3), -220.168)
  expect_equal(round(attr(probabilities$likelihood, "gradient"), 3), matrix(c(-76.815, -72.813, -80.579), ncol = 3))
  expect_equal(round(attr(probabilities$likelihood, "hessian")[1,], 3), c(-20.055, -18.670, -19.012))
}) 


test_that("model is GPCM, 3 dimensions, 4 categories, estimator is maximum_likelihood, deriv is true", {
  theta <- c(-1, 2, 3)
  model <- "GPCM"
  administered <- 1:50
  number_dimensions <- 3
  estimator <- "maximum_likelihood"
  responses <- c(rep(c(0, 1, 2), 16), 1, 1)
  prior <- NULL
  
  number_items <- 50
  number_answer_categories <- 4
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  eta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (number_answer_categories - 1))))
  beta <- t(apply(eta, 1, cumsum))
  guessing <- NULL
  
  probabilities <- probabilities_and_likelihood(theta, responses, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior = prior, output = "both")
  
  expect_equal(round(probabilities$probabilities[1,], 3), c(.000, .025, .342, .633))
  expect_equal(round(probabilities$probabilities[40,], 3), c(.000, .005, .178, .816))
  expect_equal(dim(probabilities$probabilities), c(50, 4))
  expect_equal(length(probabilities$likelihood), 1)
  expect_equal(dim(attr(probabilities$likelihood, "gradient")), c(1, 3))
  expect_equal(dim(attr(probabilities$likelihood, "hessian")), c(3, 3))
  expect_equal(length(probabilities), 2)
}) 

test_that("model is GPCM, 3 dimensions, 4 categories, estimator is expected_aposteriori, deriv is false", {
  theta <- c(-2, -1.3, -2.1)
  model <- "GPCM"
  administered <- 1:50
  number_dimensions <- 3
  estimator <- "expected_aposteriori"
  responses <- c(rep(c(0, 1, 2), 16), 1, 1)
  prior <- NULL
  
  number_items <- 50
  number_answer_categories <- 4
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  eta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (number_answer_categories - 1))))
  beta <- t(apply(eta, 1, cumsum))
  guessing <- NULL
  
  probabilities <- probabilities_and_likelihood(theta, responses, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior = prior, output = "probs")
  
  expect_equal(round(probabilities[1,], 3), c(.380, .521, .097, .002))
  expect_equal(round(probabilities[40,], 3), c(.775, .217, .008, .000))
  expect_equal(dim(probabilities), c(50, 4))
}) 

context("GRM model")

test_that("model is GRM, 3 dimensions, 4 categories, estimator is expected_aposteriori, deriv is true", {
  theta <- c(.2, 2, 2.5)
  model <- "GRM"
  administered <- 1:50
  number_dimensions <- 3
  estimator <- "expected_aposteriori"
  responses <- c(rep(c(0, 1, 2), 16), 0, 1)
  prior <- matrix(c(1.2, 1.5, 1.7, 1.5, .9, 1.5, 1.7, 1.5, 1.1), ncol = 3)
  
  number_items <- 50
  number_answer_categories <- 4
  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  beta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (number_answer_categories - 1))))
  guessing <- NULL
  
  probabilities <- probabilities_and_likelihood(theta, responses, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior = prior, output = "both")
  
  expect_equal(round(probabilities$probabilities[1,], 3), c(.007, .042, .226, .725))
  expect_equal(round(probabilities$probabilities[40,], 3), c(.003, .021, .130, .846))
  expect_equal(dim(probabilities$probabilities), c(50, 4))
  expect_equal(dim(probabilities$likelihood), c(1, 1))
  expect_equal(dim(attr(probabilities$likelihood, "gradient")), c(1, 3))
  expect_equal(dim(attr(probabilities$likelihood, "hessian")), c(3, 3))
  expect_equal(length(probabilities), 2)
  
  expect_equal(round(probabilities$likelihood[1], 3), -218.875)
  expect_equal(round(attr(probabilities$likelihood, "gradient"), 3), matrix(c(-43.099, -38.046, -40.282), ncol = 3))
  expect_equal(round(attr(probabilities$likelihood, "hessian")[1,], 3), c(-1.864, -3.085, -2.745))
}) 

test_that("model is GRM, 3 dimensions, 4 categories, estimator is maximum_likelihood, deriv is true", {
  theta <- c(2, 1, 3)
  model <- "GRM"
  administered <- 1:50
  number_dimensions <- 3
  estimator <- "maximum_likelihood"
  responses <- c(rep(c(0, 1, 2), 16), 1, 1)
  prior <- NULL
  
  number_items <- 50
  number_answer_categories <- 4
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  beta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (number_answer_categories - 1))))
  guessing <- NULL
  
  probabilities <- probabilities_and_likelihood(theta, responses, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior = prior, output = "both")
  
  expect_equal(round(probabilities$probabilities[1,], 3), c(.003, .018, .113, .867))
  expect_equal(round(probabilities$probabilities[40,], 3), c(.002, .012, .081, .905))
  expect_equal(dim(probabilities$probabilities), c(50, 4))
  expect_equal(length(probabilities$likelihood), 1)
  expect_equal(dim(attr(probabilities$likelihood, "gradient")), c(1, 3))
  expect_equal(dim(attr(probabilities$likelihood, "hessian")), c(3, 3))
  expect_equal(length(probabilities), 2)
}) 

test_that("model is GRM, 3 dimensions, varying number of categories, estimator is expected_aposteriori, deriv is true", {
  theta <- c(-2, -1.1, 2)
  model <- "GRM"
  administered <- 1:50
  number_dimensions <- 3
  estimator <- "expected_aposteriori"
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
  
  probabilities <- probabilities_and_likelihood(theta, responses, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior = prior, output = "both")
  
  expect_equal(round(probabilities$probabilities[1,], 3), c(.069, .000, .985, NA, NA))
  expect_equal(round(probabilities$probabilities[40,], 3), c(.113, .000, .023, .928, NA))
  expect_equal(round(probabilities$probabilities[50,], 3), c(.261, .000, .000, .036, .933))
  expect_equal(dim(probabilities$probabilities), c(50, 5))
  expect_equal(dim(probabilities$likelihood), c(1, 1))
  expect_equal(dim(attr(probabilities$likelihood, "gradient")), c(1, 3))
  expect_equal(dim(attr(probabilities$likelihood, "hessian")), c(3, 3))
  expect_equal(length(probabilities), 2)
  
  expect_equal(round(probabilities$likelihood[1], 3), -410.588)
  expect_equal(round(attr(probabilities$likelihood, "gradient"), 3), matrix(c(-15.236, -15.309, -13.404), ncol = 3))
  expect_equal(round(attr(probabilities$likelihood, "hessian")[1,], 3), c(-11.989, -11.077, -11.042))
}) 

test_that("model is GRM, 3 dimensions, 4 categories, estimator is expected_aposteriori, deriv is false", {
  theta <- c(.1, 2.2, 1.8)
  model <- "GRM"
  administered <- 1:50
  number_dimensions <- 3
  estimator <- "expected_aposteriori"
  responses <- c(rep(c(0, 1, 2), 16), 1, 1)
  prior <- NULL
  
  number_items <- 50
  number_answer_categories <- 4
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  beta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (number_answer_categories - 1))))
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  
  probabilities <- probabilities_and_likelihood(theta, responses, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior = prior, output = "probs")
  
  expect_equal(round(probabilities[1,], 3), c(.010, .059, .285, .646))
  expect_equal(round(probabilities[40,], 3), c(.005, .032, .185, .778))
  expect_equal(dim(probabilities), c(50, 4))
}) 

context("SM model")

test_that("model is SM, 3 dimensions, 4 categories, estimator is expected_aposteriori, deriv is true", {
  theta <- c(2, 1.2, -2.3)
  model <- "SM"
  administered <- 1:50
  number_dimensions <- 3
  estimator <- "expected_aposteriori"
  responses <- c(rep(c(0, 1, 2), 16), 0, 1)
  prior <- matrix(c(1.2, 1.5, 1.7, 1.5, .9, 1.5, 1.7, 1.5, 1.1), ncol = 3)
  
  number_items <- 50
  number_answer_categories <- 4
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  beta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (number_answer_categories - 1))))
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  
  probabilities <- probabilities_and_likelihood(theta, responses, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior = prior, output = "both")
  
  expect_equal(round(probabilities$probabilities[1,], 3), c(.045, .245, .510, .201))
  expect_equal(round(probabilities$probabilities[40,], 3), c(.094, .393, .436, .077))
  expect_equal(dim(probabilities$probabilities), c(50, 4))
  expect_equal(dim(probabilities$likelihood), c(1, 1))
  expect_equal(dim(attr(probabilities$likelihood, "gradient")), c(1, 3))
  expect_equal(dim(attr(probabilities$likelihood, "hessian")), c(3, 3))
  expect_equal(length(probabilities), 2)
  
  expect_equal(round(probabilities$likelihood[1], 3), -93.544)
  expect_equal(round(attr(probabilities$likelihood, "gradient"), 3), matrix(c(-16.913, -14.88, -18.074), ncol = 3))
  expect_equal(round(attr(probabilities$likelihood, "hessian")[1,], 3), c(-9.229, -9.184, -10.147))
}) 

test_that("model is SM, 3 dimensions, varying number of categories, estimator is expected_aposteriori, deriv is true", {
  theta <- c(1.8, -1.4, 2.5)
  model <- "SM"
  administered <- 1:50
  number_dimensions <- 3
  estimator <- "expected_aposteriori"
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
  
  probabilities <- probabilities_and_likelihood(theta, responses, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior = prior, output = "both")
  
  expect_equal(round(probabilities$probabilities[1,], 3), c(.009, .002, .990, NA, NA))
  expect_equal(round(probabilities$probabilities[40,], 3), c(.017, .007, .010, .966, NA))
  expect_equal(round(probabilities$probabilities[50,], 3), c(.003, .000, .000, .001, .996))
  expect_equal(dim(probabilities$probabilities), c(50, 5))
  expect_equal(dim(probabilities$likelihood), c(1, 1))
  expect_equal(dim(attr(probabilities$likelihood, "gradient")), c(1, 3))
  expect_equal(dim(attr(probabilities$likelihood, "hessian")), c(3, 3))
  expect_equal(length(probabilities), 2)
  
  expect_equal(round(probabilities$likelihood[1], 3), -223.143)
  expect_equal(round(attr(probabilities$likelihood, "gradient"), 3), matrix(c(-36.651, -40.743, -36.271), ncol = 3))
  expect_equal(round(attr(probabilities$likelihood, "hessian")[1,], 3), c(-1.398, -3.486, -2.896))
}) 

test_that("model is SM, 3 dimensions, 4 categories, estimator is maximum_likelihood, deriv is true", {
  theta <- c(-2.4, 1, 2.7)
  model <- "SM"
  administered <- 1:50
  number_dimensions <- 3
  estimator <- "maximum_likelihood"
  responses <- c(rep(c(0, 1, 2), 16), 1, 1)
  prior <- NULL
  
  number_items <- 50
  number_answer_categories <- 4
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  beta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (number_answer_categories - 1))))
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  
  probabilities <- probabilities_and_likelihood(theta, responses, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior = prior, output = "both")
  
  expect_equal(round(probabilities$probabilities[1,], 3), c(.032, .189, .500, .278))
  expect_equal(round(probabilities$probabilities[40,], 3), c(.020, .127, .447, .407))
  expect_equal(dim(probabilities$probabilities), c(50, 4))
  expect_equal(length(probabilities$likelihood), 1)
  expect_equal(dim(attr(probabilities$likelihood, "gradient")), c(1, 3))
  expect_equal(dim(attr(probabilities$likelihood, "hessian")), c(3, 3))
  expect_equal(length(probabilities), 2)
}) 

test_that("model is SM, 3 dimensions, 4 categories, estimator is maximum_likelihood, deriv is false", {
  theta <- c(.2, -.4, 1.1)
  model <- "SM"
  administered <- 1:50
  number_dimensions <- 3
  estimator <- "maximum_likelihood"
  responses <- c(rep(c(0, 1, 2), 16), 1, 1)
  prior <- NULL
  
  number_items <- 50
  number_answer_categories <- 4
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  beta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (number_answer_categories - 1))))
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  
  probabilities <- probabilities_and_likelihood(theta, responses, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior = prior, output = "probs")
    
  expect_equal(round(probabilities[1,], 3), c(.030, .181, .496, .293))
  expect_equal(round(probabilities[40,], 3), c(.051, .271, .506, .171))
  expect_equal(dim(probabilities), c(50, 4))
}) 
