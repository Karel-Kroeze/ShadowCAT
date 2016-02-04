# only for whithin R:
'
library(testthat)
'

make_random_seed_exist <- rnorm(1)

context("3PLM model")

test_that("model is 3PLM, 1 dimension, 2 categories", {
  theta <- .5
  model <- "3PLM"
  administered <- 1:50
  number_dimensions <- 1 
  number_items <- 50
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  
  probabilities <- get_probabilities(theta, model, administered, alpha, beta, guessing)
  
  expect_equal(round(probabilities[1,], 3), c(.293, .707))
  expect_equal(round(probabilities[40,], 3), c(.329, .671))
  expect_equal(dim(probabilities), c(50, 2))
})

test_that("model is 3PLM, 2 dimensions, 2 categories", {
  theta <- c(.3, .5)
  model <- "3PLM"
  administered <- 1:50
  number_dimensions <- 2  
  number_items <- 50
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  
  probabilities <- get_probabilities(theta, model, administered, alpha, beta, guessing)
  
  expect_equal(round(probabilities[1,], 3), c(.232, .768))
  expect_equal(round(probabilities[40,], 3), c(.25, .75))
  expect_equal(dim(probabilities), c(50, 2))
})

test_that("not all items administered", {
  theta <- c(.3, .5)
  model <- "3PLM"
  administered <- c(1, 4, 19)
  number_dimensions <- 2  
  number_items <- 50
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  
  probabilities <- get_probabilities(theta, model, administered, alpha, beta, guessing)
  
  expect_equal(round(probabilities[1,], 3), c(.232, .768))
  expect_equal(round(probabilities[3,], 3), c(.637, .363))
  expect_equal(dim(probabilities), c(3, 2))
})

context("model is GPCM")

test_that("model is GPCM, 1 dimensions, 2 categories", {
  theta <- .6
  model <- "GPCM"
  administered <- 1:50
  number_dimensions <- 1
  number_items <- 50
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  
  probabilities <- get_probabilities(theta, model, administered, alpha, beta, guessing)
  
  expect_equal(round(probabilities[1,], 3), c(.23, .77))
  expect_equal(round(probabilities[40,], 3), c(.369, .631))
  expect_equal(dim(probabilities), c(50, 2))
}) 

test_that("model is GPCM, 3 dimensions, 4 categories", {
  theta <- c(-2, 1, .4)
  model <- "GPCM"
  administered <- 1:50
  number_dimensions <- 3  
  number_items <- 50
  number_answer_categories <- 4
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  eta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (number_answer_categories - 1))))
  beta <- t(apply(eta, 1, cumsum))
  guessing <- NULL
  
  probabilities <- get_probabilities(theta, model, administered, alpha, beta, guessing)
  
  expect_equal(round(probabilities[1,], 3), c(.033, .352, .514, .102))
  expect_equal(round(probabilities[40,], 3), c(.037, .370, .501, .092))
  expect_equal(dim(probabilities), c(50, 4))
}) 


test_that("model is GPCM, 3 dimensions, varying number of categories", {
  theta <- c(.7, -1.2, 2)
  model <- "GPCM"
  administered <- 1:50
  number_dimensions <- 3  
  number_items <- 50
  max_number_answer_categories <- 5
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  eta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (max_number_answer_categories - 1))))
  eta[c(1, 5:10), 3:4] <- NA
  eta[c(40:45), 4] <- NA
  beta <- row_cumsum(eta)
  guessing <- NULL
  
  probabilities <- get_probabilities(theta, model, administered, alpha, beta, guessing)
  
  expect_equal(round(probabilities[1,], 3), c(.001, .066, .932, NA,  NA))
  expect_equal(round(probabilities[40,], 3), c(.002, .047, .332, .620,  NA))
  expect_equal(round(probabilities[50,], 3), c(.000, .003, .066, .374, .556))
  expect_equal(dim(probabilities), c(50, 5))
})

context("GRM model")

test_that("model is GRM, 1 dimensions, 2 categories, estimator is expected_aposteriori, deriv is true", {
  theta <- .5
  model <- "GRM"
  administered <- 1:50
  number_dimensions <- 1 
  number_items <- 50
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  
  probabilities <- get_probabilities(theta, model, administered, alpha, beta, guessing)
  
  expect_equal(round(probabilities[1,], 3), c(.239, .761))
  expect_equal(round(probabilities[40,], 3), c(.381, .619))
  expect_equal(dim(probabilities), c(50, 2))
}) 

test_that("model is GRM, 3 dimensions, 4 categories, estimator is expected_aposteriori, deriv is true", {
  theta <- c(.2, 2, 2.5)
  model <- "GRM"
  administered <- 1:50
  number_dimensions <- 3
  number_items <- 50
  number_answer_categories <- 4
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  beta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (number_answer_categories - 1))))
  guessing <- NULL
  
  probabilities <- get_probabilities(theta, model, administered, alpha, beta, guessing)
  
  expect_equal(round(probabilities[1,], 3), c(.007, .042, .226, .725))
  expect_equal(round(probabilities[40,], 3), c(.003, .021, .130, .846))
  expect_equal(dim(probabilities), c(50, 4))
}) 

test_that("model is GRM, 3 dimensions, varying number of categories, estimator is expected_aposteriori, deriv is true", {
  theta <- c(-2, -1.1, 2)
  model <- "GRM"
  administered <- 1:50
  number_dimensions <- 3
  number_items <- 50
  max_number_answer_categories <- 5
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  eta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (max_number_answer_categories - 1))))
  eta[c(1, 5:10), 3:4] <- NA
  eta[c(40:45), 4] <- NA
  beta <- row_cumsum(eta)
  guessing <- NULL
  
  probabilities <- get_probabilities(theta, model, administered, alpha, beta, guessing)
  
  expect_equal(round(probabilities[1,], 3), c(.069, .000, .985, NA, NA))
  expect_equal(round(probabilities[40,], 3), c(.113, .000, .023, .928, NA))
  expect_equal(round(probabilities[50,], 3), c(.261, .000, .000, .036, .933))
  expect_equal(dim(probabilities), c(50, 5))
}) 


context("SM model")

test_that("model is SM, 1 dimension, 2 categories", {
  theta <- -2.3
  model <- "SM"
  administered <- 1:5
  number_dimensions <- 3
  number_items <- 50
  number_answer_categories <- 2
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  
  probabilities <- get_probabilities(theta, model, administered, alpha, beta, guessing)
  
  expect_equal(round(probabilities[1,], 3), c(.575, .425))
  expect_equal(round(probabilities[3,], 3), c(.979, .021))
  expect_equal(dim(probabilities), c(5, 2))
}) 

test_that("model is SM, 3 dimensions, 4 categories", {
  theta <- c(2, 1.2, -2.3)
  model <- "SM"
  administered <- 1:50
  number_dimensions <- 3
  number_items <- 50
  number_answer_categories <- 4
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  beta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (number_answer_categories - 1))))
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  
  probabilities <- get_probabilities(theta, model, administered, alpha, beta, guessing)
  
  expect_equal(round(probabilities[1,], 3), c(.045, .245, .510, .201))
  expect_equal(round(probabilities[40,], 3), c(.094, .393, .436, .077))
  expect_equal(dim(probabilities), c(50, 4))
}) 

test_that("model is SM, 3 dimensions, varying number of categories, estimator is expected_aposteriori, deriv is true", {
  theta <- c(1.8, -1.4, 2.5)
  model <- "SM"
  administered <- 1:50
  number_dimensions <- 3
  number_items <- 50
  max_number_answer_categories <- 5
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  eta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (max_number_answer_categories - 1))))
  eta[c(1, 5:10), 3:4] <- NA
  eta[c(40:45), 4] <- NA
  beta <- row_cumsum(eta)
  guessing <- NULL
  
  probabilities <- get_probabilities(theta, model, administered, alpha, beta, guessing)
  
  expect_equal(round(probabilities[1,], 3), c(.009, .002, .990, NA, NA))
  expect_equal(round(probabilities[40,], 3), c(.017, .007, .010, .966, NA))
  expect_equal(round(probabilities[50,], 3), c(.003, .000, .000, .001, .996))
  expect_equal(dim(probabilities), c(50, 5))
}) 


