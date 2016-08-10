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
  answers <- rep(c(1, 0), 25)
  
  probs_and_likelihoods_per_item <- get_probs_and_likelihoods_per_item(theta, model, alpha, beta, guessing, answers, with_likelihoods = TRUE)
    
  expect_equal(round(probs_and_likelihoods_per_item$P[1,], 3), c(.293, .707))
  expect_equal(round(probs_and_likelihoods_per_item$P[50,], 3), c(.083, .917))
  expect_equal(dim(probs_and_likelihoods_per_item$P), c(50, 2))
  expect_equal(round(probs_and_likelihoods_per_item$l[5:7], 3), c(.727, .334, .529))
  expect_equal(round(probs_and_likelihoods_per_item$d[5:7], 3), c(.262, -.629, .425))
  expect_equal(round(probs_and_likelihoods_per_item$D[5:7], 3), c(-.171, -.233, -.160))
  expect_equal(length(probs_and_likelihoods_per_item$l), 50)
  expect_equal(length(probs_and_likelihoods_per_item$d), 50)
  expect_equal(length(probs_and_likelihoods_per_item$D), 50)
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
  answers <- rep(c(1, 0), 25)
  
  probs_and_likelihoods_per_item <- get_probs_and_likelihoods_per_item(theta, model, alpha, beta, guessing, answers, with_likelihoods = TRUE)
  
  expect_equal(round(probs_and_likelihoods_per_item$P[1,], 3), c(.232, .768))
  expect_equal(round(probs_and_likelihoods_per_item$P[50,], 3), c(.028, .972))
  expect_equal(dim(probs_and_likelihoods_per_item$P), c(50, 2))
  expect_equal(round(probs_and_likelihoods_per_item$l[5:7], 3), c(.742, .297, .452))
  expect_equal(round(probs_and_likelihoods_per_item$d[5:7], 3), c(.248, -.670, .474))
  expect_equal(round(probs_and_likelihoods_per_item$D[5:7], 3), c(-.167, -.221, -.121))
  expect_equal(length(probs_and_likelihoods_per_item$l), 50)
  expect_equal(length(probs_and_likelihoods_per_item$d), 50)
  expect_equal(length(probs_and_likelihoods_per_item$D), 50)
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
  answers <- rep(c(1, 0), 25)
  
  probs_and_likelihoods_per_item <- get_probs_and_likelihoods_per_item(theta, model, alpha, beta, guessing, answers, with_likelihoods = TRUE)
  
  expect_equal(round(probs_and_likelihoods_per_item$P[1,], 3), c(.23, .77))
  expect_equal(round(probs_and_likelihoods_per_item$P[50,], 3), c(.123, .877))
  expect_equal(dim(probs_and_likelihoods_per_item$P), c(50, 2))
  expect_equal(round(probs_and_likelihoods_per_item$l[5:7], 3), c(.719, .326, .393))
  expect_equal(round(probs_and_likelihoods_per_item$d[5:7], 3), c(.281, -.674, .607))
  expect_equal(round(probs_and_likelihoods_per_item$D[5:7], 3), c(-.202, -.220, -.239))
  expect_equal(length(probs_and_likelihoods_per_item$l), 50)
  expect_equal(length(probs_and_likelihoods_per_item$d), 50)
  expect_equal(length(probs_and_likelihoods_per_item$D), 50)
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
  answers <- c(rep(c(0, 1, 2), 16), 3, 3)
  
  probs_and_likelihoods_per_item <- get_probs_and_likelihoods_per_item(theta, model, alpha, beta, guessing, answers, with_likelihoods = TRUE)
  
  expect_equal(round(probs_and_likelihoods_per_item$P[1,], 3), c(.033, .352, .514, .102))
  expect_equal(round(probs_and_likelihoods_per_item$P[50,], 3), c(.086, .493, .381, .040))
  expect_equal(dim(probs_and_likelihoods_per_item$P), c(50, 4))
  expect_equal(round(probs_and_likelihoods_per_item$l[33:35], 3), c(.010, .218, .491))
  expect_equal(round(probs_and_likelihoods_per_item$d[33:35], 3), c(1.743, -1.004, -.380))
  expect_equal(round(probs_and_likelihoods_per_item$D[33:35], 3), c(-.211, -.459, -.487))
  expect_equal(length(probs_and_likelihoods_per_item$l), 50)
  expect_equal(length(probs_and_likelihoods_per_item$d), 50)
  expect_equal(length(probs_and_likelihoods_per_item$D), 50)
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
  answers <- c(rep(0,5), rep(1,5), rep(c(1:3),9), rep(c(0:2), 4), 3)
  
  probs_and_likelihoods_per_item <- get_probs_and_likelihoods_per_item(theta, model, alpha, beta, guessing, answers, with_likelihoods = TRUE)
  
  expect_equal(round(probs_and_likelihoods_per_item$P[1,], 3), c(.001, .066, .932, NA, NA))
  expect_equal(round(probs_and_likelihoods_per_item$P[50,], 3), c(.000, .003, .066, .374, .556))
  expect_equal(dim(probs_and_likelihoods_per_item$P), c(50, 5))
  expect_equal(round(probs_and_likelihoods_per_item$l[3:5], 3), c(.006, .000, .002))
  expect_equal(round(probs_and_likelihoods_per_item$d[3:5], 3), c(-2.515, -3.474, -1.903))
  expect_equal(round(probs_and_likelihoods_per_item$D[3:5], 3), c(-.703, -.407, -.093))
  expect_equal(length(probs_and_likelihoods_per_item$l), 50)
  expect_equal(length(probs_and_likelihoods_per_item$d), 50)
  expect_equal(length(probs_and_likelihoods_per_item$D), 50)
})

context("GRM model")

test_that("model is GRM, 1 dimensions, 2 categories", {
  theta <- .5
  model <- "GRM"
  administered <- 1:50
  number_dimensions <- 1 
  number_items <- 50
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  answers <- rep(c(1, 0), 25)
  
  probs_and_likelihoods_per_item <- get_probs_and_likelihoods_per_item(theta, model, alpha, beta, guessing, answers, with_likelihoods = TRUE)
  
  expect_equal(round(probs_and_likelihoods_per_item$P[1,], 3), c(.239, .761))
  expect_equal(round(probs_and_likelihoods_per_item$P[50,], 3), c(.138, .862))
  expect_equal(dim(probs_and_likelihoods_per_item$P), c(50, 2))
  expect_equal(round(probs_and_likelihoods_per_item$l[3:5], 3), c(.251, .201, .689))
  expect_equal(round(probs_and_likelihoods_per_item$d[3:5], 3), c(.749, -0.799, .311))
  expect_equal(round(probs_and_likelihoods_per_item$D[3:5], 3), c(-.188, -.160, -.214))
  expect_equal(length(probs_and_likelihoods_per_item$l), 50)
  expect_equal(length(probs_and_likelihoods_per_item$d), 50)
  expect_equal(length(probs_and_likelihoods_per_item$D), 50)
}) 

test_that("model is GRM, 3 dimensions, 4 categories", {
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
  answers <- c(rep(c(0, 1, 2), 16), 3, 3)
  
  probs_and_likelihoods_per_item <- get_probs_and_likelihoods_per_item(theta, model, alpha, beta, guessing, answers, with_likelihoods = TRUE)
  
  expect_equal(round(probs_and_likelihoods_per_item$P[1,], 3), c(.007, .042, .226, .725))
  expect_equal(round(probs_and_likelihoods_per_item$P[50,], 3), c(.001, .008, .052, .939))
  expect_equal(dim(probs_and_likelihoods_per_item$P), c(50, 4))
  expect_equal(round(probs_and_likelihoods_per_item$l[3:5], 3), c(.062, .000, .032))
  expect_equal(round(probs_and_likelihoods_per_item$d[3:5], 3), c(-.916, -1.000, -.957))
  expect_equal(round(probs_and_likelihoods_per_item$D[3:5], 3), c(-.078, .000, -.042))
  expect_equal(length(probs_and_likelihoods_per_item$l), 50)
  expect_equal(length(probs_and_likelihoods_per_item$d), 50)
  expect_equal(length(probs_and_likelihoods_per_item$D), 50)
}) 

test_that("model is GRM, 3 dimensions, varying number of categories", {
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
  answers <- c(rep(0,5), rep(1,5), rep(c(1:3),9), rep(c(0:2), 4), 3)
  
  probs_and_likelihoods_per_item <- get_probs_and_likelihoods_per_item(theta, model, alpha, beta, guessing, answers, with_likelihoods = TRUE)
  
  expect_equal(round(probs_and_likelihoods_per_item$P[1,], 3), c(.069, .000, .985, NA, NA))
  expect_equal(round(probs_and_likelihoods_per_item$P[50,], 3), c(.261, .000, .000, .036, .933))
  expect_equal(dim(probs_and_likelihoods_per_item$P), c(50, 5))
  expect_equal(round(probs_and_likelihoods_per_item$l[3:5], 3), c(.461, .041, .547))
  expect_equal(round(probs_and_likelihoods_per_item$d[3:5], 3), c(-.539, -.959, -.453))
  expect_equal(round(probs_and_likelihoods_per_item$D[3:5], 3), c(-.248, -.039, -.248))
  expect_equal(length(probs_and_likelihoods_per_item$l), 50)
  expect_equal(length(probs_and_likelihoods_per_item$d), 50)
  expect_equal(length(probs_and_likelihoods_per_item$D), 50)
}) 


context("SM model")

test_that("model is SM, 1 dimension, 2 categories", {
  theta <- -2.3
  model <- "SM"
  administered <- 1:50
  number_dimensions <- 1
  number_items <- 50
  number_answer_categories <- 2
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  answers <- rep(c(1, 0), 25)
  
  probs_and_likelihoods_per_item <- get_probs_and_likelihoods_per_item(theta, model, alpha, beta, guessing, answers, with_likelihoods = TRUE)
  
  expect_equal(round(probs_and_likelihoods_per_item$P[1,], 3), c(.575, .425))
  expect_equal(round(probs_and_likelihoods_per_item$P[50,], 3), c(.849, .151))
  expect_equal(dim(probs_and_likelihoods_per_item$P), c(50, 2))
  expect_equal(round(probs_and_likelihoods_per_item$l[5:7], 3), c(.039, .969, .147))
  expect_equal(round(probs_and_likelihoods_per_item$d[5:7], 3), c(.961, -.031, .853))
  expect_equal(round(probs_and_likelihoods_per_item$D[5:7], 3), c(-.037, -.030, -.126))
  expect_equal(length(probs_and_likelihoods_per_item$l), 50)
  expect_equal(length(probs_and_likelihoods_per_item$d), 50)
  expect_equal(length(probs_and_likelihoods_per_item$D), 50)
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
  
  answers <- c(rep(c(0, 1, 2), 16), 3, 3)
  
  probs_and_likelihoods_per_item <- get_probs_and_likelihoods_per_item(theta, model, alpha, beta, guessing, answers, with_likelihoods = TRUE)
  
  expect_equal(round(probs_and_likelihoods_per_item$P[1,], 3), c(.045, .245, .510, .201))
  expect_equal(round(probs_and_likelihoods_per_item$P[50,], 3), c(.005, .037, .214, .743))
  expect_equal(dim(probs_and_likelihoods_per_item$P), c(50, 4))
  expect_equal(round(probs_and_likelihoods_per_item$l[3:5], 3), c(.091, .054, .097))
  expect_equal(round(probs_and_likelihoods_per_item$d[3:5], 3), c(1.233, -.946, -.888))
  expect_equal(round(probs_and_likelihoods_per_item$D[3:5], 3), c(-.401, -.051, -.103))
  expect_equal(length(probs_and_likelihoods_per_item$l), 50)
  expect_equal(length(probs_and_likelihoods_per_item$d), 50)
  expect_equal(length(probs_and_likelihoods_per_item$D), 50)
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
  answers <- c(rep(0,5), rep(1,5), rep(c(1:3),9), rep(c(0:2), 4), 3)
  
  probs_and_likelihoods_per_item <- get_probs_and_likelihoods_per_item(theta, model, alpha, beta, guessing, answers, with_likelihoods = TRUE)
  
  expect_equal(round(probs_and_likelihoods_per_item$P[1,], 3), c(.009, .002, .990, NA, NA))
  expect_equal(round(probs_and_likelihoods_per_item$P[50,], 3), c(.003, .000, .000, .001, .996))
  expect_equal(dim(probs_and_likelihoods_per_item$P), c(50, 5))
  expect_equal(round(probs_and_likelihoods_per_item$l[3:5], 3), c(.013, .005, .005))
  expect_equal(round(probs_and_likelihoods_per_item$d[3:5], 3), c(-0.987, -0.995, -0.995))
  expect_equal(round(probs_and_likelihoods_per_item$D[3:5], 3), c(-.013, -.005, -.005))
  expect_equal(length(probs_and_likelihoods_per_item$l), 50)
  expect_equal(length(probs_and_likelihoods_per_item$d), 50)
  expect_equal(length(probs_and_likelihoods_per_item$D), 50)
}) 


