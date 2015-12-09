# only for whithin R:
'
library(testthat)
'

make_random_seed_exist <- rnorm(1)

context("3PLM model")

test_that("model is 3PLM, 1 dimensions, 2 categories", {
  number_dimensions <- 1
  theta <- .5
  model <- "3PLM"
  estimator <- "ML"
  number_items <- 50
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  fisher_information <- get_fisher_information(theta, model, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item)
  
  expect_equal(dim(fisher_information), c(1, 1, 50))
  expect_equal(round(fisher_information[,,c(1, 3, 40, 50)], 3), c(.051, .129, .039, .117))
})

test_that("model is 3PLM, 3 dimensions, 2 categories", {
  number_dimensions <- 3
  theta <- c(2.1, -1.1, .7 )
  model <- "3PLM"
  estimator <- "ML"
  number_items <- 50
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  fisher_information <- get_fisher_information(theta, model, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item)
  
  expect_equal(dim(fisher_information), c(3, 3, 50))
  expect_equal(dim(fisher_information[,,2]), c(3, 3))
  expect_equal(round(fisher_information[,1,2], 3), c(.104, .029, .074))
  expect_equal(round(fisher_information[,2,3], 3), c(.004, .004, .006))
  expect_equal(round(fisher_information[1,,24], 3), c(.002, .002, .006))
  expect_equal(round(fisher_information[3,,40], 3), c(.049, .072, .079))
})

test_that("model is GPCM, 1 dimensions, 2 categories", {
  number_dimensions <- 1
  theta <- -.7
  model <- "GPCM"
  estimator <- "MAP"
  number_items <- 50
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  fisher_information <- get_fisher_information(theta, model, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item)
  
  expect_equal(dim(fisher_information), c(1, 1, 50))
  expect_equal(round(fisher_information[,,c(1,3, 40, 50)], 3), c(.063, .082, .058, .395))
}) 


test_that("model is GPCM, 3 dimensions, varying numbers of categories", {
  number_dimensions <- 3
  theta <- c(2, 1, -2)
  model <- "GPCM"
  estimator <- "ML"
  number_items <- 50
  max_number_answer_categories <- 5
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  eta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (max_number_answer_categories - 1))))
  eta[c(1, 5:10), 3:4] <- NA
  eta[c(40:45), 4] <- NA
  beta <- row_cumsum(eta)
  guessing <- NULL
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  fisher_information <- get_fisher_information(theta, model, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item)
  
  expect_equal(dim(fisher_information), c(3, 3, 50))
  expect_equal(dim(fisher_information[,,2]), c(3, 3))
  expect_equal(round(fisher_information[,1,2], 3), c(.910, .253, .647))
  expect_equal(round(fisher_information[,2,3], 3), c(.673, .763, 1.006))
  expect_equal(round(fisher_information[1,,24], 3), c(.054, .058, .153))
  expect_equal(round(fisher_information[3,,40], 3), c(.209, .307, .338))
})

context("SM model")

test_that("model is SM 1 dimensions, 2 categories", {
  number_dimensions <- 1
  theta <- .9
  model <- "SM"
  estimator <- "ML"
  number_items <- 50
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  fisher_information <- get_fisher_information(theta, model, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item)
  
  expect_equal(dim(fisher_information), c(1, 1, 50))
  expect_equal(round(fisher_information[,,c(1,3, 40, 50)], 3), c(.044, .217, .052, .129))
})

test_that("model is SM, 3 dimensions, 3 dimensions, 4 categories", {
  number_dimensions <- 3
  theta <- c(-.5, .6, 1.2)
  model <- "SM"
  estimator <- "ML"
  number_items <- 50
  number_answer_categories <- 4
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  beta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (number_answer_categories - 1))))
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  fisher_information <- get_fisher_information(theta, model, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item)
  
  expect_equal(dim(fisher_information), c(3, 3, 50))
  expect_equal(dim(fisher_information[,,2]), c(3, 3))
  expect_equal(round(fisher_information[,1,2], 3), c(.404, .112, .288))
  expect_equal(round(fisher_information[,2,3], 3), c(.444, .504, .664))
  expect_equal(round(fisher_information[1,,24], 3), c(.084, .091, .237))
  expect_equal(round(fisher_information[3,,40], 3), c(.137, .202, .222))
})

test_that("model is SM, 3 dimensions, varying number of categories", {
  number_dimensions <- 3
  theta <- c(-2, 1, 2)
  model <- "SM"
  estimator <- "ML"
  number_items <- 50
  max_number_answer_categories <- 5
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  eta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (max_number_answer_categories - 1))))
  eta[c(1, 5:10), 3:4] <- NA
  eta[c(40:45), 4] <- NA
  beta <- row_cumsum(eta)
  guessing <- NULL
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  fisher_information <- get_fisher_information(theta, model, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item)
  
  expect_equal(dim(fisher_information), c(3, 3, 50))
  expect_equal(dim(fisher_information[,,2]), c(3, 3))
  expect_equal(round(fisher_information[,1,2], 3), c(.172, .048, .123))
  expect_equal(round(fisher_information[,2,3], 3), c(.558, .633, .834))
  expect_equal(round(fisher_information[1,,24], 3), c(.057, .062, .162))
  expect_equal(round(fisher_information[3,,40], 3), c(.022, .032, .035))
})

context("GRM model")

test_that("model is GRM 1 dimensions, 2 categories", {
  number_dimensions <- 1
  theta <- -2
  model <- "GRM"
  estimator <- "ML"
  number_items <- 50
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  fisher_information <- get_fisher_information(theta, model, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item)
  
  expect_equal(dim(fisher_information), c(1, 1, 50))
  expect_equal(round(fisher_information[,,c(1,3, 40, 50)], 3), c(.068, .026, .051, .265))
})

test_that("model is GRM, 3 dimensions, varying numbers of categories", {
  number_dimensions <- 3
  theta <- c(2, .1, -1)
  model <- "GRM"
  estimator <- "ML"
  number_items <- 50
  number_answer_categories <- 4
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  beta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (number_answer_categories - 1))))
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  fisher_information <- get_fisher_information(theta, model, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item)
  
  expect_equal(dim(fisher_information), c(3, 3, 50))
  expect_equal(dim(fisher_information[,,2]), c(3, 3))
  expect_equal(round(fisher_information[,1,2], 3), c(.326, .091, .232))
  expect_equal(round(fisher_information[,2,3], 3), c(.281, .318, .420))
  expect_equal(round(fisher_information[1,,24], 3), c(.043, .046, .121))
  expect_equal(round(fisher_information[3,,40], 3), c(.110, .162, .178))
})

test_that("model is GRM, 3 dimensions, varying number of categories", {
  number_dimensions <- 3
  theta <- c(1.4, -1.4, 2.1)
  model <- "GRM"
  estimator <- "ML"
  number_items <- 50
  max_number_answer_categories <- 5
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  eta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (max_number_answer_categories - 1))))
  eta[c(1, 5:10), 3:4] <- NA
  eta[c(40:45), 4] <- NA
  beta <- row_cumsum(eta)
  guessing <- NULL
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  fisher_information <- get_fisher_information(theta, model, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item)
  
  expect_equal(dim(fisher_information), c(3, 3, 50))
  expect_equal(dim(fisher_information[,,2]), c(3, 3))
  expect_equal(round(fisher_information[,1,2], 3), c(.335, .093, .238))
  expect_equal(round(fisher_information[,2,3], 3), c(.307, .349, .459))
  expect_equal(round(fisher_information[1,,24], 3), c(.061, .066, .173))
  expect_equal(round(fisher_information[3,,40], 3), c(.008, .012, .013))
})

