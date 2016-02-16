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
  number_items <- 50
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  fisher_information <- get_fisher_information(theta, model, number_dimensions, alpha, beta, guessing, number_itemsteps_per_item)
  
  expect_equal(dim(fisher_information), c(1, 1, 50))
  expect_equal(round(fisher_information[,,c(1, 3, 40, 50)], 3), c(.051, .129, .039, .117))
})

test_that("model is 3PLM, 3 dimensions, 2 categories", {
  number_dimensions <- 3
  theta <- c(2.1, -1.1, .7 )
  model <- "3PLM"
  number_items <- 50
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  fisher_information <- get_fisher_information(theta, model, number_dimensions, alpha, beta, guessing, number_itemsteps_per_item)
  
  expect_equal(dim(fisher_information), c(3, 3, 50))
  expect_equal(dim(fisher_information[,,2]), c(3, 3))
  expect_equal(round(fisher_information[,1,2], 3), c(.104, .029, .074))
  expect_equal(round(fisher_information[,2,3], 3), c(.004, .004, .006))
  expect_equal(round(fisher_information[1,,24], 3), c(.002, .002, .006))
  expect_equal(round(fisher_information[3,,40], 3), c(.049, .072, .079))
})

context("GPCM model")

test_that("model is GPCM, 1 dimensions, 2 categories", {
  number_dimensions <- 1
  theta <- -.7
  model <- "GPCM"
  number_items <- 50
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  fisher_information <- get_fisher_information(theta, model, number_dimensions, alpha, beta, guessing, number_itemsteps_per_item)
  
  expect_equal(dim(fisher_information), c(1, 1, 50))
  expect_equal(round(fisher_information[,,c(1,3, 40, 50)], 3), c(.063, .082, .058, .395))
}) 


test_that("model is GPCM, 3 dimensions, varying numbers of categories", {
  number_dimensions <- 3
  theta <- c(2, 1, -2)
  model <- "GPCM"
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
  
  fisher_information <- get_fisher_information(theta, model, number_dimensions, alpha, beta, guessing, number_itemsteps_per_item)
  
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
  number_items <- 50
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  fisher_information <- get_fisher_information(theta, model, number_dimensions, alpha, beta, guessing, number_itemsteps_per_item)
  
  expect_equal(dim(fisher_information), c(1, 1, 50))
  expect_equal(round(fisher_information[,,c(1,3, 40, 50)], 3), c(.044, .217, .052, .129))
})

test_that("model is SM, 3 dimensions, 3 dimensions, 4 categories", {
  number_dimensions <- 3
  theta <- c(-.5, .6, 1.2)
  model <- "SM"
  number_items <- 50
  number_answer_categories <- 4
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  beta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (number_answer_categories - 1))))
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  fisher_information <- get_fisher_information(theta, model, number_dimensions, alpha, beta, guessing, number_itemsteps_per_item)
  
  expect_equal(dim(fisher_information), c(3, 3, 50))
  expect_equal(dim(fisher_information[,,2]), c(3, 3))
  expect_equal(round(fisher_information[,1,2], 3), c(.487, .135, .347))
  expect_equal(round(fisher_information[,2,3], 3), c(.413, .468, .617))
  expect_equal(round(fisher_information[1,,24], 3), c(.085, .092, .241))
  expect_equal(round(fisher_information[3,,40], 3), c(.136, .200, .221))
})

test_that("model is SM, 3 dimensions, varying number of categories", {
  number_dimensions <- 3
  theta <- c(-2, 1, 2)
  model <- "SM"
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
  
  fisher_information <- get_fisher_information(theta, model, number_dimensions, alpha, beta, guessing, number_itemsteps_per_item)
  
  expect_equal(dim(fisher_information), c(3, 3, 50))
  expect_equal(dim(fisher_information[,,2]), c(3, 3))
  expect_equal(round(fisher_information[,1,2], 3), c(.620, .172, .441))
  expect_equal(round(fisher_information[,2,3], 3), c(.416, .471, .621))
  expect_equal(round(fisher_information[1,,24], 3), c(.081, .088, .231))
  expect_equal(round(fisher_information[3,,40], 3), c(.020, .030, .033))
})

context("GRM model")

test_that("model is GRM 1 dimensions, 2 categories", {
  number_dimensions <- 1
  theta <- -2
  model <- "GRM"
  number_items <- 50
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  fisher_information <- get_fisher_information(theta, model, number_dimensions, alpha, beta, guessing, number_itemsteps_per_item)
  
  expect_equal(dim(fisher_information), c(1, 1, 50))
  expect_equal(round(fisher_information[,,c(1,3, 40, 50)], 3), c(.068, .026, .051, .265))
})

test_that("model is GRM, 3 dimensions, varying numbers of categories", {
  number_dimensions <- 3
  theta <- c(2, .1, -1)
  model <- "GRM"
  number_items <- 50
  number_answer_categories <- 4
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  beta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (number_answer_categories - 1))))
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  fisher_information <- get_fisher_information(theta, model, number_dimensions, alpha, beta, guessing, number_itemsteps_per_item)
  
  expect_equal(dim(fisher_information), c(3, 3, 50))
  expect_equal(dim(fisher_information[,,2]), c(3, 3))
  expect_equal(round(fisher_information[,1,2], 3), c(.376, .104, .268))
  expect_equal(round(fisher_information[,2,3], 3), c(.321, .364, .480))
  expect_equal(round(fisher_information[1,,24], 3), c(.059, .064, .169))
  expect_equal(round(fisher_information[3,,40], 3), c(.110, .162, .178))
})

test_that("model is GRM, 3 dimensions, varying number of categories", {
  number_dimensions <- 3
  theta <- c(1.4, -1.4, 2.1)
  model <- "GRM"
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
  
  fisher_information <- get_fisher_information(theta, model, number_dimensions, alpha, beta, guessing, number_itemsteps_per_item)
  
  expect_equal(dim(fisher_information), c(3, 3, 50))
  expect_equal(dim(fisher_information[,,2]), c(3, 3))
  expect_equal(round(fisher_information[,1,2], 3), c(.125, .035, .089))
  expect_equal(round(fisher_information[,2,3], 3), c(.301, .341, .449))
  expect_equal(round(fisher_information[1,,24], 3), c(.061, .067, .174))
  expect_equal(round(fisher_information[3,,40], 3), c(.006, .010, .011))
})

