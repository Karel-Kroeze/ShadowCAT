# only for whithin R:
'
library(testthat)
library(lpSolve)
'

make_random_seed_exist <- rnorm(1)

context("objective is D (determinant)")

test_that("First 1 item administered, information_summary is D (determinant), one dimension, with padding", {
  model <- "3PLM"
  number_items <- 50
  number_dimensions <- 1
  estimator <- "maximum_aposteriori"
  estimate <- rep(0, number_dimensions)
  information_summary <- "determinant"
  responses <- 1
  available <- 2:50
  administered <- 1
  prior <- .4
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  lower_bound <- -3
  upper_bound <- 3
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  item_information <- get_item_information(information_summary, estimate, model, responses, prior, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, pad = TRUE)

  expect_equal(length(item_information), 50)
  expect_equal(round(item_information[c(1:5, 35:38,47:50)], 3), c(.000, .313, .139, .104, .477, .190, .236, .170, .123, .419, .142, .193, .237))
})

test_that("First 1 item administered, information_summary is D (determinant), three dimensions, with padding", {
  model <- "3PLM"
  number_items <- 50
  number_dimensions <- 3
  estimator <- "maximum_aposteriori"
  estimate <- rep(0, number_dimensions)
  information_summary <- "determinant"
  responses <- 1
  available <- 2:50
  administered <- 1
  prior <- diag(c(.4, .8, 1.5))
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  lower_bound <- rep(-3, number_dimensions)
  upper_bound <- rep(3, number_dimensions)
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  item_information <- get_item_information(information_summary, estimate, model, responses, prior, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, pad = TRUE)
  
  expect_equal(length(item_information), 50)
  expect_equal(round(item_information[c(1:5, 35:38,47:50)], 23), c(0, -4.050e-21, 0, -2.687e-20, -6.800e-22,  3.100e-22, -3.858e-20, -2.400e-22, -5.120e-21,  4.880e-21,  1.100e-22,  1.980e-21, -6.410e-21))
})

test_that("First 10 item administered, information_summary is D (determinant), three dimensions, with padding", {
  model <- "3PLM"
  number_items <- 50
  number_dimensions <- 3
  estimator <- "maximum_aposteriori"
  estimate <- rep(0, number_dimensions)
  information_summary <- "determinant"
  responses <- rep(c(1, 0), 5)
  available <- 11:50
  administered <- 1:10
  prior <- diag(c(.4, .8, 1.5))
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  lower_bound <- rep(-3, number_dimensions)
  upper_bound <- rep(3, number_dimensions)
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  item_information <- get_item_information(information_summary, estimate, model, responses, prior, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, pad = TRUE)
  
  expect_equal(length(item_information), 50)
  expect_equal(round(item_information[c(1:5, 35:38,47:50)], 3), c(.000, .000, .000, .000, .000, .036, .040, .033, .050, .045, .039, .047, .034))
})

context("information_summary is PD (posterior determinant)")

test_that("First 1 item administered, objective is PD (posterior determinant), one dimension, no padding", { 
  model <- "3PLM"
  number_items <- 50
  number_dimensions <- 1
  estimator <- "maximum_aposteriori"
  estimate <- rep(0, number_dimensions)
  information_summary <- "posterior_determinant"
  responses <- 1
  available <- 2:50
  administered <- 1
  prior <- .4
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  lower_bound <- -3
  upper_bound <- 3
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  item_information <- get_item_information(information_summary, estimate, model, responses, prior, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, pad = FALSE)
  
  expect_equal(length(item_information), 49)
  expect_equal(round(item_information[c(1:5, 35:38,47:49)], 3), c(2.813, 2.639, 2.604, 2.977, 2.963, 2.736, 2.670, 2.623, 2.754, 2.642, 2.693, 2.737))
})

test_that("First 1 item administered, information_summary is PD (posterior determinant), three dimensions, no padding", {
  model <- "3PLM"
  number_items <- 50
  number_dimensions <- 3
  estimator <- "maximum_aposteriori"
  estimate <- rep(0, number_dimensions)
  information_summary <- "posterior_determinant"
  responses <- 1
  available <- 2:50
  administered <- 1
  prior <- diag(c(.4, .8, 1.5))
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  lower_bound <- rep(-3, number_dimensions)
  upper_bound <- rep(3, number_dimensions)
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  item_information <- get_item_information(information_summary, estimate, model, responses, prior, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, pad = FALSE)
  
  expect_equal(length(item_information), 49)
  expect_equal(round(item_information[c(1:5, 35:38,47:49)], 3), c(2.906, 2.283, 2.483, 3.047, 3.685, 2.879, 2.293, 2.672, 2.754, 2.645, 2.856, 2.388))
})

test_that("First 10 item administered, information_summary is PD (posterior determinant), three dimensions, no padding", {
  model <- "3PLM"
  number_items <- 50
  number_dimensions <- 3
  estimator <- "maximum_aposteriori"
  estimate <- rep(0, number_dimensions)
  information_summary <- "posterior_determinant"
  responses <- rep(c(1, 0), 5)
  available <- 11:50
  administered <- 1:10
  prior <- diag(c(.4, .8, 1.5))
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  lower_bound <- rep(-3, number_dimensions)
  upper_bound <- rep(3, number_dimensions)
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  item_information <- get_item_information(information_summary, estimate, model, responses, prior, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, pad = FALSE)
  
  expect_equal(length(item_information), 40)
  expect_equal(round(item_information[c(1:5, 35:40)], 3), c(9.370, 8.728, 9.485, 9.083, 8.596, 8.824, 8.596, 9.845, 9.095, 9.439, 8.723))
})

test_that("10 scattered items administered, information_summary is PD (posterior determinant), 1 dimension, with padding", {
  model <- "3PLM"
  number_items <- 50
  number_dimensions <- 1
  estimator <- "maximum_aposteriori"
  estimate <- rep(0, number_dimensions)
  information_summary <- "posterior_determinant"
  responses <- rep(c(1, 0), 5)
  available <- c(1, 3, 4, 6, 7, 9, 10, 12:18, 20:24, 31:50)
  administered <- c(2, 5, 8, 11, 19, 25:30)
  prior <- .4
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  lower_bound <- -3
  upper_bound <- 3
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  item_information <- get_item_information(information_summary, estimate, model, responses, prior, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, pad = TRUE)
  
  expect_equal(length(item_information), 50)
  expect_equal(round(item_information[c(1:5, 35:40)], 3), c(4.374, .000, 4.404, 4.368, .000, 4.454, 4.501, 4.434, 4.387, 4.518, 4.358))
  expect_equal(which(round(item_information, 5) == 0), c(2, 5, 8, 11, 19, 25:30))
})

test_that("10 scattered items administered, information_summary is PD (posterior determinant), 3 dimensions, with padding", {
  model <- "3PLM"
  number_items <- 50
  number_dimensions <- 3
  estimator <- "maximum_aposteriori"
  estimate <- rep(0, number_dimensions)
  information_summary <- "posterior_determinant"
  responses <- rep(c(1, 0), 5)
  available <- c(1, 3, 4, 6, 7, 9, 10, 12:18, 20:24, 31:50)
  administered <- c(2, 5, 8, 11, 19, 25:30)
  prior <- diag(c(.4, .8, 1.5))
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  lower_bound <- rep(-3, number_dimensions)
  upper_bound <- rep(3, number_dimensions)
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  item_information <- get_item_information(information_summary, estimate, model, responses, prior, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, pad = TRUE)
  
  expect_equal(length(item_information), 50)
  expect_equal(round(item_information[c(1:5, 35:40)], 3), c(8.537, 0.000, 8.252, 8.546, 0.000, 8.605, 8.995, 8.263, 8.879, 8.873, 8.934))
  expect_equal(which(round(item_information, 5) == 0), c(2, 5, 8, 11, 19, 25:30))
})

context("information_summary is A (trace)")

test_that("First 1 item administered, objective is A (trace), one dimension, with padding", {
  model <- "3PLM"
  number_items <- 50
  number_dimensions <- 1
  estimator <- "maximum_aposteriori"
  estimate <- rep(0, number_dimensions)
  information_summary <- "trace"
  responses <- 1
  available <- 2:50
  administered <- 1
  prior <- .4
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  lower_bound <- -3
  upper_bound <- 3
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  item_information <- get_item_information(information_summary, estimate, model, responses, prior, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, pad = TRUE)
  
  expect_equal(length(item_information), 50)
  expect_equal(round(item_information[c(1:5, 35:38,47:50)], 3), c(.000, .313, .139, .104, .477, .190, .236, .170, .123, .419, .142, .193, .237))
})

test_that("First 1 item administered, information_summary is A (trace), three dimensions, with padding", {
  model <- "3PLM"
  number_items <- 50
  number_dimensions <- 3
  estimator <- "maximum_aposteriori"
  estimate <- rep(0, number_dimensions)
  information_summary <- "trace"
  responses <- 1
  available <- 2:50
  administered <- 1
  prior <- diag(c(.4, .8, 1.5))
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  lower_bound <- rep(-3, number_dimensions)
  upper_bound <- rep(3, number_dimensions)
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  item_information <- get_item_information(information_summary, estimate, model, responses, prior, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, pad = TRUE)
  expect_equal(length(item_information), 50)
  expect_equal(round(item_information[c(1:5, 35:38,47:50)], 3), c(.000, .488, .101, .190, .691, .266, .417, .107, .258, .684, .295, .401, .176))
})

test_that("First 10 item administered, information_summary is A (trace), three dimensions, with padding", {
  model <- "3PLM"
  number_items <- 50
  number_dimensions <- 3
  estimator <- "maximum_aposteriori"
  estimate <- rep(0, number_dimensions)
  information_summary <- "trace"
  responses <- rep(c(1, 0), 5)
  available <- 11:50
  administered <- 1:10
  prior <- diag(c(.4, .8, 1.5))
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  lower_bound <- rep(-3, number_dimensions)
  upper_bound <- rep(3, number_dimensions)
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  item_information <- get_item_information(information_summary, estimate, model, responses, prior, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, pad = TRUE)
  
  expect_equal(length(item_information), 50)
  expect_equal(round(item_information[c(1:5, 35:38,47:50)], 3), c(.000, .000, .000, .000, .000, 3.678, 3.829, 3.520, 3.670, 4.096, 3.707, 3.813, 3.588))
})

context("objective is PA (posterior trace)")

test_that("First 1 item administered, information_summary is PA (posterior trace), one dimension, with padding", {
  model <- "3PLM"
  number_items <- 50
  number_dimensions <- 1
  estimator <- "maximum_aposteriori"
  estimate <- rep(0, number_dimensions)
  information_summary <- "posterior_trace"
  responses <- 1
  available <- 2:50
  administered <- 1
  prior <- .4
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  lower_bound <- -3
  upper_bound <- 3
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  item_information <- get_item_information(information_summary, estimate, model, responses, prior, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, pad = TRUE)
  
  expect_equal(length(item_information), 50)
  expect_equal(round(item_information[c(1:5, 35:38,47:50)], 3), c(.000, 2.813, 2.639, 2.604, 2.977, 2.690, 2.736, 2.670, 2.623, 2.919, 2.642, 2.693, 2.737))
})

test_that("First 1 item administered, information_summary is PA (posterior trace), three dimensions, with padding", {
  model <- "3PLM"
  number_items <- 50
  number_dimensions <- 3
  estimator <- "maximum_aposteriori"
  estimate <- rep(0, number_dimensions)
  information_summary <- "posterior_trace"
  responses <- 1
  available <- 2:50
  administered <- 1
  prior <- diag(c(.4, .8, 1.5))
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  lower_bound <- rep(-3, number_dimensions)
  upper_bound <- rep(3, number_dimensions)
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  item_information <- get_item_information(information_summary, estimate, model, responses, prior, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, pad = TRUE)
  
  expect_equal(length(item_information), 50)
  expect_equal(round(item_information[c(1:5, 35:38,47:50)], 3), c(.000, 4.904, 4.518, 4.606, 5.107, 4.683, 4.834, 4.524, 4.675, 5.101, 4.711, 4.818, 4.592))
})

test_that("First 10 item administered, information_summary is PA (posterior trace), three dimensions, with padding", {
  model <- "3PLM"
  number_items <- 50
  number_dimensions <- 3
  estimator <- "maximum_aposteriori"
  estimate <- rep(0, number_dimensions)
  information_summary <- "posterior_trace"
  responses <- rep(c(1, 0), 5)
  available <- 11:50
  administered <- 1:10
  prior <- diag(c(.4, .8, 1.5))
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  lower_bound <- rep(-3, number_dimensions)
  upper_bound <- rep(3, number_dimensions)
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  item_information <- get_item_information(information_summary, estimate, model, responses, prior, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, pad = TRUE)
  
  expect_equal(length(item_information), 50)
  expect_equal(round(item_information[c(1:5, 35:38,47:50)], 3), c(.000, .000, .000, .000, .000, 8.095, 8.246, 7.936, 8.087, 8.513, 8.124, 8.230, 8.005))
})

context("information_summary is posterior_expected_kullback_leibler")

test_that("First 1 item administered, information_summary is PEKL, one dimension, with padding", {
  model <- "3PLM"
  number_items <- 50
  number_dimensions <- 1
  estimator <- "maximum_aposteriori"
  estimate <- rep(0, number_dimensions)
  attr(estimate, 'variance') <- .4
  information_summary <- "posterior_expected_kullback_leibler"
  responses <- 1
  available <- 2:50
  administered <- 1
  prior <- .4
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  lower_bound <- -3
  upper_bound <- 3
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  item_information <- get_item_information(information_summary, estimate, model, responses, prior, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, pad = TRUE)
  
  expect_equal(length(item_information), 50)
  expect_equal(round(item_information[c(1:5, 35:38,47:50)], 3), c(.000, .048, .019, .010, .074, .025, .033, .023, .013, .062, .017, .026, .034))
})


context("many item characteristics equal")

test_that("information_summary is PD", {
  model <- "3PLM"
  number_items <- 50
  number_dimensions <- 1
  estimator <- "maximum_aposteriori"
  estimate <- rep(0, number_dimensions)
  information_summary <- "posterior_determinant"
  responses <- c(1, 0)
  available <- 3:50
  administered <- 1:2
  prior <- .4
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(rep(1, number_items * number_dimensions), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(rep(1, number_items), nrow = number_items, ncol = 1)
  
  lower_bound <- -3
  upper_bound <- 3
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  item_information <- get_item_information(information_summary, estimate, model, responses, prior, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, pad = TRUE)
  
  expect_equal(length(item_information), 50)
  expect_equal(round(item_information[c(1:5, 35:38,47:50)], 3), c(.000, .000, 2.917, 2.917, 2.917, 2.880, 2.880, 2.880, 2.880, 2.880, 2.880, 2.880, 2.880))
})

test_that("information_summary is PEKL", {
  model <- "3PLM"
  number_items <- 50
  number_dimensions <- 1
  estimator <- "maximum_aposteriori"
  estimate <- rep(0, number_dimensions)
  attr(estimate, 'variance') <- .4
  information_summary <- "posterior_expected_kullback_leibler"
  responses <- c(1, 0)
  available <- 3:50
  administered <- 1:2
  prior <- .4
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(rep(1, number_items * number_dimensions), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(rep(1, number_items), nrow = number_items, ncol = 1)
  
  lower_bound <- -3
  upper_bound <- 3
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  item_information <- get_item_information(information_summary, estimate, model, responses, prior, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, pad = TRUE)
  
  expect_equal(length(item_information), 50)
  expect_equal(round(item_information[c(1:5, 35:38,47:50)], 3), c(.000, .000, .009, .009, .009, .007, .007, .007, .007, .007, .007, .007, .007))
})

context("invalid input")

test_that("information_summary is of unknown type", {
  model <- "3PLM"
  number_items <- 50
  number_dimensions <- 1
  estimator <- "maximum_aposteriori"
  estimate <- rep(0, number_dimensions)
  information_summary <- "PPP"
  responses <- c(1, 0)
  available <- 3:50
  administered <- 1:2
  prior <- .4
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(rep(1, number_items * number_dimensions), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(rep(1, number_items), nrow = number_items, ncol = 1)
  
  lower_bound <- -3
  upper_bound <- 3
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  item_information <- get_item_information(information_summary, estimate, model, responses, prior, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, pad = TRUE)
    
  expect_equal(item_information$errors$information_summary, "of unknown type")
})

