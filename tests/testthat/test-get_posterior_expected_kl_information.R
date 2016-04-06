# only for whithin R:
'
library(testthat)
'

make_random_seed_exist <- rnorm(1)

test_that("1 dimensions, 2 categories, estimator is maximum_likelihood", {
  model <- "3PLM"
  number_items <- 50
  number_dimensions <- 1
  estimator <- "maximum_likelihood"
  estimate <- 0
  responses <- rep(c(1, 0), 17)
  administered <- c(6:20, 31:49)
  available <- c(1:5, 21:30, 50)
  prior <- diag(1) * .4
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  lower_bound <- -3
  upper_bound <- 3
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  posterior_expected_kl_information <- get_posterior_expected_kl_information(estimate, model, responses, administered, available, number_dimensions, estimator, alpha, beta, guessing, prior, number_itemsteps_per_item, lower_bound, upper_bound)
  
  expect_equal(length(posterior_expected_kl_information), 16)
  expect_equal(round(posterior_expected_kl_information[c(1, 3, 16)], 15), c(1.554e-12, 8.470e-13, 6.775e-12))
})

test_that("1 dimensions, 2 categories, estimator is expected_aposteriori", {
  model <- "3PLM"
  number_items <- 50
  number_dimensions <- 1
  estimator <- "expected_aposteriori"
  estimate <- 0
  responses <- rep(c(1, 0), 17)
  administered <- c(6:20, 31:49)
  available <- c(1:5, 21:30, 50)
  prior <- diag(1) * .4
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  lower_bound <- -3
  upper_bound <- 3
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  posterior_expected_kl_information <- get_posterior_expected_kl_information(estimate, model, responses, administered, available, number_dimensions, estimator, alpha, beta, guessing, prior, number_itemsteps_per_item, lower_bound, upper_bound)
  
  expect_equal(length(posterior_expected_kl_information), 16)
  expect_equal(round(posterior_expected_kl_information[c(1, 3, 16)], 19), c(4.112e-13, 3.457e-13, 1.753e-12))
})


test_that("model is GPCM, 3 dimensions, varying numbers of categories", {
  model <- "GPCM"
  number_items <- 50
  number_dimensions <- 3
  estimator <- "expected_aposteriori"
  estimate <- 0
  responses <- rep(c(1, 0), 17)
  administered <-  c(6:20, 31:49)
  available <-  c(1:5, 21:30, 50)
  prior <- diag(number_dimensions)
  max_number_answer_categories <- 5
  guessing <- NULL
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  eta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (max_number_answer_categories - 1))))
  eta[c(1, 5:10), 3:4] <- NA
  eta[c(40:45), 4] <- NA
  beta <- row_cumsum(eta)
  lower_bound <- -3
  upper_bound <- 3
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  posterior_expected_kl_information <- get_posterior_expected_kl_information(estimate, model, responses, administered, available, number_dimensions, estimator, alpha, beta, guessing, prior, number_itemsteps_per_item, lower_bound, upper_bound)
  
  expect_equal(length(posterior_expected_kl_information), 16)
  expect_equal(round(posterior_expected_kl_information[c(1, 3, 16)], 19), c(8.690e-17, 5.240e-17, 5.737e-16))
})
