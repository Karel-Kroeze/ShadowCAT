# only for whithin R:
'
library(testthat)
library(MultiGHQuad)
'

make_random_seed_exist <- rnorm(1)

context("estimator is maximum_likelihood")

test_that("estimator is maximum_likelihood, 1 dimension, 2 categories", {
  number_dimensions <- 1
  estimate <- rep(.3, number_dimensions)
  model <- "3PLM"
  number_items <- 50
  responses <- rep(c(1, 0), 17)
  administered <- c(6:20, 31:49)
  estimator <- "maximum_likelihood"
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  lower_bound <- rep(-3, number_dimensions)
  upper_bound <- rep(3, number_dimensions)
  prior <- NULL
  prior_var_safe_nlm <- NULL
  
  estimated_latent_trait <- estimate_latent_trait(estimate, responses, prior, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, prior_var_safe_nlm = NULL)
    
  expect_equal(round(as.vector(estimated_latent_trait), 3), -.799)
  expect_equal(round(attr(estimated_latent_trait, "variance"), 3), matrix(.271))
})


test_that("estimator is maximum_likelihood, 3 dimensions, varying number of categories", {
  number_dimensions <- 3
  estimate <- c(2, .3, -1.2)
  model <- "GPCM"
  number_items <- 50
  responses <- rep(c(1, 0), 17)
  administered <- c(6:20, 31:49)
  estimator <- "maximum_likelihood"
  max_number_answer_categories <- 5
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  eta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (max_number_answer_categories - 1))))
  eta[c(1, 5:10), 3:4] <- NA
  eta[c(40:45), 4] <- NA
  beta <- row_cumsum(eta)
  guessing <- NULL
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  lower_bound <- rep(-3, number_dimensions)
  upper_bound <- rep(3, number_dimensions)
  prior <- NULL
  prior_var_safe_nlm <- NULL
  
  estimated_latent_trait <- estimate_latent_trait(estimate, responses, prior, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, prior_var_safe_nlm = NULL)
  
  expect_equal(round(as.vector(estimated_latent_trait), 3), c(-2.03, -.408, -.449))
  expect_equal(round(attr(estimated_latent_trait, "variance"), 3)[1,], c(.749, -.495, -.181))
  expect_equal(round(attr(estimated_latent_trait, "variance"), 3)[3,], c(-.181, -.281, .484))
})

test_that("test safe_maximum_likelihood", {
  number_dimensions <- 3
  estimate <- c(0, 0, 0)
  attr(estimate, "variance") <- diag(3) * 20
  model <- "SM"
  number_items <- 300
  responses <- c(0, 0, 0, 1, 1, 1, 1, 1, 1, 0)
  administered <- c(38, 27, 4, 182, 187, 200, 241, 274, 227, 298)
  estimator <- "maximum_likelihood"
  alpha <- matrix(with_random_seed(107, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(107, rnorm)(number_items), nrow = number_items, ncol = 1)
  guessing <- NULL
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  lower_bound <- rep(-3, number_dimensions)
  upper_bound <- rep(3, number_dimensions)
  prior <- NULL
  prior_var_safe_nlm <- 1
  
  estimated_latent_trait <- estimate_latent_trait(estimate, responses, prior, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, prior_var_safe_nlm)
  
  expect_equal(round(as.vector(estimated_latent_trait), 3), c(-.635, .808, 2.586))  
  expect_equal(round(as.vector(attr(estimated_latent_trait, 'variance')), 3)[1:3], c(6.870, -3.257, -5.353))
})

context("estimator is maximum_aposteriori")

test_that("estimator is maximum_aposteriori, 1 dimension, 2 categories", { 
  number_dimensions <- 1
  estimate <- rep(.3, number_dimensions)
  model <- "3PLM"
  number_items <- 50
  responses <- rep(c(1, 0), 17)
  administered <- c(6:20, 31:49)
  estimator <- "maximum_aposteriori"
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  lower_bound <- rep(-3, number_dimensions)
  upper_bound <- rep(3, number_dimensions)
  prior <- 1
  prior_var_safe_nlm <- NULL
  
  estimated_latent_trait <- estimate_latent_trait(estimate, responses, prior, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, prior_var_safe_nlm = NULL)
  
  expect_equal(round(as.vector(estimated_latent_trait), 3), -.642)
  expect_equal(round(attr(estimated_latent_trait, "variance"), 3), matrix(.201))
})


test_that("estimator is maximum_aposteriori, 3 dimensions, varying number of categories", {
  number_dimensions <- 3
  estimate <- c(2, .3, -1.2)
  model <- "GPCM"
  number_items <- 50
  responses <- rep(c(1, 0), 17)
  administered <- c(6:20, 31:49)
  estimator <- "maximum_aposteriori"
  max_number_answer_categories <- 5
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  eta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (max_number_answer_categories - 1))))
  eta[c(1, 5:10), 3:4] <- NA
  eta[c(40:45), 4] <- NA
  beta <- row_cumsum(eta)
  guessing <- NULL
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  lower_bound <- rep(-3, number_dimensions)
  upper_bound <- rep(3, number_dimensions)
  prior <- diag(number_dimensions)
  prior_var_safe_nlm <- NULL
  
  estimated_latent_trait <- estimate_latent_trait(estimate, responses, prior, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, prior_var_safe_nlm = NULL)
  
  expect_equal(round(as.vector(estimated_latent_trait), 3), c(-1.447, -.714, -.614))
  expect_equal(round(attr(estimated_latent_trait, "variance"), 3)[1,], c(.343, -.191, -.117))
  expect_equal(round(attr(estimated_latent_trait, "variance"), 3)[3,], c(-.117, -.131, .280))
})

context("estimator is expected_aposteriori")

test_that("estimator is expected_aposteriori, 1 dimension, 2 categories", {
  number_dimensions <- 1
  estimate <- rep(.3, number_dimensions)
  attr(estimate, "variance") <- 1.2
  model <- "3PLM"
  number_items <- 50
  responses <- rep(c(1, 0), 17)
  administered <- c(6:20, 31:49)
  estimator <- "expected_aposteriori"
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  eta <- NULL
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  lower_bound <- rep(-3, number_dimensions)
  upper_bound <- rep(3, number_dimensions)
  prior <- diag(number_dimensions)
  prior_var_safe_nlm <- NULL
  
  estimated_latent_trait <- estimate_latent_trait(estimate, responses, prior, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, prior_var_safe_nlm = NULL)
  
  expect_equal(round(as.vector(estimated_latent_trait), 3), -.689)
  expect_equal(round(attr(estimated_latent_trait, "variance"), 3), matrix(.205))
})

test_that("estimator is expected_aposteriori, 3 dimensions, varying number of categories", {
  number_dimensions <- 3
  estimate <- c(2, .3, -1.2)
  attr(estimate, "variance") <- diag(c(2, 1.2, 1.5))
  model <- "GPCM"
  number_items <- 50
  responses <- rep(c(1, 0), 17)
  administered <- c(6:20, 31:49)
  estimator <- "expected_aposteriori"
  max_number_answer_categories <- 5
  guessing <- NULL
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  eta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (max_number_answer_categories - 1))))
  eta[c(1, 5:10), 3:4] <- NA
  eta[c(40:45), 4] <- NA
  beta <- row_cumsum(eta)
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  lower_bound <- rep(-3, number_dimensions)
  upper_bound <- rep(3, number_dimensions)
  prior <- diag(c(1, 1.3, 1))
  prior_var_safe_nlm <- NULL
  
  estimated_latent_trait <- estimate_latent_trait(estimate, responses, prior, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, prior_var_safe_nlm = NULL)
  
  expect_equal(round(as.vector(estimated_latent_trait), 3), c(-.861, -1.486, -.604))
  expect_equal(round(attr(estimated_latent_trait, "variance"), 3)[1,], c(.348, -.214, -.042))
  expect_equal(round(attr(estimated_latent_trait, "variance"), 3)[3,], c(-.042, -.189, .238))
})

test_that("estimator is expected_aposteriori, 3 dimensions, varying number of categories, 1 administered", {
  number_dimensions <- 3
  estimate <- c(2, .3, -1.2)
  attr(estimate, "variance") <- diag(c(2, 1.2, 1.5))
  model <- "GPCM"
  number_items <- 50
  responses <- 1
  administered <- 1
  estimator <- "expected_aposteriori"
  max_number_answer_categories <- 5
  guessing <- NULL
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  eta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (max_number_answer_categories - 1))))
  eta[c(1, 5:10), 3:4] <- NA
  eta[c(40:45), 4] <- NA
  beta <- row_cumsum(eta)
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  lower_bound <- rep(-3, number_dimensions)
  upper_bound <- rep(3, number_dimensions)
  prior <- diag(c(1, 1.3, 1))
  prior_var_safe_nlm <- NULL
  
  estimated_latent_trait <- estimate_latent_trait(estimate, responses, prior, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, prior_var_safe_nlm = NULL)
  
  expect_equal(round(as.vector(estimated_latent_trait), 3), c(-.347, -.267, -.360))
  expect_equal(round(attr(estimated_latent_trait, "variance"), 3)[1,], c(.941, -.045, -.061))
  expect_equal(round(attr(estimated_latent_trait, "variance"), 3)[3,], c(-.061, -.047, .936))
})

test_that("estimates exceed boundaries", {
  number_dimensions <- 3
  estimate <- c(-10, 10, 10)
  attr(estimate, "variance") <- diag(c(.1, .1, .1))
  model <- "GPCM"
  number_items <- 50
  responses <- rep(1, 50)
  administered <- 1:50
  estimator <- "expected_aposteriori"
  max_number_answer_categories <- 5
  guessing <- NULL
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  eta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (max_number_answer_categories - 1))))
  eta[c(1, 5:10), 3:4] <- NA
  eta[c(40:45), 4] <- NA
  beta <- row_cumsum(eta)
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  lower_bound <- rep(-3, number_dimensions)
  upper_bound <- rep(3, number_dimensions)
  prior <- diag(c(5,5,5))
  prior_var_safe_nlm <- NULL
  
  estimated_latent_trait <- estimate_latent_trait(estimate, responses, prior, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, prior_var_safe_nlm = NULL)
  
  expect_equal(round(as.vector(estimated_latent_trait), 3), c(-3, 3, 3))
  expect_equal(round(attr(estimated_latent_trait, "variance"), 3)[1,], c(0, 0, 0))
  expect_equal(round(attr(estimated_latent_trait, "variance"), 3)[3,], c(0, 0, 0))
})

