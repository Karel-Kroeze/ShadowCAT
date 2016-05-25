# only for whithin R:
'
library(testthat)
library(Matrix)
library(mvtnorm)
'

answers <- rep(c(1, 0), 20)
model <- "3PLM"
administered <- 1:40
number_items <- 50
beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  

context("normal prior")

test_that("normal prior, 1 dimension, without adapt", {
  number_dimensions <- 1
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  eap <- get_eap_estimate_riemannsum(dimension = number_dimensions, 
                                     likelihood = likelihood_or_post_density, 
                                     prior_form = "normal", 
                                     prior_parameters = list(mu = rep(0, number_dimensions), Sigma = diag(number_dimensions) * 9), 
                                     adapt = NULL, 
                                     number_gridpoints = 50,
                                     answers = answers, model = model, items_to_include = administered, number_dimensions = number_dimensions, estimator = "maximum_likelihood", alpha = alpha, beta = beta, guessing = guessing, return_log_likelihood_or_post_density = FALSE)
  expect_equal(as.vector(round(eap, 3)), -.941)
  expect_equal(round(attr(eap, "variance"), 3), matrix(.731))
})

test_that("normal prior, 3 dimensions, without adapt", {
  number_dimensions <- 3
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  eap <- get_eap_estimate_riemannsum(dimension = number_dimensions, 
                                     likelihood = likelihood_or_post_density, 
                                     prior_form = "normal", 
                                     prior_parameters = list(mu = rep(0, number_dimensions), Sigma = diag(number_dimensions) * 9), 
                                     adapt = NULL, 
                                     number_gridpoints = 6,
                                     answers = answers, model = model, items_to_include = administered, number_dimensions = number_dimensions, estimator = "maximum_likelihood", alpha = alpha, beta = beta, guessing = guessing, return_log_likelihood_or_post_density = FALSE)
  expect_equal(as.vector(round(eap, 3)), c(-.203, -2.428, -5.071))
  expect_equal(round(attr(eap, "variance"), 3)[1,], c(11.843, -1.192, -7.408))
  expect_equal(round(attr(eap, "variance"), 3)[2,], c(-1.192,  4.975, -2.696))
  expect_equal(round(attr(eap, "variance"), 3)[3,], c(-7.408, -2.696, 10.284))
})

test_that("normal prior, 1 dimension, with adapt", {
  number_dimensions <- 1
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  eap <- get_eap_estimate_riemannsum(dimension = number_dimensions, 
                                     likelihood = likelihood_or_post_density, 
                                     prior_form = "normal", 
                                     prior_parameters = list(mu = rep(0, number_dimensions), Sigma = diag(number_dimensions) * 9), 
                                     adapt = list(mu = -1, Sigma = diag(1) * .8), 
                                     number_gridpoints = 50,
                                     answers = answers, model = model, items_to_include = administered, number_dimensions = number_dimensions, estimator = "maximum_likelihood", alpha = alpha, beta = beta, guessing = guessing, return_log_likelihood_or_post_density = FALSE)
  expect_equal(as.vector(round(eap, 3)), -.933)
  expect_equal(round(attr(eap, "variance"), 3), matrix(.685))
})

test_that("normal prior, 3 dimensions, with adapt", {
  number_dimensions <- 3
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  eap <- get_eap_estimate_riemannsum(dimension = number_dimensions, 
                                     likelihood = likelihood_or_post_density, 
                                     prior_form = "normal", 
                                     prior_parameters = list(mu = rep(0, number_dimensions), Sigma = diag(number_dimensions) * 9), 
                                     adapt = list(mu = c(-2, -2.5, -5), Sigma = diag(3) * c(12, 5, 10)), 
                                     number_gridpoints = 6,
                                     answers = answers, model = model, items_to_include = administered, number_dimensions = number_dimensions, estimator = "maximum_likelihood", alpha = alpha, beta = beta, guessing = guessing, return_log_likelihood_or_post_density = FALSE)
  expect_equal(as.vector(round(eap, 3)), c(3.032, -3.133, -6.609))
  expect_equal(round(attr(eap, "variance"), 3)[1,], c(19.022, -3.136, -9.397))
  expect_equal(round(attr(eap, "variance"), 3)[2,], c(-3.136, 4.303, -1.704))
  expect_equal(round(attr(eap, "variance"), 3)[3,], c(-9.397, -1.704, 9.164))
})

context("uniform prior")

test_that("uniform prior, 1 dimension, without adapt", {
  number_dimensions <- 1
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  eap <- get_eap_estimate_riemannsum(dimension = number_dimensions, 
                                     likelihood = likelihood_or_post_density, 
                                     prior_form = "uniform", 
                                     prior_parameters = list(lower_bound = -3, upper_bound = 3), 
                                     adapt = NULL, 
                                     number_gridpoints = 50,
                                     answers = answers, model = model, items_to_include = administered, number_dimensions = number_dimensions, estimator = "maximum_likelihood", alpha = alpha, beta = beta, guessing = guessing, return_log_likelihood_or_post_density = FALSE)
  expect_equal(as.vector(round(eap, 3)), -.934)
  expect_equal(round(attr(eap, "variance"), 3), matrix(.564))
})

test_that("uniform prior, 3 dimensions, without adapt", {
  number_dimensions <- 3
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  eap <- get_eap_estimate_riemannsum(dimension = number_dimensions, 
                                     likelihood = likelihood_or_post_density, 
                                     prior_form = "uniform", 
                                     prior_parameters = list(lower_bound = rep(-3, 3), upper_bound = rep(3, 3)), 
                                     adapt = NULL, 
                                     number_gridpoints = 6,
                                     answers = answers, model = model, items_to_include = administered, number_dimensions = number_dimensions, estimator = "maximum_likelihood", alpha = alpha, beta = beta, guessing = guessing, return_log_likelihood_or_post_density = FALSE)
  expect_equal(as.vector(round(eap, 3)), c(-1.329, -2.037, -2.271))
  expect_equal(round(attr(eap, "variance"), 3)[1,], c(1.930, -.193, -.120))
  expect_equal(round(attr(eap, "variance"), 3)[2,], c(-.193, .493, -.043))
  expect_equal(round(attr(eap, "variance"), 3)[3,], c(-.120, -.043, .235))
})

test_that("uniform prior, 1 dimension, with adapt", {
  number_dimensions <- 1
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  eap <- get_eap_estimate_riemannsum(dimension = number_dimensions, 
                                     likelihood = likelihood_or_post_density, 
                                     prior_form = "uniform", 
                                     prior_parameters = list(lower_bound = -3, upper_bound = 3), 
                                     adapt = list(mu = 0, Sigma = diag(1)), 
                                     number_gridpoints = 50,
                                     answers = answers, model = model, items_to_include = administered, number_dimensions = number_dimensions, estimator = "maximum_likelihood", alpha = alpha, beta = beta, guessing = guessing, return_log_likelihood_or_post_density = FALSE)
  expect_equal(as.vector(round(eap, 3)), -.934)
  expect_equal(round(attr(eap, "variance"), 3), matrix(.564))
})

test_that("uniform prior, 3 dimensions, with adapt", {
  number_dimensions <- 3
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  eap <- get_eap_estimate_riemannsum(dimension = number_dimensions, 
                                     likelihood = likelihood_or_post_density, 
                                     prior_form = "uniform", 
                                     prior_parameters = list(lower_bound = rep(-3, 3), upper_bound = rep(3, 3)), 
                                     adapt = list(mu = c(0, .5, -.3), Sigma = diag(3) * 2), 
                                     number_gridpoints = 6,
                                     answers = answers, model = model, items_to_include = administered, number_dimensions = number_dimensions, estimator = "maximum_likelihood", alpha = alpha, beta = beta, guessing = guessing, return_log_likelihood_or_post_density = FALSE)
  expect_equal(as.vector(round(eap, 3)), c(-1.385, -1.424, -2.343))
  expect_equal(round(attr(eap, "variance"), 3)[1,], c(1.184, -0.085, -0.037))
  expect_equal(round(attr(eap, "variance"), 3)[2,], c(-.085, .257, -.012))
  expect_equal(round(attr(eap, "variance"), 3)[3,], c(-.037, -.012, .106))
})

test_that("uniform prior, 3 dimensions, with adapt with large variances", {
  number_dimensions <- 3
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  eap <- get_eap_estimate_riemannsum(dimension = number_dimensions, 
                                     likelihood = likelihood_or_post_density, 
                                     prior_form = "uniform", 
                                     prior_parameters = list(lower_bound = rep(-3, 3), upper_bound = rep(3, 3)), 
                                     adapt = list(mu = c(0, .5, -.3), Sigma = diag(3) * 150), 
                                     number_gridpoints = 6,
                                     answers = answers, model = model, items_to_include = administered, number_dimensions = number_dimensions, estimator = "maximum_likelihood", alpha = alpha, beta = beta, guessing = guessing, return_log_likelihood_or_post_density = FALSE)
  expect_equal(as.vector(round(eap, 3)), c(-2.572, -2.375, -1.271))
  expect_equal(round(attr(eap, "variance"), 3)[1,], c(.944, -.041, -.011))
  expect_equal(round(attr(eap, "variance"), 3)[2,], c(-.041, .235, -.004))
  expect_equal(round(attr(eap, "variance"), 3)[3,], c(-.011, -.004, .058))
})

  