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
  answers <- rep(c(1, 0), 17)
  administered <- c(6:20, 31:49)
  estimator <- "maximum_likelihood"
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  prior_form <- NULL
  prior_parameters <- NULL
  
  estimated_latent_trait <- estimate_latent_trait(estimate, answers, prior_form, prior_parameters, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item)
    
  expect_equal(round(as.vector(estimated_latent_trait), 3), -.799)
  expect_equal(round(attr(estimated_latent_trait, "variance"), 3), matrix(.271))
})


test_that("estimator is maximum_likelihood, 3 dimensions, varying number of categories", {
  number_dimensions <- 3
  estimate <- c(2, .3, -1.2)
  model <- "GPCM"
  number_items <- 50
  answers <- rep(c(1, 0), 17)
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
  prior_form <- NULL
  prior_parameters <- NULL
  
  estimated_latent_trait <- estimate_latent_trait(estimate, answers, prior_form, prior_parameters, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item)
  
  expect_equal(round(as.vector(estimated_latent_trait), 3), c(-2.03, -.408, -.449))
  expect_equal(round(attr(estimated_latent_trait, "variance"), 3)[1,], c(.749, -.495, -.181))
  expect_equal(round(attr(estimated_latent_trait, "variance"), 3)[3,], c(-.181, -.281, .484))
})

context("estimator is maximum_aposteriori")

test_that("estimator is maximum_aposteriori, 1 dimension, 2 categories, prior is normal", { 
  number_dimensions <- 1
  estimate <- rep(.3, number_dimensions)
  model <- "3PLM"
  number_items <- 50
  answers <- rep(c(1, 0), 17)
  administered <- c(6:20, 31:49)
  estimator <- "maximum_aposteriori"
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  prior_form <- "normal"
  prior_parameters <- list(mu = 0, Sigma = diag(1))
  
  estimated_latent_trait <- estimate_latent_trait(estimate, answers, prior_form, prior_parameters, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item)
  
  expect_equal(round(as.vector(estimated_latent_trait), 3), -.642)
  expect_equal(round(attr(estimated_latent_trait, "variance"), 3), matrix(.201))
})

test_that("estimator is maximum_aposteriori, 1 dimension, 2 categories, prior is uniform", { 
  number_dimensions <- 1
  estimate <- rep(.3, number_dimensions)
  model <- "3PLM"
  number_items <- 50
  answers <- rep(c(1, 0), 17)
  administered <- c(6:20, 31:49)
  estimator <- "maximum_aposteriori"
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  prior_form <- "uniform"
  prior_parameters <- list(lower_bound = -3, upper_bound = 3)
  
  estimated_latent_trait <- estimate_latent_trait(estimate, answers, prior_form, prior_parameters, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item)
  
  expect_equal(round(as.vector(estimated_latent_trait), 3), -.799)
  expect_equal(round(attr(estimated_latent_trait, "variance"), 3), matrix(.271))
})

test_that("estimator is maximum_aposteriori, 3 dimensions, varying number of categories, prior is normal", {
  number_dimensions <- 3
  estimate <- c(2, .3, -1.2)
  model <- "GPCM"
  number_items <- 50
  answers <- rep(c(1, 0), 17)
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
  prior_form = "normal"
  prior_parameters = list(mu = rep(0, number_dimensions), Sigma = diag(number_dimensions))
  
  estimated_latent_trait <- estimate_latent_trait(estimate, answers, prior_form, prior_parameters, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item)
  
  expect_equal(round(as.vector(estimated_latent_trait), 3), c(-1.447, -.714, -.614))
  expect_equal(round(attr(estimated_latent_trait, "variance"), 3)[1,], c(.343, -.191, -.117))
  expect_equal(round(attr(estimated_latent_trait, "variance"), 3)[3,], c(-.117, -.131, .280))
})

test_that("estimator is maximum_aposteriori, 3 dimensions, varying number of categories, prior is uniform", {
  number_dimensions <- 3
  estimate <- c(2, .3, -1.2)
  model <- "GPCM"
  number_items <- 50
  answers <- rep(c(1, 0), 17)
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
  prior_form = "uniform"
  prior_parameters = list(lower_bound = c(-2, -2, -3), upper_bound = c(3, 3, 3))
  
  estimated_latent_trait <- estimate_latent_trait(estimate, answers, prior_form, prior_parameters, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item)
  
  expect_equal(round(as.vector(estimated_latent_trait), 3), c(-2.000, -.428, -.457))
  expect_equal(round(attr(estimated_latent_trait, "variance"), 3)[1,], c(.747, -.494, -.182))
  expect_equal(round(attr(estimated_latent_trait, "variance"), 3)[3,], c(-.182, -.280, .484))
})

context("estimator is expected_aposteriori")

test_that("estimator is expected_aposteriori, 1 dimension, 2 categories, normal prior", {
  number_dimensions <- 1
  estimate <- rep(.3, number_dimensions)
  attr(estimate, "variance") <- 1.2
  model <- "3PLM"
  number_items <- 50
  answers <- rep(c(1, 0), 17)
  administered <- c(6:20, 31:49)
  estimator <- "expected_aposteriori"
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  eta <- NULL
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  prior_form = "normal"
  prior_parameters = list(mu = 0, Sigma = diag(1))
  
  estimated_latent_trait_gauss_hermite <- estimate_latent_trait(estimate, answers, prior_form, prior_parameters, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, eap_estimation_procedure = "gauss_hermite_quad")
  estimated_latent_trait_riemann <- estimate_latent_trait(estimate, answers, prior_form, prior_parameters, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item,, eap_estimation_procedure = "riemannsum")
  
  expect_equal(round(as.vector(estimated_latent_trait_gauss_hermite), 3), -.689)
  expect_equal(round(as.vector(estimated_latent_trait_riemann), 3), -.689)
  expect_equal(as.vector(round(attr(estimated_latent_trait_gauss_hermite, "variance"), 3)), .205)
  expect_equal(as.vector(round(attr(estimated_latent_trait_riemann, "variance"), 3)), .205)
})

test_that("estimator is expected_aposteriori, 1 dimension, 2 categories, uniform", {
  number_dimensions <- 1
  estimate <- rep(.3, number_dimensions)
  attr(estimate, "variance") <- 1.2
  model <- "3PLM"
  number_items <- 50
  answers <- rep(c(1, 0), 17)
  administered <- c(6:20, 31:49)
  estimator <- "expected_aposteriori"
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  eta <- NULL
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  prior_form = "uniform"
  prior_parameters = list(lower_bound = -3, upper_bound = 3)
  
  estimated_latent_trait_riemann <- estimate_latent_trait(estimate, answers, prior_form, prior_parameters, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item,, eap_estimation_procedure = "riemannsum")
  
  expect_equal(round(as.vector(estimated_latent_trait_riemann), 3), -.902)
  expect_equal(as.vector(round(attr(estimated_latent_trait_riemann, "variance"), 3)), .305)
})

test_that("estimator is expected_aposteriori, 3 dimensions, varying number of categories, normal prior", {
  number_dimensions <- 3
  estimate <- c(2, .3, -1.2)
  attr(estimate, "variance") <- diag(c(2, 1.2, 1.5))
  model <- "GPCM"
  number_items <- 50
  answers <- rep(c(1, 0), 17)
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
  prior_form = "normal"
  prior_parameters = list(mu = c(0, 0, 0), Sigma = diag(c(1, 1.3, 1)))
  
  estimated_latent_trait_gauss_hermite <- estimate_latent_trait(estimate, answers, prior_form, prior_parameters, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, eap_estimation_procedure = "gauss_hermite_quad")
  estimated_latent_trait_riemann <- estimate_latent_trait(estimate, answers, prior_form, prior_parameters, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, eap_estimation_procedure = "riemannsum")
  
  expect_equal(round(as.vector(estimated_latent_trait_gauss_hermite), 3), c(-.861, -1.486, -.604))
  expect_equal(round(as.vector(estimated_latent_trait_riemann), 3), c(-1.536, -0.611, -0.182))
  expect_equal(round(attr(estimated_latent_trait_gauss_hermite, "variance"), 3)[1,], c(.348, -.214, -.042))
  expect_equal(round(attr(estimated_latent_trait_gauss_hermite, "variance"), 3)[3,], c(-.042, -.189, .238))
  expect_equal(round(attr(estimated_latent_trait_riemann, "variance"), 3)[1,], c(.002, -.002, .000))
  expect_equal(round(attr(estimated_latent_trait_riemann, "variance"), 3)[3,], c(-.000, -.005, .006))
})

test_that("estimator is expected_aposteriori, 3 dimensions, varying number of categories, uniform prior", {
  number_dimensions <- 3
  estimate <- c(2, .3, -1.2)
  attr(estimate, "variance") <- diag(c(2, 1.2, 1.5))
  model <- "GPCM"
  number_items <- 50
  answers <- rep(c(1, 0), 17)
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
  prior_form = "uniform"
  prior_parameters = list(lower_bound = c(-3, -3, -3), upper_bound = c(4, 4, 4))
  
  estimated_latent_trait_riemann <- estimate_latent_trait(estimate, answers, prior_form, prior_parameters, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, eap_estimation_procedure = "riemannsum")
  
  expect_equal(round(as.vector(estimated_latent_trait_riemann), 3), c(-1.376, -.625, -.614))
  expect_equal(round(attr(estimated_latent_trait_riemann, "variance"), 3)[1,], c(.067, -.061, -.032))
  expect_equal(round(attr(estimated_latent_trait_riemann, "variance"), 3)[3,], c(-.032, -.454, .535))
})

test_that("estimator is expected_aposteriori, 3 dimensions, varying number of categories, 1 administered", {
  number_dimensions <- 3
  estimate <- c(2, .3, -1.2)
  attr(estimate, "variance") <- diag(c(2, 1.2, 1.5))
  model <- "GPCM"
  number_items <- 50
  answers <- 1
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
  prior_form = "normal"
  prior_parameters = list(mu = c(0, 0, 0), Sigma = diag(c(1, 1.3, 1)))
  
  estimated_latent_trait_gauss_hermite <- estimate_latent_trait(estimate, answers, prior_form, prior_parameters, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, eap_estimation_procedure = "gauss_hermite_quad")
  estimated_latent_trait_riemann <- estimate_latent_trait(estimate, answers, prior_form, prior_parameters, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, eap_estimation_procedure = "riemannsum")
  
  expect_equal(round(as.vector(estimated_latent_trait_gauss_hermite), 3), c(-.347, -.267, -.360))
  expect_equal(round(as.vector(estimated_latent_trait_riemann), 3), c(-.355, -.273, -.368))
  
  expect_equal(round(attr(estimated_latent_trait_gauss_hermite, "variance"), 3)[1,], c(.941, -.045, -.061))
  expect_equal(round(attr(estimated_latent_trait_gauss_hermite, "variance"), 3)[3,], c(-.061, -.047, .936))
  expect_equal(round(attr(estimated_latent_trait_riemann, "variance"), 3)[1,], c(.946, -.046, -.060))
  expect_equal(round(attr(estimated_latent_trait_riemann, "variance"), 3)[3,], c(-.060, -.047, .940))
})

test_that("test safe_eap", {
  number_dimensions <- 3
  estimate <- c(0, 0, 0)
  attr(estimate, "variance") <- diag(c(0, 0, 0))
  model <- "GPCM"
  number_items <- 300
  answers <- rep(1, 300)
  administered <- 1:300
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
  prior_form = "uniform"
  prior_parameters = list(lower_bound = rep(-1, 3), upper_bound = rep(1, 3))
  safe_eap = TRUE
  
  estimated_latent_trait_riemann <- estimate_latent_trait(estimate, answers, prior_form, prior_parameters, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, safe_eap = safe_eap, eap_estimation_procedure = "riemannsum")
  
  expect_equal(round(as.vector(estimated_latent_trait_riemann), 3), c(-.540, -0.368, -0.748))
  expect_equal(round(attr(estimated_latent_trait_riemann, "variance"), 3)[1,], c(.034, -.016, -.015))
  expect_equal(round(attr(estimated_latent_trait_riemann, "variance"), 3)[3,], c(-.015, -.017, .034))
})

