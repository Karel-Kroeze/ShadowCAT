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
  answers <- 1
  available <- 2:50
  administered <- 1
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  item_information_normal <- get_summarized_information(information_summary, estimate, model, answers, prior_form = "normal", prior_parameters = NULL, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, pad = TRUE)
  item_information_uniform <- get_summarized_information(information_summary, estimate, model, answers, prior_form = "uniform", prior_parameters = NULL, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, pad = TRUE)
  
  expect_equal(length(item_information_normal), 50)
  expect_equal(round(item_information_normal[c(1:5, 35:38,47:50)], 3), c(.000, .313, .139, .104, .477, .190, .236, .170, .123, .419, .142, .193, .237))
  expect_equal(item_information_uniform, item_information_normal)
})

test_that("First 1 item administered, information_summary is D (determinant), three dimensions, with padding", {
  model <- "3PLM"
  number_items <- 50
  number_dimensions <- 3
  estimator <- "maximum_aposteriori"
  estimate <- rep(0, number_dimensions)
  information_summary <- "determinant"
  answers <- 1
  available <- 2:50
  administered <- 1
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  item_information_normal <- get_summarized_information(information_summary, estimate, model, answers, prior_form = "normal", prior_parameters = NULL, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, pad = TRUE)
  item_information_uniform <- get_summarized_information(information_summary, estimate, model, answers, prior_form = "uniform", prior_parameters = NULL, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, pad = TRUE)
  
  expect_equal(length(item_information_normal), 50)
  expect_equal(round(item_information_normal[c(1:5, 35:38,47:50)], 23), c(0, 2.429e-20, 0.000e+00, 4.240e-21, 5.622e-20, 3.500e-22, -3.813e-20, -4.600e-22, 8.050e-21, 3.700e-22, 1.100e-22, 1.980e-21, -2.070e-21))
  expect_equal(item_information_uniform, item_information_normal)
})

test_that("First 10 item administered, information_summary is D (determinant), three dimensions, with padding", {
  model <- "3PLM"
  number_items <- 50
  number_dimensions <- 3
  estimator <- "maximum_aposteriori"
  estimate <- rep(0, number_dimensions)
  information_summary <- "determinant"
  answers <- rep(c(1, 0), 5)
  available <- 11:50
  administered <- 1:10
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  item_information_normal <- get_summarized_information(information_summary, estimate, model, answers, prior_form = "normal", prior_parameters = NULL, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, pad = TRUE)
  item_information_uniform <- get_summarized_information(information_summary, estimate, model, answers, prior_form = "uniform", prior_parameters = NULL, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, pad = TRUE)
  
  expect_equal(length(item_information_normal), 50)
  expect_equal(round(item_information_normal[c(1:5, 35:38,47:50)], 3), c(.000, .000, .000, .000, .000, .036, .040, .033, .050, .045, .039, .047, .034))
  expect_equal(item_information_uniform, item_information_normal)
})

context("information_summary is PD (posterior determinant)")

test_that("First 1 item administered, objective is PD (posterior determinant), one dimension, no padding", { 
  model <- "3PLM"
  number_items <- 50
  number_dimensions <- 1
  estimator <- "maximum_aposteriori"
  estimate <- rep(0, number_dimensions)
  information_summary <- "posterior_determinant"
  answers <- 1
  available <- 2:50
  administered <- 1
  prior_parameters_normal <- list(mu = 0, Sigma = diag(1) * .4)
  prior_parameters_uniform <- list(lower_bound = -3, upper_bound = 3)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  item_information_normal <- get_summarized_information(information_summary, estimate, model, answers, prior_form = "normal", prior_parameters = prior_parameters_normal, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, pad = FALSE)
  item_information_uniform <- get_summarized_information(information_summary, estimate, model, answers, prior_form = "uniform", prior_parameters = prior_parameters_uniform, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, pad = FALSE)
  
  expect_equal(length(item_information_normal), 49)
  expect_equal(round(item_information_normal[c(1:5, 35:38,47:49)], 3), c(2.813, 2.639, 2.604, 2.977, 2.963, 2.736, 2.670, 2.623, 2.754, 2.642, 2.693, 2.737))
  expect_equal(round(item_information_uniform[c(1:5, 35:38,47:49)], 3), c(.313, .139, .104, .477, .463, .236, .170, .123, .254, .142, .193, .237))
})

test_that("First 1 item administered, information_summary is PD (posterior determinant), three dimensions, no padding", {
  model <- "3PLM"
  number_items <- 50
  number_dimensions <- 3
  estimator <- "maximum_aposteriori"
  estimate <- rep(0, number_dimensions)
  information_summary <- "posterior_determinant"
  answers <- 1
  available <- 2:50
  administered <- 1
  prior_parameters_normal <- list(mu = 0, Sigma = diag(c(.4, .8, 1.5)))
  prior_parameters_uniform <- list(lower_bound = rep(-3, 3), upper_bound = rep(3, 3))
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  item_information_normal <- get_summarized_information(information_summary, estimate, model, answers, prior_form = "normal", prior_parameters = prior_parameters_normal, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, pad = FALSE)
  item_information_uniform <- get_summarized_information(information_summary, estimate, model, answers, prior_form = "uniform", prior_parameters = prior_parameters_uniform, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, pad = FALSE)
  
  expect_equal(length(item_information_normal), 49)
  expect_equal(round(item_information_normal[c(1:5, 35:38,47:49)], 3), c(2.906, 2.283, 2.483, 3.047, 3.685, 2.879, 2.293, 2.672, 2.754, 2.645, 2.856, 2.388))
  expect_equal(round(item_information_uniform[c(1:5, 35:38,47:49)], 23), c(2.429e-20, 0.000e+00, 4.240e-21, 5.622e-20, -6.5983e-19, -3.813e-20, -4.600e-22, 8.050e-21, -1.2727e-19, 1.100e-22, 1.980e-21, -2.070e-21))
})

test_that("First 10 item administered, information_summary is PD (posterior determinant), three dimensions, no padding", {
  model <- "3PLM"
  number_items <- 50
  number_dimensions <- 3
  estimator <- "maximum_aposteriori"
  estimate <- rep(0, number_dimensions)
  information_summary <- "posterior_determinant"
  answers <- rep(c(1, 0), 5)
  available <- 11:50
  administered <- 1:10
  prior_parameters_normal <- list(mu = 0, Sigma = diag(c(.4, .8, 1.5)))
  prior_parameters_uniform <- list(lower_bound = rep(-3, 3), upper_bound = rep(3, 3))
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  item_information_normal <- get_summarized_information(information_summary, estimate, model, answers, prior_form = "normal", prior_parameters = prior_parameters_normal, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, pad = FALSE)
  item_information_uniform <- get_summarized_information(information_summary, estimate, model, answers, prior_form = "uniform", prior_parameters = prior_parameters_uniform, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, pad = FALSE)
  
  expect_equal(length(item_information_normal), 40)
  expect_equal(round(item_information_normal[c(1:5, 35:40)], 3), c(9.370, 8.728, 9.485, 9.083, 8.596, 8.824, 8.596, 9.845, 9.095, 9.439, 8.723))
  expect_equal(round(item_information_uniform[c(1:5, 35:40)], 3), c(.048, .035, .038, .045, .033, .034, .033, .045, .039, .047, .034))
})

test_that("10 scattered items administered, information_summary is PD (posterior determinant), 1 dimension, with padding", {
  model <- "3PLM"
  number_items <- 50
  number_dimensions <- 1
  estimator <- "maximum_aposteriori"
  estimate <- rep(0, number_dimensions)
  information_summary <- "posterior_determinant"
  answers <- rep(c(1, 0), 5)
  available <- c(1, 3, 4, 6, 7, 9, 10, 12:18, 20:24, 31:50)
  administered <- c(2, 5, 8, 11, 19, 25:30)
  prior_parameters_normal <- list(mu = 0, Sigma = diag(1) * .4)
  prior_parameters_uniform <- list(lower_bound = -3, upper_bound = 3)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  item_information_normal <- get_summarized_information(information_summary, estimate, model, answers, prior_form = "normal", prior_parameters = prior_parameters_normal, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, pad = TRUE)
  item_information_uniform <- get_summarized_information(information_summary, estimate, model, answers, prior_form = "uniform", prior_parameters = prior_parameters_uniform, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, pad = TRUE)
  
  expect_equal(length(item_information_normal), 50)
  expect_equal(round(item_information_normal[c(1:5, 35:40)], 3), c(4.374, .000, 4.404, 4.368, .000, 4.454, 4.501, 4.434, 4.387, 4.518, 4.358))
  expect_equal(round(item_information_uniform[c(1:5, 35:40)], 3), c(1.874, .000, 1.904, 1.868, .000, 1.954, 2.001, 1.934, 1.887, 2.018, 1.858))
  expect_equal(which(round(item_information_normal, 5) == 0), c(2, 5, 8, 11, 19, 25:30))
  expect_equal(which(round(item_information_uniform, 5) == 0), c(2, 5, 8, 11, 19, 25:30))
})

test_that("10 scattered items administered, information_summary is PD (posterior determinant), 3 dimensions, with padding", {
  model <- "3PLM"
  number_items <- 50
  number_dimensions <- 3
  estimator <- "maximum_aposteriori"
  estimate <- rep(0, number_dimensions)
  information_summary <- "posterior_determinant"
  answers <- rep(c(1, 0), 5)
  available <- c(1, 3, 4, 6, 7, 9, 10, 12:18, 20:24, 31:50)
  administered <- c(2, 5, 8, 11, 19, 25:30)
  prior_parameters_normal <- list(mu = 0, Sigma = diag(c(.4, .8, 1.5)))
  prior_parameters_uniform <- list(lower_bound = rep(-3, 3), upper_bound = rep(3, 3)) 
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  item_information_normal <- get_summarized_information(information_summary, estimate, model, answers, prior_form = "normal", prior_parameters = prior_parameters_normal, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, pad = TRUE)
  item_information_uniform <- get_summarized_information(information_summary, estimate, model, answers, prior_form = "uniform", prior_parameters = prior_parameters_uniform, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, pad = TRUE)
  
  expect_equal(length(item_information_normal), 50)
  expect_equal(round(item_information_normal[c(1:5, 35:40)], 3), c(8.537, 0.000, 8.252, 8.546, 0.000, 8.605, 8.995, 8.263, 8.879, 8.873, 8.934))
  expect_equal(round(item_information_uniform[c(1:5, 35:40)], 3), c(.047, .000, .042, .051, .000, .046, .050, .042, .060, .051, .057))
  expect_equal(which(round(item_information_normal, 5) == 0), c(2, 5, 8, 11, 19, 25:30))
  expect_equal(which(round(item_information_uniform, 5) == 0), c(2, 5, 8, 11, 19, 25:30))
})

context("information_summary is A (trace)")

test_that("First 1 item administered, objective is A (trace), one dimension, with padding", {
  model <- "3PLM"
  number_items <- 50
  number_dimensions <- 1
  estimator <- "maximum_aposteriori"
  estimate <- rep(0, number_dimensions)
  information_summary <- "trace"
  answers <- 1
  available <- 2:50
  administered <- 1
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  item_information_normal <- get_summarized_information(information_summary, estimate, model, answers, prior_form = "normal", prior_parameters = NULL, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, pad = TRUE)
  item_information_uniform <- get_summarized_information(information_summary, estimate, model, answers, prior_form = "uniform", prior_parameters = NULL, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, pad = TRUE)
  
  expect_equal(length(item_information_normal), 50)
  expect_equal(round(item_information_normal[c(1:5, 35:38,47:50)], 3), c(.000, .313, .139, .104, .477, .190, .236, .170, .123, .419, .142, .193, .237))
  expect_equal(item_information_uniform, item_information_normal)
})

test_that("First 1 item administered, information_summary is A (trace), three dimensions, with padding", {
  model <- "3PLM"
  number_items <- 50
  number_dimensions <- 3
  estimator <- "maximum_aposteriori"
  estimate <- rep(0, number_dimensions)
  information_summary <- "trace"
  answers <- 1
  available <- 2:50
  administered <- 1
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  item_information_normal <- get_summarized_information(information_summary, estimate, model, answers, prior_form = "normal", prior_parameters = NULL, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, pad = TRUE)
  item_information_uniform <- get_summarized_information(information_summary, estimate, model, answers, prior_form = "uniform", prior_parameters = NULL, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, pad = TRUE)
  
  expect_equal(length(item_information_normal), 50)
  expect_equal(round(item_information_normal[c(1:5, 35:38,47:50)], 3), c(.000, .488, .101, .190, .691, .266, .417, .107, .258, .684, .295, .401, .176))
  expect_equal(item_information_uniform, item_information_normal)
})

test_that("First 10 item administered, information_summary is A (trace), three dimensions, with padding", {
  model <- "3PLM"
  number_items <- 50
  number_dimensions <- 3
  estimator <- "maximum_aposteriori"
  estimate <- rep(0, number_dimensions)
  information_summary <- "trace"
  answers <- rep(c(1, 0), 5)
  available <- 11:50
  administered <- 1:10
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  item_information_normal <- get_summarized_information(information_summary, estimate, model, answers, prior_form = "normal", prior_parameters = NULL, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, pad = TRUE)
  item_information_uniform <- get_summarized_information(information_summary, estimate, model, answers, prior_form = "uniform", prior_parameters = NULL, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, pad = TRUE)
  
  expect_equal(length(item_information_normal), 50)
  expect_equal(round(item_information_normal[c(1:5, 35:38,47:50)], 3), c(.000, .000, .000, .000, .000, 3.678, 3.829, 3.520, 3.670, 4.096, 3.707, 3.813, 3.588))
  expect_equal(item_information_uniform, item_information_normal)
})

context("objective is PA (posterior trace)")

test_that("First 1 item administered, information_summary is PA (posterior trace), one dimension, with padding", {
  model <- "3PLM"
  number_items <- 50
  number_dimensions <- 1
  estimator <- "maximum_aposteriori"
  estimate <- rep(0, number_dimensions)
  information_summary <- "posterior_trace"
  answers <- 1
  available <- 2:50
  administered <- 1
  prior_parameters_normal <- list(mu = 0, Sigma = diag(1) * .4)
  prior_parameters_uniform <- list(lower_bound = -3, upper_bound = 3)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  item_information_normal <- get_summarized_information(information_summary, estimate, model, answers, prior_form = "normal", prior_parameters = prior_parameters_normal, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, pad = TRUE)
  item_information_uniform <- get_summarized_information(information_summary, estimate, model, answers, prior_form = "uniform", prior_parameters = prior_parameters_uniform, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, pad = TRUE)
  
  expect_equal(length(item_information_normal), 50)
  expect_equal(round(item_information_normal[c(1:5, 35:38,47:50)], 3), c(.000, 2.813, 2.639, 2.604, 2.977, 2.690, 2.736, 2.670, 2.623, 2.919, 2.642, 2.693, 2.737))
  expect_equal(round(item_information_uniform[c(1:5, 35:38,47:50)], 3), c(.000, .313, .139, .104, .477, .190, .236, .170, .123, .419, .142, .193, .237))
})

test_that("First 1 item administered, information_summary is PA (posterior trace), three dimensions, with padding", {
  model <- "3PLM"
  number_items <- 50
  number_dimensions <- 3
  estimator <- "maximum_aposteriori"
  estimate <- rep(0, number_dimensions)
  information_summary <- "posterior_trace"
  answers <- 1
  available <- 2:50
  administered <- 1
  prior_parameters_normal <- list(mu = 0, Sigma = diag(c(.4, .8, 1.5)))
  prior_parameters_uniform <- list(lower_bound = rep(-3, 3), upper_bound = rep(3, 3))
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  item_information_normal <- get_summarized_information(information_summary, estimate, model, answers, prior_form = "normal", prior_parameters = prior_parameters_normal, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, pad = TRUE)
  item_information_uniform <- get_summarized_information(information_summary, estimate, model, answers, prior_form = "uniform", prior_parameters = prior_parameters_uniform, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, pad = TRUE)
  
  expect_equal(length(item_information_normal), 50)
  expect_equal(round(item_information_normal[c(1:5, 35:38,47:50)], 3), c(.000, 4.904, 4.518, 4.606, 5.107, 4.683, 4.834, 4.524, 4.675, 5.101, 4.711, 4.818, 4.592))
  expect_equal(round(item_information_uniform[c(1:5, 35:38,47:50)], 3), c(.000, .488, .101, .190, .691, .266, .417, .107, .258, .684, .295, .401, .176))
})

test_that("First 10 item administered, information_summary is PA (posterior trace), three dimensions, with padding", {
  model <- "3PLM"
  number_items <- 50
  number_dimensions <- 3
  estimator <- "maximum_aposteriori"
  estimate <- rep(0, number_dimensions)
  information_summary <- "posterior_trace"
  answers <- rep(c(1, 0), 5)
  available <- 11:50
  administered <- 1:10
  prior_parameters_normal <- list(mu = 0, Sigma = diag(c(.4, .8, 1.5)))
  prior_parameters_uniform <- list(lower_bound = rep(-3, 3), upper_bound = rep(3, 3))
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  item_information_normal <- get_summarized_information(information_summary, estimate, model, answers, prior_form = "normal", prior_parameters = prior_parameters_normal, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, pad = TRUE)
  item_information_uniform <- get_summarized_information(information_summary, estimate, model, answers, prior_form = "uniform", prior_parameters = prior_parameters_uniform, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, pad = TRUE)
  
  expect_equal(length(item_information_normal), 50)
  expect_equal(round(item_information_normal[c(1:5, 35:38,47:50)], 3), c(.000, .000, .000, .000, .000, 8.095, 8.246, 7.936, 8.087, 8.513, 8.124, 8.230, 8.005))
  expect_equal(round(item_information_uniform[c(1:5, 35:38,47:50)], 3), c(.000, .000, .000, .000, .000, 3.678, 3.829, 3.520, 3.670, 4.096, 3.707, 3.813, 3.588))
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
  answers <- 1
  available <- 2:50
  administered <- 1
  prior_parameters_normal <- list(mu = 0, Sigma = diag(1) * .4)
  prior_parameters_uniform <- list(lower_bound = -3, upper_bound = 3)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  item_information_normal <- get_summarized_information(information_summary, estimate, model, answers, prior_form = "normal", prior_parameters = prior_parameters_normal, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, pad = TRUE)
  item_information_uniform <- get_summarized_information(information_summary, estimate, model, answers, prior_form = "uniform", prior_parameters = prior_parameters_uniform, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, pad = TRUE)
  
  expect_equal(length(item_information_normal), 50)
  expect_equal(round(item_information_normal[c(1:5, 35:38,47:50)], 3), c(.000, .122, .047, .025, .189, .064, .085, .057, .033, .159, .043, .066, .086))
  expect_equal(round(item_information_uniform[c(1:5, 35:38,47:50)], 3), c(.000, 3.691, 2.137, .879, 4.892, 1.949, 2.423, 1.872, 1.127, 3.799, 1.432, 2.053, 2.443))
})


context("many item characteristics equal")

test_that("information_summary is PD", {
  model <- "3PLM"
  number_items <- 50
  number_dimensions <- 1
  estimator <- "maximum_aposteriori"
  estimate <- rep(0, number_dimensions)
  information_summary <- "posterior_determinant"
  answers <- c(1, 0)
  available <- 3:50
  administered <- 1:2
  prior_parameters <- list(mu = 0, Sigma = diag(1) * .4)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(rep(1, number_items * number_dimensions), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(rep(1, number_items), nrow = number_items, ncol = 1)
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  item_information <- get_summarized_information(information_summary, estimate, model, answers, prior_form = "normal", prior_parameters = prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, pad = TRUE)
  
  expect_equal(length(item_information), 50)
  expect_equal(round(item_information[c(1:5, 35:38,47:50)], 3), c(.000, .000, 2.917, 2.917, 2.917, 2.880, 2.880, 2.880, 2.880, 2.880, 2.880, 2.880, 2.880))
})

test_that("information_summary is posterior_expected_kullback_leibler", {
  model <- "3PLM"
  number_items <- 50
  number_dimensions <- 1
  estimator <- "maximum_aposteriori"
  estimate <- rep(0, number_dimensions)
  attr(estimate, 'variance') <- .4
  information_summary <- "posterior_expected_kullback_leibler"
  answers <- c(1, 0)
  available <- 3:50
  administered <- 1:2
  prior_parameters <- list(mu = 0, Sigma = diag(1) * .4)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(rep(1, number_items * number_dimensions), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(rep(1, number_items), nrow = number_items, ncol = 1)
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  item_information <- get_summarized_information(information_summary, estimate, model, answers, prior_form = "normal", prior_parameters = prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, pad = TRUE)
  
  expect_equal(length(item_information), 50)
  expect_equal(round(item_information[c(1:5, 35:38,47:50)], 3), c(.000, .000, .022, .022, .022, .017, .017, .017, .017, .017, .017, .017, .017))
})

context("invalid input")

test_that("information_summary is of unknown type", {
  model <- "3PLM"
  number_items <- 50
  number_dimensions <- 1
  estimator <- "maximum_aposteriori"
  estimate <- rep(0, number_dimensions)
  information_summary <- "PPP"
  answers <- c(1, 0)
  available <- 3:50
  administered <- 1:2
  prior_parameters <- list(mu = 0, Sigma = diag(1) * .4)
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(rep(1, number_items * number_dimensions), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(rep(1, number_items), nrow = number_items, ncol = 1)
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  item_information <- get_summarized_information(information_summary, estimate, model, answers, prior_form = "normal", prior_parameters = prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, pad = TRUE)
    
  expect_equal(item_information$errors$information_summary, "of unknown type")
})

