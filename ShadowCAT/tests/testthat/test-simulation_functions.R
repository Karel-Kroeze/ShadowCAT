# only for whithin R:
'
library(testthat)
'

make_random_seed_exist <- rnorm(1)

context("test simulate_testbank")

test_that("one dimension, one itemstep", {
  alpha_beta <- with_random_seed(1, simulate_testbank)(model = "GPCM", number_items = 50, number_dimensions = 1, number_itemsteps = 1)
  expect_equal(dim(alpha_beta$alpha), c(50, 1))
  expect_equal(rownames(alpha_beta$alpha), str_c("item", 1:50))
  expect_equal(round(alpha_beta$alpha[4], 3), 1.39)
  expect_equal(dim(alpha_beta$beta), c(50, 1))
  expect_equal(rownames(alpha_beta$beta), str_c("item", 1:50))
  expect_equal(round(alpha_beta$beta[4], 3), -.478)
  expect_equal(round(fivenum(alpha_beta$beta), 3), c(-1.805, -0.478, -0.047, 0.611, 2.402))
  })

test_that("two dimensions, three itemsteps", {
  alpha_beta <- with_random_seed(1, simulate_testbank)(model = "GPCM", number_items = 50, number_dimensions = 2, number_itemsteps = 3)
  expect_equal(dim(alpha_beta$alpha), c(50, 2))
  expect_equal(rownames(alpha_beta$alpha), str_c("item", 1:50))
  expect_equal(round(alpha_beta$alpha[4,], 3), c(1.390, .594))
  expect_equal(dim(alpha_beta$beta), c(50, 3))
  expect_equal(rownames(alpha_beta$beta), str_c("item", 1:50))
  expect_equal(round(alpha_beta$beta[4,], 3), c(-3.129, -4.259, -3.388))
  expect_equal(round(alpha_beta$beta[18,], 3), c(-0.534, 0.931, 4.397))
  expect_equal(round(fivenum(alpha_beta$beta), 3), c(-5.610, -2.608, -1.563, 0.126, 7.205))
})

test_that("items load on one dimension and with varying number of item steps", {
  alpha_beta <- with_random_seed(1, simulate_testbank)(model = "GPCM", number_items = 51, number_dimensions = 2, number_itemsteps = 3, items_load_one_dimension = TRUE, varying_number_item_steps = TRUE)
  expect_equal(dim(alpha_beta$alpha), c(51, 2))
  expect_equal(rownames(alpha_beta$alpha), str_c("item", 1:51))
  expect_equal(round(alpha_beta$alpha[26,], 3), c(.763, 0))
  expect_equal(round(alpha_beta$alpha[27,], 3), c(0, .768))
  expect_equal(unname(alpha_beta$alpha[1:26, 2]), rep(0, 26))
  expect_equal(unname(alpha_beta$alpha[27:51, 1]), rep(0, 25))
  expect_equal(dim(alpha_beta$beta), c(51, 3))
  expect_equal(rownames(alpha_beta$beta), str_c("item", 1:51))
  expect_equal(round(alpha_beta$beta[4,], 3), c(-0.567, 0.866, 4.299))
  expect_equal(round(alpha_beta$beta[6,], 3), c(-2.367, NA, NA))
})

context("test simulate_answer")

test_that("one dimension, one itemstep, one answer", {
  alpha_beta <- with_random_seed(1, simulate_testbank)(model = "GPCM", number_items = 50, number_dimensions = 1, number_itemsteps = 1)
  answer <- with_random_seed(1, simulate_answer)(theta = 2, model = "GPCM", alpha = alpha_beta$alpha, beta = alpha_beta$beta, guessing = NULL, item_keys = "item5")
  expect_equal(answer, 0)
})

test_that("three dimensions, four itemsteps, one_answer", {
  number_items <- 50
  alpha_beta <- with_random_seed(1, simulate_testbank)(model = "GPCM", number_items = number_items, number_dimensions = 3, number_itemsteps = 4)
  answer <- with_random_seed(1, simulate_answer)(theta = c(1, 0, 2), model = "GPCM", alpha = alpha_beta$alpha, beta = alpha_beta$beta, guessing = NULL, item_keys = "item5")
  expect_equal(answer, 4)
})

test_that("three dimensions, four itemsteps, five answers", {
  alpha_beta <- with_random_seed(1, simulate_testbank)(model = "GPCM", number_items = 50, number_dimensions = 3, number_itemsteps = 4)
  answer <- with_random_seed(1, simulate_answer)(theta = c(-2, 2, 3), model = "GPCM", alpha = alpha_beta$alpha, beta = alpha_beta$beta, guessing = NULL, item_keys = c("item5", "item2", "item8", "item1", "item18"))
  expect_equal(answer, c(4, 4, 3, 4, 1))
})

test_that("with guessing", {
  number_items <- 50
  alpha_beta <- with_random_seed(1, simulate_testbank)(model = "3PLM", number_items = number_items, number_dimensions = 1, number_itemsteps = 1)
  guessing <- matrix(rep(.25, number_items), dimnames = list(str_c("item", 1:number_items)))
  answer <- with_random_seed(1, simulate_answer)(theta = 2, model = "3PLM", alpha = alpha_beta$alpha, beta = alpha_beta$beta, guessing = guessing, item_keys = "item5")
  expect_equal(answer, 1)
})

test_that("invalid input", {
  number_items <- 50
  alpha_beta <- with_random_seed(1, simulate_testbank)(model = "3PLM", number_items = number_items, number_dimensions = 1, number_itemsteps = 1)
  guessing <- matrix(rep(.25, number_items), dimnames = list(str_c("i", 1:number_items)))
  error_messages <- simulate_answer(theta = c(1, 1, 1), model = "PLM", alpha = alpha_beta$alpha, beta = alpha_beta$beta, guessing = guessing, item_keys = "ite5")
  expect_equal(error_messages$errors$theta, "should have length equal to the number of columns of alpha")
  expect_equal(error_messages$errors$model, "of unknown type")
  expect_equal(error_messages$errors$alpha_beta_eta_guessing, "should have equal row names, in same order")
  expect_equal(error_messages$errors$item_keys, "unknown")
  
  error_message_alpha <- simulate_answer(theta = c(1, 1, 1), model = "PLM", alpha = .1, beta = alpha_beta$beta, guessing = guessing, item_keys = "ite5")
  error_message_beta <- simulate_answer(theta = c(1, 1, 1), model = "PLM", alpha = alpha_beta$alpha, beta = .1, guessing = guessing, item_keys = "ite5")
  expect_equal(error_message_alpha$errors$alpha, "should be a matrix with item keys as row names")
  expect_equal(error_message_beta$errors$beta, "should be a matrix with item keys as row names")
})

context("test test_shadowcat")

test_that("invalid input", {
  # define true theta for simulation of answers
  true_theta <- c(1, 0, 2)
  
  # define item characteristics
  number_items <- 300
  number_dimensions <- 3
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  item_keys <- str_c("item", 1:number_items)
  model <- '3PLM'
  
  alpha_beta <- simulate_testbank(model = model, number_items = number_items, number_dimensions = number_dimensions, number_itemsteps = number_answer_categories - 1)
  
  start_items <- list(type = 'random', n = 3)
  stop_test <- list(max_n = 300)
  estimator <- 'maximum_aposteriori'
  information_summary <- 'posterior_determinant'
  
  # define prior
  prior_form <- "uniform"  
  prior_parameters <- list(lower_bound = c(-1, -1, -1), upper_bound = c(1, 1, 1))
  
  error_message_true_theta <- test_shadowcat(true_theta = 2, prior_form = prior_form, prior_parameters = prior_parameters, model = model, alpha = alpha_beta$alpha, beta = alpha_beta$beta, start_items = start_items, stop_test = stop_test, estimator = estimator, information_summary = information_summary)
  expect_equal(error_message_true_theta$errors$theta, "should have length equal to the number of columns of alpha")
  
  error_message_true_theta <- test_shadowcat(true_theta = true_theta, prior_form = "exponential", prior_parameters = prior_parameters, model = model, alpha = alpha_beta$alpha, beta = alpha_beta$beta, start_items = start_items, stop_test = stop_test, estimator = estimator, information_summary = information_summary)
  expect_equal(error_message_true_theta$errors$prior_form, "of unknown type")
})


