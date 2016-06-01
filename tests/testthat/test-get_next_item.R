# only for whithin R:
'
library(testthat)
'

make_random_seed_exist <- rnorm(1)

context("number of answers is at least as large as number of required starting items")

model <- "3PLM"
number_items <- 50
number_dimensions <- 1
start_items <- list(type = 'random', n = 5)
estimator <- "maximum_aposteriori"
information_summary <- "posterior_determinant"
lp_constraints <- NULL
lp_characters <- NULL
estimate <- 0
attr(estimate, "variance") <- diag(1)
prior_form = "normal"
prior_parameters <- list(mu = 0, Sigma = diag(1) * .4)
stop_test <- list(max_n = number_items)

guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)

test_that("number of answers is 5, number of required starting items 5", {
  answers <- rep(1, 5)
  available <- c(3:10, 14:50)
  administered <- c(1, 3, 11, 12, 13)
  
  item_next <- get_next_item(start_items, information_summary, lp_constraints, lp_characters, estimate, model, answers, prior_form, prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, stop_test)
  expect_equal(item_next, 5)
})

test_that("number of answers is 7, number of required starting items 5", {
  answers <- rep(1, 7)
  available <- c(3:4, 7:10, 14:50)
  administered <- c(1, 3, 5:6, 11, 12, 13)
  
  item_next <- get_next_item(start_items, information_summary, lp_constraints, lp_characters, estimate, model, answers, prior_form, prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, stop_test)
  expect_equal(item_next, 47)
})

context("start type is random")

test_that("number of answers is 0, number of required starting items 5, start type random", {
  answers <- numeric(0)
  available <- 1:50
  administered <- numeric(0)
  
  item_next <- with_random_seed(2, get_next_item)(start_items, information_summary, lp_constraints, lp_characters, estimate, model, answers, prior_form, prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, stop_test)
  expect_equal(item_next, 10)
})

test_that("number of answers is 2, number of required starting items 5, start type random", {
  answers <- 1:2
  available <- 3:50
  administered <- 1:2
  
  item_next <- with_random_seed(2, get_next_item)(start_items, information_summary, lp_constraints, lp_characters, estimate, model, answers, prior_form, prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, stop_test)
  expect_equal(item_next, 11)
})


context("start type is fixed")

start_items <- list(type = 'fixed', indices = c(2, 4, 7:9), n = 5)

test_that("number of answers is 0, number of required starting items 5, start type fixed", {
  answers <- numeric(0)
  available <- 1:50
  administered <- numeric(0)
  
  item_next <- get_next_item(start_items, information_summary, lp_constraints, lp_characters, estimate, model, answers, prior_form, prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, stop_test)
  expect_equal(item_next, 2)
})

test_that("number of answers is 4, number of required starting items 5, start type fixed", {
  answers <- c(1, 0, 1, 1)
  available <- c(1, 3, 5:6, 9:50)
  administered <- c(2, 4, 7:8)
  item_next <- get_next_item(start_items, information_summary, lp_constraints, lp_characters, estimate, model, answers, prior_form, prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, stop_test)
  expect_equal(item_next, 9)
})

context("start type is random by dimension, items ordered per dimension")

model <- "3PLM"
number_items <- 50
number_dimensions <- 3
start_items <- list(type = 'random_by_dimension', n_by_dimension = 3, n = 9) 
estimator <- "maximum_aposteriori"
information_summary <- "posterior_determinant"
lp_constraints <- NULL
lp_characters <- NULL
estimate <- 0
attr(estimate, "variance") <- diag(1)
prior_form <- "normal"
prior_parameters <- list(mu = 0, Sigma = diag(1) * .4)
stop_test <- list(max_n = number_items)

guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
alpha[1:15,2:3] <- 0
alpha[16:30,c(1,3)] <- 0
alpha[31:50,1:2] <- 0
beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)

test_that("number of answers is 0, number of required starting items 3 per dimension, random by dimension", {
  answers <- numeric(0)
  available <- 1:50
  administered <- numeric(0)
  
  item_next <- with_random_seed(2, get_next_item)(start_items, information_summary, lp_constraints, lp_characters, estimate, model, answers, prior_form, prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, stop_test)
  expect_equal(item_next, 3)
})

test_that("number of answers is 6, number of required starting items 3 per dimension, random by dimension", {
  available <- c(2, 4:5, 7:15, 17:18, 20:26, 28:50)
  answers <- rep(1, 6)
  administered <- c(1, 3, 6, 16, 19, 27)
  
  item_next <- with_random_seed(2, get_next_item)(start_items, information_summary, lp_constraints, lp_characters, estimate, model, answers, prior_form, prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, stop_test)
  expect_equal(item_next, 34)
})

test_that("number of answers is 0, number of required starting items for each dimension is 2, 4, 2, random by dimension", {
  start_items <- list(type = 'random_by_dimension', n_by_dimension = c(2, 4, 2), n = 8)
  available <- 1:50
  answers <- numeric(0)
  administered <- numeric(0)
  
  item_next <- with_random_seed(2, get_next_item)(start_items, information_summary, lp_constraints, lp_characters, estimate, model, answers, prior_form, prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, stop_test)
  expect_equal(item_next, 3)
})

test_that("number of answers is 6, number of required starting items for each dimension is 2, 5, 1, random by dimension", {
  start_items <- list(type = 'random_by_dimension', n_by_dimension = c(2, 5, 1), n = 8)
  available <- c(2, 4:5, 7:15, 17:18, 20:26, 28:50)
  answers <- rep(1, 6)
  administered <- c(1, 3, 6, 16, 19, 27)
  
  item_next <- with_random_seed(2, get_next_item)(start_items, information_summary, lp_constraints, lp_characters, estimate, model, answers, prior_form, prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, stop_test)
  expect_equal(item_next, 20)
})


test_that("items load on all three dimension while start type is n_by_dimension", {
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  available <- c(2, 4:5, 7:15, 17:18, 20:26, 28:50)
  administered <- c(1, 3, 6, 16, 19, 27)
  answers <- rep(1, 6)
  
  item_next <- with_random_seed(2, get_next_item)(start_items, information_summary, lp_constraints, lp_characters, estimate, model, answers, prior_form, prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, stop_test)
  
  expect_equal(item_next, 12)
})

context("start type is random by dimension, items not ordered per dimension")

model <- "3PLM"
number_items <- 50
number_dimensions <- 3
start_items <- list(type = 'random_by_dimension', n_by_dimension = 3, n = 9) 
estimator <- "maximum_aposteriori"
information_summary <- "posterior_determinant"
lp_constraints <- NULL
lp_characters <- NULL
estimate <- 0
attr(estimate, "variance") <- diag(1)
prior_form <- "normal"
prior_parameters <- list(mu = 0, Sigma = diag(1) * .4)
stop_test <- list(max_n = number_items)

guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
indices_shuffled <- with_random_seed(2, sample)(1:50, 50)
alpha[indices_shuffled[1:10],2:3] <- 0
alpha[indices_shuffled[11:20],c(1,3)] <- 0
alpha[indices_shuffled[21:50],1:2] <- 0
beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)

test_that("number of answers is 0, number of required starting items 3 per dimension, random by dimension", {
  answers <- numeric(0)
  available <- 1:50
  administered <- numeric(0)
  
  item_next <- with_random_seed(2, get_next_item)(start_items, information_summary, lp_constraints, lp_characters, estimate, model, answers, prior_form, prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, stop_test)
  expect_equal(item_next, 8)
})

test_that("number of answers is 6, number of required starting items 3 per dimension, random by dimension", {
  available <- c(2, 4:5, 7:15, 17:18, 20:26, 28:50)
  answers <- rep(1, 6)
  administered <- c(1, 3, 6, 16, 19, 27)
  
  item_next <- with_random_seed(2, get_next_item)(start_items, information_summary, lp_constraints, lp_characters, estimate, model, answers, prior_form, prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, stop_test)
  expect_equal(item_next, 11)
})

test_that("number of answers is 0, number of required starting items for each dimension is 2, 4, 2, random by dimension", {
  start_items <- list(type = 'random_by_dimension', n_by_dimension = c(2, 4, 2), n = 8)
  available <- 1:50
  answers <- numeric(0)
  administered <- numeric(0)
  
  item_next <- with_random_seed(2, get_next_item)(start_items, information_summary, lp_constraints, lp_characters, estimate, model, answers, prior_form, prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, stop_test)
  expect_equal(item_next, 8)
})

test_that("number of answers is 6, number of required starting items for each dimension is 2, 5, 1, random by dimension", {
  start_items <- list(type = 'random_by_dimension', n_by_dimension = c(2, 5, 1), n = 8)
  available <- c(2, 4:5, 7:15, 17:18, 20:26, 28:50)
  answers <- rep(1, 6)
  administered <- c(1, 3, 6, 16, 19, 27)
  
  item_next <- with_random_seed(2, get_next_item)(start_items, information_summary, lp_constraints, lp_characters, estimate, model, answers, prior_form, prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, stop_test)
  expect_equal(item_next, 15)
})

test_that("number of answers is 6, number of required starting items for each dimension is 2, 5, 1, random by dimension", {
  start_items <- list(type = 'random_by_dimension', n_by_dimension = c(2, 5, 1), n = 8)
  available <- c(2, 4:5, 7:14, 17:18, 20:26, 28:50)
  answers <- rep(1, 7)
  administered <- c(1, 3, 6, 16, 19, 27, 15)
  
  item_next <- with_random_seed(2, get_next_item)(start_items, information_summary, lp_constraints, lp_characters, estimate, model, answers, prior_form, prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, stop_test)
  expect_equal(item_next, 11)
})


context("invalid input")

available <- c(2, 4:5, 7:15, 17:18, 20:26, 28:50)
administered <- c(1, 3, 6, 16, 19, 27)
answers <- rep(1, 6)

test_that("sum n_by_dimension is not equal to n, same n for each dimension", {
  start_items <- list(type = 'random_by_dimension', n_by_dimension = 2, n = 9)
  
  item_next <- get_next_item(start_items, information_summary, lp_constraints, lp_characters, estimate, model, answers, prior_form, prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, stop_test)
  expect_equal(item_next$errors$start, "contains inconsistent information. Total length of start phase and sum of length per dimension do not match")
})

test_that("sum n_by_dimension is not equal to n, different n for each dimension", {
  start_items <- list(type = 'random_by_dimension', n_by_dimension = c(2, 5, 1), n = 9)
                         
  item_next <- get_next_item(start_items, information_summary, lp_constraints, lp_characters, estimate, model, answers, prior_form, prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, stop_test)
  expect_equal(item_next$errors$start, "contains inconsistent information. Total length of start phase and sum of length per dimension do not match (n != sum(n_by_dimension)")
})

