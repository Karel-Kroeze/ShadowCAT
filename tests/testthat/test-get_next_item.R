# only for whithin R:
'
library(testthat)
'

make_random_seed_exist <- rnorm(1)

context("number of responses is at least as large as number of required starting items")

model <- "3PLM"
number_items <- 50
number_dimensions <- 1
start_items <- list(type = 'random', n = 5)
estimator <- "maximum_aposteriori"
information_summary <- "posterior_determinant"
lp_constraints <- NULL
lp_characters <- NULL
estimate <- 0
prior <- .4
lower_bound <- -3
upper_bound <- 3

guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)

test_that("number of responses is 5, number of required starting items 5", {
  responses <- rep(1, 5)
  available <- c(3:10, 14:50)
  administered <- c(1, 3, 11, 12, 13)
  
  item_next <- get_next_item(start_items, information_summary, lp_constraints, lp_characters, estimate, model, responses, prior, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound)
  expect_equal(item_next, 5)
})

test_that("number of responses is 7, number of required starting items 5", {
  responses <- rep(1, 7)
  available <- c(3:4, 7:10, 14:50)
  administered <- c(1, 3, 5:6, 11, 12, 13)
  
  item_next <- get_next_item(start_items, information_summary, lp_constraints, lp_characters, estimate, model, responses, prior, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound)
  expect_equal(item_next, 47)
})

context("start type is random")

test_that("number of responses is 0, number of required starting items 5, start type random", {
  responses <- numeric(0)
  available <- 1:50
  administered <- numeric(0)
  
  item_next <- with_random_seed(2, get_next_item)(start_items, information_summary, lp_constraints, lp_characters, estimate, model, responses, prior, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound)
  expect_equal(item_next, 10)
})

test_that("number of responses is 2, number of required starting items 5, start type random", {
  responses <- 1:2
  available <- 3:50
  administered <- 1:2
  
  item_next <- with_random_seed(2, get_next_item)(start_items, information_summary, lp_constraints, lp_characters, estimate, model, responses, prior, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound)
  expect_equal(item_next, 11)
})


context("start type is fixed")

start_items <- list(type = 'fixed', indeces = c(2, 4, 7:9), n = 5)

test_that("number of responses is 0, number of required starting items 5, start type fixed", {
  responses <- numeric(0)
  available <- 1:50
  administered <- numeric(0)
  
  item_next <- get_next_item(start_items, information_summary, lp_constraints, lp_characters, estimate, model, responses, prior, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound)
  expect_equal(item_next, 2)
})

test_that("number of responses is 4, number of required starting items 5, start type fixed", {
  responses <- c(1, 0, 1, 1)
  available <- c(1, 3, 5:6, 9:50)
  administered <- c(2, 4, 7:8)
  item_next <- get_next_item(start_items, information_summary, lp_constraints, lp_characters, estimate, model, responses, prior, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound)
  expect_equal(item_next, 9)
})

context("start type is random by dimension")

model <- "3PLM"
number_items <- 50
number_dimensions <- 3
start_items <- list(type = 'random_by_dimension', n_by_dimension = 3, n = 9) 
estimator <- "maximum_aposteriori"
information_summary <- "posterior_determinant"
lp_constraints <- NULL
lp_characters <- NULL
estimate <- 0
prior <- .4
lower_bound <- -3
upper_bound <- 3

guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
alpha[1:15,2:3] <- 0
alpha[16:30,c(1,3)] <- 0
alpha[31:50,1:2] <- 0
beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)

test_that("number of responses is 0, number of required starting items 3 per dimension, random by dimension", {
  responses <- numeric(0)
  available <- 1:50
  administered <- numeric(0)
  
  item_next <- with_random_seed(2, get_next_item)(start_items, information_summary, lp_constraints, lp_characters, estimate, model, responses, prior, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound)
  expect_equal(item_next, 3)
})

test_that("number of responses is 6, number of required starting items 3 per dimension, random by dimension", {
  available <- c(2, 4:5, 7:15, 17:18, 20:26, 28:50)
  responses <- rep(1, 6)
  administered <- c(1, 3, 6, 16, 19, 27)
  
  item_next <- with_random_seed(2, get_next_item)(start_items, information_summary, lp_constraints, lp_characters, estimate, model, responses, prior, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound)
  expect_equal(item_next, 34)
})

test_that("number of responses is 0, number of required starting items for each dimension is 2, 4, 2, random by dimension", {
  available <- 1:50
  responses <- numeric(0)
  administered <- numeric(0)
  
  item_next <- with_random_seed(2, get_next_item)(start_items, information_summary, lp_constraints, lp_characters, estimate, model, responses, prior, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound)
  expect_equal(item_next, 3)
})

test_that("number of responses is 6, number of required starting items for each dimension is 2, 5, 1, random by dimension", {
  start_items <- list(type = 'random_by_dimension', n_by_dimension = c(2, 5, 1), n = 8)
  
  available <- c(2, 4:5, 7:15, 17:18, 20:26, 28:50)
  responses <- rep(1, 6)
  administered <- c(1, 3, 6, 16, 19, 27)
  
  item_next <- with_random_seed(2, get_next_item)(start_items, information_summary, lp_constraints, lp_characters, estimate, model, responses, prior, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound)
  expect_equal(item_next, 20)
})


test_that("items load on all three dimension while start type is n_by_dimension", {
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  available <- c(2, 4:5, 7:15, 17:18, 20:26, 28:50)
  administered <- c(1, 3, 6, 16, 19, 27)
  responses <- rep(1, 6)
  
  item_next <- with_random_seed(2, get_next_item)(start_items, information_summary, lp_constraints, lp_characters, estimate, model, responses, prior, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound)
  
  expect_equal(item_next, 12)
})

context("invalid input")

available <- c(2, 4:5, 7:15, 17:18, 20:26, 28:50)
administered <- c(1, 3, 6, 16, 19, 27)
responses <- rep(1, 6)

test_that("sum n_by_dimension is not equal to n, same n for each dimension", {
  start_items <- list(type = 'random_by_dimension', n_by_dimension = 2, n = 9)
  
  item_next <- get_next_item(start_items, information_summary, lp_constraints, lp_characters, estimate, model, responses, prior, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound)
  expect_equal(item_next$errors$start, "contains inconsistent information. Total length of start phase and sum of length per dimension do not match")
})

test_that("sum n_by_dimension is not equal to n, different n for each dimension", {
  start_items <- list(type = 'random_by_dimension', n_by_dimension = c(2, 5, 1), n = 9)
                         
  item_next <- get_next_item(start_items, information_summary, lp_constraints, lp_characters, estimate, model, responses, prior, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound)
  expect_equal(item_next$errors$start, "contains inconsistent information. Total length of start phase and sum of length per dimension do not match (n != sum(n_by_dimension)")
})


