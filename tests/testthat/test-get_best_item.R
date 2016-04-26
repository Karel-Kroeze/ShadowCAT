# only for whithin R:
'
library(testthat)
library(lpSolve)
'

make_random_seed_exist <- rnorm(1)

model <- "3PLM"
information_summary <- "posterior_determinant"
estimator <- "maximum_aposteriori"
max_n <- 20
stop_test <- list(max_n = max_n)
number_items <- 50
number_dimensions <- 1
estimate <- 0
attr(estimate, "variance") <- diag(1)
prior_form = "normal"
prior_parameters = list(mu = 0, Sigma = diag(1) * .4)
guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)

characteristics <- data.frame(content = with_random_seed(2, sample)(c('algebra','physics','calculus'), number_items, TRUE),
                              time = with_random_seed(2, rnorm)(number_items),
                              exclusive = sapply(1:number_items, FUN = function (x) { if (x %in% with_random_seed(2, sample)(number_items, 4)) 1 else 0 }))
constraints <- list(list(name = 'content/algebra',
                         op = '><',
                         target = c(5, 10)),
                    list(name = 'content/physics',
                         op = '><',
                         target = c(2, 5)),
                    list(name = 'time',
                         op = '<',
                         target = 20),
                    list(name = 'exclusive',
                         op = '<',
                         target = 2))

constraints_formatted <- constraints_lp_format(max_n = max_n, number_items = number_items, characteristics = characteristics, constraints = constraints)

test_that("First 4 items administered, MI", {
  answers <- rep(c(1, 0), 2)
  available <- c(5:50)
  administered <- c(1:4)
  constraints_formatted <- NULL
  
  best_item <- get_best_item(information_summary, constraints_formatted$lp_constraints, constraints_formatted$lp_chars, estimate, model, answers, prior_form, prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, stop_test)

  expect_equal(best_item, 5)
})

test_that("First 4 items administered, Shadow", {
  answers = rep(c(1, 0), 2)
  available <- c(5:50)
  administered <- c(1:4)
  
  best_item <- get_best_item(information_summary, constraints_formatted$lp_constraints, constraints_formatted$lp_chars, estimate, model, answers, prior_form, prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, stop_test)
  
  expect_equal(best_item, 5)
})

test_that("7 items administered, MI", {
  answers <- c(1, 0, 1, 0, 1, 0, 1)
  available <- c(2, 4, 6, 8, 10, 12, 14:50)
  administered <- c(1, 3, 5, 7, 9, 11, 13)
  constraints_formatted <- NULL
  
  best_item <- get_best_item(information_summary, constraints_formatted$lp_constraints, constraints_formatted$lp_chars, estimate, model, answers, prior_form, prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, stop_test)
  
  expect_equal(best_item, 6)
})

test_that("7 items administered, Shadow", {
  answers <- c(1, 0, 1, 0, 1, 0, 1)
  available <- c(2, 4, 6, 8, 10, 12, 14:50)
  administered <- c(1, 3, 5, 7, 9, 11, 13)
  
  best_item <- get_best_item(information_summary, constraints_formatted$lp_constraints, constraints_formatted$lp_chars, estimate, model, answers, prior_form, prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, stop_test)
  
  expect_equal(best_item, 6)
})

test_that("First 10 items administered, MI", {
  answers <- rep(1, 10)
  available <- c(11:50)
  administered <- c(1:10)
  
  best_item <- get_best_item(information_summary, constraints_formatted$lp_constraints, constraints_formatted$lp_chars, estimate, model, answers, prior_form, prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, stop_test)
  
  expect_equal(best_item, 47)
})

test_that("First 10 items administered, Shadow", {
  answers <- rep(1, 10)
  available <- c(11:50)
  administered <- c(1:10)
  
  best_item <- get_best_item(information_summary, constraints_formatted$lp_constraints, constraints_formatted$lp_chars, estimate, model, answers, prior_form, prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, stop_test)
  
  expect_equal(best_item, 47)
})

context("many item characteristics equal")

model <- "3PLM"
information_summary <- "posterior_determinant"
estimator <- "maximum_aposteriori"
max_n <- 20
stop_test <- list(max_n = max_n)
number_items <- 50
number_dimensions <- 1
estimate <- 0
attr(estimate, "variance") <- diag(1)
prior_form = "normal"
prior_parameters <- list(mu = 0, Sigma = diag(1) * .4)
guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
alpha <- matrix(rep(1, number_items * number_dimensions), nrow = number_items, ncol = number_dimensions)
beta <- matrix(rep(1, number_items), nrow = number_items, ncol = 1)
number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)

characteristics <- data.frame(content = with_random_seed(2, sample)(c('algebra','physics','calculus'), number_items, TRUE),
                              time = with_random_seed(2, rnorm)(number_items),
                              exclusive = sapply(1:number_items, FUN = function (x) { if (x %in% with_random_seed(2, sample)(number_items, 4)) 1 else 0 }))
constraints <- list(list(name = 'content/algebra',
                         op = '><',
                         target = c(5, 10)),
                    list(name = 'content/physics',
                         op = '><',
                         target = c(2, 5)),
                    list(name = 'time',
                         op = '<',
                         target = 20),
                    list(name = 'exclusive',
                         op = '<',
                         target = 2))
constraints_formatted <- constraints_lp_format(max_n = max_n, number_items = number_items, characteristics = characteristics, constraints = constraints)


test_that("First 10 items administered, MI", {
  answers <- rep(1, 10)
  available <- c(11:50)
  administered <- c(1:10)
  constraints_formatted <- NULL
  
  best_item <- with_random_seed(2, get_best_item)(information_summary, constraints_formatted$lp_constraints, constraints_formatted$lp_chars, estimate, model, answers, prior_form, prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, stop_test)
  
  expect_equal(best_item, 13)
})

test_that("First 10 items administered, Shadow", {
  answers <- rep(1, 10)
  available <- c(11:50)
  administered <- c(1:10)
  
  best_item <- with_random_seed(2, get_best_item)(information_summary, constraints_formatted$lp_constraints, constraints_formatted$lp_chars, estimate, model, answers, prior_form, prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, stop_test)
  
  expect_equal(best_item, 13)
})


