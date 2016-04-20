# only for whithin R:
'
library(testthat)
library(lpSolve)
'

make_random_seed_exist <- rnorm(1)

model <- '3PLM'
number_items <- 50
number_dimensions <- 1
number_items <- 50
max_n <- 20
estimator <- "maximum_aposteriori"
information_summary <- "posterior_determinant"
estimate <- 0
prior_form <- "normal"
prior_parameters <- list(mu = 0, Sigma = matrix(.4))

guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)

#create item characteristics and constraints
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

chars_constraints_lp <- constraints_lp_format(max_n, number_items, characteristics, constraints)

test_that("First 4 items administered", {
  responses <- rep(c(1, 0), 2)
  available <- c(5:50)
  administered <- c(1:4)

  item_information <- get_summarized_information(information_summary, estimate, model, responses, prior_form, prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, pad = TRUE) 
  best_item <- get_item_index_max_information_constrained(number_items, administered, available, responses, chars_constraints_lp$lp_constraints, chars_constraints_lp$lp_chars, item_information)

  expect_equal(best_item, 5)
})

test_that("7 items administered", {
  responses <- c(1, 0, 1, 0, 1, 0, 1)
  available <- c(2, 4, 6, 8, 10, 12, 14:50)
  administered <- c(1, 3, 5, 7, 9, 11, 13)
  
  item_information <- get_summarized_information(information_summary, estimate, model, responses, prior_form, prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, pad = TRUE) 
  best_item <- get_item_index_max_information_constrained(number_items, administered, available, responses, chars_constraints_lp$lp_constraints, chars_constraints_lp$lp_chars, item_information)
  
  expect_equal(best_item, 6)
})

test_that("First 10 items administered", {
  responses <- rep(1, 10)
  available <- c(11:50)
  administered <- c(1:10)
  
  item_information <- get_summarized_information(information_summary, estimate, model, responses, prior_form, prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, pad = TRUE) 
  best_item <- get_item_index_max_information_constrained(number_items, administered, available, responses, chars_constraints_lp$lp_constraints, chars_constraints_lp$lp_chars, item_information)
  
  expect_equal(best_item, 47)
})

test_that("None administered", {
  responses <- numeric(0)
  available <- 1:50
  administered <- numeric(0)
  
  item_information <- get_summarized_information(information_summary, estimate, model, responses, prior_form, prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, pad = TRUE) 
  best_item <- get_item_index_max_information_constrained(number_items, administered, available, responses, chars_constraints_lp$lp_constraints, chars_constraints_lp$lp_chars, item_information)
  
  expect_equal(best_item, 5)
})


