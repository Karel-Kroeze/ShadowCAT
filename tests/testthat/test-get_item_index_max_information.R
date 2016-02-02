# only for whithin R:
'
library(testthat)
'

make_random_seed_exist <- rnorm(1)

test_that("with padding in objective function", {
  information_summary <- "posterior_determinant"
  estimate <- 0
  model <- "3PLM"
  estimator <- "ML"
  number_items <- 50
  number_dimensions <- 1
  available <- c(6:10, 21:30, 50)
  administered <- c(1:5, 11:20, 31:49)
  responses <- rep(c(1, 0), 17)
  prior <- .4
  lower_bound <- -3
  upper_bound <- 3
  
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  item_information <- get_item_information(information_summary, estimate, model, responses, prior, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound)
  item_index_max_information <- get_item_index_max_information(available, item_information)
  
  expect_equal(item_index_max_information, 6)
})


