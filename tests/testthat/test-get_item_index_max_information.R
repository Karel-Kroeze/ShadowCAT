# only for whithin R:
'
library(testthat)
'

make_random_seed_exist <- rnorm(1)

test_that("one dimension", {
  information_summary <- "posterior_determinant"
  estimate <- 0
  attr(estimate, "variance") <- matrix(.2^2)
  model <- "3PLM"
  estimator <- "maximum_aposteriori"
  number_items <- 50
  number_dimensions <- 1
  available <- c(6:10, 21:30, 50)
  administered <- c(1:5, 11:20, 31:49)
  responses <- rep(c(1, 0), 17)
  
  stop_test <- list(max_n = number_items, target = .1^2)
  
  prior_form <- "normal"
  prior_parameters <- list(mu = 0, Sigma = matrix(.4))
  
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  item_information <- get_summarized_information(information_summary, estimate, model, responses, prior_form, prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item)
  item_index_max_information <- get_item_index_max_information(available, item_information, estimate, stop_test, alpha, number_answers = length(administered))
  
  test_item_index_max_information <- validate_and_run.test(fn = get_item_index_max_information, available = available, item_information = item_information, estimate = estimate, stop_test = stop_test, alpha = alpha, number_answers = length(administered))
  helper_get_uncompleted_dimensions <- get('get_uncompleted_dimensions', environment(test_item_index_max_information))
  helper_get_useful_item_indices <- get('get_useful_item_indices', environment(test_item_index_max_information))
  
  expect_equal(helper_get_uncompleted_dimensions(), NULL)
  expect_equal(helper_get_useful_item_indices(NULL), 1:nrow(alpha))
  expect_equal(item_index_max_information, 6)
})

test_that("three dimensions, items load all dimensions, no variance target", {
  information_summary <- "posterior_determinant"
  estimate <- c(.2, .1, 2.5)
  attr(estimate, "variance") <- diag(c(.3, .1, .2))
  model <- "3PLM"
  estimator <- "maximum_aposteriori"
  number_items <- 50
  number_dimensions <- 3
  available <- c(6:10, 21:30, 50)
  administered <- c(31:49, 5:1, 11:20)
  responses <- rep(c(1, 0), 17)
  
  stop_test <- list(max_n = number_items)
  
  prior_form <- "normal"
  prior_parameters <- list(mu = rep(0, 3), Sigma = diag(number_dimensions) * .4)
  
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  item_information <- get_summarized_information(information_summary, estimate, model, responses, prior_form, prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item)
  item_index_max_information <- get_item_index_max_information(available, item_information, estimate, stop_test, alpha, number_answers = length(administered))
  
  test_item_index_max_information <- validate_and_run.test(fn = get_item_index_max_information, available = available, item_information = item_information, estimate = estimate, stop_test = stop_test, alpha = alpha, number_answers = length(administered))
  helper_get_uncompleted_dimensions <- get('get_uncompleted_dimensions', environment(test_item_index_max_information))
  helper_get_useful_item_indices <- get('get_useful_item_indices', environment(test_item_index_max_information))
  
  expect_equal(helper_get_uncompleted_dimensions(), NULL)
  expect_equal(helper_get_useful_item_indices(NULL), 1:nrow(alpha))
  expect_equal(item_index_max_information, 29)
})

test_that("three dimensions, items load all dimensions, with variance target", {
  information_summary <- "posterior_determinant"
  estimate <- c(.2, .1, 2.5)
  attr(estimate, "variance") <- diag(c(.3, .1, .2))
  model <- "3PLM"
  estimator <- "maximum_aposteriori"
  number_items <- 50
  number_dimensions <- 3
  available <- c(6:10, 21:30, 50)
  administered <- c(31:49, 5:1, 11:20)
  responses <- rep(c(1, 0), 17)
  
  stop_test <- list(max_n = number_items, target = c(.2, .2, .2))
  
  prior_form <- "normal"
  prior_parameters <- list(mu = rep(0, 3), Sigma = diag(number_dimensions) * .4)
  
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  item_information <- get_summarized_information(information_summary, estimate, model, responses, prior_form, prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item)
  item_index_max_information <- get_item_index_max_information(available, item_information, estimate, stop_test, alpha, number_answers = length(administered))
  
  test_item_index_max_information <- validate_and_run.test(fn = get_item_index_max_information, available = available, item_information = item_information, estimate = estimate, stop_test = stop_test, alpha = alpha, number_answers = length(administered))
  helper_get_uncompleted_dimensions <- get('get_uncompleted_dimensions', environment(test_item_index_max_information))
  helper_get_useful_item_indices <- get('get_useful_item_indices', environment(test_item_index_max_information))
  
  expect_equal(helper_get_uncompleted_dimensions(), c(1, 3))
  expect_equal(helper_get_useful_item_indices(c(1, 3)), 1:nrow(alpha))
  expect_equal(item_index_max_information, 29)
})


test_that("three dimensions, items load one or two dimensions, no variance target", {
  information_summary <- "posterior_determinant"
  estimate <- c(.2, .1, 2.5)
  attr(estimate, "variance") <- diag(c(.3, .1, .2))
  model <- "3PLM"
  estimator <- "maximum_aposteriori"
  number_items <- 50
  number_dimensions <- 3
  available <- c(6:10, 21:30, 50)
  administered <- c(31:49, 5:1, 11:20)
  responses <- rep(c(1, 0), 17)
  
  stop_test <- list(max_n = number_items)
  
  prior_form <- "normal"
  prior_parameters <- list(mu = rep(0, 3), Sigma = diag(number_dimensions) * .4)
  
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  alpha[1:20, 2:3] <- 0
  alpha[21:40,c(1,3)] <- 0
  
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  item_information <- get_summarized_information(information_summary, estimate, model, responses, prior_form, prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item)
  item_index_max_information <- get_item_index_max_information(available, item_information, estimate, stop_test, alpha, number_answers = length(administered))
  
  test_item_index_max_information <- validate_and_run.test(fn = get_item_index_max_information, available = available, item_information = item_information, estimate = estimate, stop_test = stop_test, alpha = alpha, number_answers = length(administered))
  helper_get_uncompleted_dimensions <- get('get_uncompleted_dimensions', environment(test_item_index_max_information))
  helper_get_useful_item_indices <- get('get_useful_item_indices', environment(test_item_index_max_information))
  helper_get_available_and_useful_items <- get('get_available_and_useful_items', environment(test_item_index_max_information))
  
  expect_equal(helper_get_uncompleted_dimensions(), NULL)
  expect_equal(helper_get_useful_item_indices(NULL), 1:nrow(alpha))
  expect_equal(helper_get_available_and_useful_items(1:nrow(alpha)), available)
  
  expect_equal(item_index_max_information, 6)
})

test_that("three dimensions, items load one or two dimensions, with variance target", {
  information_summary <- "posterior_determinant"
  estimate <- c(.2, .1, 2.5)
  attr(estimate, "variance") <- diag(c(.3, .1, .2))
  model <- "3PLM"
  estimator <- "maximum_aposteriori"
  number_items <- 50
  number_dimensions <- 3
  available <- c(6:10, 21:30, 50)
  administered <- c(31:49, 5:1, 11:20)
  responses <- rep(c(1, 0), 17)
  
  stop_test <- list(max_n = number_items, target = c(.2, .2, .2))
  
  prior_form <- "normal"
  prior_parameters <- list(mu = rep(0, 3), Sigma = diag(number_dimensions) * .4)
  
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  alpha[1:20, 2:3] <- 0
  alpha[21:40,c(1,3)] <- 0
  
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  item_information <- get_summarized_information(information_summary, estimate, model, responses, prior_form, prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item)
  item_index_max_information <- get_item_index_max_information(available, item_information, estimate, stop_test, alpha, number_answers = length(administered))
  
  test_item_index_max_information <- validate_and_run.test(fn = get_item_index_max_information, available = available, item_information = item_information, estimate = estimate, stop_test = stop_test, alpha = alpha, number_answers = length(administered))
  helper_get_uncompleted_dimensions <- get('get_uncompleted_dimensions', environment(test_item_index_max_information))
  helper_get_useful_item_indices <- get('get_useful_item_indices', environment(test_item_index_max_information))
  helper_get_available_and_useful_items <- get('get_available_and_useful_items', environment(test_item_index_max_information))
  
  expect_equal(helper_get_uncompleted_dimensions(), c(1, 3))
  expect_equal(helper_get_useful_item_indices(c(1, 3)), c(1:20, 41:50))
  expect_equal(helper_get_available_and_useful_items(c(1:20, 41:50)), c(6:10, 50))
  
  expect_equal(item_index_max_information, 6)
})

test_that("no items left that are both available and useful", {
  information_summary <- "posterior_determinant"
  estimate <- c(.2, .1, 2.5)
  attr(estimate, "variance") <- diag(c(.3, .1, .2))
  model <- "3PLM"
  estimator <- "maximum_aposteriori"
  number_items <- 50
  number_dimensions <- 3
  available <- c(1:24)
  administered <- c(25:50)
  responses <- rep(c(1, 0), 12)
  
  stop_test <- list(max_n = number_items, target = c(.2, .2, .2))
  
  prior_form <- "normal"
  prior_parameters <- list(mu = rep(0, 3), Sigma = diag(number_dimensions) * .4)
  
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  alpha[1:30, c(1, 3)] <- 0
  
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  item_information <- get_summarized_information(information_summary, estimate, model, responses, prior_form, prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item)
  item_index_max_information <- get_item_index_max_information(available, item_information, estimate, stop_test, alpha, number_answers = length(administered))
  
  test_item_index_max_information <- validate_and_run.test(fn = get_item_index_max_information, available = available, item_information = item_information, estimate = estimate, stop_test = stop_test, alpha = alpha, number_answers = length(administered))
  helper_get_uncompleted_dimensions <- get('get_uncompleted_dimensions', environment(test_item_index_max_information))
  helper_get_useful_item_indices <- get('get_useful_item_indices', environment(test_item_index_max_information))
  helper_get_available_and_useful_items <- get('get_available_and_useful_items', environment(test_item_index_max_information))
  
  expect_equal(helper_get_uncompleted_dimensions(), c(1, 3))
  expect_equal(helper_get_useful_item_indices(c(1, 3)), c(31:50))
  expect_equal(helper_get_available_and_useful_items(c(31:50)), available)
  expect_equal(item_index_max_information, 8)
})

test_that("one a single appropriate item", {
  information_summary <- "posterior_determinant"
  estimate <- c(.2, .1, .5)
  attr(estimate, "variance") <- diag(c(.3, .1, .02))
  model <- "3PLM"
  estimator <- "maximum_aposteriori"
  number_items <- 50
  number_dimensions <- 3
  available <- c(6:10, 21:30, 50)
  administered <- c(31:49, 5:1, 11:20)
  responses <- rep(c(1, 0), 17)
  
  stop_test <- list(max_n = number_items, target = c(.2, .2, .2))
  
  prior_form <- "normal"
  prior_parameters <- list(mu = rep(0, 3), Sigma = diag(number_dimensions) * .4)
  
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  alpha[1:49, 1:2] <- 0
  
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
  
  item_information <- get_summarized_information(information_summary, estimate, model, responses, prior_form, prior_parameters, available, administered, number_items, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item)
  item_index_max_information <- get_item_index_max_information(available, item_information, estimate, stop_test, alpha, number_answers = length(administered))
  
  test_item_index_max_information <- validate_and_run.test(fn = get_item_index_max_information, available = available, item_information = item_information, estimate = estimate, stop_test = stop_test, alpha = alpha, number_answers = length(administered))
  helper_get_uncompleted_dimensions <- get('get_uncompleted_dimensions', environment(test_item_index_max_information))
  helper_get_useful_item_indices <- get('get_useful_item_indices', environment(test_item_index_max_information))
  
  expect_equal(helper_get_uncompleted_dimensions(), 1)
  expect_equal(helper_get_useful_item_indices(1), 50)
  expect_equal(item_index_max_information, 50)
})


