# only for whithin R:
'
library(testthat)
'
context("all matrices one column")

test_that("3PLM model, 1 dimension, 2 categories", {
  model <- '3PLM'
  number_items <- 50
  number_dimensions <- 1
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  eta <- NULL # only relevant for GPCM model
  subset_indices <- c(2, 4, 7, 11, 20:25, 28)
  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  
  item_characteristics_shadowcat_format <- initItembank(model = model, 
                                                        alpha = alpha,
                                                        beta = beta,
                                                        guessing = guessing, 
                                                        eta = eta,
                                                        silent = TRUE)
  item_characteristics_subset <- subset.ShadowCAT.items(item_characteristics_shadowcat_format, subset_indices)
  
  expect_equal(round(item_characteristics_subset$pars$alpha[c(1, 2, 11),], 3), c(1.143, .502, .728))
  expect_equal(dim(item_characteristics_subset$pars$alpha), c(11, 1))
  expect_equal(round(item_characteristics_subset$pars$beta[c(1, 2, 11),], 3), c(.185, -1.130, -.597))
  expect_equal(dim(item_characteristics_subset$pars$beta), c(11, 1))
  expect_equal(round(item_characteristics_subset$pars$guessing[c(1, 2, 11),], 3), c(.1, .1, .2))
  expect_equal(dim(item_characteristics_subset$pars$guessing), c(11, 1))
  expect_equal(item_characteristics_subset$pars$index, c(2, 4, 7, 11, 20:25, 28))
  expect_equal(length(item_characteristics_subset$pars$index), 11)
  expect_equal(item_characteristics_subset$pars$m, rep(1, 11))
  expect_equal(length(item_characteristics_subset$pars$index), 11)
  expect_equal(item_characteristics_subset$Q, 1)
  expect_equal(item_characteristics_subset$K, 11)
  expect_equal(item_characteristics_subset$M, 1)
  expect_equal(item_characteristics_subset$model, "3PLM")
  expect_equal(item_characteristics_subset$subset, c(2, 4, 7, 11, 20:25, 28))  
})

context("matrices three columns")

test_that("GPCM model, 3 dimension, 4 categories", {
  model <- 'GPCM'
  number_items <- 50
  number_dimensions <- 3
  number_answer_categories <- 4
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  subset_indices <- c(2, 4, 7, 11, 20:25, 28)
  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  eta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (number_answer_categories - 1))))
  
  item_characteristics_shadowcat_format <- initItembank(model = model, 
                                                        alpha = alpha,
                                                        beta = NULL,
                                                        guessing = guessing, 
                                                        eta = eta,
                                                        silent = TRUE)
  item_characteristics_subset <- subset.ShadowCAT.items(item_characteristics_shadowcat_format, subset_indices)
  
  expect_equal(round(item_characteristics_subset$pars$alpha[2,], 3), c(.502, 1.416, 1.295))
  expect_equal(round(item_characteristics_subset$pars$alpha[11,], 3), c(.728, 1.341, .666))
  expect_equal(dim(item_characteristics_subset$pars$alpha), c(11, 3))
  expect_equal(round(item_characteristics_subset$pars$beta[1,], 3), c(-1.815, -1.630, 0.555))
  expect_equal(round(item_characteristics_subset$pars$beta[10,], 3), c(-1.995, -1.990, .015))
  expect_equal(dim(item_characteristics_subset$pars$beta), c(11, 3))
  expect_equal(round(item_characteristics_subset$pars$guessing[c(1, 2, 11),], 3), c(.1, .1, .2))
  expect_equal(dim(item_characteristics_subset$pars$guessing), c(11, 1))
  expect_equal(item_characteristics_subset$pars$index, c(2, 4, 7, 11, 20:25, 28))
  expect_equal(length(item_characteristics_subset$pars$index), 11)
  expect_equal(item_characteristics_subset$pars$m, rep(3, 11))
  expect_equal(length(item_characteristics_subset$pars$index), 11)
  expect_equal(item_characteristics_subset$Q, 3)
  expect_equal(item_characteristics_subset$K, 11)
  expect_equal(item_characteristics_subset$M, 3)
  expect_equal(item_characteristics_subset$model, "GPCM")
  expect_equal(item_characteristics_subset$subset, c(2, 4, 7, 11, 20:25, 28)) 
  
})
