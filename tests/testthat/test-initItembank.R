# only for whithin R:
'
library(testthat)
'
make_random_seed_exist <- rnorm(1)

context("3PLM model")

test_that("3PLM model, 1 dimension, 2 categories, guessing .1 and .2", {
  model <- '3PLM'
  number_items <- 50
  number_dimensions <- 1
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  eta <- NULL # only relevant for GPCM model
  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
     
  item_characteristics_shadowcat_format <- initItembank(model = model, 
                                                        alpha = alpha,
                                                        beta = beta,
                                                        guessing = guessing, 
                                                        eta = eta,
                                                        silent = TRUE)
  
  expect_equal(round(item_characteristics_shadowcat_format$pars$alpha[5], 3), 1.433)
  expect_equal(round(item_characteristics_shadowcat_format$pars$alpha[50], 3), 1.272)
  expect_equal(round(item_characteristics_shadowcat_format$pars$beta[1], 3), -.897)
  expect_equal(round(item_characteristics_shadowcat_format$pars$beta[45], 3), .622)
  expect_equal(dim(item_characteristics_shadowcat_format$pars$alpha), c(50, 1))
  expect_equal(dim(item_characteristics_shadowcat_format$pars$beta), c(50, 1))
  expect_equal(item_characteristics_shadowcat_format$pars$guessing[c(2,26),], c(.1, .2))
  expect_equal(dim(item_characteristics_shadowcat_format$pars$guessing), c(50, 1))
  expect_equal(item_characteristics_shadowcat_format$pars$index, 1:50)
  expect_equal(item_characteristics_shadowcat_format$pars$m, rep(1, 50))
  expect_equal(item_characteristics_shadowcat_format$Q, 1)
  expect_equal(item_characteristics_shadowcat_format$K, 50)
  expect_equal(item_characteristics_shadowcat_format$M, 1)
  expect_equal(item_characteristics_shadowcat_format$model, "3PLM")
  expect_equal(class(item_characteristics_shadowcat_format), "ShadowCAT.items")
})

test_that("3PLM model, 3 dimensions, 2 categories, guessing .1 and .2", {
  model <- '3PLM'
  number_items <- 60
  number_dimensions <- 3
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  eta <- NULL # only relevant for GPCM model
  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  
  item_characteristics_shadowcat_format <- initItembank(model = model, 
                                                        alpha = alpha,
                                                        beta = beta,
                                                        guessing = guessing, 
                                                        eta = eta,
                                                        silent = TRUE)
  
  expect_equal(round(item_characteristics_shadowcat_format$pars$alpha[4,], 3), c(.502, .612, 1.363))
  expect_equal(round(item_characteristics_shadowcat_format$pars$alpha[50,], 3), c(1.272, .454, .947))
  expect_equal(round(item_characteristics_shadowcat_format$pars$beta[2,], 3), .185)
  expect_equal(round(item_characteristics_shadowcat_format$pars$beta[40,], 3), -.247)
  expect_equal(dim(item_characteristics_shadowcat_format$pars$alpha), c(60, 3))
  expect_equal(dim(item_characteristics_shadowcat_format$pars$beta), c(60, 1))
  expect_equal(item_characteristics_shadowcat_format$pars$guessing[c(30,31),], c(.1, .2))
  expect_equal(dim(item_characteristics_shadowcat_format$pars$guessing), c(60, 1))
  expect_equal(item_characteristics_shadowcat_format$pars$index, 1:60)
  expect_equal(item_characteristics_shadowcat_format$pars$m, rep(1, 60))
  expect_equal(item_characteristics_shadowcat_format$Q, 3)
  expect_equal(item_characteristics_shadowcat_format$K, 60)
  expect_equal(item_characteristics_shadowcat_format$M, 1)
  expect_equal(item_characteristics_shadowcat_format$model, "3PLM")
  expect_equal(class(item_characteristics_shadowcat_format), "ShadowCAT.items")
})

context("GPCM model")

test_that("GPCM model, 1 dimension, 2 categories, guessing .1 and .2", {
  model <- 'GPCM'
  number_items <- 50
  number_dimensions <- 1
  number_answer_categories <- 2
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  eta <- NULL # only relevant for GPCM model
  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  
  item_characteristics_shadowcat_format <- initItembank(model = model, 
                                                        alpha = alpha,
                                                        beta = beta,
                                                        guessing = guessing, 
                                                        eta = eta,
                                                        silent = TRUE)
  
  expect_equal(round(item_characteristics_shadowcat_format$pars$alpha[5], 3), 1.433)
  expect_equal(round(item_characteristics_shadowcat_format$pars$alpha[50], 3), 1.272)
  expect_equal(round(item_characteristics_shadowcat_format$pars$beta[1], 3), -.897)
  expect_equal(round(item_characteristics_shadowcat_format$pars$beta[45], 3), .622)
  expect_equal(dim(item_characteristics_shadowcat_format$pars$alpha), c(50, 1))
  expect_equal(dim(item_characteristics_shadowcat_format$pars$beta), c(50, 1))
  expect_equal(item_characteristics_shadowcat_format$pars$guessing[c(2,26),], c(.1, .2))
  expect_equal(dim(item_characteristics_shadowcat_format$pars$guessing), c(50, 1))
  expect_equal(item_characteristics_shadowcat_format$pars$index, 1:50)
  expect_equal(item_characteristics_shadowcat_format$pars$m, rep(1, 50))
  expect_equal(item_characteristics_shadowcat_format$Q, 1)
  expect_equal(item_characteristics_shadowcat_format$K, 50)
  expect_equal(item_characteristics_shadowcat_format$M, 1)
  expect_equal(item_characteristics_shadowcat_format$model, "GPCM")
  expect_equal(class(item_characteristics_shadowcat_format), "ShadowCAT.items")
})

test_that("GPCM model, 1 dimension, 2 categories, guessing .1 and .2", {
  model <- 'GPCM'
  number_items <- 50
  number_dimensions <- 1
  number_answer_categories <- 2
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  eta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  beta <- matrix((apply(eta, 1, cumsum)), ncol = dim(eta)[2])
  
  item_characteristics_shadowcat_format <- initItembank(model = model, 
                                                        alpha = alpha,
                                                        beta = beta,
                                                        guessing = guessing, 
                                                        eta = eta,
                                                        silent = TRUE)
  
  expect_equal(round(item_characteristics_shadowcat_format$pars$alpha[5], 3), 1.433)
  expect_equal(round(item_characteristics_shadowcat_format$pars$alpha[50], 3), 1.272)
  expect_equal(round(item_characteristics_shadowcat_format$pars$beta[1], 3), -.897)
  expect_equal(round(item_characteristics_shadowcat_format$pars$beta[45], 3), .622)
  expect_equal(dim(item_characteristics_shadowcat_format$pars$alpha), c(50, 1))
  expect_equal(dim(item_characteristics_shadowcat_format$pars$beta), c(50, 1))
  expect_equal(item_characteristics_shadowcat_format$pars$guessing[c(2,26),], c(.1, .2))
  expect_equal(dim(item_characteristics_shadowcat_format$pars$guessing), c(50, 1))
  expect_equal(item_characteristics_shadowcat_format$pars$index, 1:50)
  expect_equal(item_characteristics_shadowcat_format$pars$m, rep(1, 50))
  expect_equal(item_characteristics_shadowcat_format$Q, 1)
  expect_equal(item_characteristics_shadowcat_format$K, 50)
  expect_equal(item_characteristics_shadowcat_format$M, 1)
  expect_equal(item_characteristics_shadowcat_format$model, "GPCM")
  expect_equal(class(item_characteristics_shadowcat_format), "ShadowCAT.items")
  
})

test_that("GPCM model, 3 dimension, 4 categories, guessing .1 and .2", {
  model <- 'GPCM'
  number_items <- 50
  number_dimensions <- 3
  number_answer_categories <- 4
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  temp_vector <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  eta <- t(apply(temp_vector, 1, function(x) x + seq(-2, 2, length.out = (number_answer_categories - 1))))
  beta <- matrix((apply(eta, 1, cumsum)), ncol = dim(eta)[2])
  
  item_characteristics_shadowcat_format <- initItembank(model = model, 
                                                        alpha = alpha,
                                                        beta = beta,
                                                        guessing = guessing, 
                                                        eta = eta,
                                                        silent = TRUE)
  
  expect_equal(round(item_characteristics_shadowcat_format$pars$alpha[5,], 3), c(1.433, .630, .644))
  expect_equal(round(item_characteristics_shadowcat_format$pars$alpha[50,], 3), c(1.272, .828, .650))
  expect_equal(round(item_characteristics_shadowcat_format$pars$beta[1,], 3), c(-2.897, 2.636, -2.568))
  expect_equal(round(item_characteristics_shadowcat_format$pars$beta[45, ], 3),  c(5.347, -1.362, -2.184))
  expect_equal(dim(item_characteristics_shadowcat_format$pars$alpha), c(50, 3))
  expect_equal(dim(item_characteristics_shadowcat_format$pars$beta), c(50, 3))
  expect_equal(item_characteristics_shadowcat_format$pars$guessing[c(2,26),], c(.1, .2))
  expect_equal(dim(item_characteristics_shadowcat_format$pars$guessing), c(50, 1))
  expect_equal(item_characteristics_shadowcat_format$pars$index, 1:50)
  expect_equal(item_characteristics_shadowcat_format$pars$m, rep(3, 50))
  expect_equal(item_characteristics_shadowcat_format$Q, 3)
  expect_equal(item_characteristics_shadowcat_format$K, 50)
  expect_equal(item_characteristics_shadowcat_format$M, 3)
  expect_equal(item_characteristics_shadowcat_format$model, "GPCM")
  expect_equal(class(item_characteristics_shadowcat_format), "ShadowCAT.items")
})

