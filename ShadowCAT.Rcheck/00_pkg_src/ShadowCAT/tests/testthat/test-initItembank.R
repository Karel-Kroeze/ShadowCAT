# only for whithin R:
# library(testthat)

# Something goes wrong with with_random_seed
if (FALSE) {

test_that("3PLM model, 1 dimension, 2 categories, guessing .1 and .2, ", {
  model <- '3PLM'
  number_items <- 50
  number_dimensions <- 1
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  eta <- NULL # only relevant for GPCM model
  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
     
  item_characteristics_shadowcat_format <- initItembank(model = model, 
                                                        alpha = simulated_item_characteristics$pars$alpha,
                                                        beta = simulated_item_characteristics$pars$beta,
                                                        guessing = guessing, 
                                                        eta = eta,
                                                        silent = TRUE)
  
  expect_equal(round(item_characteristics_shadowcat_format$pars$alpha[5], 3), 0.437)
  expect_equal(round(item_characteristics_shadowcat_format$pars$alpha[50], 3), 0.611)
  expect_equal(round(item_characteristics_shadowcat_format$pars$beta[1], 3), 0.371)
  expect_equal(round(item_characteristics_shadowcat_format$pars$beta[45], 3), -1.105)
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
}