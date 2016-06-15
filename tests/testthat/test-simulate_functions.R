# only for whithin R:
'
library(testthat)
'

make_random_seed_exist <- rnorm(1)

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
