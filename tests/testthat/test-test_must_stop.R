# only for whithin R:
'
library(testthat)
'

context("stop rule is variance")

test_that("stop rule is variance, targets not reached, one dimension", {
  should_stop <- test_must_stop(number_responses = 15, variance_current_estimate = matrix(.4), min_n = NULL, max_n = 20, stop_variance_target = .2)
  expect_equal(should_stop, FALSE)
})

test_that("stop rule is variance, targets reached, one dimension", {
  should_stop <- test_must_stop(number_responses = 15, variance_current_estimate = matrix(.1), min_n = NULL, max_n = 20, stop_variance_target = .2)
  expect_equal(should_stop, TRUE)
})

test_that("stop rule is variance, targets not reached, three dimensions", {
  should_stop <- test_must_stop(number_responses = 15, variance_current_estimate = diag(c(.4, .5, .1)), min_n = NULL, max_n = 20, stop_variance_target = c(.2, .2, .2))
  expect_equal(should_stop, FALSE)
})

test_that("stop rule is variance, variance reached, three dimensions", {
  should_stop <- test_must_stop(number_responses = 15, variance_current_estimate = diag(c(.19, .1, .1)), min_n = 15, max_n = 20, stop_variance_target = c(.2, .2, .2))
  expect_equal(should_stop, TRUE)
})

test_that("stop rule is variance, maximum number of items reached, three dimensions", {
  should_stop <- test_must_stop(number_responses = 20, variance_current_estimate = diag(c(.6, .1, .1)), min_n = NULL, max_n = 20, stop_variance_target = c(.2, .2, .2))
  expect_equal(should_stop, TRUE)
})

test_that("stop rule is variance, variance reached, minimum number of items not reached, three dimensions", {
  should_stop <- test_must_stop(number_responses = 15, variance_current_estimate = diag(c(.19, .1, .1)), min_n = 16, max_n = 20, stop_variance_target = c(.2, .2, .2))
  expect_equal(should_stop, FALSE)
})


context("stop rule is number of items")

test_that("stop rule is number of items, target not reached", { 
  should_stop <- test_must_stop(number_responses = 14, variance_current_estimate = diag(c(.1, .1, .1)), min_n = NULL, max_n = 15)
  expect_equal(should_stop, FALSE)
})

test_that("stop rule is number of items, target reached", {
  should_stop <- test_must_stop(number_responses = 15, variance_current_estimate = diag(c(.1, .1, .1)), min_n = NULL, max_n = 15)
  expect_equal(should_stop, TRUE)
})

