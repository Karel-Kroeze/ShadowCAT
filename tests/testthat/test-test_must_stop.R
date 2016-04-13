# only for whithin R:
'
library(testthat)
'

context("stop rule is variance")

test_that("stop rule is variance, targets not reached, one dimension", {
  estimate <- 1
  attr(estimate, "variance") <- matrix(.4)
  should_stop <- test_must_stop(number_answers = 15, estimate, min_n = NULL, max_n = 20, stop_variance_target = .2)
  expect_equal(should_stop, FALSE)
})

test_that("stop rule is variance, targets reached, one dimension", {
  estimate <- 1
  attr(estimate, "variance") <- matrix(.1)
  should_stop <- test_must_stop(number_answers = 15, estimate, min_n = NULL, max_n = 20, stop_variance_target = .2)
  expect_equal(should_stop, TRUE)
})

test_that("stop rule is variance, targets not reached, three dimensions", {
  estimate <- 1
  attr(estimate, "variance") <- diag(c(.4, .5, .1))
  should_stop <- test_must_stop(number_answers = 15, estimate, min_n = NULL, max_n = 20, stop_variance_target = c(.2, .2, .2))
  expect_equal(should_stop, FALSE)
})

test_that("stop rule is variance, variance reached, three dimensions", {
  estimate <- 1
  attr(estimate, "variance") <- diag(c(.19, .1, .1))
  should_stop <- test_must_stop(number_answers = 15, estimate, min_n = 15, max_n = 20, stop_variance_target = c(.2, .2, .2))
  expect_equal(should_stop, TRUE)
})

test_that("stop rule is variance, maximum number of items reached, three dimensions", {
  estimate <- 1
  attr(estimate, "variance") <- diag(c(.6, .1, .1))
  should_stop <- test_must_stop(number_answers = 20, estimate, min_n = NULL, max_n = 20, stop_variance_target = c(.2, .2, .2))
  expect_equal(should_stop, TRUE)
})

test_that("stop rule is variance, variance reached, minimum number of items not reached, three dimensions", {
  estimate <- 1
  attr(estimate, "variance") <- diag(c(.19, .1, .1))
  should_stop <- test_must_stop(number_answers = 15, estimate, min_n = 16, max_n = 20, stop_variance_target = c(.2, .2, .2))
  expect_equal(should_stop, FALSE)
})


context("stop rule is number of items")

test_that("stop rule is number of items, target not reached", { 
  estimate <- 1
  attr(estimate, "variance") <- diag(c(.1, .1, .1))
  should_stop <- test_must_stop(number_answers = 14, estimate, min_n = NULL, max_n = 15)
  expect_equal(should_stop, FALSE)
})

test_that("stop rule is number of items, target reached", {
  estimate <- 1
  attr(estimate, "variance") <- diag(c(.1, .1, .1))
  should_stop <- test_must_stop(number_answers = 15, estimate, min_n = NULL, max_n = 15)
  expect_equal(should_stop, TRUE)
})

context("stop rule is cut off")

cutoffs <- with_random_seed(2, matrix)(runif(75, 1, 2), ncol = 3)

test_that("stop rule is cut off, not far enough below cut off, one dimension", {
  estimate <- .8
  attr(estimate, "variance") <- matrix(.4^2)
  should_stop <- test_must_stop(number_answers = 15, estimate, min_n = NULL, max_n = 20, stop_variance_target = .1, cutoffs = cutoffs[,1, drop = FALSE])
  expect_equal(should_stop, FALSE)
})

test_that("stop rule is cut off, far enough below cut off, one dimension", {
  estimate <- .3
  attr(estimate, "variance") <- matrix(.4^2)
  should_stop <- test_must_stop(number_answers = 15, estimate, min_n = NULL, max_n = 20, stop_variance_target = .1, cutoffs = cutoffs[,1, drop = FALSE])
  expect_equal(should_stop, TRUE)
})

test_that("stop rule is cut off, not far enough below cut off, three dimensions", {
  estimate <- c(.7, 1.1, -2.5)
  attr(estimate, "variance") <- diag(c(.4, .1, .2)^2)
  should_stop <- test_must_stop(number_answers = 15, estimate, min_n = NULL, max_n = 20, stop_variance_target = .1, cutoffs = cutoffs)
  expect_equal(should_stop, FALSE)
})

test_that("stop rule is cut off, far enough below cut off, three dimensions", {
  estimate <- c(.5, 1.1, -2.5)
  attr(estimate, "variance") <- diag(c(.4, .1, .2)^2)
  should_stop <- test_must_stop(number_answers = 15, estimate, min_n = NULL, max_n = 20, stop_variance_target = .1, cutoffs = cutoffs)
  expect_equal(should_stop, TRUE)
})

test_that("stop rule is cut off, far enough below cut off, three dimensions, but minimum number of items not reached", {
  estimate <- c(.5, 1.1, -2.5)
  attr(estimate, "variance") <- diag(c(.4, .1, .2)^2)
  should_stop <- test_must_stop(number_answers = 15, estimate, min_n = 16, max_n = 20, stop_variance_target = .1, cutoffs = cutoffs)
  expect_equal(should_stop, FALSE)
})

test_that("stop rule is cut off, not far enough below cut off, three dimensions, but max number of items reached", {
  estimate <- c(.7, 1.1, -2.5)
  attr(estimate, "variance") <- diag(c(.4, .1, .2)^2)
  should_stop <- test_must_stop(number_answers = 75, estimate, min_n = NULL, max_n = 75, stop_variance_target = .1, cutoffs = cutoffs)
  expect_equal(should_stop, TRUE)
})

