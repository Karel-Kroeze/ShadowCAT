# only for whithin R:
'
library(testthat)
library(lpSolve)
'

make_random_seed_exist <- rnorm(1)

# create item characteristics
model <- '3PLM'
number_items <- 50
number_dimensions <- 1
number_answer_categories <- 2 # can only be 2 for 3PLM model
guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
eta <- NULL # only relevant for GPCM model

alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)

item_characteristics_shadowcat_format <- initItembank(model = model, alpha = alpha, beta = beta, guessing = guessing, eta = eta, silent = TRUE)

#create item characteristics and constraints
characteristics <- data.frame(content = with_random_seed(2, sample)(c('algebra','physics','calculus'), item_characteristics_shadowcat_format$K, TRUE),
                              time = with_random_seed(2, rnorm)(item_characteristics_shadowcat_format$K),
                              exclusive = sapply(1:item_characteristics_shadowcat_format$K, FUN = function (x) { if (x %in% with_random_seed(2, sample)(item_characteristics_shadowcat_format$K, 4)) 1 else 0 }))
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

initiated_test <- initTest(item_characteristics_shadowcat_format,
                           start = list(type = 'fixed', indices = c(2, 4, 5), n = 3),
                           stop = list(type = 'variance', target = .2, n = 20),
                           max_n = 20,
                           estimator = 'MAP',
                           objective = 'PD',
                           selection = 'MI',
                           constraints = list(constraints = constraints,
                                              characteristics = characteristics),
                           exposure = NULL,
                           lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                           upperBound = rep(3, item_characteristics_shadowcat_format$Q))

# get initiated person
initiated_person <- initPerson(item_characteristics_shadowcat_format, theta = rep(.2, item_characteristics_shadowcat_format$Q), prior = .4)

test_that("First 4 items administered", {
  initiated_person$responses = rep(c(1, 0), 2)
  initiated_person$available <- c(5:50)
  initiated_person$administered <- c(1:4)

  item_information <- objective(initiated_test, initiated_person, pad = TRUE)

  best_item <- Shadow(initiated_test, initiated_person, item_information)

  expect_equal(best_item, 5)
})

test_that("7 items administered", {
  initiated_person$responses = c(1, 0, 1, 0, 1, 0, 1)
  initiated_person$available <- c(2, 4, 6, 8, 10, 12, 14:50)
  initiated_person$administered <- c(1, 3, 5, 7, 9, 11, 13)
  
  item_information <- objective(initiated_test, initiated_person, pad = TRUE)
  
  best_item <- Shadow(initiated_test, initiated_person, item_information)
  
  expect_equal(best_item, 6)
})

test_that("First 10 items administered", {
  initiated_person$responses = rep(1, 10)
  initiated_person$available <- c(11:50)
  initiated_person$administered <- c(1:10)
  
  item_information <- objective(initiated_test, initiated_person, pad = TRUE)
  
  best_item <- Shadow(initiated_test, initiated_person, item_information)
  
  expect_equal(best_item, 47)
})



