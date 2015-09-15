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
                           stop = list(type = 'variance', target = .2),
                           max_n = 50, # utter maximum
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

test_that("available items large", {
  initiated_person$responses = rep(c(1, 0), 12) 
  initiated_person$available <- c(6:10, 11:30, 50)
  initiated_person$administered <- c(1:5, 31:49)
  
  item_information <- objective(initiated_test, initiated_person, pad = TRUE)
  
  best_item <- Shadow(initiated_test, initiated_person, item_information)
  
  expect_equal(best_item, 6)  
})

test_that("available items small", {
  initiated_person$responses = rep(c(1, 0), 12) 
  initiated_person$available <- c(7:10)
  initiated_person$administered <- c(1:6, 11:50)
  
  item_information <- objective(initiated_test, initiated_person, pad = TRUE)
  
  best_item <- Shadow(initiated_test, initiated_person, item_information)
  
  expect_equal(best_item, 8)
})

test_that("available items one", {
  initiated_person$responses = rep(c(1, 0), 12) 
  initiated_person$available <- c(10)
  initiated_person$administered <- c(1:9, 11:50)
  
  item_information <- objective(initiated_test, initiated_person, pad = TRUE)
  
  best_item <- Shadow(initiated_test, initiated_person, item_information)
  
  expect_equal(best_item, 10)
})

