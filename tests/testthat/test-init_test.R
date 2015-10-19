# only for whithin R:
'
library(testthat)
'

make_random_seed_exist <- rnorm(1)

test_that("start is random, stop is n, constraints is null", {
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
  
  # get initiated test
  initiated_test <- initTest(item_characteristics_shadowcat_format, 
                             start = list(type = 'random', n = 5), 
                             stop = list(type = 'length', n = 30),
                             max_n = 50, # utter maximum
                             estimator = 'MAP',
                             objective = 'PD',
                             selection = 'MI',
                             constraints = NULL,
                             exposure = NULL,
                             lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                             upperBound = rep(3, item_characteristics_shadowcat_format$Q))
  
  expect_equal(round(initiated_test$items$pars$alpha[2, 1], 3), 1.143)
  expect_equal(initiated_test$start$type, "random")
  expect_equal(initiated_test$start$n, 5)
  expect_equal(initiated_test$stop$type, "length")
  expect_equal(initiated_test$stop$n, 30)
  expect_equal(initiated_test$max_n, 50)
  expect_equal(initiated_test$lowerBound, -3)
  expect_equal(initiated_test$upperBound, 3)
  expect_equal(initiated_test$estimator, "MAP")
  expect_equal(initiated_test$objective, "PD")
  expect_equal(initiated_test$selection, "MI")
  expect_equal(initiated_test$constraints$characteristics, data.frame(length = rep(1, 50)))
  expect_equal(initiated_test$constraints$constraints, data.frame(name = "length", op = "=", target = 30, stringsAsFactors = FALSE))
  expect_equal(initiated_test$constraints$lp_chars, data.frame(length = rep(1, 50)))  
  expect_equal(initiated_test$internal, list())  
  
})


test_that("start is fixed, stop is variance, constraints is defined", {
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
  
  #create characteristics and constraints
  characteristics <- data.frame(content = with_random_seed(2, sample)(c('algebra','physics','calculus'), item_characteristics_shadowcat_format$K, TRUE),
                                time = with_random_seed(2, rnorm)(item_characteristics_shadowcat_format$K),
                                exclusive = sapply(1:item_characteristics_shadowcat_format$K, FUN = function (x) { if (x %in% with_random_seed(2, sample)(item_characteristics_shadowcat_format$K, 4)) 1 else 0 }))
  constraints <- list(list(name = 'content/algebra',
                           op = '><',
                           target = c(5,10)),
                      list(name = 'content/physics',
                           op = '><',
                           target = c(2,5)),
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
  
  expect_equal(round(initiated_test$items$pars$alpha[2, 1], 3), 1.143)
  expect_equal(initiated_test$start$type, "fixed")
  expect_equal(initiated_test$start$indices, c(2, 4, 5))
  expect_equal(initiated_test$start$n, 3)
  expect_equal(initiated_test$stop$type, "variance")
  expect_equal(initiated_test$stop$target, .2)
  expect_equal(initiated_test$max_n, 50)
  expect_equal(initiated_test$lowerBound, -3)
  expect_equal(initiated_test$upperBound, 3)
  expect_equal(initiated_test$estimator, "MAP")
  expect_equal(initiated_test$objective, "PD")
  expect_equal(initiated_test$selection, "MI")
  expect_equal(initiated_test$constraints$characteristics$length, rep(1, 50))
  expect_equal(as.matrix(initiated_test$constraints$characteristics["content/algebra"])[1:4,], c(1, 0, 0, 1))
  expect_equal(as.matrix(initiated_test$constraints$characteristics["content/calculus"])[1:4,], c(0, 1, 0, 0))
  expect_equal(as.matrix(initiated_test$constraints$characteristics["content/physics"])[1:4,], c(0, 0, 1, 0))
  expect_equal(round(as.matrix(initiated_test$constraints$characteristics$time)[1:4,], 3), c(-.897, .185, 1.588, -1.130))
  expect_equal(round(as.matrix(initiated_test$constraints$characteristics$exclusive)[8:11,], 3), c(1, 0, 1, 0))
  expect_equal(as.matrix(initiated_test$constraints$constraints$name)[1:4,], c("length", "content/algebra", "content/algebra", "content/physics"))
  expect_equal(as.matrix(initiated_test$constraints$constraints$op)[1:4,], c("=", ">", "<", ">"))
  expect_equal(as.matrix(initiated_test$constraints$constraints$target)[1:4,], c("50", "5", "10", "2"))
  expect_equal(initiated_test$constraints$lp_chars$length, rep(1, 50))
  expect_equal(as.matrix(initiated_test$constraints$lp_chars["content/algebra"])[1:4,], c(1, 0, 0, 1))
  expect_equal(as.matrix(initiated_test$constraints$lp_chars["content/algebra.1"])[1:4,], c(1, 0, 0, 1))
  expect_equal(as.matrix(initiated_test$constraints$lp_chars["content/physics"])[1:4,], c(0, 0, 1, 0))
  expect_equal(as.matrix(initiated_test$constraints$lp_chars["content/physics.1"])[1:4,], c(0, 0, 1, 0))
  expect_equal(round(as.matrix(initiated_test$constraints$lp_chars$time)[1:4,], 3), c(-.897, .185, 1.588, -1.130))
  expect_equal(round(as.matrix(initiated_test$constraints$lp_chars$exclusive)[8:11,], 3), c(1, 0, 1, 0))
  expect_equal(initiated_test$internal, list())    
})
