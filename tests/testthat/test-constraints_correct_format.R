# only for whithin R:
'
library(testthat)
library(lpSolve)
'

make_random_seed_exist <- rnorm(1)

test_that("no charateristics and constraints", {
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
  
  initiated_test <- initTest(item_characteristics_shadowcat_format, 
                             start = list(type = 'fixed', indices = c(2, 4, 5), n = 3), 
                             stop = list(type = 'variance', target = .2),
                             max_n = 50, # utter maximum
                             estimator = 'MAP',
                             objective = 'PD',
                             selection = 'MI',
                             constraints = NULL,
                             exposure = NULL,
                             lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                             upperBound = rep(3, item_characteristics_shadowcat_format$Q))
  
  constraints_formatted <- createConstraints(initiated_test, characteristics = NULL, constraints = NULL)

  expect_equal(length(constraints_formatted), 3)
  expect_equal(constraints_formatted$characteristics, data.frame(length = rep(1, 50)))
  expect_equal(constraints_formatted$constraints, data.frame(name = "length", op = "=", target = 50, stringsAsFactors = FALSE))
  expect_equal(constraints_formatted$lp_chars, data.frame(length = rep(1, 50))) 
})

test_that("constraints 1", {
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
                             stop = list(type = 'length', n = 30),
                             max_n = 50, # utter maximum
                             estimator = 'MAP',
                             objective = 'PD',
                             selection = 'MI',
                             constraints = list(constraints = constraints, 
                                                characteristics = characteristics),
                             exposure = NULL,
                             lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                             upperBound = rep(3, item_characteristics_shadowcat_format$Q))
  
  constraints_formatted <- createConstraints(initiated_test, characteristics = characteristics, constraints = constraints)
  
  expect_equal(length(constraints_formatted), 3)
  expect_equal(names(constraints_formatted$characteristics), c("length", "content/algebra", "content/calculus", "content/physics", "time", "exclusive"))
  expect_equal(names(constraints_formatted$lp_chars), c("length", "content/algebra", "content/algebra.1", "content/physics", "content/physics.1", "time", "exclusive"))
  
  expect_equal(unname(round(constraints_formatted$characteristics[1:5,], 3)), unname(data.frame(length = rep(1, 5), 
                                                                                                algebra = c(1, 0, 0, 1, 0),
                                                                                                calculus = c(0, 1, 0, 0, 1),
                                                                                                psychics = c(0, 0, 1, 0, 0),
                                                                                                time = c(-.897, .185, 1.588, -1.130, -.080),
                                                                                                exclusive = rep(0, 5))))
  expect_equal(constraints_formatted$constraints, data.frame(name = c("length", "content/algebra", "content/algebra", "content/physics", "content/physics", "time", "exclusive"), 
                                                             op = c("=", ">", "<", ">", "<", "<", "<"), 
                                                             target = c("30", "5", "10", "2", "5", "20", "2"), 
                                                             stringsAsFactors = FALSE))
  
  expect_equal(unname(round(constraints_formatted$lp_chars[1:5,], 3)), unname(data.frame(length = rep(1, 5), 
                                                                                          algebra = c(1, 0, 0, 1, 0),
                                                                                          algebra1 = c(1, 0, 0, 1, 0),
                                                                                          psychics = c(0, 0, 1, 0, 0),
                                                                                          psychics1 = c(0, 0, 1, 0, 0),
                                                                                          time = c(-.897, .185, 1.588, -1.130, -.080),
                                                                                          exclusive = rep(0, 5)))) 
})

test_that("constraints 2", {
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
  characteristics <- data.frame(time = with_random_seed(2, rnorm)(item_characteristics_shadowcat_format$K),
                                type = with_random_seed(2, sample)(c('depression', 'insomnia', 'anxiety'), item_characteristics_shadowcat_format$K, TRUE),
                                stressful = rep(c(0, 0, 1, 0, 0), 10),
                                extra = rep(c(1, 0, 0, 0, 0), 10))
  constraints <- list(list(name = 'time',
                           op = '><',
                           target = c(10, 20)),
                      list(name = 'type/depression',
                           op = '><',
                           target = c(2, 5)),
                      list(name = 'type/anxiety',
                           op = '><',
                           target = c(5, 10)),
                      list(name = 'stressful',
                           op = '<',
                           target = 3))
  
  initiated_test <- initTest(item_characteristics_shadowcat_format, 
                             start = list(type = 'fixed', indices = c(2, 4, 5), n = 3), 
                             stop = list(type = 'length', n = 30),
                             max_n = 50, # utter maximum
                             estimator = 'MAP',
                             objective = 'PD',
                             selection = 'MI',
                             constraints = list(constraints = constraints, 
                                                characteristics = characteristics),
                             exposure = NULL,
                             lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                             upperBound = rep(3, item_characteristics_shadowcat_format$Q))
    
  constraints_formatted <- createConstraints(initiated_test, characteristics = characteristics, constraints = constraints)
  
  expect_equal(length(constraints_formatted), 3)
  expect_equal(names(constraints_formatted$characteristics), c("length", "time", "type/depression", "type/anxiety", "type/insomnia", "stressful", "extra"))
  expect_equal(names(constraints_formatted$lp_chars), c("length", "time", "time.1", "type/depression", "type/depression.1", "type/anxiety", "type/anxiety.1", "stressful"))
  
  expect_equal(unname(round(constraints_formatted$characteristics[1:5,], 3)), unname(data.frame(length = rep(1, 5), 
                                                                                                time = c(-.897, .185, 1.588, -1.130, -.080),
                                                                                                depression = c(1, 0, 0, 1, 0),
                                                                                                anxiety = c(0, 1, 0, 0, 1),
                                                                                                insomnia = c(0, 0, 1, 0, 0),
                                                                                                stressful = c(0, 0, 1, 0, 0),
                                                                                                extra = c(1, 0, 0, 0, 0))))
  expect_equal(constraints_formatted$constraints, data.frame(name = c("length", "time", "time", "type/depression", "type/depression", "type/anxiety", "type/anxiety", "stressful"), 
                                                             op = c("=", ">", "<", ">", "<", ">", "<", "<"), 
                                                             target = c("30", "10", "20", "2", "5", "5", "10", "3"), 
                                                             stringsAsFactors = FALSE))
  
  expect_equal(unname(round(constraints_formatted$lp_chars[1:5,], 3)), unname(data.frame(length = rep(1, 5), 
                                                                                         time = c(-.897, .185, 1.588, -1.130, -.080),
                                                                                         time = c(-.897, .185, 1.588, -1.130, -.080),
                                                                                         depression = c(1, 0, 0, 1, 0),
                                                                                         depression = c(1, 0, 0, 1, 0),
                                                                                         anxiety = c(0, 1, 0, 0, 1),
                                                                                         anxiety = c(0, 1, 0, 0, 1),
                                                                                         stressful = c(0, 0, 1, 0, 0)))) 
})

