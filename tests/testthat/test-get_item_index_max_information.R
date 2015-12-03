# only for whithin R:
'
library(testthat)
'

make_random_seed_exist <- rnorm(1)

test_that("with padding in objective function", {
  # create item characteristics
  model <- "3PLM"
  number_items <- 50
  number_dimensions <- 1
  number_answer_categories <- 2 # can only be 2 for 3PLM model
  guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
  eta <- NULL # only relevant for GPCM model
  
  alpha <- matrix(with_random_seed(2, runif)(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
  beta <- matrix(with_random_seed(2, rnorm)(number_items), nrow = number_items, ncol = 1)
  
  item_characteristics_shadowcat_format <- initItembank(model = model, alpha = alpha, beta = beta, guessing = guessing, eta = eta, silent = TRUE)
  
  # get initiated test: adaptive test rules
  initiated_test <- initTest(item_characteristics_shadowcat_format, 
                             start = list(type = 'random', n = 5), 
                             stop = list(type = 'length', n = 30),
                             max_n = 50, # utter maximum
                             estimator = 'ML',
                             objective = 'PD',
                             selection = 'MI',
                             constraints = NULL,
                             exposure = NULL,
                             lowerBound = rep(-3, item_characteristics_shadowcat_format$Q),
                             upperBound = rep(3, item_characteristics_shadowcat_format$Q))
  
  # get initiated person
  initiated_person <- initPerson(item_characteristics_shadowcat_format, theta = rep(.2, item_characteristics_shadowcat_format$Q), prior = .4)
  initiated_person$available <- c(6:10, 21:30, 50)
  initiated_person$administered <- c(1:5, 11:20, 31:49)
  initiated_person$responses <- rep(c(1, 0), 17)
  
  item_information <- objective(initiated_test, initiated_person, pad = TRUE)
  item_index_max_information <- get_item_index_max_information(initiated_person$available, item_information)
  
  expect_equal(item_index_max_information, 6)
})


