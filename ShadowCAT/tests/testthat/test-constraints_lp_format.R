# only for whithin R:
'
library(testthat)
library(lpSolve)
'

make_random_seed_exist <- rnorm(1)

context("valid input")

test_that("no charateristics and constraints", {  
  constraints_formatted <- constraints_lp_format(max_n = 50, number_items = 50, characteristics = NULL, constraints = NULL)

  expect_equal(length(constraints_formatted), 2)
  expect_equal(constraints_formatted$lp_constraints, data.frame(name = "length", op = "=", target = 50, stringsAsFactors = FALSE))
  expect_equal(constraints_formatted$lp_chars, data.frame(length = rep(1, 50))) 
})

test_that("constraints 1", {
  number_items <- 50
  max_n <- 30
  
  #create item characteristics and constraints
  characteristics <- data.frame(content = with_random_seed(2, sample)(c('algebra','physics','calculus'), number_items, TRUE),
                                time = with_random_seed(2, rnorm)(number_items),
                                exclusive = sapply(1:number_items, FUN = function (x) { if (x %in% with_random_seed(2, sample)(number_items, 4)) 1 else 0 }))
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
  
  constraints_formatted <- constraints_lp_format(max_n = max_n, number_items = number_items, characteristics = characteristics, constraints = constraints)
  
  expect_equal(length(constraints_formatted), 2)
  expect_equal(names(constraints_formatted$lp_chars), c("length", "content/algebra", "content/algebra.1", "content/physics", "content/physics.1", "time", "exclusive"))
  expect_equal(constraints_formatted$lp_constraints, data.frame(name = c("length", "content/algebra", "content/algebra", "content/physics", "content/physics", "time", "exclusive"), 
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
  number_items <- 50
  max_n <- 30
  
  #create item characteristics and constraints
  characteristics <- data.frame(time = with_random_seed(2, rnorm)(number_items),
                                type = with_random_seed(2, sample)(c('depression', 'insomnia', 'anxiety'), number_items, TRUE),
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
      
  constraints_formatted <- constraints_lp_format(max_n = max_n, number_items = number_items, characteristics = characteristics, constraints = constraints)
  
  expect_equal(length(constraints_formatted), 2)
  expect_equal(names(constraints_formatted$lp_chars), c("length", "time", "time.1", "type/depression", "type/depression.1", "type/anxiety", "type/anxiety.1", "stressful"))
  expect_equal(constraints_formatted$lp_constraints, data.frame(name = c("length", "time", "time", "type/depression", "type/depression", "type/anxiety", "type/anxiety", "stressful"), 
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

test_that("constraints 3", {
  number_items <- 50
  max_n <- 30
  
  #create item characteristics and constraints
  characteristics <- data.frame(type = with_random_seed(2, sample)(c('depression', 'insomnia', 'anxiety'), number_items, TRUE))
  constraints <- list(list(name = 'type/depression',
                           op = '><',
                           target = c(2, 5)),
                      list(name = 'type/anxiety',
                           op = '><',
                           target = c(2, 4)))
    
  constraints_formatted <- constraints_lp_format(max_n = max_n, number_items = number_items, characteristics = characteristics, constraints = constraints)
  
  expect_equal(length(constraints_formatted), 2)
  expect_equal(names(constraints_formatted$lp_chars), c("length", "type/depression", "type/depression.1", "type/anxiety", "type/anxiety.1"))
  expect_equal(constraints_formatted$lp_constraints, data.frame(name = c("length", "type/depression", "type/depression", "type/anxiety", "type/anxiety"), 
                                                             op = c("=", ">", "<", ">", "<"), 
                                                             target = c("30", "2", "5", "2", "4"), 
                                                             stringsAsFactors = FALSE))
  
  expect_equal(unname(round(constraints_formatted$lp_chars[1:5,], 3)), unname(data.frame(length = rep(1, 5),
                                                                                         depression = c(1, 0, 0, 1, 0),
                                                                                         depression = c(1, 0, 0, 1, 0),
                                                                                         anxiety = c(0, 1, 0, 0, 1),
                                                                                         anxiety = c(0, 1, 0, 0, 1)))) 
})

context("invalid input")

test_that("only characteristics defined", {
  number_items <- 50
  max_n <- 30
  
  characteristics <- data.frame(time = with_random_seed(2, rnorm)(number_items),
                                type = with_random_seed(2, sample)(c('depression', 'insomnia', 'anxiety'), number_items, TRUE),
                                stressful = rep(c(0, 0, 1, 0, 0), 10),
                                extra = rep(c(1, 0, 0, 0, 0), 10))
  constraints_formatted <- constraints_lp_format(max_n = max_n, number_items = number_items, characteristics = characteristics)
  expect_equal(constraints_formatted$errors$constraints_and_characteristics,"should either be defined both or not at all")
})

test_that("characteristics too few rows", {
  number_items <- 50
  max_n <- 30
  
  characteristics <- data.frame(time = with_random_seed(2, rnorm)(40),
                                type = with_random_seed(2, sample)(c('depression', 'insomnia', 'anxiety'), 40, TRUE),
                                stressful = rep(c(0, 0, 1, 0, 0), 8))
  constraints <- list(list(name = 'time',
                           op = '><',
                           target = c(10, 20)),
                      list(name = 'type/depression',
                           op = '><',
                           target = c(2, 5)),
                      list(name = 'type/insomnia',
                           op = '><',
                           target = c(2, 5)),
                      list(name = 'type/anxiety',
                           op = '><',
                           target = c(5, 10)),
                      list(name = 'stressful',
                           op = '<',
                           target = 3))
  constraints_formatted <- constraints_lp_format(max_n = max_n, number_items = number_items, characteristics = characteristics, constraints = constraints)
  expect_equal(constraints_formatted$errors$characteristics, "should be a data frame with number of rows equal to the number of items in the item bank")
})

test_that("characteristics is list", {
  number_items <- 50
  max_n <- 30
  
  characteristics <- list(time = with_random_seed(2, rnorm)(number_items),
                          type = with_random_seed(2, sample)(c('depression', 'insomnia', 'anxiety'), number_items, TRUE),
                          stressful = rep(c(0, 0, 1, 0, 0), 10))
  constraints <- list(list(name = 'time',
                           op = '><',
                           target = c(10, 20)),
                      list(name = 'type/depression',
                           op = '><',
                           target = c(2, 5)),
                      list(name = 'type/insomnia',
                           op = '><',
                           target = c(2, 5)),
                      list(name = 'stressful',
                           op = '<',
                           target = 3))
  constraints_formatted <- constraints_lp_format(max_n = max_n, number_items = number_items, characteristics = characteristics, constraints = constraints)
  expect_equal(constraints_formatted$errors$characteristics, "should be a data frame with number of rows equal to the number of items in the item bank")
})


test_that("constraints wrong structure", {
  number_items <- 50
  max_n <- 30
  
  #create item characteristics and constraints
  characteristics <- data.frame(time = with_random_seed(2, rnorm)(number_items),
                                type = with_random_seed(2, sample)(c('depression', 'insomnia', 'anxiety'), number_items, TRUE),
                                stressful = rep(c(0, 0, 1, 0, 0), 10),
                                extra = rep(c(1, 0, 0, 0, 0), 10))
  constraints <- data.frame(name = 2)
  
  constraints_formatted <- constraints_lp_format(max_n = max_n, number_items = number_items, characteristics = characteristics, constraints = constraints)
  expect_equal(constraints_formatted$errors$constraints_structure, "should be a list of length three lists, with elements named 'name', 'op', 'target'")
})

test_that("constraints wrong name elements", {
  number_items <- 50
  max_n <- 30
  
  characteristics <- data.frame(time = with_random_seed(2, rnorm)(number_items),
                                type = with_random_seed(2, sample)(c('depression', 'insomnia', 'anxiety'), number_items, TRUE),
                                stressful = rep(c(0, 0, 1, 0, 0), 10))
  constraints <- list(list(name = 'time',
                           op = '><',
                           target = c(10, 20)),
                      list(name = 'depression',
                           op = '><',
                           target = c(2, 5)),
                      list(name = 'insomnia',
                           op = '><',
                           target = c(2, 5)),
                      list(name = 'stressful',
                           op = '<',
                           target = 3))
  constraints_formatted <- constraints_lp_format(max_n = max_n, number_items = number_items, characteristics = characteristics, constraints = constraints)
  expect_equal(constraints_formatted$errors$constraints_name_elements, "should be defined as described in the details section")
})

test_that("constraints wrong operator elements", {
  number_items <- 50
  max_n <- 30
  
  characteristics <- data.frame(time = with_random_seed(2, rnorm)(number_items),
                                type = with_random_seed(2, sample)(c('depression', 'insomnia', 'anxiety'), number_items, TRUE),
                                stressful = rep(c(0, 0, 1, 0, 0), 10))
  constraints <- list(list(name = 'time',
                           op = c('>', '<'),
                           target = c(10, 20)),
                      list(name = 'type/depression',
                           op = c('>', '<'),
                           target = c(2, 5)),
                      list(name = 'type/insomnia',
                           op = c('>', '<'),
                           target = c(2, 5)),
                      list(name = 'stressful',
                           op = '==',
                           target = 3))
  constraints_formatted <- constraints_lp_format(max_n = max_n, number_items = number_items, characteristics = characteristics, constraints = constraints)
  expect_equal(constraints_formatted$errors$constraints_operator_elements, "should be defined as described in the details section")
})

test_that("constraints wrong target elements", {
  number_items <- 50
  max_n <- 30
  
  characteristics <- data.frame(time = with_random_seed(2, rnorm)(number_items),
                                type = with_random_seed(2, sample)(c('depression', 'insomnia', 'anxiety'), number_items, TRUE),
                                stressful = rep(c(0, 0, 1, 0, 0), 10))
  constraints <- list(list(name = 'time',
                           op = '><',
                           target = c(10, 20)),
                      list(name = 'type/depression',
                           op = '><',
                           target = list(2, 5)),
                      list(name = 'type/insomnia',
                           op = '><',
                           target = c(5)),
                      list(name = 'stressful',
                           op = '<',
                           target = "3"))
  constraints_formatted <- constraints_lp_format(max_n = max_n, number_items = number_items, characteristics = characteristics, constraints = constraints)
  expect_equal(constraints_formatted$errors$constraints_target_elements, "should be defined as described in the details section")
})
