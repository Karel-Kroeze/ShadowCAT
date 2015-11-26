#' Creates a constraints object which can be used in lp function from lpSolve package.
#' 
#' Convenience wrapper to create a constraints object, which can be added to a \code{test} object.
#' 
#' @section Set constraints from \code{\link{initTest}}:
#' Note that constraints can also be set with the initTest call. Simply include a list with named elements 
#' characteristics and constraints as described above to the initTest function call.
#' 
#' @section Constraint specification:
#' \code{constraints} should be specified as a list of constraints, each constraint is a list with three named values;
#' \code{name} the column name of the characteristic this constraint applies to. For categorical characteristics the level should be specified as \code{name/value}.
#' \code{op} the logoical operator to be used. Valid options are "<", "=", ">" and "><".
#' \code{target} the target value, numeric. If the operator is "><", this should be a length two vector in between which the target should fall. 
#' 
#' @section Constraints in test object:
#' The specified constraints will be stored within a list, which can/should be saved in the test object (test$constraints). 
#' Constraints in Shadow Tests are implemented through linear programming, for which the package lpSolve is used. 
#' The constraints object is a list with three named elements; characterstics, constraints and lp_chars. Characteristics is a copy of the argument given, constraints and lp_chars
#' are set up to work with lpSolve, and should not be manually edited.
#' 
#' @section Note about length:
#' Note that the maximum test length (either max_length or the length stopping parameter) is always included as an additional constraint.
#' 
#' @examples
#' # set up a simple itembank and test.
#' items <- createTestBank("GPCM")
#' test <- initTest(items, selection = "Shadow" , objective = "PEKL")
#' 
#' # set up some dummy characteristics.
#' content <- sample(c('algebra','physics','calculus'), items$K, TRUE)
#' time <- rnorm(items$K)
#' exclusive <- rep(0, items$K)
#' exclusive[sample(items$K, 4)] <- 1
#' 
#' # bind them in a data.fame
#' characteristics <- data.frame(content, time, exclusive)
#' 
#' # set up the constraints
#' constraints <- list(
#'   list(name = 'content/algebra',
#'        op = '><',
#'        target = c(5,10)),
#'   list(name = 'content/physics',
#'        op = '><',
#'        target = c(2,5)),
#'   list(name = 'time',
#'        op = '<',
#'        target = 20),
#'   list(name = 'exclusive',
#'        op = '<',
#'        target = 2))
#' 
#' # update the test object
#' test$constraints <- createConstraints(test, characteristics, constraints)
#' 
#' # or do it all at once;
#' test2 <- initTest(items, constraints = list(characteristics = characteristics, constraints = constraints))
#' 
#' # results are identical (initTest uses createConstraints internally);
#' all.equal(test$constraints, test2$constraints)
#' 
#' @param test Test object, see \code{\link{initTest}}
#' @param characteristics \code{data.frame} with characteristics, one row per item, one column per characteristic.
#' @param constraints \code{list} of constraints, see \code{details}.
#' @return Constraints object, see \code{details}.
#' @export
createConstraints <- function(test, characteristics = NULL, constraints = NULL) {
  result <- function() {
    # test stops at whichever occurs first, length stopping rule or max_n
    max_n <- ifelse(test$stop$type == 'length', min(test$stop$n, test$max_n), test$max_n)
    
    characteristics_numeric <- (if (is.null(characteristics)) 
                                  data.frame(length = rep(1, test$items$K)) 
                                else
                                  as.data.frame(cbind(data.frame(length = rep(1, test$items$K)), 
                                                     get_characteristics_numeric())))
    
    constraints_lp <- (if (is.null(constraints))
                         data.frame(name = 'length', op = '=', target = max_n, stringsAsFactors = FALSE)
                       else
                         as.data.frame(rbind(data.frame(name = 'length', op = '=', target = max_n, stringsAsFactors = FALSE),
                                             get_constraints_lp())))

    characteristics_numeric_lp <- characteristics_numeric[,constraints_lp$name, drop = FALSE]
    
    characteristics_and_constraints_lp <- list(characteristics = characteristics_numeric, constraints = constraints_lp, lp_chars = characteristics_numeric_lp)
    attr(characteristics_and_constraints_lp, 'class') <- "ShadowCAT.constraints"
    
    characteristics_and_constraints_lp
  }
  
  get_characteristics_numeric <- function() {
    numeric_characteristics_list <- lapply(colnames(characteristics), 
                                           FUN = function(key) {   
                                            if (is.character(characteristics[[key]]) || is.factor(characteristics[[key]])) {
                                              dummy_matrix <- sapply(unique(characteristics[[key]]), FUN = categorical_to_dummy, categorical_vector = characteristics[[key]])
                                              colnames(dummy_matrix) <- paste(key, unique(characteristics[[key]]), sep = '/')
                                              dummy_matrix 
                                            }  
                                            else
                                              matrix(characteristics[[key]], ncol = 1, dimnames = list(c(), key))
                                            } 
                                    )

    if (is.list(numeric_characteristics_list))
      do.call(cbind, numeric_characteristics_list)
    else
      numeric_characteristics_list 
  }
  
  get_names_numeric_characteristics_list <- function(numeric_characteristics_list) {
    unlist(sapply(numeric_characteristics_list, 
                  FUN = function(key) { colnames(key) }))
  }
  
  get_constraints_lp <- function() {
    constraints_lp_format_list <- lapply(constraints, FUN = function(constraint) { if (constraint$op == "><")  
                                                                                     rbind(c(constraint$name, ">", constraint$target[1]), 
                                                                                           c(constraint$name, "<", constraint$target[2]))
                                                                                   else
                                                                                     c(constraint$name, constraint$op, constraint$target)
                                                                                 }) 

    constraints_lp_format <- ( if (is.list(constraints_lp_format_list))
                                 do.call(rbind, constraints_lp_format_list) 
                               else
                                 constraints_lp_format_list )
      
    colnames(constraints_lp_format) <- c("name", "op", "target")
    constraints_lp_format
  }
  
  validate <- function() {
    if (!is.null(characteristics) && (!is.data.frame(characteristics) || any(is.null(colnames(characteristics)))))
      return(add_error("characteristics", "should be a data.frame with unique column names."))
    if (!is.null(constraints) && any(sapply(constraints, FUN = function(constraint){ !is.list(constraint) || length(constraint) != 3 || names(constraint) != c("name", "op", "target") })))
      return(add_error("constraints", "is of invalid structure, should be list of lists, each sublist containing three elements called name, op, and target."))
    if (!is.null(characteristics) && !is.null(constraints) && any(sapply(constraints, 
                                                                         FUN = function(constraint){ constraint$name %not_in% unlist(sapply(colnames(characteristics), 
                                                                                                                                            FUN = function(key) { 
                                                                                                                                              if (is.character(characteristics[[key]]) || is.factor(characteristics[[key]]))
                                                                                                                                                paste(key, unique(characteristics[[key]]), sep = '/')
                                                                                                                                              else
                                                                                                                                                key })) 
                                                                                                     })))
      add_error("constraints_names", "should be contained in names characteristics")
    if (!is.null(constraints) && any(sapply(constraints, FUN = function(constraint){ constraint$op %not_in% c("<", "=", ">", "><", "<=", ">=") })))
      add_error("constraints_operators", "containts invalid operator(s)")
    if (!is.null(constraints) && any(sapply(constraints, FUN = function(constraint){ !is.numeric(constraint$target) })))
      add_error("constraints_targets", "must be numeric")
  }
  
  invalid_result <- function() {
    list(characteristics = NA,
         constraints = NA,
         lp_chars = NA,
         errors = errors())
  }
  
  validate_and_run() 
}
