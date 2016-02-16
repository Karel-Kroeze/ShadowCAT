#' Creates a list with characteristics and constraints which can be used in lp function from lpSolve package.
#' 
#' @section Constraint specification:
#' \code{constraints} should be specified as a list of constraints, each constraint is a list with three named values;
#' \code{name} the column name of the characteristic this constraint applies to. For categorical characteristics the level should be specified as \code{name/value}.
#' \code{op} the logoical operator to be used. Valid options are "<", "=", ">" and "><".
#' \code{target} the target value, numeric. If the operator is "><", this should be a length two vector in between which the target should fall. 
#' 
#' @section Return object:
#' The constraints and characteristics in lp format will be stored within a list. Setting the constraints and characteristics
#' arguments to NULL results in a list with only the test length constraint.
#' Constraints in Shadow Tests are implemented through linear programming, for which the package lpSolve is used. 
#' The returned constraints object is a list with two named elements; lp_constraints and lp_chars. The lp_constraints 
#' and lp_chars are set up to work with lpSolve, and should not be manually edited.
#' 
#' @section Note about length:
#' Note that the maximum test length is always included as an additional constraint.
#' 
#' @examples
#' max_n <- 30
#' number_items <- 50
#' 
#' # set up some dummy characteristics.
#' content <- sample(c('algebra','physics','calculus'), number_items, TRUE)
#' time <- rnorm(number_items)
#' exclusive <- rep(0, number_items)
#' exclusive[sample(number_items, 4)] <- 1
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
#' # get list of characteristics and constraintrs in lp format
#' chars_constraints_lp <- constraints_lp_format(max_n, number_items, characteristics, constraints)
#' 
#' @param max_n test length at which testing should stop
#' @param number_items number of items available in the item bank
#' @param characteristics \code{data.frame} with characteristics, one row per item, one column per characteristic.
#' @param constraints \code{list} of constraints, see \code{details}.
#' @return list containing characteristics and constraints in lp format; 
#' the maximum test length is always included as an additional constraint; see \code{details}.
#' @export
constraints_lp_format <- function(max_n, number_items, characteristics = NULL, constraints = NULL) {
  result <- function() {
    characteristics_numeric <- get_characteristics_numeric()
    constraints_lp <- get_constraints_lp()
    characteristics_numeric_lp <- characteristics_numeric[,constraints_lp$name, drop = FALSE]
    
    list(lp_constraints = constraints_lp, 
         lp_chars = characteristics_numeric_lp)
  }
  
  get_characteristics_numeric <- function() {
    if (is.null(characteristics)) 
      return(data.frame(length = rep(1, number_items)))
    
    numeric_characteristics_list <- lapply(colnames(characteristics), 
                                           function(key) {   
                                             if (is.character(characteristics[[key]]) || is.factor(characteristics[[key]])) {
                                               dummy_matrix <- sapply(unique(characteristics[[key]]), categorical_to_dummy, categorical_vector = characteristics[[key]])
                                               colnames(dummy_matrix) <- paste(key, unique(characteristics[[key]]), sep = '/')
                                               dummy_matrix 
                                            }  
                                             else
                                               matrix(characteristics[[key]], ncol = 1, dimnames = list(c(), key))
                                            } 
                                    )

    numeric_characteristics <- unlist_numeric_characteristics(numeric_characteristics_list)
    
    cbind(data.frame(length = rep(1, number_items)),
          numeric_characteristics)
  }
 
  
  get_names_numeric_characteristics_list <- function(numeric_characteristics_list) {
    unlist(sapply(numeric_characteristics_list, function(key) { colnames(key) }))
  }
  
  get_constraints_lp <- function() {
    if (is.null(constraints))
      return(data.frame(name = 'length', op = '=', target = max_n, stringsAsFactors = FALSE))
    
    constraints_lp_format_list <- lapply(constraints, function(constraint) { 
                                                        if (constraint$op == "><")  
                                                          rbind(c(constraint$name, ">", constraint$target[1]), 
                                                                c(constraint$name, "<", constraint$target[2]))
                                                        else
                                                          c(constraint$name, constraint$op, constraint$target)
                                                       }) 

    constraints_lp_format <- unlist_and_name_constraints_lp_format(constraints_lp_format_list)    

    as.data.frame(rbind(data.frame(name = 'length', op = '=', target = max_n, stringsAsFactors = FALSE),
                        constraints_lp_format))
  }
  
  unlist_numeric_characteristics <- function(numeric_characteristics_list) {
    if (is.list(numeric_characteristics_list))
      do.call(cbind, numeric_characteristics_list)
    else
      numeric_characteristics_list
  }
  
  unlist_and_name_constraints_lp_format <- function(constraints_lp_format_list) {
    if (is.list(constraints_lp_format_list))
      constraints_lp_format <- do.call(rbind, constraints_lp_format_list) 
    else
      constraints_lp_format <- constraints_lp_format_list
    colnames(constraints_lp_format) <- c("name", "op", "target")
    constraints_lp_format
  }
  
  validate <- function() {
    if (is.null(characteristics) && !is.null(constraints))
      return(add_error("characteristics", "is missing while constraints is defined"))
    if (!is.null(characteristics) && is.null(constraints))
      return(add_error("constraints", "is missing while characteristics is defined"))
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
