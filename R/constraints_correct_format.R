#' Creates a constraints object.
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
#' Constraints in Shadow Tests are implemented through linear programming, for which the package \code{\link{lpSolve}} is used. 
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
#' 
#' @return constraints Constraints object, see \code{details}.
#' @importFrom lpSolve lp
#' @export
createConstraints <- function(test, characteristics = NULL, constraints = NULL) {
  result <- function() {
    # create characteristics data.frame, add administered and N.

    characteristics_numeric <- if (is.null(characteristics)) 
                                 data.frame(length = rep(1, test$items$K))
                               else
                                 as.data.frame(cbind(data.frame(length = rep(1, test$items$K)), 
                                                     get_characteristics_numeric()))
    
    
    # test stops at whichever occurs first, length stopping rule or max_n
    max_n <- ifelse(test$stop$type == 'length', min(test$stop$n, test$max_n), test$max_n)
    
    # create constraints for (max) n
    CONSTS <- data.frame(name = 'length', 
                         op = '=', 
                         target = max_n,
                         stringsAsFactors = FALSE)    
    
    ### add user constraints
    if (! is.null(constraints)){
      for(con in constraints){
        # stop if not list
        if(! is.list(con) || length(con) != 3) stop("Each constraint should be a list of three elements.")
        if(! con[[1]] %in% colnames(characteristics_numeric)) stop("Each constraint name should have a matching characteristic.")
        if(! con[[2]] %in% c("<", "=", ">", "><", "<=", ">=")) stop("Invalid operator.")
        if(! is.numeric(con[[3]])) stop("Target value must be numeric.")
        
        if (con[[2]] == "><"){
          CONSTS <- rbind(CONSTS, c(con[[1]], ">", con[[3]][1]))
          CONSTS <- rbind(CONSTS, c(con[[1]], "<", con[[3]][2]))
        } else {
          CONSTS <- rbind(CONSTS, c(con[[1]], con[[2]], con[[3]]))
        }
      }  
    }
    
    # build characteristics table for use in lpSolve
    LPCHARS <- characteristics_numeric[,CONSTS$name, drop = FALSE]
    
    # create constraints object
    constraints <- list(characteristics = characteristics_numeric, constraints = CONSTS, lp_chars = LPCHARS)
    attr(constraints, 'class') <- "ShadowCAT.constraints"
    
    constraints
  }
  
  get_characteristics_numeric <- function() {
    numeric_characteristics_list <- sapply(colnames(characteristics), 
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
    matrix(unlist(numeric_characteristics_list), 
           nrow = test$items$K, 
           dimnames = list(c(), get_names_numeric_characteristics_list(numeric_characteristics_list)))
  }
  
  get_names_numeric_characteristics_list <- function(numeric_characteristics_list) {
    unlist(sapply(numeric_characteristics_list, 
                  FUN = function(key) { colnames(key) }))
  } 
  
  validate <- function() {
    if (!is.null(characteristics) && (!is.data.frame(characteristics) || any(is.null(colnames(characteristics)))))
      add_error("Characteristics", "should be a data.frame with unique column names.")
  }
  
  invalid_result <- function() {
    list(characteristics = NA,
         constraints = NA,
         lp_chars = NA,
         errors = errors())
  }
  
  validate_and_run() 
}
