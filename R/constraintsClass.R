#' Convenience wrapper to add constraitns to a test object for use with ShadowCAT.
#' 
#' @param test ShadowCAT test object
#' @param characteristics data.frame with characteristics, one row per item, one column per characteristic.
#' @param constraints list of constraints, each constraint is a list with three named values;
#' \code{name} the column name of the characteristic this constraint applies to. For categorical characteristics the level should be specified as \code{name/value}.
#' \code{op} the logoical operator to be used. Valid options are "<", "=", ">" and "><".
#' \code{target} the target value, numeric. If the operator is "<>", this should be a length two vector in between which the target should fall.
#' @return test ShadowCAT test object, including the constraints set.
#' @export
createConstraints <- function(test, characteristics = NULL, constraints = NULL) {
  # create characteristics data.frame, add administered and N.
  CHARS <- data.frame(length = rep(1, test$items$K))
  
  # create name vector
  NAMES <- c('length')
  
  # create constraints for administered and N
  CONSTS <- data.frame(name = 'length', 
                       op = '=', 
                       target = test$stop$n,
                       stringsAsFactors = FALSE)
  
  # add provided characteristics.
  if (!is.null(characteristics)){
    # stop if not named data.frame
    if (!is.data.frame(characteristics) || any(is.null(colnames(characteristics)))) stop("Characteristics should be a data.frame with unique column names.")
    for (name in colnames(characteristics)){
      char <- characteristics[[name]]
      
      # character vectors are much easier to work with than factors.
      if (is.factor(char)) char <- as.character(char)    
      
      # add given characteristics to characteristics
      if (is.character(char)) {    
        # add each level as a binary dummy
        for (value in unique(char)){
          CHARS <- cbind(CHARS, as.numeric(char == value))
          NAMES <- c(NAMES, paste(name,value,sep='/'))
        }
      } else {
        CHARS <- cbind(CHARS, char)
        NAMES <- c(NAMES, name)
      }
    }
    
    # name characteristics
    colnames(CHARS) <- NAMES
  }
  
  
  ### add user constraints
  if (! is.null(constraints)){
    for(con in constraints){
      # stop if not list
      if(! is.list(con) || length(con) != 3) stop("Each constraint should be a list of three elements.")
      if(! con[[1]] %in% NAMES) stop("Each constraint name should have a matching characteristic.")
      if(! con[[2]] %in% c("<", "=", ">", "><")) stop("Invalid operator.")
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
  LPCHARS <- CHARS[,CONSTS$name, drop = FALSE]
  
  # create constraints object
  constraints <- list(characteristics = CHARS, constraints = CONSTS, lp_chars = LPCHARS)
  attr(constraints, 'class') <- "ShadowCAT.constraints"
  
  # attach and return
  test$constraints = constraints
  return(test)
}
