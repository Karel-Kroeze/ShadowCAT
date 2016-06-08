#' Template to structure and test larger methods
#'
#' First validate is called, if any errors are added, the value of invalid_result is returned, else the value of result.
#' Loads the following methods into the environment of the calling method:
#'   add_error(key, value)
#'   errors() #returns list(key=value, key2=value2)
#'   validate_and_runner() # Used internally, contains the logic of functions to call.
#' @examples
#' two_times_two <- function(two, also_two) {
#'   result <- function() helper()
#'
#'   helper <- function() two * also_two
#'
#'   validate <- function() {
#'     if(two != 2) add_error('two', 'should be 2')
#'     if(also_two != 2) add_error('also_two', 'should also be two')
#'   }
#'
#'   invalid_result <- function() {
#'     list(result=NA, errors=errors())
#'   }
#'
#'   ShadowCAT:::validate_and_run()
#' }
#' two_times_two(2,2) == 4 || stop('not 4')
#' two_times_two(3,2)$errors$two == 'should be 2' || stop('wrong result')
validate_and_run <- function() {
  .errors <- list()
  
  add_error <- function(key, value=TRUE) { .errors[key] <<- value }
  
  errors <- function() { .errors }
  
  # Call validate if exists.
  # Return value of result if valid, of invalid_result if not.
  # Return reference to result function when in test mode.
  validate_and_runner <- function() {
    if (exists('validate', parent.frame(), inherits = FALSE))
      do.call('validate', list(), envir=parent.frame())
    
    #if (get0('test_inner_functions', envir = parent.frame(n=2), inherits = FALSE, ifnotfound = FALSE))
    if (exists('test_inner_functions', envir=parent.frame(n=2), inherits = FALSE))
      get('result', parent.frame())
    else if (length(errors()) == 0)
      do.call('result', list(), envir=parent.frame())
    else
      do.call('invalid_result', list(), envir=parent.frame())
  }
  
  # copy all variables to the calling function.
  for(n in ls(environment())) assign(n, get(n, environment()), parent.frame())
  
  do.call('validate_and_runner', list(), envir=parent.frame())
}

# Empty definitions so the R check doesn't complain. (Seems to be the accepted solution in de community. Still feel dirty.)
add_error <- function(key, value=TRUE){}
errors <- function(){}

#' Returns the result function instead of the result of the result function
#' @param fn The function to get the function result of
#' @param ... Arguments for fn
#' @examples
#' fn <- function(a) {
#'   result <- function() helper()
#'   helper <- function() 2*a
#'   ShadowCAT:::validate_and_run()
#' }
#' tst <- ShadowCAT:::validate_and_run.test(fn, 4)
#' helper <- get('helper', environment(tst))
#' helper() == 8 || stop('error')
validate_and_run.test <- function(fn, ...) {
  test_inner_functions <- TRUE
  fn(...)
}
