#' ShadowCAT
#' 
#' Run test with a specified person, and a specified test. See initPerson and initTest.
#'  This is a simple wrapper to call the right methods, options should be defined in the test object.
#'
#' Details
#' @param person last update of person
#' @param test test object
#' @param new_response new responses from respondent
#' @param indeces indeces of newly administered items 
#' @return next item and updated person object; when test is finished, new item is NULL
#' @export
ShadowCAT_roqua <- function(person, test, new_responses, indeces) {
  result <- function() {
    if (is.null(new_responses)) # first iteration: no responses given yet
      return(list(next_item = next_item(person, test),
                  person = person))
    
    person <- update_person_estimate(person)
    if (!stop_test(person, test))
      list(next_item = next_item(person, test),
           person = person)
    else
      list(next_item = NULL,
           person = person)
  }
  
  # if inititial items have been administered (so we are in the CAT phase), update person estimate after each newly answered item
  update_person_estimate <- function(person) {
    person$responses <- c(person$responses, new_responses)
    person$administered <- c(person$administered, indeces)
    person$available <- person$available[-which(person$available %in% indeces)]
    if (length(person$responses) > test$start$n) 
      estimate(person, test)
    else
      person
  }
  
  validate <- function() {
    if (is.null(person))
      add_error("person", "is missing")
    if (is.null(test))
      add_error("test", "is missing")
  }
  
  invalid_result <- function() {
    list(errors = errors())
  }
  
  validate_and_run()
}


#' ShadowCAT Karel Kroeze
#' 
#' Run test with a specified person, and a specified test. See initPerson and initTest.
#'  This is a simple wrapper to call the right methods, options should be defined in the test object.
#'
#' Details
#' @param person
#' @param test
#' @param verbose if larger than 0, print estimate and variance of estimate after test is finished. If larger than 1, also print estimate and variance of estimate at each iteration
#' @return person
#' @export
ShadowCAT <- function(person, test, verbose = FALSE) {
  ## Start CAT
  if (verbose > 0) cat("\n")
  
  while(!stop_test(person, test)) {
    # update person with new answer
    person <- answer(person, test, indeces = next_item(person, test))
    
    # if inititial items have been administered (so we are in the CAT phase), update person estimate after each newly answered item
    if (length(person$responses) > test$start$n) 
      person <- estimate(person, test)
    
    if (verbose > 1) 
      cat("\r", paste0(round(person$estimate, 2), collapse = ', '), " | ", paste0(round(diag(attr(person$estimate,'variance')), 2), collapse = ', '), '         ')
  }
  
  if (verbose > 0) cat("\r", paste0(round(person$estimate, 2), collapse = ', '), " | ", paste0(round(diag(attr(person$estimate,'variance')), 2), collapse = ', '), '         ', "\n")
  
  invisible(person)
}


