if(getRversion() >= "2.15.1")  utils::globalVariables(c("person_updated_after_new_response", "index_new_item"))
#' Returns the next item to be administered given a new response
#' 
#' Run test with a specified person, and a specified test. See initPerson and initTest.
#'  This is a simple wrapper to call the right methods, options should be defined in the test object.
#'
#' Details
#' @param new_response new response from respondent
#' @param person initialized person
#' @param test test object
#' @return index next item; when test is finished, person object
#' @export
ShadowCAT_roqua <- function(new_response, person, test) {
  result <- function() {
    if (is.null(new_response)) { # first iteration: no responses given yet
      assign("person_updated_after_new_response", person, envir = .GlobalEnv)
      assign("index_new_item", next_item(person, test), envir = .GlobalEnv)
      return(index_new_item)
    } 
    
    assign("person_updated_after_new_response", update_person_estimate(person_updated_after_new_response), envir = .GlobalEnv)
    if (!stop_test(person_updated_after_new_response, test)) {
      assign("index_new_item", next_item(person_updated_after_new_response, test), envir = .GlobalEnv)
      index_new_item
    }
    else {
      person_updated_after_new_response
    }
  }
  
  # if inititial items have been administered (so we are in the CAT phase), update person estimate after each newly answered item
  update_person_estimate <- function(person) {
    person$responses <- c(person$responses, new_response)
    person$administered <- c(person$administered, index_new_item)
    person$available <- person$available[-which(person$available %in% index_new_item)]
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


