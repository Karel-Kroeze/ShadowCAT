#' Done
#' 
#' Control function to check if the test is completed.
#' 
#' @param test
#' @param person
#' @return boolean completed
#' @export
stop_test <- function(person, test) {
  if (test$stop$type == 'length') {
    if (length(person$responses) >= test$stop$n) return(TRUE)
  }
  
  return(FALSE)
}


#' Next item
#' 
#' Control function to select the next item, takes account of starting / stopping conditions.
#' 
#' @param person
#' @param test
#' @return integer item index
#' @export
next_item <- function(person, test) {
  if (length(person$responses) < test$start$n) {
    # Pre-CAT
    if (test$start$type == 'random'){
      out <- sample(person$available, 1)
    }
    
    if (test$start$type == 'fixed'){
      out <- test$start$indices[length(person$responses) + 1]
    }
    
  } else {
    # CAT 
    out <- best_item(person, test)
  }
  
  return(out)
}

#' ShadowCAT
#' 
#' Run test with a specified person, and a specified test. See initPerson and initTest.
#'  This is a simple wrapper to call the right methods, options should be defined in the test object.
#'
#' Details
#' @param person
#' @param test
#' @return person
#' @export
ShadowCAT <- function(person, test, verbose = FALSE) {
  ## Start CAT
  if (verbose > 0) cat("\n")
  while(! stop_test(person, test)) {
    # respones hack to always give 1,0,1,0 etc response pattern.
    # person$responses <- rep(c(1,0), length.out = length(person$responses))
    # TODO: remove 
    
    person <- answer(person, test, next_item(person, test))
    if (length(person$responses) > test$start$n) person <- estimate(person, test, check.analyticals = TRUE)
    if (verbose > 1) cat("\r", paste0(round(person$estimate, 2), collapse = ', '), " | ", paste0(round(diag(attr(person$estimate,'variance')), 2), collapse = ', '), '         ')
  }
  
  if (verbose > 0) cat("\r", paste0(round(person$estimate, 2), collapse = ', '), " | ", paste0(round(diag(attr(person$estimate,'variance')), 2), collapse = ', '), '         ', "\n")
  return(person)
}
