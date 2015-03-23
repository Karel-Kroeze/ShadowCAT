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
ShadowCAT <- function(person, test) {
  
  ### Pre-CAT
  if (test$start$type == 'random'){
    start_items <- sample(person$available, test$start$n)
    person <- answer(person, test, start_items)
  }
  
  if (test$start$type == 'fixed'){
    person <- answer(person, test, test$start$indeces)
  }
  
  
  ### Define stopping rule(s)
  go_on <- function(test, person) {
    if (test$stop$type == 'length') {
      if (length(person$responses) > test$stop$n) return(FALSE)
    }
    
    return(TRUE)
  }
  
  
  ## Start CAT
  while(go_on(test, person)) {
    # respones hack to always give 1,0,1,0 etc response pattern.
    # W <- length(person$responses)
    # person$responses <- rep(c(1,0), length.out = W)
    
    person <- estimate(person, test, check.analyticals = FALSE)
    person <- answer(person, test, next_item(test, person))
  }
  
  return(person)
}