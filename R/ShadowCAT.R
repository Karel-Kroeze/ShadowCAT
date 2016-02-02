#' Done
#' 
#' Control function to check if the test is completed.
#' 
#' @param test
#' @param person
#' @return boolean completed
#' @export
stop_test <- function(person, test) {
  # current test duration
  l <- length(person$responses)
  
  # throw stopping criteria at the test until we can make something stick
  # max length rule, for use in combination with non-length based stopping rules (<= 0 never stops for max length!)
  if (test$max_n > 0 && l >= test$max_n) return(TRUE)
  
  # length stop rule ()
  if (test$stop$type == 'length') {
    if (l >= test$stop$n) return(TRUE)
  }
  
  # variance stop rule (target (vector of length Q with 'threshold' variances))
  if (test$stop$type == 'variance') {
    variance <- diag(attr(person$estimate, "variance"))
    if (all(variance < test$stop$target)) return(TRUE)
  }
  
  # if nothing stuck, don't stop!
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
  l <- length(person$responses)
  if (l < test$start$n) {
    # Pre-CAT
    # utter random selection
    if (test$start$type == 'random'){
      out <- sample(person$available, 1)
    }
    
    # set a pre-defined list of starting items (indices)
    if (test$start$type == 'fixed'){
      out <- test$start$indices[length(person$responses) + 1]
    }
    
    # picks nByDimension starting items per dimension (or n_i if nByDimension is a length Q vector), assumes between models, 
    # if any item has a non-zero loading on a dimension, it is considered to be part of that dimension. 
    # they CAN overlap, which may cause unwanted side effects, and in within models the result is identical to 'normal' random starting.
    if (test$start$type == 'randomByDimension'){
      # check for vector or scalar argument, make vector
      if (length(test$start$nByDimension) == 1) test$start$nByDimension <- rep(test$start$nByDimension, test$items$Q)
      
      # validate n and nByDimension
      if (test$start$n != sum(test$start$nByDimension))
        stop("Total length of start phase and sum of length per dimension do not match (n != sum(nByDimension)")
      
      # create a design matrix of item loadings
      design <- test$items$pars$alpha > 0
      
      # loop over dims, see which we're in
      for (q in 1:test$items$Q) {
          if (l < sum(test$start$nByDimension[1:q]))
            break 
      }
      
      # get items for this dimension, cull non-applicable, and select at random.
      loading <- (1:test$items$K)[design[,q]]
      applicable <- intersect(loading, person$available)
      out <- sample(applicable, 1)
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
#' @param person Person object
#' @param test Test object
#' @param responses Response vector (if set, responses will be taken from this vector, if NULL, responses will be simulated)
#' @param verbose Print extra debug output? ( default = FALSE, 1/TRUE = basic output, 2 = detailed output )
#' @return person Final person object
#' @export
ShadowCAT <- function(person, test, verbose = FALSE, responses = NULL ) {
  ## Start CAT
  if (verbose > 0) cat("\n")
  while(! stop_test(person, test)) {
    next_item_index <- next_item(person, test)
    if (is.null(responses)){
      person <- answer(person, test, next_item_index)
    } else {
      person$administered <- c(person$administered, next_item_index)
      person$responses <- c(person$responses, responses[next_item_index])
      person$available <- person$available[-which(person$available %in% person$administered)]
    }
    if (length(person$responses) > test$start$n) person <- estimate(person, test)
    if (verbose > 1) cat("\r", paste0(round(person$estimate, 2), collapse = ', '), " | ", paste0(round(diag(attr(person$estimate,'variance')), 2), collapse = ', '), '         ')
  }
  
  if (verbose > 0) cat("\r", paste0(round(person$estimate, 2), collapse = ', '), " | ", paste0(round(diag(attr(person$estimate,'variance')), 2), collapse = ', '), '         ', "\n")
  return(invisible(person))
}
