if(getRversion() >= "2.15.1")  utils::globalVariables(c("person_updated_after_new_response", "index_new_item"))
#' Returns the index of the next item to be administered given a new response, and updates the global variables person_updated_after_new_response and index_new_item.
#' When test is finished, person_updated_after_new_response is returned 
#'
#' @param new_response new response from respondent, should be initialized with NULL
#' @param prior covariance matrix of the multi variate normal prior for theta; mean vector is fixed at zero; only used when estimator type is MAP or EAP, but at this point should always be defined
#' #' note that this prior should be a square matrix with number of rows and columns equal to the number of dimensions; values on the diagonal should be larger than 1
#' @param model String, one of '3PLM', 'GPCM', 'SM' or 'GRM', for the three-parameter logistic, generalized partial credit, sequential or graded response model respectively.
#' @param alpha Matrix of alpha paramteres, one column per dimension, one row per item. Note that so called within-dimensional models still use an alpha matrix, they simply 
#' have only one non-zero loading per item.
#' @param beta Matrix of beta parameters, one column per item step, one row per item. Note that ShadowCAT expects response categories to be sequential,
#' and without gaps. That is, the weight parameter in the GPCM model is assumed to be sequential, and equal to the position of the 'location' of the beta parameter in the Beta matrix.
#' The matrix will have a number of columns equal to the largest number of response categories, items with fewer response categories should be 
#' right-padded with \code{NA}. \code{NA} values between response categories are not allowed, and will lead to errors.
#' More flexibility in Beta parameters might be added in future versions.
#' @param guessing vector of guessing parameters per item. Optionally used in 3PLM model, ignored for all others.
#' @param eta Matrix of location parameters, optionally used in GPCM model, ignored for all others.
#' @param start_items items that are shown to the patient before adaptive proces starts; one of
#' list(type = 'random', n)
#' list(type = 'fixed', indices, n)
#' list(type = 'randomByDimension', nByDimension, n)
#' where n = total number of initial items, indices = vector of initial item indeces, 
#' nByDimension = scalar of number of initial items per dimension, or vector with number of initial items for each dimension
#' @param stop_test rule for when to stop providing new items to patient; one of
#' list(type = 'length', n = ...);
#' list(type = 'variance', target = ..., n = ...)
#' where n = test length at which testing should stop (even if target has not been reached yet in case of variance stopping rule), 
#' target = vector of maximum acceptable variances per dimension
#' @param estimator type of estimator to be used, one of "MAP" (Maximum a posteriori estimation), "EAP" (Expected A Posteriori Estimation), or "ML" (maximum likelihood)
#' @param information_summary called "objective" by Kroeze; how to summarize information; one of
#' "D" = determinant: compute determinant(info_sofar_QxQ + info_QxQ_k) for each yet available item k
#' "PD" = posterior determinant: compute determinant(info_sofar_QxQ_plus_prior + info_QxQ_k) for each yet available item k
#' "A" = trace: compute trace((info_sofar_QxQ + info_QxQ_k) for each yet available item k
#' "PA" = posterior trace: compute trace(info_sofar_QxQ_plus_prior + info_QxQ_k) for each yet available item k
#' "PEKL" = compute Posterior expected Kullback-Leibler Information
#' @param item_selection selection criterion; one of "MI" (maximum information) or "Shadow" (maximum information and take constraints into account)
#' @param constraints list with constraints and characteristics
#' constraints should be specified as a list of constraints, each constraint is a list with three named values;
#' name: the column name of the characteristic this constraint applies to. For categorical characteristics the level should be specified as name/value.
#' op: the logical operator to be used. Valid options are "<", "=", ">" and "><".
#' target: the target value, numeric. If the operator is "><", this should be a length two vector in between which the target should fall.
#' characteristics should be a data.frame with characteristics, one row per item, one column per characteristic.
#' See constraints_correct_format() for details
#' @param lowerbound vector with lower bounds for theta per dimension; estimated theta values smaller than the lowerbound values are truncated to the lowerbound values 
#' @param upperbound vector with upper bounds for theta per dimension; estimated theta values larger than the upperbound values are truncated to the upperbound values
#' @return index next item; when test is finished, "stop_test"
#' @export
shadowcat_roqua <- function(new_response, prior, model, alpha, beta, guessing, eta = NULL, start_items, stop_test, estimator, information_summary, item_selection = "MI", constraints = NULL, lowerbound = rep(-3, ncol(alpha)), upperbound = rep(3, ncol(alpha))) {
  item_characteristics_shadowcat_format <- initItembank(model = model, alpha = alpha, beta = beta, guessing = guessing, eta = eta, silent = TRUE)
  
  person <- initPerson(items = item_characteristics_shadowcat_format, 
                       prior = prior)
  
  test <- initTest(items = item_characteristics_shadowcat_format, 
                   start = start_items, 
                   stop = stop_test,
                   max_n = stop_test$n,
                   estimator = estimator,
                   objective = information_summary,
                   selection = item_selection,
                   constraints = constraints,
                   exposure = NULL,
                   lowerBound = lowerbound,
                   upperBound = upperbound)
   
  result <- function() {
    if (is.null(new_response)) { # first iteration: no responses given yet
      assign("person_updated_after_new_response", person, envir = .GlobalEnv)
      assign("index_new_item", next_item(person, test), envir = .GlobalEnv)
      return(list(index_new_item = index_new_item,
                  person_updated_after_new_response = person_updated_after_new_response))
    } 
    
    assign("person_updated_after_new_response", update_person_estimate(person_updated_after_new_response), envir = .GlobalEnv)
    if (!stop_test(person_updated_after_new_response, test)) {
      assign("index_new_item", next_item(person_updated_after_new_response, test), envir = .GlobalEnv)
      list(index_new_item = index_new_item,
           person_updated_after_new_response = person_updated_after_new_response)
    }
    else {
      list(index_new_item = "stop_test",
           person_updated_after_new_response = person_updated_after_new_response)
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


