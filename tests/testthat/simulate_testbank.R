#' Simulate testbank
#' 
#' Quick and simple itembanks for testing purposes.
#' @param model String, one of '3PLM', 'GPCM', 'SM' or 'GRM', for the three-parameter logistic, generalized partial credit, sequential or graded response model respectively.
#' @param number_items number of items
#' @param number_dimensions number of dimensions
#' @param number_itemsteps number of item steps (number of categories minus 1); forced to 1 if model is 3PLM
#' @param items_load_one_dimension if TRUE, force items to load on one dimension each
#' @param return_testbank_properties if FALSE, a list of alpha and beta is returned; if TRUE, substract_testbank_properties() is applied on
#' the simulated testbank and the result of this is returned 
#' @param varying_number_item_steps if TRUE, some item steps are set to NA; in this case number_itemsteps is the maximum number of itemsteps
#' @return simulated itembank
simulate_testbank <- function(model, number_items = 50, number_dimensions = 1, number_itemsteps = 4, items_load_one_dimension = FALSE, return_testbank_properties = TRUE, varying_number_item_steps = FALSE){
  result <- function() {
    number_items <- get_number_items()
    number_itemsteps <- get_number_itemsteps()
    
    if (return_testbank_properties)
      substract_testbank_properties(model, alpha, beta, silent = TRUE)
    else
      list(alpha = alpha,
           beta = beta)
  }
  
  get_number_itemsteps <- function() {
    if (model == "3PLM") 
      1
    else
      number_itemsteps
  }
  
  get_number_items <- function() {
    # make sure the number of items is divisible by the number of dimensions
    ceiling(number_items / number_dimensions) * number_dimensions
  }
  
  get_alpha <- function() {
    alpha <- matrix(runif(number_items * number_dimensions, .3, 1.5), number_items, number_dimensions)
    if (items_load_one_dimension) {
      set <- number_items / number_dimensions
      for (dimension in 1:number_dimensions)
        alpha[((dimension - 1) * set + 1):(dimension * set), (1:number_dimensions)[-i]] <- 0
      alpha
    }
    else {
      alpha
    }
  }
  
  get_beta <- function() {
    beta_one_itemstep <- matrix(rnorm(number_items), number_items, 1)
    if (number_itemsteps == 1) {
      beta_one_itemstep
    }
    else {
      # spread polytomous items cats -2 to +2.
      spread <- seq(-2, 2, length.out = number_itemsteps)
      beta_multiple_itemsteps <- ( if (model == "GPCM") 
                                 row_sumsum(t(apply(beta_one_itemstep, 1, function(x) x + spread))) 
                                else
                                  t(apply(beta_one_itemstep, 1, function(x) x + spread)) )
      if (varying_number_item_steps)
        induce_varying_number_item_steps(beta_multiple_itemsteps)
      else
        beta_multiple_itemsteps
    }  
  }
  
  induce_varying_number_item_steps <- function(beta_multiple_itemsteps) {
    beta_multiple_itemsteps[sample(1:number_items, ceiling(number_items / 10)), ncol(beta_multiple_itemsteps)] <- NA
    if (ncol(beta_multiple_itemsteps) > 2)
      beta_multiple_itemsteps[sample(1:number_items, ceiling(number_items / 10)), (ncol(beta_multiple_itemsteps) - 1):ncol(beta_multiple_itemsteps)] <- NA
  }
}

#' Substract extra information about testbank
#' 
#' For all models, alpha is required (but may be set to a 1 vector/matrix for Rasch-type models).
#' Beta is required for all models, but in the case of GPCM, it will be computed from Eta if Eta is given and Beta omitted.
#' Guessing is only valid for the 3PLM model, and may safely be omitted (which implies a guessing parameter of 0 on all items).
#' 
#' Note that the GPCM expects Beta parameters to be category bounds rather than location parameters (Eta). That is, Beta_i = sum(Eta_1, ..., Eta_i).
#' 
#' TODO: model references.
#' 
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
#' @param silent if TRUE, a summary of the item bank properties is printed
#' @return list containing testbank and its properties
substract_testbank_properties <- function(model = '3PLM', alpha = NULL, beta = NULL, guessing = NULL, eta = NULL, silent = FALSE){  
  result <- function() {  
    # define output, list
    # I will change Q (number of dimensions), K (number of items), M (number of item steps; number of categories 
    # minus 1), and m (in pars: the number of item steps per item) into appropriate names later, since this will affect other functions
    item_bank <- list(pars = get_parameters_list(), 
                      Q = get_number_items_itemsteps_dimensions()$number_dimensions, 
                      K = get_number_items_itemsteps_dimensions()$number_items, 
                      M = get_number_items_itemsteps_dimensions()$number_itemsteps, 
                      model = model)
    
    # little feedback
    if (!silent) 
      cat("\nItembank for",item_bank$model,"model.",item_bank$K,"items over",item_bank$Q,"dimension(s), with up to",item_bank$M+1,"categories per item.")
    
    item_bank
  }
  
  get_beta <- function() {
    # allow calculating beta from eta.
    if (model == "GPCM" && is.null(beta) && !is.null(eta))
      row_cumsum(eta)
    else
      beta
  }
  
  # define paramater matrices. Everything should ALWAYS be a matrix, except for index and m
  get_parameters_list <- function() { 
    pars = list(alpha = as.matrix(alpha),
                beta = as.matrix(get_beta()),
                guessing = if (is.null(guessing))
                  matrix(0, nrow(as.matrix(get_beta())), 1)
                else
                  as.matrix(guessing),
                index = 1:nrow(as.matrix(get_beta())), # except for bookkeeping things...
                m = number_non_missing_cells_per_row(as.matrix(get_beta()))
    ) 
    if (!is.null(eta)) 
      pars$eta <- as.matrix(eta)  
    
    pars
  }
  
  get_number_items_itemsteps_dimensions <- function() {
    pars = get_parameters_list()
    list(number_items = nrow(pars$beta), 
         number_itemsteps = ncol(pars$beta), 
         number_dimensions = ncol(pars$alpha))
  }
  
  validate <- function() {
    if (is.null(alpha))
      add_error("alpha", "is missing")
    if (model != "GPCM" && is.null(beta))
      add_error("beta", "is missing")
    if (model == "GPCM" && is.null(beta) && is.null(eta))
      add_error("beta_and_eta", "are both missing; define at least one of them")
    if (model == "GPCM" && !is.null(beta) && !is.null(eta) && !all.equal(row_cumsum(eta), as.matrix(beta)))
      add_error("beta_and_eta", "objects do not match, see details.")
  }
  
  invalid_result <- function() {
    list(pars = NA, 
         Q = NA, 
         K = NA, 
         M = NA, 
         model = NA,
         errors = errors())
  }
  
  validate_and_run()
}

