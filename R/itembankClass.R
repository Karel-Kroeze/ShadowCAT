#' Initialize an itembank object.
#' 
#' Produce an itembank object from a list of parameters.
#' 
#' This is a convenience function to create an itembank object as expected by further ShadowCAT functions.
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
#' @return ShadowCAT.items Itembank object, as used by \code{\link{initTest}} and \code{\link{initPerson}}.
#' @export
initItembank <- function(model = '3PLM', alpha = NULL, beta = NULL, guessing = NULL, eta = NULL, silent = FALSE){  
  result <- function() {  
    # define output, list
    # I will change Q (number of dimensions), K (number of items), M (number of item steps; number of categories 
    # minus 1), and m (in pars: the number of item steps per item) into appropriate names later, since this will affect other functions
    item_bank <- list(pars = get_parameters_list(), 
                      Q = get_number_items_itemsteps_dimensions()$number_dimensions, 
                      K = get_number_items_itemsteps_dimensions()$number_items, 
                      M = get_number_items_itemsteps_dimensions()$number_itemsteps, 
                      model = model)
    attr(item_bank, 'class') <- c("ShadowCAT.items")
    
    # little feedback
    if (!silent) 
      cat("\nItembank for",item_bank$model,"model.",item_bank$K,"items over",item_bank$Q,"dimension(s), with up to",item_bank$M+1,"categories per item.")
    
    invisible(item_bank)
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


#' Subset method for ShadowCAT.items.
#' 
#' Creates a valid subset of items to work with in ShadowCAT functions. Primarily used internally.
#' @param items  class ShadowCAT.items.
#' @param subset indeces.
#' @return items Subset of provided itembank, with class ShadowCAT.items, as well as any original classes. 
#' @export
subset.ShadowCAT.items <- function(items, subset) {
  out <- items
  out$subset <- subset
  out$K <- length(subset)
  
  for (par in names(items$pars)){
    if (is.null(dim(items$pars[[par]]))) { # vector
      out$pars[[par]] <- items$pars[[par]][subset,drop = FALSE]
    } else { # matrix
      out$pars[[par]] <- items$pars[[par]][subset,,drop = FALSE]     
    }
  }
  
  return(out)
}

#' Create a test Itembank
#' 
#' Quick and simple itembanks for testing purposes.
#' @param model
#' @param K
#' @param Q
#' @param M
#' @return ShadowCAT.itembank
#' @export
createTestBank <- function(model, K = 50, Q = 1, M = 4, between = FALSE){
  # 3PLM is dichotomous by definition
  if (model == "3PLM") M <- 1 
  
  # make sure the number of items is divisible by the number of dimensions
  if (between) K <- ceiling(K/Q) * Q 
  
  # set up alpha, very rough uniform from .3 to 1.5
  alpha <- matrix(runif(K * Q, .3, 1.5), K, Q)
  
  # if between, force items to load on one dimension each.
  if (between){
    set = K / Q
    for (i in 1:Q){
      alpha[((i-1)*set+1):(i*set), (1:Q)[-i]] <- 0    
    }
  }
  
  # spread polytomous items cats -2 to +2.
  spread <- seq(-2,2,length.out=M)
  
  # base loading for items
  beta <- matrix(rnorm(K), K, 1)
  
  # apply spread for polytomous, betas are strictly monotously increasing because the spread is.
  # apply transposes the matrix...
  if (M > 1) beta <- t(apply(beta, 1, function(x) x + spread))
  
  # reparameterize GPCM
  if (model == "GPCM") {
    # make betas
    eta <- beta
    for (i in 1:M) {
      # rolling sum, apply over all items
      beta[,i] <- apply(eta[,1:i, drop=FALSE], 1, sum)
    }
  }
  
  # create Itembank object
  items <- initItembank(model, alpha, beta, silent = TRUE)  
  
  # return
  return(invisible(items))
}