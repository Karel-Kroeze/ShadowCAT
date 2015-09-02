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
#' @param beta Matrix of beta parameters, one column per response category, one row per item. Note that ShadowCAT expects response categories to be sequential,
#' and without gaps. That is, the weight parameter in the GPCM model is assumed to be sequential, and equal to the position of the 'location' of the beta parameter in the Beta matrix.
#' The matrix will have a number of columns equal to the largest number of response categories, items with fewer response categories should be 
#' right-padded with \code{NA}. \code{NA} values between response categories are not allowed, and will lead to errors.
#' More flexibility in Beta parameters might be added in future versions.
#' @param guessing vector of guessing parameters per item. Optionally used in 3PLM model, ignored for all others.
#' @param eta Matrix of location parameters, optionally used in GPCM model, ignored for all others.
#' @return ShadowCAT.items Itembank object, as used by \code{\link{initTest}} and \code{\link{initPerson}}.
#' @export
initItembank <- function(model = '3PLM', alpha = NULL, beta = NULL, guessing = NULL, eta = NULL, silent = FALSE){
  # TODO: Validate input.
  
  # check eta/beta are not mixed up
  if (model == "GPCM" && !is.null(beta) && !is.null(eta)) {
    criterion_beta <- transpose_if_ncol_and_nrow_larger_1(matrix_apply(eta, 1, cumsum))
    if (!all.equal(criterion_beta, as.matrix(beta))) stop("Beta and Eta parameters do not match, see details.")
  }
  
  # allow calculating beta from eta.
  if (model == "GPCM" && is.null(beta) && !is.null(eta))
    beta <- transpose_if_ncol_and_nrow_larger_1(matrix_apply(eta, 1, cumsum))
    
  # define paramater matrices. Everything should ALWAYS be a matrix.
  pars <- list()
  pars$alpha <- as.matrix(alpha)
  pars$beta <- as.matrix(beta)
  pars$guessing <- if (is.null(guessing))
                     matrix(0, nrow(pars$beta), 1)
                   else
                     as.matrix(guessing) 
  if (!is.null(eta)) 
    pars$eta <- as.matrix(eta)
  
  # find number of items, item steps (number of categories minus 1), dimensions
  K <- nrow(pars$beta) # row per item
  M <- ncol(pars$beta) # column per item step
  Q <- ncol(pars$alpha) # column per dimension
  
  pars$index <- 1:K # except for bookkeeping things...
  pars$m <- apply(pars$beta, 1, function(x) sum(!is.na(x))) # count number of non-NA cells per row - this is m, the number of cats per row.
   
  # define output, list
  out <- list(pars = pars, Q = Q, K = K, M = M, model = model)
  attr(out, 'class') <- c("ShadowCAT.items")
   
  # little feedback
  if ( ! silent) cat("\nItembank for",model,"model.",K,"items over",Q,"dimension(s), with up to",M+1,"categories per item.")
  
  # return itembank
  return(invisible(out))
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