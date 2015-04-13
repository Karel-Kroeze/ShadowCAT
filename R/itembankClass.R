
####### Itembank, and everything related to it.
#' initItembank
#' 
#' Produce an itembank object for a list of parameters
#' 
#' @param model
#' @param alpha
#' @param beta
#' @param guessing
#' @param eta
#' @param m
#' @return ShadowCAT.items
#' @export
initItembank <- function(model = '3PLM', alpha = NULL, beta = NULL, guessing = NULL, eta = NULL, silent = FALSE){
  # TODO: Validate input.
  
  # find number of items, categories
  if (is.null(dim(beta))) {
    K <- length(beta) # length of vector
    M <- 1 # vector, so couldn't be polytomous
  } else {
    K <- nrow(beta) # row per item
    M <- ncol(beta) # column per option
  }
  
  
  # find number of dimensions
  if (is.null(dim(alpha))) {
    Q <- 1 # this is a vector
  } else {
    Q <- ncol(alpha)
  }
  
  pars <- list()
  # define paramater matrices. Everything should ALWAYS be a matrix.
  pars$alpha <- matrix(alpha, K, Q)
  pars$beta <- matrix(beta, K, M)
  if ( ! is.null(guessing)) pars$guessing <- matrix(guessing, K, 1)
  else pars$guessing <- matrix(0, K, 1)
  if ( ! is.null(eta)) pars$eta <- matrix(eta, K, Q)
  
  pars$index <- 1:K # except for bookkeeping things...
  pars$m <- apply(pars$beta, 1, function(x) sum(! is.na(x))) # count number of non-NA cells per row - this is m, the number of cats per row.
  
  
  # define output, list
  out <- list(pars = pars, Q = Q, K = K, M = M, model = model)
  attr(out, 'class') <- c("ShadowCAT.items")
  
  
  # little feedback
  if ( ! silent) cat("\nItembank for",model,"model.",K,"items over",Q,"dimension(s), with up to",M,"categories per item.")
  
  
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