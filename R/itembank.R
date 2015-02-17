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
  if ( ! is.null(eta)) pars$eta <- matrix(eta, K, Q)
  
  pars$index <- 1:K # except for bookkeeping things...
  pars$m <- apply(pars$beta, 1, function(x) sum(! is.na(x))) # count number of non-NA cells per row - this is m, the number of cats per row.
  
  
  # define output, list
  out <- list(pars = pars, Q = Q, K = K, M = M)
  attr(out, 'class') <- c("ShadowCAT.items", paste0("ShadowCAT.items.",model))
  
  
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
      out$pars[[par]] <- items$pars[[par]][subset]
    } else { # matrix
      out$pars[[par]] <- items$pars[[par]][subset,,drop = FALSE]     
    }
  }
  
  return(out)
}