#' Create a test Itembank
#' 
#' Quick and simple itembanks for testing purposes.
#' @param model String, one of '3PLM', 'GPCM', 'SM' or 'GRM', for the three-parameter logistic, generalized partial credit, sequential or graded response model respectively.
#' @param K number of items
#' @param Q number of dimensions
#' @param M number of item steps (number of categories minus 1)
#' @param between is TRUE, force items to load on one dimension each
#' @param run_initItembank if FALSE, the simulated testbank is returned; if TRUE, initItembank() is applied on
#' the simulated testbank and the result of this is returned 
#' @return ShadowCAT.itembank
#' @export
createTestBank <- function(model, K = 50, Q = 1, M = 4, between = FALSE, run_initItembank = TRUE){
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
    beta <- row_cumsum(eta)
  } 
  
  if (run_initItembank)
    invisible(initItembank(model, alpha, beta, silent = TRUE))
  else
    list(alpha = alpha,
         beta = beta)
}

