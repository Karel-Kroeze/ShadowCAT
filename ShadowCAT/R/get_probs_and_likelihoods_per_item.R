#' Probabilities and likelihoods
#' 
#' Per item, get probabilities of scoring in each answer category given theta, likelihood for
#' the actually chosen answer category, and derivatives.
#' 
#' @param theta Vector with true or estimated theta.
#' @param answers If \code{with_likelihoods} is set to \code{TRUE}, vector with answers to all items included in \code{alpha}, 
#' \code{beta}, and \code{guessing}; else, numeric(0).
#' @param with_likelihoods If \code{FALSE}, only the probability matrix for all answer categories for all items is returned. If \code{TRUE}, the likelihoods 
#' for the actually chosen answer categories and derivatives are also returned.
#' @return If \code{with_likelihoods}, list with probability matrix, likelihoods and derivatives. Else, only the probability matrix.
#' @inheritParams shadowcat
#' @importFrom mvtnorm dmvnorm
#' @importFrom Rcpp evalCpp
#' @useDynLib ShadowCAT
get_probs_and_likelihoods_per_item <- function(theta, model, alpha, beta, guessing, answers = numeric(0), with_likelihoods) {
  probs_and_likelihoods <- switch(model,
                                  "3PLM" = PROB_3PLM(theta, alpha, beta, guessing, answers, with_likelihoods),
                                  "GRM" = PROB_GRM(theta, alpha, beta, answers, with_likelihoods),
                                  "SM" = PROB_SM(theta, alpha, beta, answers, with_likelihoods),
                                  "GPCM" = PROB_GPCM(theta, alpha, beta, answers, with_likelihoods))
  # likelihoods can never truly be zero, let alone negative
  # As far as I have seen, only the GPCM model may return very small negative likelihoods, due to
  # the code line remainder -= P(k, i+1), which sometimes is something like 1 - 2.092453e-17 - 1.243433e-08 - 1 = -2.220446e-16.
  # Other models may return likelihoods equal to 0.000000e+00. It seems safe to set these to very small positive likelihoods.
  probs_and_likelihoods$P[which(probs_and_likelihoods$P <= 0)] <- 1e-10
  if (with_likelihoods) {
    probs_and_likelihoods$l[which(probs_and_likelihoods$l <= 0)] <- 1e-10
    probs_and_likelihoods
  }
  else {
    probs_and_likelihoods$P
  }
}

