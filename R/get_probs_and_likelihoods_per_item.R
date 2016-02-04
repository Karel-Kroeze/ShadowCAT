#' get list containing probabilities of scoring in each answer category for each item, and
#' likelihood with derivative per item
#' 
#' @param theta vector with true or estimated theta
#' @param model string, one of '3PLM', 'GPCM', 'SM' or 'GRM', for the three-parameter logistic, generalized partial credit, sequential or graded response model respectively.
#' @param alpha matrix containing the alpha parameters (if with_likelihoods, only of administered items)
#' @param beta matrix containing the beta parameters of administered items (if with_likelihoods, only of administered items)
#' @param guessing matrix containing the quessing parameters of administered items (if with_likelihoods, only of administered items)
#' @param responses if with_likelihoods = TRUE, vector with person responses to the administered items; else, numeric(0) 
#' @param with_likelihoods if FALSE, only the probability matrix is returned, if TRUE, the likelihoods and derivatives are also returned
#' @return list with P = matrix with for each included item (rows) the probability of scoring in each answer category (columns), given theta,
#' and if with_likelihoods, l = vector of likelihoods per item, d = vector of first derivatives, D = vector of second derivatives
#' @importFrom mvtnorm dmvnorm
#' @importFrom Rcpp evalCpp
#' @useDynLib ShadowCAT
#' @export
get_probs_and_likelihoods_per_item <- function(theta, model, alpha, beta, guessing, responses = numeric(0), with_likelihoods) {
  probs_and_likelihoods <- switch(model,
                                  "3PLM" = PROB_3PLM(theta, alpha, beta, guessing, responses, with_likelihoods),
                                  "GRM" = PROB_GRM(theta, alpha, beta, responses, with_likelihoods),
                                  "SM" = PROB_SM(theta, alpha, beta, responses, with_likelihoods),
                                  "GPCM" = PROB_GPCM(theta, alpha, beta, responses, with_likelihoods))
  # likelihoods can never truly be zero, let alone negative
  probs_and_likelihoods$P[which(probs_and_likelihoods$P <= 0)] <- 1e-10
  if (with_likelihoods)
    probs_and_likelihoods$l[which(probs_and_likelihoods$l <= 0)] <- 1e-10
  probs_and_likelihoods
}

