#' Probabilities and likelihoods
#' 
#' Per item, get probabilities of scoring in each answer category given theta, likelihood for
#' the actually chosen answer category, and derivatives.
#' 
#' @references 
#' \itemize{
#' \item Muraki, E. (1992). A Generized Partial Credit Model: Application of an EM Algorithm. 
#'  Applied Psychological Measurement, 16(2), 159 - 176. Doi:10.1177/014662169201600206.
#'  \item Samejima, F. (1970). Estimation of latent trait ability using a response pattern of graded 
#'  scores. Psychometrika, 35(1), 139 - 139. Doi: 10.1007/BF02290599.
#' \item Segall, D. O. (2000). Principles of multidimensional adaptive testing. In W. J. van der 
#'  Linden & en C. A. W. Glas (Eds.), Computerized adaptive testing: Theory and 
#'  practice (pp. 53 - 74). Dordrecht: Kluwer Academic Publishers.
#'  \item Tutz, G. (1986). Bradley-Terry-Luce model with an ordered response. Journal of 
#'  Mathematical Psychology, 30(1), 306 - 316. doi: 10.1016/0022-2496(86)90034-9.
#' }
#' 
#' @param theta Vector with true or estimated theta.
#' @param number_dimensions Number of dimensions in theta.
#' @param number_items Number of items in alpha, beta, guessing.
#' @param number_itemsteps_per_item Vector containing the number of non missing cells per row of the beta matrix.
#' @param answers If \code{with_likelihoods} is set to \code{TRUE}, vector with answers to all items included in \code{alpha}, 
#' \code{beta}, and \code{guessing}; else, \code{NULL}.
#' @param with_likelihoods If \code{FALSE}, only the probability matrix for all answer categories for all items is returned. If \code{TRUE}, the likelihoods 
#' for the actually chosen answer categories are also returned.
#' @param with_derivatives If \code{FALSE}, only probability matrix and likelihoods are returned. If \code{TRUE}, the first and second derivatives are also returned.
#' Ignored if \code{with_likelihoods} is \code{FALSE}.
#' @return If \code{with_likelihoods} is \code{FALSE}, only the probability matrix is returned.
#' If \code{with_likelihoods} is \code{TRUE} and \code{with_derivatives} is \code{FALSE}, \code{list} with probability matrix and likelihoods.
#' If \code{with_likelihoods} is \code{TRUE} and \code{with_derivatives} is \code{TRUE}, \code{list} with probability matrix, likelihoods and derivatives.
#' Else, only the probability matrix.
#' @inheritParams shadowcat
#' @importFrom mvtnorm dmvnorm
#' @importFrom Rcpp evalCpp
#' @useDynLib ShadowCAT
get_probs_and_likelihoods_per_item <- function(theta, model, alpha, beta, guessing, number_dimensions, number_items, number_itemsteps_per_item, answers = NULL, with_likelihoods = FALSE, with_derivatives = FALSE) {
  max_number_itemsteps <- max(number_itemsteps_per_item)
  
  result <- function() {
    probs_and_likelihoods <- switch(model,
                                    "3PLM" = result_3plm(),
                                    "GRM" = result_grm(),
                                    "SM" = result_sm(),
                                    "GPCM" = result_gpcm())
   }
   
   # likelihoods can never truly be zero, let alone negative
   # As far as I have seen, only the GPCM model may return very small negative likelihoods, due to
   # the code line 1 - sum(probs), which sometimes is something like 1 - 2.092453e-17 - 1.243433e-08 - 1 = -2.220446e-16.
   # Other models may return likelihoods equal to 0.000000e+00. It seems safe to set these to very small positive likelihoods.
  result_3plm <- function() {
    probabilities <- get_probabilities_3plm()
    probabilities[which(probabilities <= 0)] <- 1e-10
    if (!with_likelihoods)
      probabilities
    else if (!with_derivatives)
      list(probabilities = probabilities,
           likelihoods = get_likelihoods(probabilities))
    else
      list(probabilities = probabilities,
           likelihoods = get_likelihoods(probabilities),
           derivatives1 = get_derivatives1_3plm(probabilities),
           derivatives2 = get_derivatives2_3plm(probabilities))
  }
  
  result_grm <- function() {
    psi <- get_psi()
    probabilities <- get_probabilities_grm(psi)
    probabilities[which(probabilities <= 0)] <- 1e-10
    if (!with_likelihoods)
      probabilities
    else if (!with_derivatives)
      list(probabilities = probabilities,
           likelihoods = get_likelihoods(probabilities))
    else
    list(probabilities = probabilities,
         likelihoods = get_likelihoods(probabilities),
         derivatives1 = get_derivatives1_grm(psi),
         derivatives2 = get_derivatives2_grm(psi))
  }
  
  result_sm <- function() {
    psi <- get_psi()
    probabilities <- get_probabilities_sm(psi)
    probabilities[which(probabilities <= 0)] <- 1e-10
    if (!with_likelihoods)
      probabilities
    else if (!with_derivatives)
      list(probabilities = probabilities,
           likelihoods = get_likelihoods(probabilities))
    else
    list(probabilities = probabilities,
         likelihoods = get_likelihoods(probabilities),
         derivatives1 = get_derivatives1_sm(psi),
         derivatives2 = get_derivatives2_sm(psi))
  }
  
  result_gpcm <- function() {
    probabilities <- get_probabilities_gpcm()
    probabilities[which(probabilities <= 0)] <- 1e-10
    if (!with_likelihoods) {
      probabilities
    }
    else if (!with_derivatives) {
      list(probabilities = probabilities,
           likelihoods = get_likelihoods(probabilities))
    }
    else {
      one_to_number_itemsteps_times_probs <- get_one_to_number_itemsteps_times_probs(probabilities_gpcm = probabilities, max_number_itemsteps = max_number_itemsteps)
      colsums_one_to_number_itemsteps_times_probs <- colSums(one_to_number_itemsteps_times_probs, na.rm = TRUE)
      list(probabilities = probabilities,
           likelihoods = get_likelihoods(probabilities),
           derivatives1 = get_derivatives1_gpcm(colsums_one_to_number_itemsteps_times_probs = colsums_one_to_number_itemsteps_times_probs),
           derivatives2 = get_derivatives2_gpcm(one_to_number_itemsteps_times_probs = one_to_number_itemsteps_times_probs, 
                                                colsums_one_to_number_itemsteps_times_probs = colsums_one_to_number_itemsteps_times_probs,
                                                max_number_itemsteps = max_number_itemsteps))
    }
  }
  
  logistic_function <- function(x) {
    exp(x) / (1 + exp(x))
  }
  
  get_likelihoods <- function(probabilities) {
    likelihoods <- sapply(1:number_items, function(item) { probabilities[item, answers[item] + 1] })
    likelihoods[which(likelihoods <= 0)] <- 1e-10
    likelihoods
  }
  
  psi_single_item_itemstep <- function(itemstep, item) { 
    logistic_function(alpha[item, ] %*% theta - beta[item, itemstep]) 
  }
  
  psi_single_item <- function(item) {
    psi_item <- sapply(1:number_itemsteps_per_item[item], psi_single_item_itemstep, item = item)
    number_zeros <- max_number_itemsteps - number_itemsteps_per_item[item] + 1
    c(1, psi_item, 0, rep(NA, max_number_itemsteps - length(psi_item)))
  }
  
  get_psi <- function() {
    t(sapply(1:number_items, psi_single_item))
  } 
  
  # Three Paramater Logistic Model (MultiDim) (Segall, 1996)
  get_probabilities_single_item_3plm <- function(item) {
    p_1 <- guessing[item] + (1 - guessing[item]) / (1 + exp(-sum(alpha[item, ] * (theta - beta[item]))))
    p_0 <- 1 - p_1
    c(p_0, p_1)
  }
  
  get_probabilities_3plm <- function() {
    t(sapply(1:number_items, get_probabilities_single_item_3plm))
  }
  
  get_derivatives1_3plm <- function(probabilities_3plm) {
    ((probabilities_3plm[, 2] - guessing) * (answers - probabilities_3plm[, 2])) / ((1 - guessing) * probabilities_3plm[, 2])
  }
  
  get_derivatives2_3plm <- function(probabilities_3plm) {
    (probabilities_3plm[, 1] * (probabilities_3plm[, 2] - guessing) * (guessing * answers - probabilities_3plm[, 2]^2)) / (probabilities_3plm[, 2]^2 * (1 - guessing)^2)
  }
  
  # Graded Response Model (Samejima, 1969)
  get_probabilities_grm <- function(psi) {
    t(apply(psi, 1, function(psi_one_item) { -diff(psi_one_item) }))
  }
  
  get_derivatives1_grm <- function(psi) {
    sapply(1:number_items,
           function(item) { 1 - psi[item, answers[item] + 1] - psi[item, answers[item] + 2] })
  }
  
  get_derivatives2_grm <- function(psi) {
    sapply(1:number_items,
           function(item) { -(psi[item, answers[item] + 1] * (1 - psi[item, answers[item] + 1]) + 
                              psi[item, answers[item] + 2] * (1 - psi[item, answers[item] + 2])) }) 
  }
  
  # Sequential Model (Tutz, 1990)
  get_probabilities_single_item_sm <- function(item, psi, psi_cumprod) { 
    probs <- (1 - psi[item, 2:(number_itemsteps_per_item[item] + 2)]) * psi_cumprod[item, 1:(number_itemsteps_per_item[item] + 1)]
    c(probs, rep(NA, max_number_itemsteps - length(probs) + 1))
  }

  get_probabilities_sm <- function(psi) {
    psi_cumprod <- t(apply(psi, 1, cumprod))
    t(sapply(1:number_items, get_probabilities_single_item_sm, psi = psi, psi_cumprod = psi_cumprod))
  }
  
  get_derivatives1_sm <- function(psi) {
    sapply(1:number_items,
           function(item) { -psi[item, answers[item] + 2] + (answers[item] + 1) - sum(psi[item, 1:(answers[item] + 1)]) })
  }
  
  get_derivatives2_sm <- function(psi) {
    sapply(1:number_items,
           function(item) { -psi[item, answers[item] + 2] * (1 - psi[item, answers[item] + 2]) - sum(psi[item, 1:(answers[item] + 1)] * (1 - psi[item, 1:(answers[item] + 1)])) })
  }
  
  
  # Generalized Partial Credit Model (Muraki, 1992)
  get_probabilities_single_item_gpcm <- function(item) {
    aux <- exp(1:number_itemsteps_per_item[item] * (alpha[item, ] %*% theta) - beta[item, 1:number_itemsteps_per_item[item]])
    probs <- aux / (1 + sum(aux))
    c(1 - sum(probs), probs, rep(NA, max_number_itemsteps - length(probs)))
  }
  
  get_probabilities_gpcm <- function() {
    t(sapply(1:number_items, get_probabilities_single_item_gpcm))
  }
  
  get_one_to_number_itemsteps_times_probs <- function(probabilities_gpcm, max_number_itemsteps) {
    t(probabilities_gpcm[, 2:ncol(probabilities_gpcm)]) * 1:max_number_itemsteps
  }
  
  get_itemsteps_minus_colsums_one_to_number_itemsteps_times_probs <- function(colsums_one_to_number_itemsteps_times_probs) {
    sapply(1:max_number_itemsteps, 
           function(itemstep, colsums_one_to_number_itemsteps_times_probs) { itemstep - colsums_one_to_number_itemsteps_times_probs },
           colsums_one_to_number_itemsteps_times_probs = colsums_one_to_number_itemsteps_times_probs)
  }
  
  get_derivatives1_gpcm <- function(colsums_one_to_number_itemsteps_times_probs) {
    answers - colsums_one_to_number_itemsteps_times_probs
  } 
  
  get_derivatives2_gpcm <- function(one_to_number_itemsteps_times_probs, colsums_one_to_number_itemsteps_times_probs, max_number_itemsteps) {
    itemsteps_minus_colsums_one_to_number_itemsteps_times_probs <- t(get_itemsteps_minus_colsums_one_to_number_itemsteps_times_probs(colsums_one_to_number_itemsteps_times_probs))
    -colSums(one_to_number_itemsteps_times_probs * itemsteps_minus_colsums_one_to_number_itemsteps_times_probs, na.rm = TRUE)
  }
  
  result()
}

