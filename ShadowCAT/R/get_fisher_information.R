#' Fisher Information
#' 
#' Fisher Information (expected information) per item.
#' 
#' @details
#' Fisher Information is given as:
#' \deqn{\mathcal{I}(\theta) = - \operatorname{E} \left[\left. \frac{\partial^2}{\partial\theta^2} \log f(X;\theta)\right|\theta \right]}{minus expectation of second derivative of the Log-Likelihood of f(theta)}
#' and is calculated as the weighted sum of second derivatives for all response categories. Information for multiple items is simply the sum of 
#' the individual information matrices.
#' 
#' Note: get_fisher_information always returns the 'raw' information; information given by prior distributions is added by the calling functions.
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
#' @param number_dimensions Number of dimensions of theta.
#' @param number_itemsteps_per_item Vector containing the number of non missing cells per row of the beta matrix.
#' @inheritParams shadowcat
#' @return Three dimensional array of information matrices, where dimensions one and two run along the Q dimensions of the model, and three runs along items.
get_fisher_information <- function(estimate, model, number_dimensions, alpha, beta, guessing, number_itemsteps_per_item) {
  # Fisher Information
  # minus the expectation of the second derivative of the log-likelihood
  # Expectation -> sum of derivatives for each category, 'weighted' by their probability.
  number_items <- nrow(alpha)
  probabilities <- get_probs_and_likelihoods_per_item(estimate, model, alpha, beta, guessing, with_likelihoods = FALSE)
    
  result <- function() {
    second_derivatives <- get_second_derivatives(model)
    
    # just dump everything back, how to deal with information is up to caller function.  
    lapply_return_array(1:number_items, 
                        dim = c(number_dimensions, number_dimensions, number_items), 
                        FUN = function(item) { 
                          (alpha[item,] %*% t(alpha[item,])) * second_derivatives[item] 
                        } )
  }
  
  get_second_derivatives <- function(model) {
    switch(model,
           "3PLM" = get_second_derivatives_3plm(),
           "GRM" = get_second_derivatives_grm(),
           "SM" = get_second_derivatives_sm(),
           "GPCM" = get_second_derivatives_gpcm())
  }
  
  logistic_function <- function(x) {
    exp(x) / (1 + exp(x))
  }
  
  get_second_derivatives_3plm <- function() {
    # exact form Segall 1997, CAT book, p. 72
    # p[,1] = q, p[, 2] = p
    (probabilities[,1] / probabilities[,2]) * ((probabilities[,2] - guessing)/(1 - guessing))^2
  }
    
  get_second_derivatives_grm <- function() {
    # Graded Response Model (Glas & Dagohoy, 2007)
    inner_product_alpha_theta <- as.vector(alpha %*% drop(estimate))
    sapply(1:number_items, 
           function(item) {
             psi <- c(1, logistic_function(inner_product_alpha_theta[item] - beta[item, 1:number_itemsteps_per_item[item]]), 0)
             sum(probabilities[item, 1:(number_itemsteps_per_item[item] + 1)] * (psi[1:(number_itemsteps_per_item[item] + 1)] * (1 - psi[1:(number_itemsteps_per_item[item] + 1)]) + psi[2:(number_itemsteps_per_item[item] + 2)] * (1 - psi[2:(number_itemsteps_per_item[item] + 2)])))
            })
  }

  get_second_derivatives_sm <- function() {
    # Sequential Model (Tutz, 1986)
    inner_product_alpha_theta <- as.vector(alpha %*% drop(estimate))
    sapply(1:number_items,
           function(item) {
             psi <- c(1, logistic_function(inner_product_alpha_theta[item] - beta[item, 1:number_itemsteps_per_item[item]]), 0) # basically: 1, logistic_function(inner_product_alpha_theta - b), 0.
             psi_sums <- sapply(1:(number_itemsteps_per_item[item] + 1),
                                function(item_step) {
                                  sum(psi[2:(item_step + 1)] * (1 - psi[2:(item_step + 1)]))
                                })
             sum(probabilities[item, 1:(number_itemsteps_per_item[item] + 1)] * psi_sums[1:(number_itemsteps_per_item[item] + 1)])
             
           })
  }
  
  get_second_derivatives_gpcm <- function() {
    # Generalized Partial Credit Model (Muraki, 1992)
    sapply(1:number_items,
           function(item) {
             mi <- 1:number_itemsteps_per_item[item]
             pi <- probabilities[item, mi + 1] # remove j = 0, index is now also correct.
             mp <- sum(mi*pi)  
             sum((mi * pi) * (mi - mp))
           })
  }
  
  result()
}
