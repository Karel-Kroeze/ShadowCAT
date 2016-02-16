#' Fisher Information
#' 
#' Fisher Information (expected information) per item.
#' 
#' Calculates the Fisher Information (expected information) of the given test, and returns a three dimensional array
#' of information matrices, where dimensions one and two run along the Q dimensions of the model, and three runs along items.
#' 
#' Fisher Information is given as;
#' \deqn{\mathcal{I}(\theta) = - \operatorname{E} \left[\left. \frac{\partial^2}{\partial\theta^2} \log f(X;\theta)\right|\theta \right]}{minus expectation of second derivative of the Log-Likelihood of f(theta)}
#' 
#' And is calculated as the weighted sum of second derivatives for all response categories. Information for multiple items is simply the sum of 
#' the individual information matrices.
#' 
#' Note: get_fisher_information always returns the 'raw' information, information given by prior distributions is added by the calling functions, 
#' and get_fisher_information(..) is normally called internally.
#' 
#' @param estimate vector containing theta estimate
#' @param model string, one of '3PLM', 'GPCM', 'SM' or 'GRM', for the three-parameter logistic, generalized partial credit, sequential or graded response model respectively.
#' @param number_dimensions number of dimensions
#' @param alpha matrix containing the alpha parameters
#' @param beta matrix containing the beta parameters
#' @param guessing matrix containing the quessing parameters
#' @param number_itemsteps_per_item vector containing the number of non missing cells per row of the beta matrix
#' @return array with an information matrix for each item (QxQxK).
#' @export
get_fisher_information <- function(estimate, model, number_dimensions, alpha, beta, guessing, number_itemsteps_per_item) {
  # Fisher Information
  # minus the expectation of the second derivative of the log-likelihood
  # Expectation -> sum of derivatives for each category, 'weighted' by their probability.
  number_items <- nrow(alpha)
  probabilities <- get_probs_and_likelihoods_per_item(estimate, model, alpha, beta, guessing, with_likelihoods = FALSE)$P
    
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
    # Sequential Model (Tutz, xxxx)
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