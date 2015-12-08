#' Simulates responses on items indicated by indeces, given true theta
#' 
#' @examples 
#' items <- simulate_testbank("GPCM")
#' 
#' # simulates responses on items indicated by indeces, given true theta
#' simulate_answer(.3, "GPCM", 1, "MAP", items$pars$alpha, items$pars$beta, items$pars$guessing, items$M, 3)
#' 
#' @param theta true theta
#' @param model String, one of '3PLM', 'GPCM', 'SM' or 'GRM', for the three-parameter logistic, generalized partial credit, sequential or graded response model respectively.
#' @param number_dimensions number of dimensions
#' @param estimator type of estimator to be used, one of "MAP" (Maximum a posteriori estimation) or "ML" (maximum likelihood); 
#' "EAP" (Expected A Posteriori Estimation) is currently not working due to problems with the MultiGHQuad package
#' @param alpha matrix of alpha parameters
#' @param beta matrix of beta parameters
#' @param guessing vector of guessing parameters
#' @param number_itemsteps number of itemsteps
#' @param indeces answer questions with the given indeces
#' @return vector responses, or updated person object if indeces is set.
simulate_answer <- function(theta, model, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps, indeces) {
  if (is.null(guessing))
    guessing <- matrix(0, nrow(as.matrix(beta)), 1)

  # probabilities, generated with TRUE theta.
  Pij <- probabilities_and_likelihoods(theta, responses = NULL, model, indeces, number_dimensions, estimator, alpha, beta, guessing, output = "probs")
  
  # cumulative probabilities
  cp <- Pij 
  for (i in 1:(number_itemsteps + 1)) cp[,i] <- apply(matrix(Pij[,1:i],ncol=i),1,sum)
   
  # rand ~ unif(0,1)
  rand <- runif(length(indeces))
  
  # answer is the number of categories that have a cumulative probability smaller than rand
  apply(rand > cp, 1, sum, na.rm=TRUE)
}