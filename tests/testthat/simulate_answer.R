#' Simulates responses on items indicated by indeces, given true theta
#' 
#' @examples 
#' items <- simulate_testbank("GPCM")
#' 
#' # simulates responses on items indicated by indeces, given true theta
#' simulate_answer(.3, "GPCM", 1, "maximum_aposteriori", items$pars$alpha, items$pars$beta, items$pars$guessing, items$M, 3)
#' 
#' @param theta vector with true theta
#' @param model String, one of '3PLM', 'GPCM', 'SM' or 'GRM', for the three-parameter logistic, generalized partial credit, sequential or graded response model respectively.
#' @param number_dimensions number of dimensions
#' @param estimator type of estimator to be used, one of "maximum_aposteriori", "maximum_likelihood", or "expected_aposteriori"
#' @param alpha matrix of alpha parameters
#' @param beta matrix of beta parameters
#' @param guessing vector of guessing parameters
#' @param number_itemsteps number of itemsteps
#' @param indeces vector of indeces for which answers should be simulated
#' @return vector with responses
simulate_answer <- function(theta, model, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps, indeces) {
  result <- function() {
    guessing <- get_guessing()
    # probabilities, generated with true theta.
    probabilities <- get_probabilities(theta, model, indeces, alpha, beta, guessing)
    cumulative_probabilities <- row_cumsum(probabilities) 
    random_numbers <- runif(length(indeces))
    
    # answer is the number of categories that have a cumulative probability smaller than random_numbers
    apply(random_numbers > cumulative_probabilities, 1, sum, na.rm=TRUE)
  }

  get_guessing <- function() {
    if (is.null(guessing))
      matrix(0, nrow(as.matrix(beta)), 1)
    else
      guessing
  }
  
  result()
}