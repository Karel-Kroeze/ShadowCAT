#' Get matrix with for each included item the probability of scoring in each answer category, given theta
#' 
#' @param theta vector with true or estimated theta
#' @param model string, one of '3PLM', 'GPCM', 'SM' or 'GRM', for the three-parameter logistic, generalized partial credit, sequential or graded response model respectively.
#' @param items_to_include vector with indeces of items to include
#' @param alpha matrix containing the alpha parameters (for complete test bank)
#' @param beta matrix containing the beta parameters (for complete test bank)
#' @param guessing matrix containing the quessing parameters (for complete test bank)
#' @return matrix with for each included item (rows) the probability of scoring in each answer category (columns), given theta
#' @export
get_probabilities <- function(theta, model, items_to_include, alpha, beta, guessing) {
  # TODO: Check input.
  # TODO priors: mean? 
  # priors: Alleen variabele deel van multivariaat normaal verdeling (exp).
  alpha <- get_subset(alpha, items_to_include)
  beta <- get_subset(beta, items_to_include)
  guessing <- get_subset(guessing, items_to_include)
  get_probs_and_likelihoods_per_item(theta, model, alpha, beta, guessing, numeric(0), FALSE)$P
}
