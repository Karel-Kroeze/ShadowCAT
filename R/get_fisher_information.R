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
#' Note: FI always returns the 'raw' information, information given by prior distributions is added by the calling functions, 
#' and FI(..) is normally called internally.
#' 
#' @param test Test object, see \code{\link{initTest}}.
#' @param person Person object, see \code{\link{initPerson}}.
#' @return array with an information matrix for each item (QxQxK).
#' @export
FI <- function(test, person) {
  # Fisher Information
  # minus the expectation of the second derivative of the log-likelihood
  # Expectation -> sum of derivatives for each category, 'weighted' by their probability.
  probabilities <- prob(test, person)$P
  
  result <- function() {
    D = get_D(test$items$model)
    
    # just dump everything back, how to deal with information is up to caller function.  
    lapply_return_array(1:test$items$K, 
                        dim = c(test$items$Q, test$items$Q, test$items$K), 
                        FUN = function(item) { 
                          (test$items$pars$alpha[item,] %*% t(test$items$pars$alpha[item,])) * D[item] 
                        } )
  }
  
  get_D <- function(model) {
    switch(model,
           "3PLM" = get_D_3PLM(),
           "GRM" = get_D_GRM(),
           "SM" = get_D_SM(),
           "GPCM" = get_D_GPCM())
  }
  
  logistic_function <- function(x) {
    exp(x)/(1+exp(x))
  }
  
  get_D_3PLM <- function() {
    # exact form Segall 1997, CAT book, p. 72
    # p[,1] = q, p[, 2] = p
    (probabilities[,1] / probabilities[,2]) * ((probabilities[,2] - test$items$pars$guessing)/(1 - test$items$pars$guessing))^2
  }
  
  get_D_GRM <- function() {
    # Graded Response Model (Glas & Dagohoy, 2007)
    inner_product_alpha_theta <- apply(test$items$pars$alpha * drop(person$estimate), 1, sum)
    D <- numeric(test$items$K)
    for (item in 1:test$items$K) {
      for (item_step in 1:(test$items$pars$m[item] + 1)) {
        Psi <- c(1, logistic_function(inner_product_alpha_theta[item] - test$items$pars$beta[item, 1:test$items$pars$m[item]]), 0)
        D[item] <- D[item] + probabilities[item, item_step] * (Psi[item_step] * (1 - Psi[item_step]) + Psi[item_step + 1] * (1 - Psi[item_step + 1]))
      }
    }
    D
  }

  get_D_SM <- function() {
    # Sequential Model (Tutz, xxxx)
    # TODO: triple check this.
    inner_product_alpha_theta <- apply(test$items$pars$alpha * drop(person$estimate), 1, sum)
    D <- numeric(test$items$K)
    for (item in 1:test$items$K) {
      for (item_step in 1:(test$items$pars$m[item] + 1)) {
        Psi <- c(1, logistic_function(inner_product_alpha_theta[item] - test$items$pars$beta[item, 1:test$items$pars$m[item]]), 0) # basically: 1, logistic_function(inner_product_alpha_theta - b), 0.
        D[item] <- D[item] + probabilities[item, item_step] * sum(Psi[2:(item_step + 1)] * (1 - Psi[2:(item_step+1)]))
      }
    }
    D
  }
  
  get_D_GPCM <- function() {
    # Generalized Partial Credit Model (Muraki, 1992)
    D <- numeric(test$items$K)
    for (item in 1:test$items$K) {
      mi <- 1:test$items$pars$m[item]
      pi <- probabilities[item, mi + 1] # remove j = 0, index is now also correct.
      mp <- sum(mi*pi)  
      D[item] <- sum((mi * pi) * (mi - mp))
    }
    D
  }
  
  validate <- function() {
    if (is.null(test))
      add_error("test", "is missing")
    if (is.null(person)) 
      add_error("person", "is missing")
  }
  
  invalid_result <- function() {
    list(errors = errors())
  }
  
  validate_and_run()
}