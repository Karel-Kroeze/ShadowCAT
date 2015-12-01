#' Simulate response(s)
#' 
#' Given a person and (a subset of) an itembank, simulate a response pattern.
#' 
#' @examples 
#' items <- createTestBank("GPCM")
#' test <- initTest(items)
#' person <- initPerson(items)
#' 
#' # simulates responses to all questions, and returns a vector response pattern.
#' answer(person, test)
#' 
#' # simulates responses to the specified question indeces, and returns an updated person object.
#' answer(person, test, sample(test$items$K, 5))
#' 
#' 
#' @param person Person object, see \code{\link{initPerson}}.
#' @param test Test object, see \code{\link{initPerson}}.
#' @param indeces If not NULL, answer questions with the given indeces, and update the person object. If NULL, answers the full test.
#' @return vector responses, or updated person object if indeces is set.
#' @export
answer <- function(person, test, indeces) {
  # probabilities, generated with TRUE theta.
  Pij <- probabilities_and_likelihoods(person$theta, test$items$model, indeces, test$items$Q, test$estimator, test$items$pars$alpha, test$items$pars$beta, test$items$pars$guessing, output = "probs")
  
  # cumulative probabilities
  cp <- Pij 
  for (i in 1:(test$items$M+1)) cp[,i] <- apply(matrix(Pij[,1:i],ncol=i),1,sum)
  
  # rand ~ unif(0,1)
  rand <- runif(length(indeces))
  
  # answer is the number of categories that have a cumulative probability smaller than rand
  responses <- apply(rand > cp, 1, sum, na.rm=TRUE)
  
  # updated person
  person$responses <- c(person$responses, responses)
  person$administered <- c(person$administered, indeces)
  person$available <- person$available[-which(person$available %in% indeces)]
  person
}