#' Simulate response(s)
#' 
#' Given a person and (a subset of) an itembank, simulate a response pattern.
#' 
#' 
#' @param person
#' @param test
#' @param indeces If not NULL, answer questions with the given indeces, and update the person object.
#' @return vector responses, or updated person object if indeces is set.
#' @export
answer <- function(person, test, indeces = NULL) {
  # subset items to relevant part
  # TODO: this is kind of awkward, we're updating the parent object with a temporary version of items. 
  # Doesn't hurt, but not pretty
  if ( ! is.null(indeces)) test$items <- subset(test$items, indeces)
  
  # attach items directly
  items <- test$items
  
  # probabilities, generated with TRUE theta.
  Pij <- prob(test,theta=person$theta)$P
  
  # cumulative probabilities
  cp <- Pij 
  for (i in 1:(items$M+1)) cp[,i] <- apply(matrix(Pij[,1:i],ncol=i),1,sum)
  
  # rand ~ unif(0,1)
  rand <- runif(items$K)
  
  # answer is the number of categories that have a cumulative probability smaller than rand
  responses <- apply(rand > cp, 1, sum, na.rm=TRUE)
  
  # what to return?
  if (is.null(indeces)){
    # vector of responses
    out <- responses
  } else {
    # updated person
    person$responses <- c(person$responses, responses)
    person$administered <- c(person$administered, indeces)
    person$available <- person$available[-which(person$available %in% indeces)]
    out <- person
  }
  
  #return
  return(out)
}