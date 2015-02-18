#' Simulate response(s)
#' 
#' Given a person and (a subset of) an itembank, simulate a response pattern.
#' 
#' 
#' @param person
#' @param items
#' @return vector responses
#' @export
answer <- function(person, items) {
  # probabilities, generated with TRUE theta.
  Pij <- prob(items,theta=person$theta)$P
  
  # cumulative probabilities
  cp <- Pij 
  for (i in 1:(items$M+1)) cp[,i] <- apply(matrix(Pij[,1:i],ncol=i),1,sum)
  
  # rand ~ unif(0,1)
  rand <- runif(items$K)
  
  # answer is the number of categories that have a cumulative probability smaller than rand
  out <- apply(rand > cp, 1, sum, na.rm=TRUE)
  
  #return
  return(out)
}