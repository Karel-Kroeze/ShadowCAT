#' Next item
#' 
#' Control function to select the next item, takes account of starting / stopping conditions.
#' 
#' @param person
#' @param test
#' @return integer item index
#' @export
next_item <- function(person, test) {
  l <- length(person$responses)
  if (l < test$start$n) {
    # Pre-CAT
    # utter random selection
    if (test$start$type == 'random'){
      out <- sample(person$available, 1)
    }
    
    # set a pre-defined list of starting items (indices)
    if (test$start$type == 'fixed'){
      out <- test$start$indices[length(person$responses) + 1]
    }
    
    # picks nByDimension starting items per dimension (or n_i if nByDimension is a length Q vector), assumes between models, 
    # if any item has a non-zero loading on a dimension, it is considered to be part of that dimension. 
    # they CAN overlap, which may cause unwanted side effects, and in within models the result is identical to 'normal' random starting.
    if (test$start$type == 'randomByDimension'){
      # check for vector or scalar argument, make vector
      if (length(test$start$nByDimension) == 1) 
        test$start$nByDimension <- rep(test$start$nByDimension, test$items$Q)
      
      # validate n and nByDimension
      if (test$start$n != sum(test$start$nByDimension))
        stop("Total length of start phase and sum of length per dimension do not match (n != sum(nByDimension)")
      
      # create a design matrix of item loadings
      design <- test$items$pars$alpha > 0
      
      # loop over dims, see which we're in
      for (q in 1:test$items$Q) {
          if (l < sum(test$start$nByDimension[1:q]))
            break 
      }
      
      # get items for this dimension, cull non-applicable, and select at random.
      loading <- (1:test$items$K)[design[,q]]
      applicable <- intersect(loading, person$available)
      out <- sample(applicable, 1)
    }
    
  } else {
    # CAT 
    out <- best_item(person, test)
  }
  
  return(out)
}
