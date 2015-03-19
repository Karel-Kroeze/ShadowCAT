#' Maximum Information item selection
#' 
#' Naive item selection based on maximum information only.
#' 
#' Selects the item with the highest value of a specified objective function.
#' @param test
#' @param person
#' @return item index of item with the highest value of the objective function.
#' @export
MI <- function(test, person, objective){
    
    # fetch highest information item
    max_info <- 0
    
    # in case nothing gets flagged, preselect a random item
    max_item <- sample(person$available, 1)
    
    # loop over items, flag best matches.
    for (i in 1:length(objective)) { 
      if (objective[i] > max_info){
        max_info <- objective[i]
        max_item <- person$available[i]
      }
    }
    
    # return
    return(max_item)
  }


