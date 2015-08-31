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
    out <- which(objective == max(objective[person$available]))
    
    # return
    return(out)
  }


