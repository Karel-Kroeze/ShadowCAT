#' Maximum Information item selection
#' 
#' Naive item selection based on maximum information only.
#' 
#' Selects the item with the highest value of a specified objective function.
#' @param test
#' @param person
#' @return item index of item with the highest value of the objective function.
#' @export
MI <- function(test, person){
  if (test$objective == "FI" | test$objective == "PFI") {
    # fetch info for all items
    info <- FI(test, person)
    
    # sum over administered items, if there are any. If none, set this to 0 matrix. 
    if (length(person$administered > 0)){
      so_far <- apply(info[,,person$administered], 3, sum)  
    } else {
      so_far <- matrix(0,test$items$Q,test$items$Q)
    }
    
    # include prior?
    if (test$objective == "PFI") so_far <- so_far + solve(person$prior)
    
    # fetch highest information item
    max_info <- 0
    for (i in person$available) { 
      item_info <- det(so_far + info[,,i])
      if (item_info > max_info){
        max_info <- item_info
        max_item <- i
      }
    }
    
    # return
    return(max_item)
  }
}