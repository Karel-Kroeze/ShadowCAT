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
    
    # fetch highest information available item
    # TODO: This is a very hacky way to stop itemselection on domains that have met variance targets.
    # Needs simplification, optimization and generalization.
    if (test$stop$type == "variance" && test$items$Q > 1){
      # check which domains are completed (should really be held in person object)
      completed_domains <- which(diag(attr(person$estimate, "variance")) < test$stop$target)
      
      # check which items load on which domains (should really be held in test object)
      item_loadings <- apply( test$items$pars$alpha, 1, function(x) abs(x) > 1e-5 )
      
      # cross-reference item loadings with completed domains, throwing out items where are loadings are on completed domains
      useful_items <- which( apply(item_loadings, 2, function(x) !all(which(x) %in% completed_domains)))
      
      # return max item(s) that have not yet been administered, and have 'useful' domain loadings
      out <- which(objective == max(objective[intersect(useful_items, person$available)]))
    } else {
      out <- which(objective == max(objective[person$available]))
    }
    
    # return
    return(out)
  }


