#' Subset method for ShadowCAT.items.
#' 
#' Creates a valid subset of items to work with in ShadowCAT functions. Primarily used internally.
#' @param items  item characteristics class ShadowCAT.items.
#' @param select indeces of items to be selected
#' @return Item subset of provided itembank, with class ShadowCAT.items, as well as any original classes. 
#' @export
subset.ShadowCAT.items <- function(items, select) {
  items$pars <- lapply(items$pars, 
                       FUN = function(x) {
                         if (is.vector(x)) 
                           x[select]
                         else
                           x[select,,drop = FALSE]
                       })
                        
  #items$pars <- lapply(items$pars, subset, subset = sapply(1:items$K, FUN = function (x) { if (x %in% select) TRUE else FALSE }))
  items$subset <- select
  items$K <- length(select)
  items
}
