#' Subset method for ShadowCAT.items.
#' 
#' Creates a valid subset of items to work with in ShadowCAT functions. Primarily used internally.
#' @param x item characteristics class ShadowCAT.items (as returned by initItembank()).
#' @param subset indeces of items to be selected
#' @param ... further arguments to be passed to or from other methods.
#' @return Item subset of provided itembank, with class ShadowCAT.items, as well as any original classes. 
#' @export
subset.ShadowCAT.items <- function(x, subset, ...) {
  x$pars <- lapply(x$pars, 
                       FUN = function(x) {
                         if (is.vector(x)) 
                           x[subset]
                         else
                           x[subset,,drop = FALSE]
                       })
                        
  x$subset <- subset
  x$K <- length(subset)
  x
}

