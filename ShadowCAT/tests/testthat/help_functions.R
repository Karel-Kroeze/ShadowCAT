#' Run function with set.seed
#' 
#' Wraption given function to seeds random with given seed. Resets seed after.
#' 
#' @param seed The seed to use.
#' @param fn The function to call using the given seed.
#' @return Function that wraps the given function.
#' @examples
#' foo <- with_random_seed(2, function(a) { runif(a) })
#' foo(2) == foo(2)
with_random_seed <- function(seed, fn) {
  function(...) {
    old <- .Random.seed
    set.seed(seed)
    res <- fn(...)
    .Random.seed <<- old
    res
  }
}