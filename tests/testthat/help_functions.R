#' Wraption given function to seeds random with given seed. Resets seed after.
#' @param seed the seed to use
#' @param fn the function to call using the given seed
#' @return function that wraps the given function
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