#' get the number of non-NA cells per row of a matrix
#' 
#' @param X a matrix
#' @return number of non-NA cells per row
#' @examples number_non_missing_cells_per_row(matrix(c(1,2,NA,NA,2,3,NA,2,3), ncol = 3)) == c(1, 3, 2) || stop("wrong")
#' @export
number_non_missing_cells_per_row <- function(X) {
  apply(X, 1, function(x) sum(!is.na(x)))
}

#' apply with input and output converted to matrix
#' 
#' @param X a vector or matrix; if vector, the vector will be converted to a matrix
#' @param margin a vector giving the subscripts which the function will be applied over, see ?apply for details
#' @param FUN the funcion to be applied
#' @param ... any additional arguments to FUN
#' @return same a apply, but always in matrix format
#' @examples matrix_apply(c(1:5, NA), 2, sum, na.rm = TRUE) == 15 || is.matrix(matrix_apply(c(1:5, NA), 2, sum, na.rm = TRUE)) || stop("wrong")
#' @export
matrix_apply <- function(X, margin, FUN, ...) {
  as.matrix(apply(as.matrix(X), margin, FUN, ...))
}

#' compute cumulative sums for each row of a matrix
#' 
#' @param X a vector or matrix; if vector, the vector will be converted to a matrix with one column
#' @return matrix containing cumulative sums for each row
#' @examples row_cumsum(matrix(1:12, ncol = 3))[2,] == c(2, 6, 10) || stop("wrong");
#' row_cumsum(1:12) == matrix(1:12, ncol = 1) || stop("wrong")
#' @export
row_cumsum <- function(X) {
  transpose_if_ncol_and_nrow_larger_1(matrix_apply(X, 1, cumsum))
}

#' Recursive list apply, given a bunch of vectors, creates lists of lists with the elements as keys
#'  and values determined by fn.
#' @param ... vectors to loop over
#' @param fn function to execute for each leaf
#' @examples
#' rsapply(c('a', 'b'), c('f', 'g'), fn=function(x,y) c(x, y))
#' #> list(a=list(f=c('a', 'f'), g=c('a', 'g')), b=list(f=c('b', 'f'), g=c('b', 'g')))
#' @export
rsapply <- function(..., fn) {
  .rsapply <- function(el=NULL, ..., args=c(), fn) {
    if(is.null(el)) return(do.call(fn, as.list(args)))
    result <- list()
    for(x in el) {
      result[[x]] <- rsapply(..., args=c(args, x), fn=fn)
    }
    result
  }
  .rsapply(..., fn=fn)
}

#' get sum of objects returned by a loop, added to a starting object
#' 
#' @param start_object the object to which the sum is to be added
#' @param loop_vector the vector to be looped over
#' @param FUN a function of the loop_vector values, returning the object to be added to start_object at each iteration. First argument should be loop values argument
#' @param ... any additional arguments to FUN
#' @return sum of objects returned by a loop, added to starting object
#' @examples X = matrix(1:6, ncol = 2);
#' sum_loop_outputs(matrix(0, 2, 2), 1:3, FUN = function(item, X) { X[item,] %*% t(X[item,]) }, X = X) == (X[1,] %*% t(X[1,]) + X[2,] %*% t(X[2,]) + X[3,] %*% t(X[3,])) || stop("wrong")
#' @export
sum_loop_outputs <- function(start_object, loop_vector, FUN, ...) {
  for (i in loop_vector) {
    start_object <- start_object + FUN(i, ...)
  }
  start_object
}

#' transpose matrix if number of columns is larger than 1
#' 
#' @param X a matrix
#' @return transpose of X if number of columns and number of rows are larger than 1, X otherwise
#' @examples transpose_if_ncol_and_nrow_larger_1(matrix(1:50, ncol = 2)) == t(matrix(1:50, ncol = 2)) || stop("wrong")
#' @export
transpose_if_ncol_and_nrow_larger_1 <- function(X) {
  if (ncol(X) > 1 && ncol(X) > 1)
    t(X)
  else
    X
}


