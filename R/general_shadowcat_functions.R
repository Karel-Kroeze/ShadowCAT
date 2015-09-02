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

#' apply with input and output converted to matrix
#' 
#' @param a vector or matrix; if vector, the vector will be converted to a matrix
#' @param margin a vector giving the subscripts which the function will be applied over, see ?apply for details
#' @param FUN the funcion to be applied
#' @param ... any additional arguments to FUN
#' @return same a apply, but always in matrix format
#' @examples matrix_apply(c(1:5, NA), 2, sum, na.rm = TRUE) == 15 || is.matrix(matrix_apply(c(1:5, NA), 2, sum, na.rm = TRUE)) || stop("wrong")
#' @export
matrix_apply <- function(X, margin, FUN, ...) {
  as.matrix(apply(as.matrix(X), margin, FUN, ...))
}