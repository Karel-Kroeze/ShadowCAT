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
#' row_cumsum(1:12) == 1:12 || stop("wrong")
#' @export
row_cumsum <- function(X) {
  transpose_if_ncol_and_nrow_larger_1(matrix_apply(X, 1, cumsum))
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




