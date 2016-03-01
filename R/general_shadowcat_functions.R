#' Return as a scalar
#' NULL -> NA
#'
#' @param obj Any object
#' @return NA (if NULL) or obj as a scalar
#' @export
as.scalar2 <- function(obj){
  if (is.null(obj)) obj <- NA
  return(structure(obj, class=c("scalar",class(obj))))
}

#' transform categorical vector into dummy variables
#' 
#' @param unique_values a vector containing the unique values of the categorical vector
#' @param categorical_vector the categorical vector to be transformed into dummy variables
#' @return a matrix containing the dummy variables, number of columns (dummies) equals the number of unique values
#' @examples expect_equal <- categorical_to_dummy(c("a", "b", "c"), rep(c("a", "b", "c"), 3))[1,] == rep(c(1, 0, 0), 3);
#' sum(expect_equal) == 9 || stop("wrong")
#' @export
categorical_to_dummy <- function(unique_values, categorical_vector) {
  sapply(unique_values, FUN = function(value) { as.numeric(categorical_vector == value) })
}

#' Get subset based on indeces
#' 
#' @param x vector or matrix from which a subset is to be taken
#' @param subset indeces of elements/rows to be selected
#' @return subset of x 
#' @examples all(get_subset(c(1, 4, 2, 6, 3), c(2, 5, 3)) == c(4, 3, 2)) || stop("wrong")
#' all(get_subset(matrix(c(1, 4, 2, 6, 3, 8, 6, 9, 0, 1), ncol = 2), c(2, 5, 3)) == matrix(c(4, 3, 2, 6, 1, 9), ncol = 2)) || stop("wrong")
#' @export
get_subset <- function(x, subset) {
  if (is.null(x))
    return(NULL)
  if (is.vector(x)) 
    x[subset]
  else
    x[subset,,drop = FALSE]
}

#' lapply with array as output
#' 
#'@param x a vector or list; see ? lapply for details
#'@param dim vector with dimensions of the returned array
#'@param FUN function to be applied to each element of X; see ? lapply for details
#'@param ... any additional arguments to FUN
#'@return lapply output, with output converted to array
#'@examples matrix_example <- matrix(1:9, ncol = 3); 
#' lapply_return_array(1:3, 
#'                     c(3, 3, 3), 
#'                     FUN = function(i, matrix_example) { matrix_example[i,] %*% t(matrix_example[i,]) }, 
#'                     matrix_example = matrix_example)[1, 1, 1] == 1 || stop("wrong")
#'@export
lapply_return_array <- function(x, dim, FUN, ...) {
  lapply_output <- lapply(x, FUN, ...)
  array(unlist(lapply_output), dim = dim)
}

#' apply with input and output converted to matrix
#' 
#' @param x a vector or matrix; if vector, the vector will be converted to a matrix
#' @param margin a vector giving the subscripts which the function will be applied over, see ?apply for details
#' @param FUN the funcion to be applied
#' @param ... any additional arguments to FUN
#' @return same a apply, but always in matrix format
#' @examples ( matrix_apply(c(1:5, NA), 2, sum, na.rm = TRUE) == 15 && is.matrix(matrix_apply(c(1:5, NA), 2, sum, na.rm = TRUE)) ) || stop("wrong")
#' @export
matrix_apply <- function(x, margin, FUN, ...) {
  as.matrix(apply(as.matrix(x), margin, FUN, ...))
}

#' not in
#'
#'@param x vector or value for which it is too be checked if its elements are in vector y
#'@param y vector for which it is too be checked if it contains the values in x
#'@return a vector with elements equal to TRUE for each element of x that is not in y
#'@examples expect_equal <- c("a", "b", "g", "k", "b") %not_in% c("a", "b", "c") == c(FALSE, FALSE, TRUE, TRUE, FALSE);
#' sum(expect_equal) == 5 || stop("wrong")
#'@export
`%not_in%` <- function(x,y) { !(x %in% y) }

#' get the number of non-NA cells per row of a matrix
#' 
#' @param x a matrix
#' @return number of non-NA cells per row
#' @examples expect_equal <- number_non_missing_cells_per_row(matrix(c(1,2,NA,NA,2,3,NA,2,3), ncol = 3)) == c(1, 3, 2);
#' sum(expect_equal) == 3 || stop("wrong")
#' @export
number_non_missing_cells_per_row <- function(x) {
  apply(x, 1, function(x_row) sum(!is.na(x_row)))
}

#' compute cumulative sums for each row of a matrix
#' 
#' @param x a vector or matrix; if vector, the vector will be converted to a matrix with one column
#' @return matrix containing cumulative sums for each row
#' @examples expect_equal1 <- row_cumsum(matrix(1:12, ncol = 3))[2,] == c(2, 8, 18);
#' expect_equal2 <- row_cumsum(1:12) == matrix(1:12, ncol = 1);
#' expect_equal3 <- row_cumsum(matrix(1:4, nrow = 1)) == matrix(c(1, 3, 6, 10), nrow = 1);
#' sum(expect_equal1) == 3 || stop("wrong")
#' sum(expect_equal2) == 12 || stop("wrong")
#' sum(expect_equal3) == 4 || stop("wrong")
#' @export
row_cumsum <- function(x) {
  x_matrix <- as.matrix(x)
  r_cumsum <- matrix_apply(x_matrix, 1, cumsum)
  if ((nrow(x_matrix) > 1 && ncol(x_matrix) > 1) || nrow(x_matrix) == 1)
    t(r_cumsum)
  else
    r_cumsum
}

#' get row sums of matrix or sum of vector
#' 
#' @param x matrix of which row sums are to be obtained, or vector of which sum is to be obtained
#' @return row sums or sum of x
#' @examples all(row_or_vector_sums(matrix(c(1:10), nrow = 2)) == c(25, 30)) || stop("wrong");
#' row_or_vector_sums(1:5) == 15 || stop("wrong")
#' @export
row_or_vector_sums <- function(x) {
  if (is.matrix(x))
    rowSums(x)
  else
    sum(x)  
}

#' Check whether all matrices have equal row names, in same order
#' 
#' @param row_names the correct character vector of row names
#' @param list_of_matrices_to_check list containing the matrices for which the row names should be checked. NULL elements in the list are ignored.
#' @return TRUE if the row names of all the matrices are equal to row_names, including the order; FALSE otherwise
#' @examples row_names_are_equal(c("a", "b", "c"), list("matrix1" = matrix(1:3, ncol = 1, dimnames = list(c("a", "b", "c"), NULL)),
#'                                                      "matrix2" = matrix(2:4, ncol = 1, dimnames = list(c("a", "b", "c"), NULL)),
#'                                                      "matrix3" = matrix(3:5, ncol = 1, dimnames = list(c("a", "b", "c"), NULL)),
#'                                                      "matrix4" = NULL)) == TRUE || stop("wrong")
#' row_names_are_equal(c("a", "b", "c"), list("matrix1" = matrix(1:3, ncol = 1, dimnames = list(c("a", "b", "c"), NULL)),
#'                                            "matrix2" = matrix(2:4, ncol = 1, dimnames = list(c("a", "c", "b"), NULL)),
#'                                            "matrix3" = matrix(3:5, ncol = 1, dimnames = list(c("a", "b", "c"), NULL)))) == FALSE || stop("wrong")
#' @export
row_names_are_equal <- function(row_names, list_of_matrices_to_check) {
  row_names_are_equal_per_matrix <- sapply(1:length(list_of_matrices_to_check), 
                                           function(x) { 
                                             if (is.null(x))
                                               return(TRUE)
                                             all(rownames(list_of_matrices_to_check[[x]]) == row_names) })
  all(row_names_are_equal_per_matrix)
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
#' @examples matrix_example <- matrix(1:6, ncol = 2);
#' expect_equal <- sum_loop_outputs(matrix(0, 2, 2), 1:3, FUN = function(item, matrix_example) { matrix_example[item,] %*% t(matrix_example[item,]) }, matrix_example = matrix_example) == 
#'                 (matrix_example[1,] %*% t(matrix_example[1,]) + matrix_example[2,] %*% t(matrix_example[2,]) + matrix_example[3,] %*% t(matrix_example[3,]));
#' sum(expect_equal) == 4 || stop("wrong")
#' @export
sum_loop_outputs <- function(start_object, loop_vector, FUN, ...) {
  for (i in loop_vector) {
    start_object <- start_object + FUN(i, ...)
  }
  start_object
}
