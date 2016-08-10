#' Return as a scalar
#' 
#' Same as \code{opencpu.encode::as.scalar}, except that it also turns \code{NULL} into \code{NA}. 
#'
#' @param obj Any object.
#' @return \code{obj} as a scalar; \code{NA} if \code{obj} is \code{NULL}.
as.scalar2 <- function(obj){
  if (is.null(obj)) obj <- NA
  return(structure(obj, class=c("scalar",class(obj))))
}

#' Transform categorical vector into dummy variables
#' 
#' Transform categorical vector into dummy variables, with number of dummy
#' variables equal to the number of unique values in the categorical vector. 
#' 
#' @param x Categorical vector to be transformed into dummy variables.
#' @return Matrix containing the dummy variables.
#' @examples x <- rep(c("a", "b", "c"), 3)
#' ShadowCAT:::categorical_to_dummy(x)
#' \dontshow{
#' all(ShadowCAT:::categorical_to_dummy(x)[1,] == rep(c(1, 0, 0), 3)) || stop("wrong")
#' }
categorical_to_dummy <- function(x) {
  sapply(unique(x), FUN = function(value) { as.numeric(x == value) })
}

#' Get subset based on indices
#' 
#' Get a subset of a vector or matrix based on a vector of (row) indices.
#' 
#' @param x Vector or matrix from which a subset is to be taken.
#' @param subset Vector with indices of elements/rows to be selected.
#' @return Subset of x.
#' @examples ShadowCAT:::get_subset(c(1, 4, 2, 6, 3), subset = c(2, 5, 3))
#' ShadowCAT:::get_subset(matrix(c(1, 4, 2, 6, 3, 8, 6, 9, 0, 1), ncol = 2), subset = c(2, 5, 3))
#' \dontshow{
#' all(ShadowCAT:::get_subset(c(1, 4, 2, 6, 3), subset = c(2, 5, 3)) == c(4, 3, 2)) || stop("wrong")
#' all(ShadowCAT:::get_subset(matrix(c(1, 4, 2, 6, 3, 8, 6, 9, 0, 1), ncol = 2), subset = c(2, 5, 3)) == matrix(c(4, 3, 2, 6, 1, 9), ncol = 2)) || stop("wrong")
#' }
get_subset <- function(x, subset) {
  if (is.null(x))
    return(NULL)
  if (is.vector(x)) 
    x[subset]
  else
    x[subset,,drop = FALSE]
}

#' Lapply with array output
#' 
#' Same as lapply, but output of lapply is converted to array.
#' 
#' @param x A vector or list; see \code{\link{lapply}} for details.
#' @param dim Vector with dimensions of the returned array.
#' @param FUN Function to be applied to each element of x; see \code{\link{lapply}} for details.
#' @param ... Any additional arguments to FUN.
#' @return As for \code{\link{lapply}}, but with output converted to array.
#' @examples matrix_example <- matrix(1:9, ncol = 3)
#' ShadowCAT:::lapply_return_array(x = 1:3, 
#'                                 dim = c(3, 3, 3), 
#'                                 FUN = function(i, matrix_example) { 
#'                                         matrix_example[i,] %*% t(matrix_example[i,]) 
#'                                       }, 
#'                                 matrix_example = matrix_example)
#' \dontshow{
#' array_out <- ShadowCAT:::lapply_return_array(x = 1:3, 
#'                                             dim = c(3, 3, 3), 
#'                                             FUN = function(i, matrix_example) { 
#'                                                     matrix_example[i,] %*% t(matrix_example[i,]) }, 
#'                                             matrix_example = matrix_example)
#' list_out <- lapply(X = 1:3,
#'                    FUN = function(i, matrix_example) { 
#'                            matrix_example[i,] %*% t(matrix_example[i,]) },
#'                    matrix_example = matrix_example)
#' all(array_out[,,1] == list_out[[1]] && array_out[,,2] == list_out[[2]] && array_out[,,3] == list_out[[3]]) || stop("wrong")
#' }
lapply_return_array <- function(x, dim, FUN, ...) {
  lapply_output <- lapply(x, FUN, ...)
  array(unlist(lapply_output), dim = dim)
}

#' Matrix apply
#' 
#' Apply with input and output converted to matrix
#' 
#' @param x Vector or matrix; if vector, the vector will be converted to a matrix.
#' @param margin Vector giving the subscripts the function will be applied over, see \code{\link{apply}} for details.
#' @param FUN The funcion to be applied, see \code{\link{apply}} for details.
#' @param ... Any additional arguments to FUN.
#' @return Output of apply, but always in matrix format.
#' @examples ShadowCAT:::matrix_apply(c(1:5, NA), 2, sum, na.rm = TRUE)
#' \dontshow{
#' ShadowCAT:::matrix_apply(c(1:5, NA), 2, sum, na.rm = TRUE) == matrix(15) || stop("wrong")
#' }
matrix_apply <- function(x, margin, FUN, ...) {
  as.matrix(apply(as.matrix(x), margin, FUN, ...))
}

#' Move values towards means
#'
#' Move a vector of values in the direction of a vector of means.
#' 
#' @param values The vector with values to move to the mean.
#' @param means The vector with means to which values should be moved.
#' @param amount_change Vector with the amounts by which each value should be moved to its corresponding mean.
#' @return \code{values} plus or minus \code{amount_change}, or \code{values} if \code{values} is equal to \code{means}.
#' @examples ShadowCAT:::move_values_to_means(values = 1, means = 2, amount_change = .1)
#' ShadowCAT:::move_values_to_means(values = 3, means = 2, amount_change = .1)
#' ShadowCAT:::move_values_to_means(values = 2, means = 2, amount_change = .1)
#' ShadowCAT:::move_values_to_means(values = c(1, 2, 3), 
#'                                  means = c(1, 5, -1), 
#'                                  amount_change = rep(.1, 3))
#' 
#' \dontshow{
#' ShadowCAT:::move_values_to_means(values = 1, means = 2, amount_change = .1) == 1.1 || stop("wrong") 
#' ShadowCAT:::move_values_to_means(values = 3, means = 2, amount_change = .1) == 2.9 || stop("wrong")
#' ShadowCAT:::move_values_to_means(values = 2, means = 2, amount_change = .1) == 2 || stop("wrong")
#' all(ShadowCAT:::move_values_to_means(values = c(1, 2, 3), means = c(1, 5, -1), amount_change = rep(.1, 3)) == c(1, 2.1, 2.9)) || stop("wrong")
#' }
move_values_to_means <- function(values, means, amount_change) {
  sapply(1:length(values), function(i) {
    if (values[i] == means[i])
      return(values[i])
    if (values[i] < means[i])
      values[i] + amount_change[i]
    else
      values[i] - amount_change[i]
  } )
}


#' Check \code{NA} pattern in vector
#' 
#' Check whether \code{NA}'s exist only at the end of a vector, that is, 
#' whether there are no values after the first \code{NA} in the vector.
#' 
#' @param x Vector.
#' @return \code{TRUE} if there are only \code{NA}'s at the end of the vector or no \code{NA} at all; \code{FALSE} otherwise.
#' @examples ShadowCAT:::na_only_end_vector(c(1:3, NA, NA))
#' ShadowCAT:::na_only_end_vector(c(1:3, NA, NA, 6))
#' 
#' \dontshow{
#' ShadowCAT:::na_only_end_vector(c(1:3, NA, NA)) || stop("wrong");
#' ShadowCAT:::na_only_end_vector(c(1:3, NA, NA, 6)) == FALSE || stop("wrong")
#' }
na_only_end_vector <- function(x) {
   all(diff(is.na(x)) >= 0)
}

#' Check \code{NA} pattern in matrix
#'
#' Check whether \code{NA}'s exist only at the end of each row, that is, 
#' whether there are no values after the first \code{NA} in the row.
#' 
#' @param x Matrix.
#' @return \code{TRUE} if there are only \code{NA}'s at the end of each row or no \code{NA} at all; \code{FALSE} otherwise.
#' @examples ShadowCAT:::na_only_end_rows(matrix(c(1:3, NA, NA, 1:5, 1:4, NA), ncol = 5, byrow = TRUE))
#' ShadowCAT:::na_only_end_rows(matrix(c(1:3, NA, NA, 1:2, NA, 4:5, 1:4, NA), ncol = 5, byrow = TRUE))
#' 
#' \dontshow{
#' ShadowCAT:::na_only_end_rows(matrix(c(1:3, NA, NA, 1:5, 1:4, NA), ncol = 5, byrow = TRUE)) || stop("wrong");
#' ShadowCAT:::na_only_end_rows(matrix(c(1:3, NA, NA, 1:2, NA, 4:5, 1:4, NA), ncol = 5, byrow = TRUE))  == FALSE || stop("wrong")
#' }
na_only_end_rows <- function(x) {
  all(apply(x, 1, na_only_end_vector))
}

#' Not in
#' 
#' Get logical vector indicating if there no match for the left operant in the right operant.
#'
#' @param x vector or value for which it is too be checked if its elements are in vector y
#' @param y vector for which it is too be checked if it contains the values in x
#' @return Vector of same length as that of \code{x}, with elements equal to \code{TRUE} for each element of x that is not in y.
#' @examples ShadowCAT:::`%not_in%`(c("a", "b", "g", "k", "b"), c("a", "b", "c"))
#' \dontshow{
#' all(ShadowCAT:::`%not_in%`(c("a", "b", "g", "k", "b"), c("a", "b", "c")) == c(FALSE, FALSE, TRUE, TRUE, FALSE)) || stop("wrong")
#' }
`%not_in%` <- function(x,y) { !(x %in% y) }


#' Count non-\code{NA} per row
#'
#' Get the number of non-\code{NA} cells per row of a matrix.
#' 
#' @param x Matrix.
#' @return Vector with the number of non-\code{NA} cells per row.
#' @examples ShadowCAT:::number_non_missing_cells_per_row(matrix(c(1,2,NA,NA,2,3,NA,2,3), ncol = 3))
#' 
#' \dontshow{
#' all(ShadowCAT:::number_non_missing_cells_per_row(matrix(c(1,2,NA,NA,2,3,NA,2,3), ncol = 3)) == c(1, 3, 2)) || stop("wrong")
#' }
number_non_missing_cells_per_row <- function(x) {
  apply(x, 1, function(x_row) sum(!is.na(x_row)))
}

#' Remove rows outside bounds
#' 
#' Remove rows of a matrix that contain one or more values outside specified bounds.
#' 
#' @param matrix_to_evaluate Matrix for which the rows with values outside the bounds should be removed.
#' @param lower_bound Vector with lower bounds for each column; length should be equal to the number of columns of matrix_to_evaluate.
#' @param upper_bound Vector with upper bounds for each column; length should be equal to the number of columns of matrix_to_evaluate.
#' @return \code{matrix_to_evaluate} without the rows that contain values outside the specified bounds.
#' @examples ShadowCAT:::remove_rows_outside_bounds(matrix(-3:3), lower_bound = -2, upper_bound = 2)
#' ShadowCAT:::remove_rows_outside_bounds(expand.grid(list(-3:3, -2:2)), 
#'                                        lower_bound = c(-2, -3), 
#'                                        upper_bound = c(2, 1))
#' 
#' \dontshow{
#' all(ShadowCAT:::remove_rows_outside_bounds(matrix(-3:3), lower_bound = -2, upper_bound = 2) == matrix(-2:2)) || stop("wrong")
#' matrix_to_evaluate <- expand.grid(list(-3:3, -2:2))
#' all(ShadowCAT:::remove_rows_outside_bounds(matrix_to_evaluate, lower_bound = c(-2, -3), upper_bound = c(2, 1)) == matrix_to_evaluate[c(2:6, 9:13, 16:20, 23:27), ]) || stop("wrong")
#' }
remove_rows_outside_bounds <- function(matrix_to_evaluate, lower_bound, upper_bound) {
  matrix_to_evaluate_transpose <- t(matrix_to_evaluate)
  inside_bounds <- colSums(matrix_to_evaluate_transpose <= upper_bound & matrix_to_evaluate_transpose >= lower_bound) == ncol(matrix_to_evaluate)
  matrix_to_evaluate[inside_bounds, , drop = FALSE]
}


#' Cumulative row sums
#' 
#' Compute cumulative sums for each row of a matrix.
#' 
#' @param x Vector or matrix. If vector, the vector will be converted to a matrix with one column.
#' @return Matrix containing cumulative sums for each row.
#' @examples ShadowCAT:::row_cumsum(matrix(1:12, ncol = 3))
#' ShadowCAT:::row_cumsum(1:12)
#' ShadowCAT:::row_cumsum(matrix(1:4, nrow = 1))
#' 
#' \dontshow{
#' all(ShadowCAT:::row_cumsum(matrix(1:12, ncol = 3))[2,] == c(2, 8, 18)) || stop("wrong")
#' all(should_be_equal2 <- ShadowCAT:::row_cumsum(1:12) == matrix(1:12, ncol = 1)) || stop("wrong")
#' all(should_be_equal3 <- ShadowCAT:::row_cumsum(matrix(1:4, nrow = 1)) == matrix(c(1, 3, 6, 10), nrow = 1)) || stop("wrong")
#' }
row_cumsum <- function(x) {
  x_matrix <- as.matrix(x)
  r_cumsum <- matrix_apply(x_matrix, 1, cumsum)
  if ((nrow(x_matrix) > 1 && ncol(x_matrix) > 1) || nrow(x_matrix) == 1)
    t(r_cumsum)
  else
    r_cumsum
}

#' Row sum or vector sum
#'
#' Get row sums of matrix or sum of vector.
#' 
#' @param x Matrix of which row sums are to be obtained, or vector of which sum is to be obtained.
#' @return Row sums or sum of x.
#' @examples ShadowCAT:::row_or_vector_sums(matrix(c(1:10), nrow = 2))
#' ShadowCAT:::row_or_vector_sums(1:5)
#' 
#' \dontshow{
#' all(ShadowCAT:::row_or_vector_sums(matrix(c(1:10), nrow = 2)) == c(25, 30)) || stop("wrong");
#' ShadowCAT:::row_or_vector_sums(1:5) == 15 || stop("wrong")
#' }
row_or_vector_sums <- function(x) {
  if (is.matrix(x))
    rowSums(x)
  else
    sum(x)  
}

#' Check row names
#'
#' Check whether all matrices have row names equal to the standard.
#' 
#' @param row_names Character vector of row names to use as the standard.
#' @param list_of_matrices_to_check List containing the matrices for which the row names should be checked. NULL elements in the list are ignored.
#' @return \code{TRUE} if the row names of all the matrices are equal to \code{row_names}, including the order; \code{FALSE} otherwise.
#' @examples 
#' ShadowCAT:::row_names_are_equal(c("a", "b", "c"), 
#'                                 list("matrix1" = matrix(1:3, 
#'                                                         ncol = 1, 
#'                                                         dimnames = list(c("a", "b", "c"), NULL)),
#'                                      "matrix2" = matrix(2:4, 
#'                                                         ncol = 1, 
#'                                                         dimnames = list(c("a", "b", "c"), NULL)),
#'                                      "matrix3" = matrix(3:5, 
#'                                                         ncol = 1, 
#'                                                         dimnames = list(c("a", "b", "c"), NULL)),
#'                                      "matrix4" = NULL))
#' ShadowCAT:::row_names_are_equal(c("a", "b", "c"), 
#'                                 list("matrix1" = matrix(1:3, 
#'                                                         ncol = 1, 
#'                                                         dimnames = list(c("a", "b", "c"), NULL)),
#'                                      "matrix2" = matrix(2:4, 
#'                                                         ncol = 1, 
#'                                                         dimnames = list(c("a", "b", "c"), NULL)),
#'                                      "matrix3" = matrix(3:4, 
#'                                                         ncol = 1, 
#'                                                         dimnames = list(c("a", "b"), NULL)),
#'                                      "matrix4" = NULL))
#' ShadowCAT:::row_names_are_equal(c("a", "b", "c"), 
#'                                 list("matrix1" = matrix(1:3, 
#'                                                         ncol = 1, 
#'                                                         dimnames = list(c("a", "b", "c"), NULL)),
#'                                      "matrix2" = matrix(2:4, 
#'                                                         ncol = 1, 
#'                                                         dimnames = list(c("a", "c", "b"), NULL)),
#'                                      "matrix3" = matrix(3:5, 
#'                                                         ncol = 1, 
#'                                                         dimnames = list(c("a", "b", "c"), NULL))))
#'
#' 
#' \dontshow{
#' ShadowCAT:::row_names_are_equal(c("a", "b", "c"), list("matrix1" = matrix(1:3, ncol = 1, dimnames = list(c("a", "b", "c"), NULL)),
#'                                                                  "matrix2" = matrix(2:4, ncol = 1, dimnames = list(c("a", "b", "c"), NULL)),
#'                                                                  "matrix3" = matrix(3:5, ncol = 1, dimnames = list(c("a", "b", "c"), NULL)),
#'                                                                  "matrix4" = NULL)) == TRUE || stop("wrong")
#' ShadowCAT:::row_names_are_equal(c("a", "b", "c"), list("matrix1" = matrix(1:3, ncol = 1, dimnames = list(c("a", "b", "c"), NULL)),
#'                                                        "matrix2" = matrix(2:4, ncol = 1, dimnames = list(c("a", "b", "c"), NULL)),
#'                                                        "matrix3" = matrix(3:4, ncol = 1, dimnames = list(c("a", "b"), NULL)),
#'                                                        "matrix4" = NULL)) == FALSE || stop("wrong")
#' ShadowCAT:::row_names_are_equal(c("a", "b", "c"), list("matrix1" = matrix(1:3, ncol = 1, dimnames = list(c("a", "b", "c"), NULL)),
#'                                                        "matrix2" = matrix(2:4, ncol = 1, dimnames = list(c("a", "c", "b"), NULL)),
#'                                                        "matrix3" = matrix(3:5, ncol = 1, dimnames = list(c("a", "b", "c"), NULL)))) == FALSE || stop("wrong")
#' }
row_names_are_equal <- function(row_names, list_of_matrices_to_check) {
  row_names_are_equal_per_matrix <- sapply(list_of_matrices_to_check, 
                                           function(matrix_to_check) { 
                                             if (is.null(matrix_to_check))
                                               return(TRUE)
                                             nrow(matrix_to_check) == length(row_names) && all(rownames(matrix_to_check) == row_names) })
  all(row_names_are_equal_per_matrix)
}

#' Recursive list apply
#' 
#' Given a bunch of vectors, creates lists of lists with the elements as keys and values determined by fn.
#' 
#' @param ... Vectors to loop over.
#' @param fn Function to execute for each leaf.
#' @return Lists of lists with the elements as keys and values determined by fn.
#' @examples
#' ShadowCAT:::rsapply(c('a', 'b'), c('f', 'g'), fn=function(x,y) c(x, y))
#' #> list(a=list(f=c('a', 'f'), g=c('a', 'g')), b=list(f=c('b', 'f'), g=c('b', 'g')))
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

#' Sum loop outputs
#' 
#' Get sum of objects returned by a loop, added to a starting object.
#' 
#' @param start_object The object to which the sum is to be added.
#' @param loop_vector The vector to be looped over.
#' @param FUN Function of the loop_vector values, returning the object to be added to \code{start_object} at each iteration. First argument should be loop values argument.
#' @param ... Any additional arguments to FUN.
#' @return Sum of objects returned by a loop, added to starting object.
#' @examples matrix_example <- matrix(1:6, ncol = 2)
#' ShadowCAT:::sum_loop_outputs(start_object = matrix(0, 2, 2), loop_vector = 1:3, 
#'                              FUN = function(item, matrix_example) { 
#'                                      matrix_example[item,] %*% t(matrix_example[item,]) }, 
#'                              matrix_example = matrix_example)
#' 
#'  \dontshow{
#' all(ShadowCAT:::sum_loop_outputs(start_object = matrix(0, 2, 2), loop_vector = 1:3, FUN = function(item, matrix_example) { matrix_example[item,] %*% t(matrix_example[item,]) }, matrix_example = matrix_example) == 
#'                    (matrix_example[1,] %*% t(matrix_example[1,]) + matrix_example[2,] %*% t(matrix_example[2,]) + matrix_example[3,] %*% t(matrix_example[3,]))) || stop("wrong")
#' }
sum_loop_outputs <- function(start_object, loop_vector, FUN, ...) {
  for (i in loop_vector) {
    start_object <- start_object + FUN(i, ...)
  }
  start_object
}

# validate functions constraints and characteristics

#' Check both \code{NULL} or not \code{NULL}
#' 
#' Check whether characteristics and constraints are either both \code{NULL} or both not \code{NULL}.
#' 
#' @param characteristics Data frame with characteristics.
#' @param constraints List of constraints.
#' @return \code{TRUE} if \code{characteristics} and \code{constraints} are either both \code{NULL} or both not \code{NULL}, \code{FALSE} otherwise.
no_missing_information <- function(characteristics, constraints) {
  (is.null(characteristics) && is.null(constraints)) || (!is.null(characteristics) && !is.null(constraints))
}

#' Check class and number of rows of characteristics
#'
#' Check whether characteristics is data frame with correct number of rows.
#' 
#' @param characteristics Data frame with characteristics.
#' @param number_items Number of items in the item bank.
#' @return \code{TRUE} if characteristics is data frame with correct number of rows, \code{FALSE} otherwise.
characteristics_correct_format <- function(characteristics, number_items) {
  if (is.null(characteristics))
    return(TRUE)
  is.data.frame(characteristics) && nrow(characteristics) == number_items
}

#' Check structure of constraints
#'
#' Check whether elements of constraints have correct structure.
#' 
#' @param constraints List of constraints.
#' @return \code{TRUE} if elements are lists of length three with correct names, \code{FALSE} otherwise.
constraints_correct_structure <- function(constraints) {
  if (is.null(constraints))
    return(TRUE)
  all(sapply(constraints, 
             function(constraint) { 
               is.list(constraint) && length(constraint) == 3 && names(constraint) == c("name", "op", "target") }))
}

#' Check names in constraints
#'
#' Check whether the name elements within the constraints elements are correct.
#' 
#' @param characteristics Data frame with characteristics.
#' @param constraints List of constraints.
#' @return \code{TRUE} if the name elements in constraints are correct, \code{FALSE} otherwise.
constraints_correct_names <- function(constraints, characteristics) {
  if (is.null(constraints) || is.null(characteristics) || !constraints_correct_structure(constraints) || !is.data.frame(characteristics))
    return(TRUE)
  characteristic_names <- unlist(sapply(colnames(characteristics), 
                                        function(key) { 
                                          if (is.character(characteristics[[key]]) || is.factor(characteristics[[key]]))
                                            paste(key, unique(characteristics[[key]]), sep = '/')
                                          else
                                            key }))
  
  all(sapply(constraints, 
             function(constraint) {  length(constraint$name) == 1 && constraint$name %in% characteristic_names }))
}
#' Check operator in constraints
#'
#' Check whether the operator elements within the constraints elements are correct.
#' 
#' @param constraints List of constraints.
#' @return \code{TRUE} if the operator elements in constraints are correct, \code{FALSE} otherwise.
constraints_correct_operators <- function(constraints) {
  if (is.null(constraints) || !constraints_correct_structure(constraints))
    return(TRUE)
  all(sapply(constraints, 
             function(constraint) { length(constraint$op) == 1 && constraint$op %in% c("<", "=", ">", "><") }))
}


#' Check target in constraints
#'
#' Check whether the target elements within the constraints elements are correct.
#' 
#' @param constraints List of constraints.
#' @return \code{TRUE} if the target elements in constraints are correct, \code{FALSE} otherwise.
constraints_correct_targets <- function(constraints) {
  if (is.null(constraints) || !constraints_correct_structure(constraints))
    return(TRUE)
  all(sapply(constraints, 
             function(constraint) { 
               is.numeric(constraint$target) && is.vector(constraint$target) && length(constraint$target) == length(strsplit(constraint$op, split = NULL)[[1]]) }))
}
