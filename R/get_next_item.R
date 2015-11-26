#' Next item
#' 
#' Control function to select the next item, takes account of starting / stopping conditions.
#' 
#' @param person person object
#' @param test test object
#' @return integer item index
#' @export
next_item <- function(person, test) {
  result <- function() {
    if (length(person$responses) < test$start$n)
      get_start_item(test$start$type)
    else
      best_item(person, test)
  }
  
  get_start_item <- function(start_type) {
    switch(start_type,
           "random" = get_start_item_random(),
           "fixed" = get_start_item_fixed(),
           "randomByDimension" = get_start_item_randomByDimension())
  }
  
  get_start_item_random <- function() {
    sample(person$available, 1)
  }
  
  get_start_item_fixed <- function() {
    test$start$indices[length(person$responses) + 1]
  }
  
  # picks nByDimension starting items per dimension (or n_i if nByDimension is a length Q vector), assumes between models, 
  # if any item has a non-zero loading on a dimension, it is considered to be part of that dimension. 
  # they CAN overlap, which may cause unwanted side effects, and in within models the result is identical to 'normal' random starting.
  get_start_item_randomByDimension <- function() {
    n_by_dimension_vector <- get_n_by_dimension_vector()
    design_matrix_item_loadings <- test$items$pars$alpha > 0
    
    # if enough items from first dimension are drawn, items from next dimension are drawn, etc.
    dimension <- find_dimension_to_draw_from(n_by_dimension_vector)
     
    # get items for this dimension, cull non-applicable, and select at random.
    items_with_positive_loading_on_dimension <- (1:test$items$K)[design_matrix_item_loadings[,dimension]]
    positive_loading_and_available <- intersect(items_with_positive_loading_on_dimension, person$available)
    sample(positive_loading_and_available, 1)
  }
  
  get_n_by_dimension_vector <- function() {
    if (length(test$start$nByDimension) == 1) 
      rep(test$start$nByDimension, test$items$Q)
    else
      test$start$nByDimension
  }
  
  find_dimension_to_draw_from <- function(n_by_dimension_vector) {
    sum(length(person$responses) >= cumsum(n_by_dimension_vector)) + 1
  }
  
  validate <- function() {
    if (is.null(person))
      return(add_error("person", "is missing"))
    if (is.null(test))
      return(add_error("test", "is missing"))
    if (test$start$type == "randomByDimension" && length(test$start$nByDimension) > 1 && test$start$n != sum(test$start$nByDimension))
      return(add_error("start", "contains inconsistent information. Total length of start phase and sum of length per dimension do not match (n != sum(nByDimension)"))
    if (test$start$type == "randomByDimension" && length(test$start$nByDimension) == 1 && test$start$n != sum(rep(test$start$nByDimension, test$items$Q)))
      return(add_error("start", "contains inconsistent information. Total length of start phase and sum of length per dimension do not match"))
  }
  
  invalid_result <- function() {
    list(errors = errors())
  }
  
  validate_and_run()
}
