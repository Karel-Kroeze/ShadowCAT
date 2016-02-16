#' Simulate alpha and beta matrices
#' 
#' Quick and simple itembanks for testing purposes.
#' @param model String, one of '3PLM', 'GPCM', 'SM' or 'GRM', for the three-parameter logistic, generalized partial credit, sequential or graded response model respectively.
#' @param number_items number of items. Can be increased if number_items / number_dimensions is not an integer
#' @param number_dimensions number of dimensions
#' @param number_itemsteps number of item steps (number of categories minus 1); forced to 1 if model is 3PLM
#' @param items_load_one_dimension if TRUE, force items to load on one dimension each
#' @param varying_number_item_steps if TRUE, some item steps are set to NA; in this case number_itemsteps is the maximum number of itemsteps
#' @return list containing simulated alpha and beta matrix
simulate_testbank <- function(model, number_items = 50, number_dimensions = 1, number_itemsteps = 4, items_load_one_dimension = FALSE, varying_number_item_steps = FALSE){
  result <- function() {
    number_items <- get_number_items()
    number_itemsteps <- get_number_itemsteps()
    alpha <- get_alpha(number_items, number_itemsteps)
    rownames(alpha) <- str_c("item", 1:number_items)
    beta <- get_beta(number_items, number_itemsteps)
    rownames(beta) <- str_c("item", 1:number_items)
    list(alpha = alpha, beta = beta)
  }
  
  get_number_itemsteps <- function() {
    if (model == "3PLM") 
      1
    else
      number_itemsteps
  }
  
  get_number_items <- function() {
    # make sure the number of items is divisible by the number of dimensions
    ceiling(number_items / number_dimensions) * number_dimensions
  }
  
  get_alpha <- function(number_items, number_itemsteps) {
    alpha <- matrix(runif(number_items * number_dimensions, .3, 1.5), number_items, number_dimensions)
    if (number_dimensions > 1 && items_load_one_dimension) {
      set <- number_items / number_dimensions
      for (dimension in 1:number_dimensions)
        alpha[((dimension - 1) * set + 1):(dimension * set), (1:number_dimensions)[-dimension]] <- 0
      alpha
    }
    else {
      alpha
    }
  }
  
  get_beta <- function(number_items, number_itemsteps) {
    beta_one_itemstep <- matrix(rnorm(number_items), number_items, 1)
    if (number_itemsteps == 1) {
      beta_one_itemstep
    }
    else {
      # spread polytomous items cats -2 to +2.
      spread <- seq(-2, 2, length.out = number_itemsteps)
      beta_multiple_itemsteps <- ( if (model == "GPCM") 
                                     row_cumsum(t(apply(beta_one_itemstep, 1, function(x) x + spread))) 
                                   else
                                     t(apply(beta_one_itemstep, 1, function(x) x + spread)) )
      if (varying_number_item_steps)
        induce_varying_number_item_steps(beta_multiple_itemsteps, number_items, number_itemsteps)
      else
        beta_multiple_itemsteps
    }  
  }
  
  induce_varying_number_item_steps <- function(beta_multiple_itemsteps, number_items, number_itemsteps) {
    beta_multiple_itemsteps[sample(1:number_items, ceiling(number_items / 10)), number_itemsteps] <- NA
    if (number_itemsteps > 2)
      beta_multiple_itemsteps[sample(1:number_items, ceiling(number_items / 10)), (number_itemsteps - 1):number_itemsteps] <- NA
    beta_multiple_itemsteps
  }
  
  result()
}
