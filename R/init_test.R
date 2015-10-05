#' Initiate ShadowCAT test object
#' 
#' requires and returns information about adaptive test rules
#' 
#' @param items item characteristics as returned by init_item_bank()
#' @param start items that are shown to the patient before adaptive proces starts; one of
#' list(type = 'random', n)
#' list(type = 'fixed', indices, n)
#' list(type = 'randomByDimension', nByDimension, n)
#' where n = total number of initial items, indices = vector of initial item indices, 
#' nByDimension = scalar of number of initial items per dimension, or vector with number of initial items for each dimension,
#' @param stop rule for when to stop providing new items to patient; one of
#' list(type = 'length', n = 30); in this case, max_n is actually superfluous
#' list(type = 'variance', target)
#' where n is test length at which testing should stop,  target = vector of maximum acceptable variances per dimension
#' @param max_n maximum test length alowed
#' @param estimator type of estimator to be used, one of "MAP" (Maximum a posteriori estimation), "EAP" (Expected A Posteriori Estimation), or "ML" (maximum likelihood)
#' @param objective how to summarize information; one of
#' "D" = compute determinant(info_sofar_QxQ + info_QxQ_k) for each yet available item k
#  "PD" = compute determinant(info_sofar_QxQ_plus_prior + info_QxQ_k) for each yet available item k
#  "A" = compute trace((info_sofar_QxQ + info_QxQ_k) for each yet available item k
#  "PA" = compute trace(info_sofar_QxQ_plus_prior + info_QxQ_k) for each yet available item k
#  "PEKL" = compute Posterior expected Kullback-Leibler Information 
#' @param selection selection criterion; one of "MI" (maximum information) or "Shadow" (maximum information and take constraints into account)
#' @param constraints list with constraints and characteristics
#' constraints should be specified as a list of constraints, each constraint is a list with three named values;
#' \code{name} the column name of the characteristic this constraint applies to. For categorical characteristics the level should be specified as \code{name/value}.
#' \code{op} the logical operator to be used. Valid options are "<", "=", ">" and "><".
#' \code{target} the target value, numeric. If the operator is "><", this should be a length two vector in between which the target should fall.
#' characteristics should be a data.frame with characteristics, one row per item, one column per characteristic.
#' @param exposure vector, need to figure out what this is yet
#' @param lowerBound vector with lower bounds for theta per dimension 
#' @param upperBound vector with upper bounds for theta per dimension
#' @param ...
#' @return ShadowCAT.test object, containing list of input, where constraints are restructured 
#' @export
initTest <- function(items, 
                     start = list(type = 'random', n = 5), 
                     stop = list(type = 'length', n = 30),
                     max_n = 50, # utter maximum
                     estimator = 'MAP',
                     objective = 'PD',
                     selection = 'MI',
                     constraints = NULL,
                     exposure = NULL,
                     lowerBound = rep(-3, items$Q),
                     upperBound = rep(3, items$Q),
                     ...)
  {

  # attach everything
  out <- list(items = items,
              start = start,
              stop = stop,
              max_n = max_n,
              lowerBound = lowerBound,
              upperBound = upperBound,
              estimator = estimator,
              objective = objective,
              selection = selection,
              constraints = constraints,
              internal = list(...))
  
  attr(out, 'class') <- c("ShadowCAT.test")  
  
  # set up default constraints
  out$constraints <- createConstraints(out, constraints$characteristics, constraints$constraints)
  
  # maybe for the future, not finished yet:
  # apply ineligibility constraints (if requested).
  #if(!is.null(exposure))
  #  out <- with(exposure, createExposureConstraint(out, exposure, eligible, feasible, total, target))
  
  invisible(out)
}