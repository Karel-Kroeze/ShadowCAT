#' Latent trait estimation
#'
#' ML, MAP and EAP estimates.
#'
#' Obtains a latent trait estimate and variance of the estimate.
#'
#' @section ML and MAP:
#' Maximum Likelihood and Maximum A-Posteriori estimates are based on a Newton-type non-linear minimization algorithm,
#' and handled with package \code{\link{nlm}}.
#'
#' @section EAP:
#' Estimated A-Posteriori estimates require the repeated evaluation of Q nested integrals, where Q is the dimensionality of the test.
#' This is performed with adaptive multidimensional Gauss-Hermite quadrature, and handled by package \code{\link{MultiGHQuad}}, see the documentation there for further details.
#' Note that the number of quadrature points used rises exponentially with the dimensionality of the test - use of EAP estimates with
#' a 3+ dimensional test may not be a good idea.
#'
#' @section WML:
#' TODO: UPDATE WITH REFERENCES - MORE PRECISE DETAILS.
#' Note that WML estimation is not included. There is no satisfying solution to multidimensional Weighted Maximum Likelihood Estimation,
#' current WML estimators as used in other sources do not account for the covariance between dimensions.
#'
#' @section Variance:
#' Variance of the estimate is added to the estimate as an attribute.
#'
#' @examples
#' # create a basic test + person
#' items <- createTestBank("GPCM")
#' test <- initTest(items, estimator = "ML")
#' person <- initPerson(items)
#'
#' # answer a few items
#' person <- answer(person, test, sample(items$K, 10))
#'
#' # obtain estimates
#' ML <- estimate(person, test)$estimate
#' test$estimator <- "MAP"
#' MAP <- estimate(person, test)$estimate
#' test$estimator <- "EAP"
#' EAP <- estimate(person, test)$estimate
#' ML; MAP; EAP
#'
#' # access variance
#' attr(ML, "variance")
#'
#' # Note that EAP takes considerably more time when dimensionality is higher...
#' items5d <- createTestBank("GPCM", Q=5)
#' test5dEAP <- initTest(items, estimator = "EAP")
#' test5dMAP <- initTest(items, estimator = "MAP")
#' person5d <- answer(initPerson(items), test, sample(items5d$K, 10))
#'
#' system.time(estimate(person5d, test5dEAP))
#' system.time(estimate(person5d, test5dMAP))
#'
#' @param person Person object, see \code{\link{initPerson}}.
#' @param test Test object, see \code{\link{initTest}}.
#' @param ... Additional parameters passed on to the underlying estimation functions. 1
#' @return person Person object, amended with the new estimate.
#' @importFrom MultiGHQuad init.quad eval.quad
#' @export
estimate <- function(person, test, ...) {
  # catch additional arguments
  args <- list(...)
  
  if (test$estimator %in% c("ML", "MAP")) {
    # LL is the target function, test, person and minimize need to be passed on. We also want the value of the hessian at the final estimate.
    # note that prior is applied in LL if posterior = true
    # suppress warnings and errors and do EAP instead.
    
    person$estimate <- withCallingHandlers(
      tryCatch(
        nlm( f = LL, p = person$estimate, test = test, person = person, minimize = TRUE # We want a maximum, but nlm produces minima -> reverse function call.
        )$estimate,
        error = function(e) {
          message(paste0( test$estimator, " failed, trying EAP estimate. \n", e$message ) )
          test$estimator <- "EAP"
          return(estimate(person, test)$estimate)
        }
      ),
      warning = function(w) {
        message(w$message)
        invokeRestart("muffleWarning")
      }
    )
    
    # variance
    # get FI
    # TODO: We should really store info somewhere so we don't have to redo this (when using FI based selection criteria).
    info <- FI(test, person)
    
    # information so far is simply the sum of individual item information
    so_far <-
      apply(info[, , person$administered, drop = FALSE], c(1, 2), sum)
    
    #add prior
    if (test$estimator == "MAP")
      so_far <- so_far + solve(person$prior)
    
    # inverse
    attr(person$estimate, "variance") <- solve(so_far)
  }
  
  if (test$estimator == "EAP") {
    # Multidimensional Gauss-Hermite Quadrature
    Q <- test$items$Q
    # TODO: prior mean is currently fixed at zero, update when/if possible.
    # TODO: allow setting ip through internals argument(s)
    adapt <- NULL
    if (length(person$responses) > 5 &
        !is.null(attr(person$estimate, 'variance')))
      adapt <-
      list(mu = person$estimate, Sigma = as.matrix(attr(person$estimate, "variance")))
    
    QP <- init.quad(
      Q = Q,
      prior = list(mu = rep(0, test$items$Q), Sigma = person$prior),
      adapt = adapt,
      ip = switch(Q, 50, 15, 6, 4, 3)
    )
    person$estimate <-
      eval.quad(
        FUN = LL,
        X = QP,
        test = test,
        person = person,
        ...
      )
  }
  
  # enforce boundaries.
  # TODO: make debug output toggleable
  # if (any(person$estimate > test$upperBound | person$estimate < test$lowerBound)) cat("Estimate outside boundaries (k =", length(person$responses), "estimate =", paste0(round(person$estimate, 2), collapse = ", "), ").\n")
  person$estimate[which(person$estimate > test$upperBound)] <-
    test$upperBound[which(person$estimate > test$upperBound)]
  person$estimate[which(person$estimate < test$lowerBound)] <-
    test$lowerBound[which(person$estimate < test$lowerBound)]
  
  return(invisible(person))
}
