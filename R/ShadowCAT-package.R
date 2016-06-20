#' Package for computerized adaptive testing
#' 
#' ShadowCAT finds the best next item given a set of answers to 
#' already adminstered items, and updates the estimate of the latent trait theta after each new response.
#' Optionally with Shadow Testing Procedure (Van der Linden, 2000). Item bank parameters are assumed to be known.
#' 
#' @section Package options:
#' \describe{
#' \item{shadowcat}{\code{\link{shadowcat}} is the main function of the package. It finds the
#' key of the best next item given a set of answers to already adminstered items, and an update of the
#' latent trait estimate. It also returns an indicator of whether the test should be continued and a list of answers 
#' to the already administered items. The function requires a fixed set of alpha and beta parameters.}
#' \item{simulation}{Three functions useful for simulations are included: \code{\link{simulate_testbank}} for simulation of alpha and beta matrices,
#' \code{\link{simulate_answer}} for simulation of answers to specific questions from the item bank, and
#' \code{\link{test_shadowcat}} for simulation of a testing routine with \code{shadowcat}}
#' }
#' 
#' 
#' @references
#' \itemize{
#' \item Glas, C. A. W., & Dagohoy, A. V. T. (2006). A Person Fit Test For Irt Models For Polytomous Items. Psychometrika, 72(2), 159-180.
#' \item Van der Linden, W. J. (2000). Constrained adaptive testing with shadow tests. In W. J. van der Linden & C. A. W. Glas (Eds.), Computerized adaptive testing: Theory and practice (pp. 27-52). Dordrecht,
#' the Netherlands: Kluwer Academic Publishers. 
#' }
"_PACKAGE"
#> [1] "_PACKAGE"
