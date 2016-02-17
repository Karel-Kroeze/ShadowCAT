#' Latent trait estimation
#' 
#' maximum likelihood, maximum a posteriori and expected a posteriori estimates.
#' 
#' Obtains a latent trait estimate and variance of the estimate.
#' 
#' @section Maximum Likelihood and Maximum A Posteriori:
#' Maximum Likelihood and Maximum A-Posteriori estimates are based on a Newton-type non-linear minimization algorithm,
#' and handled with package \code{\link{nlm}}.
#'  
#' @section Expected A Posteriori:
#' Expected A-Posteriori estimates require the repeated evaluation of Q nested integrals, where Q is the dimensionality of the test.
#' This is performed with adaptive multidimensional Gauss-Hermite quadrature, and handled by package MultiGHQuad, see the documentation there for further details.
#' Note that the number of quadrature points used rises exponentially with the dimensionality of the test - use of EAP estimates with 
#' a 3+ dimensional test may not be a good idea.
#' 
#' @section Weighted Maximum Likelihood:
#' TODO: UPDATE WITH REFERENCES - MORE PRECISE DETAILS.
#' Note that WML estimation is not included. There is no satisfying solution to multidimensional Weighted Maximum Likelihood Estimation,
#' current WML estimators as used in other sources do not account for the covariance between dimensions. 
#' 
#' @section Variance:
#' Covariance matrix of the estimate is added to the estimate as an attribute.
#' 
#' @examples 
#' number_dimensions <- 1
#' estimate <- rep(.3, number_dimensions)
#' model <- "3PLM"
#' number_items <- 50
#' responses <- rep(c(1, 0), 17)
#' administered <- c(6:20, 31:49)
#' alpha <- matrix(runif(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
#' beta <- matrix(rnorm(number_items), nrow = number_items, ncol = 1)
#' guessing <- c(rep(.1, number_items / 2), rep(.2, number_items / 2))
#' number_itemsteps_per_item <- number_non_missing_cells_per_row(beta)
#' lower_bound <- rep(-3, number_dimensions)
#' upper_bound <- rep(3, number_dimensions)
#' prior <- diag(1)
#' prior_var_safe_nlm <- diag(number_dimensions)
#'
#' # obtain estimates
#' estimator <- "maximum_likelihood"
#' ML <- estimate_latent_trait(estimate, responses, prior, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, prior_var_safe_nlm)
#' estimator <- "maximum_aposteriori"
#' MAP <- estimate_latent_trait(estimate, responses, prior, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, prior_var_safe_nlm)
#' estimator <- "expected_aposteriori"
#' EAP <- estimate_latent_trait(estimate, responses, prior, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, prior_var_safe_nlm)
#' ML; MAP; EAP
#' 
#' # access variance
#' attr(ML, "variance")
#' 
#' # Note that expected_aposteriori takes considerably more time when dimensionality is higher...
#' number_dimensions <- 5
#' estimate <- rep(.3, number_dimensions)
#' alpha <- matrix(runif(number_items * number_dimensions, .3, 1.5), nrow = number_items, ncol = number_dimensions)
#' lower_bound <- rep(-3, number_dimensions)
#' upper_bound <- rep(3, number_dimensions)
#' prior <- diag(number_dimensions) 
#' 
#' estimator <- "maximum_aposteriori"
#' system.time(estimate_latent_trait(estimate, responses, prior, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, prior_var_safe_nlm = diag(number_dimensions)))
#' estimator <- "expected_aposteriori"
#' system.time(estimate_latent_trait(estimate, responses, prior, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound))
#' 
#' @param estimate vector containing theta estimate, with covariance matrix as an attribute
#' @param responses vector with person responses
#' @param prior prior covariance matrix for theta
#' @param model string, one of '3PLM', 'GPCM', 'SM' or 'GRM', for the three-parameter logistic, generalized partial credit, sequential or graded response model respectively
#' @param administered vector with indeces of administered items
#' @param number_dimensions number of dimensions
#' @param estimator type of estimator to be used, one of "maximum_aposteriori", "maximum_likelihood", or "expected_aposteriori"
#' @param alpha matrix of alpha paramteres
#' @param beta matrix of beta paramteres
#' @param guessing matrix of guessing parameters
#' @param number_itemsteps_per_item vector containing the number of non missing cells per row of the beta matrix
#' @param lower_bound vector with lower bounds for theta per dimension; estimated theta values smaller than the lowerbound values are truncated to the lowerbound values
#' @param upper_bound vector with upper bounds for theta per dimension; estimated theta values larger than the upperbound values are truncated to the upperbound values
#' @param prior_var_safe_nlm if not NULL, expected a posteriori estimate with prior variance(s) equal to prior_var_safe_ml is computed instead of maximum_likelihood/maximum_aposteriori, if maximum_likelihood/maximum_aposteriori estimate fails. Can be a scalar 
#' (if variance for each dimension is equal) or vector
#' @return vector containing the updated estimate with the covariance matrix as attribute
#' @importFrom MultiGHQuad init.quad eval.quad
#' @export
estimate_latent_trait <- function(estimate, responses, prior, model, administered, number_dimensions, estimator, alpha, beta, guessing, number_itemsteps_per_item, lower_bound, upper_bound, prior_var_safe_nlm = NULL) {
  result <- function() {
    updated_estimate <- get_updated_estimate_and_variance_attribute(estimator)
    trim_estimate(updated_estimate)
  }
  
  get_updated_estimate_and_variance_ml <- function() {
    # for now, simple nlm (TODO: look at optim, and possible reintroducing pure N-R).
    # We want a maximum, but nlm produces minima -> reverse function call.
    estimate <- tryCatch(nlm(f = likelihood_or_post_density, p = estimate, responses, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior, inverse_likelihood_or_post_density = TRUE)$estimate,
                         error = function(e) { switch_to_eap_if_requested() },
                         warning = function(w) { switch_to_eap_if_requested() })
    # TODO: We should really store info somewhere so we don't have to redo this (when using get_fisher_information based selection criteria).
    fisher_information_items <- get_fisher_information(estimate, model, number_dimensions, alpha, beta, guessing, number_itemsteps_per_item)
    fisher_information_test_so_far <- apply(fisher_information_items[,,administered, drop = FALSE], c(1, 2), sum)
    
    attr(estimate, "variance") <- solve(fisher_information_test_so_far)
    estimate
  }
  
  get_updated_estimate_and_variance_map <- function() {
    estimate <- tryCatch(nlm(f = likelihood_or_post_density, p = estimate, responses, model, administered, number_dimensions, estimator, alpha, beta, guessing, prior, inverse_likelihood_or_post_density = TRUE)$estimate,
                         error = function(e) { switch_to_eap_if_requested() },
                         warning = function(w) { switch_to_eap_if_requested() })
    fisher_information_items <- get_fisher_information(estimate, model, number_dimensions, alpha, beta, guessing, number_itemsteps_per_item)
    fisher_information_test_so_far <- apply(fisher_information_items[,,administered, drop = FALSE], c(1, 2), sum) +
      solve(prior)
    attr(estimate, "variance") <- solve(fisher_information_test_so_far)
    estimate
  }
  
  get_updated_estimate_and_variance_eap <- function(prior) {
    # Multidimensional Gauss-Hermite Quadrature
    # TODO: prior mean is currently fixed at zero, update when/if possible.
    # TODO: allow setting ip through internals argument(s)
    adapt <- if (length(responses) > 5 & !is.null(attr(estimate, 'variance'))) list(mu = estimate, Sigma = as.matrix(attr(estimate, "variance")))
    Q_dim_grid_quad_points <- init_quad(Q = number_dimensions,
                                        prior = list(mu = rep(0, number_dimensions), Sigma = prior),
                                        adapt = adapt,
                                        ip = switch(number_dimensions, 50, 15, 6, 4, 3))
    eval.quad(FUN = likelihood_or_post_density, X = Q_dim_grid_quad_points, responses, model, administered, number_dimensions, estimator = "maximum_likelihood", alpha, beta, guessing)
  }
  
  get_updated_estimate_and_variance_attribute <- function(estimator) {
    switch(estimator,
           maximum_likelihood = get_updated_estimate_and_variance_ml(),
           maximum_aposteriori = get_updated_estimate_and_variance_map(),
           expected_aposteriori = get_updated_estimate_and_variance_eap(prior))
  }
  
  trim_estimate <- function(estimate) {
    estimate[which(estimate > upper_bound)] <- upper_bound[which(estimate > upper_bound)]
    estimate[which(estimate < lower_bound)] <- lower_bound[which(estimate < lower_bound)]
    estimate
  }
  
  switch_to_eap_if_requested <- function() {
    if (is.null(prior_var_safe_nlm))
      stop("something went wrong with nlm maximization and prior_var_safe_ml was set to NULL")
    get_updated_estimate_and_variance_eap(prior = diag(number_dimensions) * prior_var_safe_nlm)
  }
  
  result()
}

#' Q-dimensional grid of quadrature points
#' 
#' Creates a flattened, rotated grid that incorporates correlation through an eigenvalue decomposition of the covariance matrix.
#' Copied from MultiGHQuad package and bugs fixed. Should be moved back to MultiGHQuad package
#' 
#' @param Q Number of dimensions
#' @param prior List of prior mean and covariance matrix
#' @param adapt List of adaptive mean and covariance matrix; if NULL no adaptation is used
#' @param ip Number of quadrature points per dimension. Defaults to 6. Note that the total number of quadrature points is ip^Q
#' @return A list with a matrix X of ip^Q by Q quadrature points and a vector W of length ip^Q associated weights
#' @examples grid_points_and_weights1 <- init_quad(Q = 1, prior = list(mu = rep(0, 1), Sigma = diag(1)*2), adapt = list(mu = rep(1, 1), Sigma = diag(1)*5), ip = 50)
#' round(grid_points_and_weights1$X[c(2, 8, 50),], 3) == c(-25.951, -17.382,  30.037) || stop("wrong")
#' round(grid_points_and_weights1$W[c(2, 8, 50)], 3) == c(-169.013, -76.610, -225.947) || stop("wrong")
#' grid_points_and_weights2 <- init_quad(Q = 2, prior = list(mu = rep(0, 2), Sigma = diag(2)*2), adapt = list(mu = rep(1, 2), Sigma = (matrix(rep(2, 4), ncol = 2) + diag(2)*5)), ip = 20)
#' round(grid_points_and_weights2$X[c(2, 8, 50),2], 3) == c(-24.858, -14.749 , -8.557) || stop("wrong")
#' round(grid_points_and_weights2$W[c(2, 8, 50)], 3) == c(-155.016, -76.948, -40.058) || stop("wrong")
#' @export
init_quad <- function (Q, prior = list(mu = rep(0, Q), Sigma = diag(Q)), 
                            adapt = NULL, ip = 6) { 
  if (!is.null(adapt) && !is.null(attr(adapt, "variance"))) {
    adapt <- list(mu = adapt, Sigma = attr(adapt, "variance"))
  }
  if (!is.null(adapt) && (!is.list(adapt) || length(adapt$mu) != Q || dim(adapt$Sigma) != c(Q, Q))) 
    stop("Size or format of Adapt argument invalid.")
  x <- fastGHQuad::gaussHermiteData(ip)
  w <- x$w/sqrt(pi)
  x <- x$x * sqrt(2)
  X <- as.matrix(expand.grid(lapply(apply(replicate(Q, x), 
                                          2, list), unlist)))
  trans <- function(X, Sigma) {
    lambda <- with(eigen(Sigma), {
      if (any(values < 0)) 
        warning("Matrix is not positive definite.")
      if (length(values) > 1) 
        vectors %*% diag(sqrt(values))
      else vectors * sqrt(values)
    })
    t(lambda %*% t(X))
  }
  
  g <- as.matrix(expand.grid(lapply(apply(replicate(Q, w), 
                                          2, list), unlist)))
  W <- apply(log(g), 1, sum)
  
  if (is.null(adapt)) {
    X <- trans(X, prior$Sigma)
    X <- t(t(X) + prior$mu)
  }
  else {
    X <- trans(X, adapt$Sigma)
    X <- t(t(X) + adapt$mu)
    adapt$chol <- chol(adapt$Sigma)
    adapt$det <- sum(log(diag(adapt$chol)))
    adapt$aux <- diag((X -  rep(1, ip^Q) %*% t(adapt$mu)) %*% solve(adapt$Sigma) %*% t(X -  rep(1, ip^Q) %*% t(adapt$mu)))
    prior$chol <- chol(prior$Sigma)
    prior$det <- sum(log(diag(prior$chol)))
    prior$aux <- diag((X -  rep(1, ip^Q) %*% t(prior$mu)) %*% solve(prior$Sigma) %*% t(X -  rep(1, ip^Q) %*% t(prior$mu)))
    fact <- (adapt$aux - prior$aux) / 2 + adapt$det - prior$det
    W <- W + fact
  }
  
  list(X = X, W = W)
}
