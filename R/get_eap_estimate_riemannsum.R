#' EAP estimate with Riemannsum
#' 
#' Compute the expected aposteriori estimate and covariance matrix of the latent trait theta. 
#' Integration approximation occurs via a Riemannsumm, where grid points can be adapted to the location of the
#' posterior distribution.
#' 
#' @param prior_form String indicating the form of the prior; one of \code{"normal"} or \code{"uniform"}.
#' @param prior_parameters List containing mu and Sigma of the normal prior: \code{list(mu = ..., Sigma = ...)}, or 
#' the upper and lower bound of the uniform prior: \code{list(lower_bound = ..., upper_bound = ...)}.
#' The list element \code{Sigma} should always be in matrix form. List elements \code{mu}, \code{lower_bound}, and \code{upper_bound} should always be vectors.
#' The length of \code{mu}, \code{lower_bound}, and \code{upper_bound} should be equal to the number of dimensions.
#' For uniform prior, true theta should fall within \code{lower_bound} and \code{upper_bound} and be not too close to one of these bounds, in order to prevent errors. 
#' @param adapt List containing mu and Sigma for the adaptation of the grid points: list(mu = ..., Sigma = ...).
#' If \code{NULL}, adaptation with normal prior is based on the prior parameters, and no adaptation is made with uniform prior.
#' @param number_gridpoints Value indicating the number of grid points per dimension to use for the Riemannsum.
#' @return Expected aposteriori estimate of the latent trait theta, with its covariance matrix as an attribute.
#' @importFrom mvtnorm dmvnorm
#' @importFrom Matrix nearPD
get_eap_estimate_riemannsum <- function(prior_form, prior_parameters, adapt = NULL, number_gridpoints = 50, alpha, beta, guessing, answers, administered, number_dimensions, model) {
  alpha <- get_subset(alpha, administered)
  beta <- get_subset(beta, administered)
  guessing <- get_subset(guessing, administered)
  
  result <- function() {
    mid_grid_points <- get_mid_grid_points() 
    joint_distribution <- get_joint_distribution(mid_grid_points = mid_grid_points)
    sum_joint <- sum(joint_distribution)
    sum_joint_times_grid <- colSums(joint_distribution * mid_grid_points)
    eap_theta_estimate <- as.vector(sum_joint_times_grid / sum_joint)

    var_eap_theta_estimate <- matrix(0, number_dimensions, number_dimensions)
    for (i in 1:nrow(mid_grid_points)) {
      deviation <- mid_grid_points[i, ] - eap_theta_estimate
      var_eap_theta_estimate <- var_eap_theta_estimate + ( deviation %*% t(deviation) * joint_distribution[i] / sum_joint )
    }
    
    attr(eap_theta_estimate, "variance") <- var_eap_theta_estimate
    eap_theta_estimate
  }
  
  trans <- function(grid_points, Sigma) {
    sigma_positive_definite <- nearPD(Sigma)$mat
    eigen_sigma <- eigen(sigma_positive_definite)
    if (length(eigen_sigma$values) > 1)
      t((eigen_sigma$vectors %*% diag(sqrt(eigen_sigma$values))) %*% t(grid_points))
    else
      t((eigen_sigma$vectors * sqrt(eigen_sigma$values)) %*% t(grid_points))
  }
  
  get_list_gridpoints_uniform <- function() {
    lapply(1:number_dimensions,
           function(dim) {
             lower_bound <- prior_parameters$lower_bound[dim]
             upper_bound <- prior_parameters$upper_bound[dim]
             interval_length <- upper_bound - lower_bound
             seq(lower_bound + .5 * interval_length/number_gridpoints, upper_bound - .5 * interval_length/number_gridpoints, interval_length/number_gridpoints)
           })
  } 
  
  get_list_gridpoints_normal <- function() {
    lapply(1:number_dimensions,
           function(dim) {
             seq(-5 + 5/number_gridpoints, 5 - 5/number_gridpoints, 10/number_gridpoints)
           })
  }
  
  get_list_gridpoints <- function() {
    switch(prior_form,
           uniform = get_list_gridpoints_uniform(),
           normal = get_list_gridpoints_normal())
  }
  
  get_mu_sigma_for_transformation <- function() {
    if (is.null(adapt) && prior_form == "uniform")
      return(NULL)
    if (is.null(adapt) && prior_form == "normal") {
      list(mu = prior_parameters$mu, sigma = prior_parameters$Sigma)
    }
    else if (prior_form == "normal") {
      list(mu = adapt$mu, sigma = adapt$Sigma)
    }
    else if (prior_form == "uniform") {
      if (all(diag(adapt$Sigma) < 9))
        list(mu = adapt$mu, sigma = adapt$Sigma)
      else
        list(mu = adapt$mu, sigma = adapt$Sigma / (max(diag(adapt$Sigma)) / 9))
    }
  }
  
  transform_grid_points <- function(grid_points_untransformed) {
    if (is.null(adapt) && prior_form == "uniform") {
      unname(as.matrix(grid_points_untransformed))
    }
    else {
      mu_sigma <- get_mu_sigma_for_transformation()
      mid_grid_points_trans <- trans(grid_points_untransformed, mu_sigma$sigma)
      t(t(mid_grid_points_trans) + mu_sigma$mu)
    }
  }
  
  get_mid_grid_points <- function() {
    grid_points_untransformed <- expand.grid(get_list_gridpoints())
    grid_points_transformed <- transform_grid_points(grid_points_untransformed)
    if (prior_form == "uniform")
      remove_rows_outside_bounds(matrix_to_evaluate = grid_points_transformed, lower_bound = prior_parameters$lower_bound, upper_bound = prior_parameters$upper_bound)
    else
      grid_points_transformed
  }
  
  get_joint_distribution <- function(mid_grid_points) {
    mu <- if (prior_form == "normal")
            as.numeric(prior_parameters$mu)
           else
             mu <- numeric(0)
    sigma <- if (prior_form == "normal")
               as.vector(nearPD(prior_parameters$Sigma)$mat)
             else
               sigma <- numeric(0)
    model_number <- switch(model,
                           "3PLM" = 0,
                           "GPCM" = 1,
                           "SM" = 2,
                           "GRM" = 3)
           
    .Call("RjointDistribution", as.vector(t(mid_grid_points)), as.numeric(alpha), as.numeric(beta), as.numeric(guessing), as.numeric(answers), 
          nrow(mid_grid_points), as.integer(number_dimensions), length(answers), as.integer(prior_form == normal), mu, sigma, as.integer(model_number), package = "ShadowCAT")
  }  
    
  result()
}
