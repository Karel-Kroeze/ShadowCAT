% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shadowcat.R
\name{shadowcat}
\alias{shadowcat}
\title{Get new item key and update estimate}
\usage{
shadowcat(answers, estimate, variance, model, alpha, beta, start_items,
  stop_test, estimator, information_summary, prior_form = NULL,
  prior_parameters = NULL, guessing = NULL, eta = NULL,
  constraints_and_characts = NULL, lower_bound = NULL, upper_bound = NULL,
  safe_eap = FALSE, eap_estimation_procedure = "riemannsum")
}
\arguments{
\item{answers}{Named list of previous answers and new answer, with names being the item keys. Should be initialized with \code{NULL}.
Answer options should start at 0.}

\item{estimate}{Vector with current estimate of latent trait theta. Length should be equal to the number of dimensions.}

\item{variance}{Current covariance matrix of the estimate, as vector.}

\item{model}{One of \code{"3PLM"}, \code{"GPCM"}, \code{"SM"} or \code{"GRM"}, for the three-parameter logistic, generalized partial credit, sequential or graded response model, respectively.}

\item{alpha}{Matrix of alpha parameters, one column per dimension, one row per item. Row names should contain the item keys. 
Note that so called within-dimensional models still use an alpha matrix, they simply have only one non-zero loading per item.}

\item{beta}{Matrix of beta parameters, one column per item step, one row per item. Row names should contain the item keys. 
Note that \code{shadowcat} expects answer categories to be sequential, and without gaps. That is, the weight parameter in the GPCM model is assumed to be sequential, 
and equal to the position of the 'location' of the beta parameter in the beta matrix.
The matrix should have a number of columns equal to the largest number of item steps over items, items with fewer answer categories should be 
right-padded with \code{NA}. \code{NA} values between answer categories are not allowed, and will lead to errors.}

\item{start_items}{List indicating the items that should be shown to the respondent before the theta estimate will be updated
for the first time. One of
\code{list(type = "random", n = ...)},
\code{list(type = "fixed", item_keys = ..., n = ...)}, or
\code{list(type = "random_by_dimension", n_by_dimension = ..., n = ...)},
where \code{n} is the total number of burn in items, \code{item_keys} is a character vector with keys of the burn in items, 
and \code{n_by_dimension} is the number of burn in items per dimension, or a vector with the number of burn in items for each dimension.
If \code{n} is 0, only \code{n} needs to be defined.
Note that the type \code{"random_by_dimension"} assumes that items load on a single dimension; if any item has a non-zero loading on a dimension, it is considered to be part of that dimension.}

\item{stop_test}{List indicating rules for when to terminate the test. Should be a list of the form
\code{list(target = ..., max_n = ..., min_n = ..., cutoffs = ...)}, 
where \code{target} is a vector indicating the maximum acceptable variance per dimension; \code{NULL} means no variance target,
\code{max_n} is the test length at which the test should be terminated (even if the target has not been reached yet), 
\code{min_n} is the minimum test length; \code{NULL} means no mimimum test length, and
\code{cutoffs} is a matrix containing cut off values per dimension (columns) and test iteration (rows). First row contains cut off values for when no items have been
administered yet, second row for when one item has been administered, etc. If estimate + 3SE < cutoff for each dimension at a certain iteration, test terminates; 
\code{NULL} means no cut off values.}

\item{estimator}{Type of estimator to be used, one of \code{"maximum_likelihood"}, \code{"maximum_aposteriori"}, or \code{"expected_aposteriori"}; see \code{details}.}

\item{information_summary}{How to summarize Fisher information, used for item selection. One of
\code{"determinant"}, \code{"posterior_determinant"}, \code{"trace"}, \code{"posterior_trace"}, or \code{"posterior_expected_kullback_leibler"}.
Fisher Information of the test so far (including all administered items) is added to the Fsher Information of the available
item before the summary is computed.}

\item{prior_form}{String indicating the form of the prior; one of \code{"normal"} or \code{"uniform"}. Not required if estimator is maximum likelihood.}

\item{prior_parameters}{List containing mu and Sigma of the normal prior: \code{list(mu = ..., Sigma = ...)}, or 
the upper and lower bound of the uniform prior: \code{list(lower_bound = ..., upper_bound = ...)}. Not required if estimator is maximum likelihood.
The list element \code{Sigma} should always be in matrix form. List elements \code{mu}, \code{lower_bound}, and \code{upper_bound} should always be vectors.
The length of \code{mu}, \code{lower_bound}, and \code{upper_bound} should be equal to the number of dimensions.
For uniform prior in combination with expected aposteriori estimation, true theta should fall within 
\code{lower_bound} and \code{upper_bound} and be not too close to one of these bounds, in order to prevent errors. 
Setting the \code{shadowcat} argument \code{safe_eap} to \code{TRUE} ensures that the estimation switches to maximum aposteriori if the expected aposteriori estimate fails.}

\item{guessing}{Matrix with one column of guessing parameters per item. Row names should contain the item keys. Optionally used in 3PLM model, ignored for all others.}

\item{eta}{Matrix of location parameters, optionally used in GPCM model, ignored for all others. Row names should contain the item keys.
If eta is defined, the beta matrix will be derived from this eta matrix by computing the cumulative sums of the rows of eta; see
Glas and Dagohoy (2006).}

\item{constraints_and_characts}{List with constraints and characteristics for Shadow Testing; \code{NULL} means no constraints. See \code{details}.}

\item{lower_bound}{Vector with lower bounds for theta per dimension. Estimated theta values smaller than the lower bound values are truncated to the lower bound values.
Can only be defined when estimator is maximum likelihood. Setting bounds with maximum likelihood estimation is equivalent to
using maximum aposteriori estimation with a uniform prior.}

\item{upper_bound}{Vector with upper bounds for theta per dimension. Estimated theta values larger than the upper bound values are truncated to the upper bound values.
Can only be defined when estimator is maximum likelihood. Setting bounds with maximum likelihood estimation is equivalent to
using maximum aposteriori estimation with a uniform prior.}

\item{safe_eap}{Only relevant if estimator is expected aposteriori. 
Set to \code{TRUE} if estimator should switch to maximum aposteriori if the integration algorithm results in an error.
An error may occur if the prior is uniform, estimator is expected aposteriori, and the bounds of the prior do not exceed the true theta value, or are too close to it.}

\item{eap_estimation_procedure}{String indicating the estimation procedure if estimator is expected aposteriori and prior form is normal. One of \code{"riemannsum"} for integration via Riemannsum or
\code{"gauss_hermite_quad"} for integration via Gaussian Hermite Quadrature. If prior form is uniform, estimation procedure should always be \code{"riemannsum"}.}
}
\value{
List containing:
\item{key_new_item}{The key of the next item to be administered given the answers to previous items. 
Next item is the item containing the maximum information, taking constraints into account if specified (Shadow Testing).}
\item{continue_test}{\code{TRUE} if test should be continued, \code{FALSE} if test should be terminated.}
\item{estimate}{Vector containing the updated theta estimate.}
\item{variance}{Vector containing the updated covariance matrix of theta.}
\item{answers}{Named list containing the answers to the administered items.}
}
\description{
Get the key of the new item to administer and an update of the theta estimate, based on given answer set.
}
\details{
Maximum Likelihood and Maximum A-Posteriori estimates are computed using minimization algorithms
as performed by \code{\link{nlm}} and \code{\link{constrOptim}}. Expected A-Posteriori estimates require the 
repeated evaluation of Q nested integrals, where Q is the dimensionality of the test.
This is performed with an adaptive Riemannsum or multidimensional Gauss-Hermite quadrature, the latter 
handled by package \code{MultiGHQuad}, see the documentation there for further details.
Note that the number of grid points used increases strongly with the dimensionality of the test. Use of Expected A-Posteriori 
estimates with a 3+ dimensional test may not be a good idea. Note that WML estimation is not included. There is no satisfying solution to multidimensional 
Weighted Maximum Likelihood Estimation. Current WML estimators as used in other sources do not account for the covariance between dimensions.

The argument \code{constraints_and_characts} should be \code{NULL} (no constraints on item selection) or a list of characteristics and constraints (Shadow Testing; Van der Linden, 2000).
The list should consist of two elements, named \code{characteristics} and \code{constraints}.
\code{characteristics} should be specified as a data frame of characteristics. Each row indicates the characteristics of
one item. Each column indicates how all items score on a certain characteristic. Characteristics may be categorical or numeric. 
\code{constraints} should be specified as a list of constraints, each constraint is a list with three named values:
\describe{
\item{\code{name}}{The column name of the characteristic this constraint applies to. For categorical characteristics the level should be specified as \code{name/value},
where \code{name} is the column name of the characteristic and \code{value} is the specific level of the characteristic this constraint applies to.} 
\item{\code{op}}{The logical operator to be used. Valid options are \code{"<"}, \code{"="}, \code{">"} and \code{"><"}.}
\item{\code{target}}{The target value, numeric. For categorical characteristics, it indicates the number of items of the relevant characteristic that should be administered (\code{"="}), or
minimally (\code{">"}), maximally (\code{"<"}), or minimally and maximally (\code{"><"}; vector with two values required) administered. For numeric characteristics,
it indicates the minimum and/or maximum sum allowed over all administered items, e.g., maximum time allowed.}
}
}
\examples{
alpha_beta <- simulate_testbank(model = "GPCM", number_items = 100, 
                                number_dimensions = 3, number_itemsteps = 3)
model <- "GPCM"
start_items <- list(type = 'fixed', item_keys = c("item33", "item5", "item23"), n = 3)
stop_test <- list(min_n = 4, max_n = 30, target = c(.1, .1, .1))
estimator <- "maximum_aposteriori"
information_summary <- "posterior_determinant"
prior_form <- "normal"
prior_parameters <- list(mu = c(0, 0, 0), Sigma = diag(3))

# Initial call: get key of first item to adminster
call1 <- shadowcat(answers = NULL, estimate = c(0, 0, 0), variance = as.vector(diag(3) * 25), 
                   model = model, alpha = alpha_beta$alpha, beta = alpha_beta$beta, 
                   start_items = start_items, stop_test = stop_test, 
                   estimator = estimator, information_summary = information_summary,
                   prior_form = prior_form, prior_parameters = prior_parameters)
# Second to fourth call: number of start items is set to 3, so no update in theta estimate yet
call2 <- shadowcat(answers = list(item33 = 2), estimate = call1$estimate,
                   variance = call1$variance, model = model, alpha = alpha_beta$alpha, 
                   beta = alpha_beta$beta, start_items = start_items, 
                   stop_test = stop_test, estimator = estimator, 
                   information_summary = information_summary, 
                   prior_form = prior_form, prior_parameters = prior_parameters)
call3 <- shadowcat(answers = list(item33 = 2, item5 = 3), 
                   estimate = call2$estimate, variance = call2$variance, 
                   model = model, alpha = alpha_beta$alpha, beta = alpha_beta$beta, 
                   start_items = start_items, stop_test = stop_test,
                   estimator = estimator, information_summary = information_summary, 
                   prior_form = prior_form, prior_parameters = prior_parameters)
call4 <- shadowcat(answers = list(item33 = 2, item5 = 3, item23 = 3), 
                   estimate = call3$estimate, variance = call3$variance,
                   model = model, alpha = alpha_beta$alpha, beta = alpha_beta$beta, 
                   start_items = start_items, stop_test = stop_test,
                   estimator = estimator, information_summary = information_summary, 
                   prior_form = prior_form, prior_parameters = prior_parameters)
# Fifth call: first time theta estimate is updated
call5 <- shadowcat(answers = list(item33 = 2, item5 = 3, item23 = 3, item84 = 1), 
                   estimate = call4$estimate, variance = call4$variance, model = model, 
                   alpha = alpha_beta$alpha, beta = alpha_beta$beta, start_items = start_items, 
                   stop_test = stop_test, estimator = estimator, 
                   information_summary = information_summary, 
                   prior_form = prior_form, prior_parameters = prior_parameters)
# Sixth call: use the updated estimate and variance as the current values for estimate and variance
call6 <- shadowcat(answers = list(item33 = 2, item5 = 3, item23 = 3, item84 = 1, item36 = 2), 
                   estimate = call5$estimate, variance = call5$variance, model = model, 
                   alpha = alpha_beta$alpha, beta = alpha_beta$beta, start_items = start_items, 
                   stop_test = stop_test, estimator = estimator, 
                   information_summary = information_summary, 
                   prior_form = prior_form, prior_parameters = prior_parameters)

# With constraints (shadow testing)
characteristics <- data.frame(content = sample(c('algebra','physics','calculus'), 
                                               size = 100, replace = TRUE),
                              time = runif(100),
                              exclusive = sapply(1:100, 
                                                 function (x) { 
                                                   if ( x \%in\% sample(1:100, size = 4) ) 1 else 0 
                                                 } ))
constraints <- list(list(name = 'content/algebra',
                         op = '><',
                         target = c(5, 10)), # ensure number of algebra items is between 5 and 10
                    list(name = 'content/physics',
                         op = '><',
                         target = c(2, 5)), # ensure number of physics items is between 2 and 5
                    list(name = 'time',
                         op = '<',
                         target = 20), # Ensure total tests takes no longer than 20 minutes
                    list(name = 'exclusive',
                         op = '<',
                         target = 2)) # Ensure number of exclusive items equals 2
constraints_and_characteristics <- list(characteristics = characteristics,
                                        constraints = constraints)
  
shadowcat(answers = NULL, estimate = c(0, 0, 0), variance = as.vector(diag(3) * 25), 
          model = model, alpha = alpha_beta$alpha, beta = alpha_beta$beta, 
          start_items = start_items, stop_test = stop_test, estimator = estimator, 
          information_summary = information_summary, prior_form = prior_form, 
          prior_parameters = prior_parameters, 
          constraints_and_characts = constraints_and_characteristics)
}
\references{
\itemize{
\item Glas, C. A. W., & Dagohoy, A. V. T. (2006). A Person Fit Test For Irt Models For Polytomous Items. Psychometrika, 72(2), 159-180.
\item Van der Linden, W. J. (2000). Constrained adaptive testing with shadow tests. In W. J. van der Linden & C. A. W. Glas (Eds.), Computerized adaptive testing: Theory and practice (pp. 27-52). Dordrecht,
the Netherlands: Kluwer Academic Publishers. 
}
}

