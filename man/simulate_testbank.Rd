% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_functions.R
\name{simulate_testbank}
\alias{simulate_testbank}
\title{Simulate alpha and beta}
\usage{
simulate_testbank(model, number_items = 50, number_dimensions = 1,
  number_itemsteps = 4, items_load_one_dimension = FALSE,
  varying_number_item_steps = FALSE, alpha_bounds = c(0.3, 1.5))
}
\arguments{
\item{model}{Model for which the item bank should be simulated. One of \code{"3PLM"}, \code{"GPCM"}, \code{"SM"} or \code{"GRM"}, 
for the three-parameter logistic, generalized partial credit, sequential or graded response model respectively. 
If model is 3PLM, number of item steps is forced to 1.}

\item{number_items}{Number of items.}

\item{number_dimensions}{Number of dimensions of theta.}

\item{number_itemsteps}{Number of item steps (number of categories minus 1); forced to 1 if model is 3PLM.}

\item{items_load_one_dimension}{If \code{TRUE}, force items to load on one dimension each.}

\item{varying_number_item_steps}{If \code{TRUE}, some item steps are set to \code{NA}; in this case \code{number_itemsteps}
is the maximum number of itemsteps.}

\item{alpha_bounds}{Vector containing lower and upper bound, respectively, of the uniform distribution from which the alpha values are drawn.}
}
\value{
Named list containing simulated alpha and beta matrix.
}
\description{
Simulate quick and simple itembanks.
}
\examples{
simulate_testbank(model = "GPCM", number_items = 50, number_dimensions = 2, number_itemsteps = 3)
simulate_testbank(model = "GPCM", number_items = 50, number_dimensions = 3, number_itemsteps = 4, 
                  items_load_one_dimension = TRUE, varying_number_item_steps = TRUE)
}

