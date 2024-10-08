% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model.R
\name{monty_model_density}
\alias{monty_model_density}
\title{Compute log density}
\usage{
monty_model_density(model, parameters)
}
\arguments{
\item{model}{A \link{monty_model} object}

\item{parameters}{A vector or matrix of parameters}
}
\value{
A log-density value, or vector of log-density values
}
\description{
Compute log density for a model.  This is a wrapper around the
\verb{$density} property within a \link{monty_model} object.
}
\examples{
m <- monty_model_function(function(a, b) dnorm(0, a, b, log = TRUE))
monty_model_density(m, c(0, 1))
monty_model_density(m, c(0, 10))
}
\seealso{
\link{monty_model_gradient} for computing gradients and
\link{monty_model_direct_sample} for sampling from a model.
}
