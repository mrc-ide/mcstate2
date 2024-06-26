% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample.R
\name{mcstate_sample}
\alias{mcstate_sample}
\title{Sample from a model}
\usage{
mcstate_sample(
  model,
  sampler,
  n_steps,
  initial = NULL,
  n_chains = 1L,
  runner = NULL,
  observer = NULL,
  restartable = FALSE
)
}
\arguments{
\item{model}{The model to sample from; this should be a
\code{mcstate_model} for now, but we might change this in future to
test to see if things match an interface rather than a
particular class attribute.}

\item{sampler}{A sampler to use.  These will be described later,
but we hope to make these reasonably easy to implement so that
we can try out different sampling ideas.  For now, the only
sampler implemented is \code{\link[=mcstate_sampler_random_walk]{mcstate_sampler_random_walk()}}.}

\item{n_steps}{The number of steps to run the sampler for.}

\item{initial}{Optionally, initial parameter values for the
sampling.  If not given, we sample from the model (or its prior).}

\item{n_chains}{Number of chains to run.  The default is to run a
single chain, but you will likely want to run more.}

\item{runner}{A runner for your chains.  The default option is to
run chains in series (via \link{mcstate_runner_serial}).  The only
other current option is \link{mcstate_runner_parallel} which uses the
\code{parallel} package to run chains in parallel.  If you only run
one chain then this argument is best left alone.}

\item{observer}{An observer, created via \link{mcstate_observer}, which
you can use to extract additional information from your model at
points included in the chain (for example, trajectories from a
dynamical model).}

\item{restartable}{Logical, indicating if the chains should be
restartable.  This will add additional data to the chains
object.}
}
\value{
A list of parameters and densities; we'll write tools for
dealing with this later.  Elements include:
\itemize{
\item \code{pars}: A matrix with as many columns as you have parameters, and
as many rows as the total number of samples taken across all
chains (\code{n_steps * n_chains})
\item \code{density}: A vector of model log densities, one per step (length
\code{n_steps * n_chains})
\item \code{initial}: A record of the initial conditions, a matrix with as
many rows as you have parameters and \code{n_chains} columns (this is
the same format as the matrix form of the \code{initial} input
parameter)
\item \code{details}: Additional details reported by the sampler; this will
be a list of length \code{n_chains} (or \code{NULL}) and the details
depend on the sampler.  This one is subject to change.
\item \code{chain}: An integer vector indicating the chain that the samples
came from (1, 2, ..., \code{n_chains})
}
}
\description{
Sample from a model.  Uses a Monte Carlo method (or possibly
something else in future) to generate samples from your
distribution.  This is going to change a lot in future, as we add
support for distributing over workers, and for things like
parallel reproducible streams of random numbers.  For now it just
runs a single chain as a proof of concept.
}
