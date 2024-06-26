% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sampler-random-walk.R
\name{mcstate_sampler_random_walk}
\alias{mcstate_sampler_random_walk}
\title{Random Walk Sampler}
\usage{
mcstate_sampler_random_walk(vcv = NULL)
}
\arguments{
\item{vcv}{A variance covariance matrix for the proposal.}
}
\value{
A \code{mcstate_sampler} object, which can be used with
\link{mcstate_sample}
}
\description{
Create a simple random walk sampler, which uses a symmetric
proposal to move around parameter space.  This sampler supports
sampling from models where the likelihood is only computable
randomly (e.g., for pmcmc).
}
