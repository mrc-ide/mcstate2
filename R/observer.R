##' Create an observer to extract additional details from your model
##' during the sampling process.
##'
##' Sometimes you want to extract additional information from your
##' model as your chain runs.  The case we see this most is when
##' running MCMC with a particle filter (pmcmc); in this case while
##' the likelihood calculation is running we are computing lots of
##' interesting quantities such as the final state of the system
##' (required for onward simulation) and filtered trajectories through
##' time.  Because these are stochastic we can't even just rerun the
##' model with our sampled parameter sets, because the final states
##' that are recovered depend also on the random number generators
##' (practically we would not want to, as it is quite expensive to
##' compute these quantities).
##'
##' The observer mechanism allows you to carry out arbitrary
##' additional calculations with your model at the end of the step.
##'
##' @title Create observer
##'
##' @param observe A function that will run with arguments `model`
##'   (the model that you passed in to [mcstate_model]) and `rng` (an
##'   rng object).  This function should return a list.  It is best if
##'   the list returned is named, with no duplicated names, and with
##'   return values that have the same exact dimensions for every
##'   iteration.  If you do this, then you will not have to provide
##'   any of the following arguments, which are going to be hard to
##'   describe and worse to implement.
##'
##' @param finalise A function that runs after a single chain has run,
##'   and you use to simplify across all samples drawn from that
##'   chain.  Takes a single argument which is the list with one set
##'   of observations per sample.
##'
##' @param combine A function that runs after all chains have run, and
##'   you use to simplify across chains.  Takes a single argument,
##'   which is the list with one set of observations per chain.
##'
##' @param append A function that runs after a continuation of chain
##'   has run (via [mcstate_sample_continue].  Takes two arguments
##'   representing the fully simplified observations from the first
##'   and second chains.
##'
##' @return An object with class `mcstate_observer` which can be
##'   passed in to `mcstate_sample`.
##'
##' @export
mcstate_observer <- function(observe,
                             finalise = NULL,
                             combine = NULL,
                             append = NULL) {
  if (is.null(finalise)) {
    finalise <- observer_finalise_auto
  }
  if (is.null(combine)) {
    combine <- observer_finalise_auto # same impl by default
  }
  if (is.null(append)) {
    append <- observer_append_auto
  }
  ret <- list(observe = observe,
              finalise = finalise,
              combine = combine,
              append = append)
  class(ret) <- "mcstate_observer"
  ret
}


observer_finalise_auto <- function(observations) {
  nms <- names(observations[[1]])
  ineligible <- is.null(nms) || anyDuplicated(nms) ||
    !all_same(lapply(observations, names))
  if (ineligible) {
    return(observations)
  }

  ret <- set_names(vector("list", length(nms)), nms)
  for (nm in nms) {
    el <- lapply(observations, "[[", nm)
    if (all_same(lapply(el, dim2))) {
      ret[[nm]] <- array_bind(arrays = el, after = length(dim2(el[[1]])))
    } else {
      ret[[nm]] <- el
    }
  }

  ret
}


## This makes me think that we should just never return things that we
## can't aggregate this way; I wonder if we can do better with the
## weird outputs above using lists?
observer_append_auto <- function(prev, curr) {
  nms <- names(prev)
  ineligible <- is.null(nms) || anyDuplicated(nms) ||
    !identical(nms, names(curr))
  if (ineligible) {
    ## Not really sure what we should do here, it's possible that we
    ## can at least assemble the lists well.
    return(list(prev, curr))
  }

  ret <- set_names(vector("list", length(nms)), nms)
  for (nm in nms) {
    el_prev <- prev[[nm]]
    el_curr <- curr[[nm]]
    bind_on <- length(dim2(el_prev)) - 1L
    ok <- length(dim2(el_prev)) == length(dim2(el_curr)) &&
      identical(dim2(el_prev)[-bind_on], dim(el_curr)[-bind_on])
    if (ok) {
      ret[[nm]] <- array_bind(el_prev, el_curr, on = bind_on)
    } else {
      ret[[nm]] <- list(el_prev, el_curr)
    }
  }

  ret
}
