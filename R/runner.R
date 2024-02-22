##' Run MCMC chains in series (one after another).  This is the
##' simplest chain runner, and the default used by [mcstate_sample()].
##' It has nothing that can be configured (yet).
##'
##' @title Run MCMC chain in series
##'
##' @return A runner of class `mcstate_runner` that can be passed to
##'   [mcstate_sample()]
##'
##' @export
mcstate_runner_serial <- function() {
  run <- function(pars, model, sampler, n_steps, rng) {
    lapply(
      seq_along(rng),
      function(i) {
        mcstate_run_chain(pars[i, ], model, sampler, n_steps, rng[[i]])
      })
  }
  structure(list(run = run), class = "mcstate_runner")
}


##' Run MCMC chains in parallel (at the same time).  This runner uses
##' the `parallel` package to distribute your chains over a number of
##' worker processes on the same machine.  Compared with the "worker"
##' support in mcstate version 1 this is very simple and we'll improve
##' it over time.  In particular we do not report back any information
##' about progress while a chain is running on a worker or even across
##' chains.  There's also no support to warn you if your number of
##' chains do not neatly divide through by the number of workers.
##' Mostly this exists as a proof of concept for us to think about the
##' different interfaces.  Unless your chains are quite slow, the
##' parallel runner will be slower than the serial runner
##' ([mcstate_runner_serial]) due to the overhead cost of starting the
##' cluster.
##'
##' @title Run MCMC chain in parallel
##'
##' @param n_workers Number of workers to create a cluster from.  In a
##'   multi-user setting be careful not to set this to more cores than
##'   you are allowed to use.  You can use `parallel::detectCores()`
##'   to get an estimate of the number of cores you have on a single
##'   user system (but this is often an overestimate as it returns the
##'   number of logical cores, including those from "hyperthreading").
##'   Fewer cores than this will be used if you run fewer chains than
##'   you have workers.
##'
##' @return A runner of class `mcstate_runner` that can be passed to
##'   [mcstate_sample()]
##'
##' @export
mcstate_runner_parallel <- function(n_workers) {
  ## An alternative way of doing this would be to start the cluster
  ## immediately on creation of the runner object, which means that
  ## subsequent uses would reuse a cluster; that would probably be
  ## nicer, but slightly harder.  If we took that approach for a callr
  ## based backend later (similar to what mcstate1 does) then we'll
  ## get the advantage that the cluster startup happens asyncronously
  ## and may be ready by the time we actually pass any work onto it.

  run <- function(pars, model, sampler, n_steps, rng) {
    n_chains <- length(rng)
    cl <- parallel::makeCluster(min(n_chains, n_workers))
    on.exit(parallel::stopCluster(cl))
    ## There is a good chance the user will see warnings about
    ## "package X not being available", and also a chance that objects
    ## will not survive serialisation/deserialisation between
    ## processes here - but we'll face that on basically every
    ## possible parallel backend really.
    ##
    ## It would be nice to refactor the rng so that we are using
    ## something closer to the rng pointer object rather than the big
    ## heavy R6 object as that already has the logic in for doing
    ## state sync, but this can all be done transparently later.
    pars_list <- asplit(pars, MARGIN = 1)
    rng_state <- lapply(rng, function(r) r$state())

    res <- parallel::clusterMap(
      cl,
      mcstate_run_chain_parallel,
      pars = pars_list,
      rng = rng_state,
      MoreArgs = list(model = model, sampler = sampler, n_steps = n_steps))
  }
  structure(list(run = run), class = "mcstate_runner")
}


## Later we could return the mutated rng state and set it back into
## the sampler, if we needed to continue for any reason.
mcstate_run_chain_parallel <- function(pars, model, sampler, n_steps, rng) {
  rng <- mcstate_rng$new(rng)
  mcstate_run_chain(pars, model, sampler, n_steps, rng)
}


mcstate_run_chain <- function(pars, model, sampler, n_steps, rng) {
  r_rng_state <- get_r_rng_state()
  density <- model$density(pars)
  state <- list(pars = pars, density = density)
  sampler$initialise(state, model, rng)

  history_pars <- matrix(NA_real_, n_steps + 1, length(pars))
  history_pars[1, ] <- pars
  history_density <- rep(NA_real_, n_steps + 1)
  history_density[[1]] <- density

  for (i in seq_len(n_steps)) {
    state <- sampler$step(state, model, rng)
    history_pars[i + 1, ] <- state$pars
    history_density[[i + 1]] <- state$density
  }

  ## Pop the parameter names on last
  colnames(history_pars) <- model$parameters

  ## I'm not sure about the best name for this
  details <- sampler$finalise(state, model, rng)

  ## This list will hold things that we'll use internally but not
  ## surface to the user in the final object (or summarise them in
  ## some particular way with no guarantees about the format).  We
  ## might hold things like start and stop times here in future.
  internal <- list(used_r_rng = !identical(get_r_rng_state(), r_rng_state))

  list(pars = history_pars,
       density = history_density,
       details = details,
       internal = internal)
}