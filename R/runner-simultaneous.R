##' Run chains *simultaneously*.  This differs from
##' [mcstate_runner_parallel], which runs chains individually in
##' parallel by working with models that can evaluate multiple
##' densities at the same time.  There are situations where this might
##' be faster than running in parallel, but primarily this exists so
##' that we can see that samplers can work with multiple samples at
##' once.
##'
##' @title Run MCMC chains simultaneously
##'
##' @inheritParams mcstate_runner_serial
##'
##' @return A runner of class `mcstate_runner` that can be passed to
##'   [mcstate_sample()]
##'
##' @export
mcstate_runner_simultaneous <- function(progress = NULL) {
  validate_suitable <- function(model, observer) {
    if (!model$properties$allow_multiple_parameters) {
      cli::cli_abort(
        c("Using multiple parameter sets with model that does not support it",
          i = paste("The 'simultaneous' runner requires that we can",
                    "provide multiple parameters at once to your model",
                    "and get back a vector of densities, but your model",
                    "does not support this (or does not advertise that it",
                    "does), with the property 'allow_multiple_parameters'",
                    "set to FALSE")),
        call = NULL)
    }
    if (!is.null(observer)) {
      cli::cli_abort(
        "Can't yet use observers with 'mcstate_runner_simultaneous'",
        call = NULL)
    }
  }

  run <- function(pars, model, sampler, observer, n_steps, rng) {
    validate_suitable(model, observer)
    n_chains <- length(rng)
    pb <- progress_bar(n_chains, n_steps, progress, FALSE, environment())
    progress <- pb(seq_len(n_chains))
    rng_state <- lapply(rng, function(r) r$state())
    ## TODO: get the rng state back into 'rng' here, or (better) look
    ## at if we should just be using seed instead here perhaps?
    ## rng_state <- matrix(res$internal$state$rng, ncol = n_chains)
    ## for (i in seq_len(n_chains)) {
    ##   rng[[i]]$set_state(rng_state[, i]) # not supported!
    ## }
    mcstate_run_chains_simultaneous(pars, model, sampler, observer,
                                    n_steps, progress, rng_state)
  }

  continue <- function(state, model, sampler, observer, n_steps) {
    validate_suitable(model, observer)
    n_chains <- length(state)
    pb <- progress_bar(n_chains, n_steps, progress, FALSE, environment())
    progress <- pb(seq_len(n_chains))
    mcstate_continue_chains_simultaneous(state, model, sampler, observer,
                                         n_steps, progress)
  }

  structure(list(run = run, continue = continue),
            class = "mcstate_runner")
}


## There is a lot of overlap here with mcstate_run_chain, we might
## later seek to harmonise this.  However, this is really just
## intended for developing and debugging the multiple-sample versions
## of the samplers, so it does not need to be particularly lovely.
##
## * we don't allow observers at the moment
## * we take our nice array format, break it into lists, and then
##   reform into essentially the same array format (later).  This is
##   hard to avoid.
## * there's quite a lot of churn around rng state
mcstate_run_chains_simultaneous <- function(pars, model, sampler, observer,
                                            n_steps, progress, rng_state) {
  r_rng_state <- get_r_rng_state()
  n_chains <- length(rng_state)
  rng <- mcstate_rng$new(unlist(rng_state), n_chains)

  chain_state <- sampler$initialise(pars, model, observer, rng)

  mcstate_run_chains_simultaneous2(chain_state, model, sampler, observer,
                                   n_steps, progress, rng, r_rng_state)
}


mcstate_continue_chains_simultaneous <- function(state, model, sampler,
                                                 observer, n_steps, progress) {
  r_rng_state <- get_r_rng_state()
  n_chains <- length(state)
  n_pars <- length(model$parameters)

  ## NOTE this duplicates code in mcstate_run_chains_simultaneous; we
  ## could move this elsewhere if we change the interface for the
  ## sequential version too?
  rng_state <- unlist(lapply(state, "[[", "rng"))
  rng <- mcstate_rng$new(unlist(rng_state), n_chains)

  ## This is the inverse of restart_data really
  pars <- matrix(vapply(state, function(x) x$chain$pars, numeric(n_pars)),
                 n_pars, n_chains)
  density <- vnapply(state, function(x) x$chain$density)
  chain_state <- list(pars = pars, density = density, observation = NULL)

  sampler_state <- lapply(state, function(x) x$sampler)
  sampler$set_internal_state(sampler_state)

  if (model$properties$is_stochastic) {
    stop("need to work with stochastic models still...")
    model$rng_state$set(state$model_rng)
  }

  mcstate_run_chains_simultaneous2(chain_state, model, sampler, observer,
                                   n_steps, progress, rng, r_rng_state)
}


mcstate_run_chains_simultaneous2 <- function(chain_state, model, sampler,
                                             observer, n_steps, progress, rng,
                                             r_rng_state) {
  initial <- chain_state$pars
  n_pars <- length(model$parameters)
  n_chains <- length(chain_state$density)

  history_pars <- array(NA_real_, c(n_pars, n_steps, n_chains))
  history_density <- matrix(NA_real_, n_steps, n_chains)

  for (i in seq_len(n_steps)) {
    chain_state <- sampler$step(chain_state, model, observer, rng)
    history_pars[, i, ] <- chain_state$pars
    history_density[i, ] <- chain_state$density
    ## TODO: also allow observations here if enabled
    progress(i)
  }

  ## Pop the parameter names on last
  rownames(history_pars) <- model$parameters

  ## I'm not sure about the best name for this
  details <- sampler$finalise(chain_state, model, rng)

  ## This simplifies handling later; we might want to make a new
  ## version of asplit that does not leave stray attributes on later
  ## though?
  rng_state <- lapply(asplit(matrix(rng$state(), ncol = n_chains), 2),
                      as.vector)

  ## TODO: observation finalisation; this will be weird
  internal <- list(
    used_r_rng = !identical(get_r_rng_state(), r_rng_state),
    state = list(
      chain = chain_state,
      rng = rng_state,
      sampler = sampler$get_internal_state(),
      simultaneous = TRUE,
      model_rng = if (model$properties$is_stochastic) model$rng_state$get()))

  list(initial = initial,
       pars = history_pars,
       density = history_density,
       details = details,
       observations = NULL,
       internal = internal)
}
