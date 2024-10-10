##' Run MCMC chains in series (one after another).  This is the
##' simplest chain runner, and the default used by [monty_sample()].
##' It has nothing that can be configured (yet).
##'
##' @title Run MCMC chain in series
##'
##' @param progress Optional logical, indicating if we should print a
##'   progress bar while running.  If `NULL`, we use the value of the
##'   option `monty.progress` if set, otherwise we show the progress
##'   bar (as it is typically wanted).  The progress bar itself
##'   responds to cli's options; in particular
##'   `cli.progress_show_after` and `cli.progress_clear` will affect
##'   your experience.  Alternatively, you can provide a string
##'   indicating the progress bar type.  Options are `fancy`
##'   (equivalent to `TRUE`), `none` (equivalent to `FALSE`) and
##'   `simple` (a very simple text-mode progress indicator designed
##'   play nicely with logging; it does not use special codes to clear
##'   the line).
##'
##' @return A runner of class `monty_runner` that can be passed to
##'   [monty_sample()]
##'
##' @export
monty_runner_serial <- function(progress = NULL) {
  run <- function(pars, model, sampler, n_steps, rng) {
    n_chains <- length(rng)
    pb <- progress_bar(n_chains, n_steps, progress, show_overall = TRUE)
    lapply(
      seq_along(rng),
      function(i) {
        monty_run_chain(pars[, i], model, sampler, n_steps,
                        pb(i), rng[[i]])
      })
  }

  continue <- function(state, model, sampler, n_steps) {
    n_chains <- length(state)
    pb <- progress_bar(n_chains, n_steps, progress, show_overall = TRUE)
    lapply(
      seq_along(state),
      function(i) {
        monty_continue_chain(state[[i]], model, sampler, n_steps, pb(i))
      })
  }

  monty_runner("Serial",
               "monty_runner_serial",
               run,
               continue)
}


##' Run MCMC chains in parallel (at the same time).  This runner uses
##' the `parallel` package to distribute your chains over a number of
##' worker processes on the same machine.  Compared with
##' [monty_runner_callr] (Whch is similar to the "worker" support in
##' `mcstate` version 1), this is very simple.  In particular we do
##' not report back any information about progress while a chain is
##' running on a worker or even across chains.  There's also no
##' support to warn you if your number of chains do not neatly divide
##' through by the number of workers.  Mostly this exists as a proof
##' of concept for us to think about the different interfaces.  Unless
##' your chains are quite slow, the parallel runner will be slower
##' than the serial runner ([monty_runner_serial]) due to the overhead
##' cost of starting the cluster.
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
##' @return A runner of class `monty_runner` that can be passed to
##'   [monty_sample()]
##'
##' @export
monty_runner_parallel <- function(n_workers) {
  ## TODO: accept a cluster as an optional argument here, I think.
  ## There is an interface we can use easily.

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
    pars_list <- asplit(pars, MARGIN = 2)
    rng_state <- lapply(rng, function(r) r$state())

    args <- list(model = model,
                 sampler = sampler,
                 n_steps = n_steps)

    ## To debug issues in the parallel sampler, it's most efficient to
    ## replace this call with `Map` and drop the `cl` argument, then
    ## you get nice stack traces.
    parallel::clusterMap(
      cl,
      monty_run_chain_parallel,
      pars = pars_list,
      rng = rng_state,
      MoreArgs = args)
  }

  continue <- function(state, model, sampler, n_steps) {
    n_chains <- length(state)
    cl <- parallel::makeCluster(min(n_chains, n_workers))
    on.exit(parallel::stopCluster(cl))
    args <- list(model = model,
                 sampler = sampler,
                 n_steps = n_steps,
                 progress = function(i) NULL)
    parallel::clusterMap(
      cl,
      monty_continue_chain,
      state,
      MoreArgs = args)
  }

  monty_runner("Parallel",
               "monty_runner_parallel",
               run,
               continue)
}


monty_run_chain_parallel <- function(pars, model, sampler, n_steps, rng) {
  rng <- monty_rng$new(rng)
  progress <- function(i) NULL
  monty_run_chain(pars, model, sampler, n_steps, progress, rng)
}


monty_run_chain <- function(pars, model, sampler, n_steps,
                            progress, rng) {
  r_rng_state <- get_r_rng_state()
  chain_state <- sampler$initialise(pars, model, rng)

  if (!is.finite(chain_state$density)) {
    ## Ideally, we'd do slightly better than this; it might be worth
    ## validing this even earlier (in the parameter validation step)
    ## but that has its own issues (e.g., in the heirarchical sampler
    ## we can't extract out the submodel densities, and we might have
    ## issues with rng behaviour of stochastic models).  So this way
    ## around we'll fail the chain, which is sad if it's a chain that
    ## only starts running after several hours, and particularly in
    ## parallel we won't find out that anything (and therefore
    ## everything) has failed until everything has completed.
    ##
    ## Once we formalise the idea of restartable model runners though,
    ## this will come out in the wash I think as we can initialise all
    ## our chains, then "restart" them from the first point; that is
    ## hard to do with the random walk sampler as it's stateless so
    ## there's no real effect, but we'll readdress this once we get
    ## the HMC or adaptive samplers set up.
    cli::cli_abort("Chain does not have finite starting density")
  }

  monty_run_chain2(chain_state, model, sampler, n_steps, progress,
                   rng, r_rng_state)
}


monty_continue_chain <- function(state, model, sampler, n_steps,
                                 progress) {
  r_rng_state <- get_r_rng_state()
  rng <- monty_rng$new(seed = state$rng)
  sampler$set_internal_state(state$sampler)
  if (model$properties$is_stochastic) {
    model$rng_state$set(state$model_rng)
  }
  monty_run_chain2(state$chain, model, sampler, n_steps, progress,
                   rng, r_rng_state)
}


monty_run_chain2 <- function(chain_state, model, sampler, n_steps,
                             progress, rng, r_rng_state) {
  initial <- chain_state$pars
  n_pars <- length(model$parameters)
  has_observer <- model$properties$has_observer

  history_pars <- matrix(NA_real_, n_pars, n_steps)
  history_density <- rep(NA_real_, n_steps)
  history_observation <- if (has_observer) vector("list", n_steps) else NULL

  for (i in seq_len(n_steps)) {
    chain_state <- sampler$step(chain_state, model, rng)
    history_pars[, i] <- chain_state$pars
    history_density[[i]] <- chain_state$density
    if (has_observer && !is.null(chain_state$observation)) {
      history_observation[[i]] <- chain_state$observation
    }
    progress(i)
  }

  ## Pop the parameter names on last
  rownames(history_pars) <- model$parameters

  ## I'm not sure about the best name for this
  details <- sampler$finalise(chain_state, model, rng)

  if (has_observer) {
    history_observation <- model$observer$finalise(history_observation)
  }

  ## This list will hold things that we'll use internally but not
  ## surface to the user in the final object (or summarise them in
  ## some particular way with no guarantees about the format).  We
  ## might hold things like start and stop times here in future.
  internal <- list(
    used_r_rng = !identical(get_r_rng_state(), r_rng_state),
    state = list(
      chain = chain_state,
      rng = rng$state(),
      sampler = sampler$get_internal_state(),
      model_rng = if (model$properties$is_stochastic) model$rng_state$get()))

  list(initial = initial,
       pars = history_pars,
       density = history_density,
       details = details,
       observations = history_observation,
       internal = internal)
}


monty_runner <- function(name, help, run, continue) {
  ret <- list(name = name,
              help = help,
              run = run,
              continue = continue)
  class(ret) <- "monty_runner"
  ret
}


##' @export
print.monty_runner <- function(x, ...) {
  cli::cli_h1("<monty_runner: {x$name} ({x$help})>")
  cli::cli_alert_info("Use {.help monty_sample} to use this runner")
  cli::cli_alert_info("See {.help {x$help}} for more information")
  invisible(x)
}
