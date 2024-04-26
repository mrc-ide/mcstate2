##' Sample from a model.  Uses a Monte Carlo method (or possibly
##' something else in future) to generate samples from your
##' distribution.  This is going to change a lot in future, as we add
##' support for distributing over workers, and for things like
##' parallel reproducible streams of random numbers.  For now it just
##' runs a single chain as a proof of concept.
##'
##' @title Sample from a model
##'
##' @param model The model to sample from; this should be a
##'   `mcstate_model` for now, but we might change this in future to
##'   test to see if things match an interface rather than a
##'   particular class attribute.
##'
##' @param sampler A sampler to use.  These will be described later,
##'   but we hope to make these reasonably easy to implement so that
##'   we can try out different sampling ideas.  For now, the only
##'   sampler implemented is [mcstate_sampler_random_walk()].
##'
##' @param n_steps The number of steps to run the sampler for.
##'
##' @param initial Optionally, initial parameter values for the
##'   sampling.  If not given, we sample from the model (or its prior).
##'
##' @param n_chains Number of chains to run.  The default is to run a
##'   single chain, but you will likely want to run more.
##'
##' @param runner A runner for your chains.  The default option is to
##'   run chains in series (via [mcstate_runner_serial]).  The only
##'   other current option is [mcstate_runner_parallel] which uses the
##'   `parallel` package to run chains in parallel.  If you only run
##'   one chain then this argument is best left alone.
##'
##' @param observer An observer, created via [mcstate_observer], which
##'   you can use to extract additional information from your model at
##'   points included in the chain (for example, trajectories from a
##'   dynamical model).
##'
##' @param restartable Logical, indicating if the chains should be
##'   restartable.  This will add additional data to the chains
##'   object.
##'
##' @return A list of parameters and densities; we'll write tools for
##'   dealing with this later.  Elements include:
##'
##' * `pars`: A matrix with as many columns as you have parameters, and
##'   as many rows as the total number of samples taken across all
##'   chains (`n_steps * n_chains`)
##'
##' * `density`: A vector of model log densities, one per step (length
##'   `n_steps * n_chains`)
##'
##' * `initial`: A record of the initial conditions, a matrix with as
##'   many rows as you have parameters and `n_chains` columns (this is
##'   the same format as the matrix form of the `initial` input
##'   parameter)
##'
##' * `details`: Additional details reported by the sampler; this will
##'   be a list of length `n_chains` (or `NULL`) and the details
##'   depend on the sampler.  This one is subject to change.
##'
##' * `chain`: An integer vector indicating the chain that the samples
##'   came from (1, 2, ..., `n_chains`)
##'
##' @export
mcstate_sample <- function(model, sampler, n_steps, initial = NULL,
                           n_chains = 1L, runner = NULL, observer = NULL,
                           restartable = FALSE) {
  if (!inherits(model, "mcstate_model")) {
    cli::cli_abort("Expected 'model' to be an 'mcstate_model'",
                   arg = "model")
  }
  if (!inherits(sampler, "mcstate_sampler")) {
    cli::cli_abort("Expected 'sampler' to be an 'mcstate_sampler'",
                   arg = "sampler")
  }
  if (is.null(runner)) {
    runner <- mcstate_runner_serial()
  }
  if (!is.null(observer) && !inherits(observer, "mcstate_observer")) {
    cli::cli_abort("Expected 'observer' to be an 'mcstate_observer'",
                   arg = "observer")
  }
  if (!inherits(runner, "mcstate_runner")) {
    cli::cli_abort("Expected 'runner' to be an 'mcstate_runner'",
                   arg = "runner")
  }

  rng <- initial_rng(n_chains)
  pars <- initial_parameters(initial, model, rng, environment())
  res <- runner$run(pars, model, sampler, observer, n_steps, rng)

  samples <- combine_chains(res, observer)
  if (restartable) {
    samples$restart <- restart_data(res, model, sampler, observer, runner)
  }
  samples
}


##' Continue (restart) chains started by [mcstate_sample].  Requires
##' that the original chains were run with `restartable = TRUE`.
##' Running chains this way will result in the final state being
##' exactly the same as running for the total (original + continued)
##' number of steps in a single push.
##'
##' @title Continue sampling
##'
##' @param samples A `mcstate_samples` object created by
##'   [mcstate_sample()]
##'
##' @param n_steps The number of new steps to run
##'
##' @param runner Optionally, a runner for your chains.  The default
##'   is to continue with the backend that you used to start the
##'   chains via [mcstate_sample] (or on the previous restart with
##'   this function).  You can use this argument to change the runner,
##'   which might be useful if transferring a pilot run from a
##'   high-resource environment to a lower-resource environment.  If
##'   given, must be an `mcstate_runner` object such as
##'   [mcstate_runner_serial] or [mcstate_runner_parallel].  You can
##'   use this argument to change the configuration of a runner, as
##'   well as the type of runner (e.g., changing the number of
##'   allocated cores).
##'
##' @inheritParams mcstate_sample
##'
##' @return A list of parameters and densities
##' @export
mcstate_sample_continue <- function(samples, n_steps, restartable = FALSE,
                                    runner = NULL) {
  if (!inherits(samples, "mcstate_samples")) {
    cli::cli_abort("Expected 'samples' to be an 'mcstate_samples' object")
  }
  if (is.null(samples$restart)) {
    cli::cli_abort(
      c("Your chains are not restartable",
        i = paste("To work with 'mcstate_sample_continue', you must",
                  "use the argument 'restartable = TRUE' when calling",
                  "mcstate_sample()")))
  }

  if (is.null(runner)) {
    runner <- samples$restart$runner
  } else {
    assert_is(runner, "mcstate_runner")
  }
  state <- samples$restart$state
  model <- samples$restart$model
  sampler <- samples$restart$sampler
  observer <- samples$restart$observer

  res <- runner$continue(state, model, sampler, observer, n_steps)

  samples <- append_chains(samples, combine_chains(res, observer), observer)

  if (restartable) {
    samples$restart <- restart_data(res, model, sampler, observer, runner)
  }
  samples
}


mcstate_sampler <- function(name, initialise, step, finalise,
                            get_internal_state, set_internal_state) {
  ret <- list(name = name,
              initialise = initialise,
              step = step,
              finalise = finalise,
              get_internal_state = get_internal_state,
              set_internal_state = set_internal_state)
  class(ret) <- "mcstate_sampler"
  ret
}


initial_parameters <- function(initial, model, rng, call = NULL) {
  n_pars <- length(model$parameters)
  n_chains <- length(rng)
  if (is.null(initial)) {
    require_direct_sample(model,
                          "'initial' must be provided with this model",
                          arg = "initial", call = environment())
    ## Really this would just be from the prior; we can't directly
    ## sample from the posterior!
    initial <- lapply(rng, function(r) direct_sample_within_domain(model, r))
  }
  if (is.list(initial)) {
    if (length(initial) != n_chains) {
      cli::cli_abort(
        c(paste("Unexpected length for list 'initial'",
                "(given {length(initial)}, expected {n_chains})"),
          i = paste("You are running {n_chains}, so I expected a list with",
                    "that many entries, each of which has {n_pars}",
                    "element{?s}")),
        arg = "initial", call = call)
    }
    err <- lengths(initial) != n_pars
    if (any(err)) {
      ## The format here is not beautiful but it's fairly informative
      ## at least.
      len <- as.character(lengths(initial)[err])
      len[[1]] <- paste(len[[1]], ngettext(len[[1]], "value", "values"))
      detail <- sprintf("%d (%s)", which(err), len)
      cli::cli_abort(
        c(paste("Unexpected initial parameter length; expected {n_pars}"),
          i = "Incorrect length for element{?s} {detail}",
          i = paste("Your model has {n_pars} parameter{?s}, so each",
                    "element of your list must have this many elements")),
        arg = "initial", call = call)
    }
    initial <- matrix(unlist(unname(initial)), n_pars, n_chains)
  } else if (is.matrix(initial)) {
    hint_matrix <- paste(
      "Expected a matrix with {n_pars} row{?s} (1 per parameter) and",
      "{n_chains} column{?s} (1 per chain)")
    if (nrow(initial) != n_pars) {
      cli::cli_abort(
        c(paste("Unexpected number of rows in 'initial'",
                "(given {nrow(initial)}, expected {n_pars})"),
          i = hint_matrix),
        arg = "initial", call = call)
    }
    if (ncol(initial) != n_chains) {
      cli::cli_abort(
        c(paste("Unexpected number of columns in 'initial'",
                "(given {ncol(initial)}, expected {n_chains})"),
          i = hint_matrix),
        arg = "initial", call = call)
    }
  } else {
    if (length(initial) != n_pars) {
      cli::cli_abort(
        c(paste("Unexpected length for vector 'initial'",
                "(given {length(initial)}, expected {n_pars})"),
          i = "I expected a vector with {n_pars} element{?s}, 1 per parameter"),
        arg = "initial", call = call)
    }
    initial <- matrix(initial, n_pars, n_chains)
  }

  err <- initial < model$domain[, 1] | initial > model$domain[, 2]
  if (any(err)) {
    ## Reporting on this is hard because we want to know what chain
    ## starting point is the issue, and which parameters.
    err_parameters <- model$parameters[rowSums(err) > 0]
    hint_parameters <- "Issues with parameter{?s}: {squote(err_parameters)}"
    err_chain <- colSums(err) > 0
    if (all(err_chain) && n_chains > 1) {
      hint_chain <- "Issues with every chain"
    } else {
      hint_chain <-
        "Issues with {cli::qty(sum(err_chain))}chain{?s} {which(err_chain)}"
    }
    cli::cli_abort(
      c("Initial conditions do not fall within parameter domain",
        x = hint_parameters,
        x = hint_chain),
      arg = "initial", call = call)
  }
  initial
}


combine_chains <- function(res, observer = NULL) {
  if (is.null(names(res))) {
    pars <- array_bind(arrays = lapply(res, "[[", "pars"), after = 2)
    density <- array_bind(arrays = lapply(res, "[[", "density"), after = 1)
    details <- lapply(res, "[[", "details")
    details <- if (all(vlapply(details, is.null))) NULL else details

    n_pars <- nrow(pars)
    initial <- vapply(res, "[[", numeric(n_pars), "initial")
    if (n_pars == 1) {
      initial <- rbind(initial, deparse.level = 0)
    }

    if (is.null(observer)) {
      observations <- NULL
    } else {
      observations <- observer$combine(lapply(res, "[[", "observations"))
    }
    used_r_rng <- vlapply(res, function(x) x$internal$used_r_rng)
  } else {
    stopifnot(is.null(observer)) # prevented earlier
    pars <- res$pars
    density <- res$density
    details <- res$details
    observations <- res$observations
    used_r_rng <- res$internal$used_r_rng
    initial <- res$initial
  }

  rownames(initial) <- rownames(pars)

  if (any(used_r_rng)) {
    cli::cli_warn(c(
      "Detected use of R's random number generators",
      i = paste("Your model has used R's random number generators (e.g.,",
                "via rnorm, runif, sample, etc).  This means that your",
                "results will not be reproducible as you change the sample",
                "runner")))
  }

  samples <- list(pars = pars,
                  density = density,
                  initial = initial,
                  details = details,
                  observations = observations)
  class(samples) <- "mcstate_samples"
  samples
}


append_chains <- function(prev, curr, observer = NULL) {
  if (is.null(observer)) {
    observations <- NULL
  } else {
    observations <- observer$append(prev$observations, curr$observations)
  }
  samples <- list(pars = array_bind(prev$pars, curr$pars, on = 2),
                  density = array_bind(prev$density, curr$density, on = 1),
                  initial = prev$initial,
                  details = curr$details,
                  observations = observations)
  class(samples) <- "mcstate_samples"
  samples
}


initial_rng <- function(n_chains, seed = NULL) {
  lapply(mcstate_rng_distributed_state(n_nodes = n_chains, seed = seed),
         function(s) mcstate_rng$new(seed = s))
}


restart_data <- function(res, model, sampler, observer, runner) {
  if (is.null(names(res))) {
    state <- lapply(res, function(x) x$internal$state)
  } else {
    if (!is.null(res$internal$state$sampler)) {
      stop("need to spread out sampler state")
    }
    if (!is.null(res$internal$state$model_rng)) {
      stop("need to spread out model state")
    }
    n_chains <- length(res$internal$state$chain$density)
    state <- lapply(seq_len(n_chains), function(i) {
      list(chain = list(
             pars = res$internal$state$chain$pars[, i],
             density = res$internal$state$chain$density[i],
             observation = NULL),
           rng = res$internal$state$rng[[i]],
           sampler = NULL,
           model_rng = NULL)
    })
  }
  list(state = state,
       model = model,
       sampler = sampler,
       observer = observer,
       runner = runner)
}


direct_sample_within_domain <- function(model, rng, max_attempts = 100) {
  for (i in seq_len(max_attempts)) {
    x <- model$direct_sample(rng)
    if (all(x >= model$domain[, 1] & x <= model$domain[, 2])) {
      return(x)
    }
  }
  cli::cli_abort(
    c("Failed to sample initial conditions within {max_attempts} attempt{?s}",
      i = paste("Your model's 'direct_sample()' method is generating",
                "samples that fall outside your model's domain.  Probably",
                "you should fix one or both of these!")))
}
