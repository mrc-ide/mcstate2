initialise_state <- function(pars, model, rng) {
  initialise_rng_state(model, rng)
  density <- model$density(pars)
  if (model$properties$has_observer) {
    observation <- model$observer$observe()
  } else {
    observation <- NULL
  }
  list(pars = pars, density = density, observation = observation)
}


update_state <- function(state, pars, density, accept, model, rng) {
  if (any(accept)) {
    if (is.matrix(state$pars)) {
      state$pars[, accept] <- pars[, accept]
      state$density[accept] <- density[accept]
      state
    } else {
      state$pars <- pars
      state$density <- density
      if (model$properties$has_observer) {
        state$observation <- model$observer$observe()
      }
    }
  }
  state
}


initialise_rng_state <- function(model, rng) {
  if (isTRUE(model$properties$is_stochastic)) {
    model$rng_state$set(monty_rng$new(rng$state())$jump()$state())
  }
}


sampler_validate_vcv <- function(vcv, pars, call = NULL) {
  if (is.matrix(pars)) {
    n_sets <- ncol(pars)
    n_pars <- nrow(pars)
  } else {
    n_sets <- 1
    n_pars <- length(pars)
  }

  n_pars_vcv <- nrow(vcv)
  if (n_pars != n_pars_vcv) {
    cli::cli_abort(
      "Incompatible length parameters ({n_pars}) and vcv ({n_pars_vcv})",
      call = call)
  }

  if (length(dim(vcv)) == 3) {
    n_sets_vcv <- dim(vcv)[[3]]
    if (n_sets == 1 && n_sets_vcv == 1) {
      dim(vcv) <- dim(vcv)[1:2]
    } else if (n_sets_vcv == 1) {
      vcv <- array(vcv, c(dim(vcv)[1:2], n_sets))
    } else if (n_sets_vcv != n_sets) {
      cli::cli_abort(
        c(paste("Incompatible number of parameter sets ({n_sets}) and slices",
                "in vcv ({n_sets_vcv})"),
          i = paste("You configured the sampler to expect {n_sets_vcv}",
                    "parameter set{?s}, but {n_sets} parameter set{?s}",
                    "were provided when the sampler was initialised")),
        call = call)
    }
  } else if (n_sets > 1) {
    vcv <- array(vcv, c(dim(vcv)[1:2], n_sets))
  }

  vcv
}
