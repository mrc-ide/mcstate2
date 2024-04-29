initialise_state <- function(pars, model, observer, rng) {
  initialise_rng_state(model, rng)
  density <- model$density(pars)
  if (is.null(observer)) {
    observation <- NULL
  } else {
    observation <- observer$observe(model$model, rng)
  }
  list(pars = pars, density = density, observation = observation)
}


update_state <- function(state, pars, density, accept, model, observer, rng) {
  if (any(accept)) {
    if (is.matrix(pars)) {
      stopifnot(is.null(observer)) # Enforced earlier
      state$pars[, accept] <- pars[, accept]
      state$density[accept] <- density[accept]
      state
    } else {
      state$pars <- pars
      state$density <- density
      if (!is.null(observer)) {
        state$observation <- observer$observe(model$model, rng)
      }
    }
  }
  state
}


initialise_rng_state <- function(model, rng) {
  if (isTRUE(model$properties$is_stochastic)) {
    model$rng_state$set(mcstate_rng$new(rng$state())$jump()$state())
  }
}
