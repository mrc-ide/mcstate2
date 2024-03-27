initialise_state <- function(pars, model, rng) {
  initialise_rng_state(model, rng)
  list(pars = pars, density = model$density(pars))
}


initialise_rng_state <- function(model, rng) {
  if (isTRUE(model$properties$is_stochastic)) {
    model$rng_state$set(mcstate_rng$new(rng$state())$jump()$state())
  }
}
