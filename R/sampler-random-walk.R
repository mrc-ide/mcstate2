##' Create a simple random walk sampler, which uses a symmetric
##' proposal to move around parameter space.  This sampler supports
##' sampling from models where the likelihood is only computable
##' randomly (e.g., for pmcmc).
##'
##' @title Random Walk Sampler
##'
##' @param vcv A variance covariance matrix for the proposal.
##'
##' @return A `mcstate_sampler` object, which can be used with
##'   [mcstate_sample]
##'
##' @export
mcstate_sampler_random_walk <- function(vcv = NULL) {
  check_vcv(vcv, allow_3d = TRUE, call = environment())
  internal <- new.env()

  initialise <- function(pars, model, observer, rng) {
    n_pars <- length(model$parameters)
    internal$multiple_parameters <- length(dim2(pars)) > 1
    if (internal$multiple_parameters) {
      ## this is enforced elsewhere
      stopifnot(model$properties$allow_multiple_parameters)
    }

    vcv <- sampler_validate_vcv(vcv, pars)
    internal$proposal <- make_rmvnorm(vcv)
    initialise_state(pars, model, observer, rng)
  }

  step <- function(state, model, observer, rng) {
    pars_next <- internal$proposal(state$pars, rng)
    density_next <- model$density(pars_next)
    accept <- density_next - state$density > log(rng$random_real(1))
    state <- update_state(state, pars_next, density_next, accept,
                          model, observer, rng)
    state
  }

  ## These are all effectively the defaults:
  finalise <- function(state, model, rng) {
    NULL
  }

  get_internal_state <- function() {
    NULL
  }

  set_internal_state <- function(state) {
  }

  mcstate_sampler("Random walk",
                  initialise,
                  step,
                  finalise,
                  get_internal_state,
                  set_internal_state)
}
