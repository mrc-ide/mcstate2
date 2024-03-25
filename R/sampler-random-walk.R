##' Create a simple random walk sampler, which uses a symmetric
##' proposal to move around parameter space.  This sampler supports
##' sampling from models where the likelihood is only computable
##' randomly (e.g., for pmcmc).
##'
##' @title Random Walk Sampler
##'
##' @param proposal A proposal function; must take a vector of
##'   parameters and a random number generator object
##'   ([`mcstate_rng`]) and produce a new vector of proposed
##'   parameters.
##'
##' @param vcv A variance covariance matrix to generate a `proposal`
##'   function from.  If you want multivariate Gaussian proposal, this
##'   is likely simpler than supplying your own `proposal`, and
##'   generally more efficient too.
##'
##' @return A `mcstate_sampler` object, which can be used with
##'   [mcstate_sample]
##'
##' @export
mcstate_sampler_random_walk <- function(proposal = NULL, vcv = NULL) {
  if (is.null(proposal) && is.null(vcv)) {
    cli::cli_abort("One of 'proposal' or 'vcv' must be given")
  }
  if (!is.null(proposal) && !is.null(vcv)) {
    cli::cli_abort("Only one of 'proposal' or 'vcv' may be given")
  }
  if (!is.null(vcv)) { # proposal is null
    check_vcv(vcv, call = environment())
    proposal <- make_rmvnorm(vcv)
  }

  initialise <- function(pars, model, observer, rng) {
    if (!is.null(vcv)) {
      n_pars <- length(model$parameters)
      n_vcv <- nrow(vcv)
      if (n_pars != n_vcv) {
        cli::cli_abort(
          "Incompatible length parameters ({n_pars}) and vcv ({n_vcv})")
      }
    }
    initialise_state(pars, model, observer, rng)
  }

  step <- function(state, model, observer, rng) {
    pars_next <- proposal(state$pars, rng)
    density_next <- model$density(pars_next)

    if (density_next - state$density > log(rng$random_real(1))) {
      state <- update_state(state, pars_next, density_next,
                            model, observer, rng)
    }
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
