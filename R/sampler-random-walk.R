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
      n_sets <- ncol(pars)
    } else {
      n_sets <- 1L
    }
    if (n_sets > 1 && !model$properties$allow_multiple_parameters) {
      cli::cli_abort(
        "Using multiple parameter sets with model that does not support it")
    }

    n_pars_vcv <- nrow(vcv)
    if (n_pars != n_pars_vcv) {
      cli::cli_abort(
        "Incompatible length parameters ({n_pars}) and vcv ({n_pars_vcv})")
    }

    if (length(dim(vcv)) == 3) {
      n_sets_vcv <- dim(vcv)[[3]]
      if (n_sets == 1) {
        dim(vcv) <- dim(vcv)[1:2]
      } else if (n_sets_vcv == 1) {
        vcv <- array(vcv, c(dim(vcv)[1:2], n_sets))
      } else if (n_sets_vcv != n_sets) {
        cli::cli_abort(
          paste("Incompatible number of parameter sets ({n_sets}) and slices",
                "in vcv ({n_sets_vcv})"))
      }
    } else if (n_sets > 1) {
      vcv <- array(vcv, c(dim(vcv)[1:2], n_sets))
    }

    internal$proposal <- make_rmvnorm(vcv)
    initialise_state(pars, model, observer, rng)
  }

  step <- function(state, model, observer, rng) {
    pars_next <- internal$proposal(state$pars, rng)
    density_next <- model$density(pars_next)
    accept <- density_next - state$density > log(rng$random_real(1))

    if (any(accept)) {
      if (internal$multiple_parameters) {
        state <- update_state_multiple(state, accept, pars_next, density_next,
                                       model, observer, rng)
      } else {
        state <- update_state(state, pars_next, density_next,
                              model, observer, rng)
      }
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
