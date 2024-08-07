##' Create a simple random walk sampler, which uses a symmetric
##' proposal to move around parameter space.  This sampler supports
##' sampling from models where the likelihood is only computable
##' randomly (e.g., for pmcmc).
##'
##' @title Random Walk Sampler
##'
##' @param vcv A variance covariance matrix for the proposal.
##'
##' @param boundaries Control the behaviour of proposals that are
##'   outside the model domain.  The supported options are:
##'
##'   * "reflect" (the default): we reflect proposed parameters that
##'     lie outside the domain back into the domain (as many times as
##'     needed)
##'
##'   * "reject": we do not evaluate the density function, and return
##'     `-Inf` for its density instead.
##'
##'   * "ignore": evaluate the point anyway, even if it lies outside
##'     the domain.
##'
##' The initial point selected will lie within the domain, as this is
##' enforced by [mcstate_sample].
##'
##' @return A `mcstate_sampler` object, which can be used with
##'   [mcstate_sample]
##'
##' @export
mcstate_sampler_random_walk <- function(vcv = NULL, boundaries = "reflect") {
  check_vcv(vcv, allow_3d = TRUE, call = environment())
  internal <- new.env()

  boundaries <- match_value(boundaries, c("reflect", "reject", "ignore"))

  initialise <- function(pars, model, observer, rng) {
    n_pars <- length(model$parameters)
    internal$multiple_parameters <- length(dim2(pars)) > 1
    if (internal$multiple_parameters) {
      ## this is enforced elsewhere
      stopifnot(model$properties$allow_multiple_parameters)
    }

    vcv <- sampler_validate_vcv(vcv, pars)
    internal$proposal <-
      make_random_walk_proposal(vcv, model$domain, boundaries)
    initialise_state(pars, model, observer, rng)
  }

  step <- function(state, model, observer, rng) {
    browser()
    pars_next <- internal$proposal(state$pars, rng)
    reject_some <- boundaries == "reject" &&
      !all(i <- is_parameters_in_domain(pars_next, model$domain))
    if (reject_some) {
      density_next <- rep(-Inf, length(state$density))
      if (any(i)) {
        density_next[i] <- model$density(pars_next[, i, drop = FALSE])
      }
    } else {
      density_next <- model$density(pars_next)
    }
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


make_random_walk_proposal <- function(vcv, domain, boundaries) {
  mvn <- make_rmvnorm(vcv)
  if (boundaries != "reflect" || !any(is.finite(domain))) {
    return(mvn)
  }

  x_min <- domain[, 1]
  x_max <- domain[, 2]
  function(x, rng) {
    reflect_proposal(mvn(x, rng), x_min, x_max)
  }
}


## create function to reflect proposal boundaries at pars_min and pars_max
## this ensures the proposal is symmetrical and we can simplify the MH step
reflect_proposal <- function(x, x_min, x_max) {
  i <- x < x_min | x > x_max
  if (any(i)) {
    i_both <- i & is.finite(x_min) & is.finite(x_max)
    i_min <- i & is.finite(x_min) & !is.finite(x_max)
    i_max <- i & !is.finite(x_min) & is.finite(x_max)
    if (is.matrix(x)) {
      ## This will be a bit different if the model supports array
      ## specifications of the domain, which we might want to support
      ## later I guess. Replicating these helps tidy up the
      ## bookkeeping, and also means we don't have to think about this
      ## very much.
      x_min <- array(x_min, dim(x))
      x_max <- array(x_max, dim(x))
    }
    x[i_both] <- reflect_proposal_both(x[i_both], x_min[i_both], x_max[i_both])
    x[i_min] <- reflect_proposal_one(x[i_min], x_min[i_min])
    x[i_max] <- reflect_proposal_one(x[i_max], x_max[i_max])
  }
  x
}


reflect_proposal_both <- function(x, x_min, x_max) {
  x_r <- x_max - x_min
  abs((x + x_r - x_min) %% (2 * x_r) - x_r) + x_min
}


reflect_proposal_one <- function(x, x_bound) {
  2 * x_bound - x
}


is_parameters_in_domain <- function(x, domain) {
  x_min <- domain[, 1]
  x_max <- domain[, 2]
  i <- x > x_min & x < x_max
  if (all(i)) {
    if (is.matrix(x)) rep(TRUE, ncol(x)) else TRUE
  } else {
    if (is.matrix(x)) apply(i, 2, all) else FALSE
  }
}
