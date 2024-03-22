##' Create a nested random walk sampler, which uses a symmetric
##' proposal for separable sections of a model to move around in
##' parameter space.  This sampler supports sampling from models where
##' the likelihood is only computable randomly (e.g., for pmcmc).
##'
##' @title Nested Random Walk Sampler
##'
##' @param vcv A list of variance covariance matrices.  We expect this
##'   to be a list with elements `base` and `groups` corresponding to
##'   the covariance matrix for base parameters (if any) and groups.
##'
##' @return A `mcstate_sampler` object, which can be used with
##'   [mcstate_sample]
##'
##' @export
mcstate_sampler_nested_random_walk <- function(vcv) {
  if (!is.list(vcv)) {
    cli::cli_abort(
      "Expected a list for 'vcv'",
      arg = 'vcv')
  }

  internal <- new.env(parent = emptyenv())  
  proposal <- nested_proposal(vcv, environment())

  initialise <- function(pars, model, rng) {
    if (proposal$n_pars != length(model$parameters)) {
      cli::cli_abort(
        "Incompatible length parameters and vcv")
    }
    if (isTRUE(model$properties$is_stochastic)) {
      model$model$set_rng_state(rng)
    }
    density <- model$density(pars, by_group = TRUE)
    internal$density_by_group <- attr(density, "by_group")
    list(pars = pars, density = density)
  }

  ## TODO: There are probably different modes that this could run in,
  ## they'd be fairly easy to change.  This one would correspond to
  ## some sort of "full update" mode where everything is done within a
  ## step, but we could also do one where we allow for picking one
  ## update type with some schedule or probability and applying that,
  ## which would allow for faster movement of some part of the chain.
  step <- function(state, model, rng) {
    if (proposal$has_base) {
      pars_next <- proposal$base(state$pars, rng)
      density_next <- model$density(pars_next, by_group = TRUE)
      density_by_group_next <- attr(density_next, "by_group")
      accept <- density_next - state$density > log(rng$random_real(1))
      if (accept) {
        state$pars <- pars_next
        state$density <- density_next
        internal$density_by_group <- density_by_group_next
      }
    }

    pars_next <- proposal$groups(state$pars, rng)
    density_next <- model$density(pars_next, by_group = TRUE)
    density_by_group_next <- attr(density_next, "by_group")
    accept <- density_by_group_next - internal$density_by_grop >
      log(rng$random_real(length(density_by_group_next)))

    if (any(accept)) {
      if (!all(accept)) {
        ## Retain some older parameters
        i <- model$parameter_groups %in% which(!accept)
        pars_next[i] <- state$pars[i]
        density_next <- model$density(pars_next, by_group = TRUE)
        density_by_group_next <- attr(density_next, "by_group")
      }
      state$pars <- pars_next
      state$density <- density_next
      internal$density_by_group <- density_by_group_next
    }
    state
  }

  finalise <- function(state, model, rng) {
    NULL
  }

  get_internal_state <- function() {
    as.list(internal)
  }

  set_internal_state <- function(state) {
    list2env(state, internal)
  }

  mcstate_sampler("Nested random walk",
                  initialise,
                  step,
                  finalise,
                  get_internal_state,
                  set_internal_state)
}


check_parameter_groups <- function(x, n_pars, name = deparse(substitute(x)),
                                   call = NULL) {
  if (!rlang::is_integerish(x)) {
    cli::cli_abort("Expected '{name}' to be integer-like", call = call)
  }
  if (length(x) != n_pars) {
    cli::cli_abort(
      "Expected '{name}' to have length {n_pars}, but was {length(x)}",
      call = call)
  }
  n_groups <- max(n_pars)
  msg <- setdiff(seq_len(n_pars), n_groups)
  if (length(msg) > 0) {
    # TODO: Better error here that explains the situation better and
    # offers a hint as to what to do.
    cli::cli_abort("Missing groups from '{name}'", call = call)
  }
  if (min(x) < 0) {
    cli::cli_abort("Invalid negative groups in '{name}'", call = call)
  }
}


nested_proposal <- function(vcv, call = NULL) {
  if (!setequal(names(vcv), c("base", "groups"))) {
    cli::cli_abort("Expected 'vcv' to have elements 'base' and 'groups'",
                   arg = "vcv")
  }
  has_base <- !is.null(vcv$base)  
  if (has_base) {
    check_vcv(vcv$base, call = environment())
  }
  if (length(vcv$groups) < 1) {
    cli::cli_abort("Expected at least 1 group")
  }
  if (!is.list(vcv$groups)) {
    cli::cli_abort("Expected 'vcv$groups' to be a list")
  }
  for (i in seq_along(vcv$base)) {
    check_vcv(vcv$groups[[i]], call = environment())
  }
  list(
    base = if (has_base) make_rmvnorm(vcv$base),
    groups = lapply(vcv$base, make_rmvnorm),
    has_base = has_base,
    n_groups = length(vcv$groups),
    n_pars = nrow(vcv$base) + sum(vnapply(vcv$groups, nrow)))
}
