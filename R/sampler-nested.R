##' Create a nested random walk sampler, which uses a symmetric
##' proposal for separable sections of a model to move around in
##' parameter space.  This sampler supports sampling from models where
##' the likelihood is only computable randomly (e.g., for pmcmc).
##'
##' @title Nested Random Walk Sampler
##'
##' @param vcv A list of variance covariance matrices.  The format
##'   here is subject to change.
##'
##' @return A `mcstate_sampler` object, which can be used with
##'   [mcstate_sample]
##'
##' @export
mcstate_sampler_nested_random_walk <- function(vcv = NULL) {
  if (!is.list(vcv)) {
    cli::cli_abort(
      c("Expected a list for 'vcv'",
        i = "We need a list where each element is itself a vcv for a subgroup",
        i = "Later we might change this to allow a block-structured matrix"),
      arg = 'vcv')
  }

  for (i in seq_along(vcv)) {
    ## TODO: this might not create the best error message, but it's
    ## likely to be clear enough for now.

    ## TODO: we should allow the first element to be empty?
    check_vcv(vcv[[i]], call = environment())
  }

  proposal <- lapply(vcv, make_rmvnorm)

  initialise <- function(pars, model, rng) {
    ## Some checking required here for these:
    ## model$parameters
    ## model$groups
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
    pars_next <- proposal$base(state$pars, rng)
    density_next <- model$density(pars_next, by_group = TRUE)
    density_by_group_next <- attr(density_next, "by_group")
    accept <- density_next - state$density > log(rng$random_real(1))
    if (accept) {
      state$pars <- pars_next
      state$density <- density_next
      internal$density_by_group <- density_by_group_next
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
