##' Create a nested random walk sampler, which uses a symmetric
##' proposal for separable sections of a model to move around in
##' parameter space.  This sampler supports sampling from models where
##' the likelihood is only computable randomly (e.g., for pmcmc), and
##' requires that models support the `has_parameter_groups` property.
##'
##' The intended use case for this sampler is for models where the
##' density can be decomposed at least partially into chunks that are
##' independent from each other.  Our motivating example for this is a
##' model of COVID-19 transmission where some parameters
##' region-specific (e.g., patterns and rates of contact between
##' individuals), and some parameters are shared across all regions
##' (e.g., intrinsic properties of the disease such as incubation
##' period).
##'
##' The strategy is to propose all the shared parameters as a
##' deviation from the current point in parameter space as a single
##' move and accept or reject as a block. Then we generate points for
##' all the region-specific parameters, compute the density and then
##' accept or reject these updates independently.  This is possible
##' because the change in likelihood in region A is independent from
##' region B.
##'
##' We expect that this approach will be beneficial in limited
##' situations, but where it is beneficial it is likely to result in
##' fairly large speed-ups:
##'
##' * You probably need more than three regions; as the number of
##'   regions increases the benefit of independently accepting or
##'   rejecting densities increases (with 1000 separate regions your
##'   chains will mix very slowly for example).
##' * Your model is fairly computationally heavy so that the density
##'   calculation completely dominates the sampling process.
##' * You do not have access to gradient information for your model;
##'   we suspect that HMC will outperform this approach by some margin
##'   because it already includes this independence via the gradients.
##' * You can compute your independent calculations in parallel, which
##'   help this method reduce your walk time.
##'
##' @title Nested Random Walk Sampler
##'
##' @param vcv A list of variance covariance matrices.  We expect this
##'   to be a list with elements `base` and `groups` corresponding to
##'   the covariance matrix for base parameters (if any) and groups.
##'
##' @inheritParams monty_sampler_random_walk
##'
##' @return A `monty_sampler` object, which can be used with
##'   [monty_sample]
##'
##' @export
monty_sampler_nested_random_walk <- function(vcv, boundaries = "reflect",
                                             rerun_every = Inf,
                                             rerun_random = FALSE) {
  if (!is.list(vcv)) {
    cli::cli_abort(
      "Expected a list for 'vcv'",
      arg = "vcv")
  }

  if (!setequal(names(vcv), c("base", "groups"))) {
    cli::cli_abort("Expected 'vcv' to have elements 'base' and 'groups'",
                   arg = "vcv")
  }
  if (!is.null(vcv$base)) {
    check_vcv(vcv$base, allow_3d = TRUE, call = environment())
  }
  if (!is.list(vcv$groups)) {
    cli::cli_abort("Expected 'vcv$groups' to be a list")
  }
  if (length(vcv$groups) < 1) {
    cli::cli_abort("Expected 'vcv$groups' to have at least one element")
  }
  for (i in seq_along(vcv$groups)) {
    check_vcv(vcv$groups[[i]], allow_3d = TRUE,
              name = sprintf("vcv$groups[%d]", i), call = environment())
  }

  internal <- new.env(parent = emptyenv())

  boundaries <- match_value(boundaries, c("reflect", "reject", "ignore"))
  
  if (!identical(unname(rerun_every), Inf)) {
    assert_scalar_positive_integer(rerun_every)
  }
  assert_scalar_logical(rerun_random)

  initialise <- function(pars, model, rng) {
    if (!model$properties$has_parameter_groups) {
      cli::cli_abort("Your model does not have parameter groupings")
    }

    internal$multiple_parameters <- length(dim2(pars)) > 1
    if (internal$multiple_parameters) {
      ## this is enforced elsewhere
      stopifnot(model$properties$allow_multiple_parameters)
    }
    
    internal$rerun <- 
      make_rerun(rerun_every, rerun_random, model$properties$is_stochastic)
    internal$step <- 0

    vcv <- sampler_validate_nested_vcv(vcv, model$parameter_groups, pars)

    internal$proposal <- nested_proposal(vcv, model$parameter_groups, pars,
                                         model$domain, boundaries)

    initialise_rng_state(model, rng)
    density <- model$density(pars, by_group = TRUE)
    density_by_group <- attr(density, "by_group")
    n_groups <- max(model$parameter_groups)

    if (is.null(density_by_group)) {
      cli::cli_abort(
        c(paste("model$density(x, by_group = TRUE) did not produce a",
                "density with a 'by_group' attribute"),
          i = paste("I expected an attribute 'by_group' with {n_groups}",
                    "elements corresponding to parameter groups to be",
                    "included with your density")))
    }
    if (dim2(density_by_group)[1] != n_groups) {
        cli::cli_abort(
          paste("model$density(x, by_group = TRUE) produced a 'by_group'",
                "attribute with incorrect length {length(density_by_group)}",
                "but I expected length {n_groups}"))
    }

    internal$density_by_group <- density_by_group
    state <- list(pars = pars, density = c(density))
    ## TODO: we need to fix observation here; it should move into a
    ## helper as part of a setup I think; see the sampler-helpers for
    ## the single-parameter case.
    state
  }

  ## There are probably different modes that this could run in, they'd
  ## be fairly easy to change.  This one would correspond to some sort
  ## of "full update" mode where everything is done within a step, but
  ## we could also do one where we allow for picking one update type
  ## with some schedule or probability and applying that, which would
  ## allow for faster movement of some part of the chain.  We could
  ## handle this by additional arguments to the constructor, then
  ## either changing the behaviour of the step function or swapping in
  ## a different version.
  step <- function(state, model, rng) {
    internal$step <- internal$step + 1
    rerun <- internal$rerun(internal$step, rng)
    if (any(rerun)) {
      ## This is currently just setup assuming we are not using multiple
      ## parameters as currently they cannot be used with stochastic models,
      ## while the rerun is only used with stochastic models
      density <- model$density(state$pars, by_group = TRUE)
      state$density <- c(density)
      internal$density_by_group <- attr(density, "by_group")
      ## TODO: observe here
    }
    
    if (!is.null(internal$proposal$base)) {
      pars_next <- internal$proposal$base(state$pars, rng)

      reject_some <- boundaries == "reject" &&
        !all(i <- is_parameters_in_domain(pars_next, model$domain))
      if (reject_some) {
        density_next <- rep(-Inf, length(state$density))
        density_by_group_next <- array(-Inf, dim2(internal$density_by_group))
        if (any(i)) {
          density_next_i <- model$density(pars_next[, i, drop = FALSE],
                                          by_group = TRUE)
          density_next[i] <- density_next_i
          density_by_group_next[, i] <- attr(density_next_i, "by_group")
        }
      } else {
        density_next <- model$density(pars_next, by_group = TRUE)
        density_by_group_next <- attr(density_next, "by_group")
      }

      accept <- density_next - state$density > log(rng$random_real(1))
      if (any(accept)) {
        if (!all(accept)) {
          ## Retain some older parameters
          i <- which(!accept)
          pars_next[, i] <- state$pars[, i]
          density_next <- model$density(pars_next, by_group = TRUE)
          density_by_group_next <- attr(density_next, "by_group")
        }
        state$pars <- pars_next
        state$density <- density_next
        internal$density_by_group <- density_by_group_next
        ## TODO: observe here
      }
    }

    pars_next <- internal$proposal$groups(state$pars, rng)

    reject_some <- boundaries == "reject" &&
      !all(i <- is_parameters_in_domain_groups(pars_next, model$domain,
                                               model$parameter_groups))

    ## This bit is potentially inefficient - for any proposed parameters out of
    ## bounds I substitute in the current parameters, so that we can run the
    ## density on all groups. Ideally we would want to only run the density on
    ## groups with all parameters in bounds. A bit fiddly to do that in a nice
    ## way when doing simultaneous sampling
    if (reject_some) {
      density_next <- rep(-Inf, length(state$density))
      density_by_group_next <- array(-Inf, dim2(internal$density_by_group))
      if (any(i)) {
        if (internal$multiple_parameters) {
          for (j in seq_len(ncol(i))) {
            if (!all(i[, j])) {
              i_group <- model$parameter_groups %in% which(!i[, j])
              pars_next[i_group, j] <- state$pars[i_group, j]
            }
          }
        } else {
          i_group <- model$parameter_groups %in% which(!i)
          pars_next[i_group] <- state$pars[i_group]
        }
        density_next <- model$density(pars_next, by_group = TRUE)
        density_by_group_next <- attr(density_next, "by_group")
      }
    } else {
      density_next <- model$density(pars_next, by_group = TRUE)
      density_by_group_next <- attr(density_next, "by_group")
    }

    accept <- density_by_group_next - internal$density_by_group >
      log(rng$random_real(dim2(density_by_group_next)[1]))

    if (any(accept)) {
      if (!all(accept)) {
        ## Retain some older parameters
        if (internal$multiple_parameters) {
          for (j in seq_len(ncol(accept))) {
            if (!all(accept[, j])) {
              i <- model$parameter_groups %in% which(!accept[, j])
              pars_next[i, j] <- state$pars[i, j]
            }
          }
        } else {
          i <- model$parameter_groups %in% which(!accept)
          pars_next[i] <- state$pars[i]
        }
        ## If e.g. density is provided by a particle filter, would this bit
        ## mean rerunning it? Increases time cost if so, and would result in
        ## new value of density (different to that which was accepted)
        density_next <- model$density(pars_next, by_group = TRUE)
        density_by_group_next <- attr(density_next, "by_group")
      }
      state$pars <- pars_next
      state$density <- c(density_next)
      internal$density_by_group <- density_by_group_next
      ## TODO: observe here
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

  monty_sampler("Nested random walk",
                "monty_sampler_nested_random_walk",
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
      paste("Expected '{name}' to have length {n_pars}, but it had length",
            "{length(x)}"),
      call = call)
  }
  if (min(x) < 0) {
    cli::cli_abort("Invalid negative group in '{name}'", call = call)
  }
  n_groups <- max(x)
  msg <- setdiff(seq_len(n_groups), x)
  if (length(msg) > 0) {
    cli::cli_abort(
      c("Missing groups from '{name}'",
        i = paste("I expected all integers from 1 to {n_groups} to be present",
                  "in your parameter groups vector, but you are missing",
                  "{msg}")),
      call = call)
  }
}


sampler_validate_nested_vcv <- function(vcv, parameter_groups, pars, 
                                        call = NULL) {
  i_base <- parameter_groups == 0
  n_base <- sum(i_base)
  has_base <- n_base != 0
  n_groups <- max(parameter_groups)
  i_group <- lapply(seq_len(n_groups), function(i) which(parameter_groups == i))
  
  if (NROW(vcv$base) != n_base) {
    cli::cli_abort(
      c("Incompatible number of base parameters in your model and sampler",
        i = paste("Your model has {n_base} base parameters, but 'vcv$base'",
                  "implies {NROW(vcv$base)} parameters")),
      call = call)
  }
  if (length(vcv$groups) != n_groups) {
    cli::cli_abort(
      c("Incompatible number of parameter groups in your model and sampler",
        i = paste("Your model has {n_groups} parameter groups, but",
                  "'vcv$groups' has {length(vcv$groups)} groups")),
      call = call)
  }
  n_pars_by_group <- lengths(i_group)
  n_pars_by_group_vcv <- vnapply(vcv$groups, nrow)
  err <- n_pars_by_group_vcv != n_pars_by_group
  if (any(err)) {
    detail <- sprintf(
      "Group %d has %d parameters but 'vcv$groups[[%d]]' has %d",
      which(err), n_pars_by_group[err],
      which(err), n_pars_by_group_vcv[err])
    cli::cli_abort(
      c("Incompatible number of parameters within parameter group",
        set_names(detail, "i")),
      call = call)
  }
  
  if (has_base) {
    if (is.matrix(pars)) {
      vcv$base <- sampler_validate_vcv(vcv$base, pars[i_base, , drop = FALSE])
    } else {
      vcv$base <- sampler_validate_vcv(vcv$base, pars[i_base])
    }
  }
  
  for (i in seq_len(n_groups)) {
    if (is.matrix(pars)) {
      vcv$groups[[i]] <- 
        sampler_validate_vcv(vcv$groups[[i]],
                             pars[i_group[[i]], , drop = FALSE])
    } else {
      vcv$groups[[i]] <- 
        sampler_validate_vcv(vcv$groups[[i]], pars[i_group[[i]]])
    }
  }
  
  vcv
}


nested_proposal <- function(vcv, parameter_groups, pars, domain,
                            boundaries = "reflect") {
  i_base <- parameter_groups == 0
  n_base <- sum(i_base)
  has_base <- n_base != 0
  n_groups <- max(parameter_groups)
  i_group <- lapply(seq_len(n_groups), function(i) which(parameter_groups == i))
  
  if (has_base) {
    if (is.matrix(pars)) {
      mvn_base <- make_random_walk_proposal(
        vcv$base, domain[i_base, , drop = FALSE], boundaries)
      proposal_base <- function(x, rng) {
        ## This approach is likely to be a bit fragile, so we'll
        ## probably want some naming related verification here soon too.
        x[i_base] <- mvn_base(x[i_base, ], rng)
        x
      }
    } else {
      mvn_base <- make_random_walk_proposal(
        vcv$base, domain[i_base, , drop = FALSE], boundaries)
      proposal_base <- function(x, rng) {
        ## This approach is likely to be a bit fragile, so we'll
        ## probably want some naming related verification here soon too.
        x[i_base] <- mvn_base(x[i_base], rng)
        x
      }
    }

  } else {
    proposal_base <- NULL
  }

  mvn_groups <- lapply(seq_len(n_groups), function(i) {
    make_random_walk_proposal(vcv$groups[[i]],
                              domain[i_group[[i]], , drop = FALSE],
                              boundaries)
  })
  proposal_groups <- function(x, rng) {
    for (i in seq_len(n_groups)) {
      if (is.matrix(x)) {
        x[i_group[[i]], ] <- mvn_groups[[i]](x[i_group[[i]], ], rng)
      } else {
        x[i_group[[i]]] <- mvn_groups[[i]](x[i_group[[i]]], rng)
      }
    }
    x
  }

  list(base = proposal_base,
       groups = proposal_groups)
}


is_parameters_in_domain_groups <- function(x, domain, parameter_groups) {
  x_min <- domain[, 1]
  x_max <- domain[, 2]
  i <- x > x_min & x < x_max
  n_groups <- max(parameter_groups)
  i_group <- lapply(seq_len(n_groups), function(i) which(parameter_groups == i))
  if (is.matrix(x)) {
    t(vapply(i_group, function(j) apply(i[j, , drop = FALSE], 2, all),
             logical(ncol(x))))
  } else {
    vapply(i_group, function(j) all(i[j]), logical(1L))
  }
}
