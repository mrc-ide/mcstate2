##' Create a nested adaptive Metropolis-Hastings sampler, which extends the
##' adaptive sampler [monty_sampler_adaptive], tuning the variance covariance
##' matrices for proposal for the separable sections
##' of a nested model (vs the simple nested random walk sampler
##' [monty_sampler_random_walk]). This sampler requires
##' that models support the `has_parameter_groups` property.
##'
##' Much like the simple nested random walk sampler
##' [monty_sampler_random_walk], the strategy is to propose all the
##' shared parameters as a deviation from the current point in parameter space
##' as a single move and accept or reject as a block. Then we generate points
##' for all the region-specific parameters, compute the density and then
##' accept or reject these updates independently.  This is possible
##' because the change in likelihood in region A is independent from
##' region B.
##'
##' The adaptive proposal algorithm of the non-nested adaptive sampler
##' [monty_sampler_adaptive] is extended here to adaptively tune the variance
##' covariance matrix of each of these parameter chunks.
##'
##' @title Nested Adaptive Metropolis-Hastings Sampler
##'
##' @inheritParams monty_sampler_adaptive
##'
##' @return A `monty_sampler` object, which can be used with
##'   [monty_sample]
##'
##' @export
monty_sampler_nested_adaptive <- function(initial_vcv,
                                          initial_vcv_weight = 1000,
                                          initial_scaling = 1,
                                          initial_scaling_weight = NULL,
                                          min_scaling = 0,
                                          scaling_increment = NULL,
                                          log_scaling_update = TRUE,
                                          acceptance_target = 0.234,
                                          forget_rate = 0.2,
                                          forget_end = Inf,
                                          adapt_end = Inf,
                                          pre_diminish = 0,
                                          boundaries = "reflect") {
  if (!is.list(initial_vcv)) {
    cli::cli_abort(
      "Expected a list for 'initial_vcv'",
      arg = "initial_vcv")
  }

  if (!setequal(names(initial_vcv), c("base", "groups"))) {
    cli::cli_abort(
      "Expected 'initial_vcv' to have elements 'base' and 'groups'",
      arg = "initial_vcv")
  }
  if (!is.null(initial_vcv$base)) {
    check_vcv(initial_vcv$base, allow_3d = TRUE, call = environment())
  }
  if (!is.list(initial_vcv$groups)) {
    cli::cli_abort("Expected 'initial_vcv$groups' to be a list")
  }
  if (length(initial_vcv$groups) < 1) {
    cli::cli_abort(
      "Expected 'initial_vcv$groups' to have at least one element")
  }
  for (i in seq_along(initial_vcv$groups)) {
    check_vcv(initial_vcv$groups[[i]], allow_3d = TRUE,
              name = sprintf("initial_vcv$groups[%d]", i),
              call = environment())
  }

  internal <- new.env(parent = emptyenv())
  
  boundaries <- match_value(boundaries, c("reflect", "reject", "ignore"))

  initialise <- function(pars, model, rng) {
    require_deterministic(model,
                          "Can't use adaptive sampler with stochastic models")

    if (!model$properties$has_parameter_groups) {
      cli::cli_abort("Your model does not have parameter groupings")
    }
    
    internal$multiple_parameters <- length(dim2(pars)) > 1
    if (internal$multiple_parameters) {
      ## this is enforced elsewhere
      stopifnot(model$properties$allow_multiple_parameters)
    }

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
    
    initial_vcv <- 
      sampler_validate_nested_vcv(initial_vcv, model$parameter_groups, pars)
    
    if (internal$multiple_parameters) {
      internal$adaptive <-
        Map(initialise_nested_adaptive, 
            lapply(asplit(pars, 2), c),
            split_nested_vcv(initial_vcv, pars),
            MoreArgs = list(initial_vcv_weight = initial_vcv_weight,
                            initial_scaling = initial_scaling,
                            initial_scaling_weight = initial_scaling_weight,
                            min_scaling = min_scaling,
                            scaling_increment = scaling_increment,
                            log_scaling_update = log_scaling_update,
                            acceptance_target = acceptance_target,
                            forget_rate = forget_rate,
                            forget_end = forget_end,
                            adapt_end = adapt_end,
                            pre_diminish = pre_diminish,
                            parameter_groups = model$parameter_groups)
        )
    } else {
      internal$adaptive <-
        initialise_nested_adaptive(pars,
                            initial_vcv,
                            initial_vcv_weight,
                            initial_scaling,
                            initial_scaling_weight,
                            min_scaling,
                            scaling_increment,
                            log_scaling_update,
                            acceptance_target,
                            forget_rate,
                            forget_end,
                            adapt_end,
                            pre_diminish,
                            model$parameter_groups)
    }

    state <- list(pars = pars, density = c(density))
    ## TODO: mrc-5862
    if (model$properties$has_observer) {
      state$observation <- model$observer$observe()
    }
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
    proposal <- 
      nested_proposal_adaptive(internal$adaptive, internal$multiple_parameters,
                               model$parameter_groups, state$pars, model$domain,
                               boundaries)
    
    if (!is.null(proposal$base)) {
      pars_next <- proposal$base(state$pars, rng)
      
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
      
      u <- rng$random_real(1)
      accept_prob_base <- pmin(1, exp(density_next - state$density))
      accept <- u < accept_prob_base
      
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
        ## TODO: mrc-5862
        if (model$properties$has_observer) {
          state$observation <- model$observer$observe()
        }
      }
    } else {
      accept_prob_base <- NULL
    }
    
    pars_next <- proposal$groups(state$pars, rng)
    
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
    
    n_groups <- max(model$parameter_groups)
    u <- rng$random_real(n_groups)
    accept_prob_groups <- 
      pmin(array(1, dim2(density_by_group_next)),
           exp(density_by_group_next - internal$density_by_group))
    accept_prob_groups[!i] <- 0
    accept <- u < accept_prob_groups
    
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
        ##new value of density (different to that which was accepted)
        density_next <- model$density(pars_next, by_group = TRUE)
        density_by_group_next <- attr(density_next, "by_group")
      }
      state$pars <- pars_next
      state$density <- c(density_next)
      internal$density_by_group <- density_by_group_next
      ## TODO: mrc-5862
      if (model$properties$has_observer) {
        state$observation <- model$observer$observe()
      }
    }

    if (internal$multiple_parameters) {
      internal$adaptive <-
        lapply(seq_len(dim(state$pars)[2]),
               function (i) update_nested_adaptive(internal$adaptive[[i]],
                                                   state$pars[, i],
                                                   model$parameter_groups,
                                                   accept_prob_base[i],
                                                   accept_prob_groups[, i]))
    } else {
      internal$adaptive <-
        update_nested_adaptive(internal$adaptive, state$pars, 
                               model$parameter_groups, accept_prob_base,
                               accept_prob_groups)
    }
    
    state
  }

  finalise <- function(state, model, rng) {
    out <- internal$adaptive
    
    keep_adaptive <- c("autocorrelation", "mean", "vcv", "weight", "included",
                       "scaling_history", "scaling_weight", "scaling_increment")
    
    if (internal$multiple_parameters) {
      out <- lapply(out, function(x) x[keep_adaptive])
    } else {
      out <- out[keep_adaptive]
    }
    
    out
  }

  get_internal_state <- function() {
    as.list(internal)
  }

  set_internal_state <- function(state) {
    list2env(state, internal)
  }

  monty_sampler("Nested adaptive Metropolis-Hastings",
                "monty_sampler_nested_adaptive",
                initialise,
                step,
                finalise,
                get_internal_state,
                set_internal_state)
}


check_nested_adaptive <- function(x, n_groups, has_base, null_allowed = FALSE,
                                  name = deparse(substitute(x))) {
  if (null_allowed) {
    if (!is.list(x) && length(x) != 1 && !is.null(x)) {
      cli::cli_abort(
        sprintf("Expected a list, single value or NULL for %s", name))
    }
  } else {
    if (!is.list(x) && length(x) != 1) {
      cli::cli_abort(sprintf("Expected a list or single value for %s", name))
    }
  }

  if (!is.list(x)) {
    if (has_base) {
      base <- x
    } else {
      base <- NULL
    }
    ret <- list(base = base,
                groups = rep(list(x), n_groups))
  } else {
    if (!setequal(names(x), c("base", "groups"))) {
      cli::cli_abort(
        sprintf("Expected %s input as list to have elements 'base' and
                'groups'", name))
    } else {
      if (!has_base && !is.null(x$base)) {
        cli::cli_abort(
          sprintf("Expected %s to be NULL as there are no base parameters",
                  paste0(name, "$base")))
      }
      if (null_allowed) {
        if (has_base && !is.null(x$base) && length(x$base) != 1) {
          cli::cli_abort(sprintf("Expected single value or NULL for %s",
                                 paste0(name, "$base")))
        }
        if (!is.list(x$groups) && length(x$groups) != 1 && !is.null(x$groups)) {
          cli::cli_abort(
            sprintf("Expected a list, single value or NULL for %s",
                    paste0(name, "$groups")))
        }
      } else {
        if (has_base && length(x$base) != 1) {
          cli::cli_abort(sprintf("Expected single value for %s",
                                 paste0(name, "$base")))
        }
        if (!is.list(x$groups) && length(x$groups) != 1) {
          cli::cli_abort(
            sprintf("Expected a list or single value for %s",
                    paste0(name, "$groups")))
        }
      }

      if (!is.list(x$groups)) {
        x$groups <- rep(list(x$groups), n_groups)
      } else {
        if (length(x$groups) != n_groups) {
          cli::cli_abort(
            sprintf("Expected %s specified as list to have length %s",
                    paste0(name, "$groups"), n_groups))
        }
        for (i in seq_len(n_groups)) {
          if (null_allowed) {
            if (length(x$groups[[i]]) != 1 && !is.null(x$groups[[i]])) {
              cli::cli_abort(
                sprintf("Expected a single value or NULL for %s",
                        paste0(name, "$groups[[", i, "]]")))
            }
          } else {
            if (length(x$groups[[i]]) != 1) {
              cli::cli_abort(
                sprintf("Expected a single value for %s",
                        paste0(name, "$groups[[", i, "]]")))
            }
          }
        }

      }
    }
    ret <- x
  }

  ret
}


split_nested_vcv <- function(vcv, pars) {
  n_groups <- length(vcv$groups)
  
  if (!is.null(vcv$base)) {
    vcv$base <- asplit(vcv$base, 3)
  } else {
    vcv$base <- rep(list(NULL), ncol(pars))
  }
  vcv$groups <- lapply(vcv$groups, function(x) asplit(x, 3))
  
  f <- function(i) {
    list(base = vcv$base[[i]],
         groups = lapply(seq_len(n_groups), function (j) vcv$groups[[j]][[i]])
    )
  }
  
  lapply(seq_len(ncol(pars)), f)
}


initialise_nested_adaptive <- function(pars,
                                       initial_vcv,
                                       initial_vcv_weight,
                                       initial_scaling,
                                       initial_scaling_weight,
                                       min_scaling,
                                       scaling_increment,
                                       log_scaling_update,
                                       acceptance_target,
                                       forget_rate,
                                       forget_end,
                                       adapt_end,
                                       pre_diminish,
                                       parameter_groups) {
  
  i_base <- parameter_groups == 0
  n_base <- sum(i_base)
  has_base <- n_base != 0
  n_groups <- max(parameter_groups)
  i_group <- lapply(seq_len(n_groups), function(i) which(parameter_groups == i))
  
  weight <- 0
  iteration <- 0
  
  initial_vcv_weight <-
    check_nested_adaptive(initial_vcv_weight, n_groups, has_base)
  scaling <-
    check_nested_adaptive(initial_scaling, n_groups, has_base)
  scaling_weight <-
    check_nested_adaptive(initial_scaling_weight, n_groups, has_base, TRUE)
  scaling_increment <-
    check_nested_adaptive(scaling_increment, n_groups, has_base, TRUE)
  min_scaling <-
    check_nested_adaptive(min_scaling, n_groups, has_base)
  
  
  if (!has_base) {
    mean_base <- NULL
    autocorrelation_base <- NULL
    vcv_base <- NULL
    proposal_vcv_base <- NULL
  } else {
    mean_base <- pars[i_base]
    autocorrelation_base <- matrix(0, n_base, n_base)
    vcv_base <- update_vcv(mean_base, autocorrelation_base, weight)
    proposal_vcv_base <-
      calc_proposal_vcv(scaling$base, vcv_base, weight,
                        initial_vcv$base, initial_vcv_weight$base)
    scaling_increment$base <- scaling_increment$base %||%
      calc_scaling_increment(n_base, acceptance_target, log_scaling_update)
    scaling_weight$base <- scaling_weight$base %||%
      5 / (acceptance_target * (1 - acceptance_target))
  }
  
  mean_groups <- lapply(i_group, function(i) pars[i])
  autocorrelation_groups <-
    lapply(lengths(i_group), function(x) matrix(0, x, x))
  vcv_groups <- Map(update_vcv, mean_groups, autocorrelation_groups,
                    weight)
  proposal_vcv_groups <-
    Map(calc_proposal_vcv, scaling$groups, vcv_groups,
        weight, initial_vcv$groups,
        initial_vcv_weight$groups)
  
  scaling_increment$groups <-
    Map(function(x, n) {
      x %||%
        calc_scaling_increment(n, acceptance_target, log_scaling_update)
    },
    scaling_increment$groups, lengths(i_group))
  scaling_weight$groups <-
    lapply(scaling_weight$groups,
           function(x) {
             x %||%  5 / (acceptance_target * (1 - acceptance_target))})
  
  mean <- list(base = mean_base, groups = mean_groups)
  autocorrelation <-
    list(base = autocorrelation_base, groups = autocorrelation_groups)
  vcv <- list(base = vcv_base, groups = vcv_groups)
  proposal_vcv <- list(base = proposal_vcv_base, groups = proposal_vcv_groups)
  
  history_pars <- numeric()
  included <- integer()
  scaling_history <- scaling
  
  list(initial_vcv = initial_vcv,
       initial_vcv_weight = initial_vcv_weight,
       weight = weight,
       iteration = iteration,
       mean = mean,
       autocorrelation = autocorrelation,
       vcv = vcv,
       proposal_vcv = proposal_vcv,
       scaling = scaling,
       scaling_increment = scaling_increment,
       scaling_weight = scaling_weight,
       min_scaling = min_scaling,
       log_scaling_update = log_scaling_update,
       acceptance_target = acceptance_target,
       forget_rate = forget_rate,
       forget_end = forget_end,
       adapt_end = adapt_end,
       pre_diminish = pre_diminish,
       history_pars = history_pars,
       included = included,
       scaling_history = scaling_history
  )
}


update_nested_adaptive <- function(adaptive, pars, parameter_groups,
                                   accept_prob_base, accept_prob_groups) {
  adaptive$iteration <- adaptive$iteration + 1
  adaptive$history_pars <- rbind(adaptive$history_pars, pars)
  if (adaptive$iteration > adaptive$adapt_end) {
    adaptive$scaling_history$base <-
      c(adaptive$scaling_history$base, adaptive$scaling$base)
    adaptive$scaling_history$groups <-
      Map(c, adaptive$scaling_history$groups, adaptive$scaling$groups)
    return(adaptive)
  }
  
  is_replacement <- check_replacement(adaptive$iteration, adaptive$forget_rate, 
                                      adaptive$forget_end)
  if (is_replacement) {
    pars_remove <- adaptive$history_pars[adaptive$included[1L], ]
    adaptive$included <- c(adaptive$included[-1L], adaptive$iteration)
  } else {
    pars_remove <- NULL
    adaptive$included <- c(adaptive$included, adaptive$iteration)
    adaptive$weight <- adaptive$weight + 1
  }
  
  if (!is.null(accept_prob_base)) {
    i_base <- parameter_groups == 0
    pars_base <- pars[i_base]
    pars_remove_base <- pars_remove[i_base]
    
    if (adaptive$iteration > adaptive$pre_diminish) {
      adaptive$scaling_weight$base <- adaptive$scaling_weight$base + 1
    }
    
    adaptive$scaling$base <-
      update_scaling(adaptive$scaling$base, adaptive$scaling_weight$base,
                     accept_prob_base, adaptive$scaling_increment$base,
                     adaptive$min_scaling$base, adaptive$acceptance_target,
                     adaptive$log_scaling_update)
    adaptive$scaling_history$base <-
      c(adaptive$scaling_history$base, adaptive$scaling$base)
    adaptive$autocorrelation$base <-
      update_autocorrelation(pars_base, adaptive$weight,
                             adaptive$autocorrelation$base, pars_remove_base)
    adaptive$mean$base <- update_mean(pars_base, adaptive$weight,
                                      adaptive$mean$base, pars_remove_base)
    adaptive$vcv$base <-
      update_vcv(adaptive$mean$base, adaptive$autocorrelation$base,
                 adaptive$weight)
    proposal_vcv_base <-
      calc_proposal_vcv(adaptive$scaling$base, adaptive$vcv$base,
                        adaptive$weight, adaptive$initial_vcv$base,
                        adaptive$initial_vcv_weight$base)
  } else {
    proposal_vcv_base <- NULL
  }
  
  n_groups <- max(parameter_groups)
  i_group <-
    lapply(seq_len(n_groups), function(i) which(parameter_groups == i))
  pars_groups <- lapply(i_group, function(i) pars[i])
  pars_remove_groups <- lapply(i_group, function(i) pars_remove[i])
  
  if (adaptive$iteration > adaptive$pre_diminish) {
    for (i in seq_len(n_groups)) {
      adaptive$scaling_weight$groups[[i]] <-
        adaptive$scaling_weight$groups[[i]] + 1
    }
  }
  
  adaptive$scaling$groups <-
    Map(update_scaling, adaptive$scaling$groups,
        adaptive$scaling_weight$groups, accept_prob_groups,
        adaptive$scaling_increment$groups, adaptive$min_scaling$groups,
        adaptive$acceptance_target, adaptive$log_scaling_update)
  adaptive$scaling_history$groups <-
    Map(c, adaptive$scaling_history$groups, adaptive$scaling$groups)
  adaptive$autocorrelation$groups <-
    Map(update_autocorrelation, pars_groups, adaptive$weight,
        adaptive$autocorrelation$groups, pars_remove_groups)
  adaptive$mean$groups <-
    Map(update_mean, pars_groups, adaptive$weight, adaptive$mean$groups,
        pars_remove_groups)
  adaptive$vcv$groups <-
    Map(update_vcv, adaptive$mean$groups, adaptive$autocorrelation$groups,
        adaptive$weight)
  proposal_vcv_groups <-
    Map(calc_proposal_vcv, adaptive$scaling$groups, adaptive$vcv$groups,
        adaptive$weight, adaptive$initial_vcv$groups,
        adaptive$initial_vcv_weight$groups)
  
  ## Update proposal
  adaptive$proposal_vcv <- list(base = proposal_vcv_base, 
                                groups = proposal_vcv_groups)
  
  adaptive
}


nested_proposal_adaptive <- function(adaptive, multiple_parameters,
                                     parameter_groups, pars, domain,
                                     boundaries = "reflect") {
  
  i_base <- parameter_groups == 0
  n_base <- sum(i_base)
  has_base <- n_base != 0
  n_groups <- max(parameter_groups)
  i_group <- lapply(seq_len(n_groups), function(i) which(parameter_groups == i))
  n_pars_by_group <- lengths(i_group)
  
  if (multiple_parameters) {
    if (has_base) {
      proposal_vcv_base <- vapply(adaptive, function (x) x$proposal_vcv$base,
                                  array(0, c(n_base, n_base)))
      proposal_vcv_base <- array(proposal_vcv_base,
                                 c(n_base, n_base, ncol(pars)))
    } else {
      proposal_vcv_base <- NULL
    }
    
    proposal_vcv_groups <- list()
    for (i in seq_len(n_groups)) {
      proposal_vcv_groups[[i]] <- 
        vapply(adaptive, function (x) x$proposal_vcv$groups[[i]], 
               array(0, c(n_pars_by_group[i], n_pars_by_group[i])))
      proposal_vcv_groups[[i]] <-
        array(proposal_vcv_groups[[i]],
              c(n_pars_by_group[i], n_pars_by_group[i], ncol(pars)))
    }
    
    proposal_vcv <- list(base = proposal_vcv_base,
                         groups = proposal_vcv_groups)
    proposal <- nested_proposal(proposal_vcv, parameter_groups, pars,
                                domain, boundaries)
  } else {
    proposal <- nested_proposal(adaptive$proposal_vcv, parameter_groups, pars,
                                domain, boundaries)
  }
  
  proposal
}
