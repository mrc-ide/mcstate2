##' Create a nested adaptive Metropolis-Hastings sampler, which extends the 
##' adaptive sampler [mcstate_sampler_adaptive], tuning the variance covariance
##' matrices for proposal for the separable sections
##' of a nested model (vs the simple nested random walk sampler
##' [mcstate_sampler_random_walk]). This sampler requires
##' that models support the `has_parameter_groups` property.
##'
##' Much like the simple nested random walk sampler
##' [mcstate_sampler_random_walk], the strategy is to propose all the 
##' shared parameters as a deviation from the current point in parameter space
##' as a single move and accept or reject as a block. Then we generate points
##' for all the region-specific parameters, compute the density and then
##' accept or reject these updates independently.  This is possible
##' because the change in likelihood in region A is independent from
##' region B.
##'
##' The adaptive proposal algorithm of the non-nested adaptive sampler 
##' [mcstate_sampler_adaptive] is extended here to adaptively tune the variance
##' covariance matrix of each of these parameter chunks.
##'
##' @title Nested Adaptive Metropolis-Hastings Sampler
##'
##' @param initial_vcv A list of initial variance covariance matrices. We 
##'   expect this to be a list with elements `base` and `groups` corresponding
##'   to the initial covariance matrix for base parameters (if any) and groups.
##'   We'll start using these in the proposal, which will gradually become more
##'   weighted towards the empirical covariance matrices calculated from the
##'   chain.
##'
##' @param initial_vcv_weight Weight of the initial variance-covariance
##'   matrices used to build the proposal of the random-walk. Higher
##'   values translate into higher confidence of the initial
##'   variance-covariance matrix and means that update from additional
##'   samples will be slower. This can be either (a) a list with elements `base`
##'   and `groups`, where the `groups` element can in turn be a list of values
##'   for each group or a single value applying to all groups or (b) a single
##'   value which applies to the base parameters and all groups.
##'
##' @param initial_scaling The initial scalings of the variance
##'   covariance matrix to be used to generate the multivariate normal
##'   proposal for the random-walk Metropolis-Hastings algorithm. To generate
##'   the proposal matrix, the weighted variance covariance matrix is
##'   multiplied by the scaling parameter squared times 2.38^2 / n_pars (where
##'   n_pars is the number of fitted parameters in the parameter chunk). Can be
##'   either (a) a list with elements `base` and `groups`, where the `groups`
##'   element can in turn be a list of values for each group or a single value
##'   applying to all groups or (b) a single value which applies to the base
##'   parameters and all groups.
##'
##' @param initial_scaling_weight The initial weights used in the scaling
##'   update. The scaling weight will increase after the first `pre_diminish`
##'   iterations, and as the scaling weight increases the adaptation of the
##'   scaling diminishes. Can be either (a) a list with elements
##'   `base` and `groups`, where the `groups` element can in turn be a list of
##'   values for each group or a single value applying to all groups or (b) a 
##'   single value which applies to the base parameters and all groups. If
##'   `NULL` (the default) the value is 
##'   5 / (acceptance_target * (1 - acceptance_target)).
##'
##' @param min_scaling The minimum scalings of the variance covariance
##'   matrix to be used to generate the multivariate normal proposal
##'   for the random-walk Metropolis-Hastings algorithm. Can be either (a) a
##'   list with elements `base` and `groups`, where the `groups` element can in
##'   turn be a list of values for each group or a single value applying to all
##'   groups or (b) a single value which applies to the base parameters and all
##'   groups.
##'
##' @param scaling_increment The scaling increments which are added or
##'   subtracted to the scaling factor of the variance-covariance
##'   after each adaptive step. Can be either (a) a list with elements
##'   `base` and `groups`, where the `groups` element can in turn be a list of
##'   values for each group or a single value applying to all groups or (b) a 
##'   single value which applies to the base parameters and all groups. If
##'   `NULL` (the default) then an optimal value will be calculated.
##'
##' @param log_scaling_update Logical, whether or not changes to the
##'   scaling parameters are made on the log-scale.
##'
##' @param acceptance_target The target for the fraction of proposals
##'   that should be accepted (optimally) for the adaptive part of the
##'   chain.
##'
##' @param forget_rate The rate of forgetting early parameter sets from the
##'   empirical variance-covariance matrix in the MCMC chains. For example,
##'   `forget_rate = 0.2` (the default) means that once in every 5th iterations
##'   we remove the earliest parameter set included, so would remove the 1st
##'   parameter set on the 5th update, the 2nd on the 10th update, and so
##'   on. Setting `forget_rate = 0` means early parameter sets are never
##'   forgotten.
##'
##' @param forget_end The final iteration at which early parameter sets can
##'   be forgotten. Setting `forget_rate = Inf` (the default) means that the
##'   forgetting mechanism continues throughout the chains. Forgetting early
##'   parameter sets becomes less useful once the chains have settled into the
##'   posterior mode, so this parameter might be set as an estimate of how long
##'   that would take.
##'
##' @param adapt_end The final iteration at which we can adapt the multivariate
##'   normal proposal. Thereafter the empirical variance-covariance matrix, its
##'   scaling and its weight remain fixed. This allows the adaptation to be
##'   switched off at a certain point to help ensure convergence of the chain.
##'
##' @param pre_diminish The number of updates before adaptation of the scaling
##'   parameter starts to diminish. Setting `pre_diminish = 0` means there is
##'   diminishing adaptation of the scaling parameter from the offset, while
##'   `pre_diminish = Inf` would mean there is never diminishing adaptation.
##'   Diminishing adaptation should help the scaling parameter to converge
##'   better, but while the chains find the location and scale of the posterior
##'   mode it might be useful to explore with it switched off.
##'
##' @return A `mcstate_sampler` object, which can be used with
##'   [mcstate_sample]
##'
##' @export
mcstate_sampler_nested_adaptive <- function(initial_vcv,
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
                                            pre_diminish = 0) {
  if (!is.list(initial_vcv)) {
    cli::cli_abort(
      "Expected a list for 'initial_vcv'",
      arg = 'initial_vcv')
  }

  if (!setequal(names(initial_vcv), c("base", "groups"))) {
    cli::cli_abort(
      "Expected 'initial_vcv' to have elements 'base' and 'groups'",
      arg = "initial_vcv")
  }
  if (!is.null(initial_vcv$base)) {
    check_vcv(initial_vcv$base, call = environment())
  }
  if (!is.list(initial_vcv$groups)) {
    cli::cli_abort("Expected 'initial_vcv$groups' to be a list")
  }
  if (length(initial_vcv$groups) < 1) {
    cli::cli_abort(
      "Expected 'initial_vcv$groups' to have at least one element")
  }
  for (i in seq_along(initial_vcv$groups)) {
    check_vcv(initial_vcv$groups[[i]],
              name = sprintf("initial_vcv$groups[%d]", i),
              call = environment())
  }

  internal <- new.env(parent = emptyenv())

  initialise <- function(pars, model, observer, rng) {
    require_deterministic(model,
                          "Can't use adaptive sampler with stochastic models")
    
    if (!model$properties$has_parameter_groups) {
      cli::cli_abort("Your model does not have parameter groupings")
    }

    initialise_rng_state(model, rng)
    density <- model$density(pars, by_group = TRUE)
    density_by_group <- attr(density, "by_group")
    i_base <- model$parameter_groups == 0
    n_base <- sum(i_base)
    has_base <- n_base != 0
    n_groups <- max(model$parameter_groups)
    i_group <- 
      lapply(seq_len(n_groups), function(i) which(model$parameter_groups == i))

    if (is.null(density_by_group)) {
      cli::cli_abort(
        c(paste("model$density(x, by_group = TRUE) did not produce a",
                "density with a 'by_group' attribute"),
          i = paste("I expected an attribute 'by_group' with {n_groups}",
                    "elements corresponding to parameter groups to be",
                    "included with your density")))
    }
    if (length(density_by_group) != n_groups) {
        cli::cli_abort(
          paste("model$density(x, by_group = TRUE) produced a 'by_group'",
                "attribute with incorrect length {length(density_by_group)}",
                "but I expected length {n_groups}"))
    }

    internal$density_by_group <- density_by_group
    
    internal$weight <- 0
    internal$iteration <- 0
    
    internal$initial_vcv_weight <- 
      check_nested_adaptive(initial_vcv_weight, n_groups, has_base)
    internal$scaling <- 
      check_nested_adaptive(initial_scaling, n_groups, has_base)
    internal$scaling_weight <- 
      check_nested_adaptive(initial_scaling_weight, n_groups, has_base, TRUE)
    internal$scaling_increment <-
      check_nested_adaptive(scaling_increment, n_groups, has_base, TRUE)
    internal$min_scaling <-
      check_nested_adaptive(min_scaling, n_groups, has_base)
    
    
    if (!has_base) {
      mean_base <- NULL
      autocorrelation_base <- NULL
      vcv_base <- NULL
      proposal_vcv_base <- NULL
    } else {
      mean_base <- pars[i_base]
      autocorrelation_base <- matrix(0, n_base, n_base)
      vcv_base <- update_vcv(mean_base, autocorrelation_base, internal$weight)
      proposal_vcv_base <- 
        calc_proposal_vcv(internal$scaling$base, vcv_base, internal$weight,
                          initial_vcv$base, internal$initial_vcv_weight$base)
      internal$scaling_increment$base <- internal$scaling_increment$base %||%
        calc_scaling_increment(n_base, acceptance_target, log_scaling_update)
      internal$scaling_weight$base <- internal$scaling_weight$base %||%
        5 / (acceptance_target * (1 - acceptance_target))
    }
    
    mean_groups <- lapply(i_group, function(i) pars[i])
    autocorrelation_groups <-
      lapply(lengths(i_group), function (x) matrix(0, x, x))
    vcv_groups <- Map(update_vcv, mean_groups, autocorrelation_groups,
                      internal$weight)
    proposal_vcv_groups <- 
      Map(calc_proposal_vcv, internal$scaling$groups, vcv_groups,
          internal$weight, initial_vcv$groups,
          internal$initial_vcv_weight$groups)
    
    internal$scaling_increment$groups <- 
      Map(function(x, n) {x %||% 
          calc_scaling_increment(n, acceptance_target, log_scaling_update)},
          internal$scaling_increment$groups, lengths(i_group))
    internal$scaling_weight$groups <- 
      lapply(internal$scaling_weight$groups,
             function(x) {
               x %||%  5 / (acceptance_target * (1 - acceptance_target))})
    
    internal$mean <- list(base = mean_base, groups = mean_groups)
    internal$autocorrelation <- 
      list(base = autocorrelation_base, groups = autocorrelation_groups)
    internal$vcv <- list(base = vcv_base, groups = vcv_groups)
    proposal_vcv <- list(base = proposal_vcv_base, groups = proposal_vcv_groups)
    internal$proposal <- nested_proposal(proposal_vcv, model$parameter_groups)
    
    internal$history_pars <- numeric()
    internal$included <- integer()
    internal$scaling_history <- internal$scaling
    
    state <- list(pars = pars, density = c(density))
    if (!is.null(observer)) {
      state$observation <- observer$observe(model$model, rng)
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
  step <- function(state, model, observer, rng) {
    if (!is.null(internal$proposal$base)) {
      pars_next <- internal$proposal$base(state$pars, rng)
      density_next <- model$density(pars_next, by_group = TRUE)
      density_by_group_next <- attr(density_next, "by_group")
      u <- rng$random_real(1)
      accept_prob_base <- min(1, exp(density_next - state$density))
      accept <- u < accept_prob_base
      if (accept) {
        state$pars <- pars_next
        state$density <- density_next
        internal$density_by_group <- density_by_group_next
        if (!is.null(observer)) {
          state$observation <- observer$observe(model$model, rng)
        }
      }
    }

    pars_next <- internal$proposal$groups(state$pars, rng)
    density_next <- model$density(pars_next, by_group = TRUE)
    density_by_group_next <- attr(density_next, "by_group")
    u <- rng$random_real(length(density_by_group_next))
    accept_prob_groups <- 
      pmin(1, exp(density_by_group_next - internal$density_by_group))
    accept <- u < accept_prob_groups
    if (any(accept)) {
      if (!all(accept)) {
        ## Retain some older parameters
        i <- model$parameter_groups %in% which(!accept)
        pars_next[i] <- state$pars[i]
        density_next <- model$density(pars_next, by_group = TRUE)
        density_by_group_next <- attr(density_next, "by_group")
      }
      state$pars <- pars_next
      state$density <- c(density_next)
      internal$density_by_group <- density_by_group_next
      if (!is.null(observer)) {
        state$observation <- observer$observe(model$model, rng)
      }
    }
    
    internal$iteration <- internal$iteration + 1
    internal$history_pars <- rbind(internal$history_pars, state$pars)
    if (internal$iteration > adapt_end) {
      internal$scaling_history$base <- 
        c(internal$scaling_history$base, internal$scaling$base)
      internal$scaling_history$groups <- 
        Map(c, internal$scaling_history$groups, internal$scaling$groups)
      return(state)
    }
    
    is_replacement <-
      check_replacement(internal$iteration, forget_rate, forget_end)
    if (is_replacement) {
      pars_remove <- internal$history_pars[internal$included[1L], ]
      internal$included <- c(internal$included[-1L], internal$iteration)
    } else {
      pars_remove <- NULL
      internal$included <- c(internal$included, internal$iteration)
      internal$weight <- internal$weight + 1
    }
    
    if (!is.null(internal$proposal$base)) {
      i_base <- model$parameter_groups == 0
      pars_base <- state$pars[i_base]
      pars_remove_base <- pars_remove[i_base]
      
      if (internal$iteration > pre_diminish) {
        internal$scaling_weight$base <- internal$scaling_weight$base + 1
      }
      
      internal$scaling$base <-
        update_scaling(internal$scaling$base, internal$scaling_weight$base, 
                       accept_prob_base, internal$scaling_increment$base,
                       internal$min_scaling$base, acceptance_target,
                       log_scaling_update)
      internal$scaling_history$base <- 
        c(internal$scaling_history$base, internal$scaling$base)
      internal$autocorrelation$base <- 
        update_autocorrelation(pars_base, internal$weight, 
                               internal$autocorrelation$base, pars_remove_base)
      internal$mean$base <- update_mean(pars_base, internal$weight,
                                        internal$mean$base, pars_remove_base)
      internal$vcv$base <- 
        update_vcv(internal$mean$base, internal$autocorrelation$base,
                   internal$weight)
      proposal_vcv_base <- 
        calc_proposal_vcv(internal$scaling$base, internal$vcv$base,
                          internal$weight, initial_vcv$base, 
                          internal$initial_vcv_weight$base)
    } else {
      proposal_vcv_base <- NULL
    }
    
    n_groups <- max(model$parameter_groups)
    i_group <- 
      lapply(seq_len(n_groups), function(i) which(model$parameter_groups == i))
    pars_groups <- lapply(i_group, function(i) state$pars[i])
    pars_remove_groups <- lapply(i_group, function(i) pars_remove[i])
    
    if (internal$iteration > pre_diminish) {
      for (i in seq_len(n_groups)) {
        internal$scaling_weight$groups[[i]] <- 
          internal$scaling_weight$groups[[i]] + 1
      }
    }
    
    internal$scaling$groups <-
      Map(update_scaling, internal$scaling$groups, 
          internal$scaling_weight$groups, accept_prob_groups,
          internal$scaling_increment$groups, internal$min_scaling$groups,
          acceptance_target, log_scaling_update)
    internal$scaling_history$groups <- 
      Map(c, internal$scaling_history$groups, internal$scaling$groups)
    internal$autocorrelation$groups <- 
      Map(update_autocorrelation, pars_groups, internal$weight, 
          internal$autocorrelation$groups, pars_remove_groups)
    internal$mean$groups <- 
      Map(update_mean, pars_groups, internal$weight, internal$mean$groups,
          pars_remove_groups)
    internal$vcv$groups <- 
      Map(update_vcv, internal$mean$groups, internal$autocorrelation$groups, 
          internal$weight)
    proposal_vcv_groups <- 
      Map(calc_proposal_vcv, internal$scaling$groups, internal$vcv$groups,
          internal$weight, initial_vcv$groups, 
          internal$initial_vcv_weight$groups)
    
    ## Update proposal
    proposal_vcv <- list(base = proposal_vcv_base, groups = proposal_vcv_groups)
    internal$proposal <- nested_proposal(proposal_vcv, model$parameter_groups)
    
    state
  }

  finalise <- function(state, model, rng) {
    out <- as.list(internal)
    out[c("autocorrelation", "mean", "vcv", "weight", "included",
          "scaling_history", "scaling_weight", "scaling_increment")]
  }

  get_internal_state <- function() {
    as.list(internal)
  }

  set_internal_state <- function(state) {
    list2env(state, internal)
  }

  mcstate_sampler("Nested adaptive Metropolis-Hastings",
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
