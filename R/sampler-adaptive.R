##' Create an adaptive Metropolis-Hastings sampler, which will tune
##' its variance covariance matrix (vs the simple random walk
##' sampler [mcstate_sampler_random_walk]).
##'
##' Efficient exploration of the parameter space during an MCMC might
##' be difficult when the target distribution is of high
##' dimensionality, especially if the target probability distribution
##' present a high degree of correlation.  Adaptive schemes are used
##' to "learn" on the fly the correlation structure by updating the
##' proposal distribution by recalculating the empirical
##' variance-covariance matrix and rescale it at each adaptive step of
##' the MCMC.
##'
##' Our implementation of an adaptive MCMC algorithm is based on an
##' adaptation of the "accelerated shaping" algorithm in Spencer (2021).
##' The algorithm is based on a random-walk Metropolis-Hastings algorithm where
##' the proposal is a multi-variate Normal distribution centred on the current
##' point.
##'
##' Spencer SEF (2021) Accelerating adaptation in the adaptive
##' Metropolis–Hastings random walk algorithm. Australian & New Zealand Journal
##' of Statistics 63:468-484.
##'
##' @title Adaptive Metropolis-Hastings Sampler
##'
##' @param initial_vcv An initial variance covariance matrix; we'll start
##'   using this in the proposal, which will gradually become more weighted
##'   towards the empirical covariance matrix calculated from the chain.
##'
##' @param initial_vcv_weight Weight of the initial variance-covariance
##'   matrix used to build the proposal of the random-walk. Higher
##'   values translate into higher confidence of the initial
##'   variance-covariance matrix and means that update from additional
##'   samples will be slower.
##'
##' @param initial_scaling The initial scaling of the variance
##'   covariance matrix to be used to generate the multivariate normal
##'   proposal for the random-walk Metropolis-Hastings algorithm. To generate
##'   the proposal matrix, the weighted variance covariance matrix is
##'   multiplied by the scaling parameter squared times 2.38^2 / n_pars (where
##'   n_pars is the number of fitted parameters). Thus, in a Gaussian target
##'   parameter space, the optimal scaling will be around 1.
##'
##' @param initial_scaling_weight The initial weight used in the scaling update.
##'   The scaling weight will increase after the first `pre_diminish`
##'   iterations, and as the scaling weight increases the adaptation of the
##'   scaling diminishes. If `NULL` (the default) the value is
##'   5 / (acceptance_target * (1 - acceptance_target)).
##'
##' @param min_scaling The minimum scaling of the variance covariance
##'   matrix to be used to generate the multivariate normal proposal
##'   for the random-walk Metropolis-Hastings algorithm.
##'
##' @param scaling_increment The scaling increment which is added or
##'   subtracted to the scaling factor of the variance-covariance
##'   after each adaptive step. If `NULL` (the default) then an optimal
##'   value will be calculated.
##'
##' @param log_scaling_update Logical, whether or not changes to the
##'   scaling parameter are made on the log-scale.
##'
##' @param acceptance_target The target for the fraction of proposals
##'   that should be accepted (optimally) for the adaptive part of the
##'   mixture model.
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
mcstate_sampler_adaptive <- function(initial_vcv,
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
  ## This sampler is stateful; we will be updating our estimate of the
  ## mean and vcv of the target distribution, along with the our
  ## scaling factor, weight and autocorrelations.
  check_vcv(initial_vcv, allow_3d = TRUE, call = environment())
  internal <- new.env()
  
  boundaries <- match_value(boundaries, c("reflect", "reject", "ignore"))

  initialise <- function(pars, model, observer, rng) {
    require_deterministic(model,
                          "Can't use adaptive sampler with stochastic models")

    internal$multiple_parameters <- length(dim2(pars)) > 1
    if (internal$multiple_parameters) {
      ## this is enforced elsewhere
      stopifnot(model$properties$allow_multiple_parameters)
    }
    
    initial_vcv <- sampler_validate_vcv(initial_vcv, pars)
    
    if (internal$multiple_parameters) {
      d <- dim(initial_vcv)
      internal$adaptive <- 
        lapply(seq_len(d[3]),
               function(i) initialise_adaptive(
                 pars[, i],
                 array(initial_vcv[, , i], d[c(1, 2)]),
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
                 pre_diminish
               ))
    } else {
      internal$adaptive <- 
        initialise_adaptive(pars,
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
                            pre_diminish)
    }

    initialise_state(pars, model, observer, rng)
  }

  step <- function(state, model, observer, rng) {
    if (internal$multiple_parameters) {
      d <- dim(state$pars)
      proposal_vcv <- 
        vapply(seq_len(d[2]),
               function (i) 
                 calc_proposal_vcv(internal$adaptive[[i]]$scaling, 
                                   internal$adaptive[[i]]$vcv,
                                   internal$adaptive[[i]]$weight, 
                                   internal$adaptive[[i]]$initial_vcv,
                                   internal$adaptive[[i]]$initial_vcv_weight),
               array(0, c(d[1], d[1])))
      proposal_vcv <- array(proposal_vcv, c(d[1], d[1], d[2]))
    } else {
      proposal_vcv <-
        calc_proposal_vcv(internal$adaptive$scaling,
                          internal$adaptive$vcv,
                          internal$adaptive$weight, 
                          internal$adaptive$initial_vcv,
                          internal$adaptive$initial_vcv_weight)
    }

    proposal <- 
      make_random_walk_proposal(proposal_vcv, model$domain, boundaries)
    pars_next <- proposal(state$pars, rng)

    u <- rng$random_real(1)
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

    accept_prob <- pmin(1, exp(density_next - state$density))

    accept <- u < accept_prob
    state <- update_state(state, pars_next, density_next, accept,
                          model, observer, rng)

    if (internal$multiple_parameters) {
      internal$adaptive <- 
        lapply(seq_len(dim(state$pars)[2]),
               function (i) update_adaptive(internal$adaptive[[i]],
                                            state$pars[, i],
                                            accept_prob[i]))
    } else {
      internal$adaptive <-
        update_adaptive(internal$adaptive, state$pars, accept_prob)
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

  mcstate_sampler("Adaptive Metropolis-Hastings",
                  initialise,
                  step,
                  finalise,
                  get_internal_state,
                  set_internal_state)
}


initialise_adaptive <- function(pars,
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
                                pre_diminish) {
  weight <- 0
  iteration <- 0

  mean <- unname(pars)
  n_pars <- length(pars)
  autocorrelation <- array(0, dim(initial_vcv))
  vcv <- update_vcv(mean, autocorrelation, weight)

  scaling <- initial_scaling

  scaling_increment <- scaling_increment %||%
    calc_scaling_increment(n_pars, acceptance_target, log_scaling_update)
  scaling_weight <- initial_scaling_weight %||%
    5 / (acceptance_target * (1 - acceptance_target))

  history_pars <- NULL
  included <- integer()
  scaling_history <- scaling

  list(initial_vcv = initial_vcv,
       initial_vcv_weight = initial_vcv_weight,
       weight = weight,
       iteration = iteration,
       mean = mean,
       autocorrelation = autocorrelation,
       vcv = vcv,
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

update_adaptive <- function(adaptive, pars, accept_prob) {
  adaptive$iteration <- adaptive$iteration + 1
  adaptive$history_pars <- rbind(adaptive$history_pars, pars)
  if (adaptive$iteration > adaptive$adapt_end) {
    adaptive$scaling_history <- c(adaptive$scaling_history, adaptive$scaling)
    return(adaptive)
  }

  if (adaptive$iteration > adaptive$pre_diminish) {
    adaptive$scaling_weight <- adaptive$scaling_weight + 1
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

  adaptive$scaling <-
    update_scaling(adaptive$scaling, adaptive$scaling_weight, accept_prob,
                   adaptive$scaling_increment, adaptive$min_scaling,
                   adaptive$acceptance_target, adaptive$log_scaling_update)
  adaptive$scaling_history <- c(adaptive$scaling_history, adaptive$scaling)
  adaptive$autocorrelation <- update_autocorrelation(
    pars, adaptive$weight, adaptive$autocorrelation, pars_remove)
  adaptive$mean <- update_mean(pars, adaptive$weight, adaptive$mean,
                               pars_remove)
  adaptive$vcv <- update_vcv(adaptive$mean, adaptive$autocorrelation,
                             adaptive$weight)

  adaptive
}


calc_scaling_increment <- function(n_pars, acceptance_target,
                                   log_scaling_update) {
  if (log_scaling_update) {
    A <- -stats::qnorm(acceptance_target / 2)

    scaling_increment <-
      (1 - 1 / n_pars) * (sqrt(2 * pi) * exp(A^2 / 2)) / (2 * A) +
      1 / (n_pars * acceptance_target * (1 - acceptance_target))
  } else {
    scaling_increment <- 1 / 100
  }

  scaling_increment
}


qp <- function(x) {
  outer(x, x)
}


calc_proposal_vcv <- function(scaling, vcv, weight, initial_vcv,
                              initial_vcv_weight) {
  n_pars <- nrow(vcv)

  weighted_vcv <-
    ((weight - 1) * vcv + (initial_vcv_weight + n_pars + 1) * initial_vcv) /
    (weight + initial_vcv_weight + n_pars + 1)

  2.38^2 / n_pars * scaling^2 * weighted_vcv
}


check_replacement <- function(iteration, forget_rate, forget_end) {
  is_forget_step <- floor(forget_rate * iteration) >
    floor(forget_rate * (iteration - 1))
  is_before_forget_end <- iteration <= forget_end

  is_forget_step & is_before_forget_end
}


update_scaling <- function(scaling, scaling_weight, accept_prob,
                           scaling_increment, min_scaling,
                           acceptance_target, log_scaling_update) {
  scaling_change <- scaling_increment * (accept_prob - acceptance_target) /
    sqrt(scaling_weight)

  if (log_scaling_update) {
    max(min_scaling, scaling * exp(scaling_change))
  } else {
    max(min_scaling, scaling + scaling_change)
  }

}


update_autocorrelation <- function(pars, weight, autocorrelation, pars_remove) {
  if (!is.null(pars_remove)) {
    if (weight > 2) {
      autocorrelation <-
        autocorrelation + 1 / (weight - 1) * (qp(pars) - qp(pars_remove))
    } else {
      autocorrelation <- autocorrelation + qp(pars) - qp(pars_remove)
    }
  } else {
    if (weight > 2) {
      autocorrelation <-
        (1 - 1 / (weight - 1)) * autocorrelation + 1 / (weight - 1) * qp(pars)
    } else {
      autocorrelation <- autocorrelation + qp(pars)
    }
  }

  autocorrelation
}


update_mean <- function(pars, weight, mean, pars_remove) {
  if (!is.null(pars_remove)) {
    mean <- mean + 1 / weight * (pars - pars_remove)
  } else {
    mean <- (1 - 1 / weight) * mean + 1 / weight * pars
  }

  mean
}


update_vcv <- function(mean, autocorrelation, weight) {
  if (weight > 1) {
    vcv <- autocorrelation - weight / (weight - 1) * qp(mean)
  } else {
    vcv <- 0 * autocorrelation
  }

  vcv
}
