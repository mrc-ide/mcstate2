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
                                     pre_diminish = 0) {
  ## This sampler is stateful; we will be updating our estimate of the
  ## mean and vcv of the target distribution, along with the our
  ## scaling factor, weight and autocorrelations.
  internal <- new.env()

  initialise <- function(pars, model, observer, rng) {
    require_deterministic(model,
                          "Can't use adaptive sampler with stochastic models")
    internal$weight <- 0
    internal$iteration <- 0

    internal$mean <- unname(pars)
    n_pars <- length(model$parameters)
    internal$autocorrelation <- matrix(0, n_pars, n_pars)
    internal$vcv <- update_vcv(internal$mean, internal$autocorrelation,
                               internal$weight)

    internal$scaling <- initial_scaling
    internal$scaling_increment <- scaling_increment %||%
      calc_scaling_increment(n_pars, acceptance_target,
                             log_scaling_update)
    internal$scaling_weight <- initial_scaling_weight %||%
      5 / (acceptance_target * (1 - acceptance_target))

    internal$history_pars <- numeric()
    internal$included <- integer()
    internal$scaling_history <- internal$scaling

    initialise_state(pars, model, observer, rng)
  }

  step <- function(state, model, observer, rng) {
    proposal_vcv <-
      calc_proposal_vcv(internal$scaling, internal$vcv, internal$weight,
                        initial_vcv, initial_vcv_weight)

    pars_next <- rmvnorm(state$pars, proposal_vcv, rng)

    u <- rng$random_real(1)
    density_next <- model$density(pars_next)

    accept_prob <- min(1, exp(density_next - state$density))

    accept <- u < accept_prob
    if (accept) {
      state <- update_state(state, pars_next, density_next, model, observer,
                            rng)
    }

    internal$iteration <- internal$iteration + 1
    internal$history_pars <- rbind(internal$history_pars, state$pars)
    if (internal$iteration > adapt_end) {
      internal$scaling_history <- c(internal$scaling_history, internal$scaling)
      return(state)
    }

    if (internal$iteration > pre_diminish) {
      internal$scaling_weight <- internal$scaling_weight + 1
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

    internal$scaling <-
      update_scaling(internal$scaling, internal$scaling_weight, accept_prob,
                     internal$scaling_increment, min_scaling, acceptance_target,
                     log_scaling_update)
    internal$scaling_history <- c(internal$scaling_history, internal$scaling)
    internal$autocorrelation <- update_autocorrelation(
      state$pars, internal$weight, internal$autocorrelation, pars_remove)
    internal$mean <- update_mean(state$pars, internal$weight, internal$mean,
                                 pars_remove)
    internal$vcv <- update_vcv(internal$mean, internal$autocorrelation,
                               internal$weight)

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

  mcstate_sampler("Adaptive Metropolis-Hastings",
                  initialise,
                  step,
                  finalise,
                  get_internal_state,
                  set_internal_state)
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
