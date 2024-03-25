ex_simple_gamma1 <- function(shape = 1, rate = 1) {
  e <- new.env(parent = .GlobalEnv)
  e$shape <- shape
  e$rate <- rate
  with(
    e,
    mcstate_model(list(
      parameters = "gamma",
      direct_sample = function(rng) {
        rng$gamma(1, shape = shape, scale = 1 / rate)
      },
      density = function(x) dgamma(x, shape = shape, rate = rate, log = TRUE),
      gradient = function(x) (shape - 1) / x - rate,
      domain = rbind(c(0, Inf)))))
}


ex_dust_sir <- function(n_particles = 100, n_threads = 1,
                        deterministic = FALSE, save_trajectories = FALSE) {
  testthat::skip_if_not_installed("dust")
  sir <- dust::dust_example("sir")

  np <- 10
  end <- 150 * 4
  times <- seq(0, end, by = 4)
  ans <- sir$new(list(), 0, np, seed = 1L)$simulate(times)
  dat <- data.frame(time = times[-1], incidence = ans[5, 1, -1])

  ## TODO: an upshot here is that our dust models are always going to
  ## need to be initialisable; we might need to sample from the
  ## statistical parameters, or set things up to allow two-phases of
  ## initialsation (which is I think where we are heading, so that's
  ## fine).
  model <- sir$new(list(), 0, n_particles, seed = 1L, n_threads = n_threads,
                   deterministic = deterministic)
  model$set_data(dust::dust_data(dat))
  model$set_index(c(2, 4))

  prior_beta_shape <- 1
  prior_beta_rate <- 1 / 0.5
  prior_gamma_shape <- 1
  prior_gamma_rate <- 1 / 0.5

  trajectories <- NULL

  ## In the new dust wrapper we'll need to make this nicer; I think
  ## that this is pretty painful atm because we wrap via the particle
  ## filter method in mcstate1.  This version replicates most of what
  ## we need though, which is some subset of the model
  details <- function(i, j) {
    if (save_trajectories) {
      traj <- trajectories[, i, , drop = FALSE]
      dim(traj) <- dim(traj)[-2]
    } else {
      traj <- NULL
    }
    list(trajectories = traj, state = model$state()[, i])
  }

  density <- function(x) {
    beta <- x[[1]]
    gamma <- x[[2]]
    prior <- dgamma(beta, prior_beta_shape, prior_beta_rate, log = TRUE) +
      dgamma(gamma, prior_gamma_shape, prior_gamma_rate, log = TRUE)
    if (is.finite(prior)) {
      model$update_state(
        pars = list(beta = x[[1]], gamma = x[[2]]),
        time = 0,
        set_initial_state = TRUE)
      res <- model$filter(save_trajectories = save_trajectories)
      if (save_trajectories) {
        trajectories <<- res$trajectories
      }
      ll <- res$log_likelihood
    } else {
      ll <- -Inf
    }
    ll + prior
  }

  direct_sample <- function(rng) {
    c(rng$gamma(1, prior_beta_shape, 1 / prior_beta_rate),
      rng$gamma(1, prior_gamma_shape, 1 / prior_gamma_rate))
  }

  set_rng_state <- function(rng_state) {
    n_streams <- n_particles + 1
    if (length(rng_state) != 32 * n_streams) {
      ## Expand the state by short jumps; we'll make this nicer once
      ## we refactor the RNG interface and dust.
      rng_state <- mcstate_rng$new(rng_state, n_streams)$state()
    }
    model$set_rng_state(rng_state)
  }

  get_rng_state <- function() {
    model$rng_state()
  }

  mcstate_model(
    list(model = model,
         details = details,
         density = density,
         direct_sample = direct_sample,
         parameters = c("beta", "gamma"),
         domain = cbind(c(0, 0), c(Inf, Inf)),
         set_rng_state = set_rng_state,
         get_rng_state = get_rng_state),
    mcstate_model_properties(is_stochastic = !deterministic))
}


ex_simple_gaussian <- function(vcv) {
  n <- nrow(vcv)
  mcstate_model(list(
    parameters = letters[seq_len(n)],
    direct_sample = make_rmvnorm(vcv, centred = TRUE),
    density = make_ldmvnorm(vcv),
    gradient = make_deriv_ldmvnorm(vcv),
    domain = cbind(rep(-Inf, n), rep(Inf, n))))
}


ex_banana <- function(sd = 0.5) {
  mcstate_model(list(
    parameters = c("theta1", "theta2"),
    direct_sample = function(rng) {
      theta <- rng$random_normal(1)
      cbind(rng$normal(1, theta^2, sd), theta)
    },
    density = function(x) {
      dnorm(x[2], log = TRUE) + dnorm((x[1] - x[2]^2) / sd, log = TRUE)
    },
    gradient = function(x){
      c((x[2]^2 - x[1])/sd^2,
        -x[2] + 2 * x[2] * (x[1] - x[2]^2) / sd^2)
    },
    domain = cbind(rep(-Inf, 2), rep(Inf, 2))))
}


random_array <- function(dim, named = FALSE) {
  if (named) {
    dn <- lapply(seq_along(dim), function(i)
      paste0(LETTERS[[i]], letters[seq_len(dim[i])]))
    names(dn) <- paste0("d", LETTERS[seq_along(dim)])
  } else {
    dn <- NULL
  }
  array(runif(prod(dim)), dim, dimnames = dn)
}
