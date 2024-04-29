ex_simple_gamma1 <- function(shape = 1, rate = 1) {
  e <- new.env(parent = .GlobalEnv)
  e$shape <- shape
  e$rate <- rate
  with(
    e,
    mcstate_model(
      list(
        parameters = "gamma",
        direct_sample = function(rng) {
          rng$gamma(1, shape = shape, scale = 1 / rate)
        },
        density = function(x) {
          drop(dgamma(x, shape = shape, rate = rate, log = TRUE))
        },
        gradient = function(x) drop((shape - 1) / x - rate),
        domain = rbind(c(0, Inf))),
      mcstate_model_properties(allow_multiple_parameters = TRUE)))
}


ex_simple_nested <- function(n_groups) {
  e <- new.env(parent = topenv())
  e$n_groups <- n_groups
  e$i <- as.numeric(seq_len(n_groups))
  with(
    e,
    mcstate_model(list(
      parameters = paste0("mu_", i),
      direct_sample = function(rng) {
        rng$normal(n_groups, i, 1)
      },
      density = function(x, by_group = FALSE) {
        z <- dnorm(x, 2^i, 1, log = TRUE)
        value <- sum(z)
        if (by_group) structure(value, "by_group" = z) else value
      },
      parameter_groups = i)))
}


## This is example is really stupid and does not really help test much
## except that the code can run at all. Ed/Marc; can you think of a
## reasonable example here please?
ex_simple_nested_with_base <- function(n_groups) {
  e <- new.env(parent = topenv())
  e$sigma <- 5
  e$mu <- rnorm(n_groups, 0, e$sigma)
  e$n_groups <- n_groups
  with(
    e,
    mcstate_model(list(
      parameters = c("sigma", paste0("mu_", seq_len(n_groups))),
      direct_sample = function(rng) {
        sigma <- rng$uniform(1, 0, 10)
        c(sigma, rng$normal(n_groups, 0, sigma))
      },
      density = function(x, by_group = FALSE) {
        sigma <- x[[1]]
        if (sigma <= 0) {
          z <- rep(-Inf, length(x) - 1)
        } else {
          z <- dnorm(x[-1], 0, x[[1]], log = TRUE)
        }
        value <- sum(z) + dunif(x[[1]], 0, 10, log = TRUE)
        if (by_group) structure(value, "by_group" = z) else value
      },
      parameter_groups = c(0, seq_len(n_groups)),
      mu = mu)))
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
  details <- function(idx_particle) {
    if (save_trajectories) {
      traj <- trajectories[, idx_particle, , drop = FALSE]
      dim(traj) <- dim(traj)[-2]
    } else {
      traj <- NULL
    }
    list(trajectories = traj, state = model$state()[, idx_particle])
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
  mcstate_model(
    list(parameters = c("a", "b"),
         direct_sample = function(rng) {
           b <- rng$random_normal(1)
           a <- rng$normal(1, b^2, sd)
           c(a, b)
         },
         density = function(x) {
           if (length(dim2(x)) == 1) {
             a <- x[1]
             b <- x[2]
           } else {
             a <- x[1, ]
             b <- x[2, ]
           }
           dnorm(b, log = TRUE) + dnorm((a - b^2) / sd, log = TRUE)
         },
         gradient = function(x) {
           if (length(dim2(x)) == 1) {
             a <- x[1]
             b <- x[2]
             c((b^2 - a) / sd^2,
               -b + 2 * b * (a - b^2) / sd^2)
           } else {
             a <- x[1, ]
             b <- x[2, ]
             rbind((b^2 - a) / sd^2,
                   -b + 2 * b * (a - b^2) / sd^2)
           }
         },
         domain = cbind(rep(-Inf, 2), rep(Inf, 2))),
    mcstate_model_properties(allow_multiple_parameters = TRUE))
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
