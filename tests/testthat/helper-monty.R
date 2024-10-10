ex_simple_gamma1 <- function(shape = 1, rate = 1) {
  e <- new.env(parent = .GlobalEnv)
  e$shape <- shape
  e$rate <- rate
  with(
    e,
    monty_model(
      list(
        parameters = "gamma",
        direct_sample = function(rng) {
          rng$gamma_scale(1, shape = shape, scale = 1 / rate)
        },
        density = function(x) {
          drop(dgamma(x, shape = shape, rate = rate, log = TRUE))
        },
        gradient = function(x) drop((shape - 1) / x - rate),
        domain = rbind(c(0, Inf))),
      monty_model_properties(allow_multiple_parameters = TRUE)))
}


ex_simple_nested <- function(n_groups) {
  e <- new.env(parent = topenv())
  e$n_groups <- n_groups
  e$i <- as.numeric(seq_len(n_groups))
  with(
    e,
    monty_model(list(
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
    monty_model(list(
      parameters = c("sigma", paste0("mu_", seq_len(n_groups))),
      direct_sample = function(rng) {
        sigma <- rng$uniform(1, 0, 10)
        c(sigma, rng$normal(n_groups, 0, sigma))
      },
      density = function(x, by_group = FALSE) {
        density1 <- function(y) {
          sigma <- y[[1]]
          if (sigma <= 0) {
            rep(-Inf, length(y) - 1)
          } else {
            dnorm(y[-1], 0, y[[1]], log = TRUE)
          }
        }
        if (is.matrix(x)) {
          z <- vapply(seq_len(ncol(x)), function(i) density1(x[, i]),
                      numeric(nrow(x) - 1))
          value <- colSums(z) + dunif(x[1, ], 0, 10, log = TRUE)
          if (by_group) structure(value, "by_group" = z) else value
        } else {
          z <- density1(x)
          value <- sum(z) + dunif(x[[1]], 0, 10, log = TRUE)
          if (by_group) structure(value, "by_group" = z) else value
        }
      },
      parameter_groups = c(0, seq_len(n_groups)),
      mu = mu),
      monty_model_properties(allow_multiple_parameters = TRUE)))
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
    c(rng$gamma_scale(1, prior_beta_shape, 1 / prior_beta_rate),
      rng$gamma_scale(1, prior_gamma_shape, 1 / prior_gamma_rate))
  }

  set_rng_state <- function(rng_state) {
    n_streams <- n_particles + 1
    if (length(rng_state) != 32 * n_streams) {
      ## Expand the state by short jumps; we'll make this nicer once
      ## we refactor the RNG interface and dust.
      rng_state <- monty_rng$new(rng_state, n_streams)$state()
    }
    model$set_rng_state(rng_state)
  }

  get_rng_state <- function() {
    model$rng_state()
  }

  if (save_trajectories) {
    observer <- monty_observer(
      function() {
        ## TODO: It's not really clear to me (Rich) that we want the
        ## rng coming in here.  In dust2 we'll correctly use the
        ## *filter* rng.  however, even there we end up with the act
        ## of observation changing the behaviour of the sampler and I
        ## am not sure that is really desirable.  We could probably
        ## improve that in dust, but that would require that we do not
        ## pass the rng here too.  So for now we take the first
        ## particle.
        ## i <- floor(rng$random_real(1) * model$model$n_particles()) + 1L
        i <- 1L
        if (save_trajectories) {
          traj <- trajectories[, i, , drop = FALSE]
          dim(traj) <- dim(traj)[-2]
        } else {
          traj <- NULL
        }
        list(trajectories = traj, state = model$state()[, i])
      })
  } else {
    observer <- NULL
  }

  monty_model(
    list(model = model,
         density = density,
         direct_sample = direct_sample,
         parameters = c("beta", "gamma"),
         domain = cbind(c(0, 0), c(Inf, Inf)),
         observer = observer,
         set_rng_state = set_rng_state,
         get_rng_state = get_rng_state),
    monty_model_properties(is_stochastic = !deterministic))
}


random_array <- function(dim, named = FALSE) {
  if (named) {
    dn <- lapply(seq_along(dim), function(i) {
      paste0(LETTERS[[i]], letters[seq_len(dim[i])])
    })
    names(dn) <- paste0("d", LETTERS[seq_along(dim)])
  } else {
    dn <- NULL
  }
  array(runif(prod(dim)), dim, dimnames = dn)
}


ex_dust_sir_likelihood <- function(n_particles = 100, n_threads = 1,
                                   deterministic = FALSE,
                                   save_trajectories = FALSE) {
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

  trajectories <- NULL

  density <- function(x) {
    beta <- x[[1]]
    gamma <- x[[2]]
    model$update_state(
      pars = list(beta = x[[1]], gamma = x[[2]]),
      time = 0,
      set_initial_state = TRUE)
    res <- model$filter(save_trajectories = save_trajectories)
    if (save_trajectories) {
      trajectories <<- res$trajectories
    }
    res$log_likelihood
  }

  set_rng_state <- function(rng_state) {
    n_streams <- n_particles + 1
    if (length(rng_state) != 32 * n_streams) {
      ## Expand the state by short jumps; we'll make this nicer once
      ## we refactor the RNG interface and dust.
      rng_state <- monty_rng$new(rng_state, n_streams)$state()
    }
    model$set_rng_state(rng_state)
  }

  get_rng_state <- function() {
    model$rng_state()
  }

  if (save_trajectories) {
    observer <- monty_observer(
      function() {
        ## TODO: It's not really clear to me (Rich) that we want the
        ## rng coming in here.  In dust2 we'll correctly use the
        ## *filter* rng.  however, even there we end up with the act
        ## of observation changing the behaviour of the sampler and I
        ## am not sure that is really desirable.  We could probably
        ## improve that in dust, but that would require that we do not
        ## pass the rng here too.  So for now we take the first
        ## particle.
        ## i <- floor(rng$random_real(1) * model$model$n_particles()) + 1L
        i <- 1L
        if (save_trajectories) {
          traj <- trajectories[, i, , drop = FALSE]
          dim(traj) <- dim(traj)[-2]
        } else {
          traj <- NULL
        }
        list(trajectories = traj, state = model$state()[, i])
      })
  } else {
    observer <- NULL
  }

  monty_model(
    list(model = model,
         density = density,
         parameters = c("beta", "gamma"),
         observer = observer,
         set_rng_state = set_rng_state,
         get_rng_state = get_rng_state),
    monty_model_properties(is_stochastic = !deterministic))
}


scrub_manual_info <- function(x) {
  x <- sub("Manual monty sampling at '.+",
           "Manual monty sampling at '<PATH>'",
           x)
  x <- sub("Created .*[0-9]{4}.*",
           "Created <DATE>",
           x)
  x
}
