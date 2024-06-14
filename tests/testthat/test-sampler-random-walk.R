test_that("can draw samples from a trivial model", {
  m <- ex_simple_gamma1()
  sampler <- mcstate_sampler_random_walk(vcv = matrix(0.01, 1, 1))
  res <- mcstate_sample(m, sampler, 100)
  expect_equal(names(res),
               c("pars", "density", "initial", "details", "observations"))
  expect_equal(dim(res$pars), c(1, 100, 1))
})


test_that("validate sampler against model on initialisation", {
  m <- ex_simple_gamma1()

  state <- list(pars = 1, density = -Inf)
  sampler1 <- mcstate_sampler_random_walk(vcv = diag(1) * 0.01)
  sampler2 <- mcstate_sampler_random_walk(vcv = diag(2) * 0.01)
  r <- mcstate_rng$new()

  expect_no_error(sampler1$initialise(1, m, NULL, r))
  expect_error(
    sampler2$initialise(1, m, NULL, r),
    "Incompatible length parameters (1) and vcv (2)",
    fixed = TRUE)
})


test_that("can draw samples from a random model", {
  set.seed(1)
  m <- ex_dust_sir()
  vcv <- matrix(c(0.0006405, 0.0005628, 0.0005628, 0.0006641), 2, 2)
  sampler <- mcstate_sampler_random_walk(vcv = vcv)
  res <- mcstate_sample(m, sampler, 20)
  expect_setequal(names(res),
                  c("pars", "density", "initial", "details", "observations"))
})


test_that("can observe a model", {
  m <- ex_dust_sir(save_trajectories = TRUE)
  vcv <- matrix(c(0.0006405, 0.0005628, 0.0005628, 0.0006641), 2, 2)
  sampler <- mcstate_sampler_random_walk(vcv = vcv)

  observer <- mcstate_observer(
    function(model, rng) {
      i <- floor(rng$random_real(1) * model$model$n_particles()) + 1L
      model$details(i)
    })

  ## This takes quite a while, and that seems mostly to be the time
  ## taken to call the filter in dust.
  res <- mcstate_sample(m, sampler, 20, n_chains = 3, observer = observer)
  expect_setequal(names(res),
                  c("pars", "density", "initial", "details", "observations"))
  expect_equal(names(res$observations),
               c("trajectories", "state"))
  expect_equal(dim(res$observations$trajectories),
               c(2, 151, 20, 3)) # states, time steps, samples, chains
  expect_equal(dim(res$observations$state),
               c(5, 20, 3)) # states, samples, chains
  expect_equal(res$observations$state[c(2, 4), , ],
               res$observations$trajectories[, 151, , ])
})


test_that("can continue observed models", {
  m <- ex_dust_sir(save_trajectories = TRUE)
  vcv <- matrix(c(0.0006405, 0.0005628, 0.0005628, 0.0006641), 2, 2)
  sampler <- mcstate_sampler_random_walk(vcv = vcv)

  observer <- mcstate_observer(
    function(model, rng) {
      # i <- floor(rng$random_real(1) * model$model$n_particles()) + 1L
      model$details(4)
    })

  set.seed(1)
  res1 <- mcstate_sample(m, sampler, 15, n_chains = 3, observer = observer)

  set.seed(1)
  res2a <- mcstate_sample(m, sampler, 5, n_chains = 3, observer = observer,
                          restartable = TRUE)
  res2b <- mcstate_sample_continue(res2a, 10)

  expect_equal(res1$observations, res2b$observations)
})


test_that("can run multiple samples at once", {
  m <- ex_simple_gamma1()
  sampler <- mcstate_sampler_random_walk(vcv = matrix(0.01, 1, 1))
  p <- matrix(runif(5), 1)
  ## TODO: we need a much better rng support here; we'll need to make
  ## a tweak to the rng code to to a long jump between each chain.
  r <- mcstate_rng$new(n_streams = 5)
  state0 <- sampler$initialise(p, m, NULL, r)
  state1 <- sampler$step(state0, m, NULL, r)

  expect_equal(dim2(state0$pars), c(1, 5))
  expect_equal(dim2(state1$pars), c(1, 5))
  expect_equal(dim2(state0$density), 5)
  expect_equal(dim2(state1$density), 5)
})


test_that("can run random walk sampler simultaneously", {
  m <- ex_simple_gamma1()
  sampler <- mcstate_sampler_random_walk(vcv = matrix(0.01, 1, 1))

  set.seed(1)
  res1 <- mcstate_sample(m, sampler, 100, n_chains = 3)

  set.seed(1)
  runner <- mcstate_runner_simultaneous()
  res2 <- mcstate_sample(m, sampler, 100, n_chains = 3, runner = runner)
  expect_equal(res1, res2)
})


test_that("can continue a simultaneous random walk sampler", {
  m <- ex_simple_gamma1()
  sampler <- mcstate_sampler_random_walk(vcv = matrix(0.01, 1, 1))

  set.seed(1)
  res1a <- mcstate_sample(m, sampler, 30, n_chains = 3, restartable = TRUE)
  res1b <- mcstate_sample_continue(res1a, 70)

  set.seed(1)
  runner <- mcstate_runner_simultaneous()
  res2a <- mcstate_sample(m, sampler, 30, n_chains = 3, runner = runner,
                          restartable = TRUE)
  expect_equal(res2a$restart$state, res1a$restart$state)

  res2b <- mcstate_sample_continue(res2a, 70)

  expect_equal(res2b, res1b)
})


test_that("validate number of parameter sets", {
  m <- ex_simple_gamma1()
  vcv <- array(1 * 10^c(-4, -3, -2, -1, 0), c(1, 1, 5))
  sampler <- mcstate_sampler_random_walk(vcv = vcv)
  expect_error(
    mcstate_sample(m, sampler, 200, n_chains = 5),
    "Incompatible number of parameter sets (1) and slices in vcv (5)",
    fixed = TRUE)
  runner <- mcstate_runner_simultaneous()
  expect_error(
    mcstate_sample(m, sampler, 200, n_chains = 3, runner = runner),
    "Incompatible number of parameter sets (3) and slices in vcv (5)",
    fixed = TRUE)
})
