test_that("can draw samples from a trivial model", {
  m <- ex_simple_gamma1()
  sampler <- mcstate_sampler_random_walk(vcv = matrix(0.01, 1, 1))
  res <- mcstate_sample(m, sampler, 100)
  expect_equal(names(res),
               c("pars", "density", "initial", "details", "observations"))
  expect_equal(dim(res$pars), c(1, 100, 1))
})


test_that("validate construction of rw sampler", {
  expect_error(
    mcstate_sampler_random_walk(),
    "One of 'proposal' or 'vcv' must be given")
  expect_error(
    mcstate_sampler_random_walk(identity, matrix(1, 1, 1)),
    "Only one of 'proposal' or 'vcv' may be given")
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
