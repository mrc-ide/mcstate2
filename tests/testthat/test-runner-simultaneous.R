test_that("require multiple parameter support for simultaneous runner", {
  m <- mcstate_model(list(density = function(x) dnorm(x, log = TRUE),
                          parameters = "a"))
  sampler <- mcstate_sampler_random_walk(vcv = diag(1) * 0.01)
  runner <- mcstate_runner_simultaneous()
  expect_error(
    mcstate_sample(m, sampler, 100, initial = 1, runner = runner),
    "mcstate_runner_simultaneous requires support for multiple parameters")
})


test_that("require multiple parameter support for simultaneous runner", {
  m <- mcstate_model(list(density = function(x) dnorm(x, log = TRUE),
                          set_rng_state = identity,
                          get_rng_state = identity,
                          parameters = "a"),
                     mcstate_model_properties(is_stochastic = TRUE,
                                              allow_multiple_parameters = TRUE))
  sampler <- mcstate_sampler_random_walk(vcv = diag(1) * 0.01)
  runner <- mcstate_runner_simultaneous()
  expect_error(
    mcstate_sample(m, sampler, 100, initial = 1, runner = runner),
    "Can't yet use multiple parameter sets with stochastic model",
    fixed = TRUE)
})


test_that("can't provide observer yet with simultaneous runner", {
  m <- ex_simple_gamma1()
  observer <- mcstate_observer(function(model, rng) 1)
  sampler <- mcstate_sampler_random_walk(vcv = diag(1) * 0.01)
  runner <- mcstate_runner_simultaneous()
  expect_error(
    mcstate_sample(m, sampler, 100, runner = runner, observer = observer),
    "Can't yet use observers with 'mcstate_runner_simultaneous'",
    fixed = TRUE)
})


test_that("can use different covariance kernels on different chains", {
  m <- ex_simple_gamma1()
  vcv <- array(1 * 10^c(-4, -3, -2, -1, 0), c(1, 1, 5))
  sampler <- mcstate_sampler_random_walk(vcv = vcv)
  runner <- mcstate_runner_simultaneous()
  n_steps <- 30

  set.seed(1)
  res <- mcstate_sample(m, sampler, n_steps, n_chains = 5, runner = runner)

  cmp <- vapply(1:5, function(i) {
    set.seed(1)
    s <- mcstate_sampler_random_walk(vcv = matrix(vcv[, , i], 1, 1))
    mcstate_sample(m, s, n_steps, n_chains = 5)$density[, i]
  }, numeric(n_steps))
  expect_equal(res$density, cmp)
})


test_that("can use the same covariance kernels on different chains", {
  m <- ex_simple_gamma1()
  runner <- mcstate_runner_simultaneous()
  n_steps <- 30

  set.seed(1)
  sampler1 <- mcstate_sampler_random_walk(vcv = array(0.01, c(1, 1, 1)))
  res1 <- mcstate_sample(m, sampler1, n_steps, n_chains = 5, runner = runner)

  set.seed(1)
  sampler2 <- mcstate_sampler_random_walk(vcv = array(0.01, c(1, 1)))
  res2 <- mcstate_sample(m, sampler2, n_steps, n_chains = 5, runner = runner)

  expect_equal(res1, res2)
})


test_that("can simultaneously sample a single chain", {
  m <- ex_simple_gamma1()
  vcv <- array(1 * 10^c(-4, -3, -2), c(1, 1, 3))

  runner <- mcstate_runner_simultaneous()
  n_steps <- 30

  set.seed(1)
  sampler1 <- mcstate_sampler_random_walk(vcv = vcv[, , 1, drop = FALSE])
  res1 <- mcstate_sample(m, sampler1, n_steps, n_chains = 1, runner = runner)

  set.seed(1)
  sampler2 <- mcstate_sampler_random_walk(vcv = matrix(vcv[[1]], 1, 1))
  res2 <- mcstate_sample(m, sampler2, n_steps, n_chains = 1, runner = runner)

  set.seed(1)
  sampler3 <- mcstate_sampler_random_walk(vcv = vcv)
  res3 <- mcstate_sample(m, sampler3, n_steps, n_chains = 3, runner = runner)

  expect_equal(res1, res2)
  expect_equal(res3$density[, 1, drop = FALSE], res1$density)
})
