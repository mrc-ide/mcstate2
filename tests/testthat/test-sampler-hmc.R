test_that("can run hmc", {
  m <- ex_simple_gaussian(diag(2))
  sampler <- mcstate_sampler_hmc(epsilon = 0.1, n_integration_steps = 10)
  set.seed(1)
  res <- mcstate_sample(m, sampler, 2000)
  ## We'll need some long tests here at some point!
  pars <- t(array_drop(res$pars, 3))
  expect_equal(cov(unname(pars)), diag(2), tolerance = 0.2)
})


test_that("can run hmc on banana shaped model", {
  m <- ex_banana()
  sampler <- mcstate_sampler_hmc(epsilon = 0.1, n_integration_steps = 10)
  set.seed(1)
  res <- mcstate_sample(m, sampler, 2000)
  ## We'll need some long tests here at some point!
  pars <- t(array_drop(res$pars, 3))
  expect_equal(abs(mean(pars[, 2])) < 0.2, TRUE)
  expect_equal(abs(mean(pars[, 1]) - 1) < 0.4, TRUE)
  expect_equal(abs(cor(pars[, 1], pars[, 2])) < 0.3, TRUE)
})


test_that("can output debug traces", {
  m <- ex_simple_gaussian(diag(2))

  set.seed(1)
  sampler1 <- mcstate_sampler_hmc(epsilon = 0.1, n_integration_steps = 10)
  res1 <- mcstate_sample(m, sampler1, 30)

  set.seed(1)
  sampler2 <- mcstate_sampler_hmc(epsilon = 0.1, n_integration_steps = 10,
                                  debug = TRUE)
  res2 <- mcstate_sample(m, sampler2, 30)

  expect_null(res1$details)
  expect_length(res2$details, 1)
  expect_equal(names(res2$details[[1]]), "debug")
  expect_equal(res2$details[[1]]$accept, rep(TRUE, 30))

  pars <- array_drop(res2$pars, 3)
  debug_pars <- res2$details[[1]]$pars
  expect_equal(dim(debug_pars), c(2, 11, 30))
  expect_equal(dimnames(debug_pars), list(c("a", "b"), NULL, NULL))
  expect_equal(debug_pars[, 1, ], cbind(res2$initial, pars[, -30]))
  expect_equal(debug_pars[, 11, ], pars)
})


test_that("can use custom vcv", {
  m <- ex_simple_gaussian(diag(2))
  set.seed(1)
  sampler1 <- mcstate_sampler_hmc(epsilon = 0.1, n_integration_steps = 10)
  res1 <- mcstate_sample(m, sampler1, 30)

  set.seed(1)
  sampler2 <- mcstate_sampler_hmc(epsilon = 0.1, n_integration_steps = 10,
                                  vcv = diag(2))
  res2 <- mcstate_sample(m, sampler2, 30)
  expect_identical(res2, res1)

  set.seed(1)
  sampler3 <- mcstate_sampler_hmc(epsilon = 0.1, n_integration_steps = 10,
                                  vcv = diag(2) * 0.5)
  res3 <- mcstate_sample(m, sampler3, 30)
  expect_false(identical(res3, res1))
})


test_that("given vcv and parameters must be compatible", {
  m <- ex_simple_gaussian(diag(2))
  set.seed(1)
  sampler <- mcstate_sampler_hmc(epsilon = 0.1, n_integration_steps = 10,
                                 vcv = diag(3))
  expect_error(
    mcstate_sample(m, sampler, 30),
    "Incompatible length parameters (2) and vcv (3)",
    fixed = TRUE)
})


test_that("can transform model parameters to R^n and back again", {
  transform <- hmc_transform(cbind(rep(-Inf, 3), rep(Inf, 3)))
  x <- c(0.1, 0.2, 0.3)
  expect_identical(transform$model2rn(x), x)
  expect_identical(transform$rn2model(x), x)
  expect_identical(transform$deriv(x), rep(1, 3))

  transform <- hmc_transform(cbind(rep(0, 3), rep(Inf, 3)))
  x <- c(0.5, 1.5, 2.5)
  expect_identical(transform$model2rn(x), log(x))
  expect_identical(transform$rn2model(log(x)), x)
  expect_identical(transform$deriv(x), x)

  lower <- c(0, 1, 2)
  upper <- c(1, 3, 6)
  transform <- hmc_transform(cbind(lower, upper))
  x <- c(0.5, 1.5, 2.5)
  expect_identical(transform$model2rn(x), logit_bounded(x, lower, upper))
  theta <- transform$model2rn(x)
  expect_equal(transform$rn2model(theta), x)
  expect_equal(transform$deriv(x), dilogit_bounded(x, lower, upper))

  ## A mixed case:
  transform <- hmc_transform(cbind(c(-Inf, 0, 0), c(Inf, Inf, 1)))
  x <- c(0.5, 0.6, 0.7)
  theta <- transform$model2rn(x)
  expect_equal(theta, c(0.5, log(0.6), logit_bounded(0.7, 0, 1)))
  expect_equal(transform$rn2model(theta), x)
  expect_equal(transform$deriv(x),
               c(1, 0.6, dilogit_bounded(0.7, 0, 1)))
})


test_that("prevent weird distributions", {
  expect_error(
    hmc_transform(cbind(rep(0, 3), c(1, Inf, -Inf))),
    "Unhandled domain type for parameter 3")
})


test_that("can continue hmc without debug", {
  m <- ex_simple_gaussian(diag(2))

  set.seed(1)
  sampler <- mcstate_sampler_hmc(epsilon = 0.1, n_integration_steps = 10)
  res1 <- mcstate_sample(m, sampler, 30, n_chains = 3, restartable = TRUE)

  set.seed(1)
  res2a <- mcstate_sample(m, sampler, 10, n_chains = 3, restartable = TRUE)
  res2b <- mcstate_sample_continue(res2a, 20, restartable = TRUE)

  expect_equal(res2b, res1)
})


test_that("can continue hmc with debug", {
  m <- ex_simple_gaussian(diag(2))

  set.seed(1)
  sampler <- mcstate_sampler_hmc(epsilon = 0.1, n_integration_steps = 10,
                                 debug = TRUE)
  res1 <- mcstate_sample(m, sampler, 30, n_chains = 3, restartable = TRUE)

  set.seed(1)
  res2a <- mcstate_sample(m, sampler, 10, n_chains = 3, restartable = TRUE)
  res2b <- mcstate_sample_continue(res2a, 20, restartable = TRUE)

  expect_equal(res2b, res1)
})


test_that("can't use hmc with models that lack gradients", {
  m <- mcstate_model(list(density = function(x) dnorm(x, log = TRUE),
                          parameters = "a"))
  sampler <- mcstate_sampler_hmc(epsilon = 0.1, n_integration_steps = 10)
  expect_error(
    mcstate_sample(m, sampler, 30, 1, n_chains = 3),
    "Can't use HMC without a gradient")
})


test_that("can't use hmc with stochastic models", {
  set.seed(1)
  m <- ex_dust_sir()
  sampler <- mcstate_sampler_hmc(epsilon = 0.1, n_integration_steps = 10)
  expect_error(
    mcstate_sample(m, sampler, 30, n_chains = 3),
    "Can't use HMC with stochastic models")
})


test_that("can run hmc model simultaneously", {
  m <- ex_banana()
  sampler <- mcstate_sampler_hmc(epsilon = 0.1, n_integration_steps = 10)
  runner <- mcstate_runner_simultaneous()
  set.seed(1)
  res1 <- mcstate_sample(m, sampler, 200, n_chains = 3)
  set.seed(1)
  res2 <- mcstate_sample(m, sampler, 200, n_chains = 3, runner = runner)
  expect_equal(res1, res2)
})


test_that("can run hmc model simultaneously, with debug", {
  m <- ex_banana()
  sampler <- mcstate_sampler_hmc(epsilon = 0.1, n_integration_steps = 10,
                                 debug = TRUE)
  runner <- mcstate_runner_simultaneous()
  set.seed(1)
  res1 <- mcstate_sample(m, sampler, 30, n_chains = 3)
  ## We don't yet do a good job of auto squashing the details.  I'll
  ## make this change in a future PR so it's more obvious (mrc-5293)
  res1$details <- observer_finalise_auto(res1$details)
  set.seed(1)
  res2 <- mcstate_sample(m, sampler, 30, n_chains = 3, runner = runner)
  expect_equal(res1, res2)
})
