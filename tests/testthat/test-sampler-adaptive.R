test_that("Empirical VCV calculated correctly with forget_rate = 0", {
  m <- ex_simple_gaussian(vcv = rbind(c(0.02, 0.01), c(0.01, 0.03)))

  sampler <- mcstate_sampler_adaptive(initial_vcv = diag(c(0.01, 0.01)),
                                      forget_rate = 0,
                                      log_scaling_update = FALSE)
  res <- mcstate_sample(m, sampler, 1000)
  expect_equal(names(res),
               c("pars", "density", "initial", "details", "observations"))

  ## forget_rate = 0 so full chain should be included in VCV
  expect_equal(res$details[[1]]$weight, 1000)
  expect_equal(res$details[[1]]$included, seq_len(1000))
  pars <- t(array_drop(res$pars, 3))
  expect_equal(res$details[[1]]$vcv, cov(pars), ignore_attr = TRUE)
})


test_that("Empirical VCV calculated correctly with forget_rate = 0.1", {
  m <- ex_simple_gaussian(vcv = rbind(c(0.02, 0.01), c(0.01, 0.03)))

  sampler <- mcstate_sampler_adaptive(initial_vcv = diag(c(0.01, 0.01)),
                                      forget_rate = 0.1)
  res <- mcstate_sample(m, sampler, 1000)
  expect_equal(names(res),
               c("pars", "density", "initial", "details", "observations"))

  ## forget_rate = 0.1 so VCV should exclude first 100 parameter sets
  expect_equal(res$details[[1]]$weight, 900)
  expect_equal(res$details[[1]]$included, seq(101, 1000, by = 1))
  pars <- t(array_drop(res$pars, 3))
  expect_equal(res$details[[1]]$vcv, cov(pars[101:1000, ]),
               ignore_attr = TRUE)
})


test_that("Empirical VCV correct using both forget_rate and forget_end", {
  m <- ex_simple_gaussian(vcv = rbind(c(0.02, 0.01), c(0.01, 0.03)))

  sampler <- mcstate_sampler_adaptive(initial_vcv = diag(c(0.01, 0.01)),
                                      forget_rate = 0.5,
                                      forget_end = 200)
  res <- mcstate_sample(m, sampler, 1000)
  expect_equal(names(res),
               c("pars", "density", "initial", "details", "observations"))

  ## forget_rate = 0.5 and forget_end = 200 so VCV should exclude first
  ## 100 parameter sets
  expect_equal(res$details[[1]]$weight, 900)
  expect_equal(res$details[[1]]$included, seq(101, 1000, by = 1))
  pars <- t(array_drop(res$pars, 3))
  expect_equal(res$details[[1]]$vcv, cov(pars[101:1000,]),
               ignore_attr = TRUE)
})


test_that("Empirical VCV correct using forget_rate, forget_end and adapt_end", {
  m <- ex_simple_gaussian(vcv = rbind(c(0.02, 0.01), c(0.01, 0.03)))

  sampler <- mcstate_sampler_adaptive(initial_vcv = diag(c(0.01, 0.01)),
                                      forget_rate = 0.25,
                                      forget_end = 100,
                                      adapt_end = 300)
  res <- mcstate_sample(m, sampler, 1000)
  expect_equal(names(res),
               c("pars", "density", "initial", "details", "observations"))

  ## forget_rate = 0.25, forget_end = 500 and adapt_end = 300 so VCV should
  ## only include parameter sets 26 to 300
  expect_equal(res$details[[1]]$weight, 275)
  expect_equal(res$details[[1]]$included, seq(26, 300, by = 1))
  pars <- t(array_drop(res$pars, 3))
  expect_equal(res$details[[1]]$vcv, cov(pars[26:300,]),
               ignore_attr = TRUE)
})


test_that("can continue adaptive sampler", {
  m <- ex_simple_gaussian(vcv = rbind(c(0.02, 0.01), c(0.01, 0.03)))
  sampler <- mcstate_sampler_adaptive(initial_vcv = diag(c(0.01, 0.01)))

  set.seed(1)
  res1 <- mcstate_sample(m, sampler, 30, n_chains = 3, restartable = TRUE)

  set.seed(1)
  res2a <- mcstate_sample(m, sampler, 10, n_chains = 3, restartable = TRUE)
  res2b <- mcstate_sample_continue(res2a, 20, restartable = TRUE)

  expect_equal(res2b, res1)
})


test_that("can't use adaptive sampler with stochastic models", {
  set.seed(1)
  m <- ex_dust_sir()
  sampler <- mcstate_sampler_adaptive(initial_vcv = diag(c(0.01, 0.01)))
  expect_error(
    mcstate_sample(m, sampler, 30, n_chains = 3),
    "Can't use adaptive sampler with stochastic models")
})


## This can be done, but it's pretty tedious bookkeeping, and worth
## doing separately (once Ed has done the nested adaptive case too,
## which is even worse bookkeeping).
test_that("can't use adaptive sampler with simultaneous runner", {
  m <- ex_simple_gamma1()
  sampler <- mcstate_sampler_adaptive(initial_vcv = matrix(0.01, 1, 1))
  runner <- mcstate_runner_simultaneous()
  expect_error(
    mcstate_sample(m, sampler, 100, n_chains = 3, runner = runner),
    "Can't use 'mcstate_sampler_adaptive' with simultaneous chains")
})
