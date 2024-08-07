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


test_that("can run adaptive sampler simultaneously", {
  m <- ex_simple_gamma1()
  sampler <- mcstate_sampler_adaptive(initial_vcv = matrix(0.01, 1, 1))

  set.seed(1)
  res1 <- mcstate_sample(m, sampler, 100, n_chains = 3)

  set.seed(1)
  runner <- mcstate_runner_simultaneous()
  res2 <- mcstate_sample(m, sampler, 100, n_chains = 3, runner = runner)
  expect_equal(res1, res2)
})

test_that("can run sampler with reflecting boundaries", {
  model <- mcstate_model(
    list(parameters = "x",
         domain = cbind(-1, 1),
         density = function(x) {
           if (abs(x) > 1) {
             stop("parameter out of bounds")
           }
           0.5
         },
         direct_sample = function(rng) {
           rng$uniform(1, -1, 1)
         }))
  
  s1 <- mcstate_sampler_adaptive(matrix(0.5, 1, 1), boundaries = "ignore")
  s2 <- mcstate_sampler_adaptive(matrix(0.5, 1, 1), boundaries = "reflect")
  s3 <- mcstate_sampler_adaptive(matrix(0.5, 1, 1), boundaries = "reject")
  
  expect_error(mcstate_sample(model, s1, 100), "parameter out of bounds")
  
  ## Seems to be a bit of a problem here - since the density is flat, the
  ## acceptance probability is always 1. The effect on the scaling then is that
  ## it just keeps increasing! When it gets large enough, it causes problems
  ## with the reflecting, which throws some warnings about probable complete
  ## loss of accuracy in modulus, and results in values lying on the boundary
  res2 <- mcstate_sample(model, s2, 100)
  r2 <- range(drop(res2$pars))
  expect_gt(diff(r2), 0.75)
  expect_gt(r2[[1]], -1)
  expect_lt(r2[[2]], 1)
  
  res3 <- mcstate_sample(model, s3, 100)
  r3 <- range(drop(res3$pars))
  expect_gt(diff(r3), 0.75)
  expect_gt(r3[[1]], -1)
  expect_lt(r3[[2]], 1)
  
  ## Different with rejection than reflection, and more step
  ## rejections when rejection used.
  expect_true(!all(res2$pars == res3$pars))
  expect_gt(sum(diff(drop(res3$pars)) == 0),
            sum(diff(drop(res2$pars)) == 0))
})


test_that("can run sampler with rejecting boundaries", {
  model <- mcstate_model(
    list(parameters = "x",
         domain = cbind(-1, 1),
         density = function(x) {
           if (abs(x) > 1) {
             stop("parameter out of bounds")
           }
           0.5
         },
         direct_sample = function(rng) {
           rng$uniform(1, -1, 1)
         }))
  
  s1 <- mcstate_sampler_adaptive(matrix(0.5, 1, 1), boundaries = "ignore")
  s2 <- mcstate_sampler_adaptive(matrix(0.5, 1, 1), boundaries = "reject")
  
  expect_error(mcstate_sample(model, s1, 100), "parameter out of bounds")
  res <- mcstate_sample(model, s2, 100)
  r <- range(drop(res$pars))
  expect_gt(diff(r), 0.75)
  expect_gt(r[[1]], -1)
  expect_lt(r[[2]], 1)
})


test_that("can run sampler with rejecting boundaries simultaneously", {
  m <- mcstate_model(
    list(parameters = "x",
         domain = cbind(-1, 1),
         density = function(x) {
           ifelse(abs(x) > 1, NA_real_, log(0.5))
         },
         direct_sample = function(rng) {
           rng$uniform(1, -1, 1)
         }),
    mcstate_model_properties(allow_multiple_parameters = TRUE))
  
  s <- mcstate_sampler_adaptive(matrix(0.5, 1, 1), boundaries = "reject")
  runner <- mcstate_runner_simultaneous()
  
  n_steps <- 30
  
  set.seed(1)
  res <- mcstate_sample(m, s, n_steps, n_chains = 4, runner = runner)
  set.seed(1)
  cmp <- mcstate_sample(m, s, n_steps, n_chains = 4)
  
  expect_equal(res, cmp)
})
