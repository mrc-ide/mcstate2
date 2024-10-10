test_that("can draw samples from a trivial model", {
  m <- ex_simple_gamma1()
  sampler <- monty_sampler_random_walk(vcv = matrix(0.01, 1, 1))
  res <- monty_sample(m, sampler, 100)
  expect_equal(names(res),
               c("pars", "density", "initial", "details", "observations"))
  expect_equal(dim(res$pars), c(1, 100, 1))
})


test_that("validate sampler against model on initialisation", {
  m <- ex_simple_gamma1()

  state <- list(pars = 1, density = -Inf)
  sampler1 <- monty_sampler_random_walk(vcv = diag(1) * 0.01)
  sampler2 <- monty_sampler_random_walk(vcv = diag(2) * 0.01)
  r <- monty_rng$new()

  expect_no_error(sampler1$initialise(1, m, r))
  expect_error(
    sampler2$initialise(1, m, r),
    "Incompatible length parameters (1) and vcv (2)",
    fixed = TRUE)
})


test_that("can draw samples from a random model", {
  set.seed(1)
  m <- ex_dust_sir()
  vcv <- matrix(c(0.0006405, 0.0005628, 0.0005628, 0.0006641), 2, 2)
  sampler <- monty_sampler_random_walk(vcv = vcv)
  res <- monty_sample(m, sampler, 20)
  expect_setequal(names(res),
                  c("pars", "density", "initial", "details", "observations"))
})


test_that("can observe a model", {
  skip("FIXME: add model-based observer")
  m <- ex_dust_sir(save_trajectories = TRUE)
  vcv <- matrix(c(0.0006405, 0.0005628, 0.0005628, 0.0006641), 2, 2)
  sampler <- monty_sampler_random_walk(vcv = vcv)

  observer <- monty_observer(
    function(model, rng) {
      i <- floor(rng$random_real(1) * model$model$n_particles()) + 1L
      model$details(i)
    })

  ## This takes quite a while, and that seems mostly to be the time
  ## taken to call the filter in dust.
  res <- monty_sample(m, sampler, 20, n_chains = 3, observer = observer)
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
  skip("FIXME: add model-based observer")
  m <- ex_dust_sir(save_trajectories = TRUE)
  vcv <- matrix(c(0.0006405, 0.0005628, 0.0005628, 0.0006641), 2, 2)
  sampler <- monty_sampler_random_walk(vcv = vcv)

  observer <- monty_observer(
    function(model, rng) {
      ## ideally we get a random sample here, but that's not easy with
      ## current dust
      model$details(4)
    })

  set.seed(1)
  res1 <- monty_sample(m, sampler, 15, n_chains = 3, observer = observer)

  set.seed(1)
  res2a <- monty_sample(m, sampler, 5, n_chains = 3, observer = observer,
                          restartable = TRUE)
  res2b <- monty_sample_continue(res2a, 10)

  expect_equal(res1$observations, res2b$observations)
})


test_that("can run multiple samples at once", {
  m <- ex_simple_gamma1()
  sampler <- monty_sampler_random_walk(vcv = matrix(0.01, 1, 1))
  p <- matrix(runif(5), 1)
  ## TODO: we need a much better rng support here; we'll need to make
  ## a tweak to the rng code to to a long jump between each chain.
  r <- monty_rng$new(n_streams = 5)
  state0 <- sampler$initialise(p, m, r)
  state1 <- sampler$step(state0, m, r)

  expect_equal(dim2(state0$pars), c(1, 5))
  expect_equal(dim2(state1$pars), c(1, 5))
  expect_equal(dim2(state0$density), 5)
  expect_equal(dim2(state1$density), 5)
})


test_that("can run random walk sampler simultaneously", {
  m <- ex_simple_gamma1()
  sampler <- monty_sampler_random_walk(vcv = matrix(0.01, 1, 1))

  set.seed(1)
  res1 <- monty_sample(m, sampler, 100, n_chains = 3)

  set.seed(1)
  runner <- monty_runner_simultaneous()
  res2 <- monty_sample(m, sampler, 100, n_chains = 3, runner = runner)
  expect_equal(res1, res2)
})


test_that("can continue a simultaneous random walk sampler", {
  m <- ex_simple_gamma1()
  sampler <- monty_sampler_random_walk(vcv = matrix(0.01, 1, 1))

  set.seed(1)
  res1a <- monty_sample(m, sampler, 30, n_chains = 3, restartable = TRUE)
  res1b <- monty_sample_continue(res1a, 70)

  set.seed(1)
  runner <- monty_runner_simultaneous()
  res2a <- monty_sample(m, sampler, 30, n_chains = 3, runner = runner,
                          restartable = TRUE)
  expect_equal(res2a$restart$state, res1a$restart$state)

  res2b <- monty_sample_continue(res2a, 70)

  expect_equal(res2b, res1b)
})


test_that("validate number of parameter sets", {
  m <- ex_simple_gamma1()
  vcv <- array(1 * 10^c(-4, -3, -2, -1, 0), c(1, 1, 5))
  sampler <- monty_sampler_random_walk(vcv = vcv)
  expect_error(
    monty_sample(m, sampler, 200, n_chains = 5),
    "Incompatible number of parameter sets (1) and slices in vcv (5)",
    fixed = TRUE)
  runner <- monty_runner_simultaneous()
  expect_error(
    monty_sample(m, sampler, 200, n_chains = 3, runner = runner),
    "Incompatible number of parameter sets (3) and slices in vcv (5)",
    fixed = TRUE)
})


test_that("can reflect parameters", {
  expect_equal(reflect_proposal(0.4, 0, 1), 0.4)
  expect_equal(reflect_proposal(1.4, 0, 1), 0.6)
  expect_equal(reflect_proposal(2.4, 0, 1), 0.4)

  expect_equal(reflect_proposal(6, 1, 5), 4)
  expect_equal(reflect_proposal(0, 1, 5), 2)
  expect_equal(reflect_proposal(10, 1, 5), 2)

  expect_equal(reflect_proposal(0:2 + 0.4, rep(0, 3), rep(1, 3)),
               c(0.4, 0.6, 0.4))
  expect_equal(reflect_proposal(0:2 + 1.4, rep(1, 3), rep(2, 3)),
               c(1.4, 1.6, 1.4))
})


test_that("Reflect parameters: infinite", {
  p <- rnorm(10)
  expect_equal(reflect_proposal(p, rep(-Inf, 5), rep(Inf, 5)), p)
})


test_that("reflect parameters: lower", {
  expect_equal(reflect_proposal(-0.4, 0, Inf), 0.4)
  expect_equal(reflect_proposal(0.6, 1, Inf), 1.4)
  expect_equal(reflect_proposal(0.6, 0, Inf), 0.6)
})


test_that("reflect parameters: upper", {
  expect_equal(reflect_proposal(1.4, -Inf, 1), 0.6)
  expect_equal(reflect_proposal(0.4, -Inf, 0), -0.4)
  expect_equal(reflect_proposal(0.4, -Inf, 1), 0.4)
})


test_that("can reflect parameters in matrix form", {
  p <- matrix(1:12, 3, 4) + .1

  expect_equal(reflect_proposal(p, rep(-Inf, 3), rep(Inf, 3)), p)

  x_min <- c(3, -Inf, 9)
  x_max <- c(5, Inf, 10)
  expect_equal(reflect_proposal(p, x_min, x_max),
               sapply(1:4, function(i) reflect_proposal(p[, i], x_min, x_max)))

  xa <- c(3, 4, 9)
  xb <- c(Inf, Inf, Inf)
  expect_equal(reflect_proposal(p, xa, xb),
               sapply(1:4, function(i) reflect_proposal(p[, i], xa, xb)))
  expect_equal(reflect_proposal(p, -xb, xa),
               sapply(1:4, function(i) reflect_proposal(p[, i], -xb, xa)))


  mvn <- make_rmvnorm(diag(3) * c(4, 8, 16), centred = TRUE)
  r <- monty_rng$new(seed = 42)
  x <- replicate(10, mvn(r))
  x_min <- c(-2, -Inf, -1)
  x_max <- c(2, Inf, 1)
  expect_equal(reflect_proposal(x, x_min, x_max),
               apply(x, 2, reflect_proposal, x_min, x_max))
})


test_that("Can create a reflected random walk proposal", {
  vcv <- diag(3) * c(4, 8, 16)
  domain <- cbind(c(-2, -Inf, -1), c(2, Inf, 1))
  x <- rep(0, 3)

  p1 <- make_random_walk_proposal(vcv, domain, "ignore")
  p2 <- make_random_walk_proposal(vcv, domain[c(2, 2, 2), ], "reflect")
  p3 <- make_random_walk_proposal(vcv, domain, "reflect")

  expect_equal(body(p2), body(p1))

  r1 <- monty_rng$new(seed = 42)
  r2 <- monty_rng$new(seed = 42)
  r3 <- monty_rng$new(seed = 42)

  x1 <- replicate(10, p1(x, r1))
  x2 <- replicate(10, p2(x, r2))
  x3 <- replicate(10, p3(x, r3))

  expect_identical(x1, x2)
  expect_equal(x3, reflect_proposal(x1, domain[, 1], domain[, 2]))
})


test_that("can run sampler with reflecting boundaries", {
  model <- monty_model(
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

  s1 <- monty_sampler_random_walk(matrix(0.5, 1, 1), boundaries = "ignore")
  s2 <- monty_sampler_random_walk(matrix(0.5, 1, 1), boundaries = "reflect")
  s3 <- monty_sampler_random_walk(matrix(0.5, 1, 1), boundaries = "reject")

  expect_error(monty_sample(model, s1, 100), "parameter out of bounds")

  res2 <- monty_sample(model, s2, 100)
  r2 <- range(drop(res2$pars))
  expect_gt(diff(r2), 0.75)
  expect_gt(r2[[1]], -1)
  expect_lt(r2[[2]], 1)

  res3 <- monty_sample(model, s3, 100)
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
  model <- monty_model(
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

  s1 <- monty_sampler_random_walk(matrix(0.5, 1, 1), boundaries = "ignore")
  s2 <- monty_sampler_random_walk(matrix(0.5, 1, 1), boundaries = "reject")

  expect_error(monty_sample(model, s1, 100), "parameter out of bounds")
  res <- monty_sample(model, s2, 100)
  r <- range(drop(res$pars))
  expect_gt(diff(r), 0.75)
  expect_gt(r[[1]], -1)
  expect_lt(r[[2]], 1)
})


test_that("can run sampler with rejecting boundaries simultaneously", {
  m <- monty_model(
    list(parameters = "x",
         domain = cbind(-1, 1),
         density = function(x) {
           ifelse(abs(x) > 1, NA_real_, log(0.5))
         },
         direct_sample = function(rng) {
           rng$uniform(1, -1, 1)
         }),
    monty_model_properties(allow_multiple_parameters = TRUE))

  s <- monty_sampler_random_walk(matrix(0.5, 1, 1), boundaries = "reject")
  runner <- monty_runner_simultaneous()

  n_steps <- 30

  set.seed(1)
  res <- monty_sample(m, s, n_steps, n_chains = 4, runner = runner)
  set.seed(1)
  cmp <- monty_sample(m, s, n_steps, n_chains = 4)

  expect_equal(res, cmp)
})


test_that("can print a sampler object", {
  s <- monty_sampler_random_walk(diag(1))
  res <- evaluate_promise(withVisible(print(s)))
  expect_mapequal(res$result, list(value = s, visible = FALSE))
  expect_match(
    res$messages,
    "<monty_sampler: Random walk (monty_random_walk)>",
    fixed = TRUE, all = FALSE)
})


test_that("Can rerun a stochastic model", {
  set.seed(1)
  m <- ex_dust_sir()
  vcv <- matrix(c(0.0006405, 0.0005628, 0.0005628, 0.0006641), 2, 2)
  sampler1 <- monty_sampler_random_walk(vcv = vcv, rerun_every = 2)
  sampler2 <- monty_sampler_random_walk(vcv = vcv, rerun_every = 2, 
                                        rerun_random = TRUE)
  res1 <- monty_sample(m, sampler1, 20)
  res2 <- monty_sample(m, sampler2, 20)
  
  expect_gt(sum(diff(res1$density) != 0), sum(diff(res1$pars[1, , 1]) != 0))
  expect_true(all(diff(res1$density)[seq(1, 20, by = 2)] != 0))
  expect_gt(sum(diff(res2$density) != 0), sum(diff(res2$pars[1, , 1]) != 0))
})
