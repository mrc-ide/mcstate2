test_that("can run a simple nested random walk", {
  set.seed(1)
  ng <- 3
  m <- ex_simple_nested(ng)
  v <- list(base = NULL,
            groups = rep(list(matrix(50)), ng))
  s <- mcstate_sampler_nested_random_walk(v)
  res <- mcstate_sample(m, s, 100)

  ## Something we can look for is that acceptances are not equally
  ## shared; I've used a very poor vcv here so that acceptance is
  ## relatively rare and we can see this:
  accept <- t(diff(t(res$pars[, , 1])) != 0)
  expect_setequal(colSums(accept), 0:3)
  expect_setequal(colSums(accept * 2^(seq_len(ng) - 1)), 0:7)
})


test_that("validate vcv inputs on construction of sampler", {
  expect_error(
    mcstate_sampler_nested_random_walk(NULL),
    "Expected a list for 'vcv'")
  expect_error(
    mcstate_sampler_nested_random_walk(list()),
    "Expected 'vcv' to have elements 'base' and 'groups'")
  expect_error(
    mcstate_sampler_nested_random_walk(list(base = TRUE, groups = TRUE)),
    "Expected a matrix for 'vcv$base'",
    fixed = TRUE)
  expect_error(
    mcstate_sampler_nested_random_walk(list(base = NULL, groups = TRUE)),
    "Expected 'vcv$groups' to be a list",
    fixed = TRUE)
  expect_error(
    mcstate_sampler_nested_random_walk(list(base = NULL, groups = list())),
    "Expected 'vcv$groups' to have at least one element",
    fixed = TRUE)
  expect_error(
    mcstate_sampler_nested_random_walk(list(base = NULL, groups = list(TRUE))),
    "Expected a matrix for 'vcv$groups[1]'",
    fixed = TRUE)

  vcv <- list(base = diag(1), groups = list(diag(2), diag(3)))
  expect_s3_class(mcstate_sampler_nested_random_walk(vcv), "mcstate_sampler")
})


test_that("can't use nested sampler with models that are not nested", {
  m <- ex_simple_gaussian(diag(3))
  vcv <- list(base = diag(1), groups = list(diag(2), diag(3)))
  s <- mcstate_sampler_nested_random_walk(vcv)
  expect_error(
    mcstate_sample(m, s, 100),
    "Your model does not have parameter groupings")
})


test_that("Validate that the model produces correct parameter groupings", {
  m <- mcstate_model(list(
    parameters = c("a", "b", "c"),
    density = function(x, by_group = TRUE) {
      x
    },
    parameter_groups = c(0, 1, 2)))
  s <- mcstate_sampler_nested_random_walk(
    list(base = diag(1), groups = list(diag(1), diag(1))))
  expect_error(
    mcstate_sample(m, s, 100, c(0, 0, 0)),
    "model$density(x, by_group = TRUE) did not produce a density",
    fixed = TRUE)
})


test_that("Validate length of by group output", {
  m <- mcstate_model(list(
    parameters = c("a", "b", "c"),
    density = function(x, by_group = TRUE) {
      structure(x, by_group = rep(0, 5))
    },
    parameter_groups = c(0, 1, 2)))
  s <- mcstate_sampler_nested_random_walk(
    list(base = diag(1), groups = list(diag(1), diag(1))))
  expect_error(
    mcstate_sample(m, s, 100, c(0, 0, 0)),
    "model$density(x, by_group = TRUE) produced a 'by_group' attribute with",
    fixed = TRUE)
})


test_that("can build nested proposal functions", {
  v <- matrix(c(1, .9999, .9999, 1), 2, 2)
  vcv <- list(base = NULL,
              groups = list(v, v / 100))
  g <- c(1, 1, 2, 2)
  f <- nested_proposal(vcv, g)

  expect_null(f$base)
  expect_true(is.function(f$groups))

  r <- mcstate_rng$new(seed = 1)
  y <- replicate(10000, f$groups(rep(0, 4), r))
  v <- var(t(y))
  cmp <- rbind(c(1, 1, 0, 0),
               c(1, 1, 0, 0),
               c(0, 0, 0.01, 0.01),
               c(0, 0, 0.01, 0.01))
  expect_equal(v, cmp, tolerance = 1e-2)
})


test_that("can build nested proposal functions with base components", {
  v <- matrix(c(1, .5, .5, 1), 2, 2)
  vcv <- list(base = v / 10,
              groups = list(v, v / 100))
  f <- nested_proposal(vcv, c(0, 0, 1, 1, 2, 2))
  g <- nested_proposal(list(base = NULL, groups = vcv$groups), c(1, 1, 2, 2))

  expect_true(is.function(f$base))
  expect_true(is.function(f$groups))

  r <- mcstate_rng$new(seed = 1)
  y <- replicate(20000, f$base(1:6, r))
  expect_equal(y[3:6, ], matrix(3:6, 4, 20000))
  expect_equal(var(t(y[1:2, ])), vcv$base, tolerance = 0.01)

  r1 <- mcstate_rng$new(seed = 1)
  r2 <- mcstate_rng$new(seed = 1)
  expect_equal(
    replicate(100, f$groups(1:6, r1)[3:6]),
    replicate(100, g$groups(3:6, r2)))

  expect_equal(replicate(100, f$groups(1:6, r1)[1:2]),
               matrix(1:2, 2, 100))
})


test_that("validate that the proposal and model are compatible", {
  vcv <- list(base = diag(2), groups = list(diag(1), diag(2), diag(3)))
  expect_error(
    nested_proposal(vcv, c(1, 2, 2, 3, 3, 3)),
    "Incompatible number of base parameters in your model and sampler")
  expect_error(
    nested_proposal(vcv, c(0, 0, 1:5)),
    "Incompatible number of parameter groups in your model and sampler")
  expect_error(
    nested_proposal(vcv, c(0, 0, 1, 1, 2, 2, 3, 3, 3)),
    "Incompatible number of parameters within parameter group")
})


test_that("can continue nested sampler correctly", {
  set.seed(1)
  ng <- 3
  m <- ex_simple_nested(ng)
  v <- list(base = NULL,
            groups = rep(list(matrix(1)), ng))
  sampler <- mcstate_sampler_nested_random_walk(v)

  set.seed(1)
  res1 <- mcstate_sample(m, sampler, 30, n_chains = 3, restartable = TRUE)

  set.seed(1)
  res2a <- mcstate_sample(m, sampler, 10, n_chains = 3, restartable = TRUE)
  res2b <- mcstate_sample_continue(res2a, 20, restartable = TRUE)

  expect_equal(res2b, res1)
})


test_that("can run a sampler with shared parameters", {
  set.seed(1)
  ng <- 5
  m <- ex_simple_nested_with_base(ng)
  s <- mcstate_sampler_nested_random_walk(
    list(base = diag(1), groups = rep(list(diag(1)), ng)))
  res <- mcstate_sample(m, s, 100)
  ## This is not great, until we get a real example we can use; but it
  ## does test that the sampler runs, and that it mixes.
  expect_true(length(unique(res$pars[1, , ])) > 1)
})


test_that("can run an observer during a nested fit", {
  set.seed(1)
  ng <- 5
  m <- ex_simple_nested_with_base(ng)
  s <- mcstate_sampler_nested_random_walk(
    list(base = diag(1), groups = rep(list(diag(1)), ng)))
  counter <- 0
  observer <- mcstate_observer(function(...) {
    counter <<- counter + 1
    list(n = counter)
  })
  res <- mcstate_sample(m, s, 100, observer = observer)
  expect_equal(
    dim(res$observations$n),
    c(1, 100, 1))
  expect_gt(max(res$observations$n), 120) # called way more than once per step
})
