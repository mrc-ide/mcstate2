test_that("check argument matching", {
  candidate <- c("mean", "sd")
  expect_equal(match_call_candidate(list("a", "b"), candidate), 1:2)
  expect_null(match_call_candidate(list("a", "b", "c"), candidate))
  expect_equal(match_call_candidate(list("a", sd = "b"), candidate), 1:2)
  expect_equal(match_call_candidate(list(mean = "a", sd = "b"), candidate), 1:2)
  expect_null(match_call_candidate(list(sd = "a", sd = "b"), candidate), 1:2)
  expect_null(match_call_candidate(list(mu = "a", sd = "b"), candidate), 1:2)
  expect_equal(match_call_candidate(list(sd = "a", mean = "b"), candidate), 2:1)
})


test_that("match against set of candidates", {
  candidates <- list(c("a", "b"), c("shape", "scale"))
  expect_equal(match_call(list(0, 0), candidates),
               list(success = TRUE, index = 1, args = 1:2))
  expect_equal(match_call(list(0, b = 0), candidates),
               list(success = TRUE, index = 1, args = 1:2))
  expect_equal(match_call(list(a = 0, b = 0), candidates),
               list(success = TRUE, index = 1, args = 1:2))
  expect_equal(match_call(list(b = 0, a = 0), candidates),
               list(success = TRUE, index = 1, args = 2:1))
  expect_equal(match_call(list(0, scale = 0), candidates),
               list(success = TRUE, index = 2, args = 1:2))
  expect_equal(match_call(list(scale = 0, shape = 0), candidates),
               list(success = TRUE, index = 2, args = 2:1))

  error <- c(x = "Failed to match given arguments: .",
             i = "Call should match one of:",
             "*" = "a, b",
             "*" = "shape, scale")
  expect_equal(
    match_call(list(0), candidates),
    list(success = FALSE, error = error))

  error <- c(x = "Failed to match given arguments: 1, 2, 3",
             i = "Call should match:",
             "*" = "a, b")
  expect_equal(
    match_call(list(a = 1, b = 2, c = 3), candidates[1]),
    list(success = FALSE, error = error))
})


## These are all very poor tests in that they just restate the
## implementation.  I am not really clear how we can do better than
## this at the moment.  That said, there's a chance we'll swap out the
## implementation later for something symbolically differentiable so
## at that point these would become fairly useful.
test_that("can compute density for exponential distribution", {
  expect_equal(distr_exponential_rate$density(3.1, 0.4),
               dexp(3.1, 0.4, log = TRUE))
  expect_equal(distr_exponential_mean$density(3.1, 0.4),
               dexp(3.1, 1 / 0.4, log = TRUE))
})


test_that("can sample from exponential distribution", {
  r1 <- mcstate_rng$new(1)
  r2 <- mcstate_rng$new(1)
  expect_equal(distr_exponential_rate$sample(r1, 0.4),
               r2$exponential(1, 0.4))
  expect_equal(distr_exponential_mean$sample(r1, 0.4),
               r2$exponential(1, 1 / 0.4))
})


test_that("can compute density for normal distribution", {
  expect_equal(distr_normal$density(0, 1, 2),
               dnorm(0, 1, 2, log = TRUE))
})


test_that("can sample from normal distribution", {
  r1 <- mcstate_rng$new(1)
  r2 <- mcstate_rng$new(1)
  expect_equal(distr_normal$sample(r1, 1, 2),
               r2$normal(1, 1, 2))
})


test_that("can compute density for uniform distribution", {
  expect_equal(distr_uniform$density(1, 2, 3),
               dunif(1, 2, 3, log = TRUE))
})


test_that("can sample from uniform distribution", {
  r1 <- mcstate_rng$new(1)
  r2 <- mcstate_rng$new(1)
  expect_equal(distr_uniform$sample(r1, 2, 3),
               r2$uniform(1, 2, 3))
})


## Needs special testing because the primary use is currently during
## package build
test_that("can create a distribution object", {
  density <- function(x, a, b) NULL
  sample <- function(rng, a, b) NULL
  d <- distribution("Foo", density, sample)
  expect_equal(
    d,
    list(name = "Foo", args = c("a", "b"), density = density, sample = sample))
})
