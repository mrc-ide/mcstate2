test_that("repeated calls rmvnorm agrees with one-shot", {
  vcv <- matrix(c(4, 2, 2, 3), ncol = 2)
  x <- runif(4)
  g1 <- mcstate_rng$new(seed = 42)
  g2 <- mcstate_rng$new(seed = 42)
  r1 <- rmvnorm(x, vcv, g1)
  r2 <- make_rmvnorm(vcv)(x, g2)
  expect_identical(r2, r1)
})


test_that("rmvnorm has correct behaviour in trivial case", {
  vcv <- matrix(1, 1, 1)
  g1 <- mcstate_rng$new(seed = 42)
  g2 <- mcstate_rng$new(seed = 42)
  r1 <- replicate(100, rmvnorm(0, vcv, g1))
  r2 <- g2$random_normal(100)
  expect_identical(r1, r2)
})


## This is still just too slow to really push around:
test_that("rvnorm produces correct expectation", {
  vcv <- matrix(c(4, 2, 2, 3), ncol = 2)
  f <- make_rmvnorm(vcv)
  rng <- mcstate_rng$new(seed = 42)
  x <- c(1, 2)
  r <- t(replicate(1e5, f(x, rng)))
  expect_equal(colMeans(r), c(1, 2), tolerance = 0.01)
  expect_equal(var(r), vcv, tolerance = 0.01)
})


test_that("can compute the log density of a multivariate normal", {
  skip_if_not_installed("mvtnorm")
  vcv <- matrix(c(4, 2, 2, 3), ncol = 2)
  set.seed(1)
  x <- matrix(rnorm(20, sd = 10), ncol = 2)
  expect_equal(apply(x, 1, make_ldmvnorm(vcv)),
               mvtnorm::dmvnorm(x, sigma = vcv, log = TRUE))
  expect_identical(apply(x, 1, make_ldmvnorm(vcv)),
                   apply(x, 1, ldmvnorm, vcv))
})


test_that("can compute derivatives of multivariate normal log density", {
  skip_if_not_installed("numDeriv")
  vcv <- matrix(c(4, 2, 2, 3), ncol = 2)
  set.seed(1)
  x <- matrix(rnorm(20, sd = 10), ncol = 2)
  f <- make_ldmvnorm(vcv)
  g <- make_deriv_ldmvnorm(vcv)
  expect_equal(apply(x, 1, g),
               apply(x, 1, function(xi) numDeriv::grad(f, xi)))
  expect_identical(apply(x, 1, g),
                   apply(x, 1, deriv_ldmvnorm, vcv))
})


test_that("Can draw samples from many single-variable MVNs", {
  r1 <- mcstate_rng$new(seed = 42, n_streams = 4)
  r2 <- mcstate_rng$new(seed = 42, n_streams = 4)
  vcv <- array(1, c(1, 1, 4))
  x <- runif(4)
  expect_equal(rmvnorm(x, vcv, r2), r1$random_normal(1) + x)
  expect_equal(rmvnorm(x, 0.1 * vcv, r2), r1$random_normal(1) * sqrt(0.1) + x)
  expect_equal(rmvnorm(x, 0.1 * 1:4 * vcv, r2),
               r1$random_normal(1) * sqrt(0.1 * 1:4) + x)
})


test_that("Can draw samples from many centred single-variable MVNs", {
  r1 <- mcstate_rng$new(seed = 42, n_streams = 4)
  r2 <- mcstate_rng$new(seed = 42, n_streams = 4)
  vcv <- array(1, c(1, 1, 4))
  x <- runif(4)
  expect_equal(make_rmvnorm(vcv, centred = TRUE)(r2),
               r1$random_normal(1))

  expect_equal(make_rmvnorm(vcv, centred = TRUE)(r2),
               r1$random_normal(1))
  expect_equal(make_rmvnorm(0.1 * vcv, centred = TRUE)(r2),
               r1$random_normal(1) * sqrt(0.1))
  expect_equal(make_rmvnorm(0.1 * 1:4 * vcv, centred = TRUE)(r2),
               r1$random_normal(1) * sqrt(0.1 * 1:4))
})


test_that("Can draw samples from many bivariate MVNs", {
  r1 <- mcstate_rng$new(seed = 42, n_streams = 5)
  r2 <- mcstate_rng$new(seed = 42, n_streams = 5)
  r3 <- mcstate_rng$new(seed = 42, n_streams = 5)
  vcv <- array(0, c(2, 2, 5))
  set.seed(1)
  vcv[1, 1, ] <- 1:5
  vcv[2, 2, ] <- 1
  vcv[1, 2, ] <- vcv[2, 1, ] <- rnorm(5, 0, 0.1)

  x <- matrix(runif(10), 2, 5)
  y <- rmvnorm(x, vcv, r2)
  expect_equal(dim(y), dim(x))

  ## A bit of work do do these separately:
  rng_state <- matrix(r1$state(), ncol = 5)
  z <- vapply(1:5, function(i) {
    rmvnorm(x[, i], vcv[, , i], mcstate_rng$new(rng_state[, i]))
  }, numeric(2))
  expect_identical(y, z)

  expect_equal(make_rmvnorm(vcv, centred = TRUE)(r3),
               y - x)
})
