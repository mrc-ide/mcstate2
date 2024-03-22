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
