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
    "Expected 'vcv$base' to be a matrix",
    fixed = TRUE)
  expect_error(
    mcstate_sampler_nested_random_walk(list(base = NULL, groups = TRUE)),
    "Expected 'vcv$groups' to be a matrix",
    fixed = TRUE)
  expect_error(
    mcstate_sampler_nested_random_walk(list(base = NULL, groups = list())),
    "Expected 'vcv$groups' to have at least one element",
    fixed = TRUE)
  expect_error(
    mcstate_sampler_nested_random_walk(list(base = NULL, groups = list(TRUE))),
    "Expected 'vcv$groups[1]' to be a matrix",
    fixed = TRUE)

  vcv <- list(base = diag(1), groups = list(diag(2), diag(3)))
  expect_s3_class(mcstate_sampler_nested_random_walk(vcv), "mcstate_sampler")
})
