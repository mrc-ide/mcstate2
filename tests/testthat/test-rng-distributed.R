test_that("Can create a set of distributed random number states", {
  s <- monty_rng_distributed_state(1L, n_nodes = 2)
  expect_type(s, "list")
  expect_length(s, 2)
  cmp <- monty_rng$new(1)
  expect_equal(s[[1]], cmp$state())
  cmp$long_jump()
  expect_equal(s[[2]], cmp$state())
})


test_that("Can create a set of distributed random number states", {
  s <- monty_rng_distributed_state(1L, n_streams = 3, n_nodes = 2)
  expect_type(s, "list")
  expect_length(s, 2)
  cmp <- monty_rng$new(1, 3)
  expect_equal(s[[1]], cmp$state())
  cmp$long_jump()
  expect_equal(s[[2]], cmp$state())
})


test_that("can create distributed rng pointers", {
  algorithm <- "xoroshiro128starstar"
  p <- monty_rng_distributed_pointer(1L, n_streams = 3, n_nodes = 2,
                                    algorithm = algorithm)
  expect_type(p, "list")
  expect_length(p, 2)
  expect_equal(
    p[[1]]$state(),
    monty_rng_pointer$new(1, 3, 0, algorithm = algorithm)$state())
  expect_equal(
    p[[2]]$state(),
    monty_rng_pointer$new(1, 3, 1, algorithm = algorithm)$state())
})
