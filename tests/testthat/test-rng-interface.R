test_that("Can create pointer object", {
  obj <- monty_rng_pointer$new()
  expect_true(obj$is_current())
  expect_type(obj$state(), "raw")
  expect_length(obj$state(), 32)
  expect_equal(obj$algorithm, "xoshiro256plus")
  r <- obj$state()
  obj$sync()
  expect_identical(r, obj$state())
})


test_that("Using object requires sync", {
  obj <- monty_rng_pointer$new()
  expect_true(obj$is_current())
  r <- obj$state()
  test_xoshiro_run(obj)
  expect_false(obj$is_current())
  obj$sync()
  expect_true(obj$is_current())
  expect_false(identical(obj$state(), r))
})


test_that("Fetching state syncs", {
  obj <- monty_rng_pointer$new()
  expect_true(obj$is_current())
  r <- obj$state()
  test_xoshiro_run(obj)
  expect_false(obj$is_current())
  expect_false(identical(obj$state(), r))
  expect_true(obj$is_current())
})


test_that("Invalidated pointers can be rebuilt", {
  obj1 <- monty_rng_pointer$new()
  obj2 <- corrupt_pointer(obj1)
  expect_equal(
    test_xoshiro_run(obj2),
    test_xoshiro_run(obj1))

  ## This saves that the pointer is invalid:
  obj3 <- corrupt_pointer(obj2)
  expect_error(
    test_xoshiro_run(obj3),
    "Can't unserialise an rng pointer that was not synced")

  ## But if we sync things it's ok:
  obj2$sync()
  obj4 <- corrupt_pointer(obj2)
  expect_equal(
    test_xoshiro_run(obj4),
    test_xoshiro_run(obj1))
})


test_that("can't create invalid pointer types", {
  expect_error(
    monty_rng_pointer$new(algorithm = "mt19937"),
    "Unknown algorithm 'mt19937'")
})


test_that("Validate pointers on fetch", {
  obj <- monty_rng_pointer$new(algorithm = "xoshiro256starstar")
  expect_error(
    test_rng_pointer_get(obj, 1),
    "Incorrect rng type: given xoshiro256starstar, expected xoshiro256plus")
  obj <- monty_rng_pointer$new(algorithm = "xoshiro256plus", n_streams = 4)
  expect_error(
    test_rng_pointer_get(obj, 20),
    "Requested a rng with 20 streams but only have 4")
  expect_silent(
    test_rng_pointer_get(obj, 0))
  expect_silent(
    test_rng_pointer_get(obj, 1))
})


test_that("Create pointer with a long jump", {
  s0 <- monty_rng_pointer$new(1, 4, 0)$state()
  s1 <- monty_rng_pointer$new(1, 4, 1)$state()
  s2 <- monty_rng_pointer$new(1, 4, 2)$state()

  cmp <- monty_rng$new(1, 4)
  expect_equal(s0, cmp$state())
  expect_equal(s1, cmp$long_jump()$state())
  expect_equal(s2, cmp$long_jump()$state())
})


test_that("can summarise errors", {
  r <- monty_rng$new(n_streams = 10)
  err <- expect_error(
    r$binomial(1, 1, -1),
    "10 generators reported errors")
  expect_match(err$message, "and 6 more")
})
