test_that("null-or-value works", {
  expect_equal(1 %||% NULL, 1)
  expect_equal(1 %||% 2, 1)
  expect_equal(NULL %||% NULL, NULL)
  expect_equal(NULL %||% 2, 2)
})


test_that("reject inappropriate vcv matrices", {
  expect_error(check_vcv(1, "m"),
               "Expected 'm' to be a matrix")
  expect_error(check_vcv(matrix(runif(6), 2, 3), "m"),
               "Expected 'm' to be symmetric")
  expect_error(check_vcv(matrix(runif(16), 4, 4), "m"),
               "Expected 'm' to be symmetric")
  m <- var(mtcars[1:3])
  expect_no_error(check_vcv(m))
  diag(m) <- 0
  expect_error(check_vcv(m), "Expected 'm' to be positive definite")
})


test_that("can bind lists of matrices", {
  m1 <- matrix(1:6, 3, 2)
  m2 <- matrix(7:12, 3, 2)
  m3 <- matrix(13:16, 2, 2)

  expect_equal(rbind_list(list(m1)), m1)
  expect_equal(rbind_list(list(m1, m2)), rbind(m1, m2))
  expect_equal(rbind_list(list(m1, m2, m3)), rbind(m1, m2, m3))
})


test_that("can get generalised array dimensions", {
  expect_equal(dim2(1:6), 6)
  expect_equal(dim2(matrix(1:6, 2, 3)), c(2, 3))
})


test_that("can get generalised array names", {
  v <- 1:6
  m <- matrix(v, 2, 3)
  expect_null(dimnames2(v))
  expect_null(dimnames2(m))

  names(v) <- letters[1:6]
  dimnames(m) <- list(letters[1:2], letters[3:5])
  expect_equal(dimnames2(v), list(letters[1:6]))
  expect_equal(dimnames2(m), list(letters[1:2], letters[3:5]))
})


test_that("can test if things are the same", {
  expect_true(all_same(integer()))
  expect_true(all_same(1))
  expect_true(all_same(c(1, 1)))
  expect_true(all_same(rep(1, 10)))
  expect_false(all_same(c(1, 1, 2, 1)))
})
