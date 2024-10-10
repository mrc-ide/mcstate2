test_that("can't create empty packer", {
  expect_error(monty_packer(),
               "Trying to generate an empty packer")
})


test_that("trivial packer", {
  xp <- monty_packer("a")
  expect_equal(xp$parameters, "a")
  expect_equal(xp$unpack(1), list(a = 1))
  expect_equal(xp$unpack(c(a = 1)), list(a = 1))
  expect_error(xp$unpack(c(b = 1)),
               "Incorrect names in input")
  expect_error(xp$unpack(1:2),
               "Incorrect length input; expected 1 but given 2")
  expect_equal(xp$pack(list(a = 1)), 1)
  expect_equal(xp$pack(list(a = 1:2)),
               1:2)

  expect_equal(xp$index(), list(a = 1))
})


test_that("multiple scalar unpacking", {
  xp <- monty_packer(c("a", "b", "c"))
  expect_equal(xp$parameters, c("a", "b", "c"))
  expect_equal(xp$unpack(1:3), list(a = 1, b = 2, c = 3))
  expect_equal(xp$pack(list(a = 1, b = 2, c = 3)), 1:3)
  expect_equal(xp$index(), list(a = 1, b = 2, c = 3))
})


test_that("can bind data into an unpacked list", {
  xp <- monty_packer(c("a", "b"), fixed = list(x = 1:5, y = 10))
  expect_equal(xp$parameters, c("a", "b"))
  expect_equal(xp$unpack(1:2), list(a = 1, b = 2, x = 1:5, y = 10))
  expect_equal(xp$pack(list(a = 1, b = 2, x = 1:5, y = 10)), 1:2)
  expect_equal(xp$pack(list(a = 1, b = 2)), 1:2)
})


test_that("can unpack arrays", {
  xp <- monty_packer("a", list(b = 3))
  expect_equal(xp$parameters, c("a", "b[1]", "b[2]", "b[3]"))
  expect_equal(xp$unpack(1:4), list(a = 1, b = 2:4))
})


test_that("can use integer vectors for array inputs", {
  xp <- monty_packer("a", c(b = 3, c = 4))
  expect_equal(xp$parameters,
               c("a", sprintf("b[%d]", 1:3), sprintf("c[%d]", 1:4)))
  expect_equal(xp$unpack(1:8), list(a = 1, b = 2:4, c = 5:8))
  expect_equal(xp$index(), list(a = 1, b = 2:4, c = 5:8))
})


test_that("can create packer with only arrays", {
  xp <- monty_packer(array = list(a = 2, b = 3))
  expect_equal(xp$unpack(1:5), list(a = 1:2, b = 3:5))
})


test_that("can create packers with higher-level dimensionsality", {
  xp <- monty_packer(
    array = list(a = 1, b = 2, c = 2:3, d = 2:4))
  expect_equal(
    xp$parameters,
    c("a[1]", "b[1]", "b[2]",
      sprintf("c[%d,%d]", 1:2, rep(1:3, each = 2)),
      sprintf("d[%d,%d,%d]", 1:2, rep(1:3, each = 2), rep(1:4, each = 6))))
  expect_equal(
    xp$unpack(1:33),
    list(a = 1, b = 2:3, c = matrix(4:9, 2, 3), d = array(10:33, 2:4)))
  expect_equal(
    xp$pack(xp$unpack(1:33)),
    1:33)
})


test_that("names for scalar, array and fixed must be reasonable", {
  expect_error(
    monty_packer(c("a", "b", "a")),
    "Elements of 'scalar' must be unique")
  expect_error(
    monty_packer("a", list(1)),
    "'array' must be named")
  expect_error(
    monty_packer("a", list(b = 1, c = 2, b = 3)),
    "'array' must have unique names")
  expect_error(
    monty_packer("a", fixed = list(x = 1, y = 2, x = 1)),
    "'fixed' must have unique names")
  expect_error(
    monty_packer(c("a", "b", "c"), list(b = 2, d = 3)),
    "Names must be distinct between 'scalar', 'array' and 'fixed'")
  expect_error(
    monty_packer(c("a", "b", "c"), list(d = 2, e = 3), list(d = 5)),
    "Names must be distinct between 'scalar', 'array' and 'fixed'")
})


test_that("validate array inputs", {
  expect_error(
    monty_packer(array = list(a = 1, b = 2.5)),
    "Elements of 'array' must be integer-like vectors, but 'b' is not")
  expect_error(
    monty_packer(array = list(a = 1, b = c(3, 0, 2))),
    "All dimensions in 'array' must be at least 1, but 'b' violates this")
})


test_that("can pass empty array elements as scalars", {
  p <- monty_packer(array = list(a = integer(), b = 2))
  expect_equal(p$parameters, c("a", "b[1]", "b[2]"))
  expect_equal(p$unpack(1:3), list(a = 1, b = 2:3))
  expect_equal(p$pack(list(a = 1, b = 2:3)), 1:3)
})


test_that("can pass empty array elements as scalars in odd order", {
  p <- monty_packer(array = list(a = 2, b = NULL))
  expect_equal(p$parameters, c("a[1]", "a[2]", "b"))
  expect_equal(p$unpack(1:3), list(a = 1:2, b = 3))
  expect_equal(p$pack(list(a = 1:2, b = 3)), 1:3)
})


test_that("require that scalar is a character vector", {
  expect_error(
    monty_packer(1:2),
    "Expected a character vector for 'scalar'")
})


test_that("can post-process parameters", {
  p <- function(x) {
    list(d = x$a + x$b + x$c)
  }
  xp <- monty_packer(c("a", "b", "c"), process = p)
  expect_equal(xp$parameters, c("a", "b", "c"))
  expect_equal(xp$unpack(1:3), list(a = 1, b = 2, c = 3, d = 6))
  expect_equal(xp$pack(xp$unpack(1:3)), 1:3)
  expect_equal(xp$pack(list(a = 1, b = 2, c = 3)), 1:3)
})


test_that("reject non-function process arguments", {
  expect_error(
    monty_packer(c("a", "b", "c"), process = TRUE),
    "Expected a function for 'process'")
})


test_that("require that process is well-behaved", {
  p <- function(x) {
    c(x, list(d = x$a + x$b + x$c))
  }
  xp <- monty_packer(c("a", "b", "c"), process = p)
  expect_error(xp$unpack(1:3),
               "'process()' is trying to overwrite entries in parameters",
               fixed = TRUE)
})


test_that("Can print a packer", {
  p <- monty_packer(c("x", "y"))
  res <- evaluate_promise(withVisible(print(p)))
  expect_mapequal(res$result, list(value = p, visible = FALSE))
  expect_match(res$messages, "<monty_packer>",
               fixed = TRUE, all = FALSE)
  expect_match(res$messages, "Packing 2 parameters: 'x' and 'y'",
               fixed = TRUE, all = FALSE)
})


test_that("unpack a matrix", {
  p <- monty_packer("x", list(y = 5, z = c(2, 3))) # > 1 + 5 + 6 = 12
  m <- matrix(seq_len(12 * 3), 12)

  res <- p$unpack(m)
  expect_equal(res$x, c(1, 13, 25))
  expect_equal(res$y, cbind(2:6, 14:18, 26:30))
  expect_equal(res$z, array(c(7:12, 19:24, 31:36), c(2, 3, 3)))
})


test_that("error if given the wrong size input to unpack", {
  p <- monty_packer("x", list(y = 5, z = c(2, 3))) # > 1 + 5 + 6 = 12
  m <- matrix(seq_len(12 * 3), 9)
  expect_error(
    p$unpack(m),
    "Incorrect length of first dimension of input; expected 12 but given 9")
})


test_that("error if given the wrong size input to unpack", {
  p <- monty_packer("x", list(y = 5, z = c(2, 3))) # > 1 + 5 + 6 = 12
  m <- matrix(seq_len(12 * 3), 12)
  rownames(m) <- letters[1:12]
  expect_error(
    p$unpack(m),
    "Incorrect rownames in input")
})


test_that("can't used fixed with array unpacking", {
  p <- monty_packer(c("x", "y"), fixed = list(a = 10))
  expect_equal(p$unpack(1:2), list(x = 1, y = 2, a = 10))
  expect_error(
    p$unpack(matrix(1:6, 2)),
    "Can't unpack a matrix where the unpacker uses 'fixed'")
})


test_that("can't used process with array unpacking", {
  p <- monty_packer(c("x", "y"), process = function(d) list(z = d$x + d$y))
  expect_equal(p$unpack(1:2), list(x = 1, y = 2, z = 3))
  expect_error(
    p$unpack(matrix(1:6, 2)),
    "Can't unpack a matrix where the unpacker uses 'process'")
})


test_that("Properly unpack scalars stored as zero-length arrays", {
  p <- monty_packer(array = list(a = integer(0), b = 1L))
  expect_equal(p$unpack(1:2), list(a = 1, b = 2))
  expect_equal(p$unpack(cbind(1:2)), list(a = 1, b = matrix(2)))
  expect_equal(p$unpack(matrix(1:10, 2, 5)),
               list(a = seq(1, 9, by = 2),
                    b = matrix(seq(2, 10, by = 2), 1, 5)))
  expect_equal(p$unpack(array(1:30, c(2, 3, 5))),
               list(a = matrix(seq(1, 29, by = 2), c(3, 5)),
                    b = array(seq(2, 30, by = 2), c(1, 3, 5))))
})


## These tests are the inverse of the tests abvove.
test_that("Roundtrip scalars stored as zero-length arrays", {
  p <- monty_packer(array = list(a = integer(0), b = 1L))
  expect_equal(p$pack(list(a = 1, b = 2)), c(1, 2))

  ## expect_equal(p$unpack(1:2), list(a = 1, b = 2))
  expect_equal(p$pack(list(a = 1, b = 2)), c(1, 2))
  expect_equal(p$pack(list(a = 1, b = matrix(2))), cbind(c(1, 2)))

  ## expect_equal(p$pack(list(a = 1, b = 2)), 1:2)
  expect_equal(p$pack(list(a = seq(1, 9, by = 2),
                           b = matrix(seq(2, 10, by = 2), 1, 5))),
               matrix(1:10, 2, 5))
  expect_equal(p$pack(list(a = matrix(seq(1, 29, by = 2), c(3, 5)),
                           b = array(seq(2, 30, by = 2), c(1, 3, 5)))),
               array(1:30, c(2, 3, 5)))
})


test_that("all-scalar corner case", {
  p1 <- monty_packer(c("a", "b"))
  p2 <- monty_packer(array = list(a = integer(0), b = integer(0)))

  ## These inputs both map to the same output, for both ways of
  ## writing the packer:
  expect_equal(p1$unpack(1:2),
               list(a = 1, b = 2))
  expect_equal(p1$unpack(cbind(1:2)),
               list(a = 1, b = 2))
  expect_equal(p2$unpack(1:2),
               list(a = 1, b = 2))
  expect_equal(p2$unpack(cbind(1:2)),
               list(a = 1, b = 2))

  ## Which means that we can't work out how to pack this output, with
  ## either packer:
  expect_equal(p1$pack(list(a = 1, b = 2)),
               1:2)
  expect_equal(p2$pack(list(a = 1, b = 2)),
               1:2)
})


test_that("validate that we can consistently unpack things", {
  p <- monty_packer("a", list(b = 2, c = 3:4))

  expect_equal(
    p$pack(list(a = 1, b = 2:3, c = matrix(4:15, 3, 4))),
    1:15)

  ## Order does not matter:
  expect_equal(
    p$pack(rev(list(a = 1, b = 2:3, c = matrix(4:15, 3, 4)))),
    1:15)

  v <- p$unpack(matrix(1:45, 15, 3))
  expect_equal(names(v), c("a", "b", "c"))
  expect_equal(v$a, c(1, 16, 31))
  expect_equal(v$b, cbind(2:3, 17:18, 32:33))
  expect_equal(v$c, array(c(4:15, 19:30, 34:45), c(3, 4, 3)))
})


test_that("give errors when input is the wrong shape, from scalar input", {
  p <- monty_packer("a", list(b = 2, c = 3:4))

  ## Check the happy path first:
  expect_equal(
    p$pack(list(a = 1, b = 2:3, c = matrix(4:15, 3, 4))),
    1:15)
  v <- p$unpack(matrix(1:45, 15, 3))
  expect_equal(p$pack(v), matrix(1:45, 15, 3))

  err <- expect_error(
    p$pack(list(a = 1, b = 1:3, c = matrix(4:15, 3, 4))),
    "Incompatible dimensions in input for 'b'")
  expect_match(err$body[[1]], "b: expected <2>, given <3>")
  err <- expect_error(
    p$pack(list(a = numeric(3),
                b = matrix(0, 3, 3),
                c = array(0, c(3, 4, 3)))),
    "Incompatible dimensions in input for 'b'")
  expect_match(err$body[[1]], "b: expected <2>, given <3>")

  err <- expect_error(
    p$pack(list(a = 1, b = 1:3, c = matrix(4:15, 4, 3))),
    "Incompatible dimensions in input for 'b'")
  expect_match(err$body[[1]], "b: expected <2>, given <3>")
  expect_match(err$body[[2]], "c: expected <3, 4>, given <4, 3>")
  err <- expect_error(
    p$pack(list(a = 1,
                b = matrix(0, 3, 3),
                c = array(0, c(4, 3, 3)))),
    "Incompatible dimensions in input for 'b'")
  expect_match(err$body[[1]], "b: expected <2>, given <3>")
  expect_match(err$body[[2]], "c: expected <3, 4>, given <4, 3>")
})


test_that("give errors when input has incorect residual dimension", {
  p <- monty_packer("a", list(b = 2, c = 3:4))

  err <- expect_error(
    p$pack(list(a = 1, b = matrix(2:3, 2, 3), c = matrix(4:15, 3, 4))),
    "Inconsistent residual dimension in inputs")
  expect_equal(
    err$body,
    c(x = "'a', 'c': <...1>",
      x = "'b': <...3>"))

  err <- expect_error(
    p$pack(list(a = numeric(4),
                b = matrix(0, 2, 3),
                c = array(0, c(3, 4, 5)))),
    "Inconsistent residual dimension in inputs")
  expect_equal(
    err$body,
    c(x = "'a': <...4>",
      x = "'b': <...3>",
      x = "'c': <...5>"))
})


test_that("validate names to pack", {
  p <- monty_packer("a", list(b = 2, c = 3:4))
  expect_error(
    p$pack(list(a = 1, b = 2, c = 3, d = 4)),
    "Unexpected element present in input to pack: 'd'")
  expect_error(
    p$pack(list(a = 1, d = 4)),
    "Missing elements from input to pack: 'b' and 'c'")
})


test_that("fixed inputs can be present or absent", {
  p <- monty_packer(c("a", "b"), fixed = list(c = 10, d = 12))
  expect_equal(p$pack(list(a = 1, b = 2)), 1:2)
  expect_equal(p$pack(list(a = 1, b = 2, c = NA)), 1:2)
})


test_that("if process is present ignore extra names", {
  p <- monty_packer(c("a", "b"), process = identity)
  expect_equal(
    p$pack(list(a = 1, b = 2)),
    1:2)
})


test_that("can subset a packer of scalars", {
  p <- monty_packer(c("a", "b", "c", "d"))
  res <- p$subset(c("b", "c"))
  expect_equal(res$index, 2:3)
  expect_equal(res$packer$parameters, c("b", "c"))
  expect_equal(res$packer$unpack(1:2), list(b = 1, c = 2))
})


test_that("can subset a packer of arrays", {
  p <- monty_packer(array = list(a = integer(), b = 2, c = c(3, 3)))
  res <- p$subset(c("a", "c"))
  expect_equal(res$index, c(1, 4:12))

  cmp <- monty_packer(array = list(a = integer(), c = c(3, 3)))
  expect_equal(res$packer$parameters, cmp$parameters)
  expect_equal(res$packer$index(), cmp$index())
})


test_that("can reorder on subset", {
  p <- monty_packer(c("a", "b"), list(c = 2, d = integer(), e = 2, f = 3))
  res <- p$subset(c("d", "a", "c"))
  cmp <- monty_packer(array = list(d = integer(), a = integer(), c = 2))

  expect_equal(res$index, c(5, 1, 3, 4))
  expect_equal(res$packer$index(), cmp$index())
})


test_that("prevent duplicates in subset", {
  p <- monty_packer(c("a", "b"), list(c = 2, d = integer(), e = 2, f = 3))
  expect_error(
    p$subset(c("a", "a", "b", "c")),
    "Duplicated name in 'keep': 'a'")
  expect_error(
    p$subset(c("a", "a", "b", "c", "b")),
    "Duplicated names in 'keep': 'a' and 'b'")
})


test_that("prevent unknown names in subset", {
  p <- monty_packer(c("a", "b"))
  expect_error(
    p$subset(c("a", "b", "c")),
    "Unknown name in 'keep': 'c'")
  expect_error(
    p$subset(c("a", "b", "c", "d")),
    "Unknown names in 'keep': 'c' and 'd'")
})


test_that("don't allow things other than character vectors for now", {
  p <- monty_packer(c("a", "b"))
  expect_error(
    p$subset(NULL),
    "Invalid input for 'keep'; this must currently be a character vector")
  expect_error(
    p$subset(1),
    "Invalid input for 'keep'; this must currently be a character vector")
})
