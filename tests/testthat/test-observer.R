test_that("can finalise observations automatically where similar", {
  expect_equal(
    observer_finalise_auto(list(list(a = 1:5), list(a = 6:10))),
    list(a = cbind(1:5, 6:10)))
  expect_equal(
    observer_finalise_auto(list(list(a = 1:5, b = 1:3),
                                list(a = 6:10, b = 4:6))),
    list(a = cbind(1:5, 6:10),
         b = cbind(1:3, 4:6)))
})


test_that("don't combine things that can't be stuck together", {
  expect_equal(
    observer_finalise_auto(list(list(a = 1:5, b = 1:3), list(a = 6:10))),
    list(list(a = 1:5, b = 1:3), list(a = 6:10)))
  expect_equal(
    observer_finalise_auto(list(list(a = 1:5, a = 1), list(a = 6:10, a = 1))),
    list(list(a = 1:5, a = 1), list(a = 6:10, a = 1)))
  expect_equal(
    observer_finalise_auto(list(list(1:5, 1), list(6:10, 1))),
    list(list(1:5, 1), list(6:10, 1)))

  expect_equal(
    observer_finalise_auto(list(list(a = 1:5, b = 1:3),
                                list(a = 6:10, b = NULL))),
    list(a = cbind(1:5, 6:10), b = list(1:3, NULL)))
})


test_that("can append two observations automatically", {
  expect_equal(
    observer_append_auto(list(a = cbind(1:5)), list(a = cbind(6:8))),
    list(a = cbind(1:8)))
})


test_that("don't append incompatible observations", {
  expect_equal(
    observer_append_auto(list(a = cbind(1:5)),
                         list(a = cbind(6:8), b = TRUE)),
    list(list(a = cbind(1:5)), list(a = cbind(6:8), b = TRUE)))

  expect_equal(
    observer_append_auto(list(a = cbind(1:5)), list(a = cbind(6:8, 9:11))),
    list(a = list(cbind(1:5), cbind(6:8, 9:11))))
})
