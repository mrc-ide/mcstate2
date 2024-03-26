test_that("parameter groups validation catches error cases", {
  expect_error(check_parameter_groups(TRUE, 5, "groups"),
               "Expected 'groups' to be integer-like")
  expect_error(check_parameter_groups(4.5, 5, "groups"),
               "Expected 'groups' to be integer-like")
  expect_error(check_parameter_groups(1:4, 5, "groups"),
               "Expected 'groups' to have length 5, but it had length 4")
  expect_error(check_parameter_groups(c(1, 2, 5, 6), 4, "groups"),
               "Missing groups from 'groups'")
  expect_error(check_parameter_groups(c(-2, 1, 2, 3), 4, "groups"),
               "Invalid negative group in 'groups'")
})
