test_that("can create a minimal model", {
  m <- mcstate_model(list(density = function(x) dnorm(x, log = TRUE),
                          parameters = "a"))
  expect_s3_class(m, "mcstate_model")
  expect_equal(m$properties,
               mcstate_model_properties(has_gradient = FALSE,
                                        has_direct_sample = FALSE,
                                        is_stochastic = FALSE,
                                        has_parameter_groups = FALSE))
  expect_equal(m$domain, rbind(a = c(-Inf, Inf)))
  expect_equal(m$parameters, "a")
  expect_equal(m$density(0), dnorm(0, log = TRUE))
})


test_that("can create a more interesting model", {
  m <- ex_simple_gamma1()
  expect_equal(m$properties,
               mcstate_model_properties(has_gradient = TRUE,
                                        has_direct_sample = TRUE,
                                        is_stochastic = FALSE,
                                        has_parameter_groups = FALSE))
  expect_equal(m$domain, rbind(gamma = c(0, Inf)))
  expect_equal(m$parameters, "gamma")
  expect_equal(m$density(1), dgamma(1, 1, 1, log = TRUE))
})


test_that("Require parameters to be given if model does not provide them", {
  expect_error(
    mcstate_model(list(density = function(x) dnorm(x, log = TRUE))),
    "Expected 'model$parameters' to be a character vector",
    fixed = TRUE)
})


test_that("require density is a function", {
  expect_error(mcstate_model(list(parameters = "a")),
               "Expected 'model$density' to be a function",
               fixed = TRUE)
})


test_that("require gradient is a function if given", {
  expect_error(
    mcstate_model(list(density = identity, gradient = TRUE, parameters = "a")),
    "Expected 'model$gradient' to be a function if non-NULL",
    fixed = TRUE)
})


test_that("require direct sample is a function if given", {
  expect_error(
    mcstate_model(list(density = identity,
                       direct_sample = TRUE,
                       parameters = "a")),
    "Expected 'model$direct_sample' to be a function if non-NULL",
    fixed = TRUE)
})


test_that("validate domain", {
  expect_error(
    mcstate_model(list(density = identity, parameters = "a", domain = list())),
    "Expected 'model$domain' to be a matrix if non-NULL",
    fixed = TRUE)
  expect_error(
    mcstate_model(list(density = identity,
                       parameters = "a",
                       domain = matrix(1:4, 2, 2))),
    "Expected 'model$domain' to have 1 row, but it had 2",
    fixed = TRUE)
  expect_error(
    mcstate_model(list(density = identity,
                       parameters = c("a", "b", "c"),
                       domain = matrix(1:4, 2, 2))),
    "Expected 'model$domain' to have 3 rows, but it had 2",
    fixed = TRUE)
  expect_error(
    mcstate_model(list(density = identity,
                       parameters = "a",
                       domain = matrix(1:4, 1, 4))),
    "Expected 'model$domain' to have 2 columns, but it had 4",
    fixed = TRUE)
  expect_error(
    mcstate_model(list(density = identity,
                       parameters = "a",
                       domain = rbind(A = c(0, 1)))),
    "Unexpected rownames on domain",
    fixed = TRUE)
})


test_that("can use properties to guarantee that a property exists", {
  m <- list(density = function(x) dnorm(x, log = TRUE), parameters = "a")
  expect_error(
    mcstate_model(m, mcstate_model_properties(has_gradient = TRUE)),
    "Did not find a function 'gradient' within your model")
  expect_error(
    mcstate_model(m, mcstate_model_properties(has_direct_sample = TRUE)),
    "Did not find a function 'direct_sample' within your model")
  expect_no_error(
    mcstate_model(m, mcstate_model_properties(has_gradient = FALSE,
                                              has_direct_sample = FALSE)))
})


test_that("can use properties to ignore properties in a model", {
  m <- mcstate_model(
    list(density = identity, parameters = "a", gradient = TRUE),
    mcstate_model_properties(has_gradient = FALSE))
  expect_false(m$properties$has_gradient)
  expect_null(m$gradient)

  m <- mcstate_model(
    list(density = identity, parameters = "a", direct_sample = TRUE),
    mcstate_model_properties(has_direct_sample = FALSE))
  expect_false(m$properties$has_direct_sample)
  expect_null(m$direct_sample)
})


test_that("require properties are correct type", {
  expect_error(
    mcstate_model(
      list(density = identity, parameters = "a", gradient = TRUE),
      list(has_gradient = FALSE)),
    "Expected 'properties' to be a 'mcstate_model_properties' object")
})


test_that("stochastic models need an rng setting function", {
  expect_error(
    mcstate_model(
      list(density = identity, parameters = "a"),
      mcstate_model_properties(is_stochastic = TRUE)),
    "Expected 'model$set_rng_state' and 'model$get_rng_state' to be functions",
    fixed = TRUE)
  expect_error(
    mcstate_model(
      list(density = identity, parameters = "a", set_rng_state = identity),
      mcstate_model_properties(is_stochastic = TRUE)),
    "Expected 'model$get_rng_state' to be a function",
    fixed = TRUE)
  expect_error(
    mcstate_model(
      list(density = identity, parameters = "a", get_rng_state = identity),
      mcstate_model_properties(is_stochastic = TRUE)),
    "Expected 'model$set_rng_state' to be a function",
    fixed = TRUE)

  m <- list(density = identity, parameters = "a", set_rng_state = TRUE,
            get_rng_state = identity)
  expect_error(
    mcstate_model(m),
    "Expected 'model$set_rng_state' to be a function",
    fixed = TRUE)
  expect_error(
    mcstate_model(m,
                  mcstate_model_properties(is_stochastic = TRUE)),
    "Expected 'model$set_rng_state' to be a function",
    fixed = TRUE)
  expect_no_error(
    res <- mcstate_model(m,
                         mcstate_model_properties(is_stochastic = FALSE)))
  expect_null(res$set_rng_state)
})


test_that("ignore groups if requested", {
  m <- mcstate_model(
    list(density = identity, parameters = "x"),
    mcstate_model_properties(has_parameter_groups = FALSE))
  expect_false(m$properties$has_parameter_groups)
  expect_null(m$parameter_groups)

  m <- mcstate_model(
    list(density = identity, parameters = "x", parameter_groups = 1),
    mcstate_model_properties(has_parameter_groups = FALSE))
  expect_false(m$properties$has_parameter_groups)
  expect_null(m$parameter_groups)
})


test_that("If parameter groups are present, then density requires arg", {
  expect_error(
    mcstate_model(
      list(density = identity, parameters = "x", parameter_groups = 1),
      mcstate_model_properties(has_parameter_groups = TRUE)),
    "Expected 'model$density' to have an argument 'by_group'",
    fixed = TRUE)
})
