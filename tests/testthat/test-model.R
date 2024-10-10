test_that("can create a minimal model", {
  m <- monty_model(list(density = function(x) dnorm(x, log = TRUE),
                        parameters = "a"))
  expect_s3_class(m, "monty_model")
  expect_equal(m$properties,
               monty_model_properties(has_gradient = FALSE,
                                      has_direct_sample = FALSE,
                                      is_stochastic = FALSE,
                                      has_observer = FALSE,
                                      has_parameter_groups = FALSE))
  expect_equal(m$domain, rbind(a = c(-Inf, Inf)))
  expect_equal(m$parameters, "a")
  expect_equal(m$density(0), dnorm(0, log = TRUE))
})


test_that("can create a more interesting model", {
  m <- ex_simple_gamma1()
  expect_equal(m$properties,
               monty_model_properties(has_gradient = TRUE,
                                      has_direct_sample = TRUE,
                                      is_stochastic = FALSE,
                                      has_parameter_groups = FALSE,
                                      has_observer = FALSE,
                                      allow_multiple_parameters = TRUE))
  expect_equal(m$domain, rbind(gamma = c(0, Inf)))
  expect_equal(m$parameters, "gamma")
  expect_equal(m$density(1), dgamma(1, 1, 1, log = TRUE))
})


test_that("Require parameters to be given if model does not provide them", {
  expect_error(
    monty_model(list(density = function(x) dnorm(x, log = TRUE))),
    "Expected 'model$parameters' to be a character vector",
    fixed = TRUE)
})


test_that("require density is a function", {
  expect_error(monty_model(list(parameters = "a")),
               "Expected 'model$density' to be a function",
               fixed = TRUE)
})


test_that("require gradient is a function if given", {
  expect_error(
    monty_model(list(density = identity, gradient = TRUE, parameters = "a")),
    "Expected 'model$gradient' to be a function if non-NULL",
    fixed = TRUE)
})


test_that("require direct sample is a function if given", {
  expect_error(
    monty_model(list(density = identity,
                     direct_sample = TRUE,
                     parameters = "a")),
    "Expected 'model$direct_sample' to be a function if non-NULL",
    fixed = TRUE)
})


test_that("require observer is a monty_observer if given", {
  expect_error(
    monty_model(list(density = identity,
                     observer = TRUE,
                     parameters = "a")),
    "Expected 'model$observer' to be a 'monty_observer' object if non-NULL",
    fixed = TRUE)
})


test_that("observer must exist if required", {
  properties <- monty_model_properties(has_observer = TRUE)
  expect_error(
    monty_model(list(parameters = "a", density = identity),
                properties = properties),
    "Did not find a 'monty_observer' object 'observer' within your model")
})


test_that("Ignore invalid observer if properties say to ignore it", {
  properties <- monty_model_properties(has_observer = FALSE)
  res <- monty_model(list(density = identity,
                          observer = TRUE,
                          parameters = "a"),
                     properties = properties)
  expect_s3_class(res, "monty_model")
  expect_false(res$properties$has_observer)
  expect_null(res$observer)
})


test_that("validate domain", {
  expect_error(
    monty_model(list(density = identity, parameters = "a", domain = list())),
    "Expected 'model$domain' to be a matrix if non-NULL",
    fixed = TRUE)
  expect_error(
    monty_model(list(density = identity,
                     parameters = "a",
                     domain = matrix(1:4, 2, 2))),
    "Expected 'model$domain' to have 1 row, but it had 2",
    fixed = TRUE)
  expect_error(
    monty_model(list(density = identity,
                     parameters = c("a", "b", "c"),
                     domain = matrix(1:4, 2, 2))),
    "Expected 'model$domain' to have 3 rows, but it had 2",
    fixed = TRUE)
  expect_error(
    monty_model(list(density = identity,
                     parameters = "a",
                     domain = matrix(1:4, 1, 4))),
    "Expected 'model$domain' to have 2 columns, but it had 4",
    fixed = TRUE)
  expect_error(
    monty_model(list(density = identity,
                     parameters = "a",
                     domain = rbind(A = c(0, 1)))),
    "Unexpected parameters found in 'model$domain' rownames",
    fixed = TRUE)
})


test_that("fill in default domain entries where only some provided", {
  base <- list(density = sum, parameters = letters[1:3])
  default <- matrix(c(-Inf, Inf), 3, 2, byrow = TRUE,
                    dimnames = list(letters[1:3], NULL))
  expect_equal(
    monty_model(modifyList(base, list(domain = NULL)))$domain,
    default)
  expect_equal(
    monty_model(modifyList(base, list(domain = rbind(b = c(0, 1)))))$domain,
    rbind(a = c(-Inf, Inf), b = c(0, 1), c = c(-Inf, Inf)))
  change <- rbind(b = c(0, 1), a = c(-1, 2))
  expect_equal(
    monty_model(modifyList(base, list(domain = change)))$domain,
    rbind(a = c(-1, 2), b = c(0, 1), c = c(-Inf, Inf)))
})

test_that("can use properties to guarantee that a property exists", {
  m <- list(density = function(x) dnorm(x, log = TRUE), parameters = "a")
  expect_error(
    monty_model(m, monty_model_properties(has_gradient = TRUE)),
    "Did not find a function 'gradient' within your model")
  expect_error(
    monty_model(m, monty_model_properties(has_direct_sample = TRUE)),
    "Did not find a function 'direct_sample' within your model")
  expect_no_error(
    monty_model(m, monty_model_properties(has_gradient = FALSE,
                                          has_direct_sample = FALSE)))
})


test_that("can use properties to ignore properties in a model", {
  m <- monty_model(
    list(density = identity, parameters = "a", gradient = TRUE),
    monty_model_properties(has_gradient = FALSE))
  expect_false(m$properties$has_gradient)
  expect_null(m$gradient)

  m <- monty_model(
    list(density = identity, parameters = "a", direct_sample = TRUE),
    monty_model_properties(has_direct_sample = FALSE))
  expect_false(m$properties$has_direct_sample)
  expect_null(m$direct_sample)
})


test_that("require properties are correct type", {
  expect_error(
    monty_model(
      list(density = identity, parameters = "a", gradient = TRUE),
      list(has_gradient = FALSE)),
    "Expected 'properties' to be a 'monty_model_properties' object")
})


test_that("stochastic models need an rng setting function", {
  expect_error(
    monty_model(
      list(density = identity, parameters = "a"),
      monty_model_properties(is_stochastic = TRUE)),
    "Expected 'model$set_rng_state' and 'model$get_rng_state' to be functions",
    fixed = TRUE)
  expect_error(
    monty_model(
      list(density = identity, parameters = "a", set_rng_state = identity),
      monty_model_properties(is_stochastic = TRUE)),
    "Expected 'model$get_rng_state' to be a function",
    fixed = TRUE)
  expect_error(
    monty_model(
      list(density = identity, parameters = "a", get_rng_state = identity),
      monty_model_properties(is_stochastic = TRUE)),
    "Expected 'model$set_rng_state' to be a function",
    fixed = TRUE)

  m <- list(density = identity, parameters = "a", set_rng_state = TRUE,
            get_rng_state = identity)
  expect_error(
    monty_model(m),
    "Expected 'model$set_rng_state' to be a function",
    fixed = TRUE)
  expect_error(
    monty_model(m,
                monty_model_properties(is_stochastic = TRUE)),
    "Expected 'model$set_rng_state' to be a function",
    fixed = TRUE)
  expect_no_error(
    res <- monty_model(m,
                       monty_model_properties(is_stochastic = FALSE)))
  expect_null(res$set_rng_state)
})


test_that("ignore groups if requested", {
  m <- monty_model(
    list(density = identity, parameters = "x"),
    monty_model_properties(has_parameter_groups = FALSE))
  expect_false(m$properties$has_parameter_groups)
  expect_null(m$parameter_groups)

  m <- monty_model(
    list(density = identity, parameters = "x", parameter_groups = 1),
    monty_model_properties(has_parameter_groups = FALSE))
  expect_false(m$properties$has_parameter_groups)
  expect_null(m$parameter_groups)
})


test_that("If parameter groups are present, then density requires arg", {
  expect_error(
    monty_model(
      list(density = identity, parameters = "x", parameter_groups = 1),
      monty_model_properties(has_parameter_groups = TRUE)),
    "Expected 'model$density' to have an argument 'by_group'",
    fixed = TRUE)
})


test_that("can print information about simple models", {
  m <- ex_simple_gamma1()
  res <- evaluate_promise(withVisible(print(m)))
  expect_mapequal(res$result, list(value = m, visible = FALSE))
  expect_match(res$messages, "<monty_model>",
               fixed = TRUE, all = FALSE)
  expect_match(res$messages, "Model has 1 parameter: 'gamma'",
               fixed = TRUE, all = FALSE)
  expect_match(res$messages, "can compute gradients",
               fixed = TRUE, all = FALSE)
  expect_match(res$messages, "can be directly sampled from",
               fixed = TRUE, all = FALSE)
})
