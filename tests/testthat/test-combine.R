test_that("can combine two simple models", {
  a <- monty_model(list(parameters = "x",
                        density = function(x) dnorm(x, log = TRUE)))
  b <- monty_model(list(parameters = "x",
                        density = function(x) dexp(x, log = TRUE)))
  ab <- a + b
  expect_equal(ab$properties, a$properties)
  expect_equal(ab$density(2), a$density(2) + b$density(2))
})


test_that("can combine a model with direct_sample and one without", {
  a <- monty_model(list(
    parameters = "x",
    density = function(x) dnorm(x, log = TRUE),
    direct_sample = function(rng) rng$random_normal(1)))
  b <- monty_model(list(
    parameters = "x",
    density = function(x) dexp(x, log = TRUE)))
  ab <- a + b
  expect_equal(ab$properties, a$properties)

  r1 <- monty_rng$new(seed = 42)
  r2 <- monty_rng$new(seed = 42)
  expect_equal(ab$direct_sample(r1), a$direct_sample(r2))
})


test_that("can't create direct_sample if both models have method", {
  a <- ex_simple_gamma1()
  b <- ex_simple_gamma1()
  ab <- a + b
  expect_null(ab$direct_sample)
  expect_false(ab$properties$has_direct_sample)
  properties <- monty_model_properties(has_direct_sample = TRUE)
  expect_error(
    monty_model_combine(a, b, properties = properties),
    "Can't create a direct_sample from these models")
})


test_that("can't create direct_sample if neither model has method", {
  a <- monty_model(list(parameters = "x",
                        density = function(x) dnorm(x, log = TRUE)))
  properties <- monty_model_properties(has_direct_sample = TRUE)
  expect_error(
    monty_model_combine(a, a, properties = properties),
    "Can't create a direct_sample from these models")
})


test_that("can combine domains", {
  a <- monty_model(list(parameters = "y",
                        density = identity,
                        domain = rbind(c(-5, 5))))
  b <- monty_model(list(parameters = c("x", "y"),
                        density = identity,
                        domain = rbind(c(-10, 10),
                                       c(-3, 10))))
  ab <- a + b
  expect_equal(ab$domain,
               rbind(y = c(-3, 5), x = c(-10, 10)))
})


test_that("direct sampling may reorder parameters", {
  a <- monty_model(list(parameters = "y",
                        density = function(x) dnorm(x, log = TRUE)))
  b <- monty_model(list(
    parameters = c("x", "y"),
    density = function(x) sum(dnorm(x, sd = c(1, 10), log = TRUE)),
    direct_sample = function(rng) rng$normal(2, 0, c(1, 10))))
  ab <- a + b
  expect_true(ab$properties$has_direct_sample)
  r1 <- monty_rng$new(seed = 1)
  r2 <- monty_rng$new(seed = 1)
  expect_equal(ab$direct_sample(r1),
               b$direct_sample(r2)[2:1]) # reversed, to align parameters
})


test_that("direct sampling requires that only one model is sampleable", {
  a <- ex_simple_gamma1(1)
  b <- ex_simple_gamma1(2)
  err <- expect_error(
    monty_model_combine(a, b,
                        monty_model_properties(has_direct_sample = TRUE)),
    "Can't create a direct_sample from these models")
  expect_match(err$body, "Both models have a 'direct_sample'")
  ab <- monty_model_combine(
    a, b, monty_model_properties(has_direct_sample = FALSE))
  expect_false(ab$properties$has_direct_sample)
})


test_that("direct sampling requires that a model is sampleable", {
  a <- monty_model(list(parameters = "a", density = identity))
  err <- expect_error(
    monty_model_combine(a, a,
                        monty_model_properties(has_direct_sample = TRUE)),
    "Can't create a direct_sample from these models")
  expect_match(err$body, "Neither of your models have 'direct_sample' methods")
})


test_that("require that all parameters can be sampled", {
  x <- monty_model(list(
    parameters = c("a", "b"),
    density = identity))
  y <- monty_model(list(
    parameters = c("a", "b", "c"),
    density = identity,
    direct_sample = identity))
  z <- monty_model(list(
    parameters = c("a", "b", "c", "d"),
    density = identity))
  expect_true((x + y)$properties$has_direct_sample)
  expect_false((x + z)$properties$has_direct_sample)
  expect_error(
    y + z,
    "Can't create a direct_sample from these models as 'lhs' does not")
})


test_that("validate input args", {
  x <- ex_simple_gamma1()
  expect_error(
    monty_model_combine(x, NULL),
    "Expected 'b' to be a 'monty_model'")
  expect_error(
    monty_model_combine(NULL, x),
    "Expected 'a' to be a 'monty_model'")
  expect_error(
    monty_model_combine(x, x, TRUE),
    "Expected 'properties' to be a 'monty_model_properties' object")
  expect_error(
    x + NULL,
    "Addition via '+' is only defined for 'monty_model",
    fixed = TRUE)
})


test_that("can combine models with gradients", {
  a <- ex_simple_gamma1(1)
  b <- ex_simple_gamma1(2)
  ab <- a + b
  expect_true(ab$properties$has_gradient)
  expect_equal(ab$gradient(3),
               a$gradient(3) + b$gradient(3))
})


test_that("both models require gradients to create gradient function", {
  a <- monty_model(list(parameters = "x",
                        density = identity,
                        gradient = identity))
  b <- monty_model(list(parameters = "x",
                        density = function(x) dnorm(x, log = TRUE)))
  expect_false(
    monty_model_combine(a, b)$properties$has_gradient)
  p <- monty_model_properties(has_gradient = FALSE)
  expect_false(
    monty_model_combine(a, b, p)$properties$has_gradient)
  p <- monty_model_properties(has_gradient = TRUE)
  expect_error(
    monty_model_combine(a, b, p),
    "Can't create a gradient from these models")
})


test_that("can combine gradients where parameters do not agree", {
  a <- monty_model(list(
    parameters = c("y", "z"),
    density = identity,
    gradient = sqrt))
  b <- monty_model(list(
    parameters = c("x", "y"),
    density = identity,
    gradient = log))
  ab <- a + b
  expect_true(ab$properties$has_gradient)
  expect_equal(
    ab$gradient(c(2, 3, 4)),
    c(sqrt(2) + log(2), sqrt(3), log(4)))
})


test_that("can combine a stochastic and deterministic model", {
  ll <- ex_dust_sir_likelihood()
  prior <- monty_dsl({
    beta ~ Gamma(shape = 1, rate = 1 / 0.5)
    gamma ~ Gamma(shape = 1, rate = 1 / 0.5)
  })
  post <- ll + prior
  expect_true(post$properties$is_stochastic)

  expect_identical(post$rng_state$get(), ll$rng_state$get())
})


test_that("can't disable stochastic model on combination", {
  ll <- ex_dust_sir_likelihood()
  prior <- monty_dsl({
    beta ~ Gamma(shape = 1, rate = 1 / 0.5)
    gamma ~ Gamma(shape = 1, rate = 1 / 0.5)
  })
  properties <- monty_model_properties(is_stochastic = FALSE)
  expect_error(
    monty_model_combine(ll, prior, properties),
    "Refusing to create non-stochastic model from stochastic components")
})


test_that("can't create model out of two stochastic halves", {
  ll <- ex_dust_sir_likelihood()
  expect_error(
    ll + ll,
    "Can't combine two stochastic models")
})


test_that("Can't force creation of stochastic model from deterministic", {
  a <- monty_model(list(parameters = "x",
                        density = function(x) dnorm(x, log = TRUE)))
  expect_error(
    monty_model_combine(a, a, monty_model_properties(is_stochastic = TRUE)),
    "Can't create stochastic support functions for these models")
})


test_that("cominining models with observers is possible", {
  a <- monty_model(
    list(
      parameters = "x",
      density = identity,
      observer = monty_observer(identity)))
  b <- monty_model(
    list(
      parameters = "x",
      density = identity))
  c <- a + b
  expect_true(c$properties$has_observer)
  expect_identical(c$observer, a$observer)
})


test_that("Can't create observer where both models have them", {
  a <- monty_model(
    list(
      parameters = "x",
      density = identity,
      observer = monty_observer(identity)))
  b <- a + a
  expect_false(b$properties$has_observer)
  expect_null(b$observer)

  properties <- monty_model_properties(has_observer = TRUE)
  err <- expect_error(
    monty_model_combine(a, a, properties = properties),
    "Can't create an observer from these models")
  expect_match(conditionMessage(err),
               "Both models have an 'observer' object")
  err <- expect_error(
    monty_model_combine(b, b, properties = properties),
    "Can't create an observer from these models")
  expect_match(conditionMessage(err),
               "Neither of your models have 'observer' objects")

  properties <- monty_model_properties(has_observer = FALSE)
  res <- monty_model_combine(a, a, properties = properties)
  expect_false(res$properties$has_observer)
})
