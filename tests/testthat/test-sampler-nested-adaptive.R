test_that("can run a simple nested adaptive sampler", {
  set.seed(1)
  ng <- 3
  m <- ex_simple_nested(ng)
  v <- list(base = NULL,
            groups = rep(list(matrix(50)), ng))
  s <- monty_sampler_nested_adaptive(v)
  res <- monty_sample(m, s, 100)

  ## Something we can look for is that acceptances are not equally
  ## shared; I've used a very poor vcv here so that acceptance is
  ## relatively rare and we can see this:
  accept <- t(diff(t(res$pars[, , 1])) != 0)
  expect_setequal(colSums(accept), 0:3)
  expect_setequal(colSums(accept * 2^(seq_len(ng) - 1)), 0:7)
})


test_that("validate vcv inputs on construction of sampler", {
  expect_error(
    monty_sampler_nested_adaptive(NULL),
    "Expected a list for 'initial_vcv'")
  expect_error(
    monty_sampler_nested_adaptive(list()),
    "Expected 'initial_vcv' to have elements 'base' and 'groups'")
  expect_error(
    monty_sampler_nested_adaptive(list(base = TRUE, groups = TRUE)),
    "Expected a matrix or 3d array for 'initial_vcv$base'",
    fixed = TRUE)
  expect_error(
    monty_sampler_nested_adaptive(list(base = NULL, groups = TRUE)),
    "Expected 'initial_vcv$groups' to be a list",
    fixed = TRUE)
  expect_error(
    monty_sampler_nested_adaptive(list(base = NULL, groups = list())),
    "Expected 'initial_vcv$groups' to have at least one element",
    fixed = TRUE)
  expect_error(
    monty_sampler_nested_adaptive(list(base = NULL, groups = list(TRUE))),
    "Expected a matrix or 3d array for 'initial_vcv$groups[1]'",
    fixed = TRUE)

  vcv <- list(base = diag(1), groups = list(diag(2), diag(3)))
  expect_s3_class(monty_sampler_nested_adaptive(vcv), "monty_sampler")
})


test_that("can't use nested sampler with models that are not nested", {
  m <- monty_example("gaussian", diag(3))
  vcv <- list(base = diag(1), groups = list(diag(2), diag(3)))
  s <- monty_sampler_nested_adaptive(vcv)
  expect_error(
    monty_sample(m, s, 100),
    "Your model does not have parameter groupings")
})


test_that("Validate that the model produces correct parameter groupings", {
  m <- monty_model(list(
    parameters = c("a", "b", "c"),
    density = function(x, by_group = TRUE) {
      x
    },
    parameter_groups = c(0, 1, 2)))
  s <- monty_sampler_nested_adaptive(
    list(base = diag(1), groups = list(diag(1), diag(1))))
  expect_error(
    monty_sample(m, s, 100, c(0, 0, 0)),
    "model$density(x, by_group = TRUE) did not produce a density",
    fixed = TRUE)
})


test_that("Validate length of by group output", {
  m <- monty_model(list(
    parameters = c("a", "b", "c"),
    density = function(x, by_group = TRUE) {
      structure(x, by_group = rep(0, 5))
    },
    parameter_groups = c(0, 1, 2)))
  s <- monty_sampler_nested_adaptive(
    list(base = diag(1), groups = list(diag(1), diag(1))))
  expect_error(
    monty_sample(m, s, 100, c(0, 0, 0)),
    "model$density(x, by_group = TRUE) produced a 'by_group' attribute with",
    fixed = TRUE)
})


test_that("can continue nested sampler correctly", {
  set.seed(1)
  ng <- 3
  m <- ex_simple_nested(ng)
  v <- list(base = NULL,
            groups = rep(list(matrix(1)), ng))
  sampler <- monty_sampler_nested_adaptive(v)

  set.seed(1)
  res1 <- monty_sample(m, sampler, 30, n_chains = 3, restartable = TRUE)

  set.seed(1)
  res2a <- monty_sample(m, sampler, 10, n_chains = 3, restartable = TRUE)
  res2b <- monty_sample_continue(res2a, 20, restartable = TRUE)

  expect_equal(res2b, res1)
})


test_that("can run a sampler with shared parameters", {
  set.seed(1)
  ng <- 5
  m <- ex_simple_nested_with_base(ng)
  s <- monty_sampler_nested_adaptive(
    list(base = diag(1), groups = rep(list(diag(1)), ng)))
  res <- monty_sample(m, s, 100)
  ## This is not great, until we get a real example we can use; but it
  ## does test that the sampler runs, and that it mixes.
  expect_true(length(unique(res$pars[1, , ])) > 1)
})


test_that("can run an observer during a nested fit", {
  set.seed(1)
  ng <- 5
  m <- ex_simple_nested_with_base(ng)
  s <- monty_sampler_nested_adaptive(
    list(base = diag(1), groups = rep(list(diag(1)), ng)))
  counter <- 0
  ## Directly wire this in for now; we really just need better
  ## examples here.
  m$observer <- monty_observer(function(...) {
    counter <<- counter + 1
    list(n = counter)
  })
  m$properties$has_observer <- TRUE
  res <- monty_sample(m, s, 100)
  expect_equal(
    dim(res$observations$n),
    c(1, 100, 1))

  pars <- cbind(res$initial[, 1], res$pars[, , 1])
  pars_update <- apply(pars, 1, diff) != 0
  ## n_update: 1 for initial + num of base updates + num of >0 group updates
  n_update <- 1 + sum(pars_update[, "sigma"]) +
    sum(rowSums(pars_update[, paste0("mu_", seq_len(5))]) > 0)
  expect_equal(max(res$observations$n), n_update)
})


test_that("Empirical VCV calculated correctly with forget_rate = 0", {
  set.seed(1)
  ng <- 3
  m <- ex_simple_nested_with_base(ng)
  s <- monty_sampler_nested_adaptive(
    initial_vcv = list(base = diag(1), groups = rep(list(diag(1)), ng)),
    forget_rate = 0,
    log_scaling_update = FALSE)
  res <- monty_sample(m, s, 1000)
  expect_equal(names(res),
               c("pars", "density", "initial", "details", "observations"))

  ## forget_rate = 0 so full chain should be included in VCV
  expect_equal(res$details[[1]]$weight, 1000)
  expect_equal(res$details[[1]]$included, seq_len(1000))
  pars <- t(array_drop(res$pars, 3))
  vcv <- res$details[[1]]$vcv
  expect_equal(vcv$base, var(pars[, "sigma"]), ignore_attr = TRUE)
  expect_equal(vcv$groups[[1]], var(pars[, "mu_1"]), ignore_attr = TRUE)
  expect_equal(vcv$groups[[2]], var(pars[, "mu_2"]), ignore_attr = TRUE)
  expect_equal(vcv$groups[[3]], var(pars[, "mu_3"]), ignore_attr = TRUE)
})


test_that("Empirical VCV calculated correctly with forget_rate = 0.1", {
  set.seed(1)
  ng <- 3
  m <- ex_simple_nested_with_base(ng)
  s <- monty_sampler_nested_adaptive(
    initial_vcv = list(base = diag(1), groups = rep(list(diag(1)), ng)),
    forget_rate = 0.1)
  res <- monty_sample(m, s, 1000)
  expect_equal(names(res),
               c("pars", "density", "initial", "details", "observations"))

  ## forget_rate = 0.1 so VCV should exclude first 100 parameter sets
  expect_equal(res$details[[1]]$weight, 900)
  expect_equal(res$details[[1]]$included, seq(101, 1000, by = 1))
  pars <- t(array_drop(res$pars, 3))
  vcv <- res$details[[1]]$vcv
  expect_equal(vcv$base, var(pars[101:1000, "sigma"]), ignore_attr = TRUE)
  expect_equal(vcv$groups[[1]], var(pars[101:1000, "mu_1"]), ignore_attr = TRUE)
  expect_equal(vcv$groups[[2]], var(pars[101:1000, "mu_2"]), ignore_attr = TRUE)
  expect_equal(vcv$groups[[3]], var(pars[101:1000, "mu_3"]), ignore_attr = TRUE)
})


test_that("Empirical VCV correct using both forget_rate and forget_end", {
  set.seed(1)
  ng <- 3
  m <- ex_simple_nested_with_base(ng)
  s <- monty_sampler_nested_adaptive(
    initial_vcv = list(base = diag(1), groups = rep(list(diag(1)), ng)),
    forget_rate = 0.5,
    forget_end = 200,
    log_scaling_update = FALSE)
  res <- monty_sample(m, s, 1000)
  expect_equal(names(res),
               c("pars", "density", "initial", "details", "observations"))

  ## forget_rate = 0.5 and forget_end = 200 so VCV should exclude first
  ## 100 parameter sets
  expect_equal(res$details[[1]]$weight, 900)
  expect_equal(res$details[[1]]$included, seq(101, 1000, by = 1))
  pars <- t(array_drop(res$pars, 3))
  vcv <- res$details[[1]]$vcv
  expect_equal(vcv$base, var(pars[101:1000, "sigma"]), ignore_attr = TRUE)
  expect_equal(vcv$groups[[1]], var(pars[101:1000, "mu_1"]), ignore_attr = TRUE)
  expect_equal(vcv$groups[[2]], var(pars[101:1000, "mu_2"]), ignore_attr = TRUE)
  expect_equal(vcv$groups[[3]], var(pars[101:1000, "mu_3"]), ignore_attr = TRUE)
})


test_that("Empirical VCV correct using forget_rate, forget_end and adapt_end", {
  set.seed(1)
  ng <- 3
  m <- ex_simple_nested_with_base(ng)
  s <- monty_sampler_nested_adaptive(
    initial_vcv = list(base = diag(1), groups = rep(list(diag(1)), ng)),
    forget_rate = 0.25,
    forget_end = 100,
    adapt_end = 300)
  res <- monty_sample(m, s, 1000)
  expect_equal(names(res),
               c("pars", "density", "initial", "details", "observations"))

  ## forget_rate = 0.25, forget_end = 500 and adapt_end = 300 so VCV should
  ## only include parameter sets 26 to 300
  expect_equal(res$details[[1]]$weight, 275)
  expect_equal(res$details[[1]]$included, seq(26, 300, by = 1))
  pars <- t(array_drop(res$pars, 3))
  vcv <- res$details[[1]]$vcv
  expect_equal(vcv$base, var(pars[26:300, "sigma"]), ignore_attr = TRUE)
  expect_equal(vcv$groups[[1]], var(pars[26:300, "mu_1"]), ignore_attr = TRUE)
  expect_equal(vcv$groups[[2]], var(pars[26:300, "mu_2"]), ignore_attr = TRUE)
  expect_equal(vcv$groups[[3]], var(pars[26:300, "mu_3"]), ignore_attr = TRUE)
})


test_that("check nested adaptive inputs correctly", {
  input <- list(base = 100, groups = list(25, 50, 75))
  expect_equal(check_nested_adaptive(input, 3, TRUE), input)

  input <- list(base = NULL, groups = list(25, 50, 75))
  expect_equal(check_nested_adaptive(input, 3, TRUE, TRUE), input)

  input <- list(base = 100, groups = 50)
  expect_equal(check_nested_adaptive(input, 3, TRUE),
               list(base = 100, groups = rep(list(50), 3)))

  input <- list(base = 100, groups = NULL)
  expect_equal(check_nested_adaptive(input, 3, TRUE, TRUE),
               list(base = 100, groups = rep(list(NULL), 3)))

  input <- list(base = NULL, groups = list(25, 50, 75))
  expect_equal(check_nested_adaptive(input, 3, FALSE, TRUE), input)

  input <- 100
  expect_equal(check_nested_adaptive(input, 3, TRUE),
               list(base = 100, groups = rep(list(100), 3)))

  input <- 100
  expect_equal(check_nested_adaptive(input, 3, FALSE),
               list(base = NULL, groups = rep(list(100), 3)))

  input <- NULL
  expect_equal(check_nested_adaptive(input, 3, TRUE, TRUE),
               list(base = NULL, groups = rep(list(NULL), 3)))

  input <- NULL
  expect_equal(check_nested_adaptive(input, 3, FALSE, TRUE),
               list(base = NULL, groups = rep(list(NULL), 3)))

  input <- c(25, 50, 75)
  expect_error(check_nested_adaptive(input, 3, TRUE, TRUE),
               "Expected a list, single value or NULL for input", fixed = TRUE)

  input <- c(25, 50, 75)
  expect_error(check_nested_adaptive(input, 3, TRUE),
               "Expected a list or single value for input", fixed = TRUE)

  input <- list(25, 50, 75)
  expect_error(
    check_nested_adaptive(input, 3, TRUE),
    "Expected input input as list to have elements 'base' and 'groups'",
    fixed = TRUE)

  input <- list(base = 100, groups = 50)
  expect_error(
    check_nested_adaptive(input, 3, FALSE),
    "Expected input$base to be NULL as there are no base parameters",
    fixed = TRUE)

  input <- list(base = NULL, groups = list(25, 50, 75))
  expect_error(check_nested_adaptive(input, 3, TRUE),
               "Expected single value for input$base", fixed = TRUE)

  input <- list(base = 100, groups = NULL)
  expect_error(check_nested_adaptive(input, 3, TRUE),
               "Expected a list or single value for input$groups", fixed = TRUE)

  input <- list(base = c(50, 100), groups = list(25, 50, 75))
  expect_error(check_nested_adaptive(input, 3, TRUE, TRUE),
               "Expected single value or NULL for input$base", fixed = TRUE)

  input <- list(base = 100, groups = c(25, 50, 75))
  expect_error(
    check_nested_adaptive(input, 3, TRUE, TRUE),
    "Expected a list, single value or NULL for input$groups", fixed = TRUE)

  input <- list(base = 100, groups = list(25, 50))
  expect_error(
    check_nested_adaptive(input, 3, TRUE),
    "Expected input$groups specified as list to have length 3", fixed = TRUE)

  input <- list(base = 100, groups = list(25, NULL, 75))
  expect_error(
    check_nested_adaptive(input, 3, TRUE),
    "Expected a single value for input$groups[[2]]", fixed = TRUE)

  input <- list(base = 100, groups = list(25, 50, c(75, 100)))
  expect_error(
    check_nested_adaptive(input, 3, TRUE, TRUE),
    "Expected a single value or NULL for input$groups[[3]]", fixed = TRUE)
})

test_that("can run nested adaptive sampler simultaneously", {
  set.seed(1)
  ng <- 5
  m <- ex_simple_nested_with_base(ng)
  sampler <- monty_sampler_nested_adaptive(
    list(base = diag(1), groups = rep(list(diag(1)), ng)))
  
  set.seed(1)
  res1 <- monty_sample(m, sampler, 100, n_chains = 3)
  
  set.seed(1)
  runner <- monty_runner_simultaneous()
  res2 <- monty_sample(m, sampler, 100, n_chains = 3, runner = runner)
  expect_equal(res1, res2)
})


test_that("can run nested adaptive sampler with rejecting boundaries
          simultaneously", {
  set.seed(1)
  ng <- 5
  m <- ex_simple_nested_with_base(ng)
  m$domain[, 1] <- -3
  m$domain[, 2] <- 3
  sampler <- monty_sampler_nested_adaptive(
    list(base = diag(1), groups = rep(list(diag(1)), ng)),
    boundaries = "reject")
  
  set.seed(1)
  res1 <- monty_sample(m, sampler, 100, n_chains = 3)
  
  set.seed(1)
  runner <- monty_runner_simultaneous()
  res2 <- monty_sample(m, sampler, 100, n_chains = 3, runner = runner)
  expect_equal(res1, res2)
})
