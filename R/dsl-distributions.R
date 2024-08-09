##' Report information about supported distributions in the dsl.  This
##' is primarily intended for use in packages which use
##' [mcstate_dsl_parse_distribution], as this function reports
##' information about which distributions and arguments would succeed
##' there.
##'
##' @title Information about supported distributions
##'
##' @return A data.frame with columns
##'
##' * `name` the name of the distribution; each name begins with a
##'   capital letter, and there are duplicate names where different
##'   parameterisations are supported.
##' * `args` the arguments of all parameters, *except* the random
##'   variable itself which is given as the first argument to density
##'   functions.
##'
##' We may expand the output here in the future to include information on
##'   if distributions have support in C++, but we might end up
##'   supporting everything this way soon.
##'
##' @export
mcstate_dsl_distributions <- function() {
  dsl_distribution_summary
}

## We will expand this later to support differentiation; that might
## impact the density calculation depending on how we handle the
## symbolic differentiation vs the automatic differentiation.
##
## We also need the domain from the distributions, but the uniform
## makes this surprisingly hard because it is either domain [a, b]
## with constant parameters or it is (-Inf, Inf) otherwise!
distribution <- function(name, density, domain, cpp, expr,
                         sample = NULL, variant = NULL) {
  args <- names(formals(density))[-1]
  list(name = name,
       variant = variant,
       args = args,
       density = density,
       domain = domain,
       sample = sample,
       expr = expr,
       cpp = cpp)
}

distr_beta <- distribution(
  name = "Beta",
  density = function(x, a, b) dbeta(x, a, b, log = TRUE),
  domain = c(0, 1),
  sample = function(rng, a, b) rng$beta(1, a, b),
  expr = list(
    density = quote((a - 1) * log(x) + (b - 1) * log(1 - x) - lbeta(a, b)),
    mean = quote(a / (a + b))),
  cpp = list(density = "beta", sample = "beta"))

distr_binomial <- distribution(
  name = "Binomial",
  density = function(x, size, prob) dbinom(x, size, prob, log = TRUE),
  domain = c(0, Inf), # size?
  sample = function(rng, size, prob) rng$binomial(1, size, prob),
  expr = list(
    density = quote(lchoose(size, x) + x * log(prob) +
                    (size - x) * log(1 - prob)),
    mean = quote(size * prob)),
  cpp = list(density = "binomial", sample = "binomial"))

distr_exponential_rate <- distribution(
  name = "Exponential",
  variant = "rate",
  density = function(x, rate) dexp(x, rate, log = TRUE),
  domain = c(0, Inf),
  sample = function(rng, rate) rng$exponential_rate(1, rate),
  expr = list(
    density = quote(log(rate) - rate * x),
    mean = quote(1 / rate)),
  cpp = list(density = "exponential_rate", sample = "exponential_rate"))

distr_exponential_mean <- distribution(
  name = "Exponential",
  variant = "mean",
  density = function(x, mean) dexp(x, 1 / mean, log = TRUE),
  domain = c(0, Inf),
  sample = function(rng, mean) rng$exponential_mean(1, mean),
  expr = list(
    density = quote(-log(mean) - x / mean),
    mean = quote(mean)),
  cpp = list(density = "exponential_mean", sample = "exponential_mean"))

distr_gamma_rate <- distribution(
  name = "Gamma",
  variant = "rate",
  density = function(x, shape, rate) dgamma(x, shape, rate = rate, log = TRUE),
  domain = c(0, Inf),
  sample = function(rng, shape, rate) rng$gamma_rate(1, shape, rate),
  expr = list(
    density = quote((shape - 1) * log(x) - rate * x -
                    lgamma(shape) + shape * log(rate)),
    mean = quote(shape / rate)),
  cpp = list(density = "gamma_rate", sample = "gamma_rate"))

distr_gamma_scale <- distribution(
  name = "Gamma",
  variant = "rate",
  density = function(x, shape, scale) {
    dgamma(x, shape, scale = scale, log = TRUE)
  },
  domain = c(0, Inf),
  sample = function(rng, shape, scale) rng$gamma_scale(1, shape, scale),
  expr = list(
    density = quote((shape - 1) * log(x) - x / scale -
                    lgamma(shape) - shape * log(scale)),
    mean = quote(shape * scale)),
  cpp = list(density = "gamma_scale", sample = "gamma_scale"))

distr_hypergeometric <- distribution(
  name = "Hypergeometric",
  density = function(x, m, n, k) dhyper(x, n1, n2, k, log = TRUE),
  domain = c(0, Inf), # the true version here is quite complex
  sample = function(rng, n1, n2, k) rng$hypergeometric(1, n1, n2, k),
  expr = list(
    density = quote(lchoose(n1, x) + lchoose(n2, k - x) - lchoose(n1 + n2, k)),
    mean = quote(k * n1 / (n1 + n2))),
  cpp = list(density = "hypergeometric", sample = "hypergeometric"))

distr_normal <- distribution(
  name = "Normal",
  density = function(x, mean, sd) dnorm(x, mean, sd, log = TRUE),
  domain = c(-Inf, Inf),
  expr = list(
    density = quote(-(x - mean)^2 / (2 * sd^2) - log(2 * pi) / 2 - log(sd)),
    mean = quote(mean)),
  sample = function(rng, mean, sd) rng$normal(1, mean, sd),
  cpp = list(density = "normal", sample = "normal"))

distr_poisson <- distribution(
  name = "Poisson",
  density = function(x, lambda) dpois(x, lambda, log = TRUE),
  domain = c(0, Inf),
  expr = list(
    density = quote(x * log(lambda) - lambda - lfactorial(x)),
    mean = quote(lambda)),
  sample = function(rng, lambda) rng$poisson(1, lambda),
  cpp = list(density = "poisson", sample = "poisson"))

distr_uniform <- distribution(
  name = "Uniform",
  density = function(x, min, max) dunif(x, min, max, log = TRUE),
  sample = function(rng, min, max) rng$uniform(1, min, max),
  expr = list(
    density = quote(if (x < min || x > max) -Inf else -log(max - min)),
    mean = quote((max - min) / 2)),
  domain = function(min, max) {
    c(if (is.na(min)) -Inf else min, if (is.na(max)) Inf else max)
  },
  cpp = list(density = "uniform", sample = "uniform"))

dsl_distributions <- local({
  d <- list(
    distr_beta,
    distr_binomial,
    distr_exponential_rate, # preferred form, listed first
    distr_exponential_mean,
    distr_gamma_rate,
    distr_gamma_scale,
    distr_hypergeometric,
    distr_normal,
    distr_poisson,
    distr_uniform)
  split(d, vapply(d, "[[", "", "name"))
})


## Just compute this once, on build:
dsl_distribution_summary <- local({
  dat <- unlist(dsl_distributions, FALSE, FALSE)
  data.frame(
    name = vapply(dat, function(x) x$name, ""),
    args = I(lapply(dat, function(x) names(formals(x$density))[-1])),
    stringsAsFactors = FALSE,
    check.names = FALSE)
})


## Utilities to match calls against lists; this would go in util.R but
## it's really very geared around what we need for validating and
## disambiguating calls to distribution function
match_call <- function(args, candidates) {
  for (i in seq_along(candidates)) {
    res <- match_call_candidate(args, candidates[[i]])
    if (!is.null(res)) {
      return(list(success = TRUE, index = i, args = res))
    }
  }

  ## Everything below here is error handling; we've failed to match
  ## the call.
  nms <- names(args)
  unnamed <- sprintf("<%d>", seq_along(args))
  if (is.null(nms)) {
    given <- paste(unnamed, collapse = ", ")
  } else {
    given <- paste(ifelse(nzchar(args), args, unnamed), collapse = ", ")
  }
  detail <- c("x" = sprintf("Failed to match given arguments: %s", given))

  expected <- vcapply(candidates, paste, collapse = ", ")
  if (length(candidates) == 1) {
    hint <- c("i" = "Call should match:")
  } else {
    hint <- c("i" = "Call should match one of:")
  }
  list(success = FALSE, error = c(detail, hint, set_names(expected, "*")))
}


match_call_candidate <- function(args, candidate) {
  if (length(args) != length(candidate)) {
    return(NULL)
  }
  nms <- names(args)
  if (is.null(nms)) {
    return(seq_along(args))
  }
  i <- nzchar(nms)
  if (anyDuplicated(nms[i])) {
    return(NULL)
  }
  j <- match(nms[i], candidate)
  if (anyNA(j)) {
    return(NULL)
  }
  ## This is the same logic as used in R's argument matching, which slots in
  ## named arguments then assigns the remaining as unnamed arguments.
  n <- length(args)
  ret <- integer(n)
  ret[i] <- j
  ret[seq_len(n)[!i]] <- seq_len(n)[-j]
  ret
}
