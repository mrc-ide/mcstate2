##' Describe properties of a model.  Use of this function is optional,
##' but you can pass the return value of this as the `properties`
##' argument of mcstate_model to enforce that your model does actually
##' have these properties.
##'
##' @title Describe model properties
##'
##' @param has_gradient Logical, indicating if the model has a
##'   `gradient` method.  Use `NULL` (the default) to detect this from
##'   the model.
##'
##' @param has_direct_sample Logical, indicating if the model has a
##'   `direct_sample` method.  Use `NULL` (the default) to detect this
##'   from the model.
##'
##' @param is_stochastic Logical, indicating if the model is
##'   stochastic.  Stochastic models must supply `set_rng_state` and
##'   `get_rng_state` methods.
##'
##' @param has_groups Logical, indicating that the model can be
##'   decomposed into parameter groups which are independent of each
##'   other.  This is indicated by using the `parameter_groups`
##'   argument to [mcstate_model], and by the presence of a `by_group`
##'   argument to `density` and (later we may also support this in
##'   `gradient`).  Use `NULL` (the default) to detect this from the
##'   model.
##'
##' @return A list of class `mcstate_model_properties` which should
##'   not be modified.
##'
##' @export
mcstate_model_properties <- function(has_gradient = NULL,
                                     has_direct_sample = NULL,
                                     is_stochastic = NULL,
                                     has_groups = NULL) {
  ret <- list(has_gradient = has_gradient,
              has_direct_sample = has_direct_sample,
              is_stochastic = is_stochastic,
              has_groups = has_groups)
  class(ret) <- "mcstate_model_properties"
  ret
}

##' Create a basic `mcstate` model.  This takes a user-supplied object
##' that minimally can compute a probability density (via a `density`
##' function) and information about parameters; with this we can
##' sample from the model using `MCMC` using [mcstate_sample].  We
##' don't imagine that many users will call this function directly,
##' but that this will be glue used by packages.
##'
##' The `model` argument can be a list or environment (something
##' indexable by `$`) and have elements:
##'
##' * `density`: A function that will compute some probability
##'   density.  It must take an argument representing a parameter
##'   vector (a numeric vector) and return a single value.  This is
##'   the posterior probability density in Bayesian inference, but it
##'   could be anything really.  Models can return `-Inf` if things
##'   are impossible, and we'll try and cope gracefully with that
##'   wherever possible.
##'
##' * `parameters`: A character vector of parameter names.  This
##'   vector is the source of truth for the length of the parameter
##'   vector.
##'
##' * `domain`: Information on the parameter domain.  This is a two
##'   column matrix with `length(parameters)` rows representing each
##'   parameter.  The parameter minimum and maximum bounds are given
##'   as the first and second column.  Infinite values (`-Inf` or
##'   `Inf`) should be used where the parameter has infinite domain up
##'   or down.  Currently used to translate from a bounded to
##'   unbounded space for HMC, but we might also use this for
##'   reflecting proposals in MCMC too.  If not present we assume that
##'   the model is valid everywhere (i.e., that all parameters are
##'   valid from `-Inf` to `Inf`.
##'
##' * `direct_sample`: A function to sample directly from the
##'   parameter space, given an [mcstate_rng] object to sample from.
##'   In the case where a model returns a posterior (e.g., in Bayesian
##'   inference), this is assumed to be sampling from the prior.
##'   We'll use this for generating initial conditions for MCMC where
##'   those are not given, and possibly other uses.  If not given then
##'   when using [mcstate_sample()] the user will have to provide a
##'   vector of initial states.
##'
##' * `gradient`: A function to compute the gradient of `density` with
##'   respect to the parameter vector; takes a parameter vector and
##'   returns a vector the same length.  For efficiency, the model may
##'   want to be stateful so that gradients can be efficiently
##'   calculated after a density calculation, or density after
##'   gradient, where these are called with the same parameters.  This
##'   function is optional (and may not be well defined or possible to
##'   define).
##'
##' * `set_rng_state`: A function to set the state (this is in
##'   contrast to the `rng` that is passed through to `direct_sample`
##'   as that is the _sampler's_ rng stream, but we assume models will
##'   look after their own stream, and that they may need many
##'   streams).  Models that provide this method are assumed to be
##'   stochastic; however, you can use the `is_stochastic` property
##'   (via [mcstate_model_properties()]) to override this (e.g., to
##'   run a stochastic model with its deterministic expectation).
##'   This function takes a raw vector of random number state from
##'   [mcstate_rng] and uses it to set the random number state for
##'   your model; this is derived from the random number stream for a
##'   particular chain, jumped ahead.
##'
##' * `get_rng_state`: A function to get the RNG state; must be
##'   provided if `set_rng_state` is present.  Must return the random
##'   number state, which is a raw vector (potentially quite long).
##'
##' * `parameter_groups`: Optionally, an integer vector indicating
##'   parameter group membership.  The format here may change
##'   (especially if we explore more complex nestings) but at present
##'   parameters with group 0 affect everything (so are accepted or
##'   rejected as a whole), while parameters in groups 1 to `n` are
##'   indepenent (for example, changing the parameters in group 2 does
##'   not affect the density of parameters proposed in group 3).
##'
##' @title Create basic model
##'
##' @param model A list or environment with elements as described in
##'   Details.
##'
##' @param properties Optionally, a [mcstate_model_properties] object,
##'   used to enforce or clarify properties of the model.
##'
##' @return An object of class `mcstate_model`.  This will have elements:
##'
##' * `model`: The model as provided
##' * `parameters`: The parameter name vector
##' * `parameter_groups`: The parameter groups
##' * `domain`: The parameter domain matrix, named with your parameters
##' * `direct_sample`: The `direct_sample` function, if provided by the model
##' * `gradient`: The `gradient` function, if provided by the model
##' * `properties`: A list of properties of the model;
##'   see [mcstate_model_properties()].  Currently this contains:
##'     * `has_gradient`: the model can compute its gradient
##'     * `has_direct_sample`: the model can sample from parameters space
##'     * `is_stochastic`: the model will behave stochastically
##'     * `has_parameter_groups`: The model has separable parameter groups
##'
##' @export
mcstate_model <- function(model, properties = NULL) {
  call <- environment() # for nicer stack traces
  parameters <- validate_model_parameters(model, call)
  domain <- validate_model_domain(model, call)
  density <- validate_model_density(model, call)

  properties <- validate_model_properties(properties, call)
  gradient <- validate_model_gradient(model, properties, call)
  direct_sample <- validate_model_direct_sample(model, properties, call)
  rng_state <- validate_model_rng_state(model, properties, call)
  parameter_groups <- validate_model_parameter_groups(model, properties, call)

  ## Update properties based on what we found:
  properties$has_gradient <- !is.null(gradient)
  properties$has_direct_sample <- !is.null(direct_sample)
  properties$is_stochastic <- !is.null(rng_state$set)
  properties$has_parameter_groups <- !is.null(parameter_groups)

  ret <- list(model = model,
              parameters = parameters,
              parameter_groups = parameter_groups,
              domain = domain,
              density = density,
              gradient = gradient,
              direct_sample = direct_sample,
              rng_state = rng_state,
              properties = properties)
  class(ret) <- "mcstate_model"
  ret
}


validate_model_properties <- function(properties, call = NULL) {
  if (is.null(properties)) {
    return(mcstate_model_properties())
  }
  assert_is(properties, "mcstate_model_properties", call = call)
  properties
}


validate_model_parameters <- function(model, call = NULL) {
  if (!is.character(model$parameters)) {
    cli::cli_abort("Expected 'model$parameters' to be a character vector",
                   arg = "model", call = call)
  }
  model$parameters
}


validate_model_domain <- function(model, call = NULL) {
  n_pars <- length(model$parameters)
  if (is.null(model$domain)) {
    domain <- cbind(rep(-Inf, n_pars), rep(Inf, n_pars))
  } else {
    domain <- model$domain
    if (!is.matrix(domain)) {
      cli::cli_abort("Expected 'model$domain' to be a matrix if non-NULL")
    }
    if (nrow(domain) != n_pars) {
      cli::cli_abort(paste(
        "Expected 'model$domain' to have {n_pars} row{?s},",
        "but it had {nrow(domain)}"))
    }
    if (ncol(domain) != 2) {
      cli::cli_abort(paste(
        "Expected 'model$domain' to have 2 columns,",
        "but it had {ncol(domain)}"))
    }
    if (!is.null(rownames(domain))) {
      if (!identical(rownames(domain), model$parameters)) {
        cli::cli_abort("Unexpected rownames on domain")
      }
    }
  }
  rownames(domain) <- model$parameters
  domain
}


validate_model_density <- function(model, call = NULL) {
  if (!is.function(model$density)) {
    cli::cli_abort("Expected 'model$density' to be a function",
                   arg = "model", call = call)
  }
  model$density
}


validate_model_optional_method <- function(model, properties,
                                           method_name, property_name,
                                           call) {
  if (isFALSE(properties[[property_name]])) {
    return(NULL)
  }
  value <- model[[method_name]]
  if (isTRUE(properties[[property_name]]) && !is.function(value)) {
    cli::cli_abort(
      paste("Did not find a function '{method_name}' within your model, but",
            "your properties say that it should do"),
      arg = "model", call = call)
  }
  if (!is.null(value) && !is.function(value)) {
    cli::cli_abort(
      "Expected 'model${method_name}' to be a function if non-NULL",
      arg = "model", call = call)
  }
  value
}


validate_model_gradient <- function(model, properties, call) {
  validate_model_optional_method(
    model, properties, "gradient", "has_gradient", call)
}


validate_model_direct_sample <- function(model, properties, call) {
  validate_model_optional_method(
    model, properties, "direct_sample", "has_direct_sample", call)
}


validate_model_rng_state <- function(model, properties, call) {
  not_stochastic <- isFALSE(properties$is_stochastic) || (
    is.null(properties$is_stochastic) &&
    is.null(model$set_rng_state) &&
    is.null(model$get_rng_state))
  if (not_stochastic) {
    return(NULL)
  }

  has_set <- is.function(model$set_rng_state)
  has_get <- is.function(model$get_rng_state)

  if (has_set && has_get) {
    return(list(set = model$set_rng_state,
                get = model$get_rng_state))
  }

  ## Everything below here is an error, we just want to make a nice
  ## one.
  if (!has_set && !has_get) {
    msg <- paste("Expected 'model$set_rng_state' and 'model$get_rng_state'",
                 "to be functions")
  } else if (!has_set) {
    msg <- "Expected 'model$set_rng_state' to be a function"
  } else {
    msg <- "Expected 'model$get_rng_state' to be a function"
  }

  if (isTRUE(properties$is_stochastic)) {
    hint <- paste("You have specified 'is_stochastic = TRUE', so in order",
                  "to use your stochastic model we need a way of setting",
                  "its state at the start of sampling and getting it back",
                  "at the end")
  } else {
    hint <- paste("Check the set_rng_state and get_rng_state properties",
                  "of your model, or set the 'is_stochastic' property to",
                  "FALSE to ignore them")
  }
  cli::cli_abort(
    c(msg, i = hint),
    arg = "model", call = call)
}


validate_model_parameter_groups <- function(model, properties, call) {
  if (isFALSE(properties$has_parameter_groups)) {
    return(NULL)
  }
  if (is.null(properties$parameter_groups) && is.null(properties$parameter_groups)) {
    return(NULL)
  }
  ## Here, we're not reporting intent very well; we throw about the
  ## issue but not about how to circumvent it with properties
  check_parameter_groups(model$parameter_groups, length(model$parameters),
                         call)
  if (!("by_group" %in% names(formals(model$density)))) {
    cli::cli_abort("Expected 'model$density' to have an argument 'by_group'",
                   call = call)
  }
  model
}


require_direct_sample <- function(model, message, ...) {
  if (!model$properties$has_direct_sample) {
    cli::cli_abort(
      c(message,
        i = paste("This model does not provide 'direct_sample()', so we",
                  "cannot directly sample from its parameter space")),
      ...)
  }
}


require_deterministic <- function(model, message, ...) {
  if (model$properties$is_stochastic) {
    cli::cli_abort(
      c(message,
        i = paste("This model is stochastic (its 'is_stochastic' property",
                  "is TRUE) so cannot be used in contexts that require",
                  "a deterministic density")),
      ...)
  }
}


require_gradient <- function(model, message, ...) {
  if (!model$properties$has_gradient) {
    cli::cli_abort(
      c(message,
        i = paste("This model does not provide a gradient (its 'has_gradient'",
                  "property is FALSE) so cannot be used in contexts that",
                  "require one")),
      ...)
  }
}
