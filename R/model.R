##' Describe properties of a model.  Use of this function is optional,
##' but you can pass the return value of this as the `properties`
##' argument of monty_model to enforce that your model does actually
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
##' @param has_parameter_groups Logical, indicating that the model can
##'   be decomposed into parameter groups which are independent of
##'   each other.  This is indicated by using the `parameter_groups`
##'   field within the `model` object passed to [monty_model], and
##'   by the presence of a `by_group` argument to `density` and (later
##'   we may also support this in `gradient`).  Use `NULL` (the
##'   default) to detect this from the model.
##'
##' @param has_observer Logical, indicating if the model has an
##'   "observation" function, which we will describe more fully soon.
##'   An observer is a function `observe` which takes no arguments and
##'   returns arbitrary data about the previously evaluated density.
##'   Use `NULL` (the default) to detect this from the model.
##'
##' @param allow_multiple_parameters Logical, indicating if the
##'   density calculation can support being passed a matrix of
##'   parameters (with each column corresponding to a different
##'   parameter set) and return a vector of densities.  If `FALSE`, we
##'   will support some different approaches to sort this out for you
##'   if this feature is needed.  This cannot be detected from the
##'   model, and the default is `FALSE` because it's not always
##'   straightforward to implement.  However, where it is possible it
##'   may be much more efficient (via vectorisation or
##'   parallelisation) to do this yourself.
##'
##' @return A list of class `monty_model_properties` which should
##'   not be modified.
##'
##' @export
##'
##' @examples
##' # Default properties:
##' monty_model_properties()
##'
##' # Set some properties:
##' monty_model_properties(has_gradient = TRUE, is_stochastic = FALSE)
monty_model_properties <- function(has_gradient = NULL,
                                   has_direct_sample = NULL,
                                   is_stochastic = NULL,
                                   has_parameter_groups = NULL,
                                   has_observer = NULL,
                                   allow_multiple_parameters = FALSE) {
  assert_scalar_logical(has_gradient, allow_null = TRUE)
  assert_scalar_logical(has_direct_sample, allow_null = TRUE)
  assert_scalar_logical(is_stochastic, allow_null = TRUE)
  assert_scalar_logical(has_parameter_groups, allow_null = TRUE)
  assert_scalar_logical(has_observer, allow_null = TRUE)
  assert_scalar_logical(allow_multiple_parameters)

  ret <- list(has_gradient = has_gradient,
              has_direct_sample = has_direct_sample,
              is_stochastic = is_stochastic,
              has_parameter_groups = has_parameter_groups,
              has_observer = has_observer,
              ## TODO: I am not convinced on this name
              allow_multiple_parameters = allow_multiple_parameters)
  class(ret) <- "monty_model_properties"
  ret
}

##' Create a basic `monty` model.  This takes a user-supplied object
##' that minimally can compute a probability density (via a `density`
##' function) and information about parameters; with this we can
##' sample from the model using `MCMC` using [monty_sample].  We
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
##'   wherever possible.  If the property `allow_multiple_parameters`
##'   is `TRUE`, then this function must be able to handle the
##'   argument parameter being a matrix,  and return a vector
##'   of densities.
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
##'   reflecting proposals in MCMC too, as well as a fast way of
##'   avoiding calculating densities where proposals fall out of
##'   bounds.  If not present we assume that the model is valid
##'   everywhere (i.e., that all parameters are valid from `-Inf` to
##'   `Inf`.  If unnamed, you must provide a domain for all
##'   parameters.  If named, then you can provide a subset, with
##'   parameters that are not included assumed to have a domain of
##'   `(-Inf, Inf)`.
##'
##' * `direct_sample`: A function to sample directly from the
##'   parameter space, given a [monty_rng] object to sample from.
##'   In the case where a model returns a posterior (e.g., in Bayesian
##'   inference), this is assumed to be sampling from the prior.
##'   We'll use this for generating initial conditions for MCMC where
##'   those are not given, and possibly other uses.  If not given then
##'   when using [monty_sample()] the user will have to provide a
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
##'   (via [monty_model_properties()]) to override this (e.g., to
##'   run a stochastic model with its deterministic expectation).
##'   This function takes a raw vector of random number state from
##'   [monty_rng] and uses it to set the random number state for
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
##'   independent (for example, changing the parameters in group 2 does
##'   not affect the density of parameters proposed in group 3).
##'
##' @title Create basic model
##'
##' @param model A list or environment with elements as described in
##'   Details.
##'
##' @param properties Optionally, a [monty_model_properties] object,
##'   used to enforce or clarify properties of the model.
##'
##' @return An object of class `monty_model`.  This will have elements:
##'
##' * `model`: The model as provided
##' * `parameters`: The parameter name vector
##' * `parameter_groups`: The parameter groups
##' * `domain`: The parameter domain matrix, named with your parameters
##' * `direct_sample`: The `direct_sample` function, if provided by the model
##' * `gradient`: The `gradient` function, if provided by the model
##' * `properties`: A list of properties of the model;
##'   see [monty_model_properties()].  Currently this contains:
##'     * `has_gradient`: the model can compute its gradient
##'     * `has_direct_sample`: the model can sample from parameters space
##'     * `is_stochastic`: the model will behave stochastically
##'     * `has_parameter_groups`: The model has separable parameter groups
##'
##' @export
##'
##' @seealso [monty_model_function], which provides a simple interface
##'   for creating models from functions.
monty_model <- function(model, properties = NULL) {
  call <- environment() # for nicer stack traces
  parameters <- validate_model_parameters(model, call)
  domain <- validate_model_domain(model, call)
  density <- validate_model_density(model, call)

  properties <- validate_model_properties(properties, call)
  gradient <- validate_model_gradient(model, properties, call)
  direct_sample <- validate_model_direct_sample(model, properties, call)
  observer <- validate_model_observer(model, properties, call)
  rng_state <- validate_model_rng_state(model, properties, call)
  parameter_groups <- validate_model_parameter_groups(model, properties, call)

  ## Update properties based on what we found:
  properties$has_gradient <- !is.null(gradient)
  properties$has_direct_sample <- !is.null(direct_sample)
  properties$is_stochastic <- !is.null(rng_state$set)
  properties$has_parameter_groups <- !is.null(parameter_groups)
  properties$has_observer <- !is.null(observer)

  ret <- list(model = model,
              parameters = parameters,
              parameter_groups = parameter_groups,
              domain = domain,
              density = density,
              gradient = gradient,
              direct_sample = direct_sample,
              observer = observer,
              rng_state = rng_state,
              properties = properties)
  class(ret) <- "monty_model"
  ret
}


##' Compute log density for a model.  This is a wrapper around the
##' `$density` property within a [monty_model] object.
##'
##' @title Compute log density
##'
##' @param model A [monty_model] object
##'
##' @param parameters A vector or matrix of parameters
##'
##' @return A log-density value, or vector of log-density values
##'
##' @seealso [monty_model_gradient] for computing gradients and
##'   [monty_model_direct_sample] for sampling from a model.
##'
##' @export
##' @examples
##' m <- monty_model_function(function(a, b) dnorm(0, a, b, log = TRUE))
##' monty_model_density(m, c(0, 1))
##' monty_model_density(m, c(0, 10))
monty_model_density <- function(model, parameters) {
  require_monty_model(model)
  check_model_parameters(model, parameters)
  model$density(parameters)
}


##' Compute the gradient of log density (which is returned by
##' [monty_model_density]) with respect to parameters.  Not all models
##' support this, and an error will be thrown if it is not possible.
##'
##' @title Compute gradient of log density
##'
##' @inheritParams monty_model_density
##'
##' @param named Logical, indicating if the output should be named
##'   using the parameter names.
##'
##' @return A vector or matrix of gradients
##'
##' @seealso [monty_model_density] for log density, and
##'   [monty_model_direct_sample] to sample from a model
##'
##' @export
##' @examples
##' m <- monty_example("banana")
##' # Global maximum at (0, 0), and gradient is zero there:
##' monty_model_density(m, c(0, 0))
##' monty_model_gradient(m, c(0, 0))
##'
##' # Nonzero gradient away from the origin:
##' monty_model_gradient(m, c(0.4, 0.2))
monty_model_gradient <- function(model, parameters, named = FALSE) {
  require_monty_model(model)
  require_gradient(
    model,
    "Can't compute gradient, as this model does not support it",
    call = environment())
  check_model_parameters(model, parameters)
  assert_scalar_logical(named, call = environment())
  ret <- model$gradient(parameters)
  if (named) {
    if (is.matrix(ret)) {
      rownames(ret) <- model$parameters
    } else {
      names(ret) <- model$parameters
    }
  }
  ret
}


##' Directly sample from a model.  Not all models support this, and an
##' error will be thrown if it is not possible.
##'
##' @title Directly sample from a model
##'
##' @inheritParams monty_model_gradient
##'
##' @param rng Random number state, created by [monty_rng].  Use of
##'   an RNG with more than one stream may or may not work as
##'   expected; this is something we need to tidy up (`mrc-5292`)
##'
##' @return A vector or matrix of sampled parameters
##'
##' @export
##' @examples
##' m <- monty_example("banana")
##'
##' r <- monty_rng$new()
##' monty_model_direct_sample(m, r)
##' monty_model_direct_sample(m, r, named = TRUE)
##'
##' r <- monty_rng$new(n_streams = 3)
##' monty_model_direct_sample(m, r)
##' monty_model_direct_sample(m, r, named = TRUE)
monty_model_direct_sample <- function(model, rng, named = FALSE) {
  require_monty_model(model)
  require_direct_sample(
    model,
    "Can't directly sample from this model")
  assert_scalar_logical(named, call = environment())
  ret <- model$direct_sample(rng)
  if (named) {
    if (is.matrix(ret)) {
      rownames(ret) <- model$parameters
    } else {
      names(ret) <- model$parameters
    }
  }
  ret
}


validate_model_properties <- function(properties, call = NULL) {
  if (is.null(properties)) {
    return(monty_model_properties())
  }
  assert_is(properties, "monty_model_properties", call = call)
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
  domain <- model$domain
  n_pars <- length(model$parameters)

  if (is.null(domain)) {
    domain <- cbind(rep(-Inf, n_pars), rep(Inf, n_pars))
    rownames(domain) <- model$parameters
    return(domain)
  }

  if (!is.matrix(domain)) {
    cli::cli_abort("Expected 'model$domain' to be a matrix if non-NULL",
                   call = call)
  }
  if (ncol(domain) != 2) {
    cli::cli_abort(
      c(paste("Expected 'model$domain' to have 2 columns,",
              "but it had {ncol(domain)}"),
        i = paste("Because your domain is unnamed, if given it must",
                  "include all parameters in the same order as your model")),
      call = call)
  }

  nms <- rownames(domain)
  if (is.null(nms)) {
    if (nrow(domain) != n_pars) {
      cli::cli_abort(
        paste("Expected 'model$domain' to have {n_pars} row{?s},",
              "but it had {nrow(domain)}"),
        call = call)
    }
    rownames(domain) <- model$parameters
  } else {
    ## We might treat parameters that begin with '[' specially and
    ## allow these to replicate.  So if the user has a[1], a[2],
    ## a[3] then a row with 'a' will apply across all of these that
    ## are not explicitly given.
    err <- setdiff(nms, model$parameters)
    if (length(err) > 0) {
      cli::cli_abort(
        c("Unexpected parameters found in 'model$domain' rownames",
          set_names(err, "x")),
        call = call)
    }
    msg <- setdiff(model$parameters, nms)
    if (length(msg) > 0) {
      extra <- cbind(rep(-Inf, length(msg)), rep(Inf, length(msg)))
      rownames(extra) <- msg
      domain <- rbind(domain, extra)
    }
    domain <- domain[model$parameters, , drop = FALSE]
  }

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


validate_model_observer <- function(model, properties, call) {
  if (isFALSE(properties$has_observer)) {
    return(NULL)
  }
  value <- model$observer
  if (isTRUE(properties$has_observer) && !inherits(value, "monty_observer")) {
    cli::cli_abort(
      paste("Did not find a 'monty_observer' object 'observer' within",
            "your model, but your properties say that it should exist"),
      arg = "model", call = call)
  }
  if (!is.null(value) && !inherits(value, "monty_observer")) {
    cli::cli_abort(
      "Expected 'model$observer' to be a 'monty_observer' object if non-NULL",
      arg = "model", call = call)
  }
  value
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
  parameter_groups <- model$parameter_groups
  if (is.null(properties$parameter_groups) && is.null(parameter_groups)) {
    return(NULL)
  }
  ## Here, we're not reporting intent very well; we throw about the
  ## issue but not about how to circumvent it with properties
  check_parameter_groups(parameter_groups, length(model$parameters),
                         name = "model$parameter_groups", call = call)
  if (!("by_group" %in% names(formals(model$density)))) {
    cli::cli_abort("Expected 'model$density' to have an argument 'by_group'",
                   call = call)
  }
  parameter_groups
}


require_monty_model <- function(model, arg = deparse(substitute(model)),
                                call = parent.frame()) {
  if (!inherits(model, "monty_model")) {
    cli::cli_abort("Expected '{arg}' to be a 'monty_model'",
                   arg = arg, call = call)
  }
}


require_direct_sample <- function(model, message, ..., call = parent.frame()) {
  if (!model$properties$has_direct_sample) {
    cli::cli_abort(
      c(message,
        i = paste("This model does not provide 'direct_sample()', so we",
                  "cannot directly sample from its parameter space")),
      ..., call = call)
  }
}


require_deterministic <- function(model, message, ..., call = parent.frame()) {
  if (model$properties$is_stochastic) {
    cli::cli_abort(
      c(message,
        i = paste("This model is stochastic (its 'is_stochastic' property",
                  "is TRUE) so cannot be used in contexts that require",
                  "a deterministic density")),
      ..., call = call)
  }
}


require_gradient <- function(model, message, ..., call = parent.frame()) {
  if (!model$properties$has_gradient) {
    cli::cli_abort(
      c(message,
        i = paste("This model does not provide a gradient (its 'has_gradient'",
                  "property is FALSE) so cannot be used in contexts that",
                  "require one")),
      ..., call = call)
  }
}


require_multiple_parameters <- function(model, message, ...) {
  if (!model$properties$allow_multiple_parameters) {
    cli::cli_abort(
      c(message,
        i = paste("This functionality requires that we can",
                  "provide multiple parameters at once to your model",
                  "and get back a vector of densities, but your model",
                  "does not support this (or does not advertise that it",
                  "does), with the property 'allow_multiple_parameters'",
                  "set to FALSE")),
      ...)
  }
}


##' @export
print.monty_model <- function(x, ...) {
  cli::cli_h1("<monty_model>")
  cli::cli_alert_info(
    "Model has {length(x$parameters)} parameter{?s}: {squote(x$parameters)}")
  ## TODO: once the interface around multiple parameters stabilises,
  ## we should reflect information about allow_multiple_parameters and
  ## has_parameter_groups back here.
  str <- monty_model_properties_str(x$properties)
  if (length(str) > 0) {
    cli::cli_alert_info("This model:")
    cli::cli_bullets(set_names(str, "*"))
  }
  cli::cli_alert_info("See {.help monty_model} for more information")
  invisible(x)
}


monty_model_properties_str <- function(properties) {
  c(if (properties$has_gradient) "can compute gradients",
    if (properties$has_direct_sample) "can be directly sampled from",
    if (properties$is_stochastic) "is stochastic",
    if (properties$has_observer) "has an observer")
}


check_model_parameters <- function(model, parameters, call = parent.frame()) {
  n_pars <- length(model$parameters)
  if (is.matrix(parameters)) {
    require_multiple_parameters(
      model,
      "'parameters' cannot be a matrix",
      arg = "parameters",
      call = call)
    if (nrow(parameters) != n_pars) {
      cli::cli_abort(
        paste("Expected 'parameters' to have {n_pars} row{?s}, but it had",
              "{nrow(parameters)} rows"),
        arg = "parameters", call = call)
    }
    nms <- rownames(parameters)
  } else {
    if (length(parameters) != n_pars) {
      cli::cli_abort(
        paste("Expected 'parameters' to have length {n_pars}, but it had",
              "length {length(parameters)}"),
        arg = "parameters")
    }
    nms <- names(parameters)
  }

  if (!is.null(nms) && !identical(nms, model$parameters)) {
    where <- if (is.matrix(parameters)) "rownames" else "names"
    cli::cli_abort(
      c("'parameters' has {where} but these disagree with 'model$parameters'",
        i = "If this is intentional, try again with 'unname(parameters)'"),
      arg = "parameters", call = call)
  }
}


##' @export
print.monty_model_properties <- function(x, ...) {
  cli::cli_h1("<monty_model_properties>")
  unset <- vlapply(x, is.null)
  is_set <- !unset
  if (any(is_set)) {
    cli::cli_bullets(
      set_names(sprintf("%s: {.code %s}",
                        names(x)[is_set], vcapply(x[is_set], as.character)),
                "*"))
  }
  if (any(unset)) {
    cli::cli_alert_info("Unset: {squote(names(x)[unset])}")
  }
}
