dsl_generate <- function(dat) {
  density <- dsl_generate_density(dat)
  mcstate_model(
    list(parameters = dat$parameters,
         density = density))
}


dsl_generate_density <- function(dat) {
  ## Here we'll generate a series of self-contained statements, and
  ## evaluate these in a generic environment.  This disallows use of
  ## fancy functions that the user provides, but I think that the dsl
  ## does not really support that anyway, so that feels reasonable.
  ##
  ## We're using .x and .density as special input/output vectors, at
  ## least for now; I don't really know how sustainable that is but
  ## stops us needing to use a second layer of indirection (e.g., all
  ## assignments go into some additional list/env and we need to
  ## manipulate all the calling environments accordingly).  When we
  ## start generating C code later this will change anyway as there
  ## we'll need to do lookups anyway.
  np <- length(dat$parameters)
  body <- c(dsl_generate_density_unpack(dat$parameters),
            quote(.density <- numeric()),
            lapply(dat$exprs, dsl_generate_density_expr),
            quote(sum(.density)))
  as_function(alist(.x = ), body, topenv())
}


dsl_generate_density_unpack <- function(parameters) {
  lapply(seq_along(parameters), function(i) {
    call("<-", as.name(parameters[[i]]), bquote(.x[.(i)]))
  })
}


dsl_generate_density_expr <- function(expr) {
  switch(expr$type,
         assignment = dsl_generate_density_assignment(expr),
         stochastic = dsl_generate_density_stochastic(expr),
         cli::cli_abort(paste(
           "Unimplemented expression type '{expr$type}';",
           "this is an mcstate2 bug")))
}


dsl_generate_density_assignment <- function(expr) {
  expr$expr
}


dsl_generate_density_stochastic <- function(expr) {
  lhs <- bquote(.density[[.(expr$name)]])
  rhs <- rlang::call2(expr$distribution$density,
                      as.name(expr$name), !!!expr$args)
  rlang::call2("<-", lhs, rhs)
}
