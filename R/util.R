`%||%` <- function(x, y) { # nolint
  if (is.null(x)) y else x
}


check_vcv <- function(vcv, name = deparse(substitute(vcv)), call = NULL) {
  if (!is.matrix(vcv)) {
    cli::cli_abort("Expected '{name}' to be a matrix",
                   arg = name, call = call)
  }
  if (!isSymmetric(vcv)) {
    cli::cli_abort("Expected '{name}' to be symmetric",
                   arg = name, call = call)
  }
  if (!is_positive_definite(vcv)) {
    cli::cli_abort("Expected '{name}' to be positive definite",
                   arg = name, call = call)
  }
}


is_positive_definite <- function(x, tol = sqrt(.Machine$double.eps)) {
  ev <- eigen(x, symmetric = TRUE)
  all(ev$values >= -tol * abs(ev$values[1]))
}


mcstate_file <- function(path) {
  system.file(path, package = "mcstate2", mustWork = TRUE)
}


## Small helper to get R's RNG state; we can use this to check that
## the state was not updated unexpectedly.
get_r_rng_state <- function() {
  if (exists(".Random.seed", globalenv(), inherits = FALSE)) {
    get(".Random.seed", globalenv(), inherits = FALSE)
  } else {
    NULL
  }
}


vlapply <- function(...) {
  vapply(..., FUN.VALUE = TRUE)
}


vnapply <- function(...) {
  vapply(..., FUN.VALUE = 1)
}


vcapply <- function(...) {
  vapply(..., FUN.VALUE = "")
}


rbind_list <- function(x) {
  stopifnot(all(vlapply(x, is.matrix)))
  if (length(x) == 1) {
    return(x[[1]])
  }
  nc <- unique(vnapply(x, ncol))
  stopifnot(length(nc) == 1)
  do.call("rbind", x)
}


squote <- function(x) {
  sprintf("'%s'", x)
}


dim2 <- function(x) {
  dim(x) %||% length(x)
}


"dim2<-" <- function(x, value) {
  if (length(value) != 1) {
    dim(x) <- value
  }
  x
}


dimnames2 <- function(x) {
  if (!is.null(dim(x))) {
    dimnames(x)
  } else {
    nms <- names(x)
    if (is.null(nms)) NULL else list(nms)
  }
}


set_names <- function(x, nms) {
  if (length(nms) == 1 && length(x) != 1) {
    nms <- rep_len(nms, length(x))
  }
  names(x) <- nms
  x
}


all_same <- function(x) {
  if (length(x) < 2) {
    return(TRUE)
  }
  all(vlapply(x[-1], identical, x[[1]]))
}


as_function <- function(args, body, envir) {
  if (!rlang::is_call(body, "{")) {
    body <- rlang::call2("{", !!!body)
  }
  as.function(c(args, body), envir = envir)
}


near_match <- function(x, possibilities, threshold = 2, max_matches = 5) {
  i <- tolower(x) == tolower(possibilities)
  if (any(i)) {
    possibilities[i]
  } else {
    d <- set_names(drop(utils::adist(x, possibilities, ignore.case = TRUE)),
                   possibilities)
    utils::head(names(sort(d[d <= threshold])), max_matches)
  }
}


duplicate_values <- function(x) {
  if (!anyDuplicated(x)) {
    return(x[integer(0)])
  }
  unique(x[duplicated(x)])
}


collapseq <- function(x) {
  paste(squote(x), collapse = ", ")
}
