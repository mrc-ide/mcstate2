##' Build a parameter packer, which can be used in models to translate
##' between an unstructured vector of numbers (the vector being
##' updated by an MCMC for example) to a structured list of named
##' values, which is easier to program against.  We refer to the
##' process of taking a named list of scalars, vectors and arrays and
##' converting into a single vector "packing" and the inverse
##' "unpacking".
##'
##' There are several places where it is most convenient to work in an
##' unstructured vector:
##'
##' * An MCMC is typically discussed as a the updating of some
##'   vector `x` to another `x'`
##' * An optimisation algorithm will try and find a set of values for
##'   a vector `x` that minimises (or maximises) some function `f(x)`
##' * An ode solver works with a vector `x(t)` (`x` at time `t`) and
##'   considers `x(t + h)` by computing the vector of derivatives
##'   `dx(t)/dt`
##'
##' In all these cases, the algorithm that needs the vector of numbers
##' knows nothing about what they represent.  Commonly, these will be
##' a packed vector of parameters.  So our vector `x` might actually
##' represent the parameters `a`, `b` and `c` in a vector as `[a, b,
##' c]` - this is a very common pattern, and you have probably
##' implemented this yourself.
##'
##' In more complex settings, we might want our vector `x` to collect
##' more structured quantities.  Suppose that you are fitting a model
##' with an age-structured or sex-structured parameter.  Rather than
##' having a series of scalars packed into your vector `x` you might
##' have a series of values destined to be treated as a vector:
##'
##' ```
##' | 1  2  3  4  5  6  7  |
##' | a  b  c  d1 d2 d3 d4 |
##' ```
##'
##' So here we might have a vector of length 7, where the first three
##' elements will represent be the scalar values `a`, `b` and `c` but
##' the next four will be a vector `d`.
##'
##' Unpacked, this might be written as:
##'
##' ```
##' list(a = 1, b = 2, c = 3, d = 4:7)
##' ```
##'
##' The machinery here is designed to make these transformations
##' simple and standardised within monty, and should be flexible
##' enough for many situations.  We will also use these from within
##' `dust2` and `odin2` for transformations in and out of vectors of
##' ODE state.
##'
##' # When to use `process`
##'
##' The `process` function is a get-out-of-jail function designed to
##' let you do arbitrary transformations when unpacking a vector.  In
##' general, this should not be the first choice to use because it is
##' less easy to reason about by other tooling (for example, as we
##' develop automatic differentiation support for use with the HMC
##' algorithm, a `process` function will be problematic because we
##' will need to make sure we can differentiate this process).
##' However, there are cases where it will be only way to achieve some
##' results.
##'
##' Imagine that you are packing a 2x2 covariance matrix into your
##' vector in order to use within an MCMC or optimisation algorithm.
##' Ultimately, our unpacked vector will need to hold four elements
##' (`b11`, `b12`, `b21`, `b22`), but there are only three distinct
##' values as the two off-diagonal elements will be the same (i.e.,
##' `b12 == b21`).  So we might write this passing in `b_raw = 3` to
##' `array`, so that our unpacked list holds `b_raw = c(b11, b12,
##' b22)`.  We would then write `process` as something like:
##'
##' ```
##' process <- function(x) {
##'   list(b = matrix(x$b_raw[c(1, 2, 2, 3)], 2, 2))
##' }
##' ```
##'
##' which creates the symmetric 2x2 matrix `b` from `b_raw`.
##'
##' # Unpacking matrices
##'
##' If you do not use `fixed` or `process` when defining your packer,
##' then you can use `$unpack()` with a matrix or higher-dimensional
##' output.  There are two ways that you might like to unpack this
##' sort of output.  Assume you have a matrix `m` with 3 rows and 2
##' columns; this means that we have two sets of parameters or state
##' (one per column) and 3 states within each; this is the format that
##' MCMC parameters will be in for example.
##'
##' The first would to be return a list where the `i`th element is the
##' result of unpacking the `i`th parameter/state vector.  You can do
##' this by running
##'
##' ```
##' apply(m, 2, p$unpack)
##' ```
##'
##' The second would be to return a named list with three elements
##' where the `ith` element is the unpacked version of the `i`th
##' state.  In this case you can pass the matrix directly in to the
##' unpacker:
##'
##' ```
##' p$unpack(m)
##' ```
##'
##' When you do this, the elements of `m` will acquire an additional
##' dimension; scalars become vectors (one per set), vectors become
##' matrices (one column per set) and so on.
##'
##' This approach generalises to higher dimensional input, though we
##' suspect you'll spend a bit of time head-scratching if you use it.
##'
##' # Packing lists into vectors and matrices
##'
##' The unpacking operation is very common - an MCMC proceeds,
##' produces an unstructured vector, and you unpack it into a list in
##' order to be able to easily work with it.  The reverse is much less
##' common, where we take a list and convert it into a vector (or
##' matrix, or multidimensional array).  Use of this direction
##' ("packing") may be more common where using packers to work with
##' the output of state-space models (e.g. in
##' [odin2](https://mrc-ide.github.io/odin2) or
##' [dust2](https://mrc-ide.github.io/dust2), which use this
##' machinery).
##'
##' The input to `pack()` will be the shape that `unpack()` returned;
##' a named list of numerical vectors, matrices and arrays.  The names
##' must correspond to the names in your packer (i.e., `scalar` and
##' the names of `array`).  Each element has dimensions
##'
##' ```
##' <...object, ...residual>
##' ```
##'
##' where `...object` is the dimensions of the data itself and
##' `...residual` is the dimensions of the hypothetical input to
##' `pack`.
##'
##' There is an unfortunate ambiguity in R's lack of true scalar types
##' that we cannot avoid.  It is hard to tell the difference packing a
##' vector vs packing an array where all dimensions are 1.  See the
##' examples, and please let us know if the behaviour needs changing.
##'
##' @title Build a parameter packer
##'
##' @param scalar Names of scalar parameters.  This is similar for
##'   listing elements in `array` with values of 1, though elements in
##'   `scalar` will be placed ahead of those listed in `array` within
##'   the final parameter vector, and elements in `array` will have
##'   generated names that include square brackets.
##'
##' @param array A list, where names correspond to the names of array
##'   parameters and values correspond to the lengths of parameters.
##'   Multiple dimensions are allowed (so if you provide an element
##'   with two entries these represent dimensions of a matrix).
##'   Zero-length integer vectors or `NULL` values are counted as
##'   scalars, which allows you to put scalars at positions other than
##'   the front of the packing vector. In future, you may be able to
##'   use *strings* as values for the lengths, in which case these
##'   will be looked for within `fixed`.
##'
##' @param fixed A named list of fixed parameters; these will be added
##'   into the final list directly.  These typically represent
##'   additional pieces of data that your model needs to run, but
##'   which you are not performing inference on.
##'
##' @param process An arbitrary R function that will be passed the
##'   final assembled parameter list; it may create any *additional*
##'   entries, which will be concatenated onto the original list.  If
##'   you use this you should take care not to return any values with
##'   the same names as entries listed in `scalar`, `array` or
##'   `fixed`, as this is an error (this is so that `pack()` is
##'   not broken).  We will likely play around with this process in
##'   future in order to get automatic differentiation to work.
##'
##' @return An object of class `monty_packer`, which has three
##'   elements:
##'
##' * `parameters`: a character vector of computed parameter names;
##'   these are the names that your statistical model will use.
##' * `unpack`: a function that can unpack an unstructured vector
##'   (say, from your statistical model parameters) into a structured
##'   list (say, for your generative model)
##' * `pack`: a function that can pack your structured list of
##'   parameters back into a numeric vector suitable for the
##'   statistical model.  This ignores values created by a
##'   `preprocess` function.
##' * `index`: a function which produces a named list where each
##'   element has the name of a value in `parameters` and each value
##'   has the indices within an unstructured vector where these values
##'   can be found.
##' * `subset`: an experimental interface which can be used to subset a
##'   packer to a packer for a subset of contents. Documentation will
##'   be provided once the interface settles.
##'
##' @export
##'
##' @examples
##' # Here's a really simple example
##' p <- monty_packer(c("a", "b", "c"))
##' p
##'
##' p$pack(list(a = 1, b = 2, c = 3))
##' p$unpack(1:3)
##'
##' # Sometimes we have a vector embedded in our parameters:
##' p <- monty_packer(c("a", "b"), list(v = 4))
##' p$pack(list(a = 1, b = 2, v = c(6, 7, 8, 9)))
##' p$unpack(c(1, 2, 6, 7, 8, 9))
##'
##' # Or a higher dimensional structure such as a matrix:
##' p <- monty_packer(c("a", "b"), list(m = c(2, 2)))
##' p$unpack(c(1, 2, 6, 7, 8, 9))
##'
##' # You can use a packer to set "fixed" parameters that do not vary
##' # with the underlying model being fit, but are required by your model.
##' # This is simpler than the "closure" approach used previously in our
##' # mcstate package and also easier to accommodate with differentiable
##' # models:
##' p <- monty_packer(
##'   c("a", "b"),
##'   fixed = list(d = data.frame(n = 1:3, m = runif(3))))
##' p$unpack(1:2)
##' p$pack(p$unpack(1:2))
##'
##' # The example from above, where we create a symmetric 2 x 2 matrix
##' # from a 3-element vector, alongside a scalar:
##' p <- monty_packer(
##'   scalar = "a",
##'   array = list(b_flat = 3),
##'   process = function(p) list(b = matrix(p$b_flat[c(1, 2, 2, 3)], 2, 2)))
##'
##' # Unpacking we see "b_flat" is still in the list, but "b" is our
##' # symmetric matrix:
##' p$unpack(1:4)
##'
##' # The processed elements are ignored on the return pack:
##' p$pack(list(a = 1, b_flat = 2:4, b = matrix(c(2, 3, 3, 4), 2, 2)))
##' p$pack(list(a = 1, b_flat = 2:4))
##'
##' # R lacks scalars, which means that some packers will unpack
##' # different inputs to the same outputs:
##' p <- monty_packer(c("a", "b"))
##' p$unpack(1:2)
##' p$unpack(cbind(1:2))
##'
##' # This means that we can't reliably pack these inputs in a way
##' # that guarantees round-tripping is possible.  We have chosen to
##' # prioritise the case where a *single vector* is round-trippable:
##' p$pack(list(a = 1, b = 2))
##'
##' # This ambiguity goes away if unpacking matices with more than one
##' # column:
##' p$unpack(matrix(1:6, 2, 3))
monty_packer <- function(scalar = NULL, array = NULL, fixed = NULL,
                         process = NULL) {
  call <- environment()

  parameters <- character(0)
  idx <- list()
  shape <- list()

  if (!is.null(scalar)) {
    if (!is.character(scalar)) {
      cli::cli_abort(
        "Expected a character vector for 'scalar'",
        arg = "scalar")
    }
    dups <- duplicate_values(scalar)
    if (length(dups) > 0) {
      cli::cli_abort(
        c("Elements of 'scalar' must be unique",
          i = "Found {length(dups)} duplicate{?s}: {collapseq(dups)}"),
        arg = "scalar")
    }
    idx[scalar] <- set_names(as.list(seq_along(scalar)), scalar)
    parameters <- c(parameters, scalar)
  }

  len <- length(scalar) # start after scalars
  if (!is.null(array)) {
    assert_named(array, unique = TRUE, call = call)
    for (nm in names(array)) {
      tmp <- prepare_pack_array(nm, array[[nm]], call)
      parameters <- c(parameters, tmp$names)
      shape[[nm]] <- tmp$shape
      idx[[nm]] <- seq_len(tmp$n) + len
      len <- len + tmp$n
    }
  }

  if (!is.null(fixed)) {
    assert_named(fixed, unique = TRUE, call = call)
    assert_list(fixed, call = call)
  }

  err <- duplicate_values(c(scalar, names(array), names(fixed)))
  if (length(err)) {
    detail <- sprintf("'%s' appears in more than one place", err)
    cli::cli_abort(
      c("Names must be distinct between 'scalar', 'array' and 'fixed'",
        set_names(detail, "x")))
  }

  if (length(idx) == 0) {
    cli::cli_abort(
      c("Trying to generate an empty packer",
        i = paste("You have not provided any entries in 'scalar' or 'array',",
                  "which implies generating from a zero-length parameter",
                  "vector")))
  }

  if (!is.null(process)) {
    if (!is.function(process)) {
      cli::cli_abort(
        "Expected a function for 'process'",
        arg = "process")
    }
  }

  unpack <- function(x) {
    if (is.null(dim(x))) {
      unpack_vector(x, parameters, len, idx, shape, fixed, process)
    } else {
      unpack_array(x, parameters, len, idx, shape, fixed, process)
    }
  }

  pack <- function(p) {
    ## We might want a more explicit error here?
    assert_named(p, unique = TRUE)
    ## TODO: drop any processed and fixed things here, which means
    ## we're more concerned about finding things we can reshape than
    ## anything else.
    shp <- pack_check_dimensions(p, scalar, shape, names(fixed), process)
    ret <- matrix(NA_real_, len, prod(shp))
    for (nm in c(scalar, names(shape))) {
      ret[idx[[nm]], ] <- p[[nm]]
    }

    drop <- length(shp) == 0 ||
      (length(shp) == 1 && (length(shape) == 0 || all(lengths(shape) == 0)))
    if (drop) {
      dim(ret) <- NULL
    } else if (length(shp) > 1) {
      dim(ret) <- c(len, shp)
    }
    ret
  }

  subset <- function(keep) {
    ## TODO: later we will allow passing integer indexes here.  This
    ## is complicated because we should probably retain structure for
    ## any compartments that are entirely captured (contiguously and
    ## in order) and convert everything else into scalars.  Or perhaps
    ## if we take a slice out of a matrix we keep it as an array.
    ## Lots of decisions to make, so do it later.
    if (is.character(keep)) {
      if (anyDuplicated(keep)) {
        dups <- unique(keep[duplicated(keep)])
        cli::cli_abort("Duplicated name{?s} in 'keep': {squote(dups)}")
      }
      i <- match(keep, names(idx))
      if (anyNA(i)) {
        cli::cli_abort("Unknown name{?s} in 'keep': {squote(keep[is.na(i)])}")
      }
      index <- unlist(idx[i], FALSE, FALSE)
      ## Convert all scalars into array spec for now; this allows
      ## reordering.
      shape2 <- c(set_names(rep(list(integer()), length(scalar)), scalar),
                  shape)
      scalar_keep <- NULL
      array_keep <- shape2[keep]
    } else {
      cli::cli_abort(
        "Invalid input for 'keep'; this must currently be a character vector")
    }
    list(index = index,
         packer = monty_packer(scalar_keep, array_keep))
  }

  ret <- list(parameters = parameters,
              unpack = unpack,
              pack = pack,
              index = function() idx,
              subset = subset)
  class(ret) <- "monty_packer"
  ret
}


## Helper function to create array bookkeeping for the
## unpacking/packing process.
prepare_pack_array <- function(name, shape, call = NULL) {
  if (is.null(shape)) {
    shape <- integer()
  }
  if (!rlang::is_integerish(shape)) {
    cli::cli_abort(
      paste("Elements of 'array' must be integer-like vectors, but",
            "'{name}' is not"),
      arg = "array", call = call)
  }
  shape <- as.integer(shape)
  if (length(shape) == 0) {
    return(list(names = name, shape = integer(0), n = 1L))
  }
  if (any(shape <= 0)) {
    cli::cli_abort(
      paste("All dimensions in 'array' must be at least 1, but",
            "'{name}' violates this"),
      arg = "array", call = call)
  }
  if (length(shape) == 1) {
    index <- as.character(seq_len(shape))
  } else {
    index <- apply(array_indices(shape), 1, paste, collapse = ",")
  }

  ## Whole bunch of bookkeeping used above;

  ## * names contains generated parameter names with index accessors -
  ##   these will render nicely though they're gross in their own way.
  ## * shape contains generated dimensions of a resulting array
  ## * n is the total number of parameters after this expansion
  list(names = sprintf("%s[%s]", name, index),
       shape = shape,
       n = length(index))
}


array_indices <- function(shape) {
  unname(as.matrix(do.call(expand.grid, lapply(shape, seq_len))))
}


##' @export
print.monty_packer <- function(x, ...) {
  cli::cli_h1("<monty_packer>")
  cli::cli_alert_info(
    "Packing {length(x$parameters)} parameter{?s}: {squote(x$parameters)}")
  cli::cli_alert_info(
    "Use '$pack()' to convert from a list to a vector")
  cli::cli_alert_info(
    "Use '$unpack()' to convert from a vector to a list")
  cli::cli_alert_info("See {.help monty_packer} for more information")
  invisible(x)
}


unpack_vector <- function(x, parameters, len, idx, shape, fixed, process) {
  call <- parent.frame()
  if (!is.null(names(x))) {
    if (!identical(names(x), parameters)) {
      ## Here, we could do better I think with this message; we
      ## might pass thropuigh empty names, and produce some summary
      ## of different names.  Something for later though.
      cli::cli_abort("Incorrect names in input")
    }
    names(x) <- NULL
  }
  if (length(x) != len) {
    cli::cli_abort(
      "Incorrect length input; expected {len} but given {length(x)}",
      call = call)
  }
  res <- lapply(idx, function(i) x[i])
  for (nm in names(shape)) {
    dim2(res[[nm]]) <- shape[[nm]]
  }
  if (!is.null(fixed)) {
    res <- c(res, fixed)
  }
  if (!is.null(process)) {
    extra <- process(res)
    err <- intersect(names(extra), names(res))
    if (length(err) > 0) {
      cli::cli_abort(
        c("'process()' is trying to overwrite entries in parameters",
          i = paste("The 'process()' function should only create elements",
                    "that are not already present in 'scalar', 'array'",
                    "or 'fixed', as this lets us reverse the transformation",
                    "process"),
          x = "{?Entry/Entries} already present: {squote(err)}"))
    }
    ## TODO: check names?
    ##
    ## TODO: this fails the multi-region use - but I think that we
    ## might want a different interface there anyway as we'll
    ## struggle to hold all the options here.
    res <- c(res, extra)
  }
  res
}


unpack_array <- function(x, parameters, len, idx, shape, fixed, process) {
  call <- parent.frame()
  dn <- dimnames(x)
  if (!is.null(dn) && !is.null(dn[[1]]) && !identical(dn[[1]], parameters)) {
    ## See comment above about reporting on this better
    cli::cli_abort("Incorrect rownames in input")
  }

  if (nrow(x) != len) {
    cli::cli_abort(
      paste("Incorrect length of first dimension of input;",
            "expected {len} but given {nrow(x)}"),
      call = call)
  }

  shape_x <- dim(x)[-1]
  dim(x) <- c(dim(x)[[1]], prod(shape_x))

  res <- lapply(idx, function(i) x[i, ])
  for (nm in names(shape)) {
    dim2(res[[nm]]) <- c(shape[[nm]], shape_x)
  }
  if (!is.null(fixed)) {
    cli::cli_abort("Can't unpack a matrix where the unpacker uses 'fixed'",
                   call = call)
  }
  if (!is.null(process)) {
    cli::cli_abort("Can't unpack a matrix where the unpacker uses 'process'",
                   call = call)
  }

  res
}


## Now, we do some annoying calculations to make sure that what we've
## been given has the correct size, etc.
pack_check_dimensions <- function(p, scalar, shape, fixed, process,
                                  call = parent.frame()) {
  if (length(scalar) > 0) {
    shape <- c(set_names(rep(list(integer()), length(scalar)), scalar),
               shape)
  }

  msg <- setdiff(names(shape), names(p))
  extra <- if (is.null(process)) setdiff(names(p), c(names(shape), fixed))

  err <- setdiff(names(shape), names(p))
  if (length(err) > 0) {
    cli::cli_abort("Missing element{?s} from input to pack: {squote(err)}",
                   call = call)
  }
  if (length(extra) > 0) {
    cli::cli_abort(
      "Unexpected element{?s} present in input to pack: {squote(extra)}",
      call = call)
  }

  ## It's very tedious to check if we can map a set of inputs to a
  ## single pack.  What we expect is if we have an element 'x' with
  ## dimensions <a, b, c, ...> and we expect that the shape 's' of
  ## this is <a, b> we need to validate that the first dimensions of
  ## 'x' are 's' and then record the remainder - that's the shape of
  ## the rest of the eventual object (so if there is no dimension left
  ## then it was originally a vector, if there's one element left our
  ## original input was a matrix and so on).  Then we check that this
  ## remainder is consistent across all elements in the list.
  ##
  ## There are a couple of additional wrinkles due to scalar variables
  ## (these are ones where 's' has length 0).  First, we can't drop
  ## things from 'd' with d[-i] because that will drop everything.
  ## Second, when working out the residual dimensions of the packed
  ## data we might have a situation where there are some scalars for
  ## which we think the residual dimension is 1 and some non-scalars
  ## where we think there is no residual dimension.  This is actually
  ## the same situation.  This is explored a bit in the tests and the
  ## examples above.
  res <- lapply(names(shape), function(nm) {
    d <- dim2(p[[nm]])
    s <- shape[[nm]]
    i <- seq_along(s)
    if (length(d) >= length(s) && identical(d[i], s)) {
      list(success = TRUE, residual = if (length(i) == 0) d else d[-i])
    } else {
      list(success = FALSE, shape = d[i])
    }
  })

  err <- !vlapply(res, "[[", "success")
  if (any(err)) {
    err_nms <- names(shape)[err]
    detail <- sprintf(
      "%s: expected <%s>, given <%s>",
      err_nms,
      vcapply(shape[err], paste, collapse = ", "),
      vcapply(res[err], function(x) paste(x$shape, collapse = ", ")))
    cli::cli_abort(
      c("Incompatible dimensions in input for {squote(names(shape)[err])}",
        set_names(detail, "x")),
      call = call)
  }

  residual <- lapply(res, "[[", "residual")
  is_scalar <- lengths(residual) == 0
  is_vector_output <- any(is_scalar) &&
    all(lengths(residual[!is_scalar]) == 1) &&
    all(vnapply(residual[!is_scalar], identity) == 1)
  if (is_vector_output) {
    return(integer())
  }

  ret <- residual[[1]]
  ok <- vlapply(residual, identical, ret)
  if (!all(ok)) {
    residual[is_scalar] <- 1L
    hash <- vcapply(residual, rlang::hash)
    detail <- vcapply(unique(hash), function(h) {
      sprintf("%s: <...%s>",
              paste(squote(names(shape)[hash == h]), collapse = ", "),
              paste(residual[hash == h][[1]], collapse = ", "))
    })
    cli::cli_abort(
      c("Inconsistent residual dimension in inputs",
        set_names(detail, "x")),
      call = call)
  }

  ret
}
