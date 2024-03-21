check_parameter_groups <- function(x, n_pars, name = deparse(substitute(x)),
                                   call = NULL) {
  if (!rlang::is_integerish(x)) {
    cli::cli_abort("Expected '{name}' to be integer-like", call = call)
  }
  if (length(x) != n_pars) {
    cli::cli_abort(
      "Expected '{name}' to have length {n_pars}, but was {length(x)}",
      call = call)
  }
  n_groups <- max(n_pars)
  msg <- setdiff(seq_len(n_pars), n_groups)
  if (length(msg) > 0) {
    # TODO: Better error here that explains the situation better and
    # offers a hint as to what to do.
    cli::cli_abort("Missing groups from '{name}'", call = call)
  }
  if (min(x) < 0) {
    cli::cli_abort("Invalid negative groups in '{name}'", call = call)
  }
}
