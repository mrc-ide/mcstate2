progress_bar <- function(n_chains, n_steps, progress, show_overall,
                         call = NULL) {
  progress <- show_progress_bar(progress, call)
  if (progress) {
    progress_bar_detail(n_chains, n_steps, show_overall)
  } else {
    progress_bar_null()
  }
}


show_progress_bar <- function(progress, call = NULL) {
  if (is.null(progress)) {
    progress <- getOption("mcstate2.progress", TRUE)
  }
  ## Error here is not great if we get this from an option; but need
  ## to disable quoting there in some cases.
  assert_scalar_logical(progress, call = call)
}


## This will be the case to use where we can report back in detail
## about what is running (for each chain we can report about the
## progress within that chain).  The other likely option will be a
## very coarse detail which is the state of each chain only; that
## might be the best we can do for some of the other parallel
## backends, but we'll see how packages that implement progress bars
## there cope.
progress_bar_detail <- function(n_chains, n_steps, show_overall) {
  e <- new.env()
  e$n <- rep(0, n_chains)
  overall <- progress_overall(n_chains, n_steps, show_overall)
  fmt <- paste("Sampling {overall(e$n)} {cli::pb_bar} |",
               "{cli::pb_percent} ETA: {cli::pb_eta}")
  id <- cli::cli_progress_bar(
    total = n_chains * n_steps,
    format = fmt,
    .auto_close = FALSE)

  function(chain_index) {
    function(at) {
      e$n[chain_index] <- at
      cli::cli_progress_update(id = id, set = sum(e$n))
    }
  }
}


## Dummy version that can be used where no progress bar is wanted.
progress_bar_null <- function(...) {
  function(chain_index) {
    function(at) {
    }
  }
}


progress_overall <- function(n_chains, n_steps, show_overall) {
  if (n_chains == 1 || !show_overall) {
    return(function(n) "")
  }
  sym <- unlist(cli::symbol[paste0("lower_block_", 1:8)], use.names = FALSE)
  at <- seq(0, n_steps, length.out = length(sym))
  ## It would be much better to look up the theme here really, and use
  ## colours from that, but I'm not totally sure what the mechanism is
  ## for that.  Doing it the theme way will result in darl/light
  ## appropriate colours being picked automatically which would be
  ## nice.
  col_running <- cli::make_ansi_style("orange")
  col_finished <- cli::make_ansi_style("green")
  function(n) {
    ret <- sym[findInterval(n, at)]
    i_finished <- n == n_steps
    i_running <- !i_finished & n > 0
    ret[i_finished] <- col_finished(ret[i_finished])
    ret[i_running] <- col_running(ret[i_running])
    paste0(c("[", ret, "]"), collapse = "")
  }
}
