progress_bar <- function(n_chains, n_steps, progress, show_overall,
                         single_chain = FALSE, call = parent.frame()) {
  progress <- show_progress_bar(progress, call)
  switch(
    progress,
    none = progress_bar_none(),
    simple = progress_bar_simple(n_steps),
    fancy = progress_bar_fancy(n_chains, n_steps, show_overall, single_chain))
}


show_progress_bar <- function(progress, call = NULL) {
  if (is.null(progress)) {
    progress <- getOption("monty.progress", TRUE)
  }
  ## Errors here are not great if we get this from an option, probably
  ## needs its own error path.
  if (is.logical(progress)) {
    assert_scalar_logical(progress, call = call)
    if (progress) "fancy" else "none"
  } else {
    match_value(progress, c("fancy", "simple", "none"))
  }
}


## A "progress bar" that we can pass through that will produce text
## that we can easily scan for in the logs of a process from callr.
##
## Arguments here allow some tuning to get this "appropriately
## responsive"; default is print every second, or at least 20 "ticks"
## within a run.
##
## These are currently not tuneable from user-facing code.
progress_bar_simple <- function(n_steps, every_s = 1, min_updates = 20) {
  function(chain_index) {
    env <- new.env(parent = emptyenv())
    env$t_next <- Sys.time()
    freq <- ceiling(n_steps / min_updates)
    function(at) {
      now <- Sys.time()
      show_progress <- at == n_steps || at %% freq == 0 || now > env$t_next
      if (show_progress) {
        env$t_next <- now + every_s
        message(sprintf("MONTY-PROGRESS: chain: %s, step: %s",
                        chain_index, at))
      }
    }
  }
}


progress_bar_fancy <- function(n_chains, n_steps, show_overall,
                               single_chain = FALSE) {
  e <- new.env()
  e$n <- rep(0, n_chains)
  overall <- progress_overall(n_chains, n_steps, show_overall, single_chain)
  fmt <- paste("Sampling {overall(e$n)} {cli::pb_bar} |",
               "{cli::pb_percent} ETA: {cli::pb_eta}")
  n_steps_total <- if (single_chain) n_steps else n_chains * n_steps
  id <- cli::cli_progress_bar(
    total = n_steps_total,
    format = fmt,
    .auto_close = FALSE)

  function(chain_index) {
    function(at) {
      ## Avoid writing into a closed progress bar, it will cause an
      ## error.  We do this by checking to see if progress has changed
      ## from last time we tried updating.
      changed <- any(e$n[chain_index] != at, na.rm = TRUE)
      if (changed) {
        e$n[chain_index] <- at
        cli::cli_progress_update(id = id, set = sum(e$n))
      }
    }
  }
}


parse_progress_bar_simple <- function(txt) {
  re <- "^MONTY-PROGRESS: chain: ([0-9]+), step: ([0-9]+)$"
  i <- grep(re, txt)
  if (length(i) == 0) {
    NULL
  } else {
    x <- txt[[last(i)]]
    list(chain_id = as.numeric(sub(re, "\\1", x)),
         step = as.numeric(sub(re, "\\2", x)))
  }
}


## Dummy version that can be used where no progress bar is wanted.
progress_bar_none <- function(...) {
  function(chain_index) {
    function(at) {
    }
  }
}


progress_overall <- function(n_chains, n_steps, show_overall, single_chain) {
  if (n_chains == 1 || !show_overall || single_chain) {
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
