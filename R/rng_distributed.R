##' Create a set of initial random number seeds suitable for using
##' within a distributed context (over multiple processes or nodes) at
##' a level higher than a single group of synchronised threads.
##'
##' See `vignette("rng_distributed")` for a proper introduction to
##' these functions.
##'
##' @title Create a set of distributed seeds
##'
##' @param seed Initial seed to use. As for [monty::monty_rng], this can
##'   be `NULL` (create a seed using R's generators), an integer or a
##'   raw vector of appropriate length.
##'
##' @param n_streams The number of streams to create per node.
##'
##' @param n_nodes The number of separate seeds to create. Each will
##'   be separated by a "long jump" for your generator.
##'
##' @param algorithm The name of an algorithm to use.
##'
##' @return A list of either raw vectors (for
##'   `monty_rng_distributed_state`) or of [monty::monty_rng_pointer]
##'   objects (for `monty_rng_distributed_pointer`)
##'
##' @export
##' @rdname monty_rng_distributed
##' @examples
##' monty::monty_rng_distributed_state(n_nodes = 2)
##' monty::monty_rng_distributed_pointer(n_nodes = 2)
monty_rng_distributed_state <- function(seed = NULL,
                                        n_streams = 1L,
                                        n_nodes = 1L,
                                        algorithm = "xoshiro256plus") {
  p <- monty_rng_pointer$new(seed, n_streams, algorithm = algorithm)

  ret <- vector("list", n_nodes)
  for (i in seq_len(n_nodes)) {
    s <- p$state()
    ret[[i]] <- s
    if (i < n_nodes) {
      p <- monty_rng_pointer$new(s, n_streams, 1L, algorithm = algorithm)
    }
  }

  ret
}


##' @export
##' @rdname monty_rng_distributed
monty_rng_distributed_pointer <- function(seed = NULL,
                                          n_streams = 1L,
                                          n_nodes = 1L,
                                          algorithm = "xoshiro256plus") {
  state <- monty_rng_distributed_state(seed, n_streams, n_nodes, algorithm)
  lapply(state, monty_rng_pointer$new,
         n_streams = n_streams, algorithm = algorithm)
}
