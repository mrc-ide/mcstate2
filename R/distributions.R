## This whole file is probably going to get an overhaul at some point
## because it's quite gross.  Mostly it is trying to tidy up some of
## the ways that we might draw from multivariate normal distributions.
## This is complicated by wanting to cache the results of the vcv
## decomposition where possible.

## Not in base R
rmvnorm <- function(x, vcv, rng) {
  make_rmvnorm(vcv)(x, rng)
}


## This is the form of the Cholesky factorisation of a matrix we use
## in the multivariate normal sampling.
chol_pivot <- function(x) {
  r <- chol(x, pivot = TRUE)
  r[, order(attr(r, "pivot", exact = TRUE))]
}


make_rmvnorm <- function(vcv, centred = FALSE) {
  n <- nrow(vcv)
  if (n == 1) {
    ## Special case for transformations in single-dimensional case; no
    ## need to work with matrix multiplication or decompositions here.
    sd <- sqrt(drop(vcv))
    if (centred) {
      function(rng) {
        rng$random_normal(1) * sd
      }
    } else {
      function(x, rng) {
        x + rng$random_normal(1) * sd
      }
    }
  } else if (is.matrix(vcv)) {
    r <- chol_pivot(vcv)
    if (centred) {
      function(rng) {
        drop(rng$random_normal(n) %*% r)
      }
    } else {
      function(x, rng) {
        x + drop(rng$random_normal(n) %*% r)
      }
    }
  } else {
    stopifnot(length(dim(vcv)) == 3)
    m <- dim(vcv)[[3]]
    r <- vapply(seq_len(m), function(i) chol_pivot(vcv[, , i]), vcv[, , 1])
    if (centred) {
      function(rng) {
        mu <- rng$random_normal(n)
        vapply(seq_len(m), function(i) mu[, i] %*% r[, , i], numeric(n))
      }
    } else {
      function(x, rng) {
        mu <- rng$random_normal(n) # + x perhaps?
        x + vapply(seq_len(m), function(i) mu[, i] %*% r[, , i], numeric(n))
      }
    }
  }
}


## log density multivariate normal
ldmvnorm <- function(x, vcv) {
  make_ldmvnorm(vcv)(x)
}


make_ldmvnorm <- function(vcv) {
  dec <- base::chol(vcv)
  constant <- -sum(log(diag(dec))) - 0.5 * nrow(vcv) * log(2 * pi)
  function(x) {
    constant - 0.5 * sum(backsolve(dec, x, transpose = TRUE)^2)
  }
}


## Not the world's most beautiful name.
deriv_ldmvnorm <- function(x, vcv) {
  make_deriv_ldmvnorm(vcv)(x)
}


make_deriv_ldmvnorm <- function(vcv) {
  vcv_inv <- solve(vcv)
  function(x) {
    -vcv_inv %*% x
  }
}
