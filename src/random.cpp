#include <cstring>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <cpp11/doubles.hpp>
#include <cpp11/external_pointer.hpp>
#include <cpp11/integers.hpp>
#include <cpp11/raws.hpp>

#include <mcstate/r/random.hpp>
#include <mcstate/random/random.hpp>
#include <mcstate/utils.hpp>

using default_rng64 = mcstate::random::prng<mcstate::random::generator<double>>;
using default_rng32 = mcstate::random::prng<mcstate::random::generator<float>>;

template <typename T>
SEXP mcstate_rng_alloc(cpp11::sexp r_seed, int n_streams, bool deterministic) {
  auto seed = mcstate::random::r::as_rng_seed<typename T::rng_state>(r_seed);
  T *rng = new T(n_streams, seed, deterministic);
  return cpp11::external_pointer<T>(rng);
}

template <typename T>
void mcstate_rng_jump(SEXP ptr) {
  T *rng = cpp11::as_cpp<cpp11::external_pointer<T>>(ptr).get();
  rng->jump();
}

template <typename T>
void mcstate_rng_long_jump(SEXP ptr) {
  T *rng = cpp11::as_cpp<cpp11::external_pointer<T>>(ptr).get();
  rng->long_jump();
}

// Little helper for returning x as a vector (m == 1) or matrix (n *
// m) by setting the dimension attribute.
cpp11::sexp sexp_matrix(cpp11::sexp x, int n, int m) {
  if (m > 1) {
    x.attr("dim") = cpp11::writable::integers{n, m};
  }
  return x;
}

template <typename real_type, typename T>
cpp11::sexp mcstate_rng_random_real(SEXP ptr, int n, int n_threads) {
  T *rng = cpp11::as_cpp<cpp11::external_pointer<T>>(ptr).get();
  const int n_streams = rng->size();

  cpp11::writable::doubles ret = cpp11::writable::doubles(n * n_streams);
  double * y = REAL(ret);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
  for (int i = 0; i < n_streams; ++i) {
    auto &state = rng->state(i);
    auto y_i = y + n * i;
    for (size_t j = 0; j < (size_t)n; ++j) {
      y_i[j] = mcstate::random::random_real<real_type>(state);
    }
  }

  return sexp_matrix(ret, n, n_streams);
}

template <typename real_type, mcstate::random::algorithm::normal A, typename T>
cpp11::sexp mcstate_rng_random_normal(SEXP ptr, int n, int n_threads) {
  T *rng = cpp11::as_cpp<cpp11::external_pointer<T>>(ptr).get();
  const int n_streams = rng->size();

  cpp11::writable::doubles ret = cpp11::writable::doubles(n * n_streams);
  double * y = REAL(ret);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
  for (int i = 0; i < n_streams; ++i) {
    auto &state = rng->state(i);
    auto y_i = y + n * i;
    for (size_t j = 0; j < (size_t)n; ++j) {
      y_i[j] = mcstate::random::random_normal<real_type, A>(state);
    }
  }

  return sexp_matrix(ret, n, n_streams);
}

struct input_vary {
  size_t len;
  size_t offset;
  bool draw;
  bool generator;
};

// See notes in R/rng.R or ?rng
input_vary check_input_type(cpp11::doubles x, int n, int m, const char *name) {
  input_vary ret {1, 1, false, false};
  if (Rf_isMatrix(x)) {
    if (Rf_ncols(x) != m) {
      cpp11::stop("If '%s' is a matrix, it must have %d columns", name, m);
    }
    ret.generator = true;
    if (Rf_nrows(x) == n) {
      ret.draw = true;
    } else if (Rf_nrows(x) != 1) {
      cpp11::stop("If '%s' is a matrix, it must have 1 or %d rows", name, n);
    }
  } else {
    if (x.size() == n) {
      ret.draw = true;
    } else if (x.size() != 1) {
      cpp11::stop("If '%s' is a vector, it must have 1 or %d elements",
                  name, n);
    }
  }

  if (ret.draw) {
    ret.offset = n;
  }

  return ret;
}

// See notes in R/rng.R or ?rng
input_vary check_input_type2(cpp11::doubles x, int n, int m, const char *name) {
  input_vary ret {1, 1, false, false};
  cpp11::sexp r_dim = x.attr("dim");
  if (r_dim == R_NilValue) {
    ret.len = x.size();
  } else if (LENGTH(r_dim) == 2) { // matrix
    auto dim = cpp11::as_cpp<cpp11::integers>(r_dim);
    ret.len = dim[0];
    if (dim[1] == n) {
      ret.draw = true;
    } else {
      // TODO: must be n, not 1 surely?
      cpp11::stop("If '%s' is a matrix, it must have %d columns", name, n);
    }
  } else if (LENGTH(r_dim) == 3) {
    auto dim = cpp11::as_cpp<cpp11::integers>(r_dim);
    ret.len = dim[0];
    if (dim[1] == n) {
      ret.draw = true;
    } else if (dim[1] != 1) {
      cpp11::stop("If '%s' is a 3d array, it must have 1 or %d columns",
                  name, n);
    }
    if (dim[2] != m) {
      cpp11::stop("If '%s' is a 3d array, it must have %d layers", name, m);
    }
    ret.generator = true;
  } else {
    cpp11::stop("'%s' must be a vector, matrix or 3d array", name);
  }

  if (ret.len < 2) {
    cpp11::stop("Input parameters imply length of '%s' of only %d (< 2)",
                name, ret.len);
  }

  if (ret.draw) {
    ret.offset = n * ret.len;
  } else {
    ret.offset = ret.len;
  }

  return ret;
}

// Below here is very repetitive, and could probably be deduplicated
// with some clever template magic. Most of the faff is because we
// want to support 4 modes of taking 1 or 2 parameters (each varying
// or not over draws and generators)
template <typename real_type, typename T>
cpp11::sexp mcstate_rng_uniform(SEXP ptr, int n,
                             cpp11::doubles r_min,
                             cpp11::doubles r_max,
                             int n_threads) {
  T *rng = cpp11::as_cpp<cpp11::external_pointer<T>>(ptr).get();
  const int n_streams = rng->size();
  cpp11::writable::doubles ret = cpp11::writable::doubles(n * n_streams);
  double * y = REAL(ret);

  const double * min = REAL(r_min);
  const double * max = REAL(r_max);
  auto min_vary = check_input_type(r_min, n, n_streams, "min");
  auto max_vary = check_input_type(r_max, n, n_streams, "max");

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
  for (int i = 0; i < n_streams; ++i) {
    auto &state = rng->state(i);
    auto y_i = y + n * i;
    auto min_i = min_vary.generator ? min + min_vary.offset * i : min;
    auto max_i = max_vary.generator ? max + max_vary.offset * i : max;
    for (size_t j = 0; j < (size_t)n; ++j) {
      auto min_ij = min_vary.draw ? min_i[j] : min_i[0];
      auto max_ij = max_vary.draw ? max_i[j] : max_i[0];
      y_i[j] = mcstate::random::uniform<real_type>(state, min_ij, max_ij);
    }
  }

  return sexp_matrix(ret, n, n_streams);
}

template <typename real_type, typename T>
cpp11::sexp mcstate_rng_exponential_rate(SEXP ptr, int n, cpp11::doubles r_rate,
                                         int n_threads) {
  T *rng = cpp11::as_cpp<cpp11::external_pointer<T>>(ptr).get();
  const int n_streams = rng->size();
  cpp11::writable::doubles ret = cpp11::writable::doubles(n * n_streams);
  double * y = REAL(ret);

  const double * rate = REAL(r_rate);
  auto rate_vary = check_input_type(r_rate, n, n_streams, "rate");

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
  for (int i = 0; i < n_streams; ++i) {
    auto &state = rng->state(i);
    auto y_i = y + n * i;
    auto rate_i = rate_vary.generator ? rate + rate_vary.offset * i : rate;
    for (size_t j = 0; j < (size_t)n; ++j) {
      auto rate_ij = rate_vary.draw ? rate_i[j] : rate_i[0];
      y_i[j] = mcstate::random::exponential_rate<real_type>(state, rate_ij);
    }
  }

  return sexp_matrix(ret, n, n_streams);
}


template <typename real_type, typename T>
cpp11::sexp mcstate_rng_exponential_mean(SEXP ptr, int n, cpp11::doubles r_mean,
                                         int n_threads) {
  T *rng = cpp11::as_cpp<cpp11::external_pointer<T>>(ptr).get();
  const int n_streams = rng->size();
  cpp11::writable::doubles ret = cpp11::writable::doubles(n * n_streams);
  double * y = REAL(ret);

  const double * mean = REAL(r_mean);
  auto mean_vary = check_input_type(r_mean, n, n_streams, "mean");

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
  for (int i = 0; i < n_streams; ++i) {
    auto &state = rng->state(i);
    auto y_i = y + n * i;
    auto mean_i = mean_vary.generator ? mean + mean_vary.offset * i : mean;
    for (size_t j = 0; j < (size_t)n; ++j) {
      auto mean_ij = mean_vary.draw ? mean_i[j] : mean_i[0];
      y_i[j] = mcstate::random::exponential_mean<real_type>(state, mean_ij);
    }
  }

  return sexp_matrix(ret, n, n_streams);
}

template <typename real_type, mcstate::random::algorithm::normal A, typename T>
cpp11::sexp mcstate_rng_normal(SEXP ptr, int n,
                            cpp11::doubles r_mean, cpp11::doubles r_sd,
                            int n_threads) {
  T *rng = cpp11::as_cpp<cpp11::external_pointer<T>>(ptr).get();
  const int n_streams = rng->size();
  cpp11::writable::doubles ret = cpp11::writable::doubles(n * n_streams);
  double * y = REAL(ret);

  const double * mean = REAL(r_mean);
  const double * sd = REAL(r_sd);
  auto mean_vary = check_input_type(r_mean, n, n_streams, "mean");
  auto sd_vary = check_input_type(r_sd, n, n_streams, "sd");

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
  for (int i = 0; i < n_streams; ++i) {
    auto &state = rng->state(i);
    auto y_i = y + n * i;
    auto mean_i = mean_vary.generator ? mean + mean_vary.offset * i : mean;
    auto sd_i = sd_vary.generator ? sd + sd_vary.offset * i : sd;
    for (size_t j = 0; j < (size_t)n; ++j) {
      auto mean_ij = mean_vary.draw ? mean_i[j] : mean_i[0];
      auto sd_ij = sd_vary.draw ? sd_i[j] : sd_i[0];
      y_i[j] = mcstate::random::normal<real_type, A>(state, mean_ij, sd_ij);
    }
  }

  return sexp_matrix(ret, n, n_streams);
}

template <typename real_type, typename T>
cpp11::sexp mcstate_rng_binomial(SEXP ptr, int n,
                              cpp11::doubles r_size, cpp11::doubles r_prob,
                              int n_threads) {
  T *rng = cpp11::as_cpp<cpp11::external_pointer<T>>(ptr).get();
  const int n_streams = rng->size();
  cpp11::writable::doubles ret = cpp11::writable::doubles(n * n_streams);
  double * y = REAL(ret);

  const double * size = REAL(r_size);
  const double * prob = REAL(r_prob);
  auto size_vary = check_input_type(r_size, n, n_streams, "size");
  auto prob_vary = check_input_type(r_prob, n, n_streams, "prob");

  mcstate::utils::openmp_errors errors(n_streams);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
  for (int i = 0; i < n_streams; ++i) {
    try {
      auto &state = rng->state(i);
      auto y_i = y + n * i;
      auto size_i = size_vary.generator ? size + size_vary.offset * i : size;
      auto prob_i = prob_vary.generator ? prob + prob_vary.offset * i : prob;
      for (size_t j = 0; j < (size_t)n; ++j) {
        auto size_ij = size_vary.draw ? size_i[j] : size_i[0];
        auto prob_ij = prob_vary.draw ? prob_i[j] : prob_i[0];
        y_i[j] = mcstate::random::binomial<real_type>(state, size_ij, prob_ij);
      }
    } catch (std::exception const& e) {
      errors.capture(e, i);
    }
  }

  errors.report("generators", 4, true);

  return sexp_matrix(ret, n, n_streams);
}


template <typename real_type, typename T>
cpp11::sexp mcstate_rng_nbinomial(SEXP ptr, int n,
                              cpp11::doubles r_size, cpp11::doubles r_prob,
                              int n_threads) {
  T *rng = cpp11::as_cpp<cpp11::external_pointer<T>>(ptr).get();
  const int n_streams = rng->size();
  cpp11::writable::doubles ret = cpp11::writable::doubles(n * n_streams);
  double * y = REAL(ret);

  const double * size = REAL(r_size);
  const double * prob = REAL(r_prob);
  auto size_vary = check_input_type(r_size, n, n_streams, "size");
  auto prob_vary = check_input_type(r_prob, n, n_streams, "prob");

  mcstate::utils::openmp_errors errors(n_streams);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
  for (int i = 0; i < n_streams; ++i) {
    try {
      auto &state = rng->state(i);
      auto y_i = y + n * i;
      auto size_i = size_vary.generator ? size + size_vary.offset * i : size;
      auto prob_i = prob_vary.generator ? prob + prob_vary.offset * i : prob;
      for (size_t j = 0; j < (size_t)n; ++j) {
        auto size_ij = size_vary.draw ? size_i[j] : size_i[0];
        auto prob_ij = prob_vary.draw ? prob_i[j] : prob_i[0];
        y_i[j] = mcstate::random::nbinomial<real_type>(state, size_ij, prob_ij);
      }
    } catch (std::exception const& e) {
      errors.capture(e, i);
    }
  }

  errors.report("generators", 4, true);

  return sexp_matrix(ret, n, n_streams);
}

template <typename real_type, typename T>
cpp11::sexp mcstate_rng_poisson(SEXP ptr, int n, cpp11::doubles r_lambda,
                             int n_threads) {
  T *rng = cpp11::as_cpp<cpp11::external_pointer<T>>(ptr).get();
  const int n_streams = rng->size();
  cpp11::writable::doubles ret = cpp11::writable::doubles(n * n_streams);
  double * y = REAL(ret);

  const double * lambda = REAL(r_lambda);
  auto lambda_vary = check_input_type(r_lambda, n, n_streams, "lambda");

  mcstate::utils::openmp_errors errors(n_streams);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
  for (int i = 0; i < n_streams; ++i) {
    try {
      auto &state = rng->state(i);
      auto y_i = y + n * i;
      auto lambda_i = lambda_vary.generator ? lambda + lambda_vary.offset * i :
        lambda;
      for (size_t j = 0; j < (size_t)n; ++j) {
        auto lambda_ij = lambda_vary.draw ? lambda_i[j] : lambda_i[0];
        y_i[j] = mcstate::random::poisson<real_type>(state, lambda_ij);
      }
    } catch (std::exception const& e) {
      errors.capture(e, i);
    }
  }

  errors.report("generators", 4, true);

  return sexp_matrix(ret, n, n_streams);
}

template <typename real_type, typename T>
cpp11::sexp mcstate_rng_multinomial(SEXP ptr, int n,
                                 cpp11::doubles r_size,
                                 cpp11::doubles r_prob,
                                 int n_threads) {
  T *rng = cpp11::as_cpp<cpp11::external_pointer<T>>(ptr).get();
  const int n_streams = rng->size();

  const double * size = REAL(r_size);
  const double * prob = REAL(r_prob);
  auto size_vary = check_input_type(r_size, n, n_streams, "size");
  auto prob_vary = check_input_type2(r_prob, n, n_streams, "prob");
  const int len = prob_vary.len;

  // Normally we return a block of doubles with the first 'n' entries
  // being the results for the first generator, the second 'n' for the
  // second, and so on. Here the first n * len are the first generator
  // (with the first 'len' being the first sample.
  cpp11::writable::doubles ret =
    cpp11::writable::doubles(len * n * n_streams);
  double * y = REAL(ret);

  mcstate::utils::openmp_errors errors(n_streams);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
  for (int i = 0; i < n_streams; ++i) {
    try {
      auto &state = rng->state(i);
      auto y_i = y + len * n * i;
      auto size_i = size_vary.generator ? size + size_vary.offset * i : size;
      auto prob_i = prob_vary.generator ? prob + prob_vary.offset * i : prob;
      for (size_t j = 0; j < (size_t)n; ++j) {
        auto size_ij = size_vary.draw ? size_i[j]        : size_i[0];
        auto prob_ij = prob_vary.draw ? prob_i + j * len : prob_i;
        auto y_ij = y_i + j * len;
        mcstate::random::multinomial<real_type>(state, size_ij, prob_ij, len,
                                             y_ij);
      }
    } catch (std::exception const& e) {
      errors.capture(e, i);
    }
  }
  errors.report("generators", 4, true);

  if (n_streams == 1) {
    ret.attr("dim") = cpp11::writable::integers{len, n};
  } else {
    ret.attr("dim") = cpp11::writable::integers{len, n, n_streams};
  }
  return ret;
}

template <typename real_type, typename T>
cpp11::sexp mcstate_rng_hypergeometric(SEXP ptr, int n,
                                    cpp11::doubles r_n1, cpp11::doubles r_n2,
                                    cpp11::doubles r_k, int n_threads) {
  T *rng = cpp11::as_cpp<cpp11::external_pointer<T>>(ptr).get();
  const int n_streams = rng->size();
  cpp11::writable::doubles ret = cpp11::writable::doubles(n * n_streams);
  double * y = REAL(ret);

  const double * n1 = REAL(r_n1);
  const double * n2 = REAL(r_n2);
  const double * k = REAL(r_k);
  auto n1_vary = check_input_type(r_n1, n, n_streams, "n1");
  auto n2_vary = check_input_type(r_n2, n, n_streams, "n1");
  auto k_vary = check_input_type(r_k, n, n_streams, "k");

  mcstate::utils::openmp_errors errors(n_streams);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
  for (int i = 0; i < n_streams; ++i) {
    try {
      auto &state = rng->state(i);
      auto y_i = y + n * i;
      auto n1_i = n1_vary.generator ? n1 + n1_vary.offset * i : n1;
      auto n2_i = n2_vary.generator ? n2 + n2_vary.offset * i : n2;
      auto k_i = k_vary.generator ? k + k_vary.offset * i : k;
      for (size_t j = 0; j < (size_t)n; ++j) {
        auto n1_ij = n1_vary.draw ? n1_i[j] : n1_i[0];
        auto n2_ij = n2_vary.draw ? n2_i[j] : n2_i[0];
        auto k_ij = k_vary.draw ? k_i[j] : k_i[0];
        y_i[j] = mcstate::random::hypergeometric<real_type>(state, n1_ij, n2_ij, k_ij);
      }
    } catch (std::exception const& e) {
      errors.capture(e, i);
    }
  }

  errors.report("generators", 4, true);

  return sexp_matrix(ret, n, n_streams);
}


template <typename real_type, typename T>
cpp11::sexp mcstate_rng_gamma_scale(SEXP ptr, int n,
                                    cpp11::doubles r_shape,
                                    cpp11::doubles r_scale,
                                    int n_threads) {
  T *rng = cpp11::as_cpp<cpp11::external_pointer<T>>(ptr).get();
  const int n_streams = rng->size();
  cpp11::writable::doubles ret = cpp11::writable::doubles(n * n_streams);
  double * y = REAL(ret);

  const double * shape = REAL(r_shape);
  const double * scale = REAL(r_scale);
  auto shape_vary = check_input_type(r_shape, n, n_streams, "shape");
  auto scale_vary = check_input_type(r_scale, n, n_streams, "scale");

  mcstate::utils::openmp_errors errors(n_streams);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
  for (int i = 0; i < n_streams; ++i) {
    try {
      auto &state = rng->state(i);
      auto y_i = y + n * i;
      auto shape_i = shape_vary.generator ? shape + shape_vary.offset * i : shape;
      auto scale_i = scale_vary.generator ? scale + scale_vary.offset * i : scale;
      for (size_t j = 0; j < (size_t)n; ++j) {
        auto shape_ij = shape_vary.draw ? shape_i[j] : shape_i[0];
        auto scale_ij = scale_vary.draw ? scale_i[j] : scale_i[0];
        y_i[j] = mcstate::random::gamma_scale<real_type>(state, shape_ij, scale_ij);
      }
    } catch (std::exception const& e) {
      errors.capture(e, i);
    }
  }

  errors.report("generators", 4, true);

  return sexp_matrix(ret, n, n_streams);
}

template <typename real_type, typename T>
cpp11::sexp mcstate_rng_gamma_rate(SEXP ptr, int n,
                                    cpp11::doubles r_shape,
                                    cpp11::doubles r_rate,
                                    int n_threads) {
  T *rng = cpp11::as_cpp<cpp11::external_pointer<T>>(ptr).get();
  const int n_streams = rng->size();
  cpp11::writable::doubles ret = cpp11::writable::doubles(n * n_streams);
  double * y = REAL(ret);

  const double * shape = REAL(r_shape);
  const double * rate = REAL(r_rate);
  auto shape_vary = check_input_type(r_shape, n, n_streams, "shape");
  auto rate_vary = check_input_type(r_rate, n, n_streams, "rate");

  mcstate::utils::openmp_errors errors(n_streams);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
  for (int i = 0; i < n_streams; ++i) {
    try {
      auto &state = rng->state(i);
      auto y_i = y + n * i;
      auto shape_i = shape_vary.generator ? shape + shape_vary.offset * i : shape;
      auto rate_i = rate_vary.generator ? rate + rate_vary.offset * i : rate;
      for (size_t j = 0; j < (size_t)n; ++j) {
        auto shape_ij = shape_vary.draw ? shape_i[j] : shape_i[0];
        auto rate_ij = rate_vary.draw ? rate_i[j] : rate_i[0];
        y_i[j] = mcstate::random::gamma_rate<real_type>(state, shape_ij, rate_ij);
      }
    } catch (std::exception const& e) {
      errors.capture(e, i);
    }
  }

  errors.report("generators", 4, true);

  return sexp_matrix(ret, n, n_streams);
}

// Below here is very repetitive, and could probably be deduplicated
// with some clever template magic. Most of the faff is because we
// want to support 4 modes of taking 1 or 2 parameters (each varying
// or not over draws and generators)
template <typename real_type, typename T>
cpp11::sexp mcstate_rng_cauchy(SEXP ptr, int n,
                            cpp11::doubles r_location,
                            cpp11::doubles r_scale,
                            int n_threads) {
  T *rng = cpp11::as_cpp<cpp11::external_pointer<T>>(ptr).get();
  const int n_streams = rng->size();
  cpp11::writable::doubles ret = cpp11::writable::doubles(n * n_streams);
  double * y = REAL(ret);

  const double * location = REAL(r_location);
  const double * scale = REAL(r_scale);
  auto location_vary = check_input_type(r_location, n, n_streams, "location");
  auto scale_vary = check_input_type(r_scale, n, n_streams, "scale");

  mcstate::utils::openmp_errors errors(n_streams);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
  for (int i = 0; i < n_streams; ++i) {
    try {
      auto &state = rng->state(i);
      auto y_i = y + n * i;
      auto location_i = location_vary.generator ? location + location_vary.offset * i : location;
      auto scale_i = scale_vary.generator ? scale + scale_vary.offset * i : scale;
      for (size_t j = 0; j < (size_t)n; ++j) {
        auto location_ij = location_vary.draw ? location_i[j] : location_i[0];
        auto scale_ij = scale_vary.draw ? scale_i[j] : scale_i[0];
        y_i[j] = mcstate::random::cauchy<real_type>(state, location_ij, scale_ij);
      }
    } catch (std::exception const& e) {
      errors.capture(e, i);
    }
  }

  errors.report("generators", 4, true);

  return sexp_matrix(ret, n, n_streams);
}

template <typename real_type, typename T>
cpp11::sexp mcstate_rng_beta(SEXP ptr, int n,
			     cpp11::doubles r_a,
			     cpp11::doubles r_b,
			     int n_threads) {
  T *rng = cpp11::as_cpp<cpp11::external_pointer<T>>(ptr).get();
  const int n_streams = rng->size();
  cpp11::writable::doubles ret = cpp11::writable::doubles(n * n_streams);
  double * y = REAL(ret);

  const double * a = REAL(r_a);
  const double * b = REAL(r_b);
  auto a_vary = check_input_type(r_a, n, n_streams, "a");
  auto b_vary = check_input_type(r_b, n, n_streams, "b");

  mcstate::utils::openmp_errors errors(n_streams);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
  for (int i = 0; i < n_streams; ++i) {
    try {
      auto &state = rng->state(i);
      auto y_i = y + n * i;
      auto a_i = a_vary.generator ? a + a_vary.offset * i : a;
      auto b_i = b_vary.generator ? b + b_vary.offset * i : b;
      for (size_t j = 0; j < (size_t)n; ++j) {
        auto a_ij = a_vary.draw ? a_i[j] : a_i[0];
        auto b_ij = b_vary.draw ? b_i[j] : b_i[0];
        y_i[j] = mcstate::random::beta<real_type>(state, a_ij, b_ij);
      }
    } catch (std::exception const& e) {
      errors.capture(e, i);
    }
  }

  errors.report("generators", 4, true);

  return sexp_matrix(ret, n, n_streams);
}

template <typename T>
cpp11::sexp mcstate_rng_state(SEXP ptr) {
  T *rng = cpp11::as_cpp<cpp11::external_pointer<T>>(ptr).get();
  auto state = rng->export_state();
  size_t len = sizeof(typename T::int_type) * state.size();
  cpp11::writable::raws ret(len);
  std::memcpy(RAW(ret), state.data(), len);
  return ret;
}

[[cpp11::register]]
SEXP mcstate_rng_alloc(cpp11::sexp r_seed, int n_streams, bool deterministic,
                    bool is_float) {
  return is_float ?
    mcstate_rng_alloc<default_rng32>(r_seed, n_streams, deterministic) :
    mcstate_rng_alloc<default_rng64>(r_seed, n_streams, deterministic);
}

[[cpp11::register]]
void mcstate_rng_jump(SEXP ptr, bool is_float) {
  if (is_float) {
    mcstate_rng_jump<default_rng32>(ptr);
  } else {
    mcstate_rng_jump<default_rng64>(ptr);
  }
}

[[cpp11::register]]
void mcstate_rng_long_jump(SEXP ptr, bool is_float) {
  if (is_float) {
    mcstate_rng_long_jump<default_rng32>(ptr);
  } else {
    mcstate_rng_long_jump<default_rng64>(ptr);
  }
}

[[cpp11::register]]
cpp11::sexp mcstate_rng_random_real(SEXP ptr, int n, int n_threads,
                                 bool is_float) {
  return is_float ?
    mcstate_rng_random_real<float, default_rng32>(ptr, n, n_threads) :
    mcstate_rng_random_real<double, default_rng64>(ptr, n, n_threads);
}

[[cpp11::register]]
cpp11::sexp mcstate_rng_random_normal(SEXP ptr, int n, int n_threads,
                                   std::string algorithm, bool is_float) {
  cpp11::sexp ret;
  if (algorithm == "box_muller") {
    constexpr auto a = mcstate::random::algorithm::normal::box_muller;
    ret = is_float ?
      mcstate_rng_random_normal<float, a, default_rng32>(ptr, n, n_threads) :
      mcstate_rng_random_normal<double, a, default_rng64>(ptr, n, n_threads);
  } else if (algorithm == "polar") {
    constexpr auto a = mcstate::random::algorithm::normal::polar;
    ret = is_float ?
      mcstate_rng_random_normal<float, a, default_rng32>(ptr, n, n_threads) :
      mcstate_rng_random_normal<double, a, default_rng64>(ptr, n, n_threads);
  } else if (algorithm == "ziggurat") {
    constexpr auto a = mcstate::random::algorithm::normal::ziggurat;
    ret = is_float ?
      mcstate_rng_random_normal<float, a, default_rng32>(ptr, n, n_threads) :
      mcstate_rng_random_normal<double, a, default_rng64>(ptr, n, n_threads);
  } else {
    cpp11::stop("Unknown normal algorithm '%s'", algorithm.c_str());
  }
  return ret;
}

[[cpp11::register]]
cpp11::sexp mcstate_rng_uniform(SEXP ptr, int n,
                             cpp11::doubles r_min,
                             cpp11::doubles r_max,
                             int n_threads,
                             bool is_float) {
  return is_float ?
    mcstate_rng_uniform<float, default_rng32>(ptr, n, r_min, r_max, n_threads) :
    mcstate_rng_uniform<double, default_rng64>(ptr, n, r_min, r_max, n_threads);
}

[[cpp11::register]]
cpp11::sexp mcstate_rng_exponential_rate(SEXP ptr, int n, cpp11::doubles r_rate,
                                         int n_threads,
                                         bool is_float) {
  return is_float ?
    mcstate_rng_exponential_rate<float, default_rng32>(ptr, n, r_rate, n_threads) :
    mcstate_rng_exponential_rate<double, default_rng64>(ptr, n, r_rate, n_threads);
}

[[cpp11::register]]
cpp11::sexp mcstate_rng_exponential_mean(SEXP ptr, int n, cpp11::doubles r_mean,
                                         int n_threads,
                                         bool is_float) {
  return is_float ?
    mcstate_rng_exponential_mean<float, default_rng32>(ptr, n, r_mean, n_threads) :
    mcstate_rng_exponential_mean<double, default_rng64>(ptr, n, r_mean, n_threads);
}

[[cpp11::register]]
cpp11::sexp mcstate_rng_normal(SEXP ptr, int n, cpp11::doubles r_mean,
                            cpp11::doubles r_sd, int n_threads,
                            std::string algorithm, bool is_float) {
  cpp11::sexp ret;
  if (algorithm == "box_muller") {
    constexpr auto a = mcstate::random::algorithm::normal::box_muller;
    ret = is_float ?
      mcstate_rng_normal<float, a, default_rng32>(ptr, n, r_mean, r_sd, n_threads) :
      mcstate_rng_normal<double, a, default_rng64>(ptr, n, r_mean, r_sd, n_threads);
  } else if (algorithm == "polar") {
    constexpr auto a = mcstate::random::algorithm::normal::polar;
    ret = is_float ?
      mcstate_rng_normal<float, a, default_rng32>(ptr, n, r_mean, r_sd, n_threads) :
      mcstate_rng_normal<double, a, default_rng64>(ptr, n, r_mean, r_sd, n_threads);
  } else if (algorithm == "ziggurat") {
    constexpr auto a = mcstate::random::algorithm::normal::ziggurat;
    ret = is_float ?
      mcstate_rng_normal<float, a, default_rng32>(ptr, n, r_mean, r_sd, n_threads) :
      mcstate_rng_normal<double, a, default_rng64>(ptr, n, r_mean, r_sd, n_threads);
  } else {
    cpp11::stop("Unknown normal algorithm '%s'", algorithm.c_str());
  }
  return ret;
}

[[cpp11::register]]
cpp11::sexp mcstate_rng_binomial(SEXP ptr, int n,
                              cpp11::doubles r_size, cpp11::doubles r_prob,
                              int n_threads, bool is_float) {
  return is_float ?
    mcstate_rng_binomial<float, default_rng32>(ptr, n, r_size, r_prob, n_threads) :
    mcstate_rng_binomial<double, default_rng64>(ptr, n, r_size, r_prob, n_threads);
}

[[cpp11::register]]
cpp11::sexp mcstate_rng_nbinomial(SEXP ptr, int n,
                              cpp11::doubles r_size, cpp11::doubles r_prob,
                              int n_threads, bool is_float) {
  return is_float ?
    mcstate_rng_nbinomial<float, default_rng32>(ptr, n, r_size, r_prob, n_threads) :
    mcstate_rng_nbinomial<double, default_rng64>(ptr, n, r_size, r_prob, n_threads);
}

[[cpp11::register]]
cpp11::sexp mcstate_rng_hypergeometric(SEXP ptr, int n,
                                    cpp11::doubles r_n1,
                                    cpp11::doubles r_n2,
                                    cpp11::doubles r_k,
                                    int n_threads, bool is_float) {
  return is_float ?
    mcstate_rng_hypergeometric<float, default_rng32>(ptr, n, r_n1, r_n2, r_k, n_threads) :
    mcstate_rng_hypergeometric<double, default_rng64>(ptr, n, r_n1, r_n2, r_k, n_threads);
}

[[cpp11::register]]
cpp11::sexp mcstate_rng_gamma_scale(SEXP ptr, int n,
                                    cpp11::doubles r_shape,
                                    cpp11::doubles r_scale,
                                    int n_threads, bool is_float) {
  return is_float ?
    mcstate_rng_gamma_scale<float, default_rng32>(ptr, n, r_shape, r_scale, n_threads) :
    mcstate_rng_gamma_scale<double, default_rng64>(ptr, n, r_shape, r_scale, n_threads);
}

[[cpp11::register]]
cpp11::sexp mcstate_rng_gamma_rate(SEXP ptr, int n,
                                    cpp11::doubles r_shape,
                                    cpp11::doubles r_rate,
                                    int n_threads, bool is_float) {
  return is_float ?
    mcstate_rng_gamma_rate<float, default_rng32>(ptr, n, r_shape, r_rate, n_threads) :
    mcstate_rng_gamma_rate<double, default_rng64>(ptr, n, r_shape, r_rate, n_threads);
}

[[cpp11::register]]
cpp11::sexp mcstate_rng_poisson(SEXP ptr, int n,
                             cpp11::doubles r_lambda,
                             int n_threads, bool is_float) {
  return is_float ?
    mcstate_rng_poisson<float, default_rng32>(ptr, n, r_lambda, n_threads) :
    mcstate_rng_poisson<double, default_rng64>(ptr, n, r_lambda, n_threads);
}

[[cpp11::register]]
cpp11::sexp mcstate_rng_cauchy(SEXP ptr, int n,
                            cpp11::doubles r_location,
                            cpp11::doubles r_scale,
                            int n_threads,
                            bool is_float) {
  return is_float ?
    mcstate_rng_cauchy<float, default_rng32>(ptr, n, r_location, r_scale, n_threads) :
    mcstate_rng_cauchy<double, default_rng64>(ptr, n, r_location, r_scale, n_threads);
}

[[cpp11::register]]
cpp11::sexp mcstate_rng_beta(SEXP ptr, int n,
			     cpp11::doubles r_a,
			     cpp11::doubles r_b,
			     int n_threads, bool is_float) {
  return is_float ?
    mcstate_rng_beta<float, default_rng32>(ptr, n, r_a, r_b, n_threads) :
    mcstate_rng_beta<double, default_rng64>(ptr, n, r_a, r_b, n_threads);
}

[[cpp11::register]]
cpp11::sexp mcstate_rng_multinomial(SEXP ptr, int n,
                                 cpp11::doubles r_size, cpp11::doubles r_prob,
                                 int n_threads, bool is_float) {
  return is_float ?
    mcstate_rng_multinomial<float, default_rng32>(ptr, n, r_size, r_prob, n_threads) :
    mcstate_rng_multinomial<double, default_rng64>(ptr, n, r_size, r_prob, n_threads);
}

[[cpp11::register]]
cpp11::sexp mcstate_rng_state(SEXP ptr, bool is_float) {
  return is_float ?
    mcstate_rng_state<default_rng32>(ptr) :
    mcstate_rng_state<default_rng64>(ptr);
}
