#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "monty/random/generator.hpp"
#include "monty/random/numeric.hpp"
#include "monty/random/exponential.hpp"
#include "monty/random/uniform.hpp"
#include "monty/random/normal.hpp"
#include "monty/random/math.hpp"

// Algorithm from George Marsaglia and Wai Wan Tsang. 2000. "A Simple Method
// for Generating Gamma Variables" *ACM Trans. Math. Softw.* 26, 3 (September 2000),
// 363-372. DOI:[10.1145/358407.358414](https://doi.acm.org/10.1145/358407.358414)
// and follows the Rust implementation https://docs.rs/rand/0.5.0/src/rand/distributions/gamma.rs.html
// but adapted to fit our needs.
namespace monty {
namespace random {
namespace {

template <typename real_type>
void gamma_validate(real_type shape, real_type scale) {
  if (shape < 0.0 || scale < 0.0) {
    char buffer[256];
    snprintf(buffer, 256,
             "Invalid call to gamma with shape = %g, scale = %g",
             shape, scale);
    monty::utils::fatal_error(buffer);
  }
}

template <typename real_type, typename rng_state_type>
real_type gamma_large(rng_state_type& rng_state, real_type shape) {
  real_type d = shape - 1.0 / 3.0;
  real_type c = 1.0 / sqrt(9.0 * d);
  while (true) {
    real_type x = normal<real_type>(rng_state, 0, 1);
    real_type v_cbrt = 1.0 + c * x;
    if (v_cbrt <= 0.0) {
      continue;
    }
    real_type v = v_cbrt * v_cbrt * v_cbrt;
    real_type u = random_real<real_type>(rng_state);
    real_type x_sqr = x * x;
    if (u < 1.0 - 0.0331 * x_sqr * x_sqr ||
      monty::math::log(u) < 0.5 * x_sqr + d * (1.0 - v + monty::math::log(v))) {
      return d * v;
    }
  }
}

template <typename real_type, typename rng_state_type>
real_type gamma_small(rng_state_type& rng_state, real_type shape) {
  real_type inv_shape = 1 / shape;
  real_type u = random_real<real_type>(rng_state);
  return gamma_large(rng_state, shape + 1.0) * monty::math::pow(u, inv_shape);
}

template <typename real_type>
real_type gamma_deterministic(real_type shape, real_type scale) {
  return shape * scale;
}

}

/// Draw random number from the gamma distribution.
/// @tparam real_type The underlying real number type, typically
/// `double` or `float`. A compile-time error will be thrown if you
/// attempt to use a non-floating point type (based on
/// `std::is_floating_point).
///
/// @tparam rng_state_type The random number state type
///
/// @param rng_state Reference to the random number state, will be
/// modified as a side-effect
///
/// @param a Shape
///
/// @param b Scale
template <typename real_type, typename rng_state_type>
__host__ __device__
real_type gamma_scale(rng_state_type& rng_state, real_type shape, real_type scale) {
  static_assert(std::is_floating_point<real_type>::value,
                "Only valid for floating-point types; use gamma<real_type>()");

  gamma_validate(shape, scale);

#ifdef __CUDA_ARCH__
  static_assert("gamma() not implemented for GPU targets");
#endif

  if (shape == 0 || scale == 0) {
    return 0;
  }

  if (rng_state.deterministic) {
    return gamma_deterministic<real_type>(shape, scale);
  }

  if (shape < 1) {
    return gamma_small<real_type>(rng_state, shape) * scale;
  }

  if (shape == 1) {
    return exponential_mean(rng_state, scale);
  }

  return gamma_large<real_type>(rng_state, shape) * scale;
}

template <typename real_type, typename rng_state_type>
__host__ __device__
real_type gamma_rate(rng_state_type& rng_state, real_type shape, real_type rate) {
  static_assert(std::is_floating_point<real_type>::value,
                "Only valid for floating-point types; use gamma<real_type>()");
  const auto scale = 1 / rate;
  gamma_validate(shape, scale);
  return gamma_scale<real_type>(rng_state, shape, scale);
}

}
}
