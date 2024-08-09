#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "mcstate/random/gamma.hpp"

namespace mcstate {
namespace random {

template <typename real_type, typename rng_state_type>
__host__ __device__
real_type beta(rng_state_type& rng_state, real_type a, real_type b) {
  static_assert(std::is_floating_point<real_type>::value,
                "Only valid for floating-point types; use beta<real_type>()");
#ifdef __CUDA_ARCH__
  static_assert("gamma() not implemented for GPU targets");
#endif

  if (rng_state.deterministic) {
    return a / (a + b);
  }

  const auto x = gamma_scale<real_type>(rng_state, a, 1);
  const auto y = gamma_scale<real_type>(rng_state, b, 1);
  return x / (x + y);
}

}
}
