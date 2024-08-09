#pragma once

#include <cmath>

#include "mcstate/random/gamma.hpp"
#include "mcstate/random/poisson.hpp"
#include "mcstate/random/generator.hpp"

namespace mcstate {
namespace random {

namespace {

template <typename real_type>
void nbinomial_validate(real_type size, real_type prob) {
   if(!R_FINITE(size) || !R_FINITE(prob) || size <= 0 || prob <= 0 || prob > 1) {
    char buffer[256];
    snprintf(buffer, 256,
             "Invalid call to nbinomial with size = %g, prob = %g",
             size, prob);
    mcstate::utils::fatal_error(buffer);
  }
}

}

template <typename real_type, typename rng_state_type>
real_type nbinomial(rng_state_type& rng_state, real_type size, real_type prob) {
#ifdef __CUDA_ARCH__
  static_assert("nbinomial() not implemented for GPU targets");
#endif
    nbinomial_validate(size, prob);

    if (rng_state.deterministic) {
      return (1 - prob) * size / prob;
    }
    return (prob == 1) ? 0 : poisson(rng_state, gamma_scale(rng_state, size, (1 - prob) / prob));
}

}
}

