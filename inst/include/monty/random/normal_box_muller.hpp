#pragma once

#include <cmath>

// Stan prevents this constant being defined on some systems; could
// pull it in from Rmath.h too, but it seems unlikely to change...
#ifndef M_PI
#define M_PI = 3.14159265358979
#endif

#include "monty/random/generator.hpp"
#include "monty/random/math.hpp"

namespace monty {
namespace random {

__nv_exec_check_disable__
template <typename real_type, typename rng_state_type>
__host__ __device__
real_type random_normal_box_muller(rng_state_type& rng_state) {
  // This function implements the Box-Muller transform:
  // https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform#Basic_form
  // Do not send a really small number to log().
  const real_type epsilon = utils::epsilon<real_type>();
  const real_type two_pi = 2 * M_PI;

  real_type u1, u2;
  do {
    u1 = random_real<real_type>(rng_state);
    u2 = random_real<real_type>(rng_state);
  } while (u1 <= epsilon);

  SYNCWARP
  return monty::math::sqrt(-2 * monty::math::log(u1)) * std::cos(two_pi * u2);
}


}
}
