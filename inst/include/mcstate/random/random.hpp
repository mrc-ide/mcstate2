#pragma once

#include "mcstate/random/generator.hpp"
#include "mcstate/random/prng.hpp"

#include "mcstate/random/beta.hpp"
#include "mcstate/random/binomial.hpp"
#include "mcstate/random/cauchy.hpp"
#include "mcstate/random/exponential.hpp"
#include "mcstate/random/gamma.hpp"
#include "mcstate/random/hypergeometric.hpp"
#include "mcstate/random/multinomial.hpp"
#include "mcstate/random/nbinomial.hpp"
#include "mcstate/random/normal.hpp"
#include "mcstate/random/poisson.hpp"
#include "mcstate/random/uniform.hpp"

#include "mcstate/random/version.hpp"

namespace mcstate {
namespace random {

namespace {

template <typename T>
struct default_rng_helper;

template <>
struct default_rng_helper<double> {
  using type = xoshiro256plus;
};

template <>
struct default_rng_helper<float> {
  using type = xoshiro128plus;
};

}

template <typename T>
using generator = typename default_rng_helper<T>::type;

}
}
