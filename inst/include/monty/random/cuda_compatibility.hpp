#pragma once

#include <stdexcept>

// There are 3 different ways we hit this file:
//
// *  __CUDA_ARCH__ is defined: we're compiling under nvcc generating
//        device code. In this case __NVCC__ is always defined.
//
// * __NVCC__: we're compiling under nvcc, either generating host or
//        device code
//
// * Neither is defined: we're compiling under gcc/clang etc and need
//   to inject our stubs

// This is necessary due to templates which are __host__ __device__;
// whenever a HOSTDEVICE function is called from another HOSTDEVICE
// function the compiler gets confused as it can't tell which one it's
// going to use. This suppresses the warning as it is ok here.
#ifdef __NVCC__
#define __nv_exec_check_disable__ _Pragma("nv_exec_check_disable")
#else
#define __nv_exec_check_disable__

// Additional stubs used to shadow the now-unneeded directives for
// nvcc, which would cause gcc/clang to raise an error
#define __host__
#define __device__
#define __global__
#define __align__(n)
#endif

#ifdef __CUDA_ARCH__
// Compiling under nvcc for the device
#define CONSTANT __constant__
#define SYNCWARP __syncwarp();
#else
// gcc/clang or nvcc for the host
#define CONSTANT const
#define SYNCWARP
#endif

namespace monty {
namespace utils {

// We cannot throw errors in GPU code, we can only send a trap signal,
// which is unrecoverable.
__host__ __device__
inline void fatal_error(const char *message) {
#ifdef __CUDA_ARCH__
  printf(message);
  printf("\n");
  __trap();
#else
  throw std::runtime_error(message);
#endif
}

}
}
