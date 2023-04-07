#ifndef CUDA_KERNEL_H
#define CUDA_KERNEL_H

#include <cuda_runtime.h>
#include <cstdint>

__global__ 
void compute_acceleration(base_particle other);

__global__
void compute_acceleration(Container particles);

#endif
