#ifndef __GPU_FUNCTIONS_
#define __GPU_FUNCTIONS_

#include <cuda.h>


__global__ void solverInit(float *i_h_old, float *i_h_new, float *i_hu_old,
    float *i_hu_new, float *i_hv_old, float *i_hv_new,
    int i_nx, int i_ny, float *i_b, int i_iter,
    float i_scaling);

__device__ void netUpdates(float *i_height_old, float *i_height_new,
    float *i_momentum_old, float *i_momentum_new,
    int i_nx, int i_ny, float *i_b, float i_scaling,
    int idx, int i_stride);

__device__ void setGhostOutflow(float *i_height, float *i_hu, float *i_hv,
         float *i_b, int i_nx, int i_ny);

#endif
