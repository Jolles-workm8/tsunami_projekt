#ifndef __GPU_FUNCTIONS_
#define __GPU_FUNCTIONS_

#include <cuda.h>


__global__ void solverInit(float *i_h, float *l_h_UpadeteR, float *l_h_UpadeteL, float *i_hu,
    float *l_hu_UpadeteR, float *l_hu_UpadeteL, float *i_hv,
    int i_nx, int i_ny, float *i_b, int i_iter,
    float i_scaling);

__device__ void netUpdates(float *i_height, float *o_height_UpdateR, float *o_height_UpdateL,
    float *i_momentum, float *o_momentum_UpdateR, float *o_momentum_UpdateL, 
    int i_xEdges, int i_yEdges, float *i_b, float i_scaling,
    int i_nx, int i_stride);

__device__ void setGhostOutflow(float *i_height, float *i_hu, float *i_hv,
         float *i_b, int i_nx, int i_ny);

__device__ void updateValues(float *o_h, float *i_h_UpdateR, float *i_h_UpdateL, float *o_hu,
                           float *i_hu_UpdateR, float *i_hu_UpdateL, int i_nx, int i_ny);

__device__ void initUpdatesZero(float *i_h_UpdateR, float *i_h_UpdateL, float *i_hu_UpdateR, float *i_hu_UpdateL, int i_nx, int i_ny);

#endif
