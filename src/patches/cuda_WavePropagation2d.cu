#include "cuda_WavePropagation2d.h"
#include "cuda_functions.cuh"
//#include "cuda_solver.cuh"
#include <iostream>
#include <cuda.h>

tsunami_lab::patches::cuda_WavePropagation2d::cuda_WavePropagation2d(
    t_idx i_xCells, t_idx i_yCells) {
  // set number of cells without GhostCells in x and y Direction
  m_xCells = i_xCells;
  m_yCells = i_yCells;

  size = (m_xCells + 2) * (m_yCells + 2);
  // allocate memory including ghostcells on each side
  m_h = (float *)malloc(size * sizeof(float));
  m_hu = (float *)malloc(size * sizeof(float));
  m_hv = (float *)malloc(size * sizeof(float));
  m_b = (float *)malloc(size * sizeof(float));

  // init to zero
  for (unsigned long l_ceY = 0; l_ceY < (m_yCells + 2); l_ceY++) {
    for (unsigned long l_ceX = 0; l_ceX < (m_xCells + 2); l_ceX++) {
      unsigned long l_ce = l_ceX + l_ceY * (m_xCells + 2);
      m_b[l_ce] = 0;
      m_h[l_ce] = 0;
      m_hu[l_ce] = 0;
      m_hv[l_ce] = 0;
    }
  }

  // allocate memory on GPU
  cudaMalloc((void **)&h_dev_old, size * sizeof(float));
  cudaMalloc((void **)&h_dev_new, size * sizeof(float));
  cudaMalloc((void **)&hu_dev_old, size * sizeof(float));
  cudaMalloc((void **)&hu_dev_new, size * sizeof(float));
  cudaMalloc((void **)&hv_dev_old, size * sizeof(float));
  cudaMalloc((void **)&hv_dev_new, size * sizeof(float));
  cudaMalloc((void **)&b_dev, size * sizeof(float));
}

tsunami_lab::patches::cuda_WavePropagation2d::~cuda_WavePropagation2d() {
  // free memory on CPU
  free(m_h);
  free(m_hu);
  free(m_hv);
  free(m_b);

  // free memory on GPU
  cudaFree(h_dev_new);
  cudaFree(h_dev_old);
  cudaFree(hu_dev_new);
  cudaFree(hu_dev_old);
  cudaFree(hv_dev_new);
  cudaFree(hv_dev_old);
  cudaFree(b_dev);
}
void tsunami_lab::patches::cuda_WavePropagation2d::MemTransfer() {
  cudaMemcpy(h_dev_old, m_h, size * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(hu_dev_old, m_hu, size * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(hv_dev_old, m_hv, size * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(b_dev, m_b, size * sizeof(float), cudaMemcpyHostToDevice);
}
void tsunami_lab::patches::cuda_WavePropagation2d::timeStep(
    t_real i_scaling, t_idx i_computeSteps) {
  t_idx l_nx = m_xCells + 2;
  t_idx l_ny = m_yCells + 2;

  dim3 threadsPerBlock(4, 4);
  dim3 BlocksPerGrid(l_nx / 4, l_ny / 4);

  solverInit<<<BlocksPerGrid, threadsPerBlock>>>(
      h_dev_old, h_dev_new, hu_dev_old, hu_dev_new, hv_dev_old, hv_dev_new,
      l_nx, l_ny, b_dev, i_computeSteps, i_scaling);

  cudaDeviceSynchronize();

  cudaMemcpy(m_h, h_dev_old, size * sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(m_hu, hu_dev_old, size * sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(m_hv, hv_dev_old, size * sizeof(float), cudaMemcpyDeviceToHost);
}

__global__ void solverInit(float *i_h_old, float *i_h_new, float *i_hu_old,
                           float *i_hu_new, float *i_hv_old, float *i_hv_new,
                           int i_nx, int i_ny, float *i_b, int i_iter,
                           float i_scaling) {
  int l_i = blockIdx.x * blockDim.x + threadIdx.x;
  int l_j = blockIdx.y * blockDim.y + threadIdx.y;
  int idx = l_i + l_j * i_nx;
  for (int i = 0; i < i_iter; ++i) {
    setGhostOutflow(i_h_old, i_hu_old, i_hv_old, i_b, i_nx, i_ny);
    __syncthreads();

    netUpdates(i_h_old, i_h_new, i_hu_old, i_hu_new, i_nx - 1, i_ny, i_b,
               i_scaling, idx, 1);
    __syncthreads();

    netUpdates(i_h_old, i_h_new, i_hv_old, i_hv_new, i_nx, i_ny - 1, i_b,
               i_scaling, idx, i_nx);
    __syncthreads();
  }
}

__device__ void setGhostOutflow(float *i_height, float *i_hu, float *i_hv,
                                float *i_b, int i_nx, int i_ny) {
  int l_i = blockIdx.x * blockDim.x + threadIdx.x;
  int l_j = blockIdx.y * blockDim.y + threadIdx.y;
  int idx = l_i + l_j * i_nx;
  
  if (l_i < i_nx && l_j < i_ny){
    if (l_j == 0) {
      i_height[idx] = i_height[idx + i_nx];
      i_hu[idx] = i_hu[idx + i_nx];
      i_hv[idx] = i_hv[idx + i_nx];
      i_b[idx] = i_b[idx + i_nx];
    } else if (l_j == i_nx - 1) {
      i_height[idx] = i_height[idx - i_nx];
      i_hu[idx] = i_hu[idx - i_nx];
      i_hv[idx] = i_hv[idx - i_nx];
      i_b[idx] = i_b[idx - i_nx];
    }

    __syncthreads();

    if (l_i == 0) {
      i_height[idx] = i_height[idx + 1];
      i_hu[idx] = i_hu[idx + 1];
      i_hv[idx] = i_hv[idx + 1];
      i_b[idx] = i_b[idx + 1];
    } else if (l_i == i_ny - 1) {
      i_height[idx] = i_height[idx - 1];
      i_hu[idx] = i_hu[idx - 1];
      i_hv[idx] = i_hv[idx - 1];
      i_b[idx] = i_b[idx - 1];
    }
  }
}

__device__ void netUpdates(float *i_height_old, float *i_height_new,
                           float *i_momentum_old, float *i_momentum_new,
                           int i_nx, int i_ny, float *i_b, float i_scaling,
                           int idx, int i_stride) {
  int l_i = blockIdx.x * blockDim.x + threadIdx.x;
  int l_j = blockIdx.y * blockDim.y + threadIdx.y;

  float m_gSqrt = sqrtf(9.812);
  float m_g = 9.812;

  i_height_new[idx] = i_height_old[idx];
  i_momentum_new[idx] = i_momentum_old[idx];

  if (l_i < i_nx-1 && l_j < i_ny-1) {
    // compute u for left and right
    float l_uL = i_momentum_old[idx] / i_height_old[idx];
    float l_uR = i_momentum_old[idx + i_stride] / i_height_old[idx + i_stride];

    float l_hL = i_height_old[idx];
    float l_hR = i_height_old[idx + i_stride];

    float l_huL = i_momentum_old[idx];
    float l_huR = i_momentum_old[idx + i_stride];

    float l_bL = i_b[idx];
    float l_bR = i_b[idx + i_stride];


    if(l_bL >= 0){
      l_hL = l_hR;
      l_huL = l_huR;
      l_bL = l_bR;
    }
    else if(l_bR >= 0){
      l_hR = l_hL;
      l_huR = l_huL;
      l_bR = l_bL;
    }


    // compute WaveSpeed ,

    float l_hSqrtL = sqrtf(l_hL);
    float l_hSqrtR = sqrtf(l_hR);

    float l_hRoe = 0.5f * (l_hL + l_hR);
    float l_uRoe = l_hSqrtL * l_uL + l_hSqrtR * l_uR;
    l_uRoe /= l_hSqrtL + l_hSqrtR;

    float l_ghSqrtRoe = m_gSqrt * sqrtf(l_hRoe);

    float l_waveSpeedL = l_uRoe - l_ghSqrtRoe;
    float l_waveSpeedR = l_uRoe + l_ghSqrtRoe;

    float l_detInv = 1 / (l_waveSpeedR - l_waveSpeedL);

    // compute the bathymetry effect
    float l_bathEff = -m_g * (l_bR - l_bL) * (l_hL + l_hR) / 2;

    // compute jump in the flux
    float l_fJump_1 = l_huR - l_huL;
    float l_fJump_2 = l_huR * l_huR / l_hR - l_huL * l_huL / l_hL +
                      (m_g / 2) * (l_hR * l_hR - l_hL * l_hL);
    l_fJump_2 -= l_bathEff;

    // compute the alpha values
    float l_strengthL =
        -i_scaling * l_detInv * (l_waveSpeedR * l_fJump_1 - l_fJump_2);
    float l_strengthR =
        -i_scaling * l_detInv * (l_fJump_2 - l_waveSpeedL * l_fJump_1);

    if(i_b[idx] >= 0){
    l_strengthL = 0;
    l_strengthR = 0;
    }
    if(i_b[idx + i_stride] >= 0){
    l_strengthL = 0;
    l_strengthR = 0;
    }
    __syncthreads();

    if (l_waveSpeedL < 0) {
      atomicAdd(&i_height_new[idx], l_strengthL);
      atomicAdd(&i_momentum_new[idx], l_strengthL * l_waveSpeedL);
    } else {
      atomicAdd(&i_height_new[idx + i_stride], l_strengthL);
      atomicAdd(&i_momentum_new[idx + i_stride], l_strengthL * l_waveSpeedL);
    }

    if (l_waveSpeedR > 0) {
      atomicAdd(&i_height_new[idx + i_stride], l_strengthR);
      atomicAdd(&i_momentum_new[idx + i_stride], l_strengthR * l_waveSpeedR);
    } else {
      atomicAdd(&i_height_new[idx], l_strengthR);
      atomicAdd(&i_momentum_new[idx], l_strengthR * l_waveSpeedR);
    }
    __syncthreads();
  }
  i_height_old[idx] = i_height_new[idx];
  i_momentum_old[idx] = i_momentum_new[idx];
}
