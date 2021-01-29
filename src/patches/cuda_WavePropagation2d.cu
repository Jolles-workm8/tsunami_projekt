#include "cuda_WavePropagation2d.h"
#include "cuda_functions.cuh"
//#include "cuda_solver.cuh"
#include <iostream>
#include <cuda.h>
#include <cooperative_groups.h>
using namespace cooperative_groups;


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
  cudaMalloc((void **)&h_dev, size * sizeof(float));
  cudaMalloc((void **)&h_dev_UpdateL, size * sizeof(float));
  cudaMalloc((void **)&h_dev_UpdateR, size * sizeof(float));
  cudaMalloc((void **)&hu_dev, size * sizeof(float));
  cudaMalloc((void **)&hu_dev_UpdateL, size * sizeof(float));
  cudaMalloc((void **)&hu_dev_UpdateR, size * sizeof(float));
  cudaMalloc((void **)&hv_dev, size * sizeof(float));
  cudaMalloc((void **)&b_dev, size * sizeof(float));
}

tsunami_lab::patches::cuda_WavePropagation2d::~cuda_WavePropagation2d() {
  // free memory on CPU
  free(m_h);
  free(m_hu);
  free(m_hv);
  free(m_b);

  // free memory on GPU
  cudaFree(h_dev);
  cudaFree(h_dev_UpdateR);
  cudaFree(h_dev_UpdateL);
  cudaFree(hu_dev);
  cudaFree(hu_dev_UpdateR);
  cudaFree(hu_dev_UpdateL);
  cudaFree(hv_dev);
  cudaFree(b_dev);
}
void tsunami_lab::patches::cuda_WavePropagation2d::MemTransfer() {
  cudaMemcpy(h_dev, m_h, size * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(hu_dev, m_hu, size * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(hv_dev, m_hv, size * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(b_dev, m_b, size * sizeof(float), cudaMemcpyHostToDevice);
}
void tsunami_lab::patches::cuda_WavePropagation2d::timeStep(
    t_real i_scaling, t_idx i_computeSteps) {
     printf("asfgh");
  t_idx l_nx = m_xCells + 2;
  t_idx l_ny = m_yCells + 2;

  /*
  dim3 threadsPerBlock(1,1);
  dim3 BlocksPerGrid(l_nx+1, l_ny+1);

  solverInit<<<BlocksPerGrid, threadsPerBlock>>>(
      h_dev, h_dev_UpdateR, h_dev_UpdateL, hu_dev, hu_dev_UpdateR, hu_dev_UpdateL, hv_dev,
      l_nx, l_ny, b_dev, i_computeSteps, i_scaling);
*/
int supportsCoopLaunch = 0;
cudaDeviceGetAttribute(&supportsCoopLaunch, cudaDevAttrCooperativeLaunch, 0);
printf("test %d \n", supportsCoopLaunch);
if(supportsCoopLaunch){

  
  /// This will launch a grid that can maximally fill the GPU, on the default stream with kernel arguments
  int numBlocksPerSm = 0;
  // Number of threads my_kernel will be launched with
  int numThreads = 256;
  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, 0);
  cudaOccupancyMaxActiveBlocksPerMultiprocessor(&numBlocksPerSm, solverInit, numThreads, 0);
  // launch
  void *kernelArgs[] = {h_dev, h_dev_UpdateR, h_dev_UpdateL, hu_dev, hu_dev_UpdateR, hu_dev_UpdateL, hv_dev,
      &l_nx, &l_ny, b_dev, &i_computeSteps, &i_scaling };
  dim3 dimBlock(16, 16, 1);
  dim3 dimGrid(deviceProp.multiProcessorCount*numBlocksPerSm, 1, 1);
  cudaLaunchCooperativeKernel((void*)solverInit, dimGrid, dimBlock, kernelArgs);

  cudaDeviceSynchronize();

  cudaMemcpy(m_h, h_dev, size * sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(m_hu, hu_dev, size * sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(m_hv, hv_dev, size * sizeof(float), cudaMemcpyDeviceToHost);
  
  } 
}

__global__ void solverInit(float *i_h, float *l_h_UpdateR, float *l_h_UpdateL, float *i_hu,
                           float *l_hu_UpdateR, float *l_hu_UpdateL, float *i_hv,
                           int i_nx, int i_ny, float *i_b, int i_iter,
                           float i_scaling) {
  
  printf("test \n");
  grid_group g = this_grid();
  printf("test %d\n", g.is_valid());
                            
  for (int i = 0; i < i_iter; ++i) {
    initUpdatesZero(l_h_UpdateR, l_h_UpdateR, l_hu_UpdateR, l_hu_UpdateL, i_nx, i_ny);
    setGhostOutflow(i_h, i_hu, i_hv, i_b, i_nx, i_ny);
    printf("test %d\n", g.is_valid());
    g.sync();
     printf("test123\n");

    //x Sweep
    
    netUpdates(i_h, l_h_UpdateR, l_h_UpdateR, i_hu, l_hu_UpdateR, l_hu_UpdateL, i_nx - 1, i_ny, i_b,
               i_scaling, i_nx, 1);
    g.sync();
    updateValues(i_h, l_h_UpdateR, l_h_UpdateR, i_hu, l_hu_UpdateR, l_hu_UpdateL, i_nx, i_ny);
    initUpdatesZero(l_h_UpdateR, l_h_UpdateR, l_hu_UpdateR, l_hu_UpdateL, i_nx, i_ny);
    g.sync();
  
    //y Sweep
    netUpdates(i_h, l_h_UpdateR, l_h_UpdateR, i_hv, l_hu_UpdateR, l_hu_UpdateL, i_nx, i_ny - 1, i_b,
               i_scaling, i_nx, i_nx);
    g.sync();
    updateValues(i_h, l_h_UpdateR, l_h_UpdateR, i_hv, l_hu_UpdateR, l_hu_UpdateL, i_nx, i_ny);
    g.sync();
  }
}

__device__ void updateValues(float *o_h, float *i_h_UpdateR, float *i_h_UpdateL, float *o_hu,
                           float *i_hu_UpdateR, float *i_hu_UpdateL, int i_nx, int i_ny){
  
  
  int l_i = blockIdx.x * blockDim.x + threadIdx.x;
  int l_j = blockIdx.y * blockDim.y + threadIdx.y;
  int idx = l_i + l_j * i_nx;

  if (l_i < i_nx && l_j < i_ny){
  o_h[idx] =o_h[idx] + i_h_UpdateR[idx] + i_h_UpdateL[idx];
  o_hu[idx] =o_hu[idx] + i_hu_UpdateR[idx] + i_hu_UpdateL[idx];      
  }                   
}

__device__ void initUpdatesZero(float *i_h_UpdateR, float *i_h_UpdateL, float *i_hu_UpdateR, float *i_hu_UpdateL, int i_nx, int i_ny){
  int l_i = blockIdx.x * blockDim.x + threadIdx.x;
  int l_j = blockIdx.y * blockDim.y + threadIdx.y;
  int idx = l_i + l_j * i_nx;

  if (l_i < i_nx && l_j < i_ny){
    i_h_UpdateR[idx] = 0;
    i_h_UpdateL[idx] = 0;
    i_hu_UpdateR[idx] = 0;
    i_hu_UpdateL[idx] = 0;
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

__device__ void netUpdates(float *i_height, float *o_height_UpdateR, float *o_height_UpdateL,
                           float *i_momentum, float *o_momentum_UpdateR, float *o_momentum_UpdateL,
                           int i_xEdges, int i_yEdges, float *i_b, float i_scaling,
                           int i_nx, int i_stride) {
  
  printf("Edges %d and %d\n",i_xEdges, i_yEdges);
  int l_i = blockIdx.x * blockDim.x + threadIdx.x;
  int l_j = blockIdx.y * blockDim.y + threadIdx.y;
  int idx = l_i + l_j * i_nx;

  float m_gSqrt = sqrtf(9.812);
  float m_g = 9.812;

  if (l_i < i_xEdges && l_j < i_yEdges) {
    printf("Hello from block %d, thread %d with access to %d and %d and %d and %d \n", l_i, l_j, idx, idx+i_stride, i_xEdges, i_yEdges);
    // compute u for left and right
    float l_uL = i_momentum[idx] / i_height[idx];
    float l_uR = i_momentum[idx + i_stride] / i_height[idx + i_stride];

    float l_hL = i_height[idx];
    float l_hR = i_height[idx + i_stride];

    float l_huL = i_momentum[idx];
    float l_huR = i_momentum[idx + i_stride];

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


    if (l_waveSpeedL < 0) {
      o_height_UpdateL[idx] = l_strengthL;
      o_momentum_UpdateL[idx] = l_strengthL * l_waveSpeedL;
    } else {
      o_height_UpdateR[idx + i_stride] = l_strengthL;
      o_momentum_UpdateR[idx + i_stride] = l_strengthL * l_waveSpeedL;
    }

    if (l_waveSpeedR > 0) {
      o_height_UpdateR[idx + i_stride] = l_strengthR;
      o_momentum_UpdateR[idx + i_stride] = l_strengthR * l_waveSpeedR;

    } else {
      o_height_UpdateL[idx] = l_strengthR;
      o_momentum_UpdateL[idx] = l_strengthR * l_waveSpeedR;
    }

    if(i_b[idx] >= 0){
    o_height_UpdateL[idx] = 0;
    o_momentum_UpdateL[idx] = 0;
    }
    if(i_b[idx + i_stride] >= 0){
    o_height_UpdateR[idx + i_stride] = 0;
    o_momentum_UpdateR[idx + i_stride] = 0;
    }
  }
}
