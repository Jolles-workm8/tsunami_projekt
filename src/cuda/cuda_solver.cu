#include <cuda.h>

#include <iostream>

__device__ void netUpdates(float *i_height_old, float *i_height_new,
                           float *i_momentum_old, float *i_momentum_new,
                           int i_nx, int i_ny, float *i_b, float i_scaling) {
  int scaling = 4;
  int l_i = blockIdx.x * blockDim.x + threadIdx.x;
  int l_j = blockIdx.y * blockDim.y + threadIdx.y;
  int idx = l_i * scaling + l_j;

  int i_displ = 1;
  float m_gSqrt = sqrtf(9.812);
  float m_g = 9.812;

  if (l_i < i_nx - 1 && l_j < i_ny - 1) {
    // compute u for left and right
    float l_uL = i_momentum_old[idx] / i_height_old[idx];
    float l_uR = i_momentum_old[idx + i_displ] / i_height_old[idx + i_displ];

    float l_hL = i_height_old[idx];
    float l_hR = i_height_old[idx + i_displ];

    float l_huL = i_momentum_old[idx];
    float l_huR = i_momentum_old[idx + i_displ];

    float l_bL = i_b[idx];
    float l_bR = i_b[idx + i_displ];

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

    __syncthreads();

    if (l_waveSpeedL < 0) {
      atomicAdd(&i_height_new[idx], l_strengthL);
      atomicAdd(&i_momentum_new[idx], l_strengthL * l_waveSpeedL);
    } else {
      atomicAdd(&i_height_new[idx + i_displ], l_strengthL);
      atomicAdd(&i_momentum_new[idx + i_displ], l_strengthL * l_waveSpeedL);
    }

    if (l_waveSpeedR > 0) {
      atomicAdd(&i_height_new[idx + i_displ], l_strengthR);
      atomicAdd(&i_momentum_new[idx + i_displ], l_strengthR * l_waveSpeedR);
    } else {
      atomicAdd(&i_height_new[idx], l_strengthR);
      atomicAdd(&i_momentum_new[idx], l_strengthR * l_waveSpeedR);
    }
    __syncthreads();

    i_height_old[idx] = i_height_new[idx];
    i_momentum_old[idx] = i_momentum_new[idx];
  }
}

__global__ void WavePropagation2d(float *i_h_old, float *i_h_new,
                                  float *i_hu_old, float *i_hu_new,
                                  float *i_hv_new, float *i_hv_old, int i_nx,
                                  int i_ny, float *i_b, int i_iter,
                                  float i_scaling) {
  for (int i = 0; i < i_iter; ++i) {
    netUpdates(i_h_old, i_h_new, i_hu_old, i_hu_new, i_nx, i_ny, i_b,
               i_scaling);
    __syncthreads();
  }
}

int main() {
  int N = 128;
  int iteration = 32;
  int size = 16;
  int i_nx = 4;
  int i_ny = 4;
  float scaling = 0.5;
  float *h_host, *hu_host, *hv_host, *h_dev_new, *h_dev_old, *hu_dev_new,
      *hu_dev_old, *hv_dev_new, *hv_dev_old, *b_host, *b_dev;
  h_host = (float *)malloc(size * sizeof(float));
  hu_host = (float *)malloc(size * sizeof(float));
  hu_host = (float *)malloc(size * sizeof(float));
  b_host = (float *)malloc(size * sizeof(float));
  cudaMalloc((void **)&h_dev_old, size * sizeof(float));
  cudaMalloc((void **)&h_dev_new, size * sizeof(float));
  cudaMalloc((void **)&hu_dev_old, size * sizeof(float));
  cudaMalloc((void **)&hu_dev_new, size * sizeof(float));
  cudaMalloc((void **)&hv_dev_old, size * sizeof(float));
  cudaMalloc((void **)&hv_dev_new, size * sizeof(float));
  cudaMalloc((void **)&b_dev, size * sizeof(float));

  for (int i = 0; i < size; i++) {
    h_host[i] = i;
    hu_host[i] = i;
    hv_host[i] = i;
    b_host[i] = 0;
  }

  cudaMemcpy(h_dev_old, h_host, size * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(hu_dev_old, hu_host, size * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(hv_dev_old, hv_host, size * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(b_dev, b_host, size * sizeof(float), cudaMemcpyHostToDevice);

  cudaError_t err = cudaGetLastError();
  dim3 threadsPerBlock(i_nx, i_ny);
  for (int i = 0; i < N / iteration; i++) {
    WavePropagation2d<<<1, threadsPerBlock>>>(
        h_dev_old, h_dev_new, hu_dev_old, hu_dev_new, hv_dev_new, hv_dev_old,
        i_nx, i_ny, b_dev, iteration, scaling);

    cudaDeviceSynchronize();
  }
  cudaMemcpy(h_host, h_dev_old, size * sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(hu_host, hu_dev_old, size * sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(hv_host, hv_dev_old, size * sizeof(float), cudaMemcpyDeviceToHost);

  //
  // for (int i = 0; i < size; i++) {
  //   std::cout << h_host[i] << '\n';
  // }

  free(h_host);
  free(hu_host);
  free(b_host);

  cudaFree(h_dev_new);
  cudaFree(h_dev_old);
  cudaFree(hu_dev_new);
  cudaFree(hu_dev_old);
  cudaFree(b_dev);

  return 0;
}
