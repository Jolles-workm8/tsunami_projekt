#include <cuda.h>
#include <thrust\device_vector.h>

#include <iostream>

__global__ void netUpdates(float *i_h_old, float *i_h_new, float *i_hu_old,
                           float *i_hu_new, int i_nx, int i_ny, float *i_b) {
  int scaling = 4;
  int l_i = blockIdx.x * blockDim.x + threadIdx.x;
  int l_j = blockIdx.y * blockDim.y + threadIdx.y;
  int tid = threadIdx.x + threadIdx.y * blockDim.x;
  int idx = l_i * scaling + l_j;

  int i_displ = 1;
  float m_gSqrt = sqrtf(9.812);
  float m_g = 9.812;

  i_h_new[idx] = 0;

  printf("%f\n %f", i_h_new[idx], i_h_old[idx]);
  printf("%d\n", idx);

  if (l_i < i_nx && l_j < i_ny) {
    printf("%d, %d\n", l_i, l_j);

    // compute u for left and right
    float l_uL = i_hu_old[idx] / i_h_old[idx];
    float l_uR = i_hu_old[idx + i_displ] / i_h_old[idx + i_displ];

    float l_hL = i_h_old[idx];
    float l_hR = i_h_old[idx + i_displ];

    float l_huL = i_hu_old[idx];
    float l_huR = i_hu_old[idx + i_displ];

    float l_bL = i_b[idx];
    float l_bR = i_b[idx + i_displ];

    __syncthreads();

    // compute WaveSpeed ,

    float l_hSqrtL = sqrtf(l_hL);
    float l_hSqrtR = sqrtf(l_hR);

    float l_hRoe = 0.5f * (l_hL + l_hR);
    float l_uRoe = l_hSqrtL * l_uL + l_hSqrtR * l_uR;
    l_uRoe /= l_hSqrtL + l_hSqrtR;

    float l_ghSqrtRoe = m_gSqrt * sqrtf(l_hRoe);

    __syncthreads();

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
    float l_strengthL = l_detInv * (l_waveSpeedR * l_fJump_1 - l_fJump_2);
    float l_strengthR = l_detInv * (l_fJump_2 - l_waveSpeedL * l_fJump_1);

    printf("Hello world  %f\n", l_waveSpeedL);
    __syncthreads();

    if (l_waveSpeedL < 0) {
      atomicAdd(&i_h_new[idx], l_strengthL);
      atomicAdd(&i_hu_new[idx], l_strengthL * l_waveSpeedL);
    } else {
      atomicAdd(&i_h_new[idx + i_displ], l_strengthL);
      atomicAdd(&i_hu_new[idx + i_displ], l_strengthL * l_waveSpeedL);
    }

    if (l_waveSpeedR > 0) {
      atomicAdd(&i_h_new[idx + i_displ], l_strengthR);
      atomicAdd(&i_hu_new[idx + i_displ], l_strengthR * l_waveSpeedR);
    } else {
      atomicAdd(&i_h_new[idx], l_strengthR);
      atomicAdd(&i_hu_new[idx], l_strengthR * l_waveSpeedR);
    }
    __syncthreads();

    printf("%f\n", i_h_new[idx]);

    swap(i_h_new, i_h_old);

    printf("%f\n", i_h_new[idx]);
  }
}

int main() {
  int size = 16;
  int i_nx = 4;
  int i_ny = 4;
  float *h_host, *hu_host, *h_dev_new, *h_dev_old, *hu_dev_new, *hu_dev_old,
      *b_host, *b_dev;
  h_host = (float *)malloc(size * sizeof(float));
  hu_host = (float *)malloc(size * sizeof(float));
  b_host = (float *)malloc(size * sizeof(float));
  cudaMalloc((void **)&h_dev_old, size * sizeof(float));
  cudaMalloc((void **)&h_dev_new, size * sizeof(float));
  cudaMalloc((void **)&hu_dev_old, size * sizeof(float));
  cudaMalloc((void **)&hu_dev_new, size * sizeof(float));
  cudaMalloc((void **)&b_dev, size * sizeof(float));

  for (int i = 0; i < size; i++) {
    h_host[i] = i;
    hu_host[i] = i;
    b_host[i] = 0;
  }

  cudaMemcpy(h_dev_old, h_host, size * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(hu_dev_old, hu_host, size * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(b_dev, b_host, size * sizeof(float), cudaMemcpyHostToDevice);

  cudaError_t err = cudaGetLastError();
  dim3 threadsPerBlock(4, 4);

  netUpdates<<<1, threadsPerBlock>>>(h_dev_old, h_dev_new, hu_dev_old,
                                     hu_dev_new, i_nx, i_ny, b_dev);

  cudaMemcpy(h_host, h_dev_old, size * sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(hu_host, hu_dev_old, size * sizeof(float), cudaMemcpyDeviceToHost);
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
