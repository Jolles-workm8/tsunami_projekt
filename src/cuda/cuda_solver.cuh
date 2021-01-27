#ifndef __GPU_CUDA_SOLVER_H__
#define __GPU_CUDA_SOLVER_H__

#include "cuda.h"

namespace tsunami_lab {
namespace cuda {
class cuda_solver;
}
}  // namespace tsunami_lab

class tsunami_lab::cuda::cuda_solver : public cuda_solver {
 private:
 public:
  __global__ void CudaWaveprop(t_idx i_xCells, t_idx i_yCells);

  __global__ ~CudaWaveprop();

  __global__ void timeStep(t_real *i_h_old, t_real *i_h_new, t_real *i_hu_old,
                           t_real *i_hu_new, t_real *i_hv_old, t_real *i_hv_new,
                           t_idx i_nx, t_idx i_ny, t_real *i_b, t_idx i_iter,
                           t_real i_scaling);

  __device__ void setGhostOutflow(t_real *i_height, t_real *i_hu, t_real *i_hv,
                                  t_real *i_b, t_idx i_nx, t_idx i_ny);

  __device__ void netUpdates(t_real *i_height_old, t_real *i_height_new,
                             t_real *i_momentum_old, t_real *i_momentum_new,
                             int i_nx, int i_ny, t_real *i_b, t_real i_scaling,
                             int idx, int i_stride);
}
#endif /* end of include guard: __GPU_CUDA_SOLVER_H__ */
