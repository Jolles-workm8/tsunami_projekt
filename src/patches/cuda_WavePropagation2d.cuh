#ifndef __GPU_CUDA_WAVEPROPAGATION_2D
#define __GPU_CUDA_WAVEPROPAGATION_2D


#include "WavePropagation.h"

namespace tsunami_lab {
    namespace patches {
    class cuda_WavePropagation2d;
    }
}

class tsunami_lab::patches::cuda_WavePropagation2d : public WavePropagation{
private:
    t_idx m_xCells;

    t_idx m_yCells;

    t_real *m_h;
    
    t_real *m_hv;

    t_real *m_hu;

    t_real *m_b;

    t_real *h_dev_new, *h_dev_old, *hu_dev_new,
      *hu_dev_old, *hv_dev_new, *hv_dev_old, *b_dev;
public:

    cuda_WavePropagation2d(t_idx i_xCells, t_idx i_yCells);

    ~cuda_WavePropagation2d();
  
    void timeStep(t_real i_scaling, t_idx i_computeSteps);

    /**
    * Gets the stride in y-direction. x-direction is stride-1.
    *
    * @return stride in y-direction.
    **/
    t_idx getStride() { return (m_xCells + 2); }

    /**
    * Gets cells' water heights.
    *
    * @return water heights.
    */
    t_real const *getHeight() { return m_h + 3 + m_xCells; }

    /**
    * Gets the cells' momenta in x-direction.
    *
    * @return momenta in x-direction.
    **/
     t_real const *getMomentumX() { return m_hu + 3 + m_xCells; }
 
    /**
    * Gets the cell's momentum in y-direction.
    *
    * @return momenta in y-direction.
    **/
    t_real const *getMomentumY() { return m_hv + 3 + m_xCells; }

    /**
    * Gets the cells' bathymetry in x-direction.
    *
    * @return bathymetry in x-direction.
    **/
    t_real const *getBathymetry() { return m_b + 3 + m_xCells; }

    /**
    * Sets the height of the cell to the given value.
    *
    * @param i_ix id of the cell in x-direction.
    * @param i_y id of the cell in y-direction.
    * @param i_h water height.
    **/
    void setHeight(t_idx i_ix, t_idx i_iy, t_real i_h) {
        m_h[(i_ix + 1) + ((i_iy + 1) * (m_xCells + 2))] = i_h;
    }

    /**
     * Sets the momentum in x-direction to the given value.
    *
     * @param i_ix id of the cell in x-direction.
    * @param i_hu momentum in x-direction.
    * @param i_y id of the cell in y-direction.
    **/
    void setMomentumX(t_idx i_ix, t_idx i_iy, t_real i_hu) {
        m_hu[(i_ix + 1) + ((i_iy + 1) * (m_xCells + 2))] = i_hu;
    }

    /**
    * Sets the momentum in y-direction to the given value.
    *
    * @param i_y id of the cell in y-direction.
    * @param i_ix id of the cell in x-direction.
    * @param i_hu momentum in x-direction.
    **/
    void setMomentumY(t_idx i_ix, t_idx i_iy, t_real i_hv) {
        m_hu[(i_ix + 1) + ((i_iy + 1) * (m_xCells + 2))] = i_hv;
    };

    /**
    * Sets the bathymetry value to the cell given it's id.
    *
    * @param i_ix id of the cell in x-direction.
    * @param i_iy id of the cell in y-direction.
    * @param i_b bathymetry value of the cell.
    **/
    void setBathymetry(t_idx i_ix, t_idx i_iy, t_real i_b) {
        m_b[(i_ix + 1) + ((i_iy + 1) * (m_xCells + 2))] = i_b;
    }

    
};

__global__ void solverInit(float *i_h_old, float *i_h_new,
                                  float *i_hu_old, float *i_hu_new,
                                  float *i_hv_old, float *i_hv_new, int i_nx,
                                  int i_ny, float *i_b, int i_iter,
                                  float i_scaling);

__device__ void netUpdates(float *i_height_old, float *i_height_new,
                           float *i_momentum_old, float *i_momentum_new,
                           int i_nx, int i_ny, float *i_b, float i_scaling,
                           int idx, int i_stride);

__device__ void setGhostOutflow(float *i_height, float *i_hu, float *i_hv,
                                float *i_b, int i_nx, int i_ny);
#endif