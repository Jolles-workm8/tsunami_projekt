/**
 * @author Julius Isken, Max Engel
 *
 * @section LICENSE
 * Copyright 2020, Julius Isken, Max Engel
 *
 * Redistribution and use in source and binary forms, with or without
 *modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *this list of conditions and the following disclaimer in the documentation
 *and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 *FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 *DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 *SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 *OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Two-dimensional wave propagation patch.
 **/
#ifndef TSUNAMI_LAB_PATCHES_WAVE_PROPAGATION_2D
#define TSUNAMI_LAB_PATCHES_WAVE_PROPAGATION_2D

#include "WavePropagation.h"

namespace tsunami_lab {
namespace patches {
class WavePropagation2d;
}
}  // namespace tsunami_lab

class tsunami_lab::patches::WavePropagation2d : public WavePropagation {
 private:
  //! current step which indicates the active values in the arrays below
  unsigned short m_step = 0;

  //! number of cells discretizing the computational domain in x-direction
  t_idx m_xCells = 0;

  //! number of cells discretizing the computational domain in y-direction
  t_idx m_yCells = 0;

  //! water heights for the current and next time step for all cells
  t_real *m_h[2] = {nullptr, nullptr};

  //! momenta for the current and the next time step for all cells in
  //! x-direction
  t_real *m_hu[2] = {nullptr, nullptr};

  //! momenta for the current and the next time step for all cells in
  //! y-direction
  t_real *m_hv[2] = {nullptr, nullptr};

  //! bathymetry data for all cells
  t_real *m_b = nullptr;

  //!  is left boundary reflecting
  bool m_reflBoundL = false;

  //!  is right boundary reflecting
  bool m_reflBoundR = false;

 public:
  /**
   * Constructs the 2d wave propagation solver.
   *
   * @param i_xCells number of cells in x-direction.
   * @param i_yCells number of cells in y-direction.
   **/
  WavePropagation2d(t_idx i_xCells, t_idx i_yCells);

  /**
   * Destructor which frees all allocated memory.
   **/
  ~WavePropagation2d();

  /**
   * Performs a time step.
   *
   * @param i_scaling scaling of the time step (dt / dx).
   **/
  void timeStep(t_real i_scaling);

  /**
   * Sets the values of the ghost cells according to outflow boundary
   *conditions.
   **/
  void setGhostOutflow();

  /**
   * calculate the position in one dimensional array
   *
   * @return the position.
   **/
  t_idx calculateArrayPosition(t_idx i_ix, t_idx i_iy) {
    return i_ix + (i_iy * (m_xCells + 2));
  }

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
  t_real const *getHeight() { return m_h[m_step] + 3 + m_xCells; }

  /**
   * Gets the cells' momenta in x-direction.
   *
   * @return momenta in x-direction.
   **/
  t_real const *getMomentumX() { return m_hu[m_step] + 3 + m_xCells; }

  /**
   * Gets the cell's momentum in y-direction.
   *
   * @return momenta in y-direction.
   **/
  t_real const *getMomentumY() { return m_hv[m_step] + 3 + m_xCells; }
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
    m_h[m_step][(i_ix + 1) + ((i_iy + 1) * (m_xCells + 2))] = i_h;
  }

  /**
   * Sets the momentum in x-direction to the given value.
   *
   * @param i_ix id of the cell in x-direction.
   * @param i_hu momentum in x-direction.
   * @param i_y id of the cell in y-direction.
   **/
  void setMomentumX(t_idx i_ix, t_idx i_iy, t_real i_hu) {
    m_hu[m_step][(i_ix + 1) + ((i_iy + 1) * (m_xCells + 2))] = i_hu;
  }

  /**
   * Sets the momentum in y-direction to the given value.
   *
   * @param i_y id of the cell in y-direction.
   * @param i_ix id of the cell in x-direction.
   * @param i_hu momentum in x-direction.
   **/
  void setMomentumY(t_idx i_ix, t_idx i_iy, t_real i_hv) {
    m_hu[m_step][(i_ix + 1) + ((i_iy + 1) * (m_xCells + 2))] = i_hv;
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

  /**
   * Sets the ghost cells to reflecting or not.
   *
   * @param i_reflL reflection of the left ghost cell.
   * @param i_reflR reflection of the right ghost cell.
   **/

  // TODO : überlegen wie man Top Bottom left righ gut impleemntiert
  void setReflection(t_idx, bool i_reflL, bool i_reflR) {
    m_reflBoundL = i_reflL;
    m_reflBoundR = i_reflR;
  }
};

#endif
