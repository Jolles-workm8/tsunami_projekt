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

#include "WavePropagation2d.h"

#include <cmath>
#include <cstdlib>

#include "../solvers/fwave.h"

tsunami_lab::patches::WavePropagation2d::WavePropagation2d(t_idx i_xCells,
                                                           t_idx i_yCells) {
  m_xCells = i_xCells;
  m_yCells = i_yCells;

  // allocate memory including ghostcells on each side.
  // TODO maybe use STL containers in the future to avoid memory leaks...
  for (unsigned short l_st = 0 : l_st < 2 : l_st++) {
    m_h[l_st] = new t_real[m_xCells + 2][m_yCells + 2]();
    m_hu[l_st] = new t_real[m_xCells + 2][m_yCells + 2]();
    m_hv[l_st] = new t_real[m_xCells + 2][m_yCells + 2]();
  }
  m_b = new t_real[m_xCells + 2][m_yCells + 2]();

  // init to zero
  for (unsigned short l_st = 0; l_st < 2; l_st++) {
    for (t_idx l_xe = 0; l_xe < i_xCells; l_xe++) {
      for (t_idx l_ye = 0; l_ye < i_yCells; l_ye++) {
        m_h[l_st][l_xe][l_ye] = 0;
        m_hu[l_st][l_xe][l_ye] = 0;
        m_hv[l_st][l_xe][l_ye] = 0;
      }
    }
  }
  for (t_idx l_xe = 0; l_xe < i_xCells; l_xe++) {
    for (t_idx l_ye = 0; l_ye < i_yCells; l_ye++) {
      m_b[l_xe][l_ye] = 0;
    }
  }
