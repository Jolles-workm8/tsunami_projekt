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

#include "../solvers/Roe.h"
#include "../solvers/fwave.h"

tsunami_lab::patches::WavePropagation2d::WavePropagation2d(t_idx i_xCells,
                                                           t_idx i_yCells) {
  m_xCells = i_xCells;
  m_yCells = i_yCells;

  // allocate memory including ghostcells on each side.
  for (unsigned short l_st = 0; l_st < 2; l_st++) {
    m_h[l_st] = new t_real[(m_xCells + 2) * (m_yCells + 2)];
    m_hu[l_st] = new t_real[(m_xCells + 2) * (m_yCells + 2)];
    m_hv[l_st] = new t_real[(m_xCells + 2) * (m_yCells + 2)];
  }
  m_b = new t_real[(m_xCells + 2) * (m_yCells + 2)];

  // init to zero
  for (unsigned short l_ce = 0; l_ce != (m_xCells + 2) * (m_yCells + 2);
       l_ce++) {
    m_b[l_ce] = 0;
    for (unsigned short l_st = 0; l_st < 2; l_st++) {
      m_h[l_st][l_ce] = 0;
      m_hu[l_st][l_ce] = 0;
      m_hv[l_st][l_ce] = 0;
    }
  }
}

tsunami_lab::patches::WavePropagation2d::~WavePropagation2d() {
  delete[] m_b;

  for (unsigned short l_st = 0; l_st < 2; l_st++) {
    delete[] m_h[l_st];
    delete[] m_hu[l_st];
    delete[] m_hv[l_st];
  }
}

void tsunami_lab::patches::WavePropagation2d::timeStep(t_real i_scaling,
                                                       int solver) {
  // pointers to old and new data
  t_real *l_hOld = m_h[m_step];
  t_real *l_huOld = m_hu[m_step];
  t_real *l_hvOld = m_hv[m_step];

  m_step = (m_step + 1) % 2;
  t_real *l_hNew = m_h[m_step];
  t_real *l_huNew = m_hu[m_step];
  t_real *l_hvNew = m_hv[m_step];

  // init new cell quantities
  for (t_idx l_ce = 0; l_ce < (m_xCells + 2) * (m_yCells + 2); l_ce++) {
    l_hNew[l_ce] = l_hOld[l_ce];
    l_huNew[l_ce] = l_huOld[l_ce];
    l_hvNew[l_ce] = l_hvOld[l_ce];
  }

  if (solver == 0) {
    // iterate over all collums in x direction with ghost cells
    for (t_idx l_ceY = 0; l_ceY < (m_yCells + 2); l_ceY++) {
      // iterate over edges in x direction and update with Riemann solutions
      for (t_idx l_ceX = 0; l_ceX < (m_xCells + 1); l_ceX++) {
        // determine left and right cell-id
        t_idx l_ceL = calculateArrayPosition(l_ceX, l_ceY);
        t_idx l_ceR = calculateArrayPosition(l_ceX + 1, l_ceY);

        // compute net-updates
        t_real l_netUpdates[2][2];

        solvers::Roe::netUpdates(l_hNew[l_ceL], l_hNew[l_ceR], l_huOld[l_ceL],
                                 l_huOld[l_ceR], l_netUpdates[0],
                                 l_netUpdates[1]);

        // update the cells' quantities
        l_hNew[l_ceL] -= i_scaling * l_netUpdates[0][0];
        l_huNew[l_ceL] -= i_scaling * l_netUpdates[0][1];

        l_hNew[l_ceR] -= i_scaling * l_netUpdates[1][0];
        l_huNew[l_ceR] -= i_scaling * l_netUpdates[1][1];
      }
    }

    // iterate over all rows in y direction without ghost cells
    for (t_idx l_ceX = 1; l_ceX < (m_xCells + 1); l_ceX++) {
      // iterate over edges in y direction and update with Riemann solutions
      for (t_idx l_ceY = 0; l_ceY < (m_yCells + 1); l_ceY++) {
        // determine left and right cell-id
        t_idx l_ceB = calculateArrayPosition(l_ceX, l_ceY);
        t_idx l_ceT = calculateArrayPosition(l_ceX, l_ceY + 1);

        // compute net-updates
        t_real l_netUpdates[2][2];

        solvers::Roe::netUpdates(l_hNew[l_ceB], l_hNew[l_ceT], l_hvOld[l_ceB],
                                 l_hvOld[l_ceT], l_netUpdates[0],
                                 l_netUpdates[1]);

        // update the cells' quantities
        l_hNew[l_ceB] -= i_scaling * l_netUpdates[0][0];
        l_hvNew[l_ceB] -= i_scaling * l_netUpdates[0][1];

        l_hNew[l_ceT] -= i_scaling * l_netUpdates[1][0];
        l_hvNew[l_ceT] -= i_scaling * l_netUpdates[1][1];
      }
    }
  }
}

void tsunami_lab::patches::WavePropagation2d::setGhostOutflow() {
  t_real *l_h = m_h[m_step];
  t_real *l_hu = m_hu[m_step];
  t_real *l_hv = m_hv[m_step];
  t_idx l_displacementFrom;
  t_idx l_displacementTo;

  // set ghost outflow for all outer cells
  // displacement = (m_xCells+2)*(m_yCells + 1);
  for (unsigned short l_ce = 1; l_ce <= (m_xCells); l_ce++) {
    l_displacementFrom = calculateArrayPosition(l_ce, 1);
    l_displacementTo = calculateArrayPosition(l_ce, 0);
    l_h[l_displacementTo] = l_h[l_displacementFrom];
    l_hu[l_displacementTo] = l_hu[l_displacementFrom];
    l_hv[l_displacementTo] = l_hv[l_displacementFrom];

    l_displacementFrom = calculateArrayPosition(l_ce, m_yCells);
    l_displacementTo = calculateArrayPosition(l_ce, m_yCells + 1);
    l_h[l_displacementTo] = l_h[l_displacementFrom];
    l_hu[l_displacementTo] = l_hu[l_displacementFrom];
    l_hv[l_displacementTo] = l_hv[l_displacementFrom];
  }

  for (unsigned short l_ce = 0; l_ce <= (m_yCells + 1); l_ce++) {
    l_displacementFrom = calculateArrayPosition(1, l_ce);
    l_displacementTo = calculateArrayPosition(0, l_ce);
    l_h[l_displacementTo] = l_h[l_displacementFrom];
    l_hu[l_displacementTo] = l_hu[l_displacementFrom];
    l_hv[l_displacementTo] = l_hv[l_displacementFrom];

    l_displacementFrom = calculateArrayPosition(m_xCells, l_ce);
    l_displacementTo = calculateArrayPosition(m_xCells + 1, l_ce);
    l_h[l_displacementTo] = l_h[l_displacementFrom];
    l_hu[l_displacementTo] = l_hu[l_displacementFrom];
    l_hv[l_displacementTo] = l_hv[l_displacementFrom];
  }
}
