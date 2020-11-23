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
 * Fwave solver for the shallow water equations.
 **/
#include "fwave.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>

// compute the lambdas we need
void tsunami_lab::solvers::fwave::waveSpeeds(t_real i_hL, t_real i_hR,
                                             t_real i_uL, t_real i_uR,
                                             t_real &o_waveSpeedL,
                                             t_real &o_waveSpeedR) {
  // pre-compute square-root ops
  t_real l_hSqrtL = std::sqrt(i_hL);
  t_real l_hSqrtR = std::sqrt(i_hR);

  // compute Roe averages
  t_real l_hRoe = 0.5f * (i_hL + i_hR);
  t_real l_uRoe = l_hSqrtL * i_uL + l_hSqrtR * i_uR;
  l_uRoe /= l_hSqrtL + l_hSqrtR;

  // compute wave speeds
  t_real l_ghSqrtRoe = m_gSqrt * std::sqrt(l_hRoe);
  o_waveSpeedL = l_uRoe - l_ghSqrtRoe;
  o_waveSpeedR = l_uRoe + l_ghSqrtRoe;
}

// compute the inverse of the matrix
void tsunami_lab::solvers::fwave::waveStrengths(
    t_real i_hL, t_real i_hR, t_real i_huL, t_real i_huR, t_real i_waveSpeedL,
    t_real i_waveSpeedR, t_real i_bL, t_real i_bR, t_real &o_strengthL,
    t_real &o_strengthR) {
  // compute inverse of right eigenvector-matrix
  t_real l_detInv = 1 / (i_waveSpeedR - i_waveSpeedL);

  // compute the bathymetry effect
  t_real l_bathEff = -m_g * (i_bR - i_bL) * (i_hL + i_hR) / 2;

  // compute jump in the flux
  t_real l_fJump_1 = i_huR - i_huL;
  t_real l_fJump_2 = i_huR * i_huR / i_hR - i_huL * i_huL / i_hL +
                     (m_g / 2) * (i_hR * i_hR - i_hL * i_hL);
  l_fJump_2 -= l_bathEff;

  // compute the alpha values
  o_strengthL = l_detInv * (i_waveSpeedR * l_fJump_1 - l_fJump_2);
  o_strengthR = l_detInv * (l_fJump_2 - i_waveSpeedL * l_fJump_1);
}

void tsunami_lab::solvers::fwave::netUpdates(t_real i_hL, t_real i_hR,
                                             t_real i_huL, t_real i_huR,
                                             t_real i_bL, t_real i_bR,
                                             t_real o_netUpdateL[2],
                                             t_real o_netUpdateR[2],
                                             t_real &o_speed) {
  // compute particle velocities, redundant
  t_real l_uL = i_huL / i_hL;
  t_real l_uR = i_huR / i_hR;

  // compute wave speeds
  t_real l_sL = 0;
  t_real l_sR = 0;

  waveSpeeds(i_hL, i_hR, l_uL, l_uR, l_sL, l_sR);

  // compute wave strengths
  t_real l_aL = 0;
  t_real l_aR = 0;

  waveStrengths(i_hL, i_hR, i_huL, i_huR, l_sL, l_sR, i_bL, i_bR, l_aL, l_aR);

  // compute scaled waves
  t_real l_waveL[2] = {0};
  t_real l_waveR[2] = {0};

  l_waveL[0] = l_aL;
  l_waveL[1] = l_sL * l_aL;

  l_waveR[0] = l_aR;
  l_waveR[1] = l_sR * l_aR;

  // compute the max wavespeed of the 2 Volumes
  o_speed = std::max(std::abs(l_sL), std::abs(l_sR));

  // set net-updates depending on wave speeds
  for (unsigned short l_qt = 0; l_qt < 2; l_qt++) {
    // init
    o_netUpdateL[l_qt] = 0;
    o_netUpdateR[l_qt] = 0;

    // 1st wave
    if (l_sL < 0) {
      o_netUpdateL[l_qt] += l_waveL[l_qt];
    } else {
      o_netUpdateR[l_qt] += l_waveL[l_qt];
    }

    // 2nd wave
    if (l_sR > 0) {
      o_netUpdateR[l_qt] += l_waveR[l_qt];
    } else {
      o_netUpdateL[l_qt] += l_waveR[l_qt];
    }
  }
}
