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
 * one-dimensional shock shock problem
 **/
#include "ArtificialTsunami.h"

#include <cmath>

tsunami_lab::setups::ArtificialTsunami::ArtificialTsunami() {
  l_bath = 0;
  l_mom_x = 0;
  l_mom_y = 0;
  l_height = 200;
}

tsunami_lab::t_real tsunami_lab::setups::ArtificialTsunami::getHeight(
    t_idx, t_idx) const {
  return l_height;
}

tsunami_lab::t_real tsunami_lab::setups::ArtificialTsunami::getMomentumX(
    t_idx, t_idx) const {
  return l_mom_x;
}

tsunami_lab::t_real tsunami_lab::setups::ArtificialTsunami::getMomentumY(
    t_idx, t_idx) const {
  return l_mom_y;
}

tsunami_lab::t_real tsunami_lab::setups::ArtificialTsunami::getBathymetry(
    t_idx i_x, t_idx i_y) const {
  float pi = 3.14159;
  if ((4500 <= i_x) && (i_x <= 5500) && (4500 <= i_y) && (i_y <= 5500)) {
    return 5 * (std::sin(((float)0.002 * i_x + 1) * pi) *
                (-1 * (-0.002 * i_y) * (0.002 * i_y) + 1));
  } else {
    return l_bath;
  }
}
