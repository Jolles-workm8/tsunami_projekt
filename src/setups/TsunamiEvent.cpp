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
#include "TsunamiEvent.h"

#include <cmath>

tsunami_lab::setups::ArtificialTsunami::ArtificialTsunami() {}

tsunami_lab::t_real tsunami_lab::setups::ArtificialTsunami::getHeight(
    t_real, t_real) const {
  return 100;
}

tsunami_lab::t_real tsunami_lab::setups::ArtificialTsunami::getMomentumX(
    t_real, t_real) const {
  return 0;
}

tsunami_lab::t_real tsunami_lab::setups::ArtificialTsunami::getMomentumY(
    t_real, t_real) const {
  return 0;
}

tsunami_lab::t_real tsunami_lab::setups::ArtificialTsunami::getBathymetry(
    t_real i_x, t_real i_y) const {
  i_x -= 5000;
  i_y -= 5000;
  float pi = 3.14159;
  if ((-500 <= i_x) && (i_x <= 500) && (-500 <= i_y) && (i_y <= 500)) {
    return 5 * (std::sin(((float)0.002 * i_x + 1) * pi) *
                (-1 * (0.002 * i_y) * (0.002 * i_y) + 1));
  } else {
    return 0;
  }
}
