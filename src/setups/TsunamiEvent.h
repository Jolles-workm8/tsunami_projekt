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
#ifndef TSUNAMI_EVENT_H_
#define TSUNAMI_EVENT_H_

#include "../io/NetCdf_Read.h"
#include "Setup.h"

namespace tsunami_lab {
namespace setups {
class TsunamiEvent;
}
}  // namespace tsunami_lab

/**
 * Tsunami Event
 **/

class tsunami_lab::setups::TsunamiEvent : public Setup {
 private:
  tsunami_lab::io::NetCdf_Read *l_netcdf;
  t_idx l_nx = 0;
  t_real const lambda = 50;


 public:
  /**
   *Constructor
   **/
  TsunamiEvent(t_idx i_nx, tsunami_lab::io::NetCdf_Read *i_netCdf);


  /**
   * Gets the water height at a given point.
   *
   * @param i_x x-coordinate of the queried point.
   * @param i_y y-coordinate of the queried point.
   * @return height at the given point.
   **/
  t_real getHeight(t_idx i_x, t_idx i_y) const;
  /**
   * Gets the momentum in x-direction.
   *
   * @return momentum in x-direction.
   **/
  t_real getMomentumX(t_idx, t_idx) const;

  /**
   * Gets the momentum in y-direction.
   *
   * @return momentum in y-direction.
   **/
  t_real getMomentumY(t_idx, t_idx) const;
  /**
   * Gets the bathymetry data.
   *
   * @param i_x x-coordinate of the queried point.
   * @param i_y y-coordinate of the queried point.
   *
   * @return the bathymetry data at the given point.
   **/
  t_real getBathymetry(t_idx i_x, t_idx i_y) const;
};
#endif
