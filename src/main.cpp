/**
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 *
 * @section LICENSE
 * Copyright 2020, Friedrich Schiller University Jena
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
 * Entry-point for simulations.
 **/

#include <netcdf.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>

#include "io/Csv.h"
#include "io/NetCdf.h"
#include "patches/WavePropagation1d.h"
#include "patches/WavePropagation2d.h"
#include "setups/CircularDamBreak2d.h"
#include "setups/DamBreak1d.h"
#include "setups/DamBreakNew.h"
#include "setups/RareRare1d.h"
#include "setups/ShockShock1d.h"
#include "setups/SubcriticalFlow.h"
#include "setups/SupercriticalFlow.h"
#include "setups/TsunamiEvent.h"

int main(int i_argc, char *i_argv[]) {
  // number of cells in x- and y-direction. Default for y-dimension is 1.
  tsunami_lab::t_idx l_nx = 0;
  tsunami_lab::t_idx l_ny = 0;
  int solver;

  // set cell size
  tsunami_lab::t_real l_dxy = 1;

  std::cout << "###################################" << std::endl;
  std::cout << "### Tsunami Lab                 ###" << std::endl;
  std::cout << "###                             ###" << std::endl;
  std::cout << "### http://scalable.uni-jena.de ###" << std::endl;
  std::cout << "###################################" << std::endl;

  if (i_argc != 3) {
    std::cerr << "invalid number of arguments, usage:" << std::endl;
    std::cerr << "  ./build/tsunami_lab N_CELLS_X SOLVER" << std::endl;
    std::cerr << "where N_CELLS_X is the number of cells in x-direction and "
                 "y-direction is computed automatically."
              << std::endl;
    std::cerr << "where SOLVER is the solver of one iteration. Type 0 for "
                 "roe-solver and 1 for fwave solver."
              << std::endl;
    return EXIT_FAILURE;
  } else {
    l_nx = atoi(i_argv[1]);
    if (l_nx < 1) {
      std::cerr << "invalid number of cells" << std::endl;
      return EXIT_FAILURE;
    }

    solver = atoi(i_argv[2]);
    if (!(solver == 0 || solver == 1)) {
      std::cerr << "invalid setup for solver, use 0 or 1 to set up solver"
                << std::endl;
      return EXIT_FAILURE;
    }
  }

  std::cout << "runtime configuration" << std::endl;
  std::cout << "  number of cells in x-direction: " << l_nx << std::endl;
  std::cout << "  number of cells in y-direction: " << l_ny << std::endl;
  std::cout << "  cell size:                      " << l_dxy << std::endl;

  if (solver == 0) {
    std::cout << "  using roe-solver" << std::endl;
  }
  if (solver == 1) {
    std::cout << "  using fwave-solver" << std::endl;
  }

  // construct NetCdf-reader
  tsunami_lab::io::NetCdf *l_netcdf;
  l_netcdf = new tsunami_lab::io::NetCdf(l_nx, "bathymetry_data.nc",
                                         "displacement_data.nc");

  l_ny = l_netcdf->get_amount_y();
  l_dxy = l_netcdf->get_dxy();

  // construct setup
  tsunami_lab::setups::Setup *l_setup;
  l_setup = new tsunami_lab::setups::TsunamiEvent(l_nx);

  // construct solver
  tsunami_lab::patches::WavePropagation *l_waveProp;
  l_waveProp = new tsunami_lab::patches::WavePropagation2d(l_nx, l_ny);

  // maximum observed height in the setup
  tsunami_lab::t_real l_hMax =
      std::numeric_limits<tsunami_lab::t_real>::lowest();

  // set up solver
  for (tsunami_lab::t_idx l_cy = 0; l_cy < l_ny; l_cy++) {
    tsunami_lab::t_real l_y = l_cy * l_dxy;

    for (tsunami_lab::t_idx l_cx = 0; l_cx < l_nx; l_cx++) {
      tsunami_lab::t_real l_x = l_cx * l_dxy;

      // get initial values of the setup
      tsunami_lab::t_real l_h = l_setup->getHeight(l_x, l_y);
      l_hMax = std::max(l_h, l_hMax);

      tsunami_lab::t_real l_hu = l_setup->getMomentumX(l_x, l_y);
      tsunami_lab::t_real l_hv = l_setup->getMomentumY(l_x, l_y);

      tsunami_lab::t_real l_b = l_setup->getBathymetry(l_x, l_y);

      // set initial values in wave propagation solver
      l_waveProp->setHeight(l_cx, l_cy, l_h);

      l_waveProp->setMomentumX(l_cx, l_cy, l_hu);

      l_waveProp->setMomentumY(l_cx, l_cy, l_hv);

      l_waveProp->setBathymetry(l_cx, l_cy, l_b);

      l_waveProp->setReflection(0, false, false);
    }
  }

  // set up time and print control
  tsunami_lab::t_idx l_timeStep = 0;
  tsunami_lab::t_real l_endTime = 50;
  tsunami_lab::t_real l_simTime = 0;

  // initialize the timescaling the momentum is ignored in the first step
  tsunami_lab::t_real l_speedMax = std::sqrt(9.81 * l_hMax);

  // derive constant time step; changes at simulation time are ignored
  tsunami_lab::t_real l_dt = 0.5 * l_dxy / l_speedMax;

  // derive scaling for a time step
  tsunami_lab::t_real l_scaling = l_dt / l_dxy;

  // write bathymetry data
  l_netcdf->writeBathymetry(l_waveProp->getStride(),
                            l_waveProp->getBathymetry());

  std::cout << "entering time loop" << std::endl;
  // iterate over time
  while (l_simTime < l_endTime) {
    if (l_timeStep % 25 == 0) {
      std::cout << "  simulation time / #time steps: " << l_simTime << " / "
                << l_timeStep << std::endl;

      l_netcdf->write(l_waveProp->getStride(), l_waveProp->getHeight(),
                      l_waveProp->getMomentumX(), l_waveProp->getMomentumY(),
                      l_timeStep / 25, l_simTime);
    }

    l_waveProp->setGhostOutflow();
    l_waveProp->timeStep(l_scaling, solver);

    l_timeStep++;
    l_simTime += l_dt;
  }

  std::cout << "finished time loop" << std::endl;

  // free memory
  std::cout << "freeing memory" << std::endl;
  delete l_setup;
  delete l_waveProp;
  delete l_netcdf;

  std::cout << "finished, exiting" << std::endl;
  return EXIT_SUCCESS;
}
