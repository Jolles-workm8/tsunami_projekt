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
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>

#include "io/NetCdf_Read.h"
#include "io/NetCdf_Write.h"
#include "patches/WavePropagation2d.h"
#include "patches/cuda_WavePropagation2d.h"
//#include "setups/ArtificialTsunami.h"
#include "setups/TsunamiEvent.h"

int main(int i_argc, char *i_argv[]) {
  // number of cells in x- and y-direction. Default for y-dimension is 1.
  tsunami_lab::t_idx l_nx = 0;
  tsunami_lab::t_idx l_ny = 0;
  tsunami_lab::t_idx l_rescaleFactor_input = 1;
  tsunami_lab::t_idx l_rescaleFactor_output = 1;
  tsunami_lab::t_idx l_computeSteps = 1;
  tsunami_lab::t_real l_endTime = 100;

  // set cell size
  tsunami_lab::t_real l_dxy = 1;

  std::cout << "###################################" << std::endl;
  std::cout << "### Tsunami Lab                 ###" << std::endl;
  std::cout << "###                             ###" << std::endl;
  std::cout << "### http://scalable.uni-jena.de ###" << std::endl;
  std::cout << "###################################" << std::endl;

  if (i_argc != 5) {
    std::cerr << "invalid number of arguments, usage:" << std::endl;
    return EXIT_FAILURE;
  } else {
    l_rescaleFactor_input = atoi(i_argv[1]);
    if (l_rescaleFactor_input < 1) {
      std::cerr << "invalid input rescaleFactor" << std::endl;
      return EXIT_FAILURE;
    }

    l_rescaleFactor_output = atoi(i_argv[2]);
    if (l_rescaleFactor_output < 1) {
      std::cerr << "invalid output rescaleFactor" << std::endl;
      return EXIT_FAILURE;
    }
    l_endTime = atoi(i_argv[3]);
    if(l_endTime<1){
      std::cerr<<"invalid seconds to compute"<< std::endl;
    }
    l_computeSteps = atoi(i_argv[4]);
    if (l_computeSteps < 1 || l_computeSteps>l_endTime){
      std::cerr<<"invalid computesteps"<< std::endl;
    }
  }

  // construct NetCdf-reader
  tsunami_lab::io::NetCdf_Read *l_netcdf_read;
  l_netcdf_read = new tsunami_lab::io::NetCdf_Read(l_rescaleFactor_input, "bathymetry_data.nc",
                                         "displacement_data.nc");

  l_nx = l_netcdf_read->get_nx();
  l_ny = l_netcdf_read->get_ny();
  l_dxy = l_netcdf_read->get_dxy();

  std::cout << "runtime configuration" << std::endl;
  std::cout << "  number of cells in x-direction: " << l_nx << std::endl;
  std::cout << "  number of cells in y-direction: " << l_ny << std::endl;
  std::cout << "  cell size:                      " << l_dxy << std::endl;

  tsunami_lab::io::NetCdf_Write *l_netcdf_write;
  l_netcdf_write = new tsunami_lab::io::NetCdf_Write(l_nx, l_ny, l_rescaleFactor_output, l_dxy);

  // construct setup
  tsunami_lab::setups::Setup *l_setup;
  l_setup = new tsunami_lab::setups::TsunamiEvent(l_nx, l_netcdf_read);

  // construct solver
  tsunami_lab::patches::WavePropagation *l_waveProp;
  l_waveProp = new tsunami_lab::patches::cuda_WavePropagation2d(l_nx, l_ny);

  // maximum observed height in the setup
  tsunami_lab::t_real l_hMax =
      std::numeric_limits<tsunami_lab::t_real>::lowest();

  std::cout << "start reading setup values " << std::endl;

  using namespace std::chrono;


  // set up solver
  // TODO Parallelize for First touch
  for (tsunami_lab::t_idx l_cy = 0; l_cy < l_ny; l_cy++) {
    //tsunami_lab::t_real l_y = l_cy * l_dxy;

    for (tsunami_lab::t_idx l_cx = 0; l_cx < l_nx; l_cx++) {
      //tsunami_lab::t_real l_x = l_cx * l_dxy;

      // get initial values of the setup
      tsunami_lab::t_real l_h = l_setup->getHeight(l_cx, l_cy);
      l_hMax = std::max(l_h, l_hMax);

      tsunami_lab::t_real l_hu = l_setup->getMomentumX(l_cx, l_cy);
      tsunami_lab::t_real l_hv = l_setup->getMomentumY(l_cx, l_cy);

      tsunami_lab::t_real l_b = l_setup->getBathymetry(l_cx, l_cy);

      // set initial values in wave propagation solver
      l_waveProp->setHeight(l_cx, l_cy, l_h);

      l_waveProp->setMomentumX(l_cx, l_cy, l_hu);

      l_waveProp->setMomentumY(l_cx, l_cy, l_hv);

      l_waveProp->setBathymetry(l_cx, l_cy, l_b);

      l_waveProp->setReflection(0, false, false);
    }
  }

  l_waveProp->MemTransfer();

  // set up time and print control
  tsunami_lab::t_idx l_timeStep = 0;

  tsunami_lab::t_real l_simTime = 0;

  std::cout << "l_hMax " << l_hMax << std::endl;
  // initialize the timescaling the momentum is ignored in the first step
  tsunami_lab::t_real l_speedMax = std::sqrt(9.81 * l_hMax);

  // derive constant time step; changes at simulation time are ignored
  tsunami_lab::t_real l_dt = 0.5 * l_dxy / l_speedMax;
  std::cout << "l_dt " << l_dt << std::endl;
  // derive scaling for a time step
  tsunami_lab::t_real l_scaling = l_dt / l_dxy;

  // write bathymetry data
  l_netcdf_write->writeBathymetry(l_waveProp->getStride(),
                            l_waveProp->getBathymetry());

  std::cout << "entering time loop" << std::endl;
  // iterate over time
  while (l_simTime < l_endTime) {
    
    std::cout << "  simulation time / #time steps: " << l_simTime << " / "
                << l_timeStep << std::endl;
 
    l_netcdf_write->write(l_waveProp->getStride(), l_waveProp->getHeight(),
                      l_waveProp->getMomentumX(), l_waveProp->getMomentumY(),
                      l_timeStep, l_simTime);
    

    l_waveProp->timeStep(l_scaling, l_computeSteps);
    l_timeStep++;
    l_simTime += l_dt * l_computeSteps;
  }


  // free memory
  std::cout << "freeing memory" << std::endl;
  delete l_setup;
  delete l_waveProp;
  delete l_netcdf_write;

  std::cout << "finished, exiting" << std::endl;
  return EXIT_SUCCESS;
}
