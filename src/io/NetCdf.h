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
 * Read and Write from and to NetCDF-files
 **/
#ifndef TSUNAMI_LAB_IO_NETCDF
#define TSUNAMI_LAB_IO_NETCDF

#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>

#include <cmath>
#include <cstring>
#include <iostream>
#include <string>

#include "../constants.h"

namespace tsunami_lab {
namespace io {
class NetCdf;
}
}  // namespace tsunami_lab

class tsunami_lab::io::NetCdf {
 private:
  t_idx l_nx;
  t_idx l_ny;
  t_real l_dxy;
  // size of the input data in m
  t_real l_size_x;
  t_real l_size_y;

  int w_x_varid;
  int w_y_varid;

  // variables needed for reading for the bathymetry file
  int r_bath_ncid;
  int r_bath_x_varid;
  int r_bath_y_varid;
  int r_bath_z_varid;
  int r_bath_x_dimid;
  int r_bath_y_dimid;
  size_t r_x_bath_length;
  size_t r_y_bath_length;
  int l_bath_cellsize;

  // variables for min and max in dimension x and y of bathymetry file
  float l_bath_max_value_x;
  float l_bath_min_value_x;
  float l_bath_max_value_y;
  float l_bath_min_value_y;

  // variables needed for reading for the displacement file
  int r_displ_ncid;
  int r_displ_x_varid;
  int r_displ_y_varid;
  int r_displ_z_varid;
  int r_displ_x_dimid;
  int r_displ_y_dimid;
  size_t r_x_displ_length;
  size_t r_y_displ_length;
  int l_displ_cellsize;

  // variables for min and max in dimension x and y of bathymetry file
  float l_displ_max_value_x;
  float l_displ_min_value_x;
  float l_displ_max_value_y;
  float l_displ_min_value_y;

  // variable for rounding the indes so the input matches the data in the
  // dataset
  float scaling_bath_x;
  float scaling_bath_y;
  float scaling_displ_x;
  float scaling_displ_y;

  // saves errors
  int retval;

  int ncid;
  int h_varid, hu_varid, hv_varid, time_varid, bath_varid;

  // variables for reading filename

  // saves a rescaling Value for the output array
  // 1 equals normal
  // 2 eauqls half
  // only values greater or equal than 1 are accepted
  t_idx rescaleFactor;

  //saves the length of the output array in the specific directions
  t_idx l_nx_out;
  t_idx l_ny_out;

 public:
  NetCdf(t_idx i_nx, t_idx i_rescaleFactor, const char* bathymetry_filename,
         const char* displacement_filename);

  ~NetCdf();

  /**
   *
   * @param i_dxy cell width in x- and y-direction.
   * @param i_nx number of cells in x-direction.
   * @param i_ny number of cells in y-direction.
   * @param i_stride stride of the data arrays in y-direction (x is assumed to
   *be stride-1).
   * @param i_h water height of the cells; optional: use nullptr if not
   *required.
   * @param i_hu momentum in x-direction of the cells; optional: use nullptr if
   *not required.
   * @param i_hv momentum in y-direction of the cells; optional: use nullptr if
   *not required.
   **/
  void write(t_idx i_stride, t_real const* i_h, t_real const* i_hu,
             t_real const* i_hv, t_idx i_timeStep, t_real i_simTime);

  void writeBathymetry(t_idx i_stride, t_real const* i_b);

  void writeArray(t_idx i_stride, t_real const* i_array, t_idx i_timeStep, int i_varid);
  /**
   * read the next neighbour bathymetry cell of the  simulation cell x y
   **/
  t_real read_bathymetry(t_idx i_x, t_idx i_y);

  /**
   * read the next neighbour displacement cell of the  simulation cell x y
   * if the cell is outside of displacement data return 0
   **/
  t_real read_displacement(t_idx i_x, t_idx i_y);

  /**
   * returns the number of cells in y direction of the simulation
   **/
  t_idx get_amount_y() { return l_ny; };

  /**
   *returns the size of a cell in the simulation
   **/
  t_real get_dxy() { return l_dxy; };

  /**
   * updates the maximum and minimum of bathymetry in x and y direction
   **/
  void update_max_min_bath();

  /**
   * updates the maximum and minimum of displacement in x and y direction
   **/
  void update_max_min_displ();

  /**
   * update the cellsize of bathymetry
   **/
  void update_bath_cellsize();

  /**
   * update the cellsize of displacement
   **/
  void update_displ_cellsize();

  /**
   * gets the real Position of a cell
   *
   **/
  void getCellPos(t_idx i_x, t_idx i_y, t_real& o_pos_x, t_real& o_pos_y);
};

#endif
