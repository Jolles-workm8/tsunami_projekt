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
 * IO-routines for writing a snapshot as Comma Separated Values (CSV).
 **/
#include "NetCdf.h"

#define SECOND "s"
#define METER "m"
#define METER_PER_SECOND "m/s"
#define ERR(e) \
  { printf("Error: %s\n", nc_strerror(e)); }

tsunami_lab::io::NetCdf::NetCdf(t_idx i_nx, t_idx i_ny, t_real i_dxy,
                                const char *bathymetry_filename,
                                const char *displacement_filename) {
  /////////////////////////////////////////
  /// Prepare writing data into a file ///
  ///////////////////////////////////////

  l_nx = i_nx;
  l_ny = i_ny;

  int x_dim, y_dim, time_dim;

  if ((retval = nc_create("solver.nc", NC_CLOBBER, &ncid))) ERR(retval);

  // define the dimensions.
  if ((retval = nc_def_dim(ncid, "x", l_nx, &x_dim))) ERR(retval);
  if ((retval = nc_def_dim(ncid, "y", l_ny, &y_dim))) ERR(retval);
  if ((retval = nc_def_dim(ncid, "seconds since", NC_UNLIMITED, &time_dim)))
    ERR(retval);

  // define coordinates
  if ((retval = nc_def_var(ncid, "x", NC_FLOAT, 1, &x_dim, &w_x_varid)))
    ERR(retval);
  if ((retval = nc_def_var(ncid, "y", NC_FLOAT, 1, &y_dim, &w_y_varid)))
    ERR(retval);
  if ((retval = nc_def_var(ncid, "seconds since", NC_FLOAT, 1, &time_dim,
                           &time_varid)))
    ERR(retval);

  // assign units attributes to coordinate and time variables
  if ((retval =
           nc_put_att_text(ncid, w_x_varid, "units", strlen(METER), METER)))
    ERR(retval);
  if ((retval =
           nc_put_att_text(ncid, w_y_varid, "units", strlen(METER), METER)))
    ERR(retval);
  if ((retval =
           nc_put_att_text(ncid, time_varid, "units", strlen(SECOND), SECOND)))
    ERR(retval);

  // dims array is used to pass dimension of variables
  int dims[3];
  dims[0] = time_dim;
  dims[1] = x_dim;
  dims[2] = y_dim;

  // define other 3 dim variables
  if ((retval = nc_def_var(ncid, "height", NC_FLOAT, 3, dims, &h_varid)))
    ERR(retval);
  if ((retval = nc_def_var(ncid, "momentum_x", NC_FLOAT, 3, dims, &hu_varid)))
    ERR(retval);
  if ((retval = nc_def_var(ncid, "momentum_y", NC_FLOAT, 3, dims, &hv_varid)))
    ERR(retval);

  // dim array for bathymetry dimensions
  int l_dimBath[2];
  l_dimBath[0] = x_dim;
  l_dimBath[1] = y_dim;

  // define 2 dim bathymetry
  if ((retval =
           nc_def_var(ncid, "bathymetry", NC_FLOAT, 2, l_dimBath, &bath_varid)))
    ERR(retval);

  // assign units attributes to other variables
  if ((retval = nc_put_att_text(ncid, h_varid, "units", strlen(METER), METER)))
    ERR(retval);
  if ((retval = nc_put_att_text(ncid, hu_varid, "units",
                                strlen(METER_PER_SECOND), METER_PER_SECOND)))
    ERR(retval);
  if ((retval = nc_put_att_text(ncid, hv_varid, "units",
                                strlen(METER_PER_SECOND), METER_PER_SECOND)))
    ERR(retval);
  if ((retval =
           nc_put_att_text(ncid, bath_varid, "units", strlen(METER), METER)))
    ERR(retval);

  // end define mode
  if ((retval = nc_enddef(ncid))) ERR(retval);

  // derive coordinates of cell center
  t_real *l_posX = new t_real[i_nx];
  t_real *l_posY = new t_real[i_ny];
  for (t_idx l_iy = 0; l_iy < i_ny; l_iy++) {
    for (t_idx l_ix = 0; l_ix < i_nx; l_ix++) {
      l_posX[l_ix] = (l_ix + 0.5) * i_dxy;
      l_posY[l_iy] = (l_iy + 0.5) * i_dxy;
    }
  }

  // write the coordinate variable data
  if ((retval = nc_put_var_float(ncid, w_x_varid, &l_posX[0]))) ERR(retval);
  if ((retval = nc_put_var_float(ncid, w_y_varid, &l_posY[0]))) ERR(retval);

  delete[] l_posX;
  delete[] l_posY;

  ////////////////////////////////////////////
  /// Prepare reading files from a source ///
  //////////////////////////////////////////
  // we assume that bathymetry and displacement have the same dimension size.
  // open the file we want to read.
  if ((retval = nc_open(bathymetry_filename, NC_NOWRITE, &r_bath_ncid)))
    ERR(retval);
  if ((retval = nc_open(displacement_filename, NC_NOWRITE, &r_displ_ncid)))
    ERR(retval);

  // BATHYMETRY FILE //
  // Get the variable id's of the x, y an z coordinates

  if ((retval = nc_inq_varid(r_bath_ncid, "x", &r_bath_x_varid))) ERR(retval);
  if ((retval = nc_inq_varid(r_bath_ncid, "y", &r_bath_y_varid))) ERR(retval);
  if ((retval = nc_inq_varid(r_bath_ncid, "z", &r_bath_z_varid))) ERR(retval);

  // get dim id of x and y
  if ((retval = nc_inq_dimid(r_bath_ncid, "x", &r_bath_x_dimid))) ERR(retval);
  if ((retval = nc_inq_dimid(r_bath_ncid, "y", &r_bath_y_dimid))) ERR(retval);

  // get the length in x-direction an y-direction
  if ((retval = nc_inq_dimlen(r_bath_ncid, r_bath_x_dimid, &r_x_length)))
    ERR(retval);
  if ((retval = nc_inq_dimlen(r_bath_ncid, r_bath_y_dimid, &r_y_length)))
    ERR(retval);

  // DISPLACEMENT FILE //
  // Get the variable id's of the x, y and z coordinates

  if ((retval = nc_inq_varid(r_displ_ncid, "x", &r_displ_x_varid))) ERR(retval);
  if ((retval = nc_inq_varid(r_displ_ncid, "y", &r_displ_y_varid))) ERR(retval);
  if ((retval = nc_inq_varid(r_displ_ncid, "z", &r_displ_z_varid))) ERR(retval);

  // allocate arrays for the data of x, y and the bathymetry;

  x_data = new int[r_x_length];
  y_data = new int[r_y_length];

  if ((retval = nc_get_var_int(r_bath_ncid, r_bath_x_varid, &x_data[0])))
    ERR(retval);
  if ((retval = nc_get_var_int(r_bath_ncid, r_bath_y_varid, &y_data[0])))
    ERR(retval);

  scaling_x = r_x_length / l_nx;
  scaling_y = r_y_length / l_ny;
}

tsunami_lab::io::NetCdf::~NetCdf() {
  int retval;
  // close the file
  if ((retval = nc_close(ncid))) ERR(retval);
  if ((retval = nc_close(ncid))) ERR(retval);

  delete[] x_data;
  delete[] y_data;
}
void tsunami_lab::io::NetCdf::writeBathymetry(t_idx i_stride,
                                              t_real const *i_b) {
  t_real *l_b = new t_real[l_nx * l_ny];

  // iterate over all cells
  for (t_idx l_iy = 0; l_iy < l_ny; l_iy++) {
    for (t_idx l_ix = 0; l_ix < l_nx; l_ix++) {
      t_idx l_id_from = l_iy * i_stride + l_ix;
      t_idx l_id_to = l_iy * l_nx + l_ix;
      l_b[l_id_to] = i_b[l_id_from];
    }
  }
  // write bathymetry data
  if ((retval = nc_put_var_float(ncid, bath_varid, &l_b[0]))) ERR(retval);

  delete[] l_b;
}

void tsunami_lab::io::NetCdf::write(t_idx i_stride, t_real const *i_h,
                                    t_real const *i_hu, t_real const *i_hv,
                                    t_idx i_timeStep, t_real i_simTime) {
  // arrays from which data is writen
  t_real *l_h = new t_real[l_nx * l_ny];
  t_real *l_hu = new t_real[l_nx * l_ny];
  t_real *l_hv = new t_real[l_nx * l_ny];

  // iterate over all cells
  for (t_idx l_iy = 0; l_iy < l_ny; l_iy++) {
    for (t_idx l_ix = 0; l_ix < l_nx; l_ix++) {
      t_idx l_id_from = l_iy * i_stride + l_ix;
      t_idx l_id_to = l_iy * l_nx + l_ix;
      l_h[l_id_to] = i_h[l_id_from];
      l_hu[l_id_to] = i_hu[l_id_from];
      l_hv[l_id_to] = i_hv[l_id_from];
    }
  }

  size_t start[3], count[3];

  // array count for data to write per time steps
  count[0] = 1;
  count[1] = l_nx;
  count[2] = l_ny;
  // array start for position displaceent in dimensions
  start[0] = i_timeStep;
  start[1] = 0;
  start[2] = 0;

  // write time since start
  if ((retval = nc_put_vara_float(ncid, time_varid, start, count, &i_simTime)))
    ERR(retval);

  // write the computed data.
  if ((retval = nc_put_vara_float(ncid, h_varid, start, count, &l_h[0])))
    ERR(retval);
  if ((retval = nc_put_vara_float(ncid, hu_varid, start, count, &l_hu[0])))
    ERR(retval);
  if ((retval = nc_put_vara_float(ncid, hv_varid, start, count, &l_hv[0])))
    ERR(retval);
  std::cout << i_simTime << std::endl;

  delete[] l_h;
  delete[] l_hu;
  delete[] l_hv;
}

tsunami_lab::t_real tsunami_lab::io::NetCdf::read_bathymetry(t_idx i_x,
                                                             t_idx i_y) {
  float bath_return_value;
  size_t index[2];
  index[0] = scaling_x * i_x;
  index[1] = scaling_y * i_y;
  if ((retval = nc_get_var1_float(r_bath_ncid, r_bath_z_varid, index,
                                  &bath_return_value)))
    ERR(retval);
  return (t_real)bath_return_value;
}

tsunami_lab::t_real tsunami_lab::io::NetCdf::read_displacement(t_idx i_x,
                                                               t_idx i_y) {
  float displ_return_value;
  size_t index[2];
  index[0] = scaling_x * i_x;
  index[1] = scaling_y * i_y;
  if ((retval = nc_get_var1_float(r_displ_ncid, r_displ_z_varid, index,
                                  &displ_return_value)))
    ERR(retval);
  return (t_real)displ_return_value;
}
