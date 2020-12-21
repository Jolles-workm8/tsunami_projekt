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
#include "NetCdf.h"

#define SECOND "s"
#define METER "m"
#define METER_PER_SECOND "m/s"
#define ERR(e) \
  { printf("Error: %s\n", nc_strerror(e)); }

tsunami_lab::io::NetCdf::NetCdf(t_idx i_nx, const char *bathymetry_filename,
                                const char *displacement_filename) {
  l_nx = i_nx;

  ////////////////////////////////////////////
  /// Prepare reading files from a source ///
  //////////////////////////////////////////
  // open the files we want to read.
  if ((retval = nc_open(bathymetry_filename, NC_NOWRITE, &r_bath_ncid)))
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
  if ((retval = nc_inq_dimlen(r_bath_ncid, r_bath_x_dimid, &r_x_bath_length)))
    ERR(retval);
  if ((retval = nc_inq_dimlen(r_bath_ncid, r_bath_y_dimid, &r_y_bath_length)))
    ERR(retval);

  // update the min and max value of the field in x and y direction
  update_max_min_bath();
  update_bath_cellsize();


  // Compute size_x
  l_size_x = l_bath_max_value_x + l_bath_cellsize - l_bath_min_value_x;

  // Compute size_y
  l_size_y = l_bath_max_value_y + l_bath_cellsize - l_bath_min_value_y;

  // calculate cellsize
  l_dxy = l_size_x / l_nx;
  // calculate number of cells in y direction and round
  l_ny = (tsunami_lab::t_idx)(l_size_y / l_dxy + 0.5);

  scaling_bath_x = r_x_bath_length / l_nx;
  scaling_bath_y = r_y_bath_length / l_ny;

  // DISPLACEMENT FILE //
  // open the file we want to read
  if ((retval = nc_open(displacement_filename, NC_NOWRITE, &r_displ_ncid)))
    ERR(retval);

  // Get the variable id's of the x, y and z coordinates
  if ((retval = nc_inq_varid(r_displ_ncid, "x", &r_displ_x_varid))) ERR(retval);
  if ((retval = nc_inq_varid(r_displ_ncid, "y", &r_displ_y_varid))) ERR(retval);
  if ((retval = nc_inq_varid(r_displ_ncid, "z", &r_displ_z_varid))) ERR(retval);

  // get dim id of x and y
  if ((retval = nc_inq_dimid(r_displ_ncid, "x", &r_displ_x_dimid))) ERR(retval);
  if ((retval = nc_inq_dimid(r_displ_ncid, "y", &r_displ_y_dimid))) ERR(retval);

  // get the length in x-direction an y-direction
  if ((retval =
           nc_inq_dimlen(r_displ_ncid, r_displ_x_dimid, &r_x_displ_length)))
    ERR(retval);
  if ((retval =
           nc_inq_dimlen(r_displ_ncid, r_displ_y_dimid, &r_y_displ_length)))
    ERR(retval);


  // update the min and max value of the field in x and y direction
  update_max_min_displ();
  update_displ_cellsize();


  /////////////////////////////////////////
  /// Prepare writing data into a file ///
  ///////////////////////////////////////

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
  dims[1] = y_dim;
  dims[2] = x_dim;

  // define other 3 dim variables
  if ((retval = nc_def_var(ncid, "height", NC_FLOAT, 3, dims, &h_varid)))
    ERR(retval);
  if ((retval = nc_def_var(ncid, "momentum_x", NC_FLOAT, 3, dims, &hu_varid)))
    ERR(retval);
  if ((retval = nc_def_var(ncid, "momentum_y", NC_FLOAT, 3, dims, &hv_varid)))
    ERR(retval);

  // dim array for bathymetry dimensions
  int l_dimBath[2];
  l_dimBath[0] = y_dim;
  l_dimBath[1] = x_dim;

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
  t_real *l_posX = new t_real[l_nx];
  t_real *l_posY = new t_real[l_ny];
  for (t_idx l_iy = 0; l_iy < l_ny; l_iy++) {
    l_posY[l_iy] = (l_iy + 0.5) * l_dxy;
  }
  for (t_idx l_ix = 0; l_ix < l_nx; l_ix++) {
    l_posX[l_ix] = (l_ix + 0.5) * l_dxy;
  }

  // write the coordinate variable data
  if ((retval = nc_put_var_float(ncid, w_x_varid, &l_posX[0]))) ERR(retval);
  if ((retval = nc_put_var_float(ncid, w_y_varid, &l_posY[0]))) ERR(retval);

  delete[] l_posX;
  delete[] l_posY;
}

tsunami_lab::io::NetCdf::~NetCdf() {
  int retval;
  // close the files
  if ((retval = nc_close(ncid))) ERR(retval);
  if ((retval = nc_close(r_displ_ncid))) ERR(retval);
  if ((retval = nc_close(r_bath_ncid))) ERR(retval);
}


void tsunami_lab::io::NetCdf::writeBathymetry(t_idx i_stride,
                                              t_real const *i_b) {
  size_t count[2] = {l_ny,l_nx};
  size_t start[2] = {0,0};
  ptrdiff_t imap[2] = {(ptrdiff_t)i_stride,1};
  if ((retval = nc_put_varm_float(ncid, bath_varid, start, count, NULL, imap, &i_b[0]))) ERR(retval);
}

void tsunami_lab::io::NetCdf::write(t_idx i_stride, t_real const *i_h,
                                    t_real const *i_hu, t_real const *i_hv,
                                    t_idx i_timeStep, t_real i_simTime) {

  size_t start[3], count[3];

  // array count for data to write per time steps
  count[0] = 1;
  count[1] = l_ny;
  count[2] = l_nx;
  // array start for position displaceent in dimensions
  start[0] = i_timeStep;
  start[1] = 0;
  start[2] = 0;

  //array for mapping the values of calculated arrays whitout ghost cells
  ptrdiff_t imap[3] = {1,(ptrdiff_t)i_stride,1};
  // write time since start
  if ((retval = nc_put_vara_float(ncid, time_varid, start, count, &i_simTime)))
    ERR(retval);

  // write the computed data.
  if ((retval = nc_put_varm_float(ncid, h_varid, start, count, NULL, imap, &i_h[0])))
    ERR(retval);
  if ((retval = nc_put_varm_float(ncid, hu_varid, start, count, NULL, imap, &i_hu[0])))
    ERR(retval);
  if ((retval = nc_put_varm_float(ncid, hv_varid, start, count, NULL, imap, &i_hv[0])))
    ERR(retval);
  std::cout << i_simTime << std::endl;
}

tsunami_lab::t_real tsunami_lab::io::NetCdf::read_bathymetry(t_idx i_x,
                                                             t_idx i_y) {
  float bath_return_value;
  size_t index[2];
  t_real o_pos_x;
  t_real o_pos_y;
  getCellPos(i_x, i_y, o_pos_x, o_pos_y);
  index[1] = (size_t)(((o_pos_x- l_bath_min_value_x)/l_bath_cellsize)+0.5);
  index[0] = (size_t)(((o_pos_y- l_bath_min_value_y)/l_bath_cellsize)+0.5);
  if ((retval = nc_get_var1_float(r_bath_ncid, r_bath_z_varid, index,
                                  &bath_return_value)))
    ERR(retval);
  return (t_real)bath_return_value;
}

tsunami_lab::t_real tsunami_lab::io::NetCdf::read_displacement(t_idx i_x,
                                                               t_idx i_y) {
  float displ_return_value;
  size_t index[2];
  t_real o_pos_x;
  t_real o_pos_y;
  getCellPos(i_x, i_y, o_pos_x, o_pos_y);


  if(o_pos_x > l_displ_min_value_x - 0.5 *l_displ_cellsize &&
      o_pos_x < l_displ_max_value_x + 0.5 *l_displ_cellsize &&
      o_pos_y > l_displ_min_value_y - 0.5 *l_displ_cellsize &&
      o_pos_y < l_displ_max_value_y + 0.5 *l_displ_cellsize){


    index[1] = (size_t)(((o_pos_x- l_displ_min_value_x)/l_displ_cellsize));
    index[0] = (size_t)(((o_pos_y- l_displ_min_value_y)/l_displ_cellsize));
    if ((retval = nc_get_var1_float(r_displ_ncid, r_displ_z_varid, index,
                                &displ_return_value)))
      ERR(retval);
      if(displ_return_value != displ_return_value){
        std::cout << "Not a Number in displacement imput, using displacement 0" << std::endl;

          return 0;
      }
    return (t_real)displ_return_value;
  }
  else{
    return 0;
  }

}


void tsunami_lab::io::NetCdf::update_max_min_bath(){
  size_t index_min[1];
  size_t index_max[1];
  index_min[0] = 0;

  // read x direction
  index_max[0] = r_x_bath_length - 1;

  if ((retval = nc_get_var1_float(r_bath_ncid, r_bath_x_varid, index_max,
                                  &l_bath_max_value_x)))
    ERR(retval);
  if ((retval = nc_get_var1_float(r_bath_ncid, r_bath_x_varid, index_min,
                                  &l_bath_min_value_x)))
    ERR(retval);


  // read y direction
  index_max[0] = r_y_bath_length - 1;

  if ((retval = nc_get_var1_float(r_bath_ncid, r_bath_y_varid, index_max,
                                  &l_bath_max_value_y)))
    ERR(retval);
  if ((retval = nc_get_var1_float(r_bath_ncid, r_bath_y_varid, index_min,
                                  &l_bath_min_value_y)))
    ERR(retval);

}


void tsunami_lab::io::NetCdf::update_max_min_displ(){
  size_t index_min[1];
  size_t index_max[1];
  index_min[0] = 0;

  // read x direction
  index_max[0] = r_x_displ_length - 1;

  if ((retval = nc_get_var1_float(r_displ_ncid, r_displ_x_varid, index_max,
                                  &l_displ_max_value_x)))
    ERR(retval);
  if ((retval = nc_get_var1_float(r_displ_ncid, r_displ_x_varid, index_min,
                                  &l_displ_min_value_x)))
    ERR(retval);



  // read y direction
  index_max[0] = r_y_displ_length - 1;

  if ((retval = nc_get_var1_float(r_displ_ncid, r_displ_y_varid, index_max,
                                  &l_displ_max_value_y)))
    ERR(retval);
  if ((retval = nc_get_var1_float(r_displ_ncid, r_displ_y_varid, index_min,
                                  &l_displ_min_value_y)))
    ERR(retval);
}


void tsunami_lab::io::NetCdf::update_bath_cellsize(){
  l_bath_cellsize = (l_bath_max_value_x -l_bath_min_value_x) / (r_x_bath_length-1);
}


void tsunami_lab::io::NetCdf::update_displ_cellsize(){
  l_displ_cellsize = (l_displ_max_value_x -l_displ_min_value_x) / (r_x_displ_length - 1);
}


void tsunami_lab::io::NetCdf::getCellPos(t_idx i_x, t_idx i_y,
                                          t_real &o_pos_x, t_real &o_pos_y){
  o_pos_x= i_x + l_bath_min_value_x+0.5*l_dxy;
  o_pos_y= i_y + l_bath_min_value_y+0.5*l_dxy;

}
