





#include "NetCdf_Write.h"

#define SECOND "s"
#define METER "m"
#define METER_PER_SECOND "m/s"
#define ERR(e) \
  { printf("Error: %s\n", nc_strerror(e)); }

tsunami_lab::io::NetCdf_Write::NetCdf_Write(t_idx i_nx, t_idx i_ny, t_idx i_rescaleFactor, t_real l_dxy) {

l_rescaleFactor = i_rescaleFactor;

/////////////////////////////////////////
  /// Prepare writing data into a file ///
  ///////////////////////////////////////

  l_nx_out = (t_idx)(i_nx / l_rescaleFactor);
  l_ny_out = (t_idx)(i_ny / l_rescaleFactor);
  int x_dim, y_dim, time_dim;

  if ((retval = nc_create("solver.nc", NC_CLOBBER, &ncid))) ERR(retval);

  // define the dimensions.
  if ((retval = nc_def_dim(ncid, "x", l_nx_out, &x_dim))) ERR(retval);
  if ((retval = nc_def_dim(ncid, "y", l_ny_out, &y_dim))) ERR(retval);
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
  t_real *l_posX = new t_real[l_nx_out];
  t_real *l_posY = new t_real[l_ny_out];
  for (t_idx l_iy = 0; l_iy < l_ny_out; l_iy++) {
    l_posY[l_iy] = (l_iy + 0.5) * l_dxy * l_rescaleFactor;
  }
  for (t_idx l_ix = 0; l_ix < l_nx_out; l_ix++) {
    l_posX[l_ix] = (l_ix + 0.5) * l_dxy * l_rescaleFactor;
  }

  // write the coordinate variable data
  if ((retval = nc_put_var_float(ncid, w_x_varid, &l_posX[0]))) ERR(retval);
  if ((retval = nc_put_var_float(ncid, w_y_varid, &l_posY[0]))) ERR(retval);

  delete[] l_posX;
  delete[] l_posY;
}

tsunami_lab::io::NetCdf_Write::~NetCdf_Write() {
    if ((retval = nc_close(ncid))) ERR(retval);
}

void tsunami_lab::io::NetCdf_Write::writeArray(t_idx i_stride, t_real const *i_array,
                                        t_idx i_timeStep, int i_varid) {

    size_t start[3], count[3];

    // array count for data to write per time steps
    count[0] = 1;
    count[1] = 1;
    count[2] = 1;
    // array start for position displaceent in dimensions
    start[0] = i_timeStep;
    start[1] = 0;
    start[2] = 0;

    t_real l_cell_Value;
    t_idx l_arrayPos;

    // write the computed data.
    for (t_idx l_ceX = 0; l_ceX < l_nx_out; l_ceX++) {
        for (t_idx l_ceY = 0; l_ceY < l_ny_out; l_ceY++) {

            l_cell_Value = 0;

            //iterate and average over the cells in one output cell
            for (t_idx l_ix = 0; l_ix < l_rescaleFactor; l_ix++) {
                for (t_idx l_iy = 0; l_iy < l_rescaleFactor; l_iy++) {
                    l_arrayPos = (l_ceX * l_rescaleFactor + l_ix) +
                        (l_ceY * l_rescaleFactor + l_iy)* i_stride;
                    l_cell_Value += i_array[l_arrayPos];
                }
            }

            l_cell_Value /= (l_rescaleFactor * l_rescaleFactor);
            start[1] = l_ceY;
            start[2] = l_ceX;

            if ((retval = nc_put_vara_float(ncid, i_varid, start, count, &l_cell_Value))) ERR(retval);

        }
    }
}

void tsunami_lab::io::NetCdf_Write::write(t_idx i_stride, t_real const *i_h,
                                        t_real const *i_hu, t_real const *i_hv,
                                        t_idx i_timeStep, t_real i_simTime) {

    size_t start[1], count[1];

    // array count for data to write per time steps
    count[0] = 1;
    // array start for position displaceent in dimensions
    start[0] = i_timeStep;

    // write time since start
    if ((retval = nc_put_vara_float(ncid, time_varid, start, count, &i_simTime)))
        ERR(retval);

    writeArray(i_stride, i_h, i_timeStep, h_varid);
    writeArray(i_stride, i_hu, i_timeStep, hu_varid);
    writeArray(i_stride, i_hv, i_timeStep, hv_varid);
    
}

void tsunami_lab::io::NetCdf_Write::writeBathymetry(t_idx i_stride,
                                              t_real const *i_b) {


    t_real l_cell_Value;
    t_idx l_arrayPos;
    size_t count[2] = {1,1};
    size_t start[2] = {0,0};
    //iterate over every cell in the output array
    for (t_idx l_ceX = 0; l_ceX < l_nx_out; l_ceX++) {
        for (t_idx l_ceY = 0; l_ceY < l_ny_out; l_ceY++) {

            l_cell_Value = 0;
            //iterate and average over the cells in one output cell
            for (t_idx l_ix = 0; l_ix < l_rescaleFactor; l_ix++) {
                for (t_idx l_iy = 0; l_iy < l_rescaleFactor; l_iy++) {
                    l_arrayPos = (l_ceX * l_rescaleFactor + l_ix) +
                        (l_ceY * l_rescaleFactor + l_iy)* i_stride;
                    l_cell_Value += i_b[l_arrayPos];
                }
            }

            l_cell_Value /= (l_rescaleFactor * l_rescaleFactor);
            start[0] = l_ceY;
            start[1] = l_ceX;
            if ((retval = nc_put_vara_float(ncid, bath_varid, start, count, &l_cell_Value))) ERR(retval);
        }
    }
}

