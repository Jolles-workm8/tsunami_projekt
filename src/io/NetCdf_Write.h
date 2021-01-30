#ifndef TSUNAMI_LAB_IO_NETCDF
#define TSUNAMI_LAB_IO_NETCDF

#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>

#include <cmath>
#include <iostream>
#include <cstring>



#include "../constants.h"

namespace tsunami_lab {
    namespace io {
        class NetCdf_Write;
    }
}  // namespace tsunami_lab

class tsunami_lab::io::NetCdf_Write {
 private:
    t_idx l_rescaleFactor;

    //saves the length of the output array in the specific directions
    t_idx l_nx_out;
    t_idx l_ny_out;

    // saves errors
    int retval;

    //variables for writing
    int ncid;
    int h_varid;
    int hu_varid;
    int hv_varid;
    int time_varid;
    int bath_varid;
    int w_x_varid;
    int w_y_varid;

 public:
    NetCdf_Write(t_idx i_nx, t_idx i_ny, t_idx rescale, t_real l_dxy);

    ~NetCdf_Write();

    void writeArray(t_idx i_stride, t_real const* i_array, t_idx i_timeStep, int i_varid);

    void write(t_idx i_stride, t_real const* i_h, t_real const* i_hu,
             t_real const* i_hv, t_idx i_timeStep, t_real i_simTime);

    void writeBathymetry(t_idx i_stride, t_real const* i_b);

};
#endif