#ifndef TSUNAMI_LAB_IO_NETCDF_READ
#define TSUNAMI_LAB_IO_NETCDF_READ

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
    class NetCdf_Read;
    }
}  // namespace tsunami_lab

class tsunami_lab::io::NetCdf_Read {
 private:
    //Cell Variables
    t_idx l_nx;
    t_idx l_ny;
    t_real l_dxy;

    //bathymetry array
    t_real *l_b = nullptr;

    //displacemet array
    t_real *l_d = nullptr;

    //for rescaling input Data
    t_idx rescaleFactor;

    // saves errors
    int retval;

    // variables needed for reading for the bathymetry file
    int r_bath_ncid;
    int r_bath_x_varid;
    int r_bath_y_varid;
    int r_bath_z_varid;
    int r_bath_x_dimid;
    int r_bath_y_dimid;
    size_t r_x_bath_length;
    size_t r_y_bath_length;
    float l_bath_cellsize;

    // variables for min and max in dimension x and y of bathymetry files
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
    float l_displ_cellsize;

    // variables for min and max in dimension x and y of bathymetry file
    float l_displ_max_value_x;
    float l_displ_min_value_x;
    float l_displ_max_value_y;
    float l_displ_min_value_y;

 public:
    NetCdf_Read(t_idx rescale, const char* bathymetry_filename,
         const char* displacement_filename);

    ~NetCdf_Read();
  
    
    /**
     * read complete bathymetry file and write to local array i_b
    **/
    void read_bathymetry(t_real *o_d);

    /**
     * read complete bathymetry file and write local to array i_d
    **/
    void read_displacement(t_real *o_d);   
    /**
     * 
    **/ 
    t_idx get_nx(){
        return l_nx;
    }
    /**
     * 
    **/ 
    t_idx get_ny(){
        return l_ny;
    }
    /**
     * 
    **/ 
    t_idx get_dxy(){
        return l_dxy;
    }

    t_real get_i_b(t_idx i_x, t_idx i_y){
        if(l_b[i_x+i_y*l_nx] == l_b[i_x+i_y*l_nx]){
            return l_b[i_x+i_y*l_nx];
        }
        return 0;
    }

    t_real get_i_d(t_idx i_x, t_idx i_y){
        if(l_d[i_x+i_y*l_nx] == l_d[i_x+i_y*l_nx]){
            return l_d[i_x+i_y*l_nx];
        }
        return 0;
    }
};
#endif