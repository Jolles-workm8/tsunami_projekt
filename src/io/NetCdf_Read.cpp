#include "NetCdf_Read.h"
#define ERR(e) \
  { printf("Error: %s\n", nc_strerror(e)); }

tsunami_lab::io::NetCdf_Read::NetCdf_Read(t_idx rescale,const char *bathymetry_filename,
                                const char *displacement_filename) {
                                    

    rescaleFactor =  rescale;                               

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

    // update the min and max value of the bathemetry in x and y direction
    size_t index_min = 0;
    size_t index_max = r_x_bath_length - 1;

    //x direction
    if ((retval = nc_get_var1_float(r_bath_ncid, r_bath_x_varid, &index_max,
                                    &l_bath_max_value_x)))
        ERR(retval);
    if ((retval = nc_get_var1_float(r_bath_ncid, r_bath_x_varid, &index_min,
                                    &l_bath_min_value_x)))
        ERR(retval);

    //y direction
    index_max = r_y_bath_length - 1;
    if ((retval = nc_get_var1_float(r_bath_ncid, r_bath_y_varid, &index_max,
                                    &l_bath_max_value_y)))
        ERR(retval);
    if ((retval = nc_get_var1_float(r_bath_ncid, r_bath_y_varid, &index_min,
                                    &l_bath_min_value_y)))
        ERR(retval);

    //calculate bathymetry cellSize
    l_bath_cellsize = (l_bath_max_value_x -l_bath_min_value_x) / (r_x_bath_length-1);

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
    index_max = r_x_displ_length - 1;

    if ((retval = nc_get_var1_float(r_displ_ncid, r_displ_x_varid, &index_max,
                                  &l_displ_max_value_x)))
        ERR(retval);
    if ((retval = nc_get_var1_float(r_displ_ncid, r_displ_x_varid, &index_min,
                                  &l_displ_min_value_x)))
        ERR(retval);

    index_max = r_y_displ_length - 1;

    if ((retval = nc_get_var1_float(r_displ_ncid, r_displ_y_varid, &index_max,
                                  &l_displ_max_value_y)))
        ERR(retval);
    if ((retval = nc_get_var1_float(r_displ_ncid, r_displ_y_varid, &index_min,
                                  &l_displ_min_value_y)))
        ERR(retval);

    //calculate displacement cellSize
    l_displ_cellsize = (l_displ_max_value_x -l_displ_min_value_x) / (r_x_displ_length - 1);


    l_b = new t_real[r_x_bath_length * r_y_bath_length];
    l_d = new t_real[r_x_bath_length * r_y_bath_length];

    
    read_bathymetry(l_b);
    read_displacement(l_d);

    l_dxy = l_bath_cellsize*rescaleFactor;
    l_nx = r_x_bath_length/rescaleFactor;
    l_ny = r_y_bath_length/rescaleFactor;

    if(rescaleFactor != 1){
        //save data temporarily
        t_real *l_b_temp = l_b;
        t_real *l_d_temp = l_d;

        //new smaller array
        l_b = new t_real[l_nx * l_ny];
        l_d = new t_real[l_nx * l_ny];

        for(size_t l_ceY = 0; l_ceY < l_ny; l_ceY++) {
            for (size_t l_ceX = 0; l_ceX < l_nx; l_ceX++) {
               l_b[l_ceX + l_ceY * l_nx] = l_b_temp[l_ceX*rescaleFactor + l_ceY*rescaleFactor * r_x_bath_length];
               l_d[l_ceX + l_ceY * l_nx] = l_d_temp[l_ceX*rescaleFactor + l_ceY*rescaleFactor * r_x_bath_length];
            }
        }


        delete[] l_b_temp;
        delete[] l_d_temp;

    }
}





void tsunami_lab::io::NetCdf_Read::read_bathymetry(t_real *o_b){    
        std::cout << "reading Bathymetry Data" << std::endl;
        /*for(size_t l_ce= 0; l_ce < r_x_bath_length * r_y_bath_length; l_ce++){
            i_b[l_ce] = 0;
        }
        */
        if ((retval = nc_get_var(r_bath_ncid, r_bath_z_varid, o_b)))
            ERR(retval);

        if ((retval = nc_close(r_bath_ncid))) ERR(retval); 

}

void tsunami_lab::io::NetCdf_Read::read_displacement(t_real *o_d){
    std::cout << "reading Displacement Data" << std::endl;

    //if dsipl array and bath array don't have the same size rescale displ array to bath array size
    if(l_displ_cellsize  == l_bath_cellsize && r_x_bath_length == r_x_displ_length && r_y_bath_length == r_y_displ_length) {
        if ((retval = nc_get_var_float(r_displ_ncid, r_displ_z_varid, o_d)))
            ERR(retval);
        if ((retval = nc_close(r_displ_ncid))) ERR(retval);
    }
    else{
        float *i_d_temp = new float[r_x_displ_length * r_y_displ_length];
        if ((retval = nc_get_var_float(r_displ_ncid, r_displ_z_varid, i_d_temp)))
            ERR(retval);
        if ((retval = nc_close(r_displ_ncid))) ERR(retval);


        //rescale to bathymetry size
        for(size_t l_ceY = 0; l_ceY < r_y_bath_length; l_ceY++) {
            for (size_t l_ceX = 0; l_ceX < r_x_bath_length; l_ceX++) {
                int posNew = l_ceX + l_ceY * r_x_bath_length;

                int posOldX = (l_bath_min_value_x + l_bath_cellsize*l_ceX - l_displ_min_value_x)/ l_displ_cellsize + 0.5;
                int posOldY = (l_bath_min_value_y + l_bath_cellsize*l_ceY- l_displ_min_value_y)/ l_displ_cellsize + 0.5;

                if(posOldX >= 0 && posOldX < (int) r_x_displ_length && posOldY >= 0 && posOldY < (int) r_y_displ_length){
                    int posOld = posOldX + posOldY * r_x_displ_length;
                    o_d[posNew] = i_d_temp[posOld];
                }
                else{
                    o_d[posNew] = 0; 
                }                
            }
        }
        delete[] i_d_temp;
    }
}
