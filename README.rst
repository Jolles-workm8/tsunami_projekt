###########
Tsunami Lab
###########

This is code of the Tsunami lab taught at Friedrich-Schiller-Universität Jena.
Further information is available from: http://scalable.uni-jena.de.


Installation
============

This code is build using the autotool scons. If you havn't installed scons yet run
* item text::

        pip install scons

Further you need Libnetcdf. Install it with the repository https://github.com/Unidata/netcdf-c
or with a Linux machine via:
* item text::


        sudo apt-get install libnetcdf-dev


For building and installing the code run.
* item text::

        scons

Running the code
================


This code relies on netcdf data in the COARS format. Save your bathymetry data as bathymetry_data.nc and your displacement data as displacement_data.nc.
* item text::

        mkdir data
        cd data
        cp bathymetry_data.nc ./$LOCATION_REPO$/data
        cp displacement_data.nc ./$LOCATION_REPO$/data
        ../tsunami_lab x y z w

where is:
-x scaling of input data
-y scaling of output data
-z simulation time
-w steps in which you want to have an output
