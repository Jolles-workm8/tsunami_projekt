#Tsunami propagation

This is code of the Tsunami lab taught at Friedrich-Schiller-Universität Jena.
Further information is available from: http://scalable.uni-jena.de.

This is a toolkit to compute Tsunami Events based on NetCDF-files or HDF5-files.
It takes into account the Bathymetry of the underlying Seadbed and solves the numerical problem with the shallow water equation.



![](https://github.com/Jolles-workm8/tsunami_projekt/markdown/images/tsunami2.gif)

## Installation

This code can only be run on Linux machines.

To compile this Code ypu need to install several frameworks.

##### Scons
Install Scons manually from https://github.com/SCons/scons or via

    pip install scons

##### LibNetCDF
Install LibNetCDF manually from https://github.com/Unidata/netcdf-c or via

    sudo apt-get install libnetcdf-dev

##### CUDA
To use the CUDA-Kernel you need to install it first. Follow the instructions on https://developer.nvidia.com/cuda-downloads

### Compile
To compile run

    scons

##Running the code


This code relies on netcdf data in the COARS format. Save your bathymetry data as bathymetry_data.nc and your displacement data as displacement_data.nc.

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
