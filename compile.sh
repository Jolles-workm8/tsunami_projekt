#!/bin/bash

module purge 
module load tools/python/3.8 nvidia/cuda/11 compiler/gcc/10.2.0 compiler/intel/2020-Update2

#pip install scons --user

scons

